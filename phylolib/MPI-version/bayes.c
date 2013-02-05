
#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "axml.h"

extern int Thorough;
extern int optimizeRateCategoryInvocations;
extern infoList iList;
extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern char run_id[128];
extern double masterTime;
extern double accumulatedTime;

extern partitionLengths pLengths[MAX_MODEL];

typedef enum
{
  SPR,
  stNNI,
  UPDATE_ALL_BL,
  UPDATE_MODEL,
  UPDATE_GAMMA,
}prop;

typedef struct {
  
  /* these 3 are independent of the state, can be taken out unless we want to pass a single pointer as an argument*/  
  nodeptr * list; /* list of possible re-insertion nodes */ 
  int maxradius;  /* maximum radius of re-insertion from the pruning point */
  tree * tr;
  
  double curprior;
  double newprior;
  double hastings;
  
  /*
   * these are the bits that are necessary for the topo moves
   */
  nodeptr p; /* node pruned */
  nodeptr nb;   /* p->next->back describes an edge that dissapears when p is pruned */
  double nbz[NUM_BRANCHES];
  nodeptr nnb; /* p->next->next->back describes an edge that dissapears when p is pruned */
  double nnbz[NUM_BRANCHES];
  nodeptr r; /* edge neighbour of re-insertion node, q->back */
  nodeptr q; /* re-insertion node */
  /*
   * these are the bits that are necessary for the NNI moves
   */
  int whichNNI; /* 0 same topo, 1, 2 */

  /*
   * necessary for the branch length moves
   */
  double qz[NUM_BRANCHES]; /* BL values prior to re-insertion */
  double bl_sliding_window_w;
  double bl_prior;
  double bl_prior_exp_lambda;
  /*
   * necessary for model moves
   */
  analdef * adef;
  int model;
  int nstates;
  int numSubsRates;
  double rt_sliding_window_w;
  double *curSubsRates;//used for resetting

  /*
   * necessary for gamma
   */
  double curAlpha;
  double gm_sliding_window_w;

} state;

static void print_state(state *s, double startLH)
{
   assert(startLH == s->tr->startLH);
   printBothOpen("tr LH %f, startLH %f, incr %f\n", s->tr->likelihood, startLH, s->tr->likelihood - startLH);
   printBothOpen("pruned %db%d nb %d, nnb %d, reinsert %db%d \n", s->p->number, s->p->back->number,
                                                    s->nb->number, s->nnb->number,
                                                    s->q->number, s->r->number);
      
}

static state *state_init(tree *tr, analdef * adef, int maxradius, double bl_w, double rt_w, double gm_w, double bl_p)
{
  state *curstate  =(state *)malloc(sizeof(state));
  nodeptr *list = (nodeptr *)malloc(sizeof(nodeptr) * 2 * tr->mxtips);
  curstate->list = list;
  curstate->maxradius = maxradius;
  curstate->tr = tr;
  curstate->bl_sliding_window_w = bl_w;
  curstate->bl_prior = 1.0;
  curstate->bl_prior_exp_lambda = bl_p;
  //this can be extended to more than one partition, but one for now
  curstate->model = 0;
  curstate->adef = adef;
  curstate->rt_sliding_window_w = rt_w;
  curstate->nstates = tr->partitionData[curstate->model].states; /* 4 for DNA */
  curstate->numSubsRates = (curstate->nstates * curstate->nstates - curstate->nstates) / 2; /* 6 for DNA */
  curstate->curSubsRates = (double *) malloc(curstate->numSubsRates * sizeof(double));
  curstate->gm_sliding_window_w = gm_w;
  assert(curstate != NULL);
  return curstate;
}
static void state_free(state *s)
{
  assert(s != NULL);
  free(s->list);
  free(s->curSubsRates);
  free(s);
}

static char *Tree2StringRecomREC(char *treestr, tree *tr, nodeptr q, boolean printBranchLengths)
{
  char  *nameptr;            
  double z;
  nodeptr p = q;
  int slot;

  if(isTip(p->number, tr->mxtips)) 
  {	       	  
    nameptr = tr->nameList[p->number];     
    sprintf(treestr, "%s", nameptr);
    while (*treestr) treestr++;
  }
  else 
  {                 	 
    while(!p->x)
      p = p->next;
    *treestr++ = '(';
    treestr = Tree2StringRecomREC(treestr, tr, q->next->back, printBranchLengths);
    *treestr++ = ',';
    treestr = Tree2StringRecomREC(treestr, tr, q->next->next->back, printBranchLengths);
    if(q == tr->start->back) 
    {
      *treestr++ = ',';
      treestr = Tree2StringRecomREC(treestr, tr, q->back, printBranchLengths);
    }
    *treestr++ = ')';                    
    // write innernode as nodenum_b_nodenumback
    sprintf(treestr, "%d", q->number);
    while (*treestr) treestr++;
    *treestr++ = 'b';                    
    sprintf(treestr, "%d", p->back->number);
    while (*treestr) treestr++;
  }

  if(q == tr->start->back) 
  {	      	 
    if(printBranchLengths)
      sprintf(treestr, ":0.0;\n");
    else
      sprintf(treestr, ";\n");	 	  	
  }
  else 
  {                   
    if(printBranchLengths)	    
    {
      //sprintf(treestr, ":%8.20f", getBranchLength(tr, SUMMARIZE_LH, p));	      	   
      assert(tr->fracchange != -1.0);
      z = q->z[0];
      if (z < zmin) 
        z = zmin;      	 
      sprintf(treestr, ":%8.20f", -log(z) * tr->fracchange);	      	   
    }
    else	    
      sprintf(treestr, "%s", "\0");	    
  }

  while (*treestr) treestr++;
  return  treestr;
}

static double exp_pdf(double lambda, double x)
{
  return (lambda * exp(-(lambda * x))); 
}

static void printSimpleTree(tree *tr, boolean printBranchLengths, analdef *adef)
{
  Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, 1, 0, 0, 0, SUMMARIZE_LH, 0,0);
  fprintf(stderr, "%s\n", tr->tree_string);
}

static void printRecomTree(tree *tr, boolean printBranchLengths, char *title)
{
  FILE *nwfile;
  nwfile = myfopen("tmp.nw", "w+");
  Tree2StringRecomREC(tr->tree_string, tr, tr->start->back, printBranchLengths);
  fprintf(nwfile,"%s\n", tr->tree_string);
  fclose(nwfile);
  if(title)
    printBothOpen("%s\n", title);
  if (printBranchLengths)
    printBothOpen("%s\n", tr->tree_string);
  printBothOpen("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  system("bin/nw_display tmp.nw");
}  

static void set_start_prior(state *s)
{
  //need to go through branch lengths and model to set the priors
  //right now just 1
  s->curprior = 1;
}

static void nodeVisitor(nodeptr p, int numtips, int radius, int maxradius, nodeptr *list, int *list_size)
{
  //printBothOpen("visited %db%dr%d\n", p->number, p->back->number, radius);
  if(!isTip(p->back->number,numtips))
  {
    //printBothOpen("Add %db%d\n", p->number, p->back->number);
    list[*list_size] = p;
    *list_size += 1;
  }
  if(isTip(p->number,numtips) || radius > maxradius)
  {
    //printBothOpen(" -r%d\n", radius);
    return;
  }
  else
  {
    nodeVisitor(p->next->back, numtips, radius + 1, maxradius, list, list_size);
    nodeVisitor(p->next->next->back, numtips, radius + 1, maxradius, list, list_size);
  }
  return;
}


static void naiveInsertionProposal(state *s)
{
    int list_size = 0;
    if(!isTip(s->nb->number, s->tr->mxtips))
    {
      nodeVisitor(s->nb->next->back,       s->tr->mxtips, 1, s->maxradius, s->list, &list_size);
      nodeVisitor(s->nb->next->next->back, s->tr->mxtips, 1, s->maxradius, s->list, &list_size);
    }
    if(!isTip(s->nnb->number, s->tr->mxtips))
    {
      nodeVisitor(s->nnb->next->back,      s->tr->mxtips, 1, s->maxradius, s->list, &list_size);
      nodeVisitor(s->nnb->next->next->back,s->tr->mxtips, 1, s->maxradius, s->list, &list_size);
    }
    assert(list_size > 0);
    s->q = s->list[rand() % list_size];
    //printBothOpen(" %d candidates:",list_size);
    //printBothOpen(" insert at %db%d :",s->q->number, s->q->back->number);
    assert(s->q != NULL);
}

static void recordBranchInfo(nodeptr p, double *bl, int numBranches)
{
  int i;
  for(i = 0; i < numBranches; i++)
    bl[i] = p->z[i];
}


static nodeptr selectRandomSubtree(tree *tr)
{
  nodeptr 
    p;

  do
    {
      int 
	exitDirection = rand() % 3; 
     
      p = tr->nodep[(rand() % (tr->mxtips - 2)) + 1 + tr->mxtips];
      
      switch(exitDirection)
	{
	case 0:
	  break;
	case 1:
	  p = p->next;
	  break;
	case 2:
	  p = p->next->next;
	  break;
	default:
	  assert(0);
	}
    }
  while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));

  assert(!isTip(p->number, tr->mxtips));

  return p;
}

static void doSPR(tree *tr, state *instate)
{
  nodeptr    
    p = selectRandomSubtree(tr);
  
  /* evaluateGeneric(tr, tr->start, TRUE);
     printf("%f \n", tr->likelihood);*/

  parsimonySPR(p, tr);
  
  /*evaluateGeneric(tr, tr->start, TRUE);
    printf("%f \n", tr->likelihood);*/

  instate->p = p;
  instate->nb  = p->next->back;
  instate->nnb = p->next->next->back;
  
  recordBranchInfo(instate->nb, instate->nbz, instate->tr->numBranches);
  recordBranchInfo(instate->nnb, instate->nnbz, instate->tr->numBranches);

  removeNodeBIG(tr, p,  tr->numBranches);
  instate->q = tr->insertNode;
  instate->r = instate->q->back;
  recordBranchInfo(instate->q, instate->qz, instate->tr->numBranches);

  assert(Thorough == 0);
  
  insertBIG(instate->tr, instate->p, instate->q, instate->tr->numBranches);
  evaluateGeneric(instate->tr, instate->p->next->next, FALSE); 
  /*testInsertBIG(tr, p, tr->insertNode);*/

  printf("%f \n", tr->likelihood);
}



static nodeptr selectRandomInnerSubtree(tree *tr)
{
  nodeptr p;
  int pruned_id;
  do
  {
    pruned_id = (rand() % (tr->mxtips - 2)) + 1 + tr->mxtips ;
    assert(pruned_id > tr->mxtips);
    assert(pruned_id <= 2*tr->mxtips - 2);
    p = tr->nodep[pruned_id];
  }
  while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));
  //printBothOpen("Selected %db%d to prune\n", p->number, p->back->number);
  return p;
}

static boolean simpleNodeProposal(state * instate)
{
  //prior is flat for these moves
  instate->newprior = 1;
  instate->p = selectRandomInnerSubtree(instate->tr);
  /* records info pre-pruning */
  instate->nb = instate->p->next->back;
  instate->nnb = instate->p->next->next->back;
  //printBothOpen("selected prune node %db%d bl %f \n", instate->p->number, instate->p->back->number, instate->p->z[0]);
  recordBranchInfo(instate->nb, instate->nbz, instate->tr->numBranches);
  recordBranchInfo(instate->nnb, instate->nnbz, instate->tr->numBranches);
  /* prune subtree p */
  if (removeNodeBIG(instate->tr, instate->p,  instate->tr->numBranches) == NULL) assert(FALSE);
  /* insert somewhere else, but it must not be in the pruned subtree */
  //printBothOpen("pruned %db%d \n", instate->p->number, instate->p->back->number);
  instate->q = (nodeptr) NULL;
  naiveInsertionProposal(instate);
  
  if(instate->q!=NULL)
    {
      instate->r = instate->q->back;
      recordBranchInfo(instate->q, instate->qz, instate->tr->numBranches);
      /*
	printBothOpen("inserted %db%d at %db%d where bl %f, Thorough is %d \n", 
	instate->p->number, instate->p->back->number,
	instate->q->number, instate->q->back->number, 
	instate->q->z[0], Thorough);
      */
      if (! insertBIG(instate->tr, instate->p, instate->q, instate->tr->numBranches)) 
	assert(FALSE);
      //TODO: breaks here evaluateGenericSpecial.c:1164: evaluateIterative: Assertion `partitionLikelihood < 0.0' failed.
      evaluateGeneric(instate->tr, instate->p->next->next, FALSE);    
      return TRUE;
    }
  else
    return FALSE;
}

static void resetSimpleNodeProposal(state * instate)
{
  /* prune the insertion */
  hookup(instate->q, instate->r, instate->qz, instate->tr->numBranches);
  instate->p->next->next->back = instate->p->next->back = (nodeptr) NULL;
  /* insert the pruned tree in its original node */
  hookup(instate->p->next,       instate->nb, instate->nbz, instate->tr->numBranches);
  hookup(instate->p->next->next, instate->nnb, instate->nnbz, instate->tr->numBranches);
  newviewGeneric(instate->tr, instate->p, FALSE); 
}


static void reset_branch_length(nodeptr p, int numBranches)
{
  int i;
  double new_value;
  for(i = 0; i < numBranches; i++)
  {
    assert(p->z_tmp[i] == p->back->z_tmp[i]);
    p->z[i] = p->back->z[i] = p->z_tmp[i];   /* restore saved value */
  }
}


/* start NNI topologycal prop */
static void recordNNIBranchInfo(nodeptr p, int numBranches)
{
  int i;
  nodeptr q = p->back;
  for(i = 0; i < numBranches; i++)
  {
    p->z_tmp[i] =  p->back->z_tmp[i] = p->z[i];
    p->next->z_tmp[i] =  p->next->back->z_tmp[i] = p->next->z[i];
    p->next->next->z_tmp[i] = p->next->next->back->z_tmp[i] = p->next->next->z[i];
    q->next->z_tmp[i] = q->next->back->z_tmp[i] = q->next->z[i];
    q->next->next->z_tmp[i] = q->next->next->back->z_tmp[i] = q->next->next->z[i];
  }
}
static showNNI_move(nodeptr p)
{
  printBothOpen("NNI from p %d %.6f, pnb %d %.6f, pnnb %d %.6f\n", 
      p->number, p->z[0],
      p->next->back->number, p->next->back->z[0], 
      p->next->next->back->number,  p->next->next->back->z[0]);
  nodeptr q = p->back;
  printBothOpen("NNI from q %d %.6f, qnb %d %.6f, qnnb %d %.6f\n", 
      q->number, q->z[0],
      q->next->back->number, q->next->back->z[0], 
      q->next->next->back->number,  q->next->next->back->z[0]);
}
//setting this out to allow for other types of setting
static void set_branch_length_sliding_window(nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
{
  int i;
  double new_value;
  double r,mx,mn;
  for(i = 0; i < numBranches; i++)
  {
    assert(p->z[i] == p->back->z[i]);
    if(record_tmp_bl)
      p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /* keep current value */
    r = (double)rand()/(double)RAND_MAX;
    mn = p->z[i]-(s->bl_sliding_window_w/2);
    mx = p->z[i]+(s->bl_sliding_window_w/2);
    new_value = fabs(mn + r * (mx-mn));
    /* Ensure always you stay within this range */
    if(new_value > zmax) new_value = zmax;
    if(new_value < zmin) new_value = zmin;
    assert(new_value <= zmax && new_value >= zmin);

    p->z[i] = p->back->z[i] = new_value;
    //assuming this will be visiting each node, and multiple threads won't be accessing this
     s->bl_prior += log(exp_pdf(s->bl_prior_exp_lambda,new_value));
    //s->bl_prior += 1;
  }
}
static void hookupBL(nodeptr p, nodeptr q, nodeptr bl_p, state *s)
{
   set_branch_length_sliding_window(bl_p, s->tr->numBranches, s, FALSE);
   hookup(p, q, bl_p->z, s->tr->numBranches);
}
static boolean stNNIproposal(state *s)
{
  //s->newprior = 1;
  s->bl_prior = 0;
  int attempts = 0;
  do{
    s->p = selectRandomInnerSubtree(s->tr); /* TODOFER do this ad hoc for NNI requirements*/
    if (++attempts > 500)
      return FALSE;
  }while(isTip(s->p->number, s->tr->mxtips) || isTip(s->p->back->number, s->tr->mxtips));
  assert(!isTip(s->p->number, s->tr->mxtips));
  nodeptr 
    p = s->p,
    q = s->p->back,
    pb1 = s->p->next->back,
    pb2 = s->p->next->next->back;
  assert(!isTip(q->number, s->tr->mxtips));
  nodeptr
    qb1 = q->next->back,
    qb2 = q->next->next->back;

  recordNNIBranchInfo(p, s->tr->numBranches);
  /* do only one type of NNI, nni1 */
  double randprop = (double)rand()/(double)RAND_MAX;
  boolean changeBL = TRUE;
  if (randprop < 1.0 / 3.0)
  {
    s->whichNNI = 1;
    if(!changeBL)
    {
      hookup(p, q, p->z, s->tr->numBranches);
      hookup(p->next,       qb1, q->next->z, s->tr->numBranches);
      hookup(p->next->next, pb2, p->next->next->z, s->tr->numBranches);
      hookup(q->next,       pb1, p->next->z, s->tr->numBranches);
      hookup(q->next->next, qb2, q->next->next->z, s->tr->numBranches);
    }
    else
    {
      hookupBL(p, q, p, s);
      hookupBL(p->next,       qb1, q->next, s);
      hookupBL(p->next->next, pb2, p->next->next, s);
      hookupBL(q->next,       pb1, p->next, s);
      hookupBL(q->next->next, qb2, q->next->next, s);
    }
  }
  else if (randprop < 2.0 / 3.0)
  {
    s->whichNNI = 2;
    if(!changeBL)
    {
      hookup(p, q, p->z, s->tr->numBranches);
      hookup(p->next,       pb1, p->next->z, s->tr->numBranches);
      hookup(p->next->next, qb1, q->next->z, s->tr->numBranches);
      hookup(q->next,       pb2, p->next->next->z, s->tr->numBranches);
      hookup(q->next->next, qb2, q->next->next->z, s->tr->numBranches);
    }
    else
    {
      hookupBL(p, q, p, s);
      hookupBL(p->next,       pb1, p->next, s);
      hookupBL(p->next->next, qb1, q->next, s);
      hookupBL(q->next,       pb2, p->next->next, s);
      hookupBL(q->next->next, qb2, q->next->next, s);
    }
  }
  else
  {
    /* change only the branch lengths */
    s->whichNNI = 0; 
    if(changeBL)
    {
      /* do it like this for symmetry */
      hookupBL(p, q, p, s);
      hookupBL(p->next,       pb1, p->next, s);
      hookupBL(p->next->next, pb2, p->next->next, s);
      hookupBL(q->next,       qb1, q->next, s);
      hookupBL(q->next->next, qb2, q->next->next, s);
    }
  }

  newviewGeneric(s->tr, p, FALSE);
  newviewGeneric(s->tr, p->back, FALSE);
  evaluateGeneric(s->tr, p, FALSE);
  return TRUE;
}
reset_stNNI(state *s)
{
  nodeptr p, q;
  p = s->p;
  q = p->back;
  nodeptr pb1, pb2;
  pb1 = p->next->back;
  pb2 = p->next->next->back;
  nodeptr qb1, qb2;
  qb1 = q->next->back;
  qb2 = q->next->next->back;

  /* whichNNI 1*/
  if(s->whichNNI == 1)
  {
    hookup(p, q, q->z_tmp, s->tr->numBranches);
    hookup(p->next,       qb1, p->next->z_tmp, s->tr->numBranches);
    hookup(p->next->next, pb2, p->next->next->z_tmp, s->tr->numBranches);
    hookup(q->next,       pb1, q->next->z_tmp, s->tr->numBranches);
    hookup(q->next->next, qb2, q->next->next->z_tmp, s->tr->numBranches);
  }
  else if(s->whichNNI == 2)
  {
    hookup(p, q, p->z_tmp, s->tr->numBranches);
    hookup(p->next,       pb1, p->next->z_tmp, s->tr->numBranches);
    hookup(p->next->next, qb1, p->next->next->z_tmp, s->tr->numBranches);
    hookup(q->next,       pb2, q->next->z_tmp, s->tr->numBranches);
    hookup(q->next->next, qb2, q->next->next->z_tmp, s->tr->numBranches);
  }
  else
  {
    assert(s->whichNNI == 0);
    hookup(p, q, p->z_tmp, s->tr->numBranches);
    hookup(p->next,       pb1, p->next->z_tmp, s->tr->numBranches);
    hookup(p->next->next, pb2, p->next->next->z_tmp, s->tr->numBranches);
    hookup(q->next,       qb1, q->next->z_tmp, s->tr->numBranches);
    hookup(q->next->next, qb2, q->next->next->z_tmp, s->tr->numBranches);
  }

  /*
  reset_branch_length(p, s->tr->numBranches);
  reset_branch_length(p->next, s->tr->numBranches);
  reset_branch_length(p->next->next, s->tr->numBranches);
  reset_branch_length(q->next, s->tr->numBranches);
  reset_branch_length(q->next->next, s->tr->numBranches);
  */

  newviewGeneric(s->tr, p, FALSE);
  newviewGeneric(s->tr, q, FALSE);
}

/* end NNI topologycal prop */


static void traverse_branches(nodeptr p, int *count, state * s, boolean resetBL)
{
  nodeptr q;
  //printf("current BL at %db%d: %f\n", p->number, p->back->number, p->z[0]);
  if(resetBL)
    reset_branch_length(p, s->tr->numBranches);
  else//can allow for other methods later
    set_branch_length_sliding_window(p, s->tr->numBranches, s, TRUE);
  *count += 1;


  if (! isTip(p->number, s->tr->mxtips)) 
  {                                  /*  Adjust descendants */
    q = p->next;
    while (q != p) 
    {
      traverse_branches(q->back, count, s, resetBL);
      q = q->next;
    }	
    newviewGeneric(s->tr, p, FALSE);     // not sure if we need this
  }
}

static void update_all_branches(state * s, boolean resetBL)
{
  int updated_branches = 0;
  assert(isTip(s->tr->start->number, s->tr->mxtips));
  /* visit each branch exactly once */
  traverse_branches(s->tr->start->back, &updated_branches, s, resetBL);
  assert(updated_branches == s->tr->mxtips + s->tr->mxtips - 3);
}

static void set_start_bl(state * instate)
{
  update_all_branches(instate, FALSE);
}


/*
 * should be sliding window proposal
 */

static boolean simpleBranchLengthProposal(state * instate)
{
   
  //for each branch get the current branch length
  //pull a uniform like
  //x = current, w =window
  //uniform(x-w/2,x+w/2)

  update_all_branches(instate, FALSE);
  evaluateGeneric(instate->tr, instate->tr->start, FALSE); /* update the tr->likelihood */

  //for prior, just using exponential for now
  //calculate for each branch length
  // where lambda is chosen and x is the branch length
  //lambda * exp(-lamba * x)

  //only calculate the new ones
  //
  return TRUE;
}

static void resetSimpleBranchLengthProposal(state * instate)
{
  update_all_branches(instate, TRUE);
}


static void printSubsRates(tree *tr,int model, int numSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  printBothOpen("Subs rates: ");
  for(i=0; i<numSubsRates; i++)
    printBothOpen("%d => %.3f, ", i, tr->partitionData[model].substRates[i]);
  printBothOpen("\n\n");
}
static void recordSubsRates(tree *tr, int model, int numSubsRates, double *prevSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    prevSubsRates[i] = tr->partitionData[model].substRates[i];
}
static void restoreSubsRates(tree *tr, analdef *adef, int model, int numSubsRates, double *prevSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    tr->partitionData[model].substRates[i] = prevSubsRates[i];
#ifndef _LOCAL_DISCRETIZATION
  initReversibleGTR(tr, model);
#endif
  /* TODO need to broadcast rates here for parallel version */

  evaluateGeneric(tr, tr->start, TRUE);
}
static void editSubsRates(tree *tr, int model, int subRatePos, double subRateValue)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  assert(subRateValue <= RATE_MAX && subRateValue >= RATE_MIN);
  int states = tr->partitionData[model].states; 
  int numSubsRates = (states * states - states) / 2;
  assert(subRatePos >= 0 && subRatePos < numSubsRates);
  tr->partitionData[model].substRates[subRatePos] = subRateValue;
}

/*
 * should be sliding window proposal
 */
static void simpleModelProposal(state * instate)
{
  //TODO: add safety to max and min values
  //record the old ones
  recordSubsRates(instate->tr, instate->model, instate->numSubsRates, instate->curSubsRates);
  //choose a random set of model params,
  //probably with dirichlet proposal
  //with uniform probabilities, no need to have other
  int state;
  double new_value,curv;
  double r,mx,mn;
  //using the branch length sliding window for a test
  for(state = 0;state<instate->numSubsRates ; state ++)
    {
      curv = instate->tr->partitionData[instate->model].substRates[state];
      r = (double)rand()/(double)RAND_MAX;
      mn = curv-(instate->rt_sliding_window_w/2);
      mx = curv+(instate->rt_sliding_window_w/2);
      new_value = fabs(mn + r * (mx-mn));
      /* Ensure always you stay within this range */
      if(new_value > RATE_MAX) new_value = RATE_MAX;
      if(new_value < RATE_MIN) new_value = RATE_MIN;
      //printf("%i %f %f\n", state, curv, new_value);
      editSubsRates(instate->tr,instate->model, state, new_value);
    }
  //recalculate eigens
#ifndef _LOCAL_DISCRETIZATION
  initReversibleGTR(instate->tr, instate->model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */
#endif

  /* TODO: need to broadcast rates here for parallel version ! */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */
  //TODO: without this, the run will fail after a successful model, but failing SPR
  evaluateGeneric(instate->tr, instate->tr->start, FALSE);
  //for prior, just use dirichlet
  // independent gamma distribution for each parameter
  //the pdf for this is
  // for gamma the prior is gamma

  //for statefreqs should all be uniform

  //only calculate the new ones
}

static void resetSimpleModelProposal(state * instate)
{
  restoreSubsRates(instate->tr, instate->adef, instate->model, instate->numSubsRates, instate->curSubsRates);
  evaluateGeneric(instate->tr, instate->tr->start, FALSE);
}

//simple sliding window
static void simpleGammaProposal(state * instate)
{
  //TODO: add safety to max and min values
  double newalpha, curv, r,mx,mn;
  curv = instate->tr->partitionData[instate->model].alpha;
  instate->curAlpha = curv;
  r = (double)rand()/(double)RAND_MAX;
  mn = curv-(instate->gm_sliding_window_w/2);
  mx = curv+(instate->gm_sliding_window_w/2);
  newalpha = fabs(mn + r * (mx-mn));
  /* Ensure always you stay within this range */
  if(newalpha > ALPHA_MAX) newalpha = ALPHA_MAX;
  if(newalpha < ALPHA_MIN) newalpha = ALPHA_MIN;
  instate->tr->partitionData[instate->model].alpha = newalpha;

#ifndef _LOCAL_DISCRETIZATION
  makeGammaCats(instate->tr->partitionData[instate->model].alpha, instate->tr->partitionData[instate->model].gammaRates, 4);
#endif

  /* TODO: for the parallel version: need to broadcast the gamma rates before re-evaluating !!!! 
     also note the _LOCAL_DISCRETIZATION flag that should only be used for the parallel stuff !
   */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE);
}

static void resetSimpleGammaProposal(state * instate)
{
  instate->tr->partitionData[instate->model].alpha = instate->curAlpha;
#ifndef _LOCAL_DISCRETIZATION
  makeGammaCats(instate->tr->partitionData[instate->model].alpha, instate->tr->partitionData[instate->model].gammaRates, 4);
#endif

   /* TODO: for the parallel version: need to broadcast the gamma rates before re-evaluating !!!! 
     also note the _LOCAL_DISCRETIZATION flag that should only be used for the parallel stuff !
   */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE);
}

static prop proposal(state * instate)
/* so here the idea would be to randomly choose among proposals? we can use typedef enum to label each, and return that */ 
{
  double randprop = (double)rand()/(double)RAND_MAX;
  boolean proposalSuccess;
  //double start_LH = evaluateGeneric(instate->tr, instate->tr->start); /* for validation */
  prop proposal_type;
  //simple proposal
  if(randprop < 0.25) 
  {
    if(randprop < 0.2)//TOPOLOGICAL MOVE
      {
	if(randprop > 0.1)//SPR MOVE
	  {
	    proposal_type = SPR;
	    // printBothOpen("Propose SPR\n");
	    if (randprop < 0.15)
	      instate->maxradius = 1;
	    else
	      instate->maxradius = 2;
	    
	    doSPR(instate->tr, instate);

	    proposalSuccess = TRUE;

	    /* TODO */
	    /*proposalSuccess = simpleNodeProposal(instate);*/
	  }
	else
	  {
	    proposal_type = stNNI;
	    proposalSuccess = stNNIproposal(instate); 
	    if(proposalSuccess == FALSE)
	      {
		/* TODOFER this came up with ds 20 and GTRPSR, see why */
		printBothOpen("WARNING!! stNNI proposal failed, doing SPR\n");
		proposal_type = SPR;
		instate->maxradius = 1;
		proposalSuccess = simpleNodeProposal(instate);
	      }
	  }
	if(proposalSuccess == FALSE)
	  {
	    assert(FALSE); // this should either never happen or look below and return PROPOSAL_FAILED to react accordingly
	  }
	else
	  {
	    /* A moved has been made, previous state is in instate */
	    if(proposal_type != stNNI) /*TODOFER delete this when bl are changed*/
	      assert(instate->tr->startLH != instate->tr->likelihood);
	  }
      }
    else{//MODEL
      proposal_type = UPDATE_MODEL;
      simpleModelProposal(instate);
    }
  }
  else
    {
      if(randprop < 0.95)//UPDATE_ALL_BL
	{
	  proposal_type = UPDATE_ALL_BL;
	  instate->bl_prior = 0;
	  //printBothOpen("Propose BL_UPDATE\n");
	  assert(proposal_type == UPDATE_ALL_BL);
	  proposalSuccess = simpleBranchLengthProposal(instate);
	  assert(instate->tr->startLH != instate->tr->likelihood);
	  assert(proposalSuccess);
	}
      else//GAMMA
	{
	  proposal_type = UPDATE_GAMMA;
	  simpleGammaProposal(instate);
	}
    }
  //record the curprior
  instate->newprior = instate->bl_prior;
  return proposal_type;
}

static void resetState(prop proposal_type, state * curstate)
{
  switch(proposal_type)
  {
    case SPR:
      resetSimpleNodeProposal(curstate);
      break;
    case stNNI:
      reset_stNNI(curstate);
      break;
    case UPDATE_ALL_BL:
      resetSimpleBranchLengthProposal(curstate);
      //printf("RESETBL\n");
      break;
    case UPDATE_MODEL:
      resetSimpleModelProposal(curstate);
      break;
    case UPDATE_GAMMA:
      resetSimpleGammaProposal(curstate);
      break;
    default:
      assert(FALSE);
  }
}

static void printStateFile(int iter, state * curstate)
{ 
  FILE *f = myfopen("sampled_states.txt", "ab");
  fprintf(f,"%d\t%f",iter, curstate->tr->likelihood);
  int i;
  for(i = 0;i < curstate->numSubsRates; i++)
  {
    fprintf(f,"\t%f",curstate->curSubsRates[i]);
  }
  fprintf(f,"\t%f",curstate->curAlpha);
  fprintf(f,"\n");
  fclose(f);
}

static void printStateFileHeader(state * curstate)
{ 
  //DELETE THE FILE IF IT EXISTS, obviously this should be better later
  remove("sampled_states.txt");
  FILE *f = myfopen("sampled_states.txt", "ab");
  fprintf(f,"iter\tlikelihood");
  int i;
  for(i = 0;i < curstate->numSubsRates; i++)
  {
    fprintf(f,"\trate%d",i);
  }
  fprintf(f,"\talpha");
  fprintf(f,"\n");
  fclose(f);
}

void mcmc(tree *tr, analdef *adef)
{
  int i=0;

  tr->startLH = tr->likelihood;
  printBothOpen("start minimalistic search with LH %f\n", tr->likelihood);
  printBothOpen("tr LH %f, startLH %f\n", tr->likelihood, tr->startLH);
  
  int insert_id;
  int j;

  int maxradius = 30;
  int accepted_spr = 0, accepted_nni = 0, accepted_bl = 0, accepted_model = 0, accepted_gamma = 0, inserts = 0;
  int rejected_spr = 0, rejected_nni = 0, rejected_bl = 0, rejected_model = 0, rejected_gamma = 0;
  int num_moves = 10000;
  boolean proposalAccepted;
  boolean proposalSuccess;
  prop which_proposal;
  double testr;
  double acceptance;

  srand (440);
  double totalTime = 0.0, proposalTime = 0.0, blTime = 0.0, printTime = 0.0;
  double t_start = gettime();
  double t;


  //allocate states
  double bl_prior_exp_lambda = 0.1;
  double bl_sliding_window_w = 0.005;
  double gm_sliding_window_w = 0.75;
  double rt_sliding_window_w = 0.5;
  state *curstate = state_init(tr, adef, maxradius, bl_sliding_window_w, rt_sliding_window_w, gm_sliding_window_w, bl_prior_exp_lambda);
  printStateFileHeader(curstate);
  set_start_bl(curstate);
  printf("start bl_prior: %f\n",curstate->bl_prior);
  set_start_prior(curstate);
  curstate->hastings = 1;//needs to be set by the proposal when necessary

  /* Set the starting LH with a full traversal */
  evaluateGeneric(tr, tr->start, TRUE);	 
  tr->startLH = tr->likelihood;
  printBothOpen("Starting with tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);

  /* Set reasonable model parameters */
  evaluateGeneric(curstate->tr, curstate->tr->start, FALSE); // just for validation 
  printBothOpen("tr LH before modOpt %f\n",curstate->tr->likelihood);
  printSubsRates(curstate->tr, curstate->model, curstate->numSubsRates);

  /* optimize the model with Brents method for reasonable starting points */
  modOpt(curstate->tr, curstate->adef, 5.0); /* not by proposal, just using std raxml machinery... */
  evaluateGeneric(curstate->tr, curstate->tr->start, FALSE); // just for validation 
  printBothOpen("tr LH after modOpt %f\n",curstate->tr->likelihood);
  printSubsRates(curstate->tr, curstate->model, curstate->numSubsRates);
  recordSubsRates(curstate->tr, curstate->model, curstate->numSubsRates, curstate->curSubsRates);

  int first = 1;
  /* beginning of the MCMC chain */
  for(j=0; j<num_moves; j++)
  {
    //printBothOpen("iter %d, tr LH %f, startLH %f\n",j, tr->likelihood, tr->startLH);
    //printRecomTree(tr, TRUE, "startiter");
    proposalAccepted = FALSE;
    t = gettime(); 

    /*
      evaluateGeneric(tr, tr->start); // just for validation 
      printBothOpen("before proposal, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
    */

    which_proposal = proposal(curstate);
    if (first == 1)
    {
      first = 0;
      curstate->curprior = curstate->newprior;
    }
   //printBothOpen("proposal done, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
    assert(which_proposal == SPR || which_proposal == stNNI ||
           which_proposal == UPDATE_ALL_BL || 
           which_proposal == UPDATE_MODEL || which_proposal == UPDATE_GAMMA);
    proposalTime += gettime() - t;
    /* decide upon acceptance */
    testr = (double)rand()/(double)RAND_MAX;
    //should look something like 
    acceptance = fmin(1,(curstate->hastings) * 
		      (exp(curstate->newprior-curstate->curprior)) * (exp(curstate->tr->likelihood-curstate->tr->startLH)));
    
    /*
      //printRecomTree(tr, FALSE, "after proposal");
      printBothOpen("after proposal, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
    */
    if(testr < acceptance)
    {
      proposalAccepted = TRUE;

      switch(which_proposal)
	{
	case SPR:      
	  //printRecomTree(tr, TRUE, "after accepted");
	  // printBothOpen("SPR new topology , iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
	   accepted_spr++;
	  break;
	case stNNI:	  
	  printBothOpen("NNI new topology , iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
	  accepted_nni++;
	  break;
	case UPDATE_ALL_BL:	  
	  //      printBothOpen("BL new , iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
	  accepted_bl++;
	  break;
	case UPDATE_MODEL:      
	  //	printBothOpen("Model new, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
	  accepted_model++;
	  break;
	case UPDATE_GAMMA:      
	  //	printBothOpen("Gamma new, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
	  accepted_gamma++;
	  break;
	default:
	  assert(0);
	}

      curstate->tr->startLH = curstate->tr->likelihood;  //new LH
      curstate->curprior = curstate->newprior;          
    }
    else
    {
      //printBothOpen("rejected , iter %d tr LH %f, startLH %f, %i \n", j, tr->likelihood, tr->startLH, which_proposal);
      resetState(which_proposal,curstate);
      
      switch(which_proposal)
	{
	case SPR:
	  rejected_spr++;
	  break;
	case stNNI:
	  rejected_nni++;
	  break;
	case UPDATE_ALL_BL:
	  rejected_bl++;
	  break;
	case UPDATE_MODEL:
	  rejected_model++;
	  break;
	case UPDATE_GAMMA:
	  rejected_gamma++;
	  break;
	default:
	  assert(0);
	}
      
      evaluateGeneric(tr, tr->start, FALSE); 
      
      // just for validation 

      if(fabs(curstate->tr->startLH - tr->likelihood) > 1.0E-10)
      {
        printBothOpen("WARNING: LH diff %.10f\n", curstate->tr->startLH - tr->likelihood);
      }
      //printRecomTree(tr, TRUE, "after reset");
      //printBothOpen("after reset, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
      assert(fabs(curstate->tr->startLH - tr->likelihood) < 1.0E-10);
    }       
    inserts++;
    
    /* need to print status */
    if (j % 50 == 0)
    {
      t = gettime(); 
      printBothOpen("sampled at iter %d, tr LH %f, startLH %f, prior %f, incr %f\n",j, tr->likelihood, tr->startLH, curstate->curprior, tr->likelihood - tr->startLH);
      boolean printBranchLengths = TRUE;
      /*printSimpleTree(tr, printBranchLengths, adef);*/
      //TODO: print some parameters to a file 
      printStateFile(j,curstate);
      printTime += gettime() - t;
    }
  }

  t = gettime(); 
  treeEvaluate(tr, 1);
  blTime += gettime() - t;
  printBothOpen("accepted SPR %d, accepted stNNI %d, accepted BL %d, accepted model %d, accepted gamma %d, num moves tried %d, SPRs with max radius %d\n", 
		accepted_spr, accepted_nni, accepted_bl, accepted_model, accepted_gamma, num_moves, maxradius);
  printBothOpen("rejected SPR %d, rejected stNNI %d, rejected BL %d, rejected model %d, rejected gamma %d\n",
		rejected_spr, rejected_nni, rejected_bl, rejected_model, rejected_gamma);
  printBothOpen("ratio SPR %f, ratio stNNI %f,  ratio BL %f, ratio model %f, ratio gamma %f\n",
		accepted_spr/(double)(rejected_spr+accepted_spr), accepted_nni/(double)(rejected_nni+accepted_nni), accepted_bl/(double)(rejected_bl+accepted_bl), 
		accepted_model/(double)(rejected_model+accepted_model), accepted_gamma/(double)(rejected_gamma+accepted_gamma));
  printBothOpen("total  %f, BL %f, printing %f, proposal %f\n", gettime()- t_start, blTime, printTime, proposalTime);
  assert(inserts == num_moves);
  state_free(curstate);
}

