/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

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

extern double accumulatedTime;

extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern char run_id[128];
extern double masterTime;
extern partitionLengths pLengths[MAX_MODEL];
extern char binaryCheckpointName[1024];
extern char binaryCheckpointInputName[1024];

boolean initrav (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  if (!isTip(p->number, tr->mxtips)) 
  {      
    q = p->next;

    do 
    {	   
      if (! initrav(tr, q->back))  return FALSE;		   
      q = q->next;	
    } 
    while (q != p);  

    newviewGeneric(tr, p, FALSE);
  }

  return TRUE;
} 












boolean update(tree *tr, nodeptr p)
{       
  nodeptr  q; 
  boolean smoothedPartitions[NUM_BRANCHES];
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  double _deltaz;

  q = p->back;   

  for(i = 0; i < tr->numBranches; i++)
    z0[i] = q->z[i];    

  if(tr->numBranches > 1)
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, TRUE);  
  else
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, FALSE);

  for(i = 0; i < tr->numBranches; i++)    
    smoothedPartitions[i]  = tr->partitionSmoothed[i];

  for(i = 0; i < tr->numBranches; i++)
  {         
    if(!tr->partitionConverged[i])
    {	  
      _deltaz = deltaz;

      if(ABS(z[i] - z0[i]) > _deltaz)  
      {	      
        smoothedPartitions[i] = FALSE;       
      }	 


      p->z[i] = q->z[i] = z[i];	 
    }
  }

  for(i = 0; i < tr->numBranches; i++)    
    tr->partitionSmoothed[i]  = smoothedPartitions[i];

  return TRUE;
}

boolean smooth (tree *tr, nodeptr p)
{
  nodeptr  q;

  if (! update(tr, p))               return FALSE; /*  Adjust branch */

  if (! isTip(p->number, tr->mxtips)) 
  {                                  /*  Adjust descendants */
    q = p->next;
    while (q != p) 
    {
      if (! smooth(tr, q->back))   return FALSE;
      q = q->next;
    }	

    if(tr->numBranches > 1 && !tr->useRecom)		  
      newviewGeneric(tr, p, TRUE);	
    else
      newviewGeneric(tr, p, FALSE);     
  }

  return TRUE;
} 

static boolean allSmoothed(tree *tr)
{
  int i;
  boolean result = TRUE;

  for(i = 0; i < tr->numBranches; i++)
  {
    if(tr->partitionSmoothed[i] == FALSE)
      result = FALSE;
    else
      tr->partitionConverged[i] = TRUE;
  }

  return result;
}


/* do maxtimes rounds of branch length optimization */
boolean smoothTree (tree *tr, int maxtimes)
{
	nodeptr  p, q;
	int i, count = 0;

	p = tr->start;
	for(i = 0; i < tr->numBranches; i++)
		tr->partitionConverged[i] = FALSE;

	while (--maxtimes >= 0)
	{
		for(i = 0; i < tr->numBranches; i++)
			tr->partitionSmoothed[i] = TRUE;

		if (! smooth(tr, p->back))       return FALSE;
		if (!isTip(p->number, tr->mxtips))
		{
			q = p->next;
			while (q != p)
			{
				if (! smooth(tr, q->back))   return FALSE;
				q = q->next;
			}
		}

		count++;

		if (allSmoothed(tr))
			break;
	}

	for(i = 0; i < tr->numBranches; i++)
		tr->partitionConverged[i] = FALSE;

	return TRUE;
} 



boolean localSmooth (tree *tr, nodeptr p, int maxtimes)
{ 
  nodeptr  q;
  int i;

  if (isTip(p->number, tr->mxtips)) return FALSE;

  for(i = 0; i < tr->numBranches; i++)	
    tr->partitionConverged[i] = FALSE;	

  while (--maxtimes >= 0) 
  {     
    for(i = 0; i < tr->numBranches; i++)	
      tr->partitionSmoothed[i] = TRUE;

    q = p;
    do 
    {
      if (! update(tr, q)) return FALSE;
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr)) 
      break;
  }

  for(i = 0; i < tr->numBranches; i++)
  {
    tr->partitionSmoothed[i] = FALSE; 
    tr->partitionConverged[i] = FALSE;
  }

  return TRUE;
}





static void resetInfoList(infoList *iList)
{
  int 
    i;

  iList->valid = 0;

  for(i = 0; i < iList->n; i++)    
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = unlikely;
  }    
}

static void initInfoList(infoList *iList, size_t n)
{
  int 
    i;

  iList->n = n;
  iList->valid = 0;
  iList->list = (bestInfo *)malloc(sizeof(bestInfo) * (size_t)n);

  for(i = 0; i < n; i++)
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = unlikely;
  }
}

static void freeInfoList(infoList *iList)
{ 
  free(iList->list);   
}


static void insertInfoList(nodeptr node, double likelihood, infoList *iList)
{
  int 
    i,
    min = 0;

  double 
    min_l =  iList->list[0].likelihood;

  for(i = 1; i < iList->n; i++)
  {
    if(iList->list[i].likelihood < min_l)
    {
      min = i;
      min_l = iList->list[i].likelihood;
    }
  }

  if(likelihood > min_l)
  {
    iList->list[min].likelihood = likelihood;
    iList->list[min].node = node;
    iList->valid += 1;
  }

  if(iList->valid > iList->n)
    iList->valid = iList->n;
}


boolean smoothRegion (tree *tr, nodeptr p, int region)
{ 
  nodeptr  q;

  if (! update(tr, p))               return FALSE; /*  Adjust branch */

  if(region > 0)
  {
    if (!isTip(p->number, tr->mxtips)) 
    {                                 
      q = p->next;
      while (q != p) 
      {
        if (! smoothRegion(tr, q->back, --region))   return FALSE;
        q = q->next;
      }	

      newviewGeneric(tr, p, FALSE);
    }
  }

  return TRUE;
}

boolean regionalSmooth (tree *tr, nodeptr p, int maxtimes, int region)
{
  nodeptr  q;
  int i;

  if (isTip(p->number, tr->mxtips)) return FALSE;            /* Should be an error */

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;

  while (--maxtimes >= 0) 
  {	
    for(i = 0; i < tr->numBranches; i++)	  
      tr->partitionSmoothed[i] = TRUE;

    q = p;
    do 
    {
      if (! smoothRegion(tr, q, region)) return FALSE;
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr)) 
      break;
  }

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;

  return TRUE;
} /* localSmooth */





nodeptr  removeNodeBIG (tree *tr, nodeptr p, int numBranches)
{  
  double   zqr[NUM_BRANCHES], result[NUM_BRANCHES];
  nodeptr  q, r;
  int i;

  q = p->next->back;
  r = p->next->next->back;

  for(i = 0; i < numBranches; i++)
    zqr[i] = q->z[i] * r->z[i];        

  makenewzGeneric(tr, q, r, zqr, iterations, result, FALSE);   

  for(i = 0; i < numBranches; i++)        
    tr->zqr[i] = result[i];

  hookup(q, r, result, numBranches); 

  p->next->next->back = p->next->back = (node *) NULL;

  return  q; 
}

nodeptr  removeNodeRestoreBIG (tree *tr, nodeptr p)
{
  nodeptr  q, r;

  q = p->next->back;
  r = p->next->next->back;  

  newviewGeneric(tr, q, FALSE);
  newviewGeneric(tr, r, FALSE);

  hookup(q, r, tr->currentZQR, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;

  return  q;
}


boolean insertBIG (tree *tr, nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r, s;
  int i;

  r = q->back;
  s = p->back;

  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];

  if(tr->thoroughInsertion)
  { 
    double  zqr[NUM_BRANCHES], zqs[NUM_BRANCHES], zrs[NUM_BRANCHES], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;      
    double defaultArray[NUM_BRANCHES];	
    double e1[NUM_BRANCHES], e2[NUM_BRANCHES], e3[NUM_BRANCHES];
    double *qz;

    qz = q->z;

    for(i = 0; i < numBranches; i++)
      defaultArray[i] = defaultz;

    makenewzGeneric(tr, q, r, qz, iterations, zqr, FALSE);           
    makenewzGeneric(tr, q, s, defaultArray, iterations, zqs, FALSE);                  
    makenewzGeneric(tr, r, s, defaultArray, iterations, zrs, FALSE);


    for(i = 0; i < numBranches; i++)
    {
      lzqr = (zqr[i] > zmin) ? log(zqr[i]) : log(zmin); 
      lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
      lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
      lzsum = 0.5 * (lzqr + lzqs + lzrs);

      lzq = lzsum - lzrs;
      lzr = lzsum - lzqs;
      lzs = lzsum - lzqr;
      lzmax = log(zmax);

      if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
      else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
      else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          

      e1[i] = exp(lzq);
      e2[i] = exp(lzr);
      e3[i] = exp(lzs);
    }
    hookup(p->next,       q, e1, numBranches);
    hookup(p->next->next, r, e2, numBranches);
    hookup(p,             s, e3, numBranches);      		  
  }
  else
  {       
    double  z[NUM_BRANCHES]; 

    for(i = 0; i < numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      

      if(z[i] < zmin) 
        z[i] = zmin;
      if(z[i] > zmax)
        z[i] = zmax;
    }

    hookup(p->next,       q, z, tr->numBranches);
    hookup(p->next->next, r, z, tr->numBranches);	                         
  }

  newviewGeneric(tr, p, FALSE);

  if(tr->thoroughInsertion)
  {     
    localSmooth(tr, p, MAX_LOCAL_SMOOTHING_ITERATIONS);   
    for(i = 0; i < numBranches; i++)
    {
      tr->lzq[i] = p->next->z[i];
      tr->lzr[i] = p->next->next->z[i];
      tr->lzs[i] = p->z[i];            
    }
  }           

  return  TRUE;
}

boolean insertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;

  r = q->back;
  s = p->back;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
    hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
    hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
  }
  else
  {       
    double  z[NUM_BRANCHES];
    int i;

    for(i = 0; i < tr->numBranches; i++)
    {
      double zz;
      zz = sqrt(q->z[i]);     
      if(zz < zmin) 
        zz = zmin;
      if(zz > zmax)
        zz = zmax;
      z[i] = zz;
    }

    hookup(p->next,       q, z, tr->numBranches);
    hookup(p->next->next, r, z, tr->numBranches);
  }   

  newviewGeneric(tr, p, FALSE);

  return  TRUE;
}


static void restoreTopologyOnly(tree *tr, bestlist *bt)
{ 
  nodeptr p = tr->removeNode;
  nodeptr q = tr->insertNode;
  double qz[NUM_BRANCHES], pz[NUM_BRANCHES], p1z[NUM_BRANCHES], p2z[NUM_BRANCHES];
  nodeptr p1, p2, r, s;
  double currentLH = tr->likelihood;
  int i;

  p1 = p->next->back;
  p2 = p->next->next->back;

  for(i = 0; i < tr->numBranches; i++)
  {
    p1z[i] = p1->z[i];
    p2z[i] = p2->z[i];
  }

  hookup(p1, p2, tr->currentZQR, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;             
  for(i = 0; i < tr->numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];           
  }

  r = q->back;
  s = p->back;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
    hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
    hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
  }
  else
  { 	
    double  z[NUM_BRANCHES];	
    for(i = 0; i < tr->numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      
      if(z[i] < zmin)
        z[i] = zmin;
      if(z[i] > zmax)
        z[i] = zmax;
    }
    hookup(p->next,       q, z, tr->numBranches);
    hookup(p->next->next, r, z, tr->numBranches);
  }     

  tr->likelihood = tr->bestOfNode;

  saveBestTree(bt, tr);

  tr->likelihood = currentLH;

  hookup(q, r, qz, tr->numBranches);

  p->next->next->back = p->next->back = (nodeptr) NULL;

  if(tr->thoroughInsertion)    
    hookup(p, s, pz, tr->numBranches);          

  hookup(p->next,       p1, p1z, tr->numBranches); 
  hookup(p->next->next, p2, p2z, tr->numBranches);      
}


boolean testInsertBIG (tree *tr, nodeptr p, nodeptr q)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;
  boolean doIt = TRUE;
  double startLH = tr->endLH;
  int i;

  r = q->back; 
  for(i = 0; i < tr->numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }



  if(doIt)
  {     
    if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;         

    evaluateGeneric(tr, p->next->next, FALSE);       

    if(tr->likelihood > tr->bestOfNode)
    {
      tr->bestOfNode = tr->likelihood;
      tr->insertNode = q;
      tr->removeNode = p;   
      for(i = 0; i < tr->numBranches; i++)
      {
        tr->currentZQR[i] = tr->zqr[i];           
        tr->currentLZR[i] = tr->lzr[i];
        tr->currentLZQ[i] = tr->lzq[i];
        tr->currentLZS[i] = tr->lzs[i];      
      }
    }

    if(tr->likelihood > tr->endLH)
    {			  
      tr->insertNode = q;
      tr->removeNode = p;   
      for(i = 0; i < tr->numBranches; i++)
        tr->currentZQR[i] = tr->zqr[i];      
      tr->endLH = tr->likelihood;                      
    }        

    hookup(q, r, qz, tr->numBranches);

    p->next->next->back = p->next->back = (nodeptr) NULL;

    if(tr->thoroughInsertion)
    {
      nodeptr s = p->back;
      hookup(p, s, pz, tr->numBranches);      
    } 

    if((tr->doCutoff) && (tr->likelihood < startLH))
    {
      tr->lhAVG += (startLH - tr->likelihood);
      tr->lhDEC++;
      if((startLH - tr->likelihood) >= tr->lhCutoff)
        return FALSE;	    
      else
        return TRUE;
    }
    else
      return TRUE;
  }
  else
    return TRUE;  
}




void addTraverseBIG(tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
  {              
    if (! testInsertBIG(tr, p, q))  return;

  }

  if ((!isTip(q->number, tr->mxtips)) && (--maxtrav > 0)) 
  {    
    addTraverseBIG(tr, p, q->next->back, mintrav, maxtrav);
    addTraverseBIG(tr, p, q->next->next->back, mintrav, maxtrav);    
  }
} 





int rearrangeBIG(tree *tr, nodeptr p, int mintrav, int maxtrav)   
{  
  double   p1z[NUM_BRANCHES], p2z[NUM_BRANCHES], q1z[NUM_BRANCHES], q2z[NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;  
  boolean doP = TRUE, doQ = TRUE;

  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
  q = p->back;




  if (!isTip(p->number, tr->mxtips) && doP) 
  {     
    p1 = p->next->back;
    p2 = p->next->next->back;


    if(!isTip(p1->number, tr->mxtips) || !isTip(p2->number, tr->mxtips))
    {
      for(i = 0; i < tr->numBranches; i++)
      {
        p1z[i] = p1->z[i];
        p2z[i] = p2->z[i];	   	   
      }

      if (! removeNodeBIG(tr, p,  tr->numBranches)) return badRear;

      if (!isTip(p1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, p, p1->next->back,
            mintrav, maxtrav);         

        addTraverseBIG(tr, p, p1->next->next->back,
            mintrav, maxtrav);          
      }

      if (!isTip(p2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, p, p2->next->back,
            mintrav, maxtrav);
        addTraverseBIG(tr, p, p2->next->next->back,
            mintrav, maxtrav);          
      }

      hookup(p->next,       p1, p1z, tr->numBranches); 
      hookup(p->next->next, p2, p2z, tr->numBranches);	   	    	    
      newviewGeneric(tr, p, FALSE);	   	    
    }
  }  

  if (!isTip(q->number, tr->mxtips) && maxtrav > 0 && doQ) 
  {
    q1 = q->next->back;
    q2 = q->next->next->back;

    /*if (((!q1->tip) && (!q1->next->back->tip || !q1->next->next->back->tip)) ||
      ((!q2->tip) && (!q2->next->back->tip || !q2->next->next->back->tip))) */
    if (
        (
         ! isTip(q1->number, tr->mxtips) && 
         (! isTip(q1->next->back->number, tr->mxtips) || ! isTip(q1->next->next->back->number, tr->mxtips))
        )
        ||
        (
         ! isTip(q2->number, tr->mxtips) && 
         (! isTip(q2->next->back->number, tr->mxtips) || ! isTip(q2->next->next->back->number, tr->mxtips))
        )
       )
    {

      for(i = 0; i < tr->numBranches; i++)
      {
        q1z[i] = q1->z[i];
        q2z[i] = q2->z[i];
      }

      if (! removeNodeBIG(tr, q, tr->numBranches)) return badRear;

      mintrav2 = mintrav > 2 ? mintrav : 2;

      if (/*! q1->tip*/ !isTip(q1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, q, q1->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, q, q1->next->next->back,
            mintrav2 , maxtrav);         
      }

      if (/*! q2->tip*/ ! isTip(q2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, q, q2->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, q, q2->next->next->back,
            mintrav2 , maxtrav);          
      }	   

      hookup(q->next,       q1, q1z, tr->numBranches); 
      hookup(q->next->next, q2, q2z, tr->numBranches);

      newviewGeneric(tr, q, FALSE); 	   
    }
  } 

  return  1;
} 





static double treeOptimizeRapid(tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt, infoList *iList)
{
  int i, index,
      *perm = (int*)NULL;   

  nodeRectifier(tr);



  if (maxtrav > tr->mxtips - 3)  
    maxtrav = tr->mxtips - 3;  



  resetInfoList(iList);

  resetBestTree(bt);

  tr->startLH = tr->endLH = tr->likelihood;

  if(tr->doCutoff)
  {
    if(tr->bigCutoff)
    {	  
      if(tr->itCount == 0)    
        tr->lhCutoff = 0.5 * (tr->likelihood / -1000.0);    
      else    		 
        tr->lhCutoff = 0.5 * ((tr->lhAVG) / ((double)(tr->lhDEC))); 	  
    }
    else
    {
      if(tr->itCount == 0)    
        tr->lhCutoff = tr->likelihood / -1000.0;    
      else    		 
        tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));   
    }    

    tr->itCount = tr->itCount + 1;
    tr->lhAVG = 0;
    tr->lhDEC = 0;
  }

  /*
     printf("DoCutoff: %d\n", tr->doCutoff);
     printf("%d %f %f %f\n", tr->itCount, tr->lhAVG, tr->lhDEC, tr->lhCutoff);

     printf("%d %d\n", mintrav, maxtrav);
     */

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
  {           
    tr->bestOfNode = unlikely;          

    if(adef->permuteTreeoptimize)
      index = perm[i];
    else
      index = i;     

    if(rearrangeBIG(tr, tr->nodep[index], mintrav, maxtrav))
    {    
      if(tr->thoroughInsertion)
      {
        if(tr->endLH > tr->startLH)                 	
        {			   	     
          restoreTreeFast(tr);	 	 
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr);
        }
        else
        { 		  
          if(tr->bestOfNode != unlikely)		    	     
            restoreTopologyOnly(tr, bt);		    
        }	   
      }
      else
      {
        insertInfoList(tr->nodep[index], tr->bestOfNode, iList);	    
        if(tr->endLH > tr->startLH)                 	
        {		      
          restoreTreeFast(tr);	  	      
          tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  
        }	    	  
      }
    }     
  }     

  if(!tr->thoroughInsertion)
  {           
    tr->thoroughInsertion = TRUE;  

    for(i = 0; i < iList->valid; i++)
    { 	  
      tr->bestOfNode = unlikely;

      if(rearrangeBIG(tr, iList->list[i].node, mintrav, maxtrav))
      {	  
        if(tr->endLH > tr->startLH)                 	
        {	 	     
          restoreTreeFast(tr);	 	 
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr);
        }
        else
        { 

          if(tr->bestOfNode != unlikely)
          {	     
            restoreTopologyOnly(tr, bt);
          }	
        }      
      }
    }       

    tr->thoroughInsertion = FALSE;
  }

  if(adef->permuteTreeoptimize)
    free(perm);

  return tr->startLH;     
}




boolean testInsertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{    
  if(tr->thoroughInsertion)
  {
    if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;    

    evaluateGeneric(tr, p->next->next, FALSE);               
  }
  else
  {
    if (! insertRestoreBIG(tr, p, q))       return FALSE;

    {
      nodeptr x, y;
      x = p->next->next;
      y = p->back;

      if(! isTip(x->number, tr->mxtips) && isTip(y->number, tr->mxtips))
      {
        while ((! x->x)) 
        {
          if (! (x->x))
            newviewGeneric(tr, x, FALSE);		     
        }
      }

      if(isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! y->x)) 
        {		  
          if (! (y->x))
            newviewGeneric(tr, y, FALSE);
        }
      }

      if(!isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! x->x) || (! y->x)) 
        {
          if (! (x->x))
            newviewGeneric(tr, x, FALSE);
          if (! (y->x))
            newviewGeneric(tr, y, FALSE);
        }
      }				      	

    }

    tr->likelihood = tr->endLH;
  }

  return TRUE;
} 

void restoreTreeFast(tree *tr)
{
  removeNodeRestoreBIG(tr, tr->removeNode);    
  testInsertRestoreBIG(tr, tr->removeNode, tr->insertNode);
}

static void myfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, stream);

  assert(bytes_written == nmemb);
}

static void myfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t
    bytes_read;

  bytes_read = fread(ptr, size, nmemb, stream);

  assert(bytes_read == nmemb);
}




static void writeTree(tree *tr, FILE *f)
{
  int 
    x = tr->mxtips + 3 * (tr->mxtips - 1);

  nodeptr
    base = tr->nodeBaseAddress;

  myfwrite(&(tr->start->number), sizeof(int), 1, f);
  myfwrite(&base, sizeof(nodeptr), 1, f);
  myfwrite(tr->nodeBaseAddress, sizeof(node), x, f);

}

int ckpCount = 0;

static void writeCheckpoint(tree *tr, int state)
{
  int   
    model; 

  char 
    extendedName[2048],
    buf[64];

  FILE 
    *f;

  strcpy(extendedName,  binaryCheckpointName);
  strcat(extendedName, "_");
  sprintf(buf, "%d", ckpCount);
  strcat(extendedName, buf);  

  ckpCount++;

  f = myfopen(extendedName, "w"); 

  /* cdta */   


  tr->ckp.accumulatedTime = accumulatedTime + (gettime() - masterTime);

  tr->ckp.state = state;

  tr->ckp.tr_optimizeRateCategoryInvocations = tr->optimizeRateCategoryInvocations;
  tr->ckp.tr_thoroughInsertion = tr->thoroughInsertion;
  tr->ckp.tr_startLH  = tr->startLH;
  tr->ckp.tr_endLH    = tr->endLH;
  tr->ckp.tr_likelihood = tr->likelihood;
  tr->ckp.tr_bestOfNode = tr->bestOfNode;

  tr->ckp.tr_lhCutoff = tr->lhCutoff;
  tr->ckp.tr_lhAVG    = tr->lhAVG;
  tr->ckp.tr_lhDEC    = tr->lhDEC;     
  tr->ckp.tr_itCount  = tr->itCount;
  tr->ckp.tr_doCutoff = tr->doCutoff;
  /* printf("Acc time: %f\n", tr->ckp.accumulatedTime); */

  /* user stupidity */


  tr->ckp.searchConvergenceCriterion = tr->searchConvergenceCriterion;
  tr->ckp.rateHetModel =  tr->rateHetModel;
  tr->ckp.maxCategories =  tr->maxCategories;
  tr->ckp.NumberOfModels = tr->NumberOfModels;
  tr->ckp.numBranches = tr->numBranches;
  tr->ckp.originalCrunchedLength = tr->originalCrunchedLength;
  tr->ckp.mxtips = tr->mxtips;
  strcpy(tr->ckp.seq_file, seq_file);

  /* handle user stupidity */


  myfwrite(&(tr->ckp), sizeof(checkPointState), 1, f);

  myfwrite(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfwrite(tr->tree1, sizeof(char), tr->treeStringLength, f);

  myfwrite(tr->rateCategory, sizeof(int), tr->originalCrunchedLength, f);
  myfwrite(tr->patrat, sizeof(double), tr->originalCrunchedLength, f);
  myfwrite(tr->patratStored, sizeof(double), tr->originalCrunchedLength, f);

  /* need to store this as well in checkpoints, otherwise the branch lengths 
     in the output tree files will be wrong, not the internal branch lengths though */

  myfwrite(tr->fracchanges,  sizeof(double), tr->NumberOfModels, f);
  myfwrite(&(tr->fracchange),   sizeof(double), 1, f);


  /* pInfo */

  for(model = 0; model < tr->NumberOfModels; model++)
  {
    int 
      dataType = tr->partitionData[model].dataType;

    myfwrite(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);
    myfwrite(tr->partitionData[model].perSiteRates, sizeof(double), tr->maxCategories, f);
    myfwrite(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfwrite(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfwrite(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);  

    myfwrite(tr->partitionData[model].frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfwrite(tr->partitionData[model].tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);  
    myfwrite(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);    
    myfwrite(&(tr->partitionData[model].alpha), sizeof(double), 1, f);
  }



  writeTree(tr, f);

  fclose(f); 

  printBothOpen("\nCheckpoint written to: %s likelihood: %f\n", extendedName, tr->likelihood);
}

static void readTree(tree *tr, FILE *f)
{
  int 
    nodeNumber,   
    x = tr->mxtips + 3 * (tr->mxtips - 1);





  nodeptr
    startAddress;

  myfread(&nodeNumber, sizeof(int), 1, f);

  tr->start = tr->nodep[nodeNumber];

  /*printf("Start: %d %d\n", tr->start->number, nodeNumber);*/

  myfread(&startAddress, sizeof(nodeptr), 1, f);

  /*printf("%u %u\n", (size_t)startAddress, (size_t)tr->nodeBaseAddress);*/



  myfread(tr->nodeBaseAddress, sizeof(node), x, f);

  {
    int i;    

    size_t         
      offset;

    boolean 
      addIt;

    if(startAddress > tr->nodeBaseAddress)
    {
      addIt = FALSE;
      offset = (size_t)startAddress - (size_t)tr->nodeBaseAddress;
    }
    else
    {
      addIt = TRUE;
      offset = (size_t)tr->nodeBaseAddress - (size_t)startAddress;
    }       

    for(i = 0; i < x; i++)
    {      	
      if(addIt)
      {	    
        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next + offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back + offset);
      }
      else
      {

        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next - offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back - offset);	   
      } 
    }

  }

  evaluateGeneric(tr, tr->start, TRUE);  

  printBothOpen("RAxML Restart with likelihood: %1.50f\n", tr->likelihood);
}


static void readCheckpoint(tree *tr)
{
  int  
    restartErrors = 0,
                  model; 

  FILE 
    *f = myfopen(binaryCheckpointInputName, "r");

  /* cdta */   

  myfread(&(tr->ckp), sizeof(checkPointState), 1, f);



  if(tr->ckp.searchConvergenceCriterion != tr->searchConvergenceCriterion)
  {
    printf("restart error, you are trying to re-start a run where the ML search criterion was turned %s\n", (tr->ckp.searchConvergenceCriterion)?"ON":"OFF");
    restartErrors++;
  }  

  if(tr->ckp.rateHetModel !=  tr->rateHetModel)
  {
    printf("restart error, you are trying to re-start a run with a different model of rate heterogeneity, the checkpoint was obtained under: %s\n", (tr->ckp.rateHetModel == GAMMA)?"GAMMA":"PSR");
    restartErrors++;
  }  

  if(tr->ckp.maxCategories !=  tr->maxCategories)
  {
    printf("restart error, you are trying to re-start a run with %d per-site rate categories, the checkpoint was obtained with: %d\n", tr->maxCategories, tr->ckp.maxCategories);
    restartErrors++;
  }

  if(tr->ckp.NumberOfModels != tr->NumberOfModels)
  {
    printf("restart error, you are trying to re-start a run with %d partitions, the checkpoint was obtained with: %d partitions\n", tr->NumberOfModels, tr->ckp.NumberOfModels);
    restartErrors++;      
  }

  if(tr->ckp.numBranches != tr->numBranches)
  {
    printf("restart error, you are trying to re-start a run where independent per-site branch length estimates were turned %s\n", (tr->ckp.numBranches > 1)?"ON":"OFF");
    restartErrors++;
  }

  if(tr->ckp.originalCrunchedLength != tr->originalCrunchedLength)
  {
    printf("restart error, you are trying to re-start a run with %d site patterns, the checkpoint was obtained with: %d site patterns\n", tr->ckp.originalCrunchedLength, tr->originalCrunchedLength);
    restartErrors++; 
  }

  if(tr->ckp.mxtips != tr->mxtips)
  {
    printf("restart error, you are trying to re-start a run with %d taxa, the checkpoint was obtained with: %d taxa\n", tr->mxtips, tr->ckp.mxtips);
    restartErrors++; 
  }

  if(strcmp(tr->ckp.seq_file, seq_file) != 0)
  {
    printf("restart error, you are trying to re-start from alignemnt file %s, the checkpoint was obtained with file: %s\n", tr->ckp.seq_file, seq_file);
    restartErrors++; 
  }

  printf("REstart errors: %d\n", restartErrors);

  if(restartErrors > 0)
  {
    printf("User induced errors with the restart from checkpoint, exiting ...\n");

    if(restartErrors > 4)
      printf(" ... maybe you should do field work instead of trying to use a computer ...\n");
    if(restartErrors > 6)
      printf(" ... kala eisai telios ilithios;\n");

    exit(-1);
  }

  tr->ntips = tr->mxtips;

  tr->startLH    = tr->ckp.tr_startLH;
  tr->endLH      = tr->ckp.tr_endLH;
  tr->likelihood = tr->ckp.tr_likelihood;
  tr->bestOfNode = tr->ckp.tr_bestOfNode;

  tr->lhCutoff   = tr->ckp.tr_lhCutoff;
  tr->lhAVG      = tr->ckp.tr_lhAVG;
  tr->lhDEC      = tr->ckp.tr_lhDEC;
  tr->itCount    = tr->ckp.tr_itCount;
  tr->thoroughInsertion       = tr->ckp.tr_thoroughInsertion;



  accumulatedTime = tr->ckp.accumulatedTime;

  /* printf("Accumulated time so far: %f\n", accumulatedTime); */

  tr->optimizeRateCategoryInvocations = tr->ckp.tr_optimizeRateCategoryInvocations;


  myfread(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfread(tr->tree1, sizeof(char), tr->treeStringLength, f);

  if(tr->searchConvergenceCriterion)
  {
    int bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 0) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 0))
    { 

#ifdef _DEBUG_CHECKPOINTING    
      printf("parsing Tree 0\n");
#endif

      treeReadTopologyString(tr->tree0, tr);   

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 0, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }

    bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 1) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 1))
    {

#ifdef _DEBUG_CHECKPOINTING
      printf("parsing Tree 1\n");
#endif

      treeReadTopologyString(tr->tree1, tr); 

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 1, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }
  }

  myfread(tr->rateCategory, sizeof(int), tr->originalCrunchedLength, f);
  myfread(tr->patrat, sizeof(double), tr->originalCrunchedLength, f);
  myfread(tr->patratStored, sizeof(double), tr->originalCrunchedLength, f);


  /* need to read this as well in checkpoints, otherwise the branch lengths 
     in the output tree files will be wrong, not the internal branch lengths though */

  myfread(tr->fracchanges,  sizeof(double), tr->NumberOfModels, f);
  myfread(&(tr->fracchange),   sizeof(double), 1, f);

  /* pInfo */

  for(model = 0; model < tr->NumberOfModels; model++)
  {
    int 
      dataType = tr->partitionData[model].dataType;

    myfread(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);
    myfread(tr->partitionData[model].perSiteRates, sizeof(double), tr->maxCategories, f);
    myfread(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfread(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfread(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);  

    myfread(tr->partitionData[model].frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfread(tr->partitionData[model].tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);  
    myfread(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);  
    myfread(&(tr->partitionData[model].alpha), sizeof(double), 1, f);

    makeGammaCats(tr->partitionData[model].alpha, tr->partitionData[model].gammaRates, 4, tr->useMedian);
  }

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  masterBarrier(THREAD_COPY_INIT_MODEL, tr);
#endif

  updatePerSiteRates(tr, FALSE);  

  readTree(tr, f);

  fclose(f); 

}

static void restoreTreeDataValuesFromCheckpoint(tree *tr)
{
  tr->optimizeRateCategoryInvocations = tr->ckp.tr_optimizeRateCategoryInvocations;  
  tr->thoroughInsertion = tr->ckp.tr_thoroughInsertion;
  tr->likelihood = tr->ckp.tr_likelihood;              
  tr->lhCutoff = tr->ckp.tr_lhCutoff;
  tr->lhAVG    = tr->ckp.tr_lhAVG;
  tr->lhDEC    = tr->ckp.tr_lhDEC;   	 
  tr->itCount = tr->ckp.tr_itCount;
  tr->doCutoff = tr->ckp.tr_doCutoff;
}

void restart(tree *tr)
{  
  readCheckpoint(tr);

  switch(tr->ckp.state)
  {
    case REARR_SETTING:      
      break;
    case FAST_SPRS:
      break;
    case SLOW_SPRS:
      break;
    default:
      assert(0);
  }
}

int determineRearrangementSetting(tree *tr,  analdef *adef, bestlist *bestT, bestlist *bt)
{
  const 
    int MaxFast = 26;

  int 
    i,   
    maxtrav = 5, 
    bestTrav = 5;

  double 
    startLH = tr->likelihood; 

  boolean 
    impr   = TRUE,
           cutoff = tr->doCutoff;

  if(adef->useCheckpoint)
  {
    assert(tr->ckp.state == REARR_SETTING);

    maxtrav = tr->ckp.maxtrav;
    bestTrav = tr->ckp.bestTrav;
    startLH  = tr->ckp.startLH;
    impr     = tr->ckp.impr;      
    cutoff = tr->ckp.cutoff;

    adef->useCheckpoint = FALSE;
  }

  tr->doCutoff = FALSE;      

  resetBestTree(bt);    

#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("MAXTRAV: %d\n", maxtrav);
#endif

  while(impr && maxtrav < MaxFast)
  {	
    recallBestTree(bestT, 1, tr);     
    nodeRectifier(tr);                      

    {
      tr->ckp.cutoff = cutoff;	
      tr->ckp.maxtrav = maxtrav;
      tr->ckp.bestTrav = bestTrav;
      tr->ckp.startLH  = startLH;
      tr->ckp.impr = impr;

      writeCheckpoint(tr, REARR_SETTING);    
    }

    if (maxtrav > tr->mxtips - 3)  
      maxtrav = tr->mxtips - 3;    

    tr->startLH = tr->endLH = tr->likelihood;

    for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {                	         
      tr->bestOfNode = unlikely;

      if(rearrangeBIG(tr, tr->nodep[i], 1, maxtrav))
      {	     
        if(tr->endLH > tr->startLH)                 	
        {		 	 	      
          restoreTreeFast(tr);	        	  	 	  	      
          tr->startLH = tr->endLH = tr->likelihood;		 
        }	         	       	
      }
    }

    treeEvaluate(tr, 8 ); // 32 * 0.25 
    saveBestTree(bt, tr); 

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("TRAV: %d lh %f\n", maxtrav, tr->likelihood);
#endif

    if(tr->likelihood > startLH)
    {	 
      startLH = tr->likelihood; 	  	  	  
      printLog(tr);	  
      bestTrav = maxtrav;	 
      impr = TRUE;
    }
    else	
      impr = FALSE;	



    if(tr->doCutoff)
    {
      tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));       

      tr->itCount =  tr->itCount + 1;
      tr->lhAVG = 0;
      tr->lhDEC = 0;
    }

    maxtrav += 5;


  }

  recallBestTree(bt, 1, tr);

  tr->doCutoff = cutoff; 

#ifdef _DEBUG_CHECKPOINTING
  printBothOpen("BestTrav %d\n", bestTrav);
#endif

  return bestTrav;     
}





void computeBIGRAPID (tree *tr, analdef *adef, boolean estimateModel) 
{   
  int
    i,
    impr, 
    bestTrav = 0, 
    rearrangementsMax = 0, 
    rearrangementsMin = 0,    
    thoroughIterations = 0,
    fastIterations = 0;

  double 
    lh = unlikely, 
       previousLh = unlikely, 
       difference, 
       epsilon;              

  bestlist 
    *bestT, 
    *bt;    

  infoList 
    *iList = (infoList*)malloc(sizeof(infoList));

  /* now here is the RAxML hill climbing search algorithm */


  /* initialize two lists of size 1 and size 20 that will keep track of the best 
     and 20 best tree topologies respectively */

  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);

  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  /* initialize an additional data structure used by the search algo, all of this is pretty 
     RAxML-specific and should probably not be in the library */

  initInfoList(iList, 50);

  /* some pretty atbitrary thresholds */

  difference = 10.0;
  epsilon = 0.01;    

  /* Thorough = 0 means that we will do fast SPR inbsertions without optimizing the 
     three branches adjacent to the subtree insertion position via Newton-Raphson 
     */

  tr->thoroughInsertion = FALSE;     

  /* if we are not using a checkpoint and estimateModel is set to TRUE we call the function 
     that optimizes model parameters, such as the CAT model assignment, the alpha paremeter
     or the rates in the GTR matrix. Otherwise we just optimize the branch lengths. Note that 
     the second parameter of treeEvaluate() controls how many times we will iterate over all branches 
     of the tree until we give up, provided that, the br-len opt. has not converged before.
     */

  if(!adef->useCheckpoint)
  {
    if(estimateModel)
      modOpt(tr, 10.0);
    else
      treeEvaluate(tr, 64); // 32 * 2
  }

  /* print some stuff to the RAxML_log file */

  printLog(tr); 

  /* save the current tree (which is the input tree parsed via -t in the bestT list */

  saveBestTree(bestT, tr);

  /* if the rearrangmenet radius has been set by the user ie. adef->initailSet == TRUE 
     then just set the apppropriate parameter.
     Otherwise, call the function  determineRearrangementSetting() that seeks 
     for the best radius by executing SPR moves on the initial tree with different radii
     and returns the smallest radius that yields the best log likelihood score after 
     applying one cycle of SPR moves to the tree 
     */

  if(!adef->initialSet)   
  {
    if((!adef->useCheckpoint) || (adef->useCheckpoint && tr->ckp.state == REARR_SETTING))
    {
      bestTrav = adef->bestTrav = determineRearrangementSetting(tr, adef, bestT, bt);     	  
      printBothOpen("\nBest rearrangement radius: %d\n", bestTrav);
    }
  }
  else
  {
    bestTrav = adef->bestTrav = adef->initial;       
    printBothOpen("\nUser-defined rearrangement radius: %d\n", bestTrav);
  }


  /* some checkpointing noise */
  if(!(adef->useCheckpoint && (tr->ckp.state == FAST_SPRS || tr->ckp.state == SLOW_SPRS)))
  {      

    /* optimize model params more thoroughly or just optimize branch lengths */
    if(estimateModel)
      modOpt(tr, 5.0);
    else
      treeEvaluate(tr, 32);   // 32 * 1 
  }

  /* save the current tree again, while the topology has not changed, the branch lengths have changed in the meantime, hence
     we need to store them again */

  saveBestTree(bestT, tr); 

  /* set the loop variable to TRUE */

  impr = 1;

  /* this is for the additional RAxML heuristics described imn this paper here:

     A. Stamatakis,  F. Blagojevic, C.D. Antonopoulos, D.S. Nikolopoulos: "Exploring new Search Algorithms and Hardware for Phylogenetics: RAxML meets the IBM Cell". 
     In Journal of VLSI Signal Processing Systems, 48(3):271-286, 2007.

     This is turned on by default 
     */


  if(tr->doCutoff)
    tr->itCount = 0;

  /* figure out where to continue computations if we restarted from a checkpoint */

  if(adef->useCheckpoint && tr->ckp.state == FAST_SPRS)
    goto START_FAST_SPRS;

  if(adef->useCheckpoint && tr->ckp.state == SLOW_SPRS)
    goto START_SLOW_SPRS;

  while(impr)
  {              
START_FAST_SPRS:
    /* if re-starting from checkpoint set the required variable values to the 
       values that they had when the checkpoint was written */

    if(adef->useCheckpoint && tr->ckp.state == FAST_SPRS)
    {	    
      impr = tr->ckp.impr;	  
      bestTrav = tr->ckp.bestTrav;	  
      rearrangementsMax = tr->ckp.rearrangementsMax;
      rearrangementsMin = tr->ckp.rearrangementsMin;
      thoroughIterations = tr->ckp.thoroughIterations;
      fastIterations = tr->ckp.fastIterations;  
      lh = tr->ckp.lh;
      previousLh = tr->ckp.previousLh;
      difference = tr->ckp.difference;
      epsilon    = tr->ckp.epsilon;  

      restoreTreeDataValuesFromCheckpoint(tr);	  

      adef->useCheckpoint = FALSE;
    }
    else
      /* otherwise, restore the currently best tree */
      recallBestTree(bestT, 1, tr); 

    /* save states of algorithmic/heuristic variables for printing the next checkpoint */


    tr->ckp.impr = impr;	
    tr->ckp.bestTrav = bestTrav;      
    tr->ckp.rearrangementsMax = rearrangementsMax;
    tr->ckp.rearrangementsMin = rearrangementsMin;
    tr->ckp.thoroughIterations = thoroughIterations;
    tr->ckp.fastIterations = fastIterations;  
    tr->ckp.lh = lh;
    tr->ckp.previousLh = previousLh;
    tr->ckp.difference = difference;
    tr->ckp.epsilon    = epsilon;              
    tr->ckp.bestTrav = bestTrav;       
    tr->ckp.impr = impr;                 

    /* write a binary checkpoint */
    writeCheckpoint(tr, FAST_SPRS);      

    /* this is the aforementioned convergence criterion that requires computing the RF,
       let's not worry about this right now */

    if(tr->searchConvergenceCriterion)
    {
      int bCounter = 0;	  	      	 	  	  	

      if(fastIterations > 1)
        cleanupHashTable(tr->h, (fastIterations % 2));		
      
      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, fastIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE, tr->threadID);	    

      {
        char 
          *buffer = (char*)calloc((size_t)tr->treeStringLength, sizeof(char));
#ifdef _DEBUG_CHECKPOINTING
        printf("Storing tree in slot %d\n", fastIterations % 2);
#endif

        Tree2String(buffer, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);

        if(fastIterations % 2 == 0)	      
          memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
        else
          memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    

        free(buffer);
      }


      assert(bCounter == tr->mxtips - 3);	    	   

      if(fastIterations > 0)
      {
        double rrf = convergenceCriterion(tr->h, tr->mxtips);

        if(rrf <= 0.01) /* 1% cutoff */
        {
          printBothOpen("ML fast search converged at fast SPR cycle %d with stopping criterion\n", fastIterations);
          printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
          cleanupHashTable(tr->h, 0);
          cleanupHashTable(tr->h, 1);
          goto cleanup_fast;
        }
        else		    
          printBothOpen("ML search convergence criterion fast cycle %d->%d Relative Robinson-Foulds %f\n", fastIterations - 1, fastIterations, rrf);
      }
    }


    /* count how many fast iterations with so-called fast SPR moves we have executed */

    fastIterations++;	

    /* optimize branch lengths */

    treeEvaluate(tr, 32);  // 32 * 1 = 32 

    /* save the tree with those branch lengths again */

    saveBestTree(bestT, tr);           

    /* print the log likelihood */

    printLog(tr);    

    /* print this intermediate tree to file */

    printResult(tr, adef, FALSE);    

    /* update the current best likelihood */

    lh = previousLh = tr->likelihood;

    /* in here we actually do a cycle of SPR moves */

    treeOptimizeRapid(tr, 1, bestTrav, adef, bt, iList);   

    /* set impr to 0 since in the immediately following for loop we check if the SPR moves above have generated 
       a better tree */

    impr = 0;

    /* loop over the 20 best trees generated by the fast SPR moves, and check if they improve the likelihood after all of their branch lengths
       have been optimized */

    for(i = 1; i <= bt->nvalid; i++)
    {	    	
      /* restore tree i from list generated by treeOptimizeRapid */

      recallBestTree(bt, i, tr);

      /* optimize branch lengths of this tree */

      treeEvaluate(tr, 8); // 0.25 * 32

      /* calc. the likelihood improvement */

      difference = ((tr->likelihood > previousLh)? 
          tr->likelihood - previousLh: 
          previousLh - tr->likelihood); 	    

      /* if the likelihood has improved save the current tree as best tree and continue */
      /* note that we always compre this tree to the likelihood of the previous best tree */

      if(tr->likelihood > lh && difference > epsilon)
      {
        impr = 1;	       
        lh = tr->likelihood;	       	     
        saveBestTree(bestT, tr);

      }	   	   
    }
#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("FAST LH: %f\n", lh);
#endif


  }

  /* needed for this RF-based convergence criterion that I actually describe in here:

     A. Stamatakis: "Phylogenetic Search Algorithms for Maximum Likelihood". In M. Elloumi, A.Y. Zomaya, editors. 
     Algorithms in Computational Biology: techniques, Approaches and Applications, John Wiley and Sons

     a copy of this book is in my office */

  if(tr->searchConvergenceCriterion)
  {
    cleanupHashTable(tr->h, 0);
    cleanupHashTable(tr->h, 1);
  }

cleanup_fast:  
  /*
     now we have jumped out of the loop that executes 
     fast SPRs, and next we will execute a loop that executes throough SPR cycles (with SPR moves 
     that optimize via newton-Raphson all adjacent branches to the insertion point) 
     until no through SPR move can be found that improves the likelihood further. A classic 
     hill climbing algo.
     */

  tr->thoroughInsertion = TRUE;
  impr = 1;

  /* restore the currently best tree. this si actually required, because we do not know which tree
     is actually stored in the tree data structure when the above loop exits */

  recallBestTree(bestT, 1, tr); 

  {
    /* RE-TRAVERSE THE ENTIRE TREE */

    evaluateGeneric(tr, tr->start, TRUE);
#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("After Fast SPRs Final %f\n", tr->likelihood);   
#endif
  }

  /* optimize model params (including branch lengths) or just 
     optimize branch lengths and leave the other model parameters (GTR rates, alhpa) 
     alone */

  if(estimateModel)
    modOpt(tr, 1.0);
  else
    treeEvaluate(tr, 32 ); //32 * 1

  /* start loop that executes thorough SPR cycles */

  while(1)
  {	 
    /* once again if we want to restart from a checkpoint that was written during this loop we need
       to restore the values of the variables appropriately */
START_SLOW_SPRS:
    if(adef->useCheckpoint && tr->ckp.state == SLOW_SPRS)
    {	        
      impr = tr->ckp.impr;	 
      bestTrav = tr->ckp.bestTrav;	 
      rearrangementsMax = tr->ckp.rearrangementsMax;
      rearrangementsMin = tr->ckp.rearrangementsMin;
      thoroughIterations = tr->ckp.thoroughIterations;
      fastIterations = tr->ckp.fastIterations;     
      lh = tr->ckp.lh;
      previousLh = tr->ckp.previousLh;
      difference = tr->ckp.difference;
      epsilon    = tr->ckp.epsilon;                    

      restoreTreeDataValuesFromCheckpoint(tr);	  

      adef->useCheckpoint = FALSE;
    }
    else
      /* otherwise we restore the currently best tree and load it from bestT into our tree data 
         structuire tr */
      recallBestTree(bestT, 1, tr);

    /* now, we write a checkpoint */

    tr->ckp.impr = impr;      
    tr->ckp.bestTrav = bestTrav;      
    tr->ckp.rearrangementsMax = rearrangementsMax;
    tr->ckp.rearrangementsMin = rearrangementsMin;
    tr->ckp.thoroughIterations = thoroughIterations;
    tr->ckp.fastIterations = fastIterations;	  
    tr->ckp.lh = lh;
    tr->ckp.previousLh = previousLh;
    tr->ckp.difference = difference;
    tr->ckp.epsilon    = epsilon;              
    tr->ckp.bestTrav = bestTrav;       
    tr->ckp.impr = impr;                 

    /* write binary checkpoint to file */

    writeCheckpoint(tr, SLOW_SPRS);     

    if(impr)
    {	    
      /* if the logl has improved write out some stuff and adapt the rearrangement radii */
      printResult(tr, adef, FALSE);
      /* minimum rearrangement radius */
      rearrangementsMin = 1;
      /* max radius, this is probably something I need to explain at the whiteboard */
      rearrangementsMax = adef->stepwidth;	

      /* once again the convergence criterion */

      if(tr->searchConvergenceCriterion)
      {
        int bCounter = 0;	      

        if(thoroughIterations > 1)
          cleanupHashTable(tr->h, (thoroughIterations % 2));		

        bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, thoroughIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
            &bCounter, 1, FALSE, FALSE, tr->threadID);	    


        {
          char 
            *buffer = (char*)calloc((size_t)tr->treeStringLength, sizeof(char));

#ifdef _DEBUG_CHECKPOINTING		
          printf("Storing tree in slot %d\n", thoroughIterations % 2);
#endif

          Tree2String(buffer, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, SUMMARIZE_LH, FALSE, FALSE);

          if(thoroughIterations % 2 == 0)	      
            memcpy(tr->tree0, buffer, tr->treeStringLength * sizeof(char));
          else
            memcpy(tr->tree1, buffer, tr->treeStringLength * sizeof(char));	    

          free(buffer);
        }

        assert(bCounter == tr->mxtips - 3);

        if(thoroughIterations > 0)
        {
          double rrf = convergenceCriterion(tr->h, tr->mxtips);

          if(rrf <= 0.01) /* 1% cutoff */
          {
            printBothOpen("ML search converged at thorough SPR cycle %d with stopping criterion\n", thoroughIterations);
            printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
            goto cleanup;
          }
          else		    
            printBothOpen("ML search convergence criterion thorough cycle %d->%d Relative Robinson-Foulds %f\n", thoroughIterations - 1, thoroughIterations, rrf);
        }
      }



      thoroughIterations++;	  
    }			  			
    else
    {
      /* if the lnl has not imrpved by the current SPR cycle adapt the min and max rearrangemnt radii and try again */

      rearrangementsMax += adef->stepwidth;
      rearrangementsMin += adef->stepwidth; 	        	      

      /* if we have already tried them then abandon this loop, the search has converged */
      if(rearrangementsMax > adef->max_rearrange)	     	     	 
        goto cleanup; 	   
    }

    /* optimize branch lengths of best tree */

    treeEvaluate(tr, 32 ); // 32 * 1

    /* do some bokkeeping and printouts again */
    previousLh = lh = tr->likelihood;	      
    saveBestTree(bestT, tr);     
    printLog(tr);

    /* do a cycle of thorough SPR moves with the minimum and maximum rearrangement radii */

    treeOptimizeRapid(tr, rearrangementsMin, rearrangementsMax, adef, bt, iList);

    impr = 0;			      		            

    /* once again get the best 20 trees produced by the SPR cycle, load them from the bt tree list into tr
       optimize their branch lengths and figure out if the LnL of the tree has improved */

    for(i = 1; i <= bt->nvalid; i++)
    {		 
      recallBestTree(bt, i, tr);	 	    	    	

      treeEvaluate(tr, 8 ); // 0.25	* 32    	 

      difference = ((tr->likelihood > previousLh)? 
          tr->likelihood - previousLh: 
          previousLh - tr->likelihood); 	    
      if(tr->likelihood > lh && difference > epsilon)
      {
        impr = 1;	       
        lh = tr->likelihood;	  	     
        saveBestTree(bestT, tr);
      }	   	   
    }  

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("SLOW LH: %f\n", lh);              
#endif
  }

cleanup: 

  /* do a final full tree traversal, not sure if this is required here */

  {
    evaluateGeneric(tr, tr->start, TRUE);

#ifdef _DEBUG_CHECKPOINTING
    printBothOpen("After SLOW SPRs Final %f\n", tr->likelihood);   
#endif
  }

  /* free data structures */

  if(tr->searchConvergenceCriterion)
  {
    freeBitVectors(tr->bitVectors, 2 * tr->mxtips);
    free(tr->bitVectors);
    freeHashTable(tr->h);
    free(tr->h);
  }

  freeBestTree(bestT);
  free(bestT);
  freeBestTree(bt);
  free(bt);
  freeInfoList(iList);  
  free(iList);

  printLog(tr);
  printResult(tr, adef, TRUE);

  /* and we are done, return to main() in axml.c  */

}



/* The number of maximum smoothing iterations is given explicitely */
boolean 
treeEvaluate (tree *tr, int maxSmoothIterations)       /* Evaluate a user tree */
{
  boolean result;
  result = smoothTree(tr, maxSmoothIterations); /* former (32 * smoothFactor) */
  assert(result); 

  evaluateGeneric(tr, tr->start, FALSE);   
  

  return TRUE;
}

/* Perform an NNI move. swap can be either 1 or 2 */
void NNI(tree * tr, nodeptr p, int swap)
{
  nodeptr       q, tmp;

  q = p->back;
  assert(!isTip(q->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));


  if(swap == 1)
   {
     tmp = p->next->back;
     hookup(p->next, q->next->back, q->next->z, tr->numBranches);
     hookup(q->next, tmp,           p->next->z, tr->numBranches);
   }
  else
   {
      tmp = p->next->next->back;
      hookup(p->next->next, q->next->back, q->next->z,       tr->numBranches);
      hookup(q->next,       tmp,           p->next->next->z, tr->numBranches);
   }
}
