#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"

void read_msa (tree * tr, char * filename);

/* This is the info you need to copy the vector*/
typedef struct
{
  int node_number;
  int num_partitions;
  size_t *partition_sizes;
  double **lh_values;
}likelihood_vector;

void free_likelihood_vector(likelihood_vector *v)
{
  if(v == NULL)
    return;
  int i;
  for(i=0; i < v->num_partitions; i++)
    free(v->lh_values[i]);
  free(v->lh_values);
  free(v->partition_sizes);
  free(v);
}

likelihood_vector *copy_likelihood_vectors (tree *tr, nodeptr p)
{
  assert(tr->useRecom == FALSE);
  likelihood_vector *v = (likelihood_vector *) malloc(sizeof(likelihood_vector));  
  v->node_number = p->number; 
  v->num_partitions = tr->NumberOfModels; 
  v->partition_sizes = (size_t *)malloc(tr->NumberOfModels * sizeof(size_t));
  v->lh_values = (double **)malloc(tr->NumberOfModels * sizeof(double *));

  /* Compute LH vector sizes for each partition */
  size_t rateHet, states, width, vector_size;
  rateHet = discreteRateCategories(tr->rateHetModel);
  int model;
  for(model = 0; model < tr->NumberOfModels; model++)
  {
    width  = (size_t)tr->partitionData[model].width;
    states = (size_t)tr->partitionData[model].states;
    vector_size = virtual_width( width ) * rateHet * states * sizeof(double);
    v->lh_values[model] = (double *)malloc(sizeof(double) * vector_size);
    assert (v->lh_values[model] != NULL);
    v->partition_sizes[model] = vector_size;
    double *lh_vector_src = tr->partitionData[model].xVector[p->number - tr->mxtips - 1];
    assert (lh_vector_src != NULL);
    vector_size = v->partition_sizes[model];
    memcpy(v->lh_values[model], lh_vector_src, vector_size);
  }
  return v;
}

void restore_vector(tree *tr, nodeptr p, likelihood_vector *v)
{
  int model;
  for(model = 0; model < tr->NumberOfModels; model++)
  {
    double *lh_vector_dest = tr->partitionData[model].xVector[p->number - tr->mxtips - 1];
    memcpy(lh_vector_dest, v->lh_values[model], v->partition_sizes[model]);
  }
}

boolean same_vector(tree *tr, nodeptr p, likelihood_vector *v)
{
  int i, model;
  for(model=0; model<tr->NumberOfModels; model++)
  {
    double *lh_vector_tree = tr->partitionData[model].xVector[p->number - tr->mxtips - 1];
    int len = (int)v->partition_sizes[model]/sizeof(double);
    for(i=0; i<len; i++)
    {
      if(v->lh_values[model][i] != lh_vector_tree[i])
      {
        printf("Diff entry in partition %d, site %d of %f\n", model, i, fabs(v->lh_values[model][i] - lh_vector_tree[i])); 
        return FALSE;
      }
    }
  }
  return TRUE;
}

int main(int argc, char * argv[])
{

  tree        * tr;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s [binary-alignment-file]\n", argv[0]);
     return (1);
   }
  tr = (tree *)malloc(sizeof(tree));

  /* read the binary input, setup tree, initialize model with alignment */
  read_msa(tr,argv[1]);
  tr->randomNumberSeed = 665;
  makeRandomTree(tr);
  printf("Number of taxa: %d\n", tr->mxtips);
  printf("Number of partitions: %d\n", tr->NumberOfModels);


  /* compute the LH of the full tree */
  printf ("Virtual root: %d\n", tr->start->number);
  evaluateGeneric(tr, tr->start, TRUE);
  printf("Likelihood: %f\n", tr->likelihood);

  /* 8 rounds of branch length optimization */
  smoothTree(tr, 1);
  evaluateGeneric(tr, tr->start, TRUE);
  printf("Likelihood after branch length optimization: %.20f\n", tr->likelihood);



  /* Now we show how to find a particular LH vector for a node */
  int i;
  int node_number = tr->mxtips + 1;
  nodeptr p = tr->nodep[node_number];
  printf("Pointing to  node %d\n", p->number);

  /* Fix as VR */
  newviewGeneric(tr, p, FALSE);
  newviewGeneric(tr, p->back, FALSE);
  evaluateGeneric(tr, p, FALSE);
  printf("Likelihood : %.f\n", tr->likelihood);

  printf("Make a copy of LH vector for node  %d\n", p->number);
  likelihood_vector *vector = copy_likelihood_vectors(tr, p);
  for(i=0; i<vector->num_partitions; i++)
     printf("Partition %d requires %d bytes\n", i, (int)vector->partition_sizes[i]);

  /* Check we have the same vector in both tree and copied one */
  assert(same_vector(tr, p, vector));

  /* Now force the p to get a new value (generally branch lengths are NOT updated like this) */
  /* This is just an example to show usage (for fast NNI eval), manually updating vectors is not recommended! */
  printf("bl : %.40f\n", p->next->z[0]);
  p->next->z[0] = p->next->back->z[0] = zmin;
  printf("bl : %.40f\n", p->next->z[0]);
  newviewGeneric(tr, p, FALSE);
  assert(!same_vector(tr, p, vector));
  evaluateGeneric(tr, p, FALSE);
  printf("Likelihood : %f\n", tr->likelihood);

  restore_vector(tr, p, vector);
  assert(same_vector(tr, p, vector));
  evaluateGeneric(tr, p, FALSE);
  printf("Likelihood after manually restoring the vector : %f\n", tr->likelihood);

  free_likelihood_vector(vector);

  /* Pick an inner branch */
  printf("numBranches %d \n", tr->numBranches);
  //tr->numBranches = 1;
  p = tr->nodep[tr->mxtips + 1];
  int partition_id = 0; /* single partition */
  double bl = get_branch_length(tr, p, partition_id);
  printf("z value: %f , bl value %f\n", p->z[partition_id], bl);
  /* set the bl to 2.5 */
  double new_bl = 2.5;
  set_branch_length(tr, p, partition_id, new_bl);
  printf("Changed BL to %f\n", new_bl);
  printf("new z value: %f , new bl value %f\n", p->z[partition_id], get_branch_length(tr, p, partition_id));
  /* set back to original */
  printf("Changed to previous BL\n");
  set_branch_length(tr, p, partition_id, bl);
  printf("new z value: %f , new bl value %f\n", p->z[partition_id], get_branch_length(tr, p, partition_id));

  return (0);
}
