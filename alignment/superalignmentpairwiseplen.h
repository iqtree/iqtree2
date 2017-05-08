//
//  superalignmentpairwiseplen.h
//  iqtree
//
//  Created by Olga on 04/05/17.
//
//

#ifndef iqtree_superalignmentpairwiseplen_h
#define iqtree_superalignmentpairwiseplen_h

#include "tree/phylosupertreeplen.h"
#include "superalignmentpairwise.h"



class SuperAlignmentPairwisePlen : public SuperAlignmentPairwise {
    
public:
    
    /**
     constructor
     */
    
    SuperAlignmentPairwisePlen();
    
    /**
     construct the pairwise alignment from two sequences of a multiple alignment
     @param aln input multiple alignment
     @param seq_id1 ID of the first sequence
     @param seq_id2 ID of the second sequence
     */
    SuperAlignmentPairwisePlen(PhyloSuperTreePlen *atree, int seq1, int seq2);
    
    ~SuperAlignmentPairwisePlen();
    
    /**
     compute the likelihood for a distance between two sequences. Used for the ML optimization of the distance.
     @param value x-value of the function
     @return log-likelihood
     */
    virtual double computeFunction(double value);
    
    /**
     This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
     used by Newton raphson method to minimize the function.
     @param value x-value of the function
     @param df (OUT) first derivative
     @param ddf (OUT) second derivative
     @return f(value) of function f you want to minimize
     */
    virtual void computeFuncDerv(double value, double &df, double &ddf);
    
    /**
     partition information
     */
    vector<PartitionInfo>* part_info;
    
};

#endif
