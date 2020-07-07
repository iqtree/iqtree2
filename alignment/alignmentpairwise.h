/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef ALIGNMENTPAIRWISE_H
#define ALIGNMENTPAIRWISE_H

#include "utils/optimization.h"
#include "tree/phylotree.h"

/**
Pairwise alignment

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class AlignmentPairwise : public Alignment, public Optimization
{
public:
    AlignmentPairwise();

    /**
        pairwise alignment with sequence numbers not yet set
        @param atree input multiple alignment
     */
    AlignmentPairwise(PhyloTree *atree);

    /**
		construct the pairwise alignment from two sequences of a multiple alignment
		@param atree input multiple alignment
		@param seq1 ID of the first sequence
		@param seq2 ID of the second sequence
	*/
    AlignmentPairwise(PhyloTree *atree, int seq1, int seq2);

    
    /**
        recalculate the pairwise alignment for a different pair of sequences of
            the same multiple alignment it was constructed for
        @param seq1 ID of the first sequence
        @param seq2 ID of the second sequence
    */
    void setSequenceNumbers(int seq1, int seq2);
    
    
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
		compute the ML distance and variance between two sequences
		@param initial_dist initial guess
		@param (OUT) second derivative of likelihood function evaluated at ML distance
		@return the ML distance
	*/
	double optimizeDist(double initial_dist);

	double optimizeDist(double initial_dist, double &d2l);

	/**
		add a pattern into the alignment
		@param state1
		@param state2 states of the pattern
		@param freq frequency of pattern
		@param cat category for the pattern (for the discrete model)
		@return TRUE if pattern contains only gaps or unknown char. 
				In that case, the pattern won't be added.
	*/
	bool addPattern(int state1, int state2, int freq, int cat = 0);

    /**
        calculate the distance (or branch length) between two sequences
        @param seq1
        @param seq2 states of the pattern
        @param initial_dist   previous estimate of distance
        @param d2l
        @return a new estimate of branch length
    */

    virtual double recomputeDist( int seq1, int seq2, double initial_dist, double &d2l );
    
	/**
		destructor
	*/
    virtual ~AlignmentPairwise();

    size_t pairCount;
    size_t derivativeCalculationCount;
    size_t costCalculationCount;
    
protected:
	PhyloTree* tree;          //multi-species alignment tree from which sequences
                              //to be aligned are to be drawn
    int        num_states_squared; //the square of num_states
    int        total_size;    //number of elements in pair_freq
    double*    pair_freq;     //array of frequency counts (owned by this instance)
                              //size is num_states_squared times 1 (or by the number
                              //of categories).
    int        trans_size;    //number of elements (rows x columns) in transition matrices
    double*    trans_mat;     //used in computeFunction(),
    double*    sum_trans_mat; //used in computeFunction()
    double*    trans_derv1;   //used in computeFuncDerv()
    double*    trans_derv2;   //used in computeFuncDerv()
    double*    sum_derv1;     //used in computeFuncDerv()
    double*    sum_derv2;     //used in computeFuncDerv()
    double*    sum_trans;     //used in computeFuncDerv()

    int        seq_id1;
    int        seq_id2;
protected:
    void setTree(PhyloTree* atree);
    
    
};

#endif
