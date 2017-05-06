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

#ifndef NGS_H
#define NGS_H

#include "tree/phylotree.h"
#include "alignment/alignmentpairwise.h"
#include "model/ratemeyerdiscrete.h"

/*
	collection of classes for Next-generation sequencing 
*/

/**
NGS Pairwise alignment

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class NGSAlignment : public AlignmentPairwise
{
public:

	/**
		constructor
		@param filename file in Fritz's format
	*/
	NGSAlignment(PhyloTree *atree);

    NGSAlignment(const char *filename);

	/**
		constructor
		@param nstate number of states
		@param ncat number of categories
		@param freq pair-state frequencies for all categories
	*/
	NGSAlignment(int nstate, int ncat, double *freq);

	NGSAlignment(int nstate, string &seq1, string &seq2);

	virtual char convertState(char state, SeqType seq_type);


	/**
		read file in Fritz's format
	*/
	void readFritzFile(const char *filename);

	/**
		compute empirical state frequencies from the alignment
		@param state_freq (OUT) is filled with state frequencies, assuming state_freq was allocated with 
			at least num_states entries.
	*/
	virtual void computeStateFreq(double *state_freq, size_t num_unknown_states = 0);

	/**
		compute the sum of pair state frequencies over all categories
		@param sum_pair_freq (OUT) will be filled in with num_states*num_states entries. 
			Memory has to be allocated before calling this function.
	*/
	void computeSumPairFreq (double *sum_pair_freq);

	/**
		compute empirical rates between state pairs
		@param rates (OUT) vector of size num_states*(num_states-1)/2 for the rates
	*/
	virtual void computeDivergenceMatrix (double *rates);

	/**
		compute the empirical distance for a category, used to initialize rate scaling factor
		@param cat specific category, between 0 and ncategory-1
	*/
	double computeEmpiricalDist(int cat);

	/**
		negative likelihood function for a category with a rate scaling factor
		@param cat specific category, between 0 and ncategory-1
		@param value a rate scaling factor
		@return negative log-likelihood (for minimization purpose)
	*/
	double computeFunctionCat(int cat, double value);

	/**
		negative likelihood and 1st and 2nd derivative function for a category with a rate scaling factor
		@param cat specific category, between 0 and ncategory-1
		@param value a rate scaling factor
		@param df (OUT) 1st derivative
		@param ddf (OUT) 2nd derivative
		@return negative log-likelihood (for minimization purpose)
	*/
	void computeFuncDervCat(int cat, double value, double &df, double &ddf);

	/**
		number of category
	*/
	int ncategory;

	//double *pair_freq;
};


class NGSTree : public PhyloTree {

public:

    /**
     * Constructor with given alignment
     * @param params program parameters
     * @param alignment
     */
	NGSTree(Params &params, NGSAlignment *alignment);	

    /**
            compute the tree likelihood
            @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihood(double *pattern_lh = NULL);

    /**
            optimize all branch lengths of the tree
            @param iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);

};

class NGSRate : public RateMeyerDiscrete {
public:

	/**
		@param tree must be NGSTree type
	*/
	NGSRate(PhyloTree *tree);

	/**
		get rate category of a specified site-pattern. 
		@param ptn pattern ID 
		@return the rate category of the specified site-pattern
	*/
	virtual int getPtnCat(int ptn) { return 0; }

	/**
		optimize rates of all site-patterns
		compute categorized rates from the "continuous" rate of the original Meyer & von Haeseler model.
		The current implementation uses the k-means algorithm with k-means++ package.
	*/
	virtual double optimizeParameters(double epsilon);


	/**
		This function is inherited from Optimization class for optimizting site rates 
		@param value x-value of the function
		@return f(value) of function f you want to minimize
	*/
	virtual double computeFunction(double value);

	/**
		This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
		@param value x-value of the function
		@param df (OUT) first derivative
		@param ddf (OUT) second derivative
		@return f(value) of function f you want to minimize
	*/
	virtual void computeFuncDerv(double value, double &df, double &ddf);

	/**
		classify rates into categories.
		@param tree_lh the current tree log-likelihood
	*/
	virtual double classifyRates(double tree_lh) { return tree_lh; }

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

};

class NGSRateCat : public RateMeyerDiscrete {
public:

	/**
		@param tree must be NGSTree type
	*/
	NGSRateCat(PhyloTree *tree, int ncat);

	/**
		optimize rates of all site-patterns
		compute categorized rates from the "continuous" rate of the original Meyer & von Haeseler model.
		The current implementation uses the k-means algorithm with k-means++ package.
	*/
	virtual double optimizeParameters(double epsilon);


	/**
		return the number of dimensions
	*/
	virtual int getNDim();
	

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

	/**
		proportion of position categories
	*/
	double *proportion;

	/**
		sum of pair freq from all positions
	*/
	double *sum_pair_freq;

protected:

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters 
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables);

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters 
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(double *variables);
};


class NGSTreeCat : public NGSTree {

public:

    /**
     * Constructor with given alignment
     * @param params program parameters
     * @param alignment
     */
	NGSTreeCat(Params &params, NGSAlignment *alignment);	
    /**
            compute the tree likelihood
            @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihood(double *pattern_lh = NULL);
};


class NGSRead : public NGSAlignment {
public:

	/** 
		constructor
	*/
	NGSRead(PhyloTree *atree);

	void init();

	//int orig_length;

	/**
		alignment score
	*/
	int score; //brauch ich das???

	/**
		read ID
	*/
	int id;

	//int scaff_id;

	/**
		matched position in the reference sequence
	*/
	int match_pos;

	/**
		TRUE for mapping forward strand, FALSE for backward
	*/
	bool direction;

	/**
		name of the reference sequence for which match found
	*/
	string chr;

	/**
		mapped portion of reference sequence
	*/
	string scaff;

	/**
		mapped portion of read
	*/
	string read;

	/**
		number of times that read is mapped (multiple optimal alignment scores)
	*/
	double times;

	/**
		TRUE if it is the first match, FALSE otherwise (in case of multiple hits)
	*/
	bool flag;

	/**
		alignment identity 
	*/
	float identity;

	/**
		read name
	*/
	string name;

	double homo_rate;

	void computePairFreq();

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
	
};

struct ReadInfo {
	int id;
	float identity;
	float distance;
	float logl;
	float homo_distance;
	float homo_logl;
};


class NGSReadSet : public vector<ReadInfo>  {
public:

	/**
		read in file containing mapped reads to the reference
		@param filename file name
		@param ref_ID reference sequence name to accept reads
		@param ident identity threshold to accept reads
		@param mismatches number of exact mismatches to accept reads
	*/
	void parseNextGen(string filename, string ref_ID="total",double ident=0.0,int mismatches=-1);

	/**
		this function will be called everytime a read is accepted from the parseNextGen()
		@param tempread read at current position while parsing 
	*/
	virtual void processReadWhileParsing(NGSRead &tempread);

	void writeFreqMatrix(ostream &out);

	/**
		write information
	*/
	void writeInfo();

	PhyloTree *tree;

	double homo_rate;

	vector<double*> pair_freq;

	vector<double*> state_freq;

	bool ngs_ignore_gaps;

};

/**
	Main function
	@param params input program parameters
*/
void runNGSAnalysis(Params &params);

#endif
