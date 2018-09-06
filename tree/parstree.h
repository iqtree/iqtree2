/*
 * parstree.h
 *
 *  Created on: Nov 6, 2014
 *      Author: diep
 */

#ifndef PARSTREE_H_
#define PARSTREE_H_

#include "iqtree.h"

#define BootValTypePars unsigned short // Diep added

enum CostMatrixType {CM_UNIFORM, CM_LINEAR};

class ParsTree: public IQTree {
public:
    /**************************************************************************
     * Methods
     *************************************************************************/
    /**
     * default constructor
     */
    ParsTree();

    /**
     * Constructor with given alignment
     * @param alignment
     */
    ParsTree(Alignment *alignment);


    /**
     * destructor
     */
    ~ParsTree();

    /**
     * read the cost matrix file
     * initialize for 'nstates' and 'columns'
     */
    void loadCostMatrixFile(char* file_name = NULL);

    /**
     * initialise cost_matrix as linear
     * initialize for 'nstates' and 'columns'
     */
    void initCostMatrix(CostMatrixType cost_type);

//    /**
//     * allocate for ptn_pars if needed
//     */
//    void allocatePtnPars(int nptn);

    /**
        compute the tree parsimony score
        @return parsimony score of the tree
    */
    int computeParsimony();

    /**
        compute partial parsimony score of the subtree rooted at dad
        @param dad_branch the branch leading to the subtree
        @param dad its dad, used to direct the traversal
    */
    void computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad);


    /**
        compute tree parsimony score based on a particular branch
        @param dad_branch the branch leading to the subtree
        @param dad its dad, used to direct the traversal
        @param branch_subst (OUT) if not NULL, the number of substitutions on this branch
        @return parsimony score of the tree
    */
    int computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);

    /**
        initialize partial_pars vector of all PhyloNeighbors, allocating central_partial_pars
     */
    virtual void initializeAllPartialPars();

    /**
        initialize partial_pars vector of all PhyloNeighbors, allocating central_partial_pars
        @param node the current node
        @param dad dad of the node, used to direct the search
        @param index the index
     */
    virtual void initializeAllPartialPars(int &index, PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * allocate memory enough for one partial_pars vector of one neighbor
     */
    virtual UINT* newPartialPars();

    /**
     * calculate the size of one partial_pars vector of one neighbor
     */
    size_t getParsBlockSize();

	/**
	 * to overwrite the one in PhyloTree
	 */
    virtual UINT * newBitsBlock();

    /*
     * For a leaf character corresponding to an ambiguous state
     * set elements corresponding to possible states to 0, others to UINT_MAX
     */
    void initLeafSiteParsForAmbiguousState(char state, UINT * site_partial_pars);

    /*
     * For the starting phase: phyloanalysis.cpp
     */
    virtual bool isParsimonyTree() {
        return true;
    }

    void initParsData(Params* pars_params);

	void printPatternScore();
	UINT findMstScore(int ptn); // find minimum spanning tree score of a given pattern

    /**************************************************************************
     * Data
     *************************************************************************/
//    SankoffCostMatrix* cost_matrix;
    unsigned int * cost_matrix; // Sep 2016: store cost matrix in 1D array
    int cost_nstates; // Sep 2016: # of states provided by cost matrix
    UINT tree_pars;
    /**
     * Store array of pattern parsimony computed in computeParsimonyBranch()
     */
    //BootValTypePars* _pattern_pars;
};

#endif /* PARSTREE_H_ */
