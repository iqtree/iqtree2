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

//    /**
//     * allocate for ptn_pars if needed
//     */
//    void allocatePtnPars(int nptn);

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
    int cost_nstates; // Sep 2016: # of states provided by cost matrix
    UINT tree_pars;
    /**
     * Store array of pattern parsimony computed in computeParsimonyBranch()
     */
    //BootValTypePars* _pattern_pars;
};

#endif /* PARSTREE_H_ */
