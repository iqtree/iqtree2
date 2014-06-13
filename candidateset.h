/*
 * candidateset.h
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#ifndef CANDIDATESET_H_
#define CANDIDATESET_H_
#include "tools.h"
#include "alignment.h"

struct CandidateTree {
	string topology; // tree topology without branch lengths for sorting purpose
	string tree; // with branch length
	double score; // log-likelihood under ML or parsimony score
};


/**
 * Candidate tree set
 */
class CandidateSet : public multimap<double, CandidateTree> {

public:
    /**
     * constructor
     */
	CandidateSet(int limit, int max_candidates, Alignment *aln);

	CandidateSet();

    /**
     * return tree with highest score
     */
    string getBestTree();

    /**
     * return randomly one candidate tree from max_candidate
     */
    string getCandidateTree();

    /**
     * update / insert tree into set of score is higher than lowest-scoring tree
     */
    bool update(string tree, double score);

    /**
     *  print score of max_candidates best trees
     */
    void printBestScores();

    /**
     * destroctor
     */
    virtual ~CandidateSet();

    /**
     * hard limit for number of trees (typically superset of candidate set)
     */
    int limit;

    /**
     *  maximum number of candidate trees
     */
    int max_candidates;

    /** index of tree topologies in set
     *
     */
    StringIntMap topologies;

    /**
     * pointer to alignment, just to assign correct IDs for taxa
     */
    Alignment *aln;

    bool treeTopologyExist(string topo);


    /**
     * return a unique topology (sorted by taxon names, rooted at taxon with alphabetically smallest name) without branch lengths
     */
    string getTopology(string tree);

};

#endif /* CANDIDATESET_H_ */
