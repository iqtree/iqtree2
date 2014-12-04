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
#include <stack>

struct CandidateTree {

	/**
	 * with branch lengths.
	 * empty for intermediate NNI tree
	 */
	string tree;


	/**
	 * tree topology WITHOUT branch lengths
	 * and WITH TAXON ID (instead of taxon names)
	 * for sorting purpose
	 */
	string topology;

	/**
	 * log-likelihood or parsimony score
	 */
	double score;

	/**
	 *  Indicate whether the tree is NNI locally optimal
	 */
	bool localOpt;

};


/**
 * Candidate tree set, sorted in ascending order of scores, i.e. the last element is the highest scoring tree
 */
class CandidateSet : public multimap<double, CandidateTree> {

public:
    /**
     * Initialization
     */
	void init(int maxCandidates, int maxPop, char* root, bool rooted, Alignment* aln);

	CandidateSet();

    /**
     * return randomly one candidate tree from max_candidate
     */
    string getRandCandTree();

    /**
     * return the next parent tree for reproduction.
     * Here we always maintain a list of candidate trees which have not
     * been used for reproduction. If all candidate trees have been used, we select the
     * current best trees as the new parent trees
     */
    string getNextCandTree();

    /**
     *  Replace an existing tree in the candidate set
     *  @param tree the new tree string that will replace the existing tree
     *  @param score the score of the new tree
     *  @return true if the topology of \a tree exist in the candidate set
     */
    bool replaceTree(string tree, double score);

    /**
     *  create the parent tree set containing top trees
     */
    void initParentTrees();

    /**
     * update / insert tree into set of score is higher than lowest-scoring tree
     * @param tree
     * 	The new tree string (with branch lengths)
     * @param score
     * 	The score (ML or parsimony) of the new tree
     * @param localOpt
     * 	Whether the tree is a locally optimal tree
     * @return false if the tree topology already exists
     *
     */
    bool update(string tree, double score, bool localOpt = true);

    /**
     *  print score of max_candidates best trees
     *
     *  @param numScore
     *  	Number of best scores to print out starting from the highest
     */
    vector<double> getBestScores(int numBestScores = 0);

    /**
     * Return the worst score in the candidate tree set
     * @return
     */
    double getWorstScore();

    /**
     * Return the best score in the candidate tree set
     * @return
     */
    double getBestScore();

    /**
     *  Return \a numTree best tree strings
     *  @param numTree number of best trees
     *  @return a list of tree
     */
    vector<string> getHighestScoringTrees(int numTree = 0);

    /**
     * Return \a numTree best local optimal trees
     * @param numTree
     * @return a vector of tree string
     */
    vector<string> getBestLocalOptimalTrees(int numTree = 0);

    /**
     * get tree(s) with highest score. More than one tree is
     * returned if there are multiple optima.
     * @return a vector containing optimal trees
     */
    vector<string> getEquallyOptimalTrees();

    /**
     * destructor
     */
    virtual ~CandidateSet();

    /**
     * check if tree topology WITHOUT branch length exist in the candidate set?
     */
    bool treeTopologyExist(string topo);

    /**
     * check if tree topology WITH branch length exist in the candidate set?
     */
    bool treeExist(string tree);

    /**
     * return a unique topology (sorted by taxon names, rooted at taxon with alphabetically smallest name) without branch lengths
     */
    string getTopology(string tree);

    /**
     * return the score of \a topology
     * @param topology
     * @return
     */
    double getTopologyScore(string topology);

    /**
     *  Empty the candidate set
     */
    void clear();

    /**
     *  Empty the \a topologies data structure;
     */
    void clearTopologies();

    /**
     * Compute the split support from the \a numTree best local optimal trees in the candidate sets
     * @param numTree the number of best trees used to calculate support values
     * @return number of splits with 100% support value
     */
    int computeSplitSupport(int numTree = 0);

    /**
     * Return a pointer to the \a CandidateTree that has topology equal to \a topology
     * @param topology
     * @return
     */
    CandidateSet::iterator getCandidateTree(string topology);

    /**
     * Remove the \a CandidateTree with topology equal to \a topology
     * @param topology
     */
    void removeCandidateTree(string topology);

    /* Getter and Setter function */
	void setAln(Alignment* aln);
	int getMaxCandidates() const;
	void setMaxCandidates(int maxCandidates);
	int getPopSize() const;
	void setPopSize(int popSize);
	void setIsRooted(bool isRooted);
	const StringDoubleHashMap& getTopologies() const {
		return topologies;
	}
	/* Getter and Setter function */

	/**
	 * get number of locally optimal trees in the set
	 * @return
	 */
	int getNumLocalOptTrees();

    /** Return a CandidateSet containing \a numTrees of current best candidate trees
     * @param numTrees
     * @return
     */
    CandidateSet getBestCandidateTrees(int numTrees);

private:
    /**
     * limit for number of trees (typically superset of candidate set)
     */
    int maxCandidates;

    /**
     *  maximum number of trees used for the evolving population
     */
    int popSize;

    /**
     *  Map data structure storing <topology_string, score>
     */
    StringDoubleHashMap topologies;

    /**
     *  Trees used for reproduction
     */
    stack<string> parentTrees;

    /**
     * pointer to alignment, just to assign correct IDs for taxa
     */
    Alignment *aln;


    /**
     *  Sequence name of the root
     */
    char* root;

	/**
	 *  Specify whether the tree is rooted or not
	 */
	bool isRooted;
};

#endif /* CANDIDATESET_H_ */
