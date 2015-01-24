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
#include "mtreeset.h"
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
	 *  Indicate whether the tree is NNI locally optimal.
	 *  The reason to have this variable is that if the -reduction is
	 *  enabled, we will also store non-locally optimal trees in the set.
	 *  This is done to identify trees that belong to the same basin of attraction
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
	void init(Alignment* aln, Params *params);

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
     * update/insert \a tree into the candidate set if its score is higher than the worst tree
     *
     * @param tree
     * 	The new tree string (with branch lengths)
     * @param score
     * 	The score (ML or parsimony) of \a tree
     * @param localOpt
     * 	Tells whether \a tree is a locally optimal (DEFAULT: true)
     * @return false if tree topology already exists
     *
     */
    bool update(string tree, double score, bool localOpt = true);

    /**
     *  Get the \a numBestScores best scores in the candidate set
     *
     *  @param numBestScores
     *  	Number of best scores
     *  @return
     *  	Vector containing \a numBestScore best scores
     */
    vector<double> getBestScores(int numBestScores = 0);

    /**
     * Get the worst score
     *
     * @return the worst score
     */
    double getWorstScore();

    /**
     * Get best score
     *
     * @return the best score
     */
    double getBestScore();

    /**
     *  Get \a numTree top scoring trees
     *
     *  @param numTree
     *  	Number of top scoring trees
     *  @return
     *  	Vector of current best trees
     */
    vector<string> getTopTrees(int numTree = 0);

    /**
     * 	Get \a numTree best locally optimal trees
     * 	@param numTree
     * 		Number of locally optimal trees
     * 	@return
     * 		Vector of current best locally optimal trees
     */
    vector<string> getBestLocalOptimalTrees(int numTree = 0);

    /**
     * 	Get tree(s) with the best score. There could be more than one
     * 	tree that share the best score (this happens frequently with parsimony)
     * 	@return
     * 		A vector containing trees with the best score
     */
    vector<string> getBestTrees();

    /**
     * destructor
     */
    virtual ~CandidateSet();

    /**
     * 	Check if tree topology \a topo already exists
     *
     * 	@param topo
     * 		Newick string of the tree topology
     */
    bool treeTopologyExist(string topo);

    /**
     * 	Check if tree \a tree already exists
     *
     * 	@param tree
     * 		Newick string of the tree topology
     */
    bool treeExist(string tree);

    /**
     * 	Return a unique topology (sorted by taxon names, rooted at taxon with alphabetically smallest name)
     * 	without branch lengths
     *
     * 	@param tree
     * 		The newick tree string, from which the topology string will be generated
     * 	@return
     * 		Newick string of the tree topology
     */
    string getTopology(string tree);

    /**
     * return the score of \a topology
     *
     * @param topology
     * 		Newick topology
     * @return
     * 		Score of the topology
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
     * Check whether the
     * @param sp the split to check, must have the same taxon set as the trees in CandidateSet.
     * @return true if the \a supportedSplits contain \a sp, false otherwise.
     */
    bool isStableSplit(Split& sp);

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

	/**
	 * get number of locally optimal trees in the set
	 * @return
	 */
	int getNumLocalOptTrees();

    /**
     * Return a CandidateSet containing \a numTrees of current best candidate trees
     * @param numTrees
     * @return
     */
    CandidateSet getBestCandidateTrees(int numTrees);

	SplitGraph& getStableSplits() {
		return stableSplit;
	}

private:

    /**
     *  Set of supported splits by the best trees
     */
    SplitGraph stableSplit;

    /**
     *  Shared params pointing to the global params
     */
    Params* params;

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

};

#endif /* CANDIDATESET_H_ */
