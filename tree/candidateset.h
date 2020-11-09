/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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

#ifndef CANDIDATESET_H_
#define CANDIDATESET_H_
//#include "phylotree.h"
#include "utils/tools.h"
#include "alignment/alignment.h"
#include "tree/mtreeset.h"
#include <stack>
#include "utils/checkpoint.h"


class IQTree;

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
};


/**
 * Candidate tree set, sorted in ascending order of scores, i.e. the last element is the highest scoring tree
 */
class CandidateSet : public multimap<double, CandidateTree>, public CheckpointFactory {

public:

    /**
     * Initialization
     */
	void init(Alignment* aln, int maxSize);

    CandidateSet();

	CandidateSet(int maxSize);

    /**
     *  Replace the current candidate trees by those in another candidate set
     *  @param candSet the candidate set whose trees will be copied over
     */
    void initTrees(CandidateSet& candSet);

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

    /**
     * return randomly one of the current best trees
     * @param numTopTrees [IN] Number of current best trees, from which a random tree is chosen.
     */
    string getRandTopTree(int numTopTrees);

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
//    bool replaceTree(string tree, double score);

    /**
     *  create the parent tree set containing top trees
     */
    void initParentTrees();

    /**
     *  update/insert \a tree into the candidate set if its score is higher than the worst tree
     *
     *  @param tree
     * 	    The new tree string (with branch lengths)
     *  @param score
     * 	    The score (ML or parsimony) of \a tree
     *  @return
     *      Relative position of the new tree to the current best tree.
     *      Return -1 if the tree topology already existed
     *      Return -2 if the candidate set is not updated
     */
    int update(string newTree, double newScore);

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
    vector<string> getBestTreeStrings(int numTree = 0);

    /**
     *  Get \a numTree top scoring trees for this MPI process. Also work for sequential version.
     *
     *  @param numTree
     *  	Number of top scoring trees
     *  @return
     *  	Vector of current best trees
     */
    vector<string> getBestTreeStringsForProcess(int totalNumTree);

    /**
     *  Return a set of trees and a set of scores
     *
     *  @param trees vector of trees
     *  @param scores vector of tree scores
     *  @param treeFormat the NEWICK format used for tree string (WT_TAXON_ID, WT_BR_LEN, ..)
     */
    void getAllTrees(vector<string> &trees, vector<double> &scores, int treeFormat = -1);

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
     * 	@param convertOption
     * 	    Use the same options as printTree() (WT_ID, WT_BR_LEN, ...)
     * 	@return
     * 		Newick string of the tree topology
     */
    string convertTreeString(const string tree, int format = WT_TAXON_ID | WT_SORT_TAXA);

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
     *  Collect all splits from the set of current best trees and compute for each of them the number of occurances.
     *
     *  @param supportThres
     *      a number in (0,1] representing the support value threshold for stable splits
     *  @return number of splits with 100% support value
     */
    int computeSplitOccurences(double supportThres);

   /**
    *   Get number of stable splits
    *   @param thresHold A number between (0,1.0], all splits have support values above this threshold
    *   are considered stable
    */
    int countStableSplits(double thresHold);

    void reportStableSplits();

    /**
     *  Update the set of stable split when a new tree is inserted
     *  to the set of best trees used for computing stable splits.
     *
     *  This function will remove all splits that belong to oldTree and add all
     *  splits of newTree
     *
     *  @param
     *  	oldTree tree that will be replace by \a newTree
     *  @param
     *  	newTree the new tree
     */
    void updateStableSplit(string oldTree, string newTree);

    /**
     * Return a pointer to the \a CandidateTree that has topology equal to \a topology
     * @param topology
     * @return
     */
    iterator getCandidateTree(string topology);

    /**
     * Remove candidate trees with topology equal to the specified topology
     * @param topology
     */
    void removeCandidateTree(string topology);

    /**
     *  Remove the worst tree in the candidate set
     */
    void removeWorstTree();

    /* Getter and Setter function */
	void setAln(Alignment* aln);

	const StringDoubleHashMap& getTopologies() const {
		return topologies;
	}

    /**
     * Return a CandidateSet containing \a numTrees candidate trees
     * @param numTrees
     * @return
     */
    CandidateSet getBestCandidateTrees(int numTrees = 0);

    /**
     *  Return a set of trees whose score are equal \a score
     */
    CandidateSet getCandidateTrees(double score);


	SplitIntMap& getCandSplits() {
		return candSplits;
	}

	/**
	 * @brief Get a random subset containing \a numSplit from the
	 * set of stable splits.
	 * @param
	 * 		numSplit size of the subset
	 * @param
	 * 		splits (OUT) a random subset of the stable splits
	 */
	//void getRandomStableSplits(int numSplit, SplitGraph& splits);

	/**
	 *  Add splits from \a treeString to the current candidate splits
	 *
	 *  @param tree collect splits from this tree
	 */
	void addCandidateSplits(string treeString);

	/**
	 *  Remove splits that appear from \a treeString.
	 *  If an existing split has weight > 1, their weight will be
	 *  reduced by 1.
	 */
	void removeCandidateSplits(string treeString);

    int getNumStableSplits() const {
        return numStableSplits;
    }

    /**
     *  Print candidate trees and their likelihood
     */
    void printTrees(string suffix);

    /**
     *  Recompute the log-likelihood of all trees
     *  @param treeObject the tree object which store other model parameters used
     *  to compute the log-likelihood.
     */
    void recomputeLoglOfAllTrees(IQTree &treeObject);

    int getMaxSize() const {
        return maxSize;
    }

    void setMaxSize(int max_size) {
        this->maxSize = max_size;
    }

private:
    /**
     *  Maximum number of candidate trees
     */
    int maxSize;

    /**
     *  Number of stable splits identified
     */
    int numStableSplits;

    /**
     *  Set of splits and the number of their occurences from the current best trees.
     *  The number of current best tree is parameterized.
     */
	SplitIntMap candSplits;

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
