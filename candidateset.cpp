/*
 * candidateset.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#include "phylotree.h"
#include "candidateset.h"
#include "mtreeset.h"


CandidateSet::CandidateSet(int maxCandidates, int maxPop, Alignment *aln) {
    assert(maxPop <= maxCandidates);
    assert(aln);
    this->maxCandidates = maxCandidates;
    this->popSize = maxPop;
    this->aln = aln;
    this->bestScore = -DBL_MAX;
    this->isRooted = false;
}

CandidateSet::CandidateSet() {
	aln = NULL;
	maxCandidates = 0;
	popSize = 0;
	bestScore = -DBL_MAX;
	isRooted = false;
}

vector<string> CandidateSet::getBestTree() {
	vector<string> res;
	for (reverse_iterator rit = rbegin(); rit != rend() && rit->second.score == bestScore; rit++) {
		res.push_back(rit->second.tree);
	}
	return res;
}

string CandidateSet::getRandCandTree() {
	if (empty())
		return "";
	// BQM: bug fix max -> min
	int id = random_int(min(popSize, (int)size()) );
	//int id = randint(0, min(max_candidates, (int)size()));
	//int id = 0 + (rand() % (int)( min(max_candidates, (int)size())) );
	for (reverse_iterator i = rbegin(); i != rend(); i++, id--)
		if (id == 0)
			return i->second.tree;
	return "";
}

vector<string> CandidateSet::getBestTrees(int numTree) {
	if (numTree == 0 || numTree > maxCandidates) {
		numTree = maxCandidates;
	}
	vector<string> res;
	int cnt = numTree;
	for (reverse_iterator rit = rbegin(); rit != rend() && cnt > 0; rit++, cnt--) {
		res.push_back(rit->second.tree);
	}
	return res;
}

vector<string> CandidateSet::getBestLocalOptimalTrees(int numTree) {
	assert(numTree <= maxCandidates);
	if (numTree == 0) {
		numTree = maxCandidates;
	}
	vector<string> res;
	int cnt = numTree;
	for (reverse_iterator rit = rbegin(); rit != rend() && cnt > 0; rit++) {
		if (rit->second.localOpt) {
			res.push_back(rit->second.tree);
			cnt--;
		}
	}
	return res;
}

bool CandidateSet::replaceTree(string tree, double score) {
    CandidateTree candidate;
    candidate.tree = tree;
    candidate.score = score;
    candidate.topology = getTopology(tree);
    if (treeTopologyExist(candidate.topology)) {
        topologies[candidate.topology] = score;
        for (reverse_iterator i = rbegin(); i != rend(); i++) {
            if (i->second.topology == candidate.topology) {
                erase( --(i.base()) );
                break;
            }
            insert(CandidateSet::value_type(score, candidate));
        }
    } else {
        return false;
    }
    return true;
}

string CandidateSet::getNextCandTree() {
    string tree;
    assert(!empty());
    if (parentTrees.empty()) {
        initParentTrees();
    }
    tree = parentTrees.top();
    parentTrees.pop();
    return tree;
}

void CandidateSet::initParentTrees() {
    if (parentTrees.empty()) {
        int count = this->popSize;
        for (reverse_iterator i = rbegin(); i != rend() && count >0 ; i++, count--) {
            parentTrees.push(i->second.tree);
            //cout << i->first << endl;
        }
    }
}

bool CandidateSet::update(string tree, double score, bool localOpt) {
	bool newTree = true;
	CandidateTree candidate;
	candidate.score = score;
	candidate.topology = getTopology(tree);
	candidate.localOpt = localOpt;
	candidate.tree = tree;
	if (candidate.score > bestScore)
		bestScore = candidate.score;
	if (treeTopologyExist(candidate.topology)) {
	    /* If tree topology already exist but the score is better, we replace the old one
	    by the new one (with new branch lengths) and update the score */
		if (topologies[candidate.topology] <= score) {
			topologies[candidate.topology] = score;
			for (CandidateSet::iterator i = begin(); i != end(); i++)
				if (i->second.topology == candidate.topology) {
					if (i->second.localOpt) {
						candidate.localOpt = i->second.localOpt;
					}
					erase(i);
					break;
				}
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
		}
		newTree = false;
	} else {
		if (size() < maxCandidates) {
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
			topologies[candidate.topology] = score;
		} else if (begin()->first < score){
			// remove the worst-scoring tree
			topologies.erase(begin()->second.topology);
			erase(begin());
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
			topologies[candidate.topology] = score;
		}
	}
	return newTree;
}

vector<double> CandidateSet::getBestScores(int numBestScore) {
	if (numBestScore == 0)
		numBestScore = size();
	vector<double> res;
	for (reverse_iterator rit = rbegin(); rit != rend() && numBestScore > 0; rit++, numBestScore--) {
		res.push_back(rit->first);
	}
	return res;
}

double CandidateSet::getWorstScore() {
	return begin()->first;
}

string CandidateSet::getTopology(string tree) {
	PhyloTree mtree;
	mtree.rooted = false;
	mtree.aln = aln;
	mtree.readTreeString(tree);
    mtree.root = mtree.findNodeName(aln->getSeqName(0));
	ostringstream ostr;
	mtree.printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	return ostr.str();
}

void CandidateSet::clear() {
	multimap<double, CandidateTree>::clear();
	topologies.clear();
}

CandidateSet::~CandidateSet() {
}

bool CandidateSet::treeTopologyExist(string topo) {
	return topologies.find(topo) != topologies.end();
}

bool CandidateSet::treeExist(string tree) {
	return treeTopologyExist(getTopology(tree));
}

int CandidateSet::computeSplitSupport(int numTree) {
	MTreeSet boot_trees;
	int numMaxSupport = 0;
	vector<string> trees = getBestLocalOptimalTrees(numTree);
	int maxSupport = trees.size();
	boot_trees.init(trees, isRooted);
	SplitGraph sg;
	SplitIntMap hash_ss;
	boot_trees.convertSplits(sg, hash_ss, SW_COUNT, -1);
	cout << sg.size() << " splits found" << endl;
	for (unordered_map<Split*,int>::iterator it = hash_ss.begin(); it != hash_ss.end(); it++) {
		if (it->second == maxSupport)
			numMaxSupport++;
	}
	cout << "Number of supported splits = " << numMaxSupport << endl;
	return numMaxSupport;
}

void CandidateSet::setAln(Alignment* aln) {
	this->aln = aln;
}

int CandidateSet::getMaxCandidates() const {
	return maxCandidates;
}

void CandidateSet::setMaxCandidates(int maxCandidates) {
	this->maxCandidates = maxCandidates;
}

int CandidateSet::getPopSize() const {
	return popSize;
}

void CandidateSet::setPopSize(int popSize) {
	this->popSize = popSize;
}

void CandidateSet::setIsRooted(bool isRooted) {
	this->isRooted = isRooted;
}

int CandidateSet::getNumLocalOptTrees() {
	int numLocalOptima = 0;
	for (reverse_iterator rit = rbegin(); rit != rend(); rit++) {
		if (rit->second.localOpt) {
			numLocalOptima++;
		}
	}
	return numLocalOptima;
}
