/*
 * candidateset.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#include "phylotree.h"
#include "candidateset.h"

void CandidateSet::init(Alignment* aln, Params *params) {
    this->aln = aln;
    this->params = params;
}

CandidateSet::~CandidateSet() {
}

CandidateSet::CandidateSet() {
	aln = NULL;
	params = NULL;
}

vector<string> CandidateSet::getBestTrees() {
	vector<string> res;
	double bestScore = rbegin()->first;
	for (reverse_iterator rit = rbegin(); rit != rend() && rit->second.score == bestScore; rit++) {
		res.push_back(rit->second.tree);
	}
	return res;
}

string CandidateSet::getRandCandTree() {
	assert(!empty());
	if (empty())
		return "";
	int id = random_int(min(params->popSize, (int)size()) );
	for (reverse_iterator i = rbegin(); i != rend(); i++, id--)
		if (id == 0)
			return i->second.tree;
	assert(0);
	return "";
}

vector<string> CandidateSet::getTopTrees(int numTree) {
	assert(numTree <= params->maxCandidates);
	if (numTree == 0) {
		numTree = params->maxCandidates;
	}
	vector<string> res;
	int cnt = numTree;
	for (reverse_iterator rit = rbegin(); rit != rend() && cnt > 0; rit++, cnt--) {
		res.push_back(rit->second.tree);
	}
	return res;
}

vector<string> CandidateSet::getBestLocalOptimalTrees(int numTree) {
	assert(numTree <= params->maxCandidates);
	if (numTree == 0) {
		numTree = params->maxCandidates;
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
        int count = params->popSize;
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

	if (treeTopologyExist(candidate.topology)) {
		newTree = false;
	    /* If tree topology already exist but the score is better, we replace the old one
	    by the new one (with new branch lengths) and update the score */
		if (topologies[candidate.topology] < score) {
			removeCandidateTree(candidate.topology);
			topologies[candidate.topology] = score;
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
		} else if (candidate.localOpt) {
			CandidateSet::iterator treePtr = getCandidateTree(candidate.topology);
			treePtr->second.localOpt = candidate.localOpt;
		}
	} else {
		if (getWorstScore() < score && size() >= params->maxCandidates) {
			// remove the worst-scoring tree
			topologies.erase(begin()->second.topology);
			erase(begin());
		}
		CandidateSet::iterator it = insert(CandidateSet::value_type(score, candidate));
		topologies[candidate.topology] = score;
		if (params->fix_stable_splits && getNumLocalOptTrees() >= params->numSupportTrees) {
			int it_pos = distance(it, end());
			// The new tree is one of the numSupportTrees best trees.
			// Thus recompute supported splits
			if (it_pos <= params->numSupportTrees) {
				int nSupportedSplits = computeSplitSupport(params->numSupportTrees);
				cout << ((double) nSupportedSplits / (aln->getNSeq() - 3)) * 100
						<< " % of the splits have 100% support and can be fixed." << endl;
			}
		}
	}
	assert(topologies.size() == size());
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

double CandidateSet::getBestScore() {
	if (size() == 0)
		return -DBL_MAX;
	else
		return rbegin()->first;
}

double CandidateSet::getWorstScore() {
	return begin()->first;
}

string CandidateSet::getTopology(string tree) {
	PhyloTree mtree;
//	mtree.rooted = params->is_rooted;
	mtree.aln = this->aln;
	mtree.setParams(params);

	stringstream str;
	str << tree;
	str.seekg(0, ios::beg);
//	freeNode();
	mtree.readTree(str, params->is_rooted);
	mtree.setAlignment(aln);
	mtree.setRootNode(params->root);

//	mtree.readTreeString(tree);
//	mtree.setRootNode(params->root);

	ostringstream ostr;
	mtree.printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	return ostr.str();
}

double CandidateSet::getTopologyScore(string topology) {
	assert(topologies.find(topology) != topologies.end());
	return topologies[topology];
}

void CandidateSet::clear() {
	multimap<double, CandidateTree>::clear();
	clearTopologies();
}

void CandidateSet::clearTopologies() {
	topologies.clear();
}


CandidateSet CandidateSet::getBestCandidateTrees(int numTrees) {
	CandidateSet res;
	if (numTrees >= size())
		numTrees = size();
	for (reverse_iterator rit = rbegin(); rit != rend() && numTrees > 0; rit++, numTrees--) {
		res.insert(*rit);
	}
	return res;
}

bool CandidateSet::treeTopologyExist(string topo) {
	return (topologies.find(topo) != topologies.end());
}

bool CandidateSet::treeExist(string tree) {
	return treeTopologyExist(getTopology(tree));
}

CandidateSet::iterator CandidateSet::getCandidateTree(string topology) {
	for (CandidateSet::reverse_iterator rit = rbegin(); rit != rend(); rit++) {
		if (rit->second.topology == topology)
			return --(rit.base());
	}
	return end();
}

void CandidateSet::removeCandidateTree(string topology) {
	bool removed = false;
	for (CandidateSet::reverse_iterator rit = rbegin(); rit != rend(); rit++) {
			if (rit->second.topology == topology) {
				erase( --(rit.base()) );
				topologies.erase(topology);
				removed = true;
				break;
			}
	}
	assert(removed);
}

bool CandidateSet::isStableSplit(Split& sp) {
	return stableSplit.containSplit(sp);
}

int CandidateSet::computeSplitSupport(int numTree) {
	stableSplit.clear();
	if (numTree == 0)
		numTree = getNumLocalOptTrees();
	SplitIntMap hash_ss;
	SplitGraph sg;
	MTreeSet boot_trees;
	int numMaxSupport = 0;
	vector<string> trees = getBestLocalOptimalTrees(numTree);
	assert(trees.size() > 1);
	int maxSupport = trees.size();
	boot_trees.init(trees, aln->getSeqNames(), params->is_rooted);
	boot_trees.convertSplits(aln->getSeqNames(), sg, hash_ss, SW_COUNT, -1, false);

	for (SplitIntMap::iterator it = hash_ss.begin(); it != hash_ss.end(); it++) {
		if (it->second == maxSupport && it->first->countTaxa() > 1) {
			numMaxSupport++;
			Split* supportedSplit = new Split(*(it->first));
			stableSplit.push_back(supportedSplit);
		}
	}
	//cout << "Number of supported splits = " << numMaxSupport << endl;
	return numMaxSupport;
}

void CandidateSet::setAln(Alignment* aln) {
	this->aln = aln;
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
