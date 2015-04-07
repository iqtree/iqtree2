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
	loglThreshold = -DBL_MAX;
}

//void CandidateSet::getRandomStableSplits(int numSplit, SplitGraph& randomStableSplits) {
//	/*
//	 *  Use reservoir sampling technique
//	 */
//	randomStableSplits.clear();
//	assert(numSplit < candidateSplitsHash.size());
//	randomStableSplits.insert(randomStableSplits.begin(), candidateSplitsHash.begin(), candidateSplitsHash.begin() + numSplit);
//	for (int i = numSplit; i < candidateSplitsHash.size(); i++) {
//		int j = random_int(1, i);
//		if (j <= numSplit) {
//			randomStableSplits[j] = candidateSplitsHash[i];
//		}
//	}
//}

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

//vector<string> CandidateSet::getBestLocalOptimalTrees(int numTree) {
//	assert(numTree <= params->maxCandidates);
//	if (numTree == 0) {
//		numTree = params->maxCandidates;
//	}
//	vector<string> res;
//	int cnt = numTree;
//	for (reverse_iterator rit = rbegin(); rit != rend() && cnt > 0; rit++) {
//		if (rit->second.localOpt) {
//			res.push_back(rit->second.tree);
//			cnt--;
//		}
//	}
//	return res;
//}

/*
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
*/


void CandidateSet::addCandidateSplits(string treeString) {
	vector<string> taxaNames = aln->getSeqNames();
	MTree tree(treeString, taxaNames, params->is_rooted);
	SplitGraph allSplits;
	tree.convertSplits(allSplits);
	for (SplitGraph::iterator splitIt = allSplits.begin(); splitIt != allSplits.end(); splitIt++) {
		int value;
		Split *sp = candidateSplitsHash.findSplit(*splitIt, value);
		if (sp != NULL) {
			sp->setWeight(value + 1);
			candidateSplitsHash.setValue(sp, value + 1);
		} else {
			sp = new Split(*(*splitIt));
			sp->setWeight(1);
			candidateSplitsHash.insertSplit(sp, 1);
		}
	}
	candidateSplitsHash.setMaxValue(candidateSplitsHash.getMaxValue() + 1);
}

void CandidateSet::removeCandidateSplits(string treeString) {
	vector<string> taxaNames = aln->getSeqNames();
	MTree tree(treeString, taxaNames, params->is_rooted);
	SplitGraph allSplits;
	tree.convertSplits(allSplits);
	for (SplitGraph::iterator splitIt = allSplits.begin(); splitIt != allSplits.end(); splitIt++) {
		int value;
		Split *sp = candidateSplitsHash.findSplit(*splitIt, value);
		assert(sp->getWeight() == value);
		if (sp == NULL) {
			cout << "Cannot find split: ";
			(*splitIt)->report(cout);
			exit(1);
		} else {
			assert(sp->getWeight() >= 1);
			if (sp->getWeight() > 1) {
				sp->setWeight(value - 1);
			} else {
				candidateSplitsHash.eraseSplit(*splitIt);
			}
		}
	}
	candidateSplitsHash.setMaxValue(candidateSplitsHash.getMaxValue() - 1);
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


bool CandidateSet::update(string tree, double score) {
	bool newTree = true;
	CandidateTree candidate;
	candidate.score = score;
	candidate.topology = getTopology(tree);
	candidate.tree = tree;

	if (treeTopologyExist(candidate.topology)) {

		newTree = false;
	    /* If tree topology already exist but the score is better, we replace the old one
	    by the new one (with new branch lengths) and update the score */
		if (topologies[candidate.topology] < score) {
			//TODO Inefficient
			removeCandidateTree(candidate.topology);
			topologies[candidate.topology] = score;
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
		}

	} else {

		if (getWorstScore() < score && size() >= params->maxCandidates) {
			// remove the worst-scoring tree
			topologies.erase(begin()->second.topology);
			erase(begin());
		}

		CandidateSet::iterator candidateTreeIt = insert(CandidateSet::value_type(score, candidate));
		topologies[candidate.topology] = score;

		if (!candidateSplitsHash.empty()) {
			// ranking of the inserted tree
			int it_pos = distance(candidateTreeIt, end());

			// A new tree is inserted in the stable tree set
			if (it_pos <= params->numSupportTrees) {
				/*
				addCandidateSplits(candidateTreeIt->second.tree);
				if (candidateSplitsHash.getMaxValue() > params->numSupportTrees) {
					assert(candidateSplitsHash.getMaxValue() == params->numSupportTrees + 1);
					removeCandidateSplits(getNthBestTree(candidateSplitsHash.getMaxValue()).tree);
				}
				 */
				buildTopSplits();
				double percentSS = (double) getNumStableSplits() / (aln->getNSeq() - 3) * 100;
				cout << percentSS << " % of the splits have 100% support" << endl;
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


void CandidateSet::getBestCandidateTrees(int numTrees, vector<CandidateTree>& candidateTrees) {
	candidateTrees.clear();
	if (numTrees >= size())
		numTrees = size();
	for (reverse_iterator rit = rbegin(); rit != rend() && numTrees > 0; rit++, numTrees--) {
		candidateTrees.push_back(rit->second);
	}
}

CandidateTree CandidateSet::getNthBestTree(int N) {
	if (N >= size())
		N = size();
	reverse_iterator rit;
	CandidateTree myTree;
	myTree.score = -DBL_MAX;
	myTree.tree = "";
	myTree.topology = "";
	for (rit = rbegin(); rit != rend(); rit++, N--) {
		myTree = rit->second;
		if (N == 0)
			break;
	}
	return myTree;
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

// TODO This function is not efficient for large number of tree
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

int CandidateSet::buildTopSplits() {
	candidateSplitsHash.clear();
	vector<CandidateTree> bestCandidateTrees;

	getBestCandidateTrees(params->numSupportTrees, bestCandidateTrees);
	assert(bestCandidateTrees.size() > 1);

	candidateSplitsHash.setMaxValue(bestCandidateTrees.size());
	loglThreshold = bestCandidateTrees.back().score;

	/* Store all splits in the best trees in candidateSplitsHash.
	 * Note that the weight of each split is set equal to the number of trees that have this split.
	 * The variable maxWeight in SpitInMap is the number of trees, from which the splits are converted.
	 * This is also the maximum weight each split can have.
	 */
	vector<CandidateTree>::iterator treeIt;
	vector<string> taxaNames = aln->getSeqNames();
	for (treeIt = bestCandidateTrees.begin(); treeIt != bestCandidateTrees.end(); treeIt++) {
		MTree tree(treeIt->tree, taxaNames, params->is_rooted);
		SplitGraph splits;
		tree.convertSplits(splits);
		SplitGraph::iterator itg;
		for (itg = splits.begin(); itg != splits.end(); itg++) {
			int value;
			Split *sp = candidateSplitsHash.findSplit(*itg, value);
			if (sp != NULL) {
				sp->setWeight(value + 1);
				candidateSplitsHash.setValue(sp, value + 1);
			}
			else {
				sp = new Split(*(*itg));
				sp->setWeight(1);
				candidateSplitsHash.insertSplit(sp, 1);
			}
		}
	}

	numStableSplits = countStableSplits();
	return numStableSplits;
}

int CandidateSet::countStableSplits() {
	if (candidateSplitsHash.empty())
		return 0;
	int numMaxSupport = 0;
	for (SplitIntMap::iterator it = candidateSplitsHash.begin(); it != candidateSplitsHash.end(); it++) {
		if (it->second == candidateSplitsHash.getMaxValue() && it->first->countTaxa() > 1) {
			assert(it->first->getWeight() == candidateSplitsHash.getMaxValue());
			numMaxSupport++;
		}
	}
	return numMaxSupport;
}

void CandidateSet::setAln(Alignment* aln) {
	this->aln = aln;
}

//int CandidateSet::getNumLocalOptTrees() {
//	int numLocalOptima = 0;
//	for (reverse_iterator rit = rbegin(); rit != rend(); rit++) {
//		if (rit->second.localOpt) {
//			numLocalOptima++;
//		}
//	}
//	return numLocalOptima;
//}
