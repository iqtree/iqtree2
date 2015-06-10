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
	maxSize = params->maxPopSize;
}

CandidateSet::~CandidateSet() {
}

CandidateSet::CandidateSet() {
	aln = NULL;
	params = NULL;
	loglThreshold = -DBL_MAX;
	numStableSplits = 0;
	maxSize = 1000;
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
	if (numTree == 0 || numTree > maxSize) {
		numTree = maxSize;
	}
	vector<string> res;
	int cnt = numTree;
	for (reverse_iterator rit = rbegin(); rit != rend() && cnt > 0; rit++, cnt--) {
		res.push_back(rit->second.tree);
	}
	return res;
}

//vector<string> CandidateSet::getBestLocalOptimalTrees(int numTree) {
//	assert(numTree <= params->maxPopSize);
//	if (numTree == 0) {
//		numTree = params->maxPopSize;
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
    candidate.topology = getTopologyString(tree);
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
	candidateSplitsHash.setNumTree(candidateSplitsHash.getNumTree() + 1);
}

void CandidateSet::removeCandidateSplits(string treeString) {
	vector<string> taxaNames = aln->getSeqNames();
	MTree tree(treeString, taxaNames, params->is_rooted);
	SplitGraph allSplits;
	tree.convertSplits(allSplits);
	for (SplitGraph::iterator splitIt = allSplits.begin(); splitIt != allSplits.end(); splitIt++) {
		int value = 0;
		Split *sp;
		sp = candidateSplitsHash.findSplit(*splitIt, value);
		if (value == 0) {
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
	candidateSplitsHash.setNumTree(candidateSplitsHash.getNumTree() - 1);
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


bool CandidateSet::update(string newTree, double newScore) {
	bool notExisted = true;
	CandidateTree candidate;
	candidate.score = newScore;
	candidate.topology = getTopology(newTree);
	candidate.tree = newTree;

	if (treeTopologyExist(candidate.topology)) {
		notExisted = false;
		// update new score if it is better the old score
		double oldScore = topologies[candidate.topology];
		if (oldScore < (newScore - 1e-6) || topologies[candidate.topology] == 0.0) {
			removeCandidateTree(candidate.topology, oldScore);
			topologies[candidate.topology] = newScore;
			// insert tree into candidate set
			insert(CandidateSet::value_type(newScore, candidate));
		}

	} else {

//		if (getWorstScore() < newScore && size() >= maxSize) {
//			// remove the worst-scoring tree
//			topologies.erase(begin()->second.topology);
//			erase(begin());
//		}

		CandidateSet::iterator candidateTreeIt = insert(CandidateSet::value_type(newScore, candidate));
		topologies[candidate.topology] = newScore;

		if (size() > maxSize) {
			// remove the worst-scoring tree
			topologies.erase(begin()->second.topology);
			erase(begin());
		}

		if (!candidateSplitsHash.empty()) {
			// ranking of the inserted tree
			int it_pos = distance(candidateTreeIt, end());

			// A new tree is inserted in the stable tree set
			if (it_pos <= params->numSupportTrees) {
//				addCandidateSplits(candidateTreeIt->second.tree);
//				if (candidateSplitsHash.getMaxValue() > params->numSupportTrees) {
//					assert(candidateSplitsHash.getMaxValue() == params->numSupportTrees + 1);
//					removeCandidateSplits(getNthBestTree(candidateSplitsHash.getMaxValue()).tree);
//				}
				buildTopSplits(params->stableSplitThreshold);
//				reportStableSplits();
			}
		}

	}

	assert(topologies.size() == size());
	return notExisted;
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

vector<string> CandidateSet::getBestTreeStrings(int numTrees) {
    if (numTrees == 0 || numTrees > size())
        numTrees = size();
    vector<string> bestTrees; 
	for (reverse_iterator rit = rbegin(); rit != rend() && numTrees > 0; 
            rit++, numTrees--) {
        bestTrees.push_back(rit->second.tree);
    }
    return bestTrees;
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
void CandidateSet::removeCandidateTree(string topology, double score) {
	bool removed = false;
	CandidateSet::iterator it;
	while ( (it = find(score)) != end()) {
		if (it->second.topology == topology) {
			erase(it);
			removed = true;
			break;
		}
	}
//	for (CandidateSet::reverse_iterator rit = rbegin(); rit != rend(); rit++) {
//			if (rit->second.topology == topology) {
//				erase( --(rit.base()) );
//				topologies.erase(topology);
//				removed = true;
//				break;
//			}
//	}
	assert(removed);
}

int CandidateSet::buildTopSplits(double supportThreshold) {
	candidateSplitsHash.clear();
	vector<CandidateTree> bestCandidateTrees;

	getBestCandidateTrees(params->numSupportTrees, bestCandidateTrees);
	//assert(bestCandidateTrees.size() > 1);

	candidateSplitsHash.setNumTree(bestCandidateTrees.size());
	loglThreshold = bestCandidateTrees.back().score;

	/* Store all splits in the best trees in candidateSplitsHash.
	 * The variable numTree in SpitInMap is the number of trees, from which the splits are converted.
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
				int newHashWeight = value + 1;
				double newSupport = (double) newHashWeight / (double) candidateSplitsHash.getNumTree();
				sp->setWeight(newSupport);
				candidateSplitsHash.setValue(sp, newHashWeight);
			}
			else {
				sp = new Split(*(*itg));
				sp->setWeight(1.0 / (double) candidateSplitsHash.getNumTree());
				candidateSplitsHash.insertSplit(sp, 1);
			}
		}
	}
	int newNumStableSplits = countStableSplits(params->stableSplitThreshold);
	if (newNumStableSplits > numStableSplits || verbose_mode >= VB_MED) {
		numStableSplits = newNumStableSplits;
		cout << ((double) numStableSplits / (aln->getNSeq() - 3)) * 100 ;
		cout << " % of the splits are stable (support threshold " << params->stableSplitThreshold;
		cout << " from " << candidateSplitsHash.getNumTree() << " trees)" << endl;
	}
	return numStableSplits;
}

int CandidateSet::countStableSplits(double thresHold) {
	if (thresHold > 1.0)
		thresHold = 1.0;
	if (candidateSplitsHash.empty())
		return 0;
	int numMaxSupport = 0;
	for (SplitIntMap::iterator it = candidateSplitsHash.begin(); it != candidateSplitsHash.end(); it++) {
		if (it->first->getWeight() > thresHold && it->first->countTaxa() > 1) {
			//cout << "Stable support: " << it->first->getWeight() << endl;
			numMaxSupport++;
		}
	}
	return numMaxSupport;
}

void CandidateSet::reportStableSplits() {
	if (candidateSplitsHash.empty()) {
		cout << "The set of stable splits is empty! " << endl;
		return;
	}

	int numMaxSupport = 0;
	for (SplitIntMap::iterator it = candidateSplitsHash.begin(); it != candidateSplitsHash.end(); it++) {
		if (it->second == candidateSplitsHash.getNumTree() && it->first->countTaxa() > 1) {
			cout << it->first->getWeight() << " / " << candidateSplitsHash.getNumTree() << endl;
			assert(it->first->getWeight() == candidateSplitsHash.getNumTree());
			it->first->report(cout);
		}
	}
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
