/*
 * candidateset.cpp
 *
 *  Created on: Jun 1, 2014
 *  Author: Tung Nguyen
 *  Email: nltung@gmail.com
 */

#include "tree/iqtree.h"
#include "candidateset.h"
#include "utils/MPIHelper.h"

void CandidateSet::init(Alignment *aln, int maxSize) {
    this->aln = aln;
    this->maxSize = maxSize;
}

CandidateSet::~CandidateSet() {
}

CandidateSet::CandidateSet() : CheckpointFactory() {
    aln = NULL;
    numStableSplits = 0;
    this->maxSize = Params::getInstance().maxCandidates;
}

void CandidateSet::initTrees(CandidateSet& candSet) {
    int curMaxSize = this->maxSize;
    *this = candSet;
    setMaxSize(curMaxSize);
}



void CandidateSet::saveCheckpoint() {
    checkpoint->startStruct("CandidateSet");
    int ntrees = min(Params::getInstance().numNNITrees, (int) size());
    checkpoint->startList(Params::getInstance().numNNITrees);
    for (reverse_iterator it = rbegin(); it != rend() && ntrees > 0; it++, ntrees--) {
        checkpoint->addListElement();
        stringstream ss;
        ss.precision(12);
        ss << it->second.score << " " << it->second.tree;
//        double score = it->second.score;
//        CKP_SAVE(score);
//        checkpoint->put("tree", it->second.tree);
        checkpoint->put("", ss.str());
    }
    checkpoint->endList();
    checkpoint->endStruct();
    CheckpointFactory::saveCheckpoint();
}

void CandidateSet::restoreCheckpoint() {
    CheckpointFactory::restoreCheckpoint();
    checkpoint->startStruct("CandidateSet");
    double score;
    string tree;
    checkpoint->startList(Params::getInstance().numNNITrees);
    for (int i = 0; i < Params::getInstance().numNNITrees; i++) {
        checkpoint->addListElement();
        string str;
        if (!checkpoint->getString("", str)) {
            break;
        }
        stringstream ss(str);
        ss >> score >> tree;
//        CKP_RESTORE(tree);
        update(tree, score);

    }
    checkpoint->endList();
    checkpoint->endStruct();
}


string CandidateSet::getRandTopTree(int numTopTrees) {
    ASSERT(!empty());
    if (empty())
        return "";
    int id = random_int(min(numTopTrees, (int) size()));
    for (reverse_iterator it = rbegin(); it != rend(); it++) {
        if (id == 0)
            return it->second.tree;
        id--;
    }
    ASSERT(0);
    return "";
}

vector<string> CandidateSet::getBestTreeStrings(int numTree) {
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

vector<string> CandidateSet::getBestTreeStringsForProcess(int numTree) {
    int numProc = MPIHelper::getInstance().getNumProcesses();
    int procID = MPIHelper::getInstance().getProcessID();

    if (numTree < numProc)
        numTree = numProc; // BUG FIX: make sure that each process gets at least 1 tree

    vector<string> alltrees = getBestTreeStrings(numTree);
    if (numProc == 1) return alltrees;
    
    if (numTree == 0 || numTree > alltrees.size()) {
        numTree = alltrees.size();
    }
    int cnt = 0;
    vector<string> res;
    // process will get trees indexed procID, procID+1*numProc, procID+2*numProc,...
    for (cnt = procID; cnt < numTree; cnt+=numProc) {
        res.push_back(alltrees[cnt]);
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
    MTree tree(treeString, taxaNames, Params::getInstance().is_rooted);
    SplitGraph allSplits;
    tree.convertSplits(allSplits);
    for (SplitGraph::iterator splitIt = allSplits.begin(); splitIt != allSplits.end(); splitIt++) {
        int value;
        Split *sp = candSplits.findSplit(*splitIt, value);
        if (sp != NULL) {
            sp->setWeight(value + 1);
            candSplits.setValue(sp, value + 1);
        } else {
            sp = new Split(*(*splitIt));
            sp->setWeight(1);
            candSplits.insertSplit(sp, 1);
        }
    }
    candSplits.setNumTree(candSplits.getNumTree() + 1);
}

void CandidateSet::removeCandidateSplits(string treeString) {
    vector<string> taxaNames = aln->getSeqNames();
    MTree tree(treeString, taxaNames, Params::getInstance().is_rooted);
    SplitGraph allSplits;
    tree.convertSplits(allSplits);
    for (SplitGraph::iterator splitIt = allSplits.begin(); splitIt != allSplits.end(); splitIt++) {
        int value = 0;
        Split *sp;
        sp = candSplits.findSplit(*splitIt, value);
        if (value == 0) {
            cout << "Cannot find split: ";
            (*splitIt)->report(cout);
            exit(1);
        } else {
            ASSERT(sp->getWeight() >= 1);
            if (sp->getWeight() > 1) {
                sp->setWeight(value - 1);
            } else {
                candSplits.eraseSplit(*splitIt);
            }
        }
    }
    candSplits.setNumTree(candSplits.getNumTree() - 1);
}

string CandidateSet::getNextCandTree() {
    string tree;
    ASSERT(!empty());
    if (parentTrees.empty()) {
        initParentTrees();
    }
    tree = parentTrees.top();
    parentTrees.pop();
    return tree;
}

void CandidateSet::initParentTrees() {
    if (parentTrees.empty()) {
        int count = Params::getInstance().popSize;
        for (reverse_iterator i = rbegin(); i != rend() && count > 0; i++, count--) {
            parentTrees.push(i->second.tree);
            //cout << i->first << endl;
        }
    }
}


int CandidateSet::update(string newTree, double newScore) {
    // Do not update candidate set if the new tree has worse score than the
    // worst tree in the candidate set
    auto front = begin();
    if ( size() >= maxSize && front!=end() && newScore < front->first ) {
        return -2;
    }
    CandidateTree candidate;
    candidate.score = newScore;
    candidate.topology = convertTreeString(newTree);
    candidate.tree = newTree;

    int treePos;
    CandidateSet::iterator candidateTreeIt;

    if (treeTopologyExist(candidate.topology)) {
        // update new score if it is better the old score
        double oldScore = topologies[candidate.topology];
        if (oldScore < newScore) {
            removeCandidateTree(candidate.topology);
            insert(CandidateSet::value_type(newScore, candidate));
            topologies[candidate.topology] = newScore;
        }
        ASSERT(topologies.size() == size());
        return -1;
    }

    candidateTreeIt = insert(CandidateSet::value_type(newScore, candidate));
    topologies[candidate.topology] = newScore;

    if (size() > maxSize) {
        removeWorstTree();
    }
    ASSERT(topologies.size() == size());

    treePos = distance(candidateTreeIt, end());

    return treePos;
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

string CandidateSet::convertTreeString(string treeString, int format) {
    MTree mtree;
    stringstream str;
    str << treeString;
    str.seekg(0, ios::beg);
    mtree.readTree(str, Params::getInstance().is_rooted);
    mtree.assignLeafID();
    string rootName = "0";
    mtree.root = mtree.findLeafName(rootName);

    ostringstream ostr;
    mtree.printTree(ostr, format);
    return ostr.str();
}

string CandidateSet::getTopology(string tree) {
//	PhyloTree mtree;
//	mtree.rooted = params->is_rooted;
//	mtree.aln = this->aln;
//	mtree.setParams(params);
    MTree mtree;

    stringstream str;
    str << tree;
    str.seekg(0, ios::beg);
//	freeNode();
    mtree.readTree(str, Params::getInstance().is_rooted);
//	mtree.setAlignment(aln);
//	mtree.setRootNode(params->root);
    mtree.assignLeafID();
    string x = "0";
    mtree.root = mtree.findLeafName(x);

    ostringstream ostr;
    mtree.printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
    return ostr.str();
}

double CandidateSet::getTopologyScore(string topology) {
    ASSERT(topologies.find(topology) != topologies.end());
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
    if (numTrees >= size() || numTrees == 0)
        numTrees = (int) size();

    for (reverse_iterator rit = rbegin(); rit != rend() && numTrees > 0; rit++, numTrees--) {
        res.insert(*rit);
    }
    return res;
}

void CandidateSet::getAllTrees(vector<string> &trees, vector<double> &scores, int format) {
    trees.clear();
    scores.clear();

    for (reverse_iterator rit = rbegin(); rit != rend(); rit++) {
        if (format != -1) {
            trees.push_back(convertTreeString(rit->second.tree, format));
        } else {
            trees.push_back(rit->second.tree);
        }
        scores.push_back(rit->first);
    }
}

bool CandidateSet::treeTopologyExist(string topo) {
    return (topologies.find(topo) != topologies.end());
}

bool CandidateSet::treeExist(string tree) {
    return treeTopologyExist(convertTreeString(tree));
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
    double treeScore;
    // Find the score of the topology
    treeScore = topologies[topology];
    // Remove the topology
    topologies.erase(topology);
    pair<CandidateSet::iterator, CandidateSet::iterator> treeItPair;
    // Find all trees with that score
    treeItPair = equal_range(treeScore);
    CandidateSet::iterator it;
    for (it = treeItPair.first; it != treeItPair.second; ++it) {
        if (it->second.topology == topology) {
            erase(it);
            removed = true;
            break;
        }
    }
    ASSERT(removed);
}


void CandidateSet::removeWorstTree() {
    topologies.erase(begin()->second.topology);
    erase(begin());
}

int CandidateSet::computeSplitOccurences(double supportThreshold) {
    candSplits.clear();
    candSplits.setNumTree(size());

    /* Store all splits in the best trees in candSplits.
     * The variable numTree in SpitInMap is the number of trees, from which the splits are converted.
     */
    CandidateSet::iterator treeIt;
    //vector<string> taxaNames = aln->getSeqNames();
    for (treeIt = begin(); treeIt != end(); treeIt++) {
        MTree tree(treeIt->second.tree, Params::getInstance().is_rooted);
        SplitGraph splits;
        tree.convertSplits(splits);
        SplitGraph::iterator itg;
        for (itg = splits.begin(); itg != splits.end(); itg++) {
            int value;
            Split *sp = candSplits.findSplit(*itg, value);
            if (sp != NULL) {
                int newHashWeight = value + 1;
                double newSupport = (double) newHashWeight / (double) candSplits.getNumTree();
                sp->setWeight(newSupport);
                candSplits.setValue(sp, newHashWeight);
            }
            else {
                sp = new Split(*(*itg));
                sp->setWeight(1.0 / (double) candSplits.getNumTree());
                candSplits.insertSplit(sp, 1);
            }
        }
    }
    int newNumStableSplits = countStableSplits(supportThreshold);
    if (verbose_mode >= VB_MED) {
        cout << ((double) newNumStableSplits / (aln->getNSeq() - 3)) * 100;
        cout << " % of the splits are stable (support threshold " << supportThreshold;
        cout << " from " << candSplits.getNumTree() << " trees)" << endl;
    }

    return numStableSplits;
}

int CandidateSet::countStableSplits(double thresHold) {
    if (thresHold >= 1.0)
        thresHold = 0.99;
    if (candSplits.empty())
        return 0;
    int numMaxSupport = 0;
    for (SplitIntMap::iterator it = candSplits.begin(); it != candSplits.end(); it++) {
        if (it->first->getWeight() >= thresHold && it->first->countTaxa() > 1) {
            //cout << "Stable support: " << it->first->getWeight() << endl;
            numMaxSupport++;
        }
    }
    return numMaxSupport;
}

void CandidateSet::reportStableSplits() {
    if (candSplits.empty()) {
        cout << "The set of stable splits is empty! " << endl;
        return;
    }

//    int numMaxSupport = 0;
    for (SplitIntMap::iterator it = candSplits.begin(); it != candSplits.end(); it++) {
        if (it->second == candSplits.getNumTree() && it->first->countTaxa() > 1) {
            cout << it->first->getWeight() << " / " << candSplits.getNumTree() << endl;
            ASSERT(it->first->getWeight() == candSplits.getNumTree());
            it->first->report(cout);
        }
    }
}

void CandidateSet::setAln(Alignment *aln) {
    this->aln = aln;
}

CandidateSet CandidateSet::getCandidateTrees(double score) {
    CandidateSet res;
    for (CandidateSet::iterator it = begin(); it != end(); it++) {
        if (abs(it->first - score) < 0.1) {
            res.insert(*it);
        }
    }
    return res;
}

void CandidateSet::printTrees(string suffix) {
    ofstream outTrees, outLHs;
    string outTreesFile = string(Params::getInstance().out_prefix) + "." + suffix;
    string outLHsFile = string(Params::getInstance().out_prefix) + "." + suffix + "_lh";
    outTrees.open(outTreesFile.c_str());
    outLHs.open(outLHsFile.c_str());
    outLHs.precision(15);
    for (reverse_iterator rit = rbegin(); rit != rend(); rit++) {
        outLHs << rit->first << endl;
        outTrees << rit->second.topology << endl;
    }
    outTrees.close();
    outLHs.close();
}

void CandidateSet::recomputeLoglOfAllTrees(IQTree &treeObject) {
    vector<string> allTreeStrings = getBestTreeStrings();
    for (vector<string>:: iterator it = allTreeStrings.begin(); it != allTreeStrings.end(); it++) {
        treeObject.readTreeString(*it);
        double score = treeObject.optimizeAllBranches(1);
        update(treeObject.getTreeString(), score);
    }
}









