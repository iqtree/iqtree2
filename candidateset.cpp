/*
 * candidateset.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#include "phylotree.h"
#include "candidateset.h"

CandidateSet::CandidateSet(int limit, int max_candidates, Alignment *aln) {
    assert(max_candidates <= limit);
    assert(aln);
    this->limit = limit;
    this->max_candidates = max_candidates;
    this->aln = aln;
}

CandidateSet::CandidateSet() {
	aln = NULL;
	limit = 0;
	max_candidates = 0;
}

string CandidateSet::getBestTree() {
	// TODO: What happen if there are multiple optimal trees?
	return (rbegin())->second.tree;
}

string CandidateSet::getRandCandTree() {
	if (empty())
		return "";
	// BQM: bug fix max -> min
	//int id = random_int(min(max_candidates, (int)size()) );
	int id = randint(0, min(max_candidates, (int)size()));
	//int id = 0 + (rand() % (int)( min(max_candidates, (int)size())) );
	for (reverse_iterator i = rbegin(); i != rend(); i++, id--)
		if (id == 0)
			return i->second.tree;
	return "";
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
        int count = this->max_candidates;
        for (reverse_iterator i = rbegin(); i != rend() && count >0 ; i++, count--) {
            parentTrees.push(i->second.tree);
            //cout << i->first << endl;
        }
    }
}

bool CandidateSet::update(string tree, double score) {
	CandidateTree candidate;
	candidate.tree = tree;
	candidate.score = score;
	candidate.topology = getTopology(tree);
	if (treeTopologyExist(candidate.topology)) {
	    // if tree topology already exist, we replace the old
	    // by the new one (with new branch lengths)
		if (topologies[candidate.topology] < score) {
			topologies[candidate.topology] = score;
			for (CandidateSet::iterator i = begin(); i != end(); i++)
				if (i->second.topology == candidate.topology) {
					erase(i);
					break;
				}
			// insert tree into candidate set
			insert(CandidateSet::value_type(score, candidate));
		}
		return false;
	}
	if (size() < limit) {
		// insert tree into candidate set
		insert(CandidateSet::value_type(score, candidate));
		topologies[candidate.topology] = score;
		return true;
	} else if (begin()->first < score){
		// remove the worst-scoring tree
		topologies.erase(begin()->second.topology);
		erase(begin());
		// insert tree into candidate set
		insert(CandidateSet::value_type(score, candidate));
		topologies[candidate.topology] = score;
		return true;
	}
	return false;
}

void CandidateSet::printBestScores() {
	int cnt = max_candidates;
	for (reverse_iterator rit = rbegin(); rit != rend() && cnt > 0; rit++, cnt--) {
		if (cnt < max_candidates) cout << " / ";
		cout << rit->first;
	}
	cout << endl;
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

CandidateSet::~CandidateSet() {
}

bool CandidateSet::treeTopologyExist(string topo) {
	return topologies.find(topo) != topologies.end();
}

bool CandidateSet::treeExist(string tree) {
	return treeTopologyExist(getTopology(tree));
}
