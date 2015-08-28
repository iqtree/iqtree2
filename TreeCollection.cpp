//
// Created by Tung Nguyen on 6/23/15.
//

#include "TreeCollection.h"

using namespace std;

TreeCollection::TreeCollection(CandidateSet &candidateTrees) {
    numTrees = candidateTrees.size();
    CandidateSet::reverse_iterator rit;
    for (rit = candidateTrees.rbegin(); rit != candidateTrees.rend(); rit++) {
       treeStrings.push_back(rit->second.tree);
       scores.push_back(rit->first);
    }
}

TreeCollection::TreeCollection(vector<string>& trees, vector<double>& scores) {
    assert(trees.size() == scores.size());
    this->treeStrings = trees;
    this->scores = scores;
    numTrees = trees.size();
}

pair<string, double> TreeCollection::getTree(int i) {
    return std::make_pair(treeStrings[i], scores[i]);
}

void TreeCollection::clear() {
    treeStrings.clear();
    scores.clear();
}

void TreeCollection::addTrees(TreeCollection &trees) {
    for (int i = 0; i < trees.getNumTrees(); i++) {
        treeStrings.push_back(trees.getTree(i).first);
        scores.push_back(trees.getTree(i).second);
    }
}

CandidateSet TreeCollection::convertToCandidateSet() {
    for (int i = 0; i < )
    return CandidateSet();
}
