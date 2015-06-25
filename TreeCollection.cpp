//
// Created by tung on 6/23/15.
//

#include "TreeCollection.h"

using namespace std;

TreeCollection::TreeCollection(vector<CandidateTree> &candidateTrees) {
    numTrees = candidateTrees.size();
    vector<CandidateTree>::iterator it;
    for (it = candidateTrees.begin(); it != candidateTrees.end(); it++) {
       treeStrings.push_back(it->tree);
       scores.push_back(it->score);
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
