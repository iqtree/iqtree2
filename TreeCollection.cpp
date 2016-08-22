//
// Created by Tung Nguyen on 6/23/15.
//

#include "TreeCollection.h"
#include "MPIHelper.h"

using namespace std;

TreeCollection::TreeCollection(CandidateSet &candidateTrees) {
    CandidateSet::reverse_iterator rit;
    for (rit = candidateTrees.rbegin(); rit != candidateTrees.rend(); rit++) {
       treeStrings.push_back(rit->second.tree);
       scores.push_back(rit->first);
       sourceProcID.push_back(MPIHelper::getInstance().getProcessID());
    }
}

TreeCollection::TreeCollection(vector<string>& trees, vector<double>& scores, vector<int> &sourceProcID) {
    assert(trees.size() == scores.size());
    this->treeStrings = trees;
    this->scores = scores;
    this->sourceProcID = sourceProcID;
//    this->sourceProcID.clear();
//    this->sourceProcID.insert(this->sourceProcID.end(), scores.size(), MPIHelper::getInstance().getProcessID());
}

pair<string, double> TreeCollection::getTree(int i) {
    assert(treeStrings.size() == scores.size());
    return std::make_pair(treeStrings[i], scores[i]);
}

void TreeCollection::clear() {
    treeStrings.clear();
    scores.clear();
}

void TreeCollection::addTrees(TreeCollection &trees) {
//    for (int i = 0; i < trees.getNumTrees(); i++) {
//        treeStrings.push_back(trees.getTree(i).first);
//        scores.push_back(trees.getTree(i).second);
//        
//    }
    treeStrings.insert(treeStrings.end(), trees.treeStrings.begin(), trees.treeStrings.end());
    scores.insert(scores.end(), trees.scores.begin(), trees.scores.end());
    sourceProcID.insert(sourceProcID.end(), trees.sourceProcID.begin(), trees.sourceProcID.end());
}

size_t TreeCollection::getNumTrees() {
    size_t numTrees = treeStrings.size();
    assert(numTrees == scores.size());
    return numTrees;
}
