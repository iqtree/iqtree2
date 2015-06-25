//
// Created by tung on 6/23/15.
//

#ifndef IQTREE_TREECOLLECTION_H
#define IQTREE_TREECOLLECTION_H
#include "candidateset.h"

class TreeCollection {

public:

    /**
     *  Constructor
     */
    TreeCollection(vector<CandidateTree>& candidateTrees);

    TreeCollection() {};

    TreeCollection(vector<string>& trees, vector<double>& scores);

    void addTrees(TreeCollection &trees);

    /*
     *  Get i-th tree and its score
    */
    pair<string, double> getTree(int i);

    void clear();

    void setTreeStrings(const vector<string> &treeStrings) {
        TreeCollection::treeStrings = treeStrings;
    }

    void setScores(const vector<double> &scores) {
        TreeCollection::scores = scores;
    }

    size_t getNumTrees() const {
        return numTrees;
    }

    const vector<string> &getTreeStrings() const {
        return treeStrings;
    }

    const vector<double> &getScores() const {
        return scores;
    }

private:
    size_t numTrees;
    vector<string> treeStrings;
    vector<double> scores;
};


#endif //IQTREE_TREECOLLECTION_H
