//
// Created by tung on 6/23/15.
//

#ifndef IQTREE_TREECOLLECTION_H
#define IQTREE_TREECOLLECTION_H
#include "tree/candidateset.h"

/**
 *  A container for a set of trees together with their scores
 */
class TreeCollection {
private:
    vector<string> treeStrings;
    vector<double> scores;
    vector<int> sourceProcID;
public:

    /**
     *  Constructor
     */
    TreeCollection() {};

    TreeCollection(vector<string>& trees, vector<double>& scores, vector<int> &sourceProcID);

    void addTrees(TreeCollection &trees);

    void addTrees(CandidateSet& candidateTrees);


    /*
     *  Get i-th tree and its score
    */
    pair<string, double> getTree(int i);

    void clear();

    void setTreeStrings(const vector<string> treeStrings) {
        TreeCollection::treeStrings = treeStrings;
    }

    void setScores(const vector<double> scores) {
        TreeCollection::scores = scores;
    }

    size_t getNumTrees();

    const vector<string> &getTreeStrings() const {
        return treeStrings;
    }

    const vector<double> &getScores() const {
        return scores;
    }

    const vector<int> &getSourceProcID() const {
        return sourceProcID;
    }

};


#endif //IQTREE_TREECOLLECTION_H
