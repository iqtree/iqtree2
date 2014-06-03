/*
 * candidateset.h
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#ifndef CANDIDATESET_H_
#define CANDIDATESET_H_
#include "tools.h"

class candidateset {
    struct tree_compare {
        bool operator()(const std::pair<string, double> &left, const std::pair<int,int> &right) {
            return left.second < right.second;
        }
    };
public:
    candidateset(int max_elements, int max_parents);

    string getBestTree();

    string getNextTree();

    bool insert(string treeString, double score);

    virtual ~candidateset();
private:
    int max_elements;

    /**
     *  maximum number of parents
     */
    int max_parents;

    int nextTreeIndex;

    set< pair<string, double>, tree_compare > candidateTrees;

    unordered_set<string> allTopos;

    bool treeExist(string treeString);

};

#endif /* CANDIDATESET_H_ */
