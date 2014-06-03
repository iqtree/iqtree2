/*
 * candidateset.cpp
 *
 *  Created on: Jun 1, 2014
 *      Author: Tung Nguyen
 */

#include "candidateset.h"

candidateset::candidateset(int max_elements, int max_parents) {
    this->max_parents = max_parents;
    this->max_elements = max_elements;
}

string candidateset::getBestTree() {
}


bool candidateset::insert(string treeString, double score) {
    if (candidateTrees.size() < max_parents) {

    }
}

candidateset::~candidateset() {
    // TODO Auto-generated destructor stub
}

bool candidateset::treeExist(string treeString) {
    PhyloTree tree;
    tree.readTreeString(treeString);

}
