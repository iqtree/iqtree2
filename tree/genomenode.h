//
//  genomenode.h
//  tree
//
//  Created by Nhan Ly-Trong on 09/02/2022.
//
#ifndef GENOMENODE_H
#define GENOMENODE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "utils/tools.h"
using namespace std;

/*--------------------------------------------------------------*/
class Insertion {
public:
    int pos;
    int length;
    bool is_append;
    
    /**
        constructor
     */
    Insertion()
    {
        pos = 0;
        length = 0;
        is_append = false;
    }
    
    /**
        constructor
     */
    Insertion(int n_pos, int n_length, bool n_is_append = false)
    {
        pos = n_pos;
        length = n_length;
        is_append = n_is_append;
    }
};
/*--------------------------------------------------------------*/

/**
A Genome Node to present a genome entry (a set of sites)
 */
class GenomeNode {
public:
    /**
        TRUE if this node contains gaps
     */
    bool is_gap;
    
    /**
        starting pos in the original genome
     */
    int pos_ori;

    /**
        length (the number of sites)
     */
    int length;
    
    /**
        total gaps from the left child branch
     */
    int cumulative_gaps_from_left_child;
    
    /**
        total gaps from the parent (NOTE: need to pre-compute before using)
     */
    int cumulative_gaps_from_parent;
    
    /**
        parent, left/right children of this node
     */
    GenomeNode *parent, *left_child, *right_child;

    /**
        constructor
     */
    GenomeNode();
    
    /**
        init a root genome node
     */
    GenomeNode(int length);
    
    /**
        constructor
     */
    GenomeNode(bool n_is_gap, int n_pos_ori, int n_length, int n_cumulative_gaps);
    
    /**
        deconstructor
     */
    ~GenomeNode();
    
    /**
        update relations
     */
    void updateRelation(GenomeNode* n_parent, GenomeNode* n_left_child, GenomeNode* n_right_child);
};
#endif
