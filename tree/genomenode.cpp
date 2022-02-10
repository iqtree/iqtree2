//
//  genomenode.cpp
//  tree
//
//  Created by Nhan Ly-Trong on 09/02/2022.
//

#include "genomenode.h"

/**
    constructor
 */
GenomeNode::GenomeNode()
{
    pos_ori = -1;
    length = 0;
    pos_new = 0;
    
    parent = NULL;
    left_child = NULL;
    right_child = NULL;
}

/**
    constructor
 */
GenomeNode::GenomeNode(int n_pos_ori, int n_length, int n_pos_new)
{
    pos_ori = n_pos_ori;
    length = n_length;
    pos_new = n_pos_new;
    
    parent = NULL;
    left_child = NULL;
    right_child = NULL;
}

/**
    init a root genome node
 */
GenomeNode::GenomeNode(int length):GenomeNode(0,length,0){}

/**
    deconstructor
 */
GenomeNode::~GenomeNode()
{
    if (left_child)
    {
        delete left_child;
        left_child = NULL;
    }
    
    if (right_child)
    {
        delete right_child;
        right_child = NULL;
    }
}

/**
    update relations
 */
void GenomeNode::updateRelation(GenomeNode* n_parent, GenomeNode* n_left_child, GenomeNode* n_right_child)
{
    parent = n_parent;
    left_child = n_left_child;
    right_child = n_right_child;
}
