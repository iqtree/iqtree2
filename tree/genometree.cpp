//
//  genometree.cpp
//  tree
//
//  Created by Nhan Ly-Trong on 09/02/2022.
//

#include "genometree.h"

/**
    constructor
 */
GenomeTree::GenomeTree()
{
    root = NULL;
}

/**
    init a root genome node
 */
GenomeTree::GenomeTree(int length)
{
    root = new GenomeNode(length);
}

/**
    deconstructor
 */
GenomeTree::~GenomeTree()
{
    if (root && root->left_child)
    {
        delete root->left_child;
        root->left_child = NULL;
    }
    
    if (root && root->right_child)
    {
        delete root->right_child;
        root->right_child = NULL;
    }
}

/**
    find a node that contains a given position
 */
GenomeNode* GenomeTree::findNodeByPos(GenomeNode* node, Insertion* insertion, int num_cumulative_gaps)
{
    // return NULL if not found
    if (!node)
        outError("Opps! Insertion occurs at an invalid position. There is something wrong!");
    
    // compute the position in the new genome (with gaps)
    int pos_new = node->pos_ori + num_cumulative_gaps + node->cumulative_gaps_from_left_child;
    
    // check if the current node contains the given pos
    // or if insertion occur at the end of the current genome (~append) -> return the right-bottom node of the tree
    if ((pos_new <= insertion->pos && insertion->pos < pos_new + node->length)
        || (insertion->is_append && insertion->pos == pos_new + node->length))
    {
        // update cumulative_gaps_from_parent before returning node
        node->cumulative_gaps_from_parent = num_cumulative_gaps;
        return node;
    }
    
    // otherwise, search on left/right branches
    if (insertion->pos < pos_new)
        return findNodeByPos(node->left_child, insertion, num_cumulative_gaps);
    else
        return findNodeByPos(node->right_child, insertion, num_cumulative_gaps + node->length + node->cumulative_gaps_from_left_child);
}

/**
    insert gaps into a "all-gap" node
 */
void GenomeTree::insertGapsIntoGaps(GenomeNode* node, int length)
{
    node->length = node->length + length;
}

/**
    insert gaps into a "normal" (all sites from the original genome) node => init a new cherry from a node
 */
void GenomeTree::insertGapsIntoNormalNode(GenomeNode* node, int pos, int length)
{
    // init the left child
    bool new_is_gap = false;
    int new_pos_ori = node->pos_ori;
    int new_length = pos - (node->pos_ori + node->cumulative_gaps_from_parent);
    int new_cumulative_gap_from_left_child = 0;
    GenomeNode* left_child = new GenomeNode(new_is_gap, new_pos_ori, new_length, new_cumulative_gap_from_left_child);
    
    // init the right child
    new_is_gap = false;
    new_pos_ori = left_child->pos_ori + left_child->length;
    new_length = node->length - left_child->length;
    new_cumulative_gap_from_left_child = 0;
    GenomeNode* right_child = new GenomeNode(new_is_gap, new_pos_ori, new_length, new_cumulative_gap_from_left_child);
    
    // update the current node
    node->is_gap = true;
    node->pos_ori = right_child->pos_ori;
    node->length = length;
    node->cumulative_gaps_from_left_child = 0;
    
    // update relations among those nodes
    node->left_child = left_child;
    node->right_child = right_child;
    left_child->parent = node;
    right_child->parent = node;
}

/**
    update the cumulative_gaps_from_left_child of all nodes on the path from the current node to root
 */
void GenomeTree::updateCumulativeGapsFromLeftChild(GenomeNode* node, int length)
{
    // make sure the current node is not NULL
    ASSERT(node);

    // stop at root
    if (node->parent)
    {
        // if the current node is left child of its parent -> update the cumulative_gaps_from_left_child
        if (node->parent->left_child == node)
            node->parent->cumulative_gaps_from_left_child += length;
        
        // go up the tree
        updateCumulativeGapsFromLeftChild(node->parent, length);
    }
}

/**
    update genome tree from a vector of insertions
 */
void GenomeTree::updateTree(Insertion* prev_insertion)
{
    // Do nothing if there is no insertion event
    if (!prev_insertion)
        return;
    
    // init current insertion
    Insertion* current_insertion = prev_insertion->next;
    
    // process insertions one by one
    for (; current_insertion;)
    {
        // find the genome node contains the position to insert gaps
        GenomeNode* node = findNodeByPos(root, current_insertion, 0);
        
        // add gaps into that node we found
        ASSERT(node);
        
        // if the current node contains only gaps -> simple inscrease its length, and update the pos_new of its children
        if (node->is_gap)
            insertGapsIntoGaps(node, current_insertion->length);
        // otherwise, insert gaps into a "normal" node (containing sites from the original genome) -> build a cherry from the current node
        else
            insertGapsIntoNormalNode(node, current_insertion->pos, current_insertion->length);
        
        // update the cumulative_gaps_from_left_child of all nodes on the path from the current node to root
        updateCumulativeGapsFromLeftChild(node, current_insertion->length);
        
        // move to the next insertion
        current_insertion = current_insertion->next;
    }
}

/**
    export new genome from original genome and genome tree
 */
vector<short int> GenomeTree::exportNewGenome(vector<short int> ori_seq, int seq_length, int UNKOWN_STATE)
{
    // init new genome
    vector<short int> new_seq(seq_length, UNKOWN_STATE);
    
    // traverse the genome tree to export new genome
    queue<GenomeNode*> genome_nodes;
    root->cumulative_gaps_from_parent = 0;
    genome_nodes.push(root);
    
    while (!genome_nodes.empty()) {
        GenomeNode* node = genome_nodes.front();
        genome_nodes.pop();
        
        // export sites from the current genome node (skip nodes with gaps as the default states are gaps)
        if (!node->is_gap && node->length > 0)
        {
            int pos_new = node->pos_ori + node->cumulative_gaps_from_parent + node->cumulative_gaps_from_left_child;
            ASSERT(node->is_gap || node->pos_ori + node->length <= ori_seq.size());
            ASSERT(pos_new + node->length <= new_seq.size());
            
            for (int i = 0; i < node->length; i++)
                new_seq[pos_new + i] = ori_seq[node->pos_ori + i];
        }
        
        // traverse the left and right children (if any)
        if (node->left_child)
        {
            node->left_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent;
            genome_nodes.push(node->left_child);
        }
        if (node->right_child)
        {
            node->right_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent + node->length + node->cumulative_gaps_from_left_child;
            genome_nodes.push(node->right_child);
        }
    }
    
    return new_seq;
}

/**
 export readable characters (for writing to file) from original genome and genome tree
 */
void GenomeTree::exportReadableCharacters(vector<short int> ori_seq, int num_sites_per_state, vector<string> state_mapping, string &output)
{
    queue<GenomeNode*> genome_nodes;
    root->cumulative_gaps_from_parent = 0;
    genome_nodes.push(root);
    
    while (!genome_nodes.empty()) {
        GenomeNode* node = genome_nodes.front();
        genome_nodes.pop();
        
        // export readable characters from the current genome node (skip nodes with gaps as the default states are gaps)
        if (!node->is_gap && node->length > 0)
        {
            int pos_new = node->pos_ori + node->cumulative_gaps_from_parent + node->cumulative_gaps_from_left_child;
            ASSERT(node->is_gap || node->pos_ori + node->length <= ori_seq.size());
            ASSERT((pos_new + node->length) * num_sites_per_state <= output.length());
            
            for (int i = 0; i < node->length; i++)
            {
                // convert normal data
                if (num_sites_per_state == 1)
                    output[pos_new + i] = state_mapping[ori_seq[node->pos_ori + i]][0];
                // convert CODON
                else
                {
                    string output_codon = state_mapping[ori_seq[node->pos_ori + i]];
                    output[(pos_new + i) * num_sites_per_state] = output_codon[0];
                    output[(pos_new + i) * num_sites_per_state + 1] = output_codon[1];
                    output[(pos_new + i) * num_sites_per_state + 2] = output_codon[2];
                }
            }
        }
        
        // traverse the left and right children (if any)
        if (node->left_child)
        {
            node->left_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent;
            genome_nodes.push(node->left_child);
        }
        if (node->right_child)
        {
            node->right_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent + node->length + node->cumulative_gaps_from_left_child;
            genome_nodes.push(node->right_child);
        }
    }
}
