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
GenomeNode* GenomeTree::findNodeByPos(GenomeNode* node, Insertion insertion)
{
    // return NULL if not found
    if (!node)
        outError("Opps! Insertion occurs at an invalid position. There is something wrong!");
    
    // check if the current node contains the given pos
    // or if insertion occur at the end of the current genome (~append) -> return the right-bottom node of the tree
    if ((node->pos_new <= insertion.pos && insertion.pos < node->pos_new + node->length)
        || (insertion.is_append && insertion.pos == node->pos_new + node->length))
        return node;
    
    // otherwise, search on left/right branches
    if (insertion.pos < node->pos_new)
        return findNodeByPos(node->left_child, insertion);
    else
        return findNodeByPos(node->right_child, insertion);
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
    int new_pos_ori = node->pos_ori;
    int new_length = pos - node->pos_new;
    int new_pos_new = node->pos_new;
    GenomeNode* left_child = new GenomeNode(new_pos_ori, new_length, new_pos_new);
    
    // init the right child
    new_pos_ori = node->pos_ori + left_child->length;
    new_length = node->length - left_child->length;
    new_pos_new = pos + length;
    GenomeNode* right_child = new GenomeNode(new_pos_ori, new_length, new_pos_new);
    
    // update the current node
    node->pos_ori = -1;
    node->length = length;
    node->pos_new = pos;
    
    // update relations among those nodes
    node->left_child = left_child;
    node->right_child = right_child;
    left_child->parent = node;
    right_child->parent = node;
}

/**
    recursively traverse the tree to update pos_new after insert gaps into a node
 */
void GenomeTree::updatePosNew(GenomeNode* node, int pos_new, int length)
{
    // make sure node is not NULL
    ASSERT(node);
    
    // stop traversing deeply if we found the node where gaps are inserted.
    if (pos_new == node->pos_new)
        return;
    
    // update the current node if its new pos is greater than pos
    if (pos_new < node->pos_new)
    {
        node->pos_new = node->pos_new + length;
        
        // browse the left child if it exists
        if (node->left_child)
            updatePosNew(node->left_child, pos_new, length);
    }
    
    // browse the right child if it exists
    if (node->right_child)
        updatePosNew(node->right_child, pos_new, length);
}

/**
    update genome tree from a vector of insertions
 */
void GenomeTree::updateTree(vector<Insertion> insertions)
{
    // process insertions one by one
    for (Insertion insertion:insertions)
    {
        // find the genome node contains the position to insert gaps
        GenomeNode* node = findNodeByPos(root, insertion);
        
        // add gaps into that node we found
        ASSERT(node);
        // if the current node contains only gaps -> simple inscrease its length, and update the pos_new of its children
        if (node->pos_ori == -1)
        {
            insertGapsIntoGaps(node, insertion.length);
            updatePosNew(node->right_child, node->pos_new, insertion.length);
        }
        // otherwise, insert gaps into a "normal" node (containing sites from the original genome) -> build a cherry from the current node
        else
            insertGapsIntoNormalNode(node, insertion.pos, insertion.length);
        
        // traverse the tree from root to update the positions in the new genome of all nodes that follows the current node. Note that chilren of the node were already updated thus, we will stop at the current node in this method.
        updatePosNew(root, node->pos_new, insertion.length);
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
    exportSiteFromGenomeNode(root, ori_seq, new_seq);
    
    return new_seq;
}

/**
    export sites from a genome node and the original genome
 */
void GenomeTree::exportSiteFromGenomeNode(GenomeNode* node, vector<short int> ori_seq, vector<short int> &new_seq)
{
    // export sites from the current genome node (skip nodes with gaps as the default states are gaps)
    if (node->pos_ori > -1 && node->length > 0)
    {
        ASSERT(node->pos_ori + node->length <= ori_seq.size());
        ASSERT(node->pos_new + node->length <= new_seq.size());
        
        for (int i = 0; i < node->length; i++)
            new_seq[node->pos_new + i] = ori_seq[node->pos_ori + i];
    }
    
    // traverse the left and right children (if any)
    if (node->left_child)
        exportSiteFromGenomeNode(node->left_child , ori_seq, new_seq);
    if (node->right_child)
        exportSiteFromGenomeNode(node->right_child , ori_seq, new_seq);
}

/**
 export readable characters (for writing to file) from original genome and genome tree
 */
string GenomeTree::exportReadableCharacters(vector<short int> ori_seq, int sequence_length, int num_sites_per_state, vector<string> state_mapping)
{
    // dummy variables
    std::string output (sequence_length * num_sites_per_state + 1, '-');
    output[output.length()-1] = '\n';
    
    // recursively traverse the genome tree to export readable characters from each genome node and the original genome
    exportReadableCharactersFromGenomeNode(root, ori_seq, output, num_sites_per_state, state_mapping);
    
    // return output
    return output;
}

/**
    export readable characters from a genome node and the original genome
 */
void GenomeTree::exportReadableCharactersFromGenomeNode(GenomeNode* node, vector<short int> ori_seq, string &output, int num_sites_per_state, vector<string> state_mapping)
{
    // export readable characters from the current genome node (skip nodes with gaps as the default states are gaps)
    if (node->pos_ori > -1 && node->length > 0)
    {
        ASSERT(node->pos_ori + node->length <= ori_seq.size());
        ASSERT((node->pos_new + node->length) * num_sites_per_state <= output.length());
        
        for (int i = 0; i < node->length; i++)
        {
            // convert normal data
            if (num_sites_per_state == 1)
                output[node->pos_new + i] = state_mapping[ori_seq[node->pos_ori + i]][0];
            // convert CODON
            else
            {
                string output_codon = state_mapping[ori_seq[node->pos_ori + i]];
                output[(node->pos_new + i) * num_sites_per_state] = output_codon[0];
                output[(node->pos_new + i) * num_sites_per_state + 1] = output_codon[1];
                output[(node->pos_new + i) * num_sites_per_state + 2] = output_codon[2];
            }
        }
    }
    
    // traverse the left and right children (if any)
    if (node->left_child)
        exportReadableCharactersFromGenomeNode(node->left_child, ori_seq, output, num_sites_per_state, state_mapping);
    if (node->right_child)
        exportReadableCharactersFromGenomeNode(node->right_child, ori_seq, output, num_sites_per_state, state_mapping);
}
