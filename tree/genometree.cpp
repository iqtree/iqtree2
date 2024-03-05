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
    if (root)
    {
        delete root;
        root = NULL;
    }
}

/**
    find a node that contains a given position
 */
GenomeNode* GenomeTree::findNodeByPos(GenomeNode* node, Insertion* insertion, int num_cumulative_gaps_from_parents, int num_cumulative_converts_from_parents)
{
    // return NULL if not found
    if (!node)
        outError("Opps! Insertion occurs at an invalid position. There is something wrong!");
    
    // compute the position in the new genome (with gaps)
    int pos_ori = node->pos_ori + num_cumulative_converts_from_parents + node->cumulative_converts_from_left_child;
    int pos_new = pos_ori + num_cumulative_gaps_from_parents + node->cumulative_gaps_from_left_child;
    
    // check if the current node contains the given pos
    // or if insertion occur at the end of the current genome (~append) -> return the right-bottom node of the tree
    if ((pos_new <= insertion->pos && insertion->pos < pos_new + node->length)
        || (insertion->is_append && insertion->pos == pos_new + node->length))
    {
        // update cumulative_gaps_from_parent before returning node
        node->cumulative_gaps_from_parent = num_cumulative_gaps_from_parents;
        node->cumulative_converts_from_parent = num_cumulative_converts_from_parents;
        return node;
    }
    
    // otherwise, search on left/right branches
    if (insertion->pos < pos_new)
        return findNodeByPos(node->left_child, insertion, num_cumulative_gaps_from_parents, num_cumulative_converts_from_parents);
    else
        return findNodeByPos(node->right_child, insertion, num_cumulative_gaps_from_parents + (node->type == GAP ? node->length : 0) + node->cumulative_gaps_from_left_child, num_cumulative_converts_from_parents + (node->type == GAP_CONVERTED_2_NORMAL ? node->length : 0) + node->cumulative_converts_from_left_child);
}

/**
    insert gaps into a "normal" (all sites from the original genome) node => init a new cherry from a node
 */
void GenomeTree::insertGapsIntoSequenceNode(GenomeNode* node, Insertion* insertion, bool attach_insertion)
{
    // backup the left_child and right_child of the current node
    GenomeNode* left_child = node->left_child;
    GenomeNode* right_child = node->right_child;
    Insertion* prev_insertion = node->insertion;
    
    // init the left child
    GenomeNodeType new_type = node->type;
    int new_pos_ori = node->pos_ori;
    int new_length = insertion->pos - (node->pos_ori + node->cumulative_gaps_from_parent + node->cumulative_gaps_from_left_child);
    int new_cumulative_gap_from_left_child = node->cumulative_gaps_from_left_child ;
    GenomeNode* new_left_child = new GenomeNode(new_type, new_pos_ori, new_length, new_cumulative_gap_from_left_child);
    
    // init the right child
    new_type = node->type;
    new_pos_ori = new_left_child->pos_ori + (new_left_child->type == GAP ? 0 : new_left_child->length);
    new_length = node->length - new_left_child->length;
    new_cumulative_gap_from_left_child = 0;
    GenomeNode* new_right_child = new GenomeNode(new_type, new_pos_ori, new_length, new_cumulative_gap_from_left_child);
    
    // update the current node
    node->type = GAP;
    node->pos_ori = new_right_child->pos_ori;
    node->length = insertion->length;
    node->cumulative_gaps_from_left_child = new_left_child->cumulative_gaps_from_left_child + (new_left_child->type == GAP ? new_left_child->length : 0);
    node->cumulative_converts_from_left_child = 0;
    if (attach_insertion)
    {
        node->insertion = insertion;
        insertion->genome_nodes.push_back(node);
    }
    
    // update relations among those nodes
    node->left_child = new_left_child;
    node->right_child = new_right_child;
    new_left_child->parent = node;
    new_right_child->parent = node;
    // attach left_child and right_child to the new_left_child, and new_right_child
    if (left_child && right_child)
    {
        new_left_child->left_child = left_child;
        new_left_child->right_child = new GenomeNode(NORMAL, node->pos_ori, 0, 0);
        new_left_child->left_child->parent = new_left_child;
        new_left_child->right_child->parent = new_left_child;
        
        new_right_child->left_child = new GenomeNode(NORMAL, node->pos_ori, 0, 0);
        new_right_child->right_child = right_child;
        new_right_child->left_child->parent = new_right_child;
        new_right_child->right_child->parent = new_right_child;
        
        if (attach_insertion)
        {
            // update list of all genome nodes that contains gaps inserted by prev_insertion
            // find and replace the current node in the genome_nodes of prev_insertions by the new_left_child
            auto iter = std::find(prev_insertion->genome_nodes.begin(), prev_insertion->genome_nodes.end(), node);
            if (iter != prev_insertion->genome_nodes.end())
                *iter = new_left_child;
            // add the new_right_child into the genome_nodes of prev_insertions
            prev_insertion->genome_nodes.push_back(new_right_child);
            
            // indicate that gaps in the new left/child was inserted by the prev_insertion
            new_left_child->insertion = prev_insertion;
            new_right_child->insertion = prev_insertion;
        }
    }
}

/**
    build a genome tree from an insertion forward the insertion list
 */
void GenomeTree::buildGenomeTree(Insertion* prev_insertion, int ori_seq_length, bool attach_insertion)
{
    // Do nothing if there is no insertion event
    if (!prev_insertion)
        return;
    
    // init root node if it's NULL
    if (!root)
        root = new GenomeNode(ori_seq_length);
    
    // init current insertion
    Insertion* current_insertion = prev_insertion->next;
    
    // process insertions one by one
    for (; current_insertion;)
    {
        // find the genome node contains the position to insert gaps
        GenomeNode* node = findNodeByPos(root, current_insertion, 0, 0);
        
        // insert gaps into a sequence node -> build a cherry from the current node
        insertGapsIntoSequenceNode(node, current_insertion, attach_insertion);
        
        // update the cumulative_gaps_from_left_child of all nodes on the path from the current node to root
        for (; node->parent; )
        {
            // if the current node is left child of its parent -> update the cumulative_gaps_from_left_child
            if (node->parent->left_child == node)
                node->parent->cumulative_gaps_from_left_child += current_insertion->length;
            
            // go up the tree
            node = node->parent;
        }
        
        // move to the next insertion
        current_insertion = current_insertion->next;
    }
}

/**
    update tree by accepting gaps (generated by a set of previous insertions) as normal characters
 */
void GenomeTree::updateGenomeTree(Insertion* start_insertion, Insertion* end_insertion)
{
    if (!start_insertion->next) return;
    
    // handle insertions one by one
    for (Insertion* insertion = start_insertion->next; insertion != end_insertion->next; )
    {
        // convert genome nodes (gaps generated by the current insertion) into normal characters one by one.
        for (int i = 0; i < insertion->genome_nodes.size(); i++)
            // convert gaps into normal characters
            convertGapsIntoNormal(insertion->genome_nodes[i]);
        
        // move to the next insertion
        insertion = insertion->next;
    }
}

/**
    convert gaps into normal characters
 */
void GenomeTree::convertGapsIntoNormal(GenomeNode *node)
{
    ASSERT(node->type == GAP);
    
    // change genome type
    node->type = GAP_CONVERTED_2_NORMAL;
    
    // get number of gaps that was converted into normal characters
    int num_converts = node->length;
    
    // update the cumulative_converts, and cumulative_gaps of all left parents (on the path from the current node to root)
    // stop at root
    for (; node->parent; )
    {
        // if the current node is left child of its parent -> update the cumulative_gaps_from_left_child, and cumulative_converts
        if (node->parent->left_child == node)
        {
            node->parent->cumulative_gaps_from_left_child -= num_converts;
            node->parent->cumulative_converts_from_left_child += num_converts;
        }
        
        // go upward the tree
        node = node->parent;
    }
}

/**
    export new genome from original genome and genome tree
 */
vector<short int> GenomeTree::exportNewGenome(vector<short int> &ori_seq, int seq_length, int UNKOWN_STATE)
{
    // init new genome
    vector<short int> new_seq(seq_length, UNKOWN_STATE);
    
    // traverse the genome tree to export new genome
    queue<GenomeNode*> genome_nodes;
    root->cumulative_gaps_from_parent = 0;
    root->cumulative_converts_from_parent = 0;
    genome_nodes.push(root);
    
    while (!genome_nodes.empty()) {
        GenomeNode* node = genome_nodes.front();
        genome_nodes.pop();
        
        // export sites from the current genome node (skip nodes with gaps as the default states are gaps)
        if (node->type != GAP && node->length > 0)
        {
            int pos_ori = node->pos_ori + node->cumulative_converts_from_parent + node->cumulative_converts_from_left_child;
            int pos_new = pos_ori + node->cumulative_gaps_from_parent + node->cumulative_gaps_from_left_child;
            ASSERT(pos_ori + node->length <= ori_seq.size());
            ASSERT(pos_new + node->length <= new_seq.size());
            
            for (int i = 0; i < node->length; i++)
                new_seq[pos_new + i] = ori_seq[pos_ori + i];
        }
        
        // traverse the left and right children (if any)
        if (node->left_child)
        {
            node->left_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent;
            node->left_child->cumulative_converts_from_parent = node->cumulative_converts_from_parent;
            genome_nodes.push(node->left_child);
        }
        if (node->right_child)
        {
            node->right_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent + (node->type == GAP ? node->length : 0) + node->cumulative_gaps_from_left_child;
            node->right_child->cumulative_converts_from_parent = node->cumulative_converts_from_parent + (node->type == GAP_CONVERTED_2_NORMAL ? node->length : 0) + node->cumulative_converts_from_left_child;
            genome_nodes.push(node->right_child);
        }
    }
    
    return new_seq;
}

/**
 export readable characters (for writing to file) from original genome and genome tree
 */
void GenomeTree::exportReadableCharacters(vector<short int> &ori_seq, int num_sites_per_state, vector<string> &state_mapping, string &output)
{
    queue<GenomeNode*> genome_nodes;
    root->cumulative_gaps_from_parent = 0;
    root->cumulative_converts_from_parent = 0;
    genome_nodes.push(root);
    
    while (!genome_nodes.empty()) {
        GenomeNode* node = genome_nodes.front();
        genome_nodes.pop();
        
        // export readable characters from the current genome node (skip nodes with gaps as the default states are gaps)
        if (node->type != GAP && node->length > 0)
        {
            int pos_ori = node->pos_ori + node->cumulative_converts_from_parent + node->cumulative_converts_from_left_child;
            int pos_new = pos_ori + node->cumulative_gaps_from_parent + node->cumulative_gaps_from_left_child;
            ASSERT(pos_ori + node->length <= ori_seq.size());
            ASSERT((num_sites_per_state == 1 ? (pos_new + node->length) : ((pos_new + node->length) * num_sites_per_state)) <= output.length());
            
            // NHANLT: potential improvement
            // change to use pointer to avoid accessing []
            // convert normal data
            if (num_sites_per_state == 1)
            {
                for (int i = 0; i < node->length; i++)
                {
                    output[pos_new + i] = state_mapping[ori_seq[pos_ori + i]][0];
                }
            }
            // convert CODON
            else
            {
                int index = pos_new * num_sites_per_state;
                for (int i = 0; i < node->length; i++, index += num_sites_per_state)
                {
                    string output_codon = state_mapping[ori_seq[pos_ori + i]];
                    output[index] = output_codon[0];
                    output[index + 1] = output_codon[1];
                    output[index + 2] = output_codon[2];
                }
            }
        }
        
        // traverse the left and right children (if any)
        if (node->left_child)
        {
            node->left_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent;
            node->left_child->cumulative_converts_from_parent = node->cumulative_converts_from_parent;
            genome_nodes.push(node->left_child);
        }
        if (node->right_child)
        {
            node->right_child->cumulative_gaps_from_parent = node->cumulative_gaps_from_parent + (node->type == GAP ? node->length : 0) + node->cumulative_gaps_from_left_child;
            node->right_child->cumulative_converts_from_parent = node->cumulative_converts_from_parent + (node->type == GAP_CONVERTED_2_NORMAL ? node->length : 0) + node->cumulative_converts_from_left_child;
            genome_nodes.push(node->right_child);
        }
    }
}
