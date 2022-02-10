//
//  genometree.h
//  tree
//
//  Created by Nhan Ly-Trong on 09/02/2022.
//
#ifndef GENOMETREE_H
#define GENOMETREE_H

#include "genomenode.h"

/**
A Genome Tree to present a genome by genome entry (each is a set of sites)
 */
class GenomeTree {
private:
    /**
        find a node that contains a given position
     */
    GenomeNode* findNodeByPos(GenomeNode* node, Insertion insertion);
    
    /**
        insert gaps into a "all-gap" node
     */
    void insertGapsIntoGaps(GenomeNode* node, int length);
    
    /**
        insert gaps into a "normal" (all sites from the original genome) node
     */
    void insertGapsIntoNormalNode(GenomeNode* node, int pos, int length);
    
    /**
        recursively traverse the tree to update pos_new after insert gaps into a node
     */
    void updatePosNew(GenomeNode* node, int pos, int length);
    
    /**
        export sites from a genome node and the original genome
     */
    void exportSiteFromGenomeNode(GenomeNode* node, vector<short int> ori_seq, vector<short int> &new_seq);
    
    /**
        export readable characters from a genome node and the original genome
     */
    void exportReadableCharactersFromGenomeNode(GenomeNode* node, vector<short int> ori_seq, string &output, int num_sites_per_state, vector<string> state_mapping);
    
    
public:
    /**
        starting pos in the original genome
     */
    GenomeNode* root;

    /**
        constructor
     */
    GenomeTree();
    
    /**
        init a root genome node
     */
    GenomeTree(int length);
    
    /**
        deconstructor
     */
    ~GenomeTree();
    
    /**
        update genome tree from a vector of insertions
     */
    void updateTree(vector<Insertion> insertions);
    
    /**
        export new genome from original genome and genome tree
     */
    vector<short int> exportNewGenome(vector<short int> ori_seq, int seq_length, int UNKOWN_STATE);
    
    /**
     export readable characters (for writing to file) from original genome and genome tree
     */
    string exportReadableCharacters(vector<short int> ori_seq, int seq_length, int num_sites_per_state, vector<string> state_mapping);

};
#endif
