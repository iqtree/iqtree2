//
//  phylotreemixlen.h
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#ifndef __iqtree__phylotreemixlen__
#define __iqtree__phylotreemixlen__

#include <stdio.h>
#include "iqtree.h"


/**
    Phylogenetic tree with mixture of branch lengths
    Started within joint project with Stephen Crotty
*/
class PhyloTreeMixlen : public IQTree {

public:

    /**
            default constructor
     */
    PhyloTreeMixlen();

    PhyloTreeMixlen(Alignment *aln, int mixlen);


    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name
            @return a new node
     */
    virtual Node* newNode(int node_id = -1, const char* node_name = NULL);

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name issued by an interger
            @return a new node
     */
    virtual Node* newNode(int node_id, int node_name);

    /**
        @return true if this is a tree with mixture branch lengths, default: false
    */
    virtual bool isMixlen() { return true; }

    /**
        set number of mixture branch lengths
    */
    void setMixlen(int mixlen);

    /**
     * assign branch lengths of a category to the overall branch length
     * @param category source category to assign
     */
    void assignMixBranches(int category, Node *node = NULL, Node *dad = NULL);

    void copyMixBranches(PhyloTree *tree, int category);

    /**
            print the tree to the output file in newick format
            @param out the output file.
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param brtype type of branch to print
            @return ID of the taxon with smallest ID
     */
//    virtual int printTree(ostream &out, int brtype, Node *node, Node *dad = NULL);

    /**
        initialize mixture branch lengths
    */
    void initializeMixBranches(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
            optimize all branch lengths of the tree
            @param iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);

    /** number of mixture categories */
    int mixlen;

protected:

    /** current category, for printTree */
    int cur_mixture;
};

#endif /* defined(__iqtree__phylotreemixlen__) */
