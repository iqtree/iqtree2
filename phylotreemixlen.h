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

    virtual ~PhyloTreeMixlen();

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
     * assign branch length by branch length of a category
     * @param category source category to assign
     */
    void assignMixBranches(int category, Node *node = NULL, Node *dad = NULL);

    /**
     * assign branch length as mean over all branch lengths of categories
     */
    void assignMeanMixBranches(Node *node = NULL, Node *dad = NULL);

    /**
     * copy branch lengths of a tree into a mixture category of this tree
     * @param tree source tree where branch lengths will be copied
     * @param category destination category
     */
    void copyMixBranches(PhyloTree *tree, int category);

    /**
     *  internal function called by printTree to print branch length
     *  @param out output stream
     *  @param length_nei target Neighbor to print
     */
    virtual void printBranchLength(ostream &out, int brtype, bool print_slash, Neighbor *length_nei);


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
     * IMPORTANT: semantic change: this function does not return score anymore, for efficiency purpose
            optimize one branch length by ML
            @param node1 1st end node of the branch
            @param node2 2nd end node of the branch
            @param clearLH true to clear the partial likelihood, otherwise false
            @param maxNRStep maximum number of Newton-Raphson steps
            @return likelihood score
     */
    virtual void optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH = true, int maxNRStep = 100);

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
    
    /** relative rate, used to initialize branch lengths */
    RateHeterogeneity *relative_rate;

    /** true to print mixture branch lengths when calling printTree */
    bool print_mix_brlen;

    /** category tree used internally to optimize parameters with EM algorithm */
    PhyloTree *cat_tree;

};

#endif /* defined(__iqtree__phylotreemixlen__) */
