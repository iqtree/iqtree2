//
//  phylosupertreeunlinked.h
//  tree
//
//  Created by Minh Bui on 2/5/18.
//

#ifndef phylosupertreeunlinked_h
#define phylosupertreeunlinked_h

#include "phylosupertree.h"

/**
    Super-tree with separate partition tree topologies
 */
class PhyloSuperTreeUnlinked : public PhyloSuperTree {
public:
    /**
     constructors
     */
    PhyloSuperTreeUnlinked(SuperAlignment *alignment);

    virtual bool isSuperTreeUnlinked() {
        return true;
    }

    /**
     read the tree from the ifstream in newick format
     @param in the input stream.
     @param is_rooted (IN/OUT) true if tree is rooted
     */
    virtual void readTree(istream &in, bool &is_rooted);
    
    
    /**
     Set the alignment, important to compute parsimony or likelihood score
     Assing taxa ids according to their position in the alignment
     @param alignment associated alignment
     */
    virtual void setAlignment(Alignment *alignment);

    /**
     * setup all necessary parameters  (declared as virtual needed for phylosupertree)
     */
    virtual void initSettings(Params& params);

    /**
     create sub-trees T|Y_1,...,T|Y_k of the current super-tree T
     and map F={f_1,...,f_k} the edges of supertree T to edges of subtrees T|Y_i
     */
    virtual void mapTrees();

    /**
     * FAST VERSION: compute parsimony tree by step-wise addition
     * @param out_prefix prefix for .parstree file
     * @param alignment input alignment
     * @param rand_stream random stream
     * @return parsimony score
     */
    virtual int computeParsimonyTree(const char *out_prefix, Alignment *alignment, int *rand_stream);

    /**
     * Assign branch lengths for branch that has no or negative length
     * With single model branch lengths are assigned using parsimony. With partition model
     * branch lengths are assigned randomly
     * @param force_change if true then force fixing also positive branch lengths
     * @return number of branches fixed
     */
    virtual int wrapperFixNegativeBranch(bool force_change);

    /** @return true if tree is bifurcating, false otherwise */
    virtual bool isBifurcating(Node *node = NULL, Node *dad = NULL);

    /**
     Read the tree saved with Taxon IDs and branch lengths.
     @param tree_string tree string to read from
     @param updatePLL if true, tree is read into PLL
     */
    virtual void readTreeString(const string &tree_string);

    /*
     * Return the tree string contining taxon IDs and branch lengths
     * @return
     * @param format (WT_TAXON_ID, WT_BR_LEN, ...)
     * @return the tree string with the specified format
     */
    virtual string getTreeString();

    /**
     save object into the checkpoint
     */
    virtual void saveCheckpoint();
    
    /**
     restore object from the checkpoint
     */
    virtual void restoreCheckpoint();
    
    /**
     * save branch lengths into a vector
     */
    virtual void saveBranchLengths(DoubleVector &lenvec, int startid = 0, PhyloNode *node = NULL, PhyloNode *dad = NULL);
    /**
     * restore branch lengths from a vector previously called with saveBranchLengths
     */
    virtual void restoreBranchLengths(DoubleVector &lenvec, int startid = 0, PhyloNode *node = NULL, PhyloNode *dad = NULL);
    
    /** set the root by name
     @param my_root root node name
     @param multi_taxa TRUE if my_root is a comma-separated list of nodes
     */
    virtual void setRootNode(const char *my_root, bool multi_taxa = false);
    
    /**
     compute the weighted average of branch lengths over partitions
     */
    virtual void computeBranchLengths();

    /**
     print the tree to the output file in newick format
     @param out the output stream.
     @param brtype type of branch to print
     */
    virtual void printTree(ostream & out, int brtype = WT_BR_LEN);

    /**
     print tree to .treefile
     @param params program parameters, field root is taken
     */
    virtual void printResultTree(string suffix = "");

    /**
     @return sum of all branch lengths
     @param node the starting node, NULL to start from the root
     @param dad dad of the node, used to direct the search
     */
    virtual double treeLength(Node *node = NULL, Node *dad = NULL);

    
    /**
     @return sum length of all internal branches
     @param node the starting node, NULL to start from the root
     @param dad dad of the node, used to direct the search
     */
    virtual double treeLengthInternal(double epsilon, Node *node = NULL, Node *dad = NULL);

    /**
     *         @brief Perform NNI search on the current tree topology
     *         @return <number_of_NNIs, number_of_NNI_steps>
     *         This function will automatically use the selected kernel (either PLL or IQ-TREE)
     */
    virtual pair<int, int> doNNISearch(bool write_info = false);

    /**
     perform tree search
     @return best likelihood found
     */
    virtual double doTreeSearch();

    /** summarize bootstrap trees */
    virtual void summarizeBootstrap(Params &params);

    /**
     write .ufboot trees file
     */
    virtual void writeUFBootTrees(Params &params);

    /**
     Test all branches of the tree with aLRT SH-like interpretation
     */
    virtual int testAllBranches(int threshold, double best_score, double *pattern_lh,
                                int reps, int lbp_reps, bool aLRT_test, bool aBayes_test,
                                PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     test the best number of threads
     */
    virtual int testNumThreads();

};

#endif /* phylosupertreeunlinked_h */
