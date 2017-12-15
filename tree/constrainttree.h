//
// C++ Interface: phylotree.h
//
// Description:
//
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CONSTRAINTTREE_H
#define CONSTRAINTTREE_H

#include "mtree.h"
#include "alignment/alignment.h"

struct NNIMove;

/**
    ConstraintTree used to guide tree search.
    Note that constraint tree may contain only a subset of taxa from a full tree.
*/
class ConstraintTree : public MTree, public SplitIntMap {
public:

    /** constructor */
    ConstraintTree();

    /** destructor */
    virtual ~ConstraintTree();

    /**
        internal function to initialize splits from tree structure
    */
    void initFromTree();

    /**
            initialize constraint tree
            @param constraint_file the name of the constraint tree file
            @param fulltaxname the full list of all taxa names
     */
    void readConstraint(const char *constraint_file, StrVector &fulltaxname);

    /**
        initialize from another constraint tree
        @param src source constraint tree
    */
    void readConstraint(MTree &src_tree);

	/** remove some taxa from the tree
	 * @param taxa_names names of taxa that will be removed
     * @return number of taxa actually removed
	 */
	virtual int removeTaxa(StrVector &taxa_names);

    /** 
        check if a "partial" split defined by two taxa name sets is compatible with the constraint tree.
        The union of 2 taxa set do not need to comprise all taxa in the constraint tree.
        @param[in] tax1 names of taxa in one side of split
        @param[in] tax2 names of taxa in other side of split
        @return true if the split is compatible with all splits in the constraint tree, false otherwise.
     */ 
    bool isCompatible(StrVector &tax1, StrVector &tax2);

    /**
        check if a branch defined by two nodes in any tree is compatible or not
        @param node1 one end node of the branch
        @param node2 the other end node of the same branch
        @return TRUE if the branch is compatible, FALSE otherwise 
    */
    bool isCompatible(Node *node1, Node *node2);

    /**
        @param tree input tree
        @return TRUE if input tree is compatible with constraint, FALSE otherwise
    */
    bool isCompatible (MTree *tree);


    /** 
        check if an NNI is compatible with the constraint tree or not
        @param nni an NNIMove
        @return TRUE if the NNI is compatible, FALSE otherwise
    */
    bool isCompatible(NNIMove &nni);

    /**
        @param taxname taxon name to search for
        @return TRUE if constraint tree has a taxon, FALSE otherwise
    */
    bool hasTaxon(string &taxname) {
        return taxname_index.find(taxname) != taxname_index.end();
    }

protected:

    /* map from taxon name to its index, used for quick taxon name search */
    StringIntMap taxname_index;

};

#endif
