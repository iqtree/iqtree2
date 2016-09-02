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
#include "alignment.h"

/**
    ConstraintTree used to guide tree search.
    Note that constraint tree may contain only a subset of taxa from a full tree.
*/
class ConstraintTree : public MTree, public SplitIntMap {
public:

    ConstraintTree();

    /**
            initialize constraint tree
            @param constraint_file the name of the constraint tree file
            @param fulltaxname the full list of all taxa names
     */
    void initConstraint(const char *constraint_file, StrVector &fulltaxname);

    /** 
        check if a "partial" split defined by two taxa name sets is compatible with the constraint tree.
        The union of 2 taxa set do not need to comprise all taxa in the constraint tree.
        @param[in] tax1 names of taxa in one side of split
        @param[in] tax2 names of taxa in other side of split
        @return true if the split is compatible with all splits in the constraint tree, false otherwise.
     */ 
    bool isCompatible(StrVector &tax1, StrVector &tax2);

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