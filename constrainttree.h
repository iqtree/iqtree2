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

/**
    ConstraintTree used to guide tree search.
    Note that constraint tree may contain only a subset of taxa from a full tree.
*/
class ConstraintTree : public MTree {
public:

    /**
            constructor, read constraint tree from user file
            @param full_tree the full tree with all taxa
            @param constraint_file the name of the constraint tree file
     */
    ConstraintTree(MTree *full_tree, const char *constraint_file);


protected:

    /** full tree that contains all taxa */
    MTree *full_tree;

};

#endif