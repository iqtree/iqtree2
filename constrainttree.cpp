//
// C++ Implementation: constrainttree.cpp
//
// Description: ConstraintTree class used to guide tree search
//
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "constrainttree.h"


ConstraintTree::ConstraintTree(MTree *full_tree, const char *constraint_file) : MTree() {
    bool is_rooted = false;
    init(constraint_file, is_rooted);
    if (is_rooted != full_tree->rooted)
        outError("Constraint tree and full tree do not agree on rooting");
    this->full_tree = full_tree;
    // check that constraint tree has a subset of taxa
    
}


