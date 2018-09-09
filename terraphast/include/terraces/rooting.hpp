#ifndef TERRACES_ROOTING_HPP
#define TERRACES_ROOTING_HPP

#include "trees.hpp"

namespace terraces {

/**
 * Returns the root split of the given tree.
 * The root split is a bitvector of size #leaves that is 1
 * for all leaves in the right subtree of the root node.
 */
std::vector<bool> root_split(const tree& t, index num_leaves);

/**
 * Re-roots the given tree in-place.
 * The node at the given \p node_idx will be placed below the new root node.
 */
tree reroot_at_node(const tree& t, index node_idx);

/**
 * Re-roots the given tree in-place.
 * The leaf corresponding to the given \p comp_taxon will be placed to the right of the root.
 * child of our new root, with the rest of the tree being the left subtree.
 */
void reroot_at_taxon_inplace(tree& t, index comp_taxon);

} // namespace terraces

#endif
