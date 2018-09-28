#ifndef SUBTREE_EXTRACTION_H
#define SUBTREE_EXTRACTION_H

#include <vector>

#include "bitmatrix.hpp"
#include "trees.hpp"

namespace terraces {

/**
 * Extracts all (pruned) subtrees from the given input tree and missing data matrix.
 * For every column of the missing data matrix, a corresponding tree is constructed,
 * containing only leaves for which this matrix contains a 1.
 */
std::vector<tree> subtrees(const tree& t, const bitmatrix& occ);

} // namespace terraces

#endif // SUBTREE_EXTRACTION_H
