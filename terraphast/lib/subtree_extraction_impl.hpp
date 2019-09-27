#ifndef SUBTREE_EXTRACTION_IMPL_HPP
#define SUBTREE_EXTRACTION_IMPL_HPP

#include <terraces/subtree_extraction.hpp>

namespace terraces {

std::pair<bitmatrix, std::vector<index>> compute_node_occ(const tree& t, const bitmatrix& occ);

index induced_lca(const tree& t, const bitmatrix& node_occ, index column);

tree subtree(const tree& t, const bitmatrix& node_occ,
             const std::vector<index>& num_leaves_per_site, index site);

} // namespace terraces

#endif // SUBTREE_EXTRACTION_IMPL_HPP
