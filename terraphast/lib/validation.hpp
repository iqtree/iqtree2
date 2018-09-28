#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include <terraces/parser.hpp>
#include <terraces/trees.hpp>

#include "bitvector.hpp"

namespace terraces {

bool is_isomorphic_unrooted(const tree& fst, const tree& snd);

bool is_isomorphic_rooted(const tree& fst, const tree& snd);

std::vector<simple_bitvector> tree_bipartitions(const tree& t);

} // namespace terraces

#endif // TREE_BIPARTITIONS_HPP
