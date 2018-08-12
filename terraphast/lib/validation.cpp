#include <numeric>
#include <unordered_map>

#include <algorithm>

#include "ranked_bitvector.hpp"
#include "supertree_helpers.hpp"
#include "trees_impl.hpp"
#include "validation.hpp"

namespace terraces {

std::vector<simple_bitvector> tree_bipartitions(const tree& t) {
	std::vector<simple_bitvector> bips(t.size(), {0, {}});
	std::vector<simple_bitvector> subtrees(t.size(), {(t.size() + 1) / 2, {}});
	foreach_postorder(t, [&](index i) {
		auto n = t[i];
		if (is_leaf(n)) {
			subtrees[i].set(n.taxon());
		} else {
			subtrees[i].set_bitwise_or(subtrees[n.lchild()], subtrees[n.rchild()]);
		}
	});
	foreach_preorder(t, [&](index i) {
		auto n = t[i];
		bool at_root = is_root(n);
		bool at_rhs_of_root =
		        !at_root && (is_root(t[n.parent()]) && t[n.parent()].rchild() == i);
		if (!(at_root || at_rhs_of_root)) {
			if (subtrees[i].get(0)) {
				subtrees[i].invert();
			}
			bips[i] = std::move(subtrees[i]);
		}
	});
	bips[t[0].rchild()] = std::move(subtrees[t[0].rchild()]);
	bips[t[0].rchild()].blank();
	bips[0] = std::move(subtrees[0]);
	bips[0].blank();
	std::sort(bips.begin(), bips.end());
	return bips;
}

bool is_isomorphic_unrooted(const tree& fst, const tree& snd) {
	assert(fst.size() == snd.size());

	auto fst_bip = tree_bipartitions(fst);
	auto snd_bip = tree_bipartitions(snd);

	return fst_bip == snd_bip;
}

bool is_isomorphic_rooted_impl(const tree& fst, const tree& snd, index fst_idx, index snd_idx) {
	auto fst_node = fst[fst_idx];
	auto snd_node = snd[snd_idx];
	if (is_leaf(fst_node) != is_leaf(snd_node)) {
		return false;
	} else if (is_leaf(fst_node)) {
		return fst_node.taxon() == snd_node.taxon();
	}

	return (is_isomorphic_rooted_impl(fst, snd, fst_node.lchild(), snd_node.lchild()) &&
	        is_isomorphic_rooted_impl(fst, snd, fst_node.rchild(), snd_node.rchild())) ||
	       (is_isomorphic_rooted_impl(fst, snd, fst_node.lchild(), snd_node.rchild()) &&
	        is_isomorphic_rooted_impl(fst, snd, fst_node.rchild(), snd_node.lchild()));
}

bool is_isomorphic_rooted(const tree& fst, const tree& snd) {
	assert(fst.size() == snd.size());
	return is_isomorphic_rooted_impl(fst, snd, 0, 0);
}

} // namespace terraces
