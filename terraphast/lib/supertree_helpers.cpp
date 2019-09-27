#include "supertree_helpers.hpp"
#include "trees_impl.hpp"

namespace terraces {

bitvector filter_constraints(const ranked_bitvector& leaves, const bitvector& c_occ,
                             const constraints& c, utils::stack_allocator<index> a) {
	bitvector result{c_occ.size(), a};
	for (auto c_i = c_occ.first_set(); c_i < c_occ.last_set(); c_i = c_occ.next_set(c_i)) {
		if (leaves.get(c[c_i].left) && leaves.get(c[c_i].shared) &&
		    leaves.get(c[c_i].right)) {
			result.set(c_i);
		}
	}
	return result;
}

union_find apply_constraints(const ranked_bitvector& leaves, const bitvector& c_occ,
                             const constraints& c, utils::stack_allocator<index> a) {
	auto sets = union_find(leaves.count(), a);
	for (auto c_i = c_occ.first_set(); c_i < c_occ.last_set(); c_i = c_occ.next_set(c_i)) {
		auto& cons = c[c_i];
		sets.merge(leaves.rank(cons.left), leaves.rank(cons.shared));
	}
	sets.compress();
	return sets;
}

constraints map_constraints(const ranked_bitvector& leaves, const constraints& cs) {
	auto result = cs;
	for (auto& c : result) {
		c.left = leaves.rank(c.left);
		c.shared = leaves.rank(c.shared);
		c.right = leaves.rank(c.right);
	}
	return result;
}

ranked_bitvector leaf_occ(const tree& tree, utils::stack_allocator<index> a) {
	ranked_bitvector leaves{tree.size(), a};
	for (index i = 0; i < tree.size(); i++) {
		if (is_leaf(tree[i])) {
			leaves.set(i);
		}
	}
	leaves.update_ranks();
	return leaves;
}

} // namespace terraces
