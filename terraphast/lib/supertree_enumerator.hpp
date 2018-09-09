#ifndef SUPERTREE_ENUMERATOR_HPP
#define SUPERTREE_ENUMERATOR_HPP

#include "bipartitions.hpp"
#include "stack_allocator.hpp"
#include "union_find.hpp"

#include "supertree_helpers.hpp"
#include "utils.hpp"

namespace terraces {

template <typename Callback>
class tree_enumerator {
	using result_type = typename Callback::result_type;

private:
	Callback m_cb;

	utils::free_list m_fl1;
	utils::free_list m_fl2;
	utils::free_list m_fl3;
	index m_fl1_allocsize;
	index m_fl2_allocsize;
	index m_fl3_allocsize;

	const constraints* m_constraints;

	result_type run(const ranked_bitvector& leaves, const bitvector& constraint_occ);
	result_type iterate(bipartitions& bip_it, const bitvector& new_constraint_occ);

	void init_freelists(index leaf_count, index constraint_count);
	utils::stack_allocator<index> leaf_allocator();
	utils::stack_allocator<index> c_occ_allocator();
	utils::stack_allocator<index> union_find_allocator();

public:
	explicit tree_enumerator(Callback cb);
	result_type run(index num_leaves, const constraints& constraints,
	                const std::vector<bool>& root_split);
	result_type run(index num_leaves, const constraints& constraints, index root_leaf);
	result_type run(index num_leaves, const constraints& constraints);
	const Callback& callback() const { return m_cb; }
};

template <typename Callback>
tree_enumerator<Callback>::tree_enumerator(Callback cb)
        : m_cb{std::move(cb)}, m_constraints{nullptr} {}

template <typename Callback>
auto tree_enumerator<Callback>::run(index num_leaves, const constraints& constraints)
        -> result_type {
	init_freelists(num_leaves, constraints.size());
	auto leaves = full_ranked_set(num_leaves, leaf_allocator());
	auto c_occ = full_set(constraints.size(), c_occ_allocator());
	assert(filter_constraints(leaves, c_occ, constraints, c_occ_allocator()) == c_occ);
	m_constraints = &constraints;
	return run(leaves, c_occ);
}

template <typename Callback>
auto tree_enumerator<Callback>::run(index num_leaves, const constraints& constraints,
                                    const std::vector<bool>& root_split) -> result_type {
	init_freelists(num_leaves, constraints.size());
	auto leaves = full_ranked_set(num_leaves, leaf_allocator());
	auto c_occ = full_set(constraints.size(), c_occ_allocator());
	assert(filter_constraints(leaves, c_occ, constraints, c_occ_allocator()) == c_occ);
	assert(root_split.size() == num_leaves);
	// enter the call
	m_cb.enter(leaves);
	// no base cases
	assert(num_leaves > 2);
	// assert(!constraints.empty()); is not necessary, since is returns correct values anyway
	// build bipartition iterator:
	auto sets = union_find::make_bipartition(root_split, union_find_allocator());
	m_constraints = &constraints;
	auto bip_it = bipartitions{leaves, sets, leaf_allocator()};
	return m_cb.exit(iterate(bip_it, c_occ));
}

template <typename Callback>
auto tree_enumerator<Callback>::run(index num_leaves, const constraints& constraints,
                                    index root_leaf) -> result_type {
	std::vector<bool> root_split(num_leaves);
	root_split[root_leaf] = true;
	return run(num_leaves, constraints, root_split);
}

template <typename Callback>
auto tree_enumerator<Callback>::run(const ranked_bitvector& leaves, const bitvector& constraint_occ)
        -> result_type {
	m_cb.enter(leaves);

	// base cases: only a few leaves
	assert(leaves.count() > 0);
	if (leaves.count() == 1) {
		return m_cb.exit(m_cb.base_one_leaf(leaves.first_set()));
	}

	if (leaves.count() == 2) {
		auto fst = leaves.first_set();
		auto snd = leaves.next_set(fst);
		return m_cb.exit(m_cb.base_two_leaves(fst, snd));
	}

	bitvector new_constraint_occ =
	        filter_constraints(leaves, constraint_occ, *m_constraints, c_occ_allocator());
	// base case: no constraints left
	if (new_constraint_occ.empty()) {
		return m_cb.exit(m_cb.base_unconstrained(leaves));
	}

	union_find sets = apply_constraints(leaves, new_constraint_occ, *m_constraints,
	                                    union_find_allocator());
	bipartitions bip_it(leaves, sets, leaf_allocator());

	return m_cb.exit(iterate(bip_it, new_constraint_occ));
}

template <typename Callback>
void tree_enumerator<Callback>::init_freelists(index leaf_count, index constraint_count) {
	m_fl1_allocsize = ranked_bitvector::alloc_size(leaf_count);
	m_fl2_allocsize = ranked_bitvector::alloc_size(constraint_count);
	m_fl3_allocsize = leaf_count;
	m_fl1 = {};
	m_fl2 = {};
	m_fl3 = {};
}

template <typename Callback>
utils::stack_allocator<index> tree_enumerator<Callback>::leaf_allocator() {
	return {m_fl1, m_fl1_allocsize};
}

template <typename Callback>
utils::stack_allocator<index> tree_enumerator<Callback>::c_occ_allocator() {
	return {m_fl2, m_fl2_allocsize};
}

template <typename Callback>
utils::stack_allocator<index> tree_enumerator<Callback>::union_find_allocator() {
	return {m_fl3, m_fl3_allocsize};
}

template <typename Callback>
auto tree_enumerator<Callback>::iterate(bipartitions& bip_it, const bitvector& new_constraint_occ)
        -> result_type {
	if (m_cb.fast_return(bip_it)) {
		return m_cb.fast_return_value(bip_it);
	}

	auto result = m_cb.begin_iteration(bip_it, new_constraint_occ, *m_constraints);
	// iterate over all possible bipartitions
	for (auto bip = bip_it.begin_bip();
	     bip < bip_it.end_bip() && m_cb.continue_iteration(result); ++bip) {
		m_cb.step_iteration(bip_it, bip);
		auto set = bip_it.get_first_set(bip, leaf_allocator());
		m_cb.left_subcall();
		auto left_result = run(set, new_constraint_occ);
		bip_it.flip_set(set);
		m_cb.right_subcall();
		auto right_result = run(set, new_constraint_occ);
		// accumulate result
		result = m_cb.accumulate(result, m_cb.combine(left_result, right_result));
	}
	m_cb.finish_iteration();

	return result;
}

} // namespace terraces

#endif // SUPERTREE_ENUMERATOR_HPP
