#include "bipartitions.hpp"

#include <cassert>
#include <ostream>

#include <terraces/errors.hpp>

#include "bits.hpp"
#include "utils.hpp"

namespace terraces {

bipartitions::bipartitions(const ranked_bitvector& leaves, const union_find& sets,
                           utils::stack_allocator<index> a)
        : m_alloc{a}, m_leaves{leaves}, m_sets{sets}, m_set_rep{find_set_reps()},
          m_end{(index(1) << (m_set_rep.count() - 1))} {
	utils::ensure<tree_count_overflow_error>(m_set_rep.count() < bits::word_bits,
	                                         "Huge terrace encountered");
}

bool bipartitions::in_left_partition(index bip, index i) const {
	return (bip & (index(1) << ((i - 1) % bits::word_bits))) != 0;
}

ranked_bitvector bipartitions::find_set_reps() const {
	ranked_bitvector set_rep(m_leaves.count(), m_alloc);
	for (index i = 0; i < m_leaves.count(); ++i) {
		set_rep.set(m_sets.simple_find(i));
	}
	set_rep.update_ranks();
	return set_rep;
}

ranked_bitvector bipartitions::get_first_set(index bip, utils::stack_allocator<index> alloc) const {
	ranked_bitvector subleaves(m_leaves.size(), alloc);
	index ii = 0;
	for (auto i = m_leaves.first_set(); i < m_leaves.last_set(); i = m_leaves.next_set(i)) {
		if (in_left_partition(bip, m_set_rep.rank(m_sets.simple_find(ii)))) {
			subleaves.set(i);
		}
		++ii;
	}
	subleaves.update_ranks();
	return subleaves;
}

void bipartitions::flip_set(ranked_bitvector& set) const {
	set.bitwise_xor(m_leaves);
	set.update_ranks();
}

std::pair<ranked_bitvector, ranked_bitvector>
bipartitions::get_both_sets(index bip, utils::stack_allocator<index> alloc) const {
	auto first_set = get_first_set(bip, alloc);
	auto second_set = first_set;
	flip_set(second_set);
	return {std::move(first_set), std::move(second_set)};
}

std::pair<ranked_bitvector, ranked_bitvector> bipartitions::get_both_sets_unsafe(index bip) const {
	return get_both_sets(bip, m_alloc);
}

} // namespace terraces
