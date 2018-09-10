#include "union_find.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cassert>

namespace terraces {

union_find::union_find(index n, utils::stack_allocator<index> a) : m_parent(n, n, a) {
#ifndef NDEBUG
	m_compressed = true;
#endif // NDEBUG
}

index union_find::find(index x) {
	assert(x < m_parent.size());
	index root = x;
	while (!is_representative(root)) {
		root = m_parent[root];
	}
	while (x != root) {
		x = utils::exchange(m_parent[x], root);
	}
	assert(is_representative(root) && root < m_parent.size());
	return root;
}

void union_find::compress() {
	for (index i = 0; i < m_parent.size(); ++i) {
		find(i);
	}
#ifndef NDEBUG
	m_compressed = true;
#endif // NDEBUG
}

void union_find::merge(index x, index y) {
#ifndef NDEBUG
	m_compressed = false;
#endif // NDEBUG
	index i = find(x);
	index j = find(y);
	if (i == j) {
		return;
	}
	if (m_parent[i] < m_parent[j]) {
		// link the smaller group to the larger one
		m_parent[i] = j;
	} else if (m_parent[i] > m_parent[j]) {
		// link the smaller group to the larger one
		m_parent[j] = i;
	} else {
		// equal rank: link arbitrarily and increase rank
		m_parent[j] = i;
		++m_parent[i];
	}
}

union_find union_find::make_bipartition(const std::vector<bool>& split,
                                        utils::stack_allocator<index> alloc) {
	union_find result(split.size(), alloc);
	std::array<index, 2> fst{{none, none}};
	for (index i = 0; i < split.size(); ++i) {
		auto& repr = fst[split[i]];
		repr = repr == none ? i : repr;
		result.merge(repr, i);
	}
	result.compress();
	return result;
}

} // namespace terraces
