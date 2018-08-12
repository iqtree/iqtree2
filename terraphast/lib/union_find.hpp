#ifndef TERRACES_UNION_FIND_HPP
#define TERRACES_UNION_FIND_HPP

#include <algorithm>

#include <terraces/trees.hpp>

#include "stack_allocator.hpp"

namespace terraces {

class union_find {
public:
	using value_type = index;

private:
	std::vector<index, utils::stack_allocator<index>> m_parent;
#ifndef NDEBUG
	bool m_compressed;
#endif // NDEBUG

public:
	union_find(index, utils::stack_allocator<index> a);
	index find(index);
	index simple_find(index x) const {
		assert(m_compressed);
		return is_representative(x) ? x : m_parent[x];
	}
	index size() const { return m_parent.size(); }
	void compress();
	void merge(index, index);
	bool is_representative(index x) const { return m_parent[x] >= m_parent.size(); }

	static union_find make_bipartition(const std::vector<bool>& split,
	                                   utils::stack_allocator<index> alloc);
};

} // namespace terraces

#endif // TERRACES_UNION_FIND_HPP
