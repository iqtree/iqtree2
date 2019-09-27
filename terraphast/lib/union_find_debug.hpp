#ifndef UNION_FIND_DEBUG_HPP
#define UNION_FIND_DEBUG_HPP

#include "union_find.hpp"

namespace terraces {
namespace debug {

class union_find_sets_iterator {
private:
	const union_find& m_sets;
	index m_i;

	index next_rep(index i) const {
		index j = i;
		for (; j < m_sets.size(); ++j) {
			if (m_sets.is_representative(j)) {
				break;
			}
		}
		return j;
	}

public:
	union_find_sets_iterator(const union_find& sets, index i = 0)
	        : m_sets{sets}, m_i{next_rep(i)} {}

	index operator*() const { return m_i; }
	union_find_sets_iterator& operator++() {
		m_i = next_rep(m_i + 1);
		return *this;
	}
	bool operator==(const union_find_sets_iterator& other) const { return m_i == other.m_i; }
	bool operator!=(const union_find_sets_iterator& other) const { return !(*this == other); }
};

struct union_find_iterable_sets {
	const union_find& sets;
	union_find_sets_iterator begin() const { return {sets}; }
	union_find_sets_iterator end() const { return {sets, sets.size()}; }
};

class union_find_set_iterator {
private:
	const union_find& m_sets;
	const index m_rep;
	index m_i;

	index next(index i) const {
		index j = i;
		for (; j < m_sets.size(); ++j) {
			if (m_sets.simple_find(j) == m_rep) {
				break;
			}
		}
		return j;
	}

public:
	union_find_set_iterator(const union_find& sets, index rep, index i = 0)
	        : m_sets{sets}, m_rep{rep}, m_i{next(i)} {}

	index operator*() const { return m_i; }
	union_find_set_iterator& operator++() {
		m_i = next(m_i + 1);
		return *this;
	}
	bool operator==(const union_find_set_iterator& other) const { return m_i == other.m_i; }
	bool operator!=(const union_find_set_iterator& other) const { return !(*this == other); }
};

struct union_find_iterable_set {
	const union_find& sets;
	index representative;
	union_find_set_iterator begin() const { return {sets, representative}; }
	union_find_set_iterator end() const { return {sets, representative, sets.size()}; }
};

} // namespace debug
} // namespace terraces

#endif // UNION_FIND_DEBUG_HPP
