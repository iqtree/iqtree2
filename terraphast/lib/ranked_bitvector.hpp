#ifndef RANKED_BITVECTOR_HPP
#define RANKED_BITVECTOR_HPP

#include "bitvector.hpp"

namespace terraces {

template <typename Alloc>
class basic_ranked_bitvector : public basic_bitvector<Alloc> {
public:
	using value_type = typename basic_bitvector<Alloc>::value_type;

private:
	using base = basic_bitvector<Alloc>;

	std::vector<value_type, Alloc> m_ranks;
	index m_count;
#ifndef NDEBUG
	bool m_ranks_dirty;
#endif // NDEBUG

public:
	basic_ranked_bitvector(index size, Alloc alloc)
	        : basic_bitvector<Alloc>{size, alloc}, m_ranks(base::m_blocks.size(), 0, alloc) {
		base::add_sentinel();
#ifndef NDEBUG
		m_ranks_dirty = true;
#endif // NDEBUG
	}

	/** Sets a bit in the bitvector. */
	void set(index i) {
		base::set(i);
#ifndef NDEBUG
		m_ranks_dirty = true;
#endif // NDEBUG
	}
	/** Clears a bit in the bitvector. */
	void clr(index i) {
		base::clr(i);
#ifndef NDEBUG
		m_ranks_dirty = true;
#endif // NDEBUG
	}
	/** Flips a bit in the bitvector. */
	void flip(index i) {
		base::flip(i);
#ifndef NDEBUG
		m_ranks_dirty = true;
#endif // NDEBUG
	}

	/** Clears all bits in the bitvector. */
	void blank();
	/** Inverts all bits in the bitvector. */
	void invert();
	/** Applies element-wise xor from another bitvector. */
	void bitwise_xor(const basic_bitvector<Alloc>& other);
	/** Sets the values of this bitvector to the bitwise or of two bitvectors. */
	void set_bitwise_or(const basic_bitvector<Alloc>& fst, const basic_bitvector<Alloc>& snd);

	/** Returns the number of set bits. */
	index count() const {
		assert(!m_ranks_dirty);
		return m_count - 1;
	}

	/** Updates the internal data structures after editing the vector. */
	void update_ranks();
	/** Returns the rank of an index, i.e. the number of set bits in the range [0..i) */
	index rank(index i) const;
	index select(index i) const;
};

using ranked_bitvector = basic_ranked_bitvector<utils::stack_allocator<index>>;

template <typename Alloc>
index basic_ranked_bitvector<Alloc>::rank(index i) const {
	assert(!m_ranks_dirty);
	assert(i <= basic_bitvector<Alloc>::m_size);
	index b = bits::block_index(i);
	return m_ranks[b] + bits::partial_popcount(base::m_blocks[b], bits::shift_index(i));
}

template <typename Alloc>
index basic_ranked_bitvector<Alloc>::select(index i) const {
	assert(i <= count());
	auto it = base::begin();
	for (index j = 0; j < i; ++j) {
		++it;
	}
	return *it;
}

template <typename Alloc>
void basic_ranked_bitvector<Alloc>::blank() {
	base::blank();
#ifndef NDEBUG
	m_ranks_dirty = true;
#endif // NDEBUG
}

template <typename Alloc>
void basic_ranked_bitvector<Alloc>::bitwise_xor(const basic_bitvector<Alloc>& other) {
	base::bitwise_xor(other);
#ifndef NDEBUG
	m_ranks_dirty = true;
#endif // NDEBUG
}

template <typename Alloc>
void basic_ranked_bitvector<Alloc>::invert() {
	base::invert();
#ifndef NDEBUG
	m_ranks_dirty = true;
#endif // NDEBUG
}

template <typename Alloc>
void basic_ranked_bitvector<Alloc>::set_bitwise_or(const basic_bitvector<Alloc>& fst,
                                                   const basic_bitvector<Alloc>& snd) {
	base::set_bitwise_or(fst, snd);
#ifndef NDEBUG
	m_ranks_dirty = true;
#endif // NDEBUG
}

template <typename Alloc>
void basic_ranked_bitvector<Alloc>::update_ranks() {
	m_count = 0;
	for (index b = 0; b < base::m_blocks.size(); ++b) {
		m_ranks[b] = m_count;
		m_count += bits::popcount(base::m_blocks[b]);
	}
	assert(m_count > 0);
#ifndef NDEBUG
	m_ranks_dirty = false;
#endif // NDEBUG
}

/** Returns a basic_ranked_bitvector<Alloc> containing size elements. */
template <typename Alloc>
basic_ranked_bitvector<Alloc> full_ranked_set(index size, Alloc a) {
	basic_ranked_bitvector<Alloc> set{size, a};
	set.invert();
	set.update_ranks();
	return set;
}

} // namespace terraces

#endif // RANKED_BITVECTOR_HPP
