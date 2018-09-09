#ifndef BITVECTOR_H
#define BITVECTOR_H

#include <cstdint>
#include <vector>

#include <terraces/trees.hpp>

#include "bits.hpp"
#include "stack_allocator.hpp"

namespace terraces {

template <typename Bitvector>
class bitvector_iterator;

template <typename Allocator>
class basic_bitvector {
public:
	using value_type = index;
	using iterator = bitvector_iterator<Allocator>;

	static index alloc_size(index size) { return size / bits::word_bits + 1; }

protected:
	index m_size;
	std::vector<value_type, Allocator> m_blocks;

	void add_sentinel() {
		// add sentinel bit for iteration
		m_blocks[bits::block_index(m_size)] |= bits::set_mask(m_size);
	}

public:
	/** Initializes a bitvector with given size. */
	basic_bitvector(index size, Allocator alloc)
	        : m_size{size}, m_blocks(alloc_size(size), 0, alloc) {
		add_sentinel();
	}
	/** Sets a bit in the bitvector. */
	void set(index i) {
		assert(i < m_size);
		m_blocks[bits::block_index(i)] |= bits::set_mask(i);
	}
	/** Clears a bit in the bitvector. */
	void clr(index i) {
		assert(i < m_size);
		m_blocks[bits::block_index(i)] &= bits::clear_mask(i);
	}
	/** Flips a bit in the bitvector. */
	void flip(index i) {
		assert(i < m_size);
		m_blocks[bits::block_index(i)] ^= bits::set_mask(i);
	}
	/** Returns a bit from the bitvector. */
	bool get(index i) const {
		assert(i < m_size);
		return ((m_blocks[bits::block_index(i)] >> bits::shift_index(i)) & 1) != 0u;
	}
	/** Returns the size of the bitvector. */
	index size() const { return m_size; }

	/** Returns true if and only if no bit is set. */
	bool empty() const;

	/** Clears all bits in the bitvector. */
	void blank();
	/** Inverts all bits in the bitvector. */
	void invert();
	/** Applies element-wise xor from another bitvector. */
	void bitwise_xor(const basic_bitvector<Allocator>& other);
	/** Sets the values of this bitvector to the bitwise or of two bitvectors. */
	void set_bitwise_or(const basic_bitvector<Allocator>& fst,
	                    const basic_bitvector<Allocator>& snd);

	/** Returns the index of the first set bit or size() if no bit is set. */
	index first_set() const;
	/** Returns the index of the next set bit after the index or size() if no bit is set. */
	index next_set(index i) const;
	/** Returns the index one past the last element. */
	index last_set() const { return m_size; }

	iterator begin() const;
	iterator end() const;

	Allocator get_allocator() const { return m_blocks.get_allocator(); }

	bool operator<(const basic_bitvector<Allocator>& other) const {
		assert(size() == other.size());
		return m_blocks < other.m_blocks;
	}
	bool operator==(const basic_bitvector<Allocator>& other) const {
		assert(size() == other.size());
		return m_blocks == other.m_blocks;
	}
	bool operator!=(const basic_bitvector<Allocator>& other) const { return !(*this == other); }
};

using bitvector = basic_bitvector<utils::stack_allocator<index>>;
using simple_bitvector = basic_bitvector<std::allocator<bitvector::value_type>>;

template <typename Alloc>
class bitvector_iterator {
public:
	using value_type = typename Alloc::value_type;

private:
	const basic_bitvector<Alloc>& m_set;
	index m_index;

public:
	bitvector_iterator(const basic_bitvector<Alloc>& set, index i) : m_set{set}, m_index{i} {}
	bitvector_iterator& operator++() {
		m_index = m_set.next_set(m_index);
		return *this;
	}
	bool operator==(const bitvector_iterator& other) const { return m_index == other.m_index; }
	bool operator!=(const bitvector_iterator& other) const { return !(*this == other); }
	const index& operator*() const { return m_index; }
};

template <typename Alloc>
auto basic_bitvector<Alloc>::begin() const -> iterator {
	return {*this, first_set()};
}

template <typename Alloc>
auto basic_bitvector<Alloc>::end() const -> iterator {
	return {*this, last_set()};
}

template <typename Alloc>
bool basic_bitvector<Alloc>::empty() const {
	for (index b = 0; b < m_blocks.size() - 1; ++b) {
		if (m_blocks[b]) {
			return false;
		}
	}
	return !(m_blocks[m_blocks.size() - 1] & bits::prefix_mask(bits::shift_index(m_size)));
}

template <typename Alloc>
void basic_bitvector<Alloc>::blank() {
	for (auto& el : m_blocks) {
		el = 0;
	}
	add_sentinel();
}

template <typename Alloc>
void basic_bitvector<Alloc>::bitwise_xor(const basic_bitvector<Alloc>& other) {
	assert(size() == other.size());
	for (index b = 0; b < m_blocks.size(); ++b) {
		m_blocks[b] ^= other.m_blocks[b];
	}
	add_sentinel();
}

template <typename Alloc>
void basic_bitvector<Alloc>::invert() {
	for (index b = 0; b < m_blocks.size() - 1; ++b) {
		m_blocks[b] = ~m_blocks[b];
	}
	m_blocks[m_blocks.size() - 1] ^= bits::prefix_mask(bits::shift_index(m_size));
}

template <typename Alloc>
void basic_bitvector<Alloc>::set_bitwise_or(const basic_bitvector<Alloc>& fst,
                                            const basic_bitvector<Alloc>& snd) {
	assert(size() == fst.size() && size() == snd.size());
	for (index b = 0; b < m_blocks.size(); ++b) {
		m_blocks[b] = fst.m_blocks[b] | snd.m_blocks[b];
	}
}

template <typename Alloc>
index basic_bitvector<Alloc>::first_set() const {
	index b = 0;
	while (!bits::has_next_bit0(m_blocks[b])) {
		++b;
	}
	return bits::next_bit0(m_blocks[b], bits::base_index(b));
}

template <typename Alloc>
index basic_bitvector<Alloc>::next_set(index i) const {
	++i;
	index b = bits::block_index(i);
	if (bits::has_next_bit(m_blocks[b], i)) {
		// the next bit is in the current block
		return bits::next_bit(m_blocks[b], i);
	}
	// the next bit is in a far-away block
	do {
		++b;
	} while (!bits::has_next_bit0(m_blocks[b]));
	return bits::next_bit0(m_blocks[b], bits::base_index(b));
}

/** Returns a bitvector containing size elements. */
template <typename Alloc>
basic_bitvector<Alloc> full_set(index size, Alloc a) {
	basic_bitvector<Alloc> set{size, a};
	set.invert();
	return set;
}

} // namespace terraces

#endif // BITVECTOR_H
