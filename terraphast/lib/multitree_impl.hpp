#ifndef MULTITREE_IMPL_HPP
#define MULTITREE_IMPL_HPP

#include "multitree.hpp"

#include <memory>

namespace terraces {
namespace multitree_impl {

template <typename T>
struct array_deleter {
	void operator()(T* ptr) { delete[] ptr; }
};

template <typename T>
std::unique_ptr<T[], array_deleter<T>> make_unique_array(std::size_t size) {
	// this may leak memory if allocation or construction of a T fails
	return std::unique_ptr<T[], array_deleter<T>>{new T[size]};
}

template <typename T>
class storage_block {
private:
	std::unique_ptr<T[], array_deleter<T>> begin;
	index size;
	index max_size;

public:
	storage_block(index max_size)
	        : begin{make_unique_array<T>(max_size)}, size{0}, max_size{max_size} {}

	bool has_space(index required = 1) { return size + required <= max_size; }

	T* get() {
		assert(has_space());
		return begin.get() + (size++);
	}

	T* get_range(index required) {
		assert(has_space(required));
		auto result = begin.get() + size;
		size += required;
		return result;
	}
};

template <typename T>
class storage_blocks {
private:
	std::vector<storage_block<T>> m_blocks;
	index m_block_size;

public:
	storage_blocks(index block_size = 1024) : m_blocks{}, m_block_size{block_size} {
		m_blocks.emplace_back(m_block_size);
	}
	storage_blocks(const storage_blocks<T>& other) : storage_blocks{other.m_block_size} {}
	storage_blocks(storage_blocks<T>&& other) = default;
	storage_blocks<T>& operator=(const storage_blocks<T>& other) {
		m_block_size = other.m_block_size;
		return *this;
	}
	storage_blocks<T>& operator=(storage_blocks<T>&& other) = default;

	T* get() {
		if (!m_blocks.back().has_space()) {
			m_blocks.emplace_back(m_block_size);
		}
		return m_blocks.back().get();
	}

	T* get_range(index required) {
		if (!m_blocks.back().has_space(required)) {
			m_blocks.emplace_back(required);
			auto result = m_blocks.back().get_range(required);
			auto last_it = --m_blocks.end();
			auto prev_it = --(--m_blocks.end());
			std::iter_swap(
			        last_it,
			        prev_it); // TODO this might lead to some bad worst-case behaviour
			return result;
		}
		return m_blocks.back().get_range(required);
	}
};

inline multitree_node* make_single_leaf(multitree_node* n, index i) {
	n->type = multitree_node_type::base_single_leaf;
	n->single_leaf = i;
	n->num_leaves = 1;
	n->num_trees = 1;
	return n;
}

inline multitree_node* make_two_leaves(multitree_node* n, index i, index j) {
	n->type = multitree_node_type::base_two_leaves;
	n->two_leaves = {i, j};
	n->num_leaves = 2;
	n->num_trees = 1;
	return n;
}

inline multitree_node* make_unconstrained(multitree_node* n, std::pair<index*, index*> range) {
	auto begin = range.first;
	auto end = range.second;
	n->type = multitree_node_type::base_unconstrained;
	n->unconstrained = {begin, end};
	n->num_leaves = (index)(end - begin);
	n->num_trees = count_unrooted_trees<index>(n->num_leaves);
	return n;
}

inline multitree_node* make_inner_node(multitree_node* n, multitree_node* left,
                                       multitree_node* right) {
	n->type = multitree_node_type::inner_node;
	n->inner_node = {left, right};
	n->num_leaves = left->num_leaves + right->num_leaves;
	n->num_trees = left->num_trees * right->num_trees;
	return n;
}

inline multitree_node* make_alternative_array(multitree_node* n, multitree_node* begin,
                                              index leaves) {
	n->type = multitree_node_type::alternative_array;
	n->alternative_array = {begin, begin};
	n->num_leaves = leaves;
	n->num_trees = 0;
	return n;
}

inline multitree_node* make_unexplored(multitree_node* n, std::pair<index*, index*> range) {
	auto begin = range.first;
	auto end = range.second;
	n->type = multitree_node_type::unexplored;
	n->unexplored = {begin, end};
	n->num_leaves = (index)(end - begin);
	n->num_trees = 0;
	return n;
}

} // namespace multitree_impl
} // namespace terraces

#endif // MULTITREE_IMPL_HPP
