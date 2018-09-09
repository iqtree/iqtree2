#ifndef MULTITREE_HPP
#define MULTITREE_HPP

#include "trees_impl.hpp"

#include <terraces/bigint.hpp>

namespace terraces {

enum class multitree_node_type {
	base_single_leaf,
	base_two_leaves,
	base_unconstrained,
	inner_node,
	alternative_array,
	unexplored,
};

struct multitree_node;

namespace multitree_nodes {
struct two_leaves {
	index left_leaf;
	index right_leaf;
};
struct unconstrained {
	index* begin;
	index* end;
	index num_leaves() const;
};
struct inner_node {
	multitree_node* left;
	multitree_node* right;
};
struct alternative_array {
	multitree_node* begin;
	multitree_node* end;
	index num_alternatives() const;
};
struct unexplored {
	index* begin;
	index* end;
	index num_leaves() const;
};
}

struct multitree_node {
	multitree_node_type type;
	index num_leaves;
	big_integer num_trees;
	union {
		index single_leaf;
		multitree_nodes::two_leaves two_leaves;
		multitree_nodes::unconstrained unconstrained;
		multitree_nodes::inner_node inner_node;
		multitree_nodes::alternative_array alternative_array;
		multitree_nodes::unexplored unexplored;
	};
};

struct newick_multitree_t {
	const multitree_node* root;
	const name_map* names;
};

std::ostream& operator<<(std::ostream& stream, newick_multitree_t tree);

inline newick_multitree_t as_newick(const multitree_node* root, const name_map& names) {
	return {root, &names};
}

inline index multitree_nodes::alternative_array::num_alternatives() const {
	return (index)(end - begin);
}
inline index multitree_nodes::unexplored::num_leaves() const { return (index)(end - begin); }
inline index multitree_nodes::unconstrained::num_leaves() const { return (index)(end - begin); }

} // namespace terraces

#endif // MULTITREE_HPP
