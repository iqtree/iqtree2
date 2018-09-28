#include "multitree.hpp"
#include "io_utils.hpp"

namespace terraces {

struct index_array_view {
	index* _begin;
	index* _end;
	index* begin() const { return _begin; }
	index* end() const { return _end; }
};

std::ostream& print_multitree_node(std::ostream& stream, const multitree_node* node,
                                   const name_map& names) {
	switch (node->type) {
	case multitree_node_type::base_single_leaf:
		return stream << names[node->single_leaf];
	case multitree_node_type::base_two_leaves: {
		auto& tl = node->two_leaves;
		return stream << '(' << names[tl.left_leaf] << ',' << names[tl.right_leaf] << ')';
	}
	case multitree_node_type::base_unconstrained: {
		auto& u = node->unconstrained;
		return stream << '{' << utils::as_comma_separated_output(
		                                index_array_view{u.begin, u.end}, names)
		              << '}';
	}
	case multitree_node_type::inner_node: {
		auto& in = node->inner_node;
		stream << '(';
		print_multitree_node(stream, in.left, names);
		stream << ',';
		print_multitree_node(stream, in.right, names);
		stream << ')';
		return stream;
	}
	case multitree_node_type::alternative_array: {
		auto& aa = node->alternative_array;
		for (auto it = aa.begin; it != aa.end; ++it) {
			if (it != aa.begin) {
				stream << '|';
			}
			print_multitree_node(stream, it, names);
		}
		return stream;
	}
	case multitree_node_type::unexplored: {
		auto& u = node->unexplored;
		return stream << '[' << utils::as_comma_separated_output(
		                                index_array_view{u.begin, u.end}, names)
		              << ']';
	}
	default:
		assert(false);
		return stream;
	}
}

std::ostream& operator<<(std::ostream& stream, newick_multitree_t tree) {
	auto node = tree.root;
	auto& names = *tree.names;
	return print_multitree_node(stream, node, names);
}

} // namespace terraces
