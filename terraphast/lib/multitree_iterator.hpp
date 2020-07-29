#ifndef MULTITREE_ITERATOR_H
#define MULTITREE_ITERATOR_H

#include "multitree.hpp"
#include "small_bipartition.hpp"

namespace terraces {

struct multitree_iterator_choicepoint {
	const multitree_node* alternatives;
	const multitree_node* current;

	multitree_iterator_choicepoint() : alternatives{nullptr}, current{nullptr} {}

	multitree_iterator_choicepoint(const multitree_node* node) {
		if (node->type == multitree_node_type::alternative_array) {
			const auto aa = node->alternative_array;
			alternatives = aa.num_alternatives() > 1 ? node : nullptr;
			current = aa.begin;
		} else {
			alternatives = nullptr;
			current = node;
		}
	}

	bool has_choices() const { return alternatives != nullptr; }

	bool is_unconstrained() const {
		return current->type == multitree_node_type::base_unconstrained;
	}

	bool is_valid() const {
		return has_choices() && current != alternatives->alternative_array.end;
	}

	bool next() {
		assert(is_valid());
		++current;
		return is_valid();
	}

	void reset() { current = alternatives->alternative_array.begin; }
};

class multitree_iterator {
private:
	terraces::tree m_tree;
	std::vector<multitree_iterator_choicepoint> m_choices;
	std::vector<small_bipartition> m_unconstrained_choices;

	bool init_subtree(index subtree_root);
	bool init_subtree(index subtree_root, index single_leaf);
	bool init_subtree(index subtree_root, multitree_nodes::two_leaves two_leaves);
	bool init_subtree(index subtree_root, multitree_nodes::inner_node inner);
	bool init_subtree(index subtree_root, multitree_nodes::unconstrained unconstrained);
	bool init_subtree_unconstrained(index subtree_root, multitree_nodes::unconstrained data);

	bool next(index root);
	bool next_unconstrained(index root, multitree_nodes::unconstrained unconstrained);
	bool reset(index root);
	bool reset_unconstrained(index root, multitree_nodes::unconstrained unconstrained);

public:
	multitree_iterator(const multitree_node* root);
	bool next();
	const terraces::tree& tree() const { return m_tree; }
};

} // namespace terraces

#endif // MULTITREE_ITERATOR_H
