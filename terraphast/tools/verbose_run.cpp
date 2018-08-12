#include "../lib/bitvector.hpp"
#include "../lib/constraints_impl.hpp"
#include "../lib/io_utils.hpp"
#include "../lib/multitree_iterator.hpp"
#include "../lib/subtree_extraction_impl.hpp"
#include "../lib/supertree_enumerator.hpp"
#include "../lib/supertree_variants_debug.hpp"
#include "../lib/supertree_variants_multitree.hpp"
#include "../lib/validation.hpp"
#include <iostream>
#include <terraces/constraints.hpp>
#include <terraces/parser.hpp>
#include <terraces/rooting.hpp>

enum class mode { count, enumerate };

using namespace terraces;

bool is_comprehensive(const tree& t, const bitmatrix& matrix) {
	const auto node_occ = compute_node_occ(t, matrix).first;
	for (terraces::index i = 0; i < matrix.cols(); ++i) {
		if (induced_lca(t, node_occ, i) != 0) {
			return false;
		}
	}
	return true;
}

void count_and_check_trees(const multitree_node* multitree, const bitmatrix& matrix,
                           const std::vector<tree>& ref_trees, const name_map& names) {
	std::cerr << "Enumerating all trees" << std::endl;

	// store all supertrees
	std::vector<tree> supertrees;
	multitree_iterator mit(multitree);
	do {
		supertrees.push_back(mit.tree());
	} while (mit.next());
	if (big_integer{supertrees.size()} != multitree->num_trees) {
		std::cerr << "PANIC! Inconsistent tree count" << std::endl;
	}

	// check for supertree property
	for (const auto& supertree : supertrees) {
		auto trees = subtrees(supertree, matrix);
		for (terraces::index i = 0; i < trees.size(); ++i) {
			if (!is_isomorphic_rooted(ref_trees[i], trees[i])) {
				std::cerr << "PANIC! Invalid supertree:\nPartition " << i
				          << "\nSupertree: " << as_newick(supertree, names)
				          << "\nExpected: " << as_newick(ref_trees[i], names)
				          << "\nReceived: " << as_newick(trees[i], names)
				          << std::endl;
			}
		}
	}

	// check for isomorphic trees
	std::vector<std::pair<std::vector<simple_bitvector>, terraces::index>> tree_bips;
	for (terraces::index i = 0; i < supertrees.size(); ++i) {
		const auto& supertree = supertrees[i];
		tree_bips.emplace_back(tree_bipartitions(supertree), i);
	}
	std::sort(tree_bips.begin(), tree_bips.end());
	int dupes = 0;
	for (terraces::index i = 1; i < supertrees.size(); ++i) {
		if (tree_bips[i - 1].first == tree_bips[i].first) {
			++dupes;
			std::cerr << "Duplicate tree found:\n"
			          << as_newick(supertrees[tree_bips[i].second], names) << '\n'
			          << as_newick(supertrees[tree_bips[i - 1].second], names)
			          << std::endl;
		}
	}
	std::cout << (supertrees.size() - dupes) << "(" << supertrees.size() << ")";
}

int main(int argc, char** argv) {
	if (argc < 5) {
		std::cerr << "Usage: " << argv[0]
		          << " [treefile] [datafile] [root_node] [use_root_split]\n";
		return -1;
	}

	terraces::index root_node = std::stoul(argv[3]);
	bool use_root_split = bool(std::stoi(argv[4]));

	try {
		auto data_stream = utils::open_ifstream(argv[2]);
		auto data = parse_bitmatrix(data_stream);
		auto tree = parse_nwk(utils::read_file_full(argv[1]), data.indices);
		auto num_leaves = data.matrix.rows();
		name_map number_names;
		for (terraces::index i = 0; i < num_leaves; ++i) {
			number_names.push_back(std::to_string(i));
		}

		std::cout << "Rooting at node " << root_node << "\n";
		auto rerooted_tree = reroot_at_node(tree, root_node);
		std::cout << as_newick(rerooted_tree, number_names) << std::endl;
		auto subtrees = terraces::subtrees(rerooted_tree, data.matrix);
		std::cout << "Per-partition trees\n";
		for (const auto& subtree : subtrees) {
			std::cout << as_newick(subtree, number_names) << "\n";
		}
		auto constraints = compute_constraints(subtrees);
		std::cout << "Constraints\n";
		for (auto c : constraints) {
			std::cout << utils::named_output<terraces::constraints,
			                                 name_map>{c, &number_names}
			          << "\n";
		}
		std::cout << "Removed " << deduplicate_constraints(constraints)
		          << " duplicate constraints\n";
		auto root_split = terraces::root_split(rerooted_tree, num_leaves);
		std::cout << "Root split\n";
		for (auto b : root_split) {
			std::cout << b;
		}
		std::cout << "\n";
		using callback = debug::variants::logging_decorator<variants::multitree_callback>;
		tree_enumerator<callback> e{callback{{}, std::cout, number_names}};
		multitree_node* result;
		if (use_root_split) {
			result = e.run(num_leaves, constraints, root_split);
		} else {
			result = e.run(num_leaves, constraints);
		}
		multitree_iterator mit(result);
		do {
			std::cout << as_newick(mit.tree(), number_names) << "\n";
		} while (mit.next());
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
	}
}
