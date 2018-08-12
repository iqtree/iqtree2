#include "../lib/bitvector.hpp"
#include "../lib/io_utils.hpp"
#include "../lib/multitree_iterator.hpp"
#include "../lib/subtree_extraction_impl.hpp"
#include "../lib/supertree_enumerator.hpp"
#include "../lib/supertree_variants.hpp"
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
	std::string usage = argc < 3 ? "" : argv[1];
	mode m;
	if (usage == "count") {
		m = mode::count;
	} else if (usage == "enumerate") {
		m = mode::enumerate;
	} else {
		std::cerr << "Usage: " << argv[0] << " {count|enumerate} [treefile] [datafile]\n";
		return -1;
	}

	try {
		auto data_stream = utils::open_ifstream(argv[3]);
		auto data = parse_bitmatrix(data_stream);
		auto tree = parse_nwk(utils::read_file_full(argv[2]), data.indices);
		auto num_leaves = data.matrix.rows();

		for (terraces::index i = 1; i < tree.size(); ++i) {
			std::cerr << "Rooting at node " << i << ":\n";
			auto rerooted_tree = reroot_at_node(tree, i);
			std::cerr << as_newick(rerooted_tree, data.names) << std::endl;
			auto subtrees = terraces::subtrees(rerooted_tree, data.matrix);
			auto constraints = compute_constraints(subtrees);
			std::cerr << "Extracted subtrees and constraints\n";
			std::cerr << "Removed " << deduplicate_constraints(constraints)
			          << " duplicate constraints" << std::endl;
			auto root_split = terraces::root_split(rerooted_tree, num_leaves);
			if (m == mode::count) {
				tree_enumerator<variants::count_callback<big_integer>> e{{}};
				auto c1 = e.run(num_leaves, constraints, root_split);
				auto c2 = e.run(num_leaves, constraints);
				std::cerr << "Counted supertrees" << std::endl;
				std::cout << c1 << "\t" << c2 << std::endl;
			} else {
				tree_enumerator<variants::multitree_callback> e{{}};
				auto t1 = e.run(num_leaves, constraints, root_split);
				auto t2 = e.run(num_leaves, constraints);
				std::cerr << "Checking supertrees" << std::endl;
				std::cout << is_comprehensive(rerooted_tree, data.matrix) << "\t";
				count_and_check_trees(t1, data.matrix, subtrees, data.names);
				std::cout << "\t";
				count_and_check_trees(t2, data.matrix, subtrees, data.names);
				std::cout << std::endl;
			}
		}
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
	}
}