#include <catch.hpp>

#include <algorithm>
#include <sstream>

#include <terraces/advanced.hpp>
#include <terraces/parser.hpp>
#include <terraces/rooting.hpp>
#include <terraces/subtree_extraction.hpp>

#include "../lib/multitree_iterator.hpp"
#include "../lib/supertree_enumerator.hpp"
#include "../lib/supertree_variants_multitree.hpp"
#include "../lib/validation.hpp"

namespace terraces {
namespace tests {

void check_unique_trees(multitree_node* root, index num_trees) {
	multitree_iterator it(root);

	std::vector<std::vector<simple_bitvector>> bipartitions;
	do {
		bipartitions.push_back(tree_bipartitions(it.tree()));
	} while (it.next());
	// check that all trees are unique
	std::sort(bipartitions.begin(), bipartitions.end());
	CHECK(std::adjacent_find(bipartitions.begin(), bipartitions.end()) == bipartitions.end());
	// check that there is the right number of trees
	CHECK(bipartitions.size() == num_trees);
}

TEST_CASE("multitree_iterator init simple", "[multitree]") {
	auto data_stream = std::istringstream{"7 4\n1 1 1 1 s1\n0 0 1 0 s2\n1 1 0 0 s3\n1 1 1 0 "
	                                      "s4\n1 1 0 1 s5\n1 0 0 1 s7\n0 0 0 1 s13"};
	auto data = terraces::parse_bitmatrix(data_stream);
	auto tree = terraces::parse_nwk("((((s2,s4),((s13,s1),s7)),s3),s5);", data.indices);

	auto supertree_data = terraces::create_supertree_data(tree, data.matrix);
	tree_enumerator<variants::multitree_callback> enumerator{{}};
	auto result = enumerator.run(supertree_data.num_leaves, supertree_data.constraints,
	                             supertree_data.root);

	check_unique_trees(result, 9);
}

TEST_CASE("multitree_iterator init unconstrained", "[multitree]") {
	name_map names{"1", "2", "3", "4", "5", "6", "7", "8"};
	constraints constraints{{0, 1, 2}};
	index root_species = 0;
	tree_enumerator<variants::multitree_callback> enumerator{{}};
	auto result = enumerator.run(names.size(), constraints, root_species);

	check_unique_trees(result, count_unrooted_trees<index>(7));
}

} // namespace tests
} // namespace terraces
