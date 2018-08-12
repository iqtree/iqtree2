#include <catch.hpp>

#include <algorithm>

#include "../lib/trees_impl.hpp"

namespace terraces {
namespace tests {

TEST_CASE("is_root(root)", "[trees]") {
	auto root_node = terraces::node{none, 1, 2, none};
	CHECK(terraces::is_root(root_node));
}

TEST_CASE("is_root(non_root)", "[trees]") {
	auto non_root_node = terraces::node{1, 2, 3, none};
	CHECK(!terraces::is_root(non_root_node));
}

TEST_CASE("is_leaf(leaf)", "[trees]") {
	auto leaf_node = terraces::node{0, none, none, none};
	CHECK(terraces::is_leaf(leaf_node));
}

TEST_CASE("is_leaf(non_leaf)", "[trees]") {
	auto non_leaf_node = terraces::node{0, 1, 2, none};
	CHECK(!terraces::is_leaf(non_leaf_node));
}

TEST_CASE("is_rooted_tree(valid)", "[trees]") {
	tree t{
	        {none, 1, 2, none},    {0, 9, 6, none},       {0, 3, 4, none},
	        {2, 5, 8, none},       {2, none, none, none}, {3, none, none, none},
	        {1, none, none, none}, {9, none, none, none}, {3, none, none, none},
	        {1, 10, 7, none},      {9, none, none, none},
	};
	check_rooted_tree(t);
}

TEST_CASE("foreach_postorder(example)", "[trees]") {
	tree t{
	        {none, 1, 2, none},    {0, 9, 6, none},       {0, 3, 4, none},
	        {2, 5, 8, none},       {2, none, none, none}, {3, none, none, none},
	        {1, none, none, none}, {9, none, none, none}, {3, none, none, none},
	        {1, 10, 7, none},      {9, none, none, none},
	};
	std::vector<index> expected{10, 7, 9, 6, 1, 5, 8, 3, 4, 2, 0};
	std::vector<index> result;
	foreach_postorder(t, [&](index i) { result.push_back(i); });
	CHECK(result == expected);
}

TEST_CASE("foreach_preorder(example)", "[trees]") {
	tree t{
	        {none, 1, 2, none},    {0, 9, 6, none},       {0, 3, 4, none},
	        {2, 5, 8, none},       {2, none, none, none}, {3, none, none, none},
	        {1, none, none, none}, {9, none, none, none}, {3, none, none, none},
	        {1, 10, 7, none},      {9, none, none, none},
	};
	std::vector<index> expected{0, 1, 9, 10, 7, 6, 2, 3, 5, 8, 4};
	std::vector<index> result;
	foreach_preorder(t, [&](index i) { result.push_back(i); });
	CHECK(result == expected);
}

TEST_CASE("foreach_postorder(trivial)", "[trees]") {
	tree t{{}};
	std::vector<index> result;
	foreach_postorder(t, [&](index i) { result.push_back(i); });
	CHECK(result.size() == 1);
	CHECK(result[0] == 0);
}

TEST_CASE("foreach_preorder(trivial)", "[trees]") {
	tree t{{}};
	std::vector<index> result;
	foreach_preorder(t, [&](index i) { result.push_back(i); });
	CHECK(result.size() == 1);
	CHECK(result[0] == 0);
}

TEST_CASE("tree_printing", "[trees][tree-printing]") {
	auto t = tree{{none, 1, 2, none},
	              {0, none, none, 1},
	              {0, 3, 4, none},
	              {2, none, none, 3},
	              {2, none, none, 4}};
	const auto names = name_map{"root", "foo", "", "bar", "baz"};
	std::ostringstream stream;
	stream << as_newick(t, names);
	CHECK(stream.str() == "(foo,(bar,baz));");
}

} // namespace tests
} // namespace terraces
