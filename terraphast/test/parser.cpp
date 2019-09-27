#include <catch.hpp>

#include <iostream>
#include <sstream>

#include <terraces/advanced.hpp>
#include <terraces/errors.hpp>
#include <terraces/parser.hpp>

namespace terraces {
namespace tests {

TEST_CASE("parsing a tree with just a root-node", "[parser]") {
	const auto results = parse_new_nwk("foo");
	const auto& tree = results.tree;
	const auto& names = results.names;
	const auto& indices = results.indices;
	REQUIRE(tree.size() == 1);
	CHECK(tree[0].parent() == none);
	CHECK(tree[0].lchild() == none);
	CHECK(tree[0].rchild() == none);
	CHECK(names[0] == "foo");
	CHECK(indices.at("foo") == 0);
}

TEST_CASE("parsing a tree with three leaves and two inner nodes", "[parser]") {
	const auto results = parse_new_nwk("((foo,bar)inner, baz)outer");
	const auto& tree = results.tree;
	const auto& names = results.names;
	const auto& indices = results.indices;
	REQUIRE(tree.size() == 5);
	CHECK(tree[0].parent() == none);
	CHECK(tree[0].lchild() == 1);
	CHECK(tree[0].rchild() == 4);
	CHECK(tree[1].parent() == 0);
	CHECK(tree[1].lchild() == 2);
	CHECK(tree[1].rchild() == 3);
	CHECK(tree[2].parent() == 1);
	CHECK(tree[2].lchild() == none);
	CHECK(tree[2].rchild() == none);
	CHECK(tree[3].parent() == 1);
	CHECK(tree[3].lchild() == none);
	CHECK(tree[3].rchild() == none);
	CHECK(tree[4].parent() == 0);
	CHECK(tree[4].lchild() == none);
	CHECK(tree[4].rchild() == none);
	CHECK(names[0] == "foo");
	CHECK(names[1] == "bar");
	CHECK(names[2] == "baz");
	CHECK(indices.at("foo") == 0);
	CHECK(indices.at("bar") == 1);
	CHECK(indices.at("baz") == 2);
	CHECK(indices.size() == 3);
}

TEST_CASE("parsing a tree with quotes, three leaves and two inner nodes", "[parser]") {
	const auto results = parse_new_nwk("((foo,bar)'inner', 'baz')outer");
	const auto& tree = results.tree;
	const auto& names = results.names;
	const auto& indices = results.indices;
	REQUIRE(tree.size() == 5);
	CHECK(tree[0].parent() == none);
	CHECK(tree[0].lchild() == 1);
	CHECK(tree[0].rchild() == 4);
	CHECK(tree[1].parent() == 0);
	CHECK(tree[1].lchild() == 2);
	CHECK(tree[1].rchild() == 3);
	CHECK(tree[2].parent() == 1);
	CHECK(tree[2].lchild() == none);
	CHECK(tree[2].rchild() == none);
	CHECK(tree[3].parent() == 1);
	CHECK(tree[3].lchild() == none);
	CHECK(tree[3].rchild() == none);
	CHECK(tree[4].parent() == 0);
	CHECK(tree[4].lchild() == none);
	CHECK(tree[4].rchild() == none);
	CHECK(names[0] == "foo");
	CHECK(names[1] == "bar");
	CHECK(names[2] == "baz");
	CHECK(indices.at("foo") == 0);
	CHECK(indices.at("bar") == 1);
	CHECK(indices.at("baz") == 2);
	CHECK(indices.size() == 3);
}

TEST_CASE("parsing an unrooted tree with three leaves", "[parser]") {
	const auto results = parse_new_nwk("(foo,bar,'baz')");
	auto exp_names = name_map{"foo", "bar", "baz"};
	auto exp_tree = tree{{none, 1, 3, none},
	                     {0, none, none, 0},
	                     {3, none, none, 1},
	                     {0, 2, 4, none},
	                     {3, none, none, 2}};
	CHECK(results.names == exp_names);
	CHECK(results.tree == exp_tree);
}

TEST_CASE("parsing tree with unnecessary parentheses", "[parser]") {
	CHECK_THROWS_AS(parse_new_nwk("()"), bad_input_error);
}

TEST_CASE("parsing trees with mismatching parentheses", "[parser]") {
	// too many closing parentheses
	CHECK_THROWS_AS(parse_new_nwk("((a,b),c))"), bad_input_error);
	// too many opening parentheses
	CHECK_THROWS_AS(parse_new_nwk("((a,b)"), bad_input_error);
	// too many opening parentheses
	CHECK_THROWS_AS(parse_new_nwk("((a,b),c"), bad_input_error);
	// ternary nodes (simple)
	CHECK_THROWS_AS(parse_new_nwk("((a,b,c),a)"), bad_input_error);
	// ternary nodes (complex)
	CHECK_THROWS_AS(parse_new_nwk("(a,((b,c),((d,e),f),g))"), bad_input_error);
}

TEST_CASE("parsing trees invalid format", "[parser]") {
	// inner node names must come after their children.
	CHECK_THROWS_AS(parse_new_nwk("a(,)"), bad_input_error);
}

TEST_CASE("parsing trees with unclosed quotes", "[parser]") {
	CHECK_THROWS_AS(parse_new_nwk("(('a',''),(('c,),),)"), bad_input_error);
}

TEST_CASE("parsing trees with duplicate taxa", "[parser]") {
	CHECK_THROWS_AS(parse_new_nwk("(a,a)"), bad_input_error);
	index_map inds{{"a", 0}};
	CHECK_THROWS_AS(parse_nwk("(a,a)", inds), bad_input_error);
}

TEST_CASE("parsing a datafile with three species and two cols", "[parser],[data-parser]") {
	auto stream = std::istringstream{"3 2\n0 1 foo\n1 1 bar\n1 1 baz\n"};
	const auto res = parse_bitmatrix(stream);
	const auto& mat = res.matrix;

	REQUIRE(mat.cols() == 2);
	REQUIRE(mat.rows() == 3);

	CHECK(!mat.get(0, 0));
	CHECK(mat.get(0, 1));
	CHECK(mat.get(1, 0));
	CHECK(mat.get(1, 1));
	CHECK(mat.get(2, 0));
	CHECK(mat.get(2, 1));

	CHECK(find_comprehensive_taxon(mat) == 1);
}

TEST_CASE("parsing a datafile with duplicate species", "[parser],[data-parser]") {
	auto stream = std::istringstream{"3 2\n1 0 foo\n1 1 bar\n1 1 bar\n"};
	CHECK_THROWS_AS(terraces::parse_bitmatrix(stream), bad_input_error);
}

TEST_CASE("parsing a complex datafile", "[parser],[data-parser]") {
	auto stream = std::istringstream{"5 3\n"
	                                 "0 1 0 foo\n"
	                                 "1 1 1 bar\n"
	                                 "0 1 0 bla\n"
	                                 "1 0 1 blub\n"
	                                 "1 1 0 gaehn\n"};
	const auto res = terraces::parse_bitmatrix(stream);
	const auto& mat = res.matrix;

	REQUIRE(mat.cols() == 3);
	REQUIRE(mat.rows() == 5);

	CHECK(!mat.get(0, 0));
	CHECK(mat.get(0, 1));
	CHECK(!mat.get(0, 2));
	CHECK(mat.get(1, 0));
	CHECK(mat.get(1, 1));
	CHECK(mat.get(1, 2));
	CHECK(!mat.get(2, 0));
	CHECK(mat.get(2, 1));
	CHECK(!mat.get(2, 2));
	CHECK(mat.get(3, 0));
	CHECK(!mat.get(3, 1));
	CHECK(mat.get(3, 2));
	CHECK(mat.get(4, 0));
	CHECK(mat.get(4, 1));
	CHECK(!mat.get(4, 2));

	CHECK(find_comprehensive_taxon(mat) == 1);
}

} // namespace tests
} // namespace terraces
