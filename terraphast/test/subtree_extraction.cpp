#include <catch.hpp>

#include <terraces/subtree_extraction.hpp>

#include "../lib/trees_impl.hpp"
#include "../lib/validation.hpp"

namespace terraces {
namespace tests {

using std::vector;

TEST_CASE("subtree extraction: full data", "[subtree_extraction]") {
	tree t{{none, 4, 5, none}, {2, none, none, 0}, {4, 6, 1, none},   {4, none, none, 1},
	       {0, 2, 3, none},    {0, none, none, 2}, {2, none, none, 3}};

	bitmatrix bm{4, 1};
	for (index row = 0; row < bm.rows(); ++row) {
		bm.set(row, 0, true);
	}

	auto t2 = subtrees(t, bm)[0];
	check_rooted_tree(t);
	check_rooted_tree(t2);
	vector<index> exp_pre{0, 1, 2, 3, 4, 5, 6};
	vector<index> exp_post{3, 4, 2, 5, 1, 6, 0};
	auto res_pre = preorder(t2);
	auto res_post = postorder(t2);
	CHECK(t2[3].taxon() == 3);
	CHECK(t2[4].taxon() == 0);
	CHECK(t2[5].taxon() == 1);
	CHECK(t2[6].taxon() == 2);
	CHECK(is_isomorphic_unrooted(t, t2));
	CHECK(exp_pre == res_pre);
	CHECK(exp_post == res_post);
}

TEST_CASE("subtree extraction: example", "[subtree_extraction]") {
	tree t{{none, 4, 5, none}, {2, none, none, 0}, {4, 6, 1, none},   {4, none, none, 1},
	       {0, 2, 3, none},    {0, none, none, 2}, {2, none, none, 3}};

	bitmatrix bm{4, 2};
	bm.set(0, 0, true);
	bm.set(0, 1, true);
	bm.set(1, 0, true);
	bm.set(2, 0, true);
	bm.set(2, 1, true);
	bm.set(3, 1, true);

	auto trees = subtrees(t, bm);
	auto t1 = trees[0];
	auto t2 = trees[1];

	vector<index> exp_pre1{0, 1, 2, 3, 4};
	vector<index> exp_pre2{0, 1, 2, 3, 4};
	vector<index> exp_post1{2, 3, 1, 4, 0};
	vector<index> exp_post2{2, 3, 1, 4, 0};
	CHECK(t1[2].taxon() == 0);
	CHECK(t1[3].taxon() == 1);
	CHECK(t1[4].taxon() == 2);
	CHECK(t2[2].taxon() == 3);
	CHECK(t2[3].taxon() == 0);
	CHECK(t2[4].taxon() == 2);
	CHECK(exp_pre1 == preorder(trees[0]));
	CHECK(exp_pre2 == preorder(trees[1]));
	CHECK(exp_post1 == postorder(trees[0]));
	CHECK(exp_post2 == postorder(trees[1]));
}

TEST_CASE("subtree extraction: moved root", "[subtree_extraction]") {
	std::stringstream mss{"5 1\n1 1\n1 2\n1 3\n0 4\n0 5"};
	auto matrix = parse_bitmatrix(mss);
	auto to_str = [&](const tree& t) {
		std::stringstream ss;
		ss << as_newick(t, matrix.names);
		return ss.str();
	};
	auto t1 = parse_nwk("((4,5),(1,(2,3)));", matrix.indices);
	auto t2 = parse_nwk("((((1,2),3),4),5);", matrix.indices);
	auto t3 = parse_nwk("(4,((1,5),(2,3)));", matrix.indices);
	auto t4 = parse_nwk("(4,((2,3),(1,5)));", matrix.indices);
	auto t5 = parse_nwk("(5,((3,(1,4)),2));", matrix.indices);
	auto st1 = subtrees(t1, matrix.matrix)[0];
	auto st2 = subtrees(t2, matrix.matrix)[0];
	auto st3 = subtrees(t3, matrix.matrix)[0];
	auto st4 = subtrees(t4, matrix.matrix)[0];
	auto st5 = subtrees(t5, matrix.matrix)[0];
	CHECK(to_str(st1) == "(1,(2,3));");
	CHECK(to_str(st2) == "((1,2),3);");
	CHECK(to_str(st3) == "(1,(2,3));");
	CHECK(to_str(st4) == "((2,3),1);");
	CHECK(to_str(st5) == "((3,1),2);");
}

} // namespace tests
} // namespace terraces
