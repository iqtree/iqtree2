#include <catch.hpp>

#include <terraces/advanced.hpp>
#include <terraces/constraints.hpp>
#include <terraces/parser.hpp>
#include <terraces/rooting.hpp>
#include <terraces/subtree_extraction.hpp>
namespace terraces {

namespace tests {

TEST_CASE("full_run_disjoint", "[rerooting],[tree_extraction],[constraints],[supertree]") {
	tree t{{none, 1, 2, none}, {0, 3, 4, none},    {0, none, none, 0},
	       {1, 5, 6, none},    {1, 7, 8, none},    {3, none, none, 1},
	       {3, none, none, 2}, {4, none, none, 3}, {4, none, none, 4}};
	bitmatrix b{5, 2};
	b.set(0, 0, true);
	b.set(0, 1, true);
	b.set(1, 1, true);
	b.set(2, 1, true);
	b.set(3, 0, true);
	b.set(4, 0, true);
	auto ts = subtrees(t, b);
	auto cs = compute_constraints(ts);
	auto names = name_map(t.size());
	index root = 0;
	CHECK(count_terrace({cs, 5, root}) == 15);
}

} // namespace tests
} // namespace terraces
