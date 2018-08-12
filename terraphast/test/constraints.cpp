#include <catch.hpp>

#include <terraces/constraints.hpp>
#include <terraces/subtree_extraction.hpp>

#include <algorithm>

#include "../lib/trees_impl.hpp"

namespace terraces {
namespace tests {

TEST_CASE("constraint extraction: full data", "[subtree_extraction],[constraints]") {
	tree t{{none, 4, 5, none}, {2, none, none, 0}, {4, 6, 1, none},   {4, none, none, 1},
	       {0, 2, 3, none},    {0, none, none, 2}, {2, none, none, 3}};

	bitmatrix bm{4, 1};
	for (index row = 0; row < bm.rows(); ++row) {
		bm.set(row, 0, true);
	}

	auto ts = subtrees(t, bm);
	auto result = compute_constraints(ts);
	auto required = constraints{{1, 3, 2}, {0, 3, 1}};
	CHECK(result == required);
}

TEST_CASE("constraint extraction: example", "[subtree_extraction],[constraints]") {
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

	auto result = compute_constraints(trees);
	auto required = constraints{{1, 0, 2}, {0, 3, 2}};
	CHECK(result == required);
}

TEST_CASE("constraint deduplication", "[deduplication], [contraints]") {
	auto dup = constraints{{0, 1, 2}, {0, 1, 2}, {3, 4, 5}, {4, 3, 5}, {7, 6, 8}, {6, 7, 8}};
	auto num = deduplicate_constraints(dup);
	CHECK(num == 3);
	CHECK(dup == (constraints{{0, 1, 2}, {3, 4, 5}, {6, 7, 8}}));
}

} // namespace tests
} // namespace terraces
