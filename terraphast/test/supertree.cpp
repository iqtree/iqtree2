#include <catch.hpp>

#include <iostream>

#include "../lib/supertree_enumerator.hpp"
#include "../lib/supertree_variants.hpp"

namespace terraces {
namespace tests {

uint64_t count_supertree(index num_leaves, const constraints& constraints) {
	tree_enumerator<variants::count_callback<uint64_t>> e{{}};
	return e.run(num_leaves, constraints);
}

bool check_supertree(index num_leaves, const constraints& constraints) {
	tree_enumerator<variants::check_callback> e{{}};
	return e.run(num_leaves, constraints) > 1;
}

TEST_CASE("count_supertree1", "[supertree]") {
	constraints c = {};
	CHECK(count_supertree(2, c) == 1);
	CHECK(!check_supertree(2, c));
}

TEST_CASE("count_supertree2", "[supertree]") {
	constraints c = {};
	CHECK(count_supertree(3, c) == 3);
	CHECK(check_supertree(3, c));
}

TEST_CASE("count_supertree3", "[supertree]") {
	constraints c = {};
	CHECK(count_supertree(7, c) == 10395);
	CHECK(check_supertree(7, c));
}

TEST_CASE("count_supertree4", "[supertree]") {
	constraints c = {{0, 1, 2}};
	CHECK(count_supertree(3, c) == 1);
	CHECK(!check_supertree(3, c));
}

TEST_CASE("count_supertree5", "[supertree]") {
	constraints c = {{0, 1, 2}, {2, 3, 4}};
	CHECK(count_supertree(5, c) == 9);
	CHECK(check_supertree(5, c));
}

TEST_CASE("count_supertree6", "[supertree]") {
	constraints c = {{1, 0, 2}, {3, 4, 1}};
	CHECK(count_supertree(5, c) == 9);
	CHECK(check_supertree(5, c));
}

TEST_CASE("count_supertree7", "[supertree]") {
	constraints c = {{0, 1, 3}, {3, 2, 0}, {4, 5, 6}, {6, 3, 4}, {2, 3, 6}, {2, 6, 7}};
	CHECK(count_supertree(8, c) == 173);
	CHECK(check_supertree(8, c));
}

TEST_CASE("count_supertree_none", "[supertree]") {
	constraints c = {{0, 1, 2}, {2, 1, 0}};
	CHECK(count_supertree(3, c) == 0);
	CHECK(!check_supertree(3, c));
}

} // namespace tests
} // namespace terraces
