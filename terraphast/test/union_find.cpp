#include <catch.hpp>

#include <iostream>

#include "../lib/union_find.hpp"

namespace terraces {
namespace tests {

TEST_CASE("union_find1", "[union_find]") {
	auto fl = utils::free_list{};
	auto alloc = utils::stack_allocator<index>{fl, 3};
	union_find leaves(3, alloc);
	leaves.merge(0, 1);
	CHECK(leaves.find(0) == leaves.find(1));
	CHECK(leaves.find(2) == 2);
}

TEST_CASE("union_find2", "[union_find]") {
	auto fl = utils::free_list{};
	auto alloc = utils::stack_allocator<index>{fl, 5};
	union_find leaves(5, alloc);
	leaves.merge(0, 1);
	leaves.merge(1, 4);
	CHECK(leaves.find(0) == leaves.find(1));
	CHECK(leaves.find(1) == leaves.find(4));
	CHECK(leaves.find(2) == 2);
	CHECK(leaves.find(3) == 3);
}

TEST_CASE("union_find::make_bipartition", "[union_find]") {
	auto fl = utils::free_list{};
	auto alloc = utils::stack_allocator<index>{fl, 8};
	std::vector<bool> b1{1, 0, 1, 1, 0, 1, 0, 0};
	std::vector<bool> b2{0, 1, 1, 1, 0, 0, 0, 0};
	std::vector<bool> b3{1, 1, 1, 1, 1, 1, 1, 1};
	auto check = [&](const std::vector<bool>& vec) {
		auto uf = union_find::make_bipartition(vec, alloc);
		auto i0 = std::distance(vec.begin(), std::find(vec.begin(), vec.end(), 0));
		auto i1 = std::distance(vec.begin(), std::find(vec.begin(), vec.end(), 1));
		for (index i = 0; i < uf.size(); ++i) {
			CHECK(uf.find(i) == uf.find(vec[i] ? i1 : i0));
		}
	};
	check(b1);
	check(b2);
	check(b3);
}

} // namespace tests
} // namespace terraces
