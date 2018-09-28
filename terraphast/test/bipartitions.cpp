#include <catch.hpp>

#include <iostream>

#include "../lib/bipartitions.hpp"
#include "../lib/ranked_bitvector.hpp"
#include "../lib/union_find.hpp"

namespace terraces {
namespace tests {

TEST_CASE("bipartition1", "[bipartition]") {
	auto fl = utils::free_list{};
	auto alloc = utils::stack_allocator<index>{fl, 4};
	union_find u(4, alloc);
	ranked_bitvector s{4, alloc};
	s.set(0);
	s.set(1);
	s.set(2);
	s.set(3);
	s.update_ranks();
	u.merge(0, 1);
	u.compress();
	bipartitions bip_it(s, u, alloc);
	CHECK(bip_it.end_bip() == 4);
	CHECK(bip_it.num_bip() == 3);
	CHECK(bip_it.begin_bip() == 1);
	auto set = bip_it.get_first_set(1, alloc);
	CHECK(!set.get(0));
	CHECK(!set.get(1));
	CHECK(set.get(2));
	CHECK(!set.get(3));
	set = bip_it.get_first_set(2, alloc);
	CHECK(!set.get(0));
	CHECK(!set.get(1));
	CHECK(!set.get(2));
	CHECK(set.get(3));
	set = bip_it.get_first_set(3, alloc);
	CHECK(!set.get(0));
	CHECK(!set.get(1));
	CHECK(set.get(2));
	CHECK(set.get(3));
	CHECK(bip_it.end_bip() == 4);
	CHECK(bip_it.num_bip() == 3);
}

} // namespace tests
} // namespace terraces
