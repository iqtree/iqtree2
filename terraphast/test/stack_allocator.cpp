#include <catch.hpp>

#include "../lib/stack_allocator.hpp"

namespace terraces {
namespace tests {

TEST_CASE("stack_allocator basic", "[utils][utils::stack_allocator]") {
	utils::free_list fl;
	utils::stack_allocator<int> alloc{fl, 10};
	std::vector<int, utils::stack_allocator<int>> vec(alloc);
	CHECK_NOTHROW(vec.reserve(10));
}

} // namespace tests
} // namespace terraces
