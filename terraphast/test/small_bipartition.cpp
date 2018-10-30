#include <catch.hpp>

#include "../lib/small_bipartition.hpp"
#include "../lib/validation.hpp"

namespace terraces {
namespace tests {

TEST_CASE("small_bipartition", "[unconstrained]") {
	small_bipartition bip{0b00101001011100};
	auto step = [&](index i) {
		CHECK(bip.is_valid());
		CHECK(bip.left_mask() == i);
		bip.next();
	};
	CHECK(bip.leftmost_leaf() == 2);
	CHECK(bip.rightmost_leaf() == 11);
	CHECK(bip.num_leaves() == 6);

	step(0b00000000000100);
	step(0b00000000001000);
	step(0b00000000001100);
	step(0b00000000010000);
	step(0b00000000010100);
	step(0b00000000011000);
	step(0b00000000011100);
	step(0b00000001000000);
	step(0b00000001000100);
	step(0b00000001001000);
	step(0b00000001001100);
	step(0b00000001010000);
	step(0b00000001010100);
	step(0b00000001011000);
	step(0b00000001011100);
	step(0b00001000000000);
	step(0b00001000000100);
	step(0b00001000001000);
	step(0b00001000001100);
	step(0b00001000010000);
	step(0b00001000010100);
	step(0b00001000011000);
	step(0b00001000011100);
	step(0b00001001000000);
	step(0b00001001000100);
	step(0b00001001001000);
	step(0b00001001001100);
	step(0b00001001010000);
	step(0b00001001010100);
	step(0b00001001011000);
	step(0b00001001011100);
	CHECK(!bip.is_valid());
}

} // namespace tests
} // namespace terraces
