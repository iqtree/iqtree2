#include <catch.hpp>

#include "../lib/bits.hpp"

namespace terraces {
namespace tests {

using bits::popcount;
using bits::bitscan;
using bits::rbitscan;
using bits::partial_popcount;
using bits::next_bit;
using bits::has_next_bit;

TEST_CASE("popcount", "[bits]") {
	CHECK(popcount(0b00000000000000000000000000000000ll) == 0);
	CHECK(popcount(0b00001001000010001010000010001100ll) == 8);
}

TEST_CASE("bitscan", "[bits]") {
	CHECK(bitscan(0b00000000000000000000000000000100ll) == 2);
	CHECK(bitscan(0b00000000000000000000100000000000ll) == 11);
	CHECK(bitscan(0b00000100000000000000000000000000ll) == 26);

	CHECK(rbitscan(0b00000000000000000000000000000100ll) == 2);
	CHECK(rbitscan(0b00000000000000000000100000000000ll) == 11);
	CHECK(rbitscan(0b10000000000000000000000000000000ll) == 31);
}

TEST_CASE("prefix_mask", "[bits]") {
	CHECK(bits::prefix_mask(23) == 0b00000000011111111111111111111111ll);
	CHECK(bits::prefix_mask(1) == 0b00000000000000000000000000000001ll);
	CHECK(bits::prefix_mask(10) == 0b00000000000000000000001111111111ll);
}

TEST_CASE("partial_popcount", "[bits]") {
	//                             33222222222211111111110000000000
	//                             10987654321098765432109876543210
	CHECK(bits::partial_popcount(0b01000011000010111100000101001000ll, 3) == 0);
	CHECK(bits::partial_popcount(0b01000011000010111100000101001000ll, 4) == 1);
	CHECK(bits::partial_popcount(0b01000011000010111100000101001000ll, 16) == 5);
	CHECK(bits::partial_popcount(0b01000011000010111100000101001000ll, 30) == 10);
}

TEST_CASE("next_bit", "[bits]") {
	CHECK(bits::has_next_bit(0b00010000000000000000000000000000ll, 28));
	CHECK(!bits::has_next_bit(0b00010000000000000000000000000000ll, 29));
	//                     33222222222211111111110000000000
	//                     10987654321098765432109876543210
	CHECK(bits::next_bit(0b10101110101011010100010101000001ll, 31) == 31);
	CHECK(bits::next_bit(0b10101110101011010100010101000001ll, 0) == 0);
	CHECK(bits::next_bit(0b10101110101011010100010101000001ll, 15) == 16);
	CHECK(bits::next_bit(0b10101110101011010100010101000001ll, 20) == 21);
	CHECK(bits::next_bit(0b10101110101011010100010101000001ll, 26) == 26);
}

} // namespace tests
} // namespace terraces
