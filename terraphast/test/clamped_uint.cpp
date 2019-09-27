#include <catch.hpp>

#include <terraces/clamped_uint.hpp>

#include <limits>
#include <terraces/errors.hpp>

namespace terraces {
namespace tests {

TEST_CASE("clamped_uint", "[checked_uint]") {
	auto max = std::numeric_limits<index>::max();
	CHECK((clamped_uint{10} + clamped_uint{417}).value() == 10 + 417);
	CHECK((clamped_uint{10} * clamped_uint{417}).value() == 10 * 417);
	CHECK((clamped_uint{max} + clamped_uint{1}).is_clamped());
	CHECK((clamped_uint{max} + clamped_uint{1}).value() == max);
	CHECK((clamped_uint{max / 2} * clamped_uint{3}).is_clamped());
	CHECK((clamped_uint{max / 2} * clamped_uint{3}).value() == max);
}

TEST_CASE("overflow_except_uint", "[checked_uint]") {
	auto max = std::numeric_limits<index>::max();
	CHECK((overflow_except_uint{10} + overflow_except_uint{417}).value() == 10 + 417);
	CHECK((overflow_except_uint{10} * overflow_except_uint{417}).value() == 10 * 417);
	CHECK_THROWS_AS(overflow_except_uint{max} + overflow_except_uint{1},
	                terraces::tree_count_overflow_error);
	CHECK_THROWS_AS(overflow_except_uint{max / 2} * overflow_except_uint{3},
	                terraces::tree_count_overflow_error);
}

} // namespace tests
} // namespace terraces
