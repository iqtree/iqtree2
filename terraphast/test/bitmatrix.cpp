#include <catch.hpp>

#include <terraces/bitmatrix.hpp>

namespace terraces {
namespace tests {

TEST_CASE("bitmatrix-construction", "[bitmatrix]") {
	auto mat = bitmatrix{10, 5};
	CHECK(mat.rows() == 10);
	CHECK(mat.cols() == 5);
}

TEST_CASE("bitmatrix set/get", "[bitmatrix]") {
	auto mat = bitmatrix{3, 4};

	CHECK(!mat.get(1, 2));
	mat.set(1, 2, true);
	CHECK(mat.get(1, 2));
}

} // namespace tests
} // namespace terraces
