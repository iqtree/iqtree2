#include <catch.hpp>

#include <iostream>
#include <string>

#include "../lib/utils.hpp"

namespace terraces {
namespace tests {

TEST_CASE("skip_ws", "[utils][utils::skip_ws]") {
	const auto str1 = std::string{""};
	CHECK(utils::skip_ws(str1.begin(), str1.end()) == str1.end());

	const auto str2 = std::string{"  \n\t  "};
	CHECK(utils::skip_ws(str2.begin(), str2.end()) == str2.end());

	const auto str3 = std::string{"foo"};
	CHECK(utils::skip_ws(str3.begin(), str3.end()) == str3.begin());

	const auto str4 = std::string{"   \t bar"};
	CHECK(utils::skip_ws(str4.begin(), str4.end()) == str4.end() - 3);
}

TEST_CASE("reverse_skip_ws", "[utils][utils::reverse_skip_ws]") {
	const auto str1 = std::string{""};
	CHECK(utils::reverse_skip_ws(str1.begin(), str1.end()) == str1.begin());

	const auto str2 = std::string{"  \n\t  "};
	CHECK(utils::reverse_skip_ws(str2.begin(), str2.end()) == str2.begin());

	const auto str3 = std::string{"foo"};
	CHECK(utils::reverse_skip_ws(str3.begin(), str3.end()) == str3.end());

	const auto str4 = std::string{"bar   \t "};
	CHECK(utils::reverse_skip_ws(str4.begin(), str4.end()) == str4.begin() + 3);
}

TEST_CASE("ensure(false)", "[utils][utils::ensure]") {
	CHECK_THROWS_AS(utils::ensure<std::exception>(false), std::exception);
}

TEST_CASE("ensure(true)", "[utils][utils::ensure]") {
	CHECK_NOTHROW(utils::ensure<std::exception>(true));
}

} // namespace tests
} // namespace terraces
