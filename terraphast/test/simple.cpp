#include <catch.hpp>
#include <sstream>
#include <terraces/errors.hpp>
#include <terraces/simple.hpp>

using namespace terraces::simple;

namespace terraces {
namespace tests {

TEST_CASE("simple_errors") {
	// malformed nwk
	CHECK_THROWS_AS(is_on_terrace("(s1,", "1 1\n1 s1"), bad_input_error);
	// mismatching sizes
	CHECK_THROWS_AS(is_on_terrace("(s1,s2)", "1 1\n1 s1"), bad_input_error);
	// malformed matrix
	CHECK_THROWS_AS(is_on_terrace("(s1,s2)", "1 1\na s1"), bad_input_error);
	// multifurcating newick tree
	CHECK_THROWS_AS(is_on_terrace("(s1,s2,s3,s4)", "4 1\n1 s1\n1 s2\n1 s3\n1 s4"),
	                bad_input_error);
	// unknown/empty name
	CHECK_THROWS_AS(is_on_terrace("((s1,s2),(,s5))", "4 2\n1 0 s1\n1 0 s2\n0 1 s3\n0 1 s4"),
	                bad_input_error);
	// no root
	CHECK_THROWS_AS(is_on_terrace("((s1,s2),(s3,s4))", "4 2\n1 0 s1\n1 0 s2\n0 1 s3\n0 1 s4"),
	                no_usable_root_error);
}

TEST_CASE("simple_results") {
	CHECK(!is_on_terrace("(s1, (s2, (s3, s4)))", "4 2\n1 0 s1\n1 1 s2\n1 1 s3\n1 1 s4"));
	CHECK(is_on_terrace("(s3, ((s1, s2), (s4, s5)))",
	                    "5 2\n1 0 s1\n1 0 s2\n1 1 s3\n0 1 s4\n0 1 s5"));
	CHECK(get_terrace_size("((s4, (s3, (s2, (s1, s6)))), s5)",
	                       "6 3\n1 0 0 s1\n1 0 0 s2\n0 0 1 s3\n0 1 1 s4\n1 1 1 s5\n0 1 1 s6") ==
	      35);
	std::stringstream ss;
	print_terrace_compressed("((s4, (s3, (s2, (s1, s6)))), s5)",
	                         "6 3\n1 0 0 s1\n1 0 0 s2\n0 0 1 s3\n0 1 1 s4\n1 1 1 s5\n0 1 1 s6",
	                         ss);
	CHECK(ss.str() == "(s5,(s2,((s3,s6),(s1,s4))|(s4,{s1,s3,s6})|((s4,(s3,s6)),s1))|((s3,s6),{"
	                  "s1,s2,s4})|({s2,s3,s6},(s1,s4))|(s4,{s1,s2,s3,s6})|((s2,s4),{s1,s3,s6})|"
	                  "((s4,(s3,s6)),(s1,s2))|(((s3,s6),(s2,s4))|(s4,{s2,s3,s6})|((s4,(s3,s6)),"
	                  "s2),s1))");
	ss.str("");
	print_terrace("((s4, (s3, (s2, (s1, s6)))), s5)",
	              "6 3\n1 0 0 s1\n1 0 0 s2\n0 0 1 s3\n0 1 1 s4\n1 1 1 s5\n0 1 1 s6", ss);
	CHECK(ss.str() == "(s5,(s2,((s3,s6),(s1,s4))));\n"
	                  "(s5,(s2,(s4,(s1,(s3,s6)))));\n"
	                  "(s5,(s2,(s4,(s3,(s1,s6)))));\n"
	                  "(s5,(s2,(s4,((s1,s3),s6))));\n"
	                  "(s5,(s2,((s4,(s3,s6)),s1)));\n"
	                  "(s5,((s3,s6),(s1,(s2,s4))));\n"
	                  "(s5,((s3,s6),(s2,(s1,s4))));\n"
	                  "(s5,((s3,s6),((s1,s2),s4)));\n"
	                  "(s5,((s2,(s3,s6)),(s1,s4)));\n"
	                  "(s5,((s3,(s2,s6)),(s1,s4)));\n"
	                  "(s5,(((s2,s3),s6),(s1,s4)));\n"
	                  "(s5,(s4,(s1,(s2,(s3,s6)))));\n"
	                  "(s5,(s4,(s1,(s3,(s2,s6)))));\n"
	                  "(s5,(s4,(s1,((s2,s3),s6))));\n"
	                  "(s5,(s4,(s2,(s1,(s3,s6)))));\n"
	                  "(s5,(s4,(s2,(s3,(s1,s6)))));\n"
	                  "(s5,(s4,(s2,((s1,s3),s6))));\n"
	                  "(s5,(s4,((s1,s2),(s3,s6))));\n"
	                  "(s5,(s4,(s3,(s1,(s2,s6)))));\n"
	                  "(s5,(s4,(s3,(s2,(s1,s6)))));\n"
	                  "(s5,(s4,(s3,((s1,s2),s6))));\n"
	                  "(s5,(s4,((s1,s3),(s2,s6))));\n"
	                  "(s5,(s4,((s2,s3),(s1,s6))));\n"
	                  "(s5,(s4,((s1,(s2,s3)),s6)));\n"
	                  "(s5,(s4,((s2,(s1,s3)),s6)));\n"
	                  "(s5,(s4,(((s1,s2),s3),s6)));\n"
	                  "(s5,((s2,s4),(s1,(s3,s6))));\n"
	                  "(s5,((s2,s4),(s3,(s1,s6))));\n"
	                  "(s5,((s2,s4),((s1,s3),s6)));\n"
	                  "(s5,((s4,(s3,s6)),(s1,s2)));\n"
	                  "(s5,(((s3,s6),(s2,s4)),s1));\n"
	                  "(s5,((s4,(s2,(s3,s6))),s1));\n"
	                  "(s5,((s4,(s3,(s2,s6))),s1));\n"
	                  "(s5,((s4,((s2,s3),s6)),s1));\n"
	                  "(s5,(((s4,(s3,s6)),s2),s1));\n");
}

TEST_CASE("simple_results_force") {
	CHECK(get_terrace_size("((s4, (s3, (s2, (s1, s6)))), s5)", "6 5\n0 1 0 0 0 s1\n0 1 0 0 0 "
	                                                           "s2\n0 0 0 0 1 s3\n0 0 1 0 1 "
	                                                           "s4\n0 1 1 0 1 s5\n0 0 1 0 1 s6",
	                       true) == 35);
}

} // namespace tests
} // namespace terraces
