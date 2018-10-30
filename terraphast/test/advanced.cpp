#include <catch.hpp>

#include <terraces/advanced.hpp>
#include <terraces/errors.hpp>
#include <terraces/parser.hpp>

namespace terraces {
namespace tests {

TEST_CASE("maximum_comprehensive_columnset", "[advanced-api]") {
	auto matrix_full = bitmatrix{4, 4};
	matrix_full.set(0, 2, 1);
	matrix_full.set(1, 2, 1);
	matrix_full.set(1, 3, 1);
	matrix_full.set(2, 2, 1);
	matrix_full.set(3, 0, 1);
	auto matrix_reduced = bitmatrix{4, 2};
	matrix_reduced.set(0, 0, 1);
	matrix_reduced.set(1, 0, 1);
	matrix_reduced.set(1, 1, 1);
	matrix_reduced.set(2, 0, 1);
	auto matrix_result = maximum_comprehensive_columnset(matrix_full);
	CHECK(find_comprehensive_taxon(matrix_result) == 1);
	CHECK(matrix_result == matrix_reduced);
}

TEST_CASE("advanced_exceptions", "[advanced-api]") {
	auto magic = 65u;
	name_map nums{"root"};
	index_map indx{{"root", 0}};
	/**
	 * Construct 65 nested three-taxon trees
	 * that create at least 65 sets in the first recursion step
	 * and thus at least 2^64 bipartitions.
	 */
	std::stringstream nwk;
	nwk << "(root,(";
	for (unsigned i = 0; i < 3 * magic; i += 3) {
		nums.push_back(std::to_string(i));
		indx.emplace(nums.back(), i + 1);
		nums.push_back(std::to_string(i + 1));
		indx.emplace(nums.back(), i + 2);
		nums.push_back(std::to_string(i + 2));
		indx.emplace(nums.back(), i + 3);
		nwk << "((" << i << ',' << (i + 1) << ")," << (i + 2) << ')';
		nwk << (i < 3 * (magic - 2) ? ",(" : (i < 3 * (magic - 1) ? "," : ""));
	}
	for (unsigned i = 0; i < magic; ++i) {
		nwk << ')';
	}
	auto tree = parse_nwk(nwk.str(), indx);
	bitmatrix matrix{3 * magic + 1, magic};
	for (unsigned i = 0; i < 3 * magic; ++i) {
		matrix.set(i + 1, i / 3, 1);
	}
	for (unsigned i = 0; i < magic; ++i) {
		matrix.set(0, i, 1);
	}
	auto d = create_supertree_data(tree, matrix);
	// all should return a large value or throw an overflow error
	CHECK(check_terrace(d));
	CHECK(fast_count_terrace(d) > 1);
	CHECK(count_terrace(d) == std::numeric_limits<decltype(count_terrace(d))>::max());
	CHECK_THROWS_AS(count_terrace_bigint(d), tree_count_overflow_error);

	// check mismatching sizes
	bitmatrix m_1x1{1, 1};
	CHECK_THROWS_AS(create_supertree_data(tree, m_1x1), bad_input_error);

	// check tiny input data
	auto trivial = parse_new_nwk("(s1,s2)").tree;
	bitmatrix mat_trivial{2, 1};
	bitmatrix mat_trivial_one{2, 1};
	mat_trivial_one.set(0, 0, 1);
	CHECK_THROWS_AS(create_supertree_data(trivial, mat_trivial_one), bad_input_error);
	CHECK_THROWS_AS(create_supertree_data(trivial, mat_trivial), no_usable_root_error);

	auto t23 = parse_new_nwk("(((((((((((((((((((((((((x,(y,z)),1),2),3),4),5),6),7),8),9),a),"
	                         "b),c),d),e),f),g),h),i),j),k),l),m),n),o)")
	                   .tree;
	bitmatrix t23_empty{27, 1};
	t23_empty.set(0, 0, 1);
	auto d23 = create_supertree_data(t23, t23_empty);
#ifndef USE_GMP
	CHECK_THROWS_AS(count_terrace_bigint(d23), tree_count_overflow_error);
#endif
	CHECK(count_terrace(d23) == std::numeric_limits<decltype(count_terrace(d23))>::max());
}

occurrence_data parse_bitmatrix_str(const std::string& str) {
	std::stringstream ss(str);
	return parse_bitmatrix(ss);
}

TEST_CASE("advanced_results", "[advanced-api]") {
	auto m1 = parse_bitmatrix_str("4 2\n1 0 s1\n1 1 s2\n1 1 s3\n1 1 s4");
	auto t1 = parse_nwk("(s1, (s2, (s3, s4)))", m1.indices);
	auto m2 = parse_bitmatrix_str("5 2\n1 0 s1\n1 0 s2\n1 1 s3\n0 1 s4\n0 1 s5");
	auto t2 = parse_nwk("(s3, ((s1, s2), (s4, s5)))", m2.indices);
	auto m3 = parse_bitmatrix_str(
	        "6 3\n1 0 0 s1\n1 0 0 s2\n0 0 1 s3\n0 1 1 s4\n1 1 1 s5\n0 1 1 s6");
	auto t3 = parse_nwk("((s4, (s3, (s2, (s1, s6)))), s5)", m3.indices);
	auto d1 = create_supertree_data(t1, m1.matrix);
	auto d2 = create_supertree_data(t2, m2.matrix);
	auto d3 = create_supertree_data(t3, m3.matrix);
	CHECK(!check_terrace(d1));
	CHECK(fast_count_terrace(d1) == 1);
	CHECK(count_terrace(d1) == 1);
	CHECK(count_terrace_bigint(d1).value() == 1);
	CHECK(check_terrace(d2));
	CHECK(fast_count_terrace(d2) > 1);
	CHECK(count_terrace(d2) > 1);
	CHECK(count_terrace_bigint(d2).value() > 1);
	CHECK(check_terrace(d3));
	CHECK(fast_count_terrace(d3) > 1);
	CHECK(count_terrace(d3) == 35);
	CHECK(count_terrace_bigint(d3).value() == 35);
	std::stringstream ss;
	print_terrace_compressed(d3, m3.names, ss);
	CHECK(ss.str() == "(s5,(s2,((s3,s6),(s1,s4))|(s4,{s1,s3,s6})|((s4,(s3,s6)),s1))|((s3,"
	                  "s6),{s1,s2,s4})|({s2,s3,s6},(s1,s4))|(s4,{s1,s2,s3,s6})|((s2,s4),{"
	                  "s1,s3,s6})|((s4,(s3,s6)),(s1,s2))|(((s3,s6),(s2,s4))|(s4,{s2,s3,s6})"
	                  "|((s4,(s3,s6)),s2),s1))");
	ss.str("");
	print_terrace(d3, m3.names, ss);
	REQUIRE(ss.str() == "(s5,(s2,((s3,s6),(s1,s4))));\n"
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
	std::stringstream ss2;
	enumerate_terrace(d3, [&](const tree& t) { ss2 << as_newick(t, m3.names) << '\n'; });
	CHECK(ss.str() == ss2.str());
}

} // namespace tests
} // namespace terraces
