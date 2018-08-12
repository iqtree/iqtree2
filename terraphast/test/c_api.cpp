#include <catch.hpp>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <terraces/terraces.h>

TEST_CASE("c_api_error_nwk", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 3, 0, 1);
	setDataMatrix(data, 3, 1, 1);
	REQUIRE(terraceAnalysis(data, "(s1,", TA_DETECT, nullptr, result) == TERRACE_NEWICK_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_nwk_unknown", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 3, 0, 1);
	setDataMatrix(data, 3, 1, 1);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s5,s4)", TA_DETECT, nullptr, result) ==
	        TERRACE_NEWICK_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_species_count", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 3, 0, 1);
	setDataMatrix(data, 3, 1, 1);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3)", TA_DETECT, nullptr, result) ==
	        TERRACE_SPECIES_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_malformed_matrix", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 4);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3)", TA_DETECT, nullptr, result) ==
	        TERRACE_MATRIX_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_num_species", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3"};
	auto data = initializeMissingData(3, 2, names);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3)", TA_DETECT, nullptr, result) ==
	        TERRACE_NUM_SPECIES_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_overflow", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	auto magic = 65u;
	std::vector<std::string> nums{"root"};
	/**
	 * Construct 65 nested three-taxon trees
	 * that create at least 65 sets in the first recursion step
	 * and thus at least 2^64 bipartitions.
	 */
	std::stringstream nwk;
	nwk << "(root,(";
	for (unsigned i = 0; i < 3 * magic; i += 3) {
		nums.push_back(std::to_string(i));
		nums.push_back(std::to_string(i + 1));
		nums.push_back(std::to_string(i + 2));
		nwk << "((" << i << ',' << (i + 1) << ")," << (i + 2) << ')';
		nwk << (i < 3 * (magic - 2) ? ",(" : (i < 3 * (magic - 1) ? "," : ""));
	}
	for (unsigned i = 0; i < magic; ++i) {
		nwk << ')';
	}
	std::vector<const char*> names;
	for (const auto& str : nums) {
		names.push_back(str.c_str());
	}
	auto data = initializeMissingData(3 * magic + 1, magic, names.data());
	for (unsigned i = 0; i < 3 * magic; ++i) {
		setDataMatrix(data, i + 1, i / 3, 1);
	}
	for (unsigned i = 0; i < magic; ++i) {
		setDataMatrix(data, 0, i, 1);
	}

	REQUIRE(terraceAnalysis(data, nwk.str().c_str(), TA_COUNT, nullptr, result) ==
	        TERRACE_SPLIT_COUNT_OVERFLOW_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_num_partitions", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3"};
	auto data = initializeMissingData(3, 0, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 2, 1, 1);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3)", TA_DETECT, nullptr, result) ==
	        TERRACE_NUM_PARTITIONS_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_multifurcating", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 3, 0, 1);
	setDataMatrix(data, 3, 1, 1);
	REQUIRE(terraceAnalysis(data, "(s1,s2,s3,s4)", TA_DETECT, nullptr, result) ==
	        TERRACE_TREE_NOT_BINARY_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_no_root", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 3, 0, 1);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3,s4)", TA_DETECT, nullptr, result) ==
	        TERRACE_NO_ROOT_SPECIES_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_file", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 3, 0, 1);
	setDataMatrix(data, 3, 1, 1);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3,s4)", TA_ENUMERATE, "", result) ==
	        TERRACE_OUTPUT_FILE_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_empty_species", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3,s4)", TA_DETECT, nullptr, result) ==
	        TERRACE_SPECIES_WITHOUT_PARTITION_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_error_flag_conflict", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3,s4)", TA_DETECT | TA_COUNT, nullptr, result) ==
	        TERRACE_FLAG_CONFLICT_ERROR);
	REQUIRE(terraceAnalysis(data, "((s1,s2),s3,s4)", TA_ENUMERATE_COMPRESS, nullptr, result) ==
	        TERRACE_FLAG_CONFLICT_ERROR);
	freeMissingData(data);
}

TEST_CASE("c_api_non_terrace_data", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4"};
	auto data = initializeMissingData(4, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 1, 1, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 2, 1, 1);
	setDataMatrix(data, 3, 0, 1);
	setDataMatrix(data, 3, 1, 1);
	REQUIRE(terraceAnalysis(data, "(s1, (s2, (s3, s4)))", TA_DETECT, nullptr, result) ==
	        TERRACE_SUCCESS);
	CHECK(mpz_cmp_ui(result, 1) == 0);
	freeMissingData(data);
}

TEST_CASE("c_api_terrace_data", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4", "s5"};
	auto data = initializeMissingData(5, 2, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 0, 1);
	setDataMatrix(data, 2, 1, 1);
	setDataMatrix(data, 3, 1, 1);
	setDataMatrix(data, 4, 1, 1);
	REQUIRE(terraceAnalysis(data, "(s3, ((s1, s2), (s4, s5)))", TA_DETECT, nullptr, result) ==
	        TERRACE_SUCCESS);
	CHECK(mpz_cmp_ui(result, 1) > 0);
	freeMissingData(data);
}

TEST_CASE("c_api_count_data", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4", "s5", "s6"};
	auto data = initializeMissingData(6, 3, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 2, 1);
	setDataMatrix(data, 3, 2, 1);
	setDataMatrix(data, 3, 1, 1);
	setDataMatrix(data, 4, 0, 1);
	setDataMatrix(data, 4, 1, 1);
	setDataMatrix(data, 4, 2, 1);
	setDataMatrix(data, 5, 1, 1);
	setDataMatrix(data, 5, 2, 1);
	REQUIRE(terraceAnalysis(data, "((s4, (s3, (s2, (s1, s6)))), s5)", TA_COUNT, nullptr,
	                        result) == TERRACE_SUCCESS);
	CHECK(mpz_cmp_ui(result, 35) == 0);
	freeMissingData(data);
}

TEST_CASE("c_api_enumerate_compressed_data", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4", "s5", "s6"};
	auto data = initializeMissingData(6, 3, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 2, 1);
	setDataMatrix(data, 3, 2, 1);
	setDataMatrix(data, 3, 1, 1);
	setDataMatrix(data, 4, 0, 1);
	setDataMatrix(data, 4, 1, 1);
	setDataMatrix(data, 4, 2, 1);
	setDataMatrix(data, 5, 1, 1);
	setDataMatrix(data, 5, 2, 1);
	auto filename = "terraphast_capi_compressed_testfile.nwk";
	REQUIRE(terraceAnalysis(data, "((s4, (s3, (s2, (s1, s6)))), s5)",
	                        TA_COUNT | TA_ENUMERATE | TA_ENUMERATE_COMPRESS, filename,
	                        result) == TERRACE_SUCCESS);
	CHECK(mpz_cmp_ui(result, 35) == 0);
	std::stringstream buffer;
	{
		std::ifstream s(filename);
		buffer << s.rdbuf();
	}
	std::remove(filename);
	CHECK(buffer.str() == "(s5,(s2,((s3,s6),(s1,s4))|(s4,{s1,s3,s6})|((s4,(s3,s6)),s1))|((s3,"
	                      "s6),{s1,s2,s4})|({s2,s3,s6},(s1,s4))|(s4,{s1,s2,s3,s6})|((s2,s4),{"
	                      "s1,s3,s6})|((s4,(s3,s6)),(s1,s2))|(((s3,s6),(s2,s4))|(s4,{s2,s3,s6})"
	                      "|((s4,(s3,s6)),s2),s1))");
	freeMissingData(data);
}

TEST_CASE("c_api_enumerate_data", "[c-api]") {
	mpz_t result;
	mpz_init(result);
	const char* names[] = {"s1", "s2", "s3", "s4", "s5", "s6"};
	auto data = initializeMissingData(6, 3, names);
	setDataMatrix(data, 0, 0, 1);
	setDataMatrix(data, 1, 0, 1);
	setDataMatrix(data, 2, 2, 1);
	setDataMatrix(data, 3, 2, 1);
	setDataMatrix(data, 3, 1, 1);
	setDataMatrix(data, 4, 0, 1);
	setDataMatrix(data, 4, 1, 1);
	setDataMatrix(data, 4, 2, 1);
	setDataMatrix(data, 5, 1, 1);
	setDataMatrix(data, 5, 2, 1);
	auto filename = "terraphast_capi_testfile.nwk";
	REQUIRE(terraceAnalysis(data, "((s4, (s3, (s2, (s1, s6)))), s5)", TA_COUNT | TA_ENUMERATE,
	                        filename, result) == TERRACE_SUCCESS);
	CHECK(mpz_cmp_ui(result, 35) == 0);
	std::stringstream buffer;
	{
		std::ifstream s(filename);
		buffer << s.rdbuf();
	}
	std::remove(filename);
	CHECK(buffer.str() == "(s5,(s2,((s3,s6),(s1,s4))));\n"
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
	freeMissingData(data);
}