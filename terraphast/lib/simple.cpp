#include <terraces/advanced.hpp>
#include <terraces/errors.hpp>
#include <terraces/parser.hpp>
#include <terraces/rooting.hpp>
#include <terraces/simple.hpp>

#include <fstream>
#include <sstream>

#include "io_utils.hpp"
#include "utils.hpp"

namespace terraces {
namespace simple {

using utils::open_ifstream;
using utils::read_file_full;
using utils::read_ifstream_full;

namespace {

std::pair<supertree_data, name_map> parse_data(const std::string& nwk_string,
                                               std::istream& matrix_stream, bool force) {
	auto occ_data = parse_bitmatrix(matrix_stream);
	auto tree = parse_nwk(nwk_string, occ_data.indices);
	if (force) {
		occ_data.matrix = maximum_comprehensive_columnset(occ_data.matrix);
	}
	auto data = create_supertree_data(tree, occ_data.matrix);
	return {data, occ_data.names};
}

} // anonymous namespace

bool is_on_terrace(const std::string& nwk_string, std::istream& matrix_stream, bool force) {
	return check_terrace(parse_data(nwk_string, matrix_stream, force).first);
}

bool is_on_terrace(std::istream& nwk_stream, std::istream& matrix_stream, bool force) {
	return is_on_terrace(read_ifstream_full(nwk_stream), matrix_stream, force);
}

bool is_on_terrace(std::istream& nwk_stream, const std::string& matrix_string, bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return is_on_terrace(read_ifstream_full(nwk_stream), matrix_stream, force);
}

bool is_on_terrace(const std::string& nwk_string, const std::string& matrix_string, bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return is_on_terrace(nwk_string, matrix_stream, force);
}

bool is_on_terrace_from_file(const std::string& nwk_filename, const std::string& matrix_filename,
                             bool force) {
	auto nwk_string = read_file_full(nwk_filename);
	auto matrix_stream = open_ifstream(matrix_filename);
	return is_on_terrace(nwk_string, matrix_stream, force);
}

index get_terrace_size(const std::string& nwk_string, std::istream& matrix_stream, bool force) {
	return count_terrace(parse_data(nwk_string, matrix_stream, force).first);
}

index get_terrace_size(std::istream& nwk_stream, std::istream& matrix_stream, bool force) {
	return get_terrace_size(read_ifstream_full(nwk_stream), matrix_stream, force);
}

index get_terrace_size(std::istream& nwk_stream, const std::string& matrix_string, bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return get_terrace_size(read_ifstream_full(nwk_stream), matrix_stream, force);
}

index get_terrace_size(const std::string& nwk_string, const std::string& matrix_string,
                       bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return get_terrace_size(nwk_string, matrix_stream, force);
}
index get_terrace_size_from_file(const std::string& nwk_filename,
                                 const std::string& matrix_filename, bool force) {
	auto nwk_string = read_file_full(nwk_filename);
	auto matrix_stream = open_ifstream(matrix_filename);
	return get_terrace_size(nwk_string, matrix_stream, force);
}

big_integer get_terrace_size_bigint(const std::string& nwk_string, std::istream& matrix_stream,
                                    bool force) {
	return count_terrace_bigint(parse_data(nwk_string, matrix_stream, force).first);
}

big_integer get_terrace_size_bigint(std::istream& nwk_stream, std::istream& matrix_stream,
                                    bool force) {
	return get_terrace_size_bigint(read_ifstream_full(nwk_stream), matrix_stream, force);
}

big_integer get_terrace_size_bigint(std::istream& nwk_stream, const std::string& matrix_string,
                                    bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return get_terrace_size_bigint(read_ifstream_full(nwk_stream), matrix_stream, force);
}

big_integer get_terrace_size_bigint(const std::string& nwk_string, const std::string& matrix_string,
                                    bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return get_terrace_size_bigint(nwk_string, matrix_stream, force);
}
big_integer get_terrace_size_bigint_from_file(const std::string& nwk_filename,
                                              const std::string& matrix_filename, bool force) {
	auto nwk_string = read_file_full(nwk_filename);
	auto matrix_stream = open_ifstream(matrix_filename);
	return get_terrace_size_bigint(nwk_string, matrix_stream, force);
}

big_integer print_terrace(const std::string& nwk_string, std::istream& matrix_stream,
                          std::ostream& output, bool force) {

	auto data = parse_data(nwk_string, matrix_stream, force);
	return print_terrace(data.first, data.second, output);
}

big_integer print_terrace(std::istream& nwk_stream, const std::string& matrix_string,
                          std::ostream& out, bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return print_terrace(read_ifstream_full(nwk_stream), matrix_stream, out, force);
}

big_integer print_terrace(std::istream& nwk_stream, std::istream& matrix_stream, std::ostream& out,
                          bool force) {
	return print_terrace(read_ifstream_full(nwk_stream), matrix_stream, out, force);
}

big_integer print_terrace(const std::string& nwk_string, const std::string& matrix_string,
                          std::ostream& output, bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return print_terrace(nwk_string, matrix_stream, output, force);
}

big_integer print_terrace_from_file(const std::string& nwk_filename,
                                    const std::string& matrix_filename, std::ostream& output,
                                    bool force) {
	auto nwk_string = read_file_full(nwk_filename);
	auto matrix_stream = open_ifstream(matrix_filename);
	return print_terrace(nwk_string, matrix_stream, output, force);
}

big_integer print_terrace_compressed(const std::string& nwk_string, std::istream& matrix_stream,
                                     std::ostream& output, bool force) {

	auto data = parse_data(nwk_string, matrix_stream, force);
	return print_terrace_compressed(data.first, data.second, output);
}

big_integer print_terrace_compressed(std::istream& nwk_stream, const std::string& matrix_string,
                                     std::ostream& out, bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return print_terrace_compressed(read_ifstream_full(nwk_stream), matrix_stream, out, force);
}

big_integer print_terrace_compressed(std::istream& nwk_stream, std::istream& matrix_stream,
                                     std::ostream& out, bool force) {
	return print_terrace_compressed(read_ifstream_full(nwk_stream), matrix_stream, out, force);
}

big_integer print_terrace_compressed(const std::string& nwk_string,
                                     const std::string& matrix_string, std::ostream& output,
                                     bool force) {
	auto matrix_stream = std::istringstream{matrix_string};
	return print_terrace_compressed(nwk_string, matrix_stream, output, force);
}

big_integer print_terrace_compressed_from_file(const std::string& nwk_filename,
                                               const std::string& matrix_filename,
                                               std::ostream& output, bool force) {
	auto nwk_string = read_file_full(nwk_filename);
	auto matrix_stream = open_ifstream(matrix_filename);
	return print_terrace_compressed(nwk_string, matrix_stream, output, force);
}

} // namespace simple
} // namespace terraces
