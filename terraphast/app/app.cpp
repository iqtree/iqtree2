#include <bitset>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>

#include <terraces/simple.hpp>

int main(int argc, char** argv) try {
	auto tree_file_name = std::string{};
	auto data_file_name = std::string{};
	if (argc == 3) {
		tree_file_name = argv[1];
		data_file_name = argv[2];
	} else {
		std::cerr << "Usage: \n" << argv[0] << " <tree-file> <occurrence file>\n";
		return 1;
	}
	auto trees = std::ostringstream{};
	const auto terraces_count =
	        terraces::simple::print_terrace_from_file(tree_file_name, data_file_name, trees);

	std::cout << "There are " << terraces_count
	          << " trees on the terrace.\n\nThe trees in question are:\n"
	          << trees.str() << '\n';
} catch (std::exception& e) {
	std::cerr << "Error: " << e.what() << "\n";
}
