#include "../lib/io_utils.hpp"
#include "../lib/trees_impl.hpp"

#include <iostream>
#include <terraces/parser.hpp>

int main(int argc, char** argv) {
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " [treefile] [rooted]";
		return -1;
	}
	auto tree = terraces::parse_new_nwk(terraces::utils::read_file_full(argv[1]));
	bool rooted = bool(std::stoi(argv[2]));
	terraces::print_tree_dot(tree.tree, tree.names, std::cout, rooted);
}