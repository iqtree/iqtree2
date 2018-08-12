#include <fstream>
#include <iostream>

#include <terraces/advanced.hpp>
#include <terraces/parser.hpp>
#include <terraces/rooting.hpp>
#include <terraces/trees.hpp>

int main(int argc, char** argv) try {
	if (argc != 3) {
		std::cerr << "Usage: " << argv[0] << " <tree-file> <occurrence file>" << std::endl;
		return 1;
	}
	auto data_file = std::ifstream{argv[2]};
	const auto data = terraces::parse_bitmatrix(data_file);

	auto tree_file = std::ifstream{argv[1]};
	auto tree_string = std::string{};
	std::getline(tree_file, tree_string);
	auto tree = terraces::parse_nwk(tree_string, data.indices);
	auto comp_taxon = terraces::find_comprehensive_taxon(data.matrix);

	if (comp_taxon == terraces::none) {
		std::cerr << "Cannot find a comprehensive taxon for the tree\n";
		return 1;
	}

	terraces::reroot_at_taxon_inplace(tree, comp_taxon);

	std::cout << as_newick(tree, data.names);
} catch (std::exception& e) {
	std::cerr << "Error: " << e.what() << "\n";
}
