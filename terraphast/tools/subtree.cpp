#include <fstream>
#include <iostream>
#include <iterator>

#include <set>

#include <terraces/parser.hpp>
#include <terraces/rooting.hpp>
#include <terraces/subtree_extraction.hpp>
#include <terraces/trees.hpp>

int main(int argc, char** argv) try {
	if (argc != 3) {
		std::cerr << "Usage: " << argv[0] << " <tree-file> <species-file>" << std::endl;
		return 1;
	}
	auto tree_file = std::ifstream{argv[1]};
	auto tree_string = std::string{};
	std::getline(tree_file, tree_string);
	auto data = terraces::parse_new_nwk(tree_string);

	auto names_file = std::ifstream{argv[2]};
	std::set<std::string> names;
	std::copy(std::istream_iterator<std::string>(names_file),
	          std::istream_iterator<std::string>(), std::inserter(names, names.begin()));
	terraces::bitmatrix occ{data.tree.size(), 1};
	for (auto& species_pair : data.indices) {
		if (names.find(species_pair.first) != names.end()) {
			occ.set(species_pair.second, 0, true);
		}
	}
	terraces::reroot_at_taxon_inplace(data.tree, data.indices[*names.begin()]);
	auto subtree = terraces::subtrees(data.tree, occ)[0];
	std::cout << as_newick(subtree, data.names);
} catch (std::exception& e) {
	std::cerr << "Error: " << e.what() << "\n";
}
