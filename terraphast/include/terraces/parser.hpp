#ifndef TERRACES_PARSER_HPP
#define TERRACES_PARSER_HPP

#include <istream>
#include <stdexcept>
#include <string>
#include <utility>

#include "bitmatrix.hpp"
#include "trees.hpp"

namespace terraces {

struct occurrence_data {
	bitmatrix matrix;
	name_map names;
	index_map indices;
};

struct named_tree {
	terraces::tree tree;
	name_map names;
	index_map indices;
};

/**
 * Parses a string in Newick format with given taxon IDs and returns the corresponding rooted tree.
 * \throws bad_input_error if the tree is malformed or an unknown taxon is encountered.
 * \returns the \ref tree corresponding to the input string
 */
tree parse_nwk(const std::string& input, const index_map& taxa);

/**
 * Parses a string in Newick format and returns the corresponding rooted tree and taxon names.
 * \throws bad_input_error if the tree is malformed or a duplicate taxon is encountered.
 * \returns the \ref named_tree corresponding to the input string
 */
named_tree parse_new_nwk(const std::string& input);

/**
 * Parses a data-file and returns the associated bit-matrix with taxon names as well as a
 * comprehensive taxon (or \ref none if none exists).
 * \throws bad_input_error if the data file is malformed or a duplicate taxon is encountered.
 * \returns the \ref occurrence_data corresponding to the input stream
 */
occurrence_data parse_bitmatrix(std::istream& input);

} // namespace terraces

#endif
