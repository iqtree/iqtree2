#include <terraces/errors.hpp>

namespace terraces {

namespace {

std::string build_error_message(bad_input_error_type type) {
	switch (type) {
	case bad_input_error_type::nwk_mismatched_quotes:
		return "Mismatching quotes in nwk tree";
	case bad_input_error_type::nwk_mismatched_parentheses:
		return "Mismatching parentheses in nwk tree";
	case bad_input_error_type::nwk_taxon_unknown:
		return "Unknown taxon in nwk tree";
	case bad_input_error_type::nwk_taxon_duplicate:
		return "Duplicate taxon in nwk tree";
	case bad_input_error_type::nwk_multifurcating:
		return "Only bifurcating trees are supported";
	case bad_input_error_type::nwk_malformed:
		return "Malformed nwk tree";
	case bad_input_error_type::nwk_tree_trivial:
		return "Less than 4 taxa in nwk tree";
	case bad_input_error_type::bitmatrix_name_duplicate:
		return "Duplicate taxon in bitmatrix ";
	case bad_input_error_type::bitmatrix_name_empty:
		return "Empty taxon name in bitmatrix";
	case bad_input_error_type::bitmatrix_size_invalid:
		return "Mismatching number of rows/columns "
		       "between bitmatrix header and content";
	case bad_input_error_type::bitmatrix_malformed:
		return "Malformed bitmatrix";
	case bad_input_error_type::tree_mismatching_size:
		return "Mismatching size between tree and bitmatrix";
	case bad_input_error_type::tree_unnamed_leaf:
		return "Unnamed leaf found in tree";
	}
	return "Unknown error";
}

} // anonymous namespace

bad_input_error::bad_input_error(bad_input_error_type type)
        : std::runtime_error{build_error_message(type)}, m_type{type} {}

bad_input_error::bad_input_error(bad_input_error_type type, std::string msg)
        : std::runtime_error{build_error_message(type) + "\n" + msg}, m_type{type} {}

} // namespace terraces
