
#ifndef TERRACES_ERRORS_HPP
#define TERRACES_ERRORS_HPP

#include <stdexcept>
#include <string>

namespace terraces {

/** Error codes for different input errors. */
enum class bad_input_error_type {
	/** Mismatched quotes in a nwk string. */
	nwk_mismatched_quotes,
	/** Mismatched parentheses in a nwk string. */
	nwk_mismatched_parentheses,
	/** Unknown taxon when parsing a nwk string using an \ref index_map. */
	nwk_taxon_unknown,
	/** Duplicate taxon in a nwk string. */
	nwk_taxon_duplicate,
	/** Multifurcating nwk tree. */
	nwk_multifurcating,
	/** Otherwise malformed nwk tree. */
	nwk_malformed,
	/** Less than 4 taxa in the tree. */
	nwk_tree_trivial,
	/** Duplicate taxon in a bitmatrix. */
	bitmatrix_name_duplicate,
	/** Empty taxon name in a bitmatrix. */
	bitmatrix_name_empty,
	/** Mismatching number of rows/columns between header and content. */
	bitmatrix_size_invalid,
	/** Otherwise malformed bitmatrix. */
	bitmatrix_malformed,
	/** Mismatching size between tree and bitmatrix. */
	tree_mismatching_size,
	/** Unnamed leaf found in a tree. */
	tree_unnamed_leaf,
};

/** This error is thrown if the input to a function is malformed. */
class bad_input_error : public std::runtime_error {
public:
	explicit bad_input_error(bad_input_error_type type);
	explicit bad_input_error(bad_input_error_type type, std::string msg);
	bad_input_error_type type() const { return m_type; }

private:
	bad_input_error_type m_type;
};

/** This error is thrown if there is no comprehensive taxon in a dataset. */
class no_usable_root_error : public std::runtime_error {
	using std::runtime_error::runtime_error;
};

/** This error is thrown if a file couldn't be opened. */
class file_open_error : public std::runtime_error {
	using std::runtime_error::runtime_error;
};

/**
 * This error is thrown if the terrace search ran into a huge terrace.
 * More specifically, this means that on a recursive subcalls,
 * more than 2^64 possible bipartitions would have to be enumerated.
 * This would take quite a lot of time.
 * It will also be thrown during a computation overflow
 * when GMP integration is disabled.
 */
class tree_count_overflow_error : public std::overflow_error {
	using std::overflow_error::overflow_error;
};

} // namespace terraces

#endif
