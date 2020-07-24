#ifndef TERRACES_CONSTRAINTS_HPP
#define TERRACES_CONSTRAINTS_HPP

#include <iosfwd>
#include <tuple>
#include <vector>

#include "trees.hpp"

namespace terraces {

/**
 * Represents a LCA constraint on leaves of a tree.
 * It is of the form \code
 * lca(left, shared) < lca(shared, right)
 * where the LCA's are compared by their height in the tree.
 */
struct constraint {
	index left;
	index shared;
	index right;

	constraint(index leftIndex, index sharedIndex, index rightIndex)
	        : left{leftIndex}, shared{sharedIndex}, right{rightIndex} {}

	bool operator==(const constraint& o) const {
		return std::tie(left, shared, right) == std::tie(o.left, o.shared, o.right);
	}

	bool operator!=(const constraint& o) const { return !(o == *this); }
};

using constraints = std::vector<constraint>;

std::ostream& operator<<(std::ostream& s, const constraint& c);

/**
 * Extracts all LCA constraints from the input trees.
 * \param trees The input trees.
 * \returns All LCA constraints for the input trees.
 * For every inner edge, one LCA constraint based on the leftmost and rightmost descendant
 * of the endpoints are is extracted.
 */
constraints compute_constraints(const std::vector<tree>& trees);

/**
 * Removes duplicate constraints from the input vector
 * \param in_c The input constraints.
 * \returns The input constraints without duplicates and possibly reordered.
 */
index deduplicate_constraints(constraints& in_c);

} // namespace terraces

#endif // TERRACES_CONSTRAINTS_HPP
