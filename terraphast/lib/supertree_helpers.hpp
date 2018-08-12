#ifndef SUPERTREE_HELPERS_HPP
#define SUPERTREE_HELPERS_HPP

#include <terraces/constraints.hpp>
#include <terraces/trees.hpp>

#include "ranked_bitvector.hpp"
#include "stack_allocator.hpp"
#include "union_find.hpp"

namespace terraces {

/**
 * Filters the given constraints using the given leaves.
 * \param leaves The leaf set.
 * \param c_occ The set containing all constraint indices that are in the constraint set.
 * \param c The constraints themselves.
 * \param a The allocator used to construct the result bitvector.
 * \returns The bitvector representing the subset of \p c_occ that contains all constraints
 *          whose entries are all within \p leaves.
 */
bitvector filter_constraints(const ranked_bitvector& leaves, const bitvector& c_occ,
                             const constraints& c, utils::stack_allocator<index> a);

/**
 * Applies the given constraints to the given leaves to
 * \param leaves The leaf set.
 * \param c_occ The indices of all constraints to be applied.
 *              They must have been filtered with \p leaves
 *              using \ref filter_constraints.
 * \param c The constraints themselves.
 * \param a The allocator used to construct the result.
 * \returns A \ref union_find representation of the sets created by applying the constraints.
 */
union_find apply_constraints(const ranked_bitvector& leaves, const bitvector& c_occ,
                             const constraints& c, utils::stack_allocator<index> a);

} // namespace terraces

#endif // SUPERTREE_HELPERS_HPP
