#ifndef ADVANCED_HPP
#define ADVANCED_HPP

#include <functional>
#include <iosfwd>

#include "bigint.hpp"
#include "bitmatrix.hpp"
#include "constraints.hpp"
#include "trees.hpp"

namespace terraces {

/**
 * Contains the necessary data to describe all possible supertrees equivalent to a certain tree.
 */
struct supertree_data {
	/** The lca constraints describing the structure of the possible supertrees */
	terraces::constraints constraints;
	/** The number of leaves of the supertree. */
	index num_leaves;
	/** The index of the 'root leaf' in leaf-based numbering.
	 * It will be placed as a leaf below the root of all supertrees.
	 * In terms of phylogenetic trees, this means the index of a comprehensive taxon. */
	index root;
};

/**
 * Returns the index of the first comprehensive taxon in a occurrence bitmatrix.
 * @param data The occurrence bitmatrix.
 * @return The (row) index of the first comprehensive taxon in the input matrix
 *         or @ref none if none exists.
 */
index find_comprehensive_taxon(const bitmatrix& data);

/**
 * Extracts a maximum size subset of the columns
 * so the resulting matrix contains a comprehensive taxon.
 * @param data The occurrence matrix
 * @return The output matrix containing a subset of the columns of the input matrix.
 */
bitmatrix maximum_comprehensive_columnset(const bitmatrix& data);

/**
 * Computes the necessary data to enumerate the supertrees of the given tree and missing data
 * matrix.
 * \param tree The phylogenetic tree.
 * \param data The missing data matrix. It must contain only data for the leaves!
 * \returns \ref supertree_data object describing all possible supertrees equivalent to the input
 * tree.
 */
supertree_data create_supertree_data(const tree& tree, const bitmatrix& data);

/**
 * Checks if a phylogenetic tree lies on a phylogenetic terrace.
 * \param data The constraints extracted from the tree and missing data matrix describing all
 * possible supertrees.
 * \return True if and only if the tree lies on a phylogenetic terrace,
 *         i.e. there are at least 2 different trees (called supertrees)
 *         that are equivalent with respect to the missing data matrix.
 */
bool check_terrace(const supertree_data& data);

/**
 * Checks if a phylogenetic tree lies on a phylogenetic terrace by computing a lower bound to the
 * terrace size.
 * \param data The constraints extracted from the tree and missing data matrix describing all
 * possible supertrees.
 * \return A lower bound to the number of trees on the phylogenetic terrace.
 */
index fast_count_terrace(const supertree_data& data);

/**
 * Counts all trees on a terrace around a phylogenetic tree.
 * Note that this number might not be representable as a 32/64 bit integer, and might thus be
 * clamped.
 * \param data The constraints extracted from the tree and missing data matrix describing all
 * possible supertrees.
 * \return The number of trees on the phylogenetic terrace containing the input tree.
 *         Note that if this result is UINT32/64_MAX = 2^32/64 - 1, the computations resulted in an
 * overflow,
 *         i.e. the result is only a lower bound on the number of trees on this terrace.
 */
index count_terrace(const supertree_data& data);

/**
 * Counts all trees on a terrace around a phylogenetic tree.
 * \param data The constraints extracted from the tree and missing data matrix describing all
 * possible supertrees.
 * \return The number of trees on the phylogenetic terrace containing the input tree.
 */
big_integer count_terrace_bigint(const supertree_data& data);

/**
 * Enumerates all trees on a terrace around a phylogenetic tree.
 * The trees will be printed in a compressed <b>multitree format</b>,
 * which is an extension of the normal Newick format:
 * In place of a subtree <i>T</i>, there may be two additional types of nodes:
 * <ul>
 * <li>An <b>alternative list</b> of all possible sub(-multi-)trees that appear on the terrace:<br>
 *     <i>T1</i>|<i>T2</i>|<i>T3</i>|...</li>
 * <li>An <b>unconstrained subtree</b> where all possible subtrees with the given leaves are part of
 * the terrace:<br>
 * {<i>leaf1</i>,<i>leaf2</i>,<i>leaf3</i>,...}
 * </li>
 * </ul>
 * \param data The constraints extracted from the tree and missing data matrix describing all
 * possible supertrees.
 * \param names The name map containing only leaf names. It will be used to output the multitree.
 * \param output The output stream into which the multitree will be written.
 * \return The number of trees on the phylogenetic terrace containing the input tree.
 */
big_integer print_terrace_compressed(const supertree_data& data, const name_map& names,
                                     std::ostream& output);

/**
 * Enumerates all trees on a terrace around a phylogenetic tree.
 * The trees will be printed in Newick format, one tree per line
 * \param data The constraints extracted from the tree and missing data matrix describing all
 * possible supertrees.
 * \param names The name map containing only leaf names. It will be used to output the multitree.
 * \param output The output stream into which the trees will be written.
 * \return The number of trees on the phylogenetic terrace containing the input tree.
 */
big_integer print_terrace(const supertree_data& data, const name_map& names, std::ostream& output);

/**
 * Enumerates all trees on a terrace around a phylogenetic tree.
 * The given callback function will be called with every tree on the terrace as a parameter.
 * \param data The constraints extracted from the tree and missing data matrix describing all
 * possible supertrees.
 * \param callback The callback function taking a tree as a parameter.
 */
void enumerate_terrace(const supertree_data& data, std::function<void(const tree&)> callback);

} // namespace terraces

#endif // ADVANCED_HPP
