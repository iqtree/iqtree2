#ifndef SUPERTREE_VARIANTS_HPP
#define SUPERTREE_VARIANTS_HPP

#include <terraces/constraints.hpp>

#include <terraces/bigint.hpp>
#include <terraces/clamped_uint.hpp>

#include "bipartitions.hpp"
#include "trees_impl.hpp"

namespace terraces {
namespace variants {

template <typename Result>
/**
 * An abstract callback used to control the execution of \ref tree_enumerator.
 * Every callback must declare a return type whose values are aggregated
 * (sum over different bipartitions) and combined (results from left and right subtree).
 * <p>Additionally, the enumeration can be terminated early using \ref fast_return
 * or \ref continue_iteration, if one does not want to iterate over the whole recursion tree.</p>
 */
class abstract_callback {
public:
	/** The result type returned by the @ref tree_enumerator. */
	using result_type = Result;

	/**
	 * Called when a (sub)call begins.
	 * \param leaves The leaf set inspected at the current recursive call.
	 */
	void enter(const ranked_bitvector& leaves) { (void)leaves; }
	/**
	 * Called when a (sub)call finishes.
	 * \param val The result that was either returned from \ref fast_return_value
	 *            or accumulated from the recursive subcalls.
	 * \returns The result that should be passed to the calling function.
	 */
	Result exit(Result val) { return val; }

	/** Returns the result for a single leaf. */
	Result base_one_leaf(index);
	/** Returns the result for two leaves. */
	Result base_two_leaves(index, index);
	/** Returns the result for multiple leaves without constraints. */
	Result base_unconstrained(const ranked_bitvector&);
	/** Returns an empty result. */
	Result null_result() const;

	/**
	 * Called before iterating over the bipartitions to check if we can return early.
	 * This may be used to avoid searching subtrees that are 'not interesting'.
	 * \see fast_return_value
	 * \returns true if and only if we want to skip iterating over the current bipartitions.
	 *          (default: false)
	 */
	bool fast_return(const bipartitions&) { return false; }
	/** Returns the result to be returned in case \ref fast_return is true. */
	Result fast_return_value(const bipartitions&);

	/**
	 * Called when we begin iterating over the bipartitions.
	 * \param
	 * \returns The initial value for the result used for accumulation
	 *          - something equivalent to 0.
	 *          (default: value-initialization of @p Result)
	 */
	Result begin_iteration(const bipartitions& bip_it, const bitvector& c_occ,
	                       const constraints& c) {
		(void)bip_it;
		(void)c_occ;
		(void)c;
		return Result{};
	}
	/**
	 * Returns true iff the iteration should continue.
	 * This may be used to stop iterating after a sufficient number of results have been
	 * accumulated.
	 * \returns true if and only if we want to continue iterating over
	 */
	bool continue_iteration(Result) { return true; }
	/** Called when an iteration step begins. */
	void step_iteration(const bipartitions&, index) {}
	/** Called when the last iteration step has finished. */
	void finish_iteration() {}
	/** Called before descending into the left subset. */
	void left_subcall() {}
	/** Called before descending into the right subset. */
	void right_subcall() {}
	/**
	 * Accumulates the results from multiple bipartitions.
	 * \param accumulator The current value for the accumulator.
	 * \param value The value to be added to the accumulator.
	 * \returns The new accumulator value.
	 */
	Result accumulate(Result accumulator, Result value);
	/**
	 * Combines the results from two subcalls.
	 * \param left The result from the left subcall.
	 * \param right The result from the right subcall.
	 * \returns The combined result.
	 */
	Result combine(Result left, Result right);
};

template <typename Number>
/**
 * A callback implementation that counts all trees using the \p Number type.
 * It gives correct results if the number type is integer and doesn't overflow.
 */
class count_callback : public abstract_callback<Number> {
public:
	using return_type = typename abstract_callback<Number>::result_type;
	// only one choice for a single leaf
	return_type base_one_leaf(index) { return 1; }
	// only one choice for two leaves
	return_type base_two_leaves(index, index) { return 1; }
	// (#unrooted trees) choices for leaves with no constraints
	return_type base_unconstrained(const ranked_bitvector& leaves) {
		return count_unrooted_trees<return_type>(leaves.count());
	}
	return_type null_result() const { return 0; }

	// The number of bipartitions gives a lower bound on the number of trees.
	index fast_return_value(const bipartitions& bip_it) { return bip_it.num_bip(); }

	// Multiple choices are counted independently
	return_type accumulate(return_type acc, return_type val) { return acc + val; }
	// Choices from two the subtrees can be combined in any possible way
	return_type combine(return_type left, return_type right) { return left * right; }
};

/**
 * A callback implementation that counts all trees using \ref clamped_uint.
 * Thus, the result will be clamped at the maximal value for uint64_t,
 * eliminating the need to fully walk through the recursion.
 */
class clamped_count_callback : public count_callback<clamped_uint> {
public:
	using return_type = typename count_callback<clamped_uint>::result_type;
	// No need to keep counting if we already overflowed.
	bool continue_iteration(return_type result) { return !result.is_clamped(); }
};

/**
 * A callback implementation that returns a simple lower bound to the number
 * of trees compatible with the given constraints.
 * This allows to quickly check whether a tree lies on a phylogenetic terrace.
 * It will always run faster than callbacks that enumerate all possible trees.
 */
class check_callback : public abstract_callback<index> {
public:
	using return_type = index;

	// No choices for a single leaf
	index base_one_leaf(index) { return 1; }
	// No choices for two leaves
	index base_two_leaves(index, index) { return 1; }
	// There are at least two possible trees if we have at least three unconstrained leaves
	index base_unconstrained(const ranked_bitvector&) { return 2; }
	return_type null_result() const { return 0; }

	/* Since we assume there are no incompatible constraints,
	 * every subcall will return at least 1, so the number of bipartitions
	 * gives a simple lower bound for the number of trees. */
	bool fast_return(const bipartitions& bip_it) { return bip_it.num_bip() > 1; }
	index fast_return_value(const bipartitions& bip_it) { return bip_it.num_bip(); }

	/* No need to keep on counting if we already know we are on a terrace. */
	bool continue_iteration(index acc) { return acc < 2; }

	index accumulate(index acc, index val) { return acc + val; }
	index combine(index left, index right) { return left * right; }
};

} // namespace variants
} // namespace terraces

#endif // SUPERTREE_VARIANTS_HPP
