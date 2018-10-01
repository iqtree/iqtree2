// This file is only used to instantiate some configurations of tree_enumerators

#include "../lib/supertree_enumerator.hpp"
#include "../lib/supertree_variants.hpp"
#include "../lib/supertree_variants_debug.hpp"
#include "../lib/supertree_variants_multitree.hpp"
#include <terraces/bigint.hpp>
#include <terraces/clamped_uint.hpp>

namespace terraces {
using namespace terraces::variants;
using namespace terraces::debug::variants;

template class tree_enumerator<check_callback>;

#ifdef USE_GMP

template class tree_enumerator<count_callback<mpz_class>>;

template class tree_enumerator<count_callback<big_integer>>;

#else
template class tree_enumerator<count_callback<clamped_uint>>;
#endif

template class tree_enumerator<clamped_count_callback>;

template class tree_enumerator<multitree_callback>;

template class tree_enumerator<logging_decorator<check_callback>>;

#ifdef USE_GMP

template class tree_enumerator<logging_decorator<count_callback<mpz_class>>>;

template class tree_enumerator<logging_decorator<count_callback<big_integer>>>;

#else
template class tree_enumerator<logging_decorator<count_callback<clamped_uint>>>;
#endif

template class tree_enumerator<logging_decorator<clamped_count_callback>>;

template class tree_enumerator<logging_decorator<multitree_callback>>;

template class tree_enumerator<stack_state_decorator<check_callback>>;

#ifdef USE_GMP

template class tree_enumerator<stack_state_decorator<count_callback<mpz_class>>>;

template class tree_enumerator<stack_state_decorator<count_callback<big_integer>>>;

#else
template class tree_enumerator<stack_state_decorator<count_callback<clamped_uint>>>;
#endif

template class tree_enumerator<stack_state_decorator<clamped_count_callback>>;

template class tree_enumerator<stack_state_decorator<multitree_callback>>;
} // namespace terraces
