#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <limits>

namespace terraces {

/** Our index type. */
using index = std::size_t;

/** Index value used to represent a non-index. */
#ifdef CLANG_UNDER_VS
#undef max
#endif
constexpr static auto none = std::numeric_limits<index>::max();

} // namespace terraces

#endif // DEFINITIONS_HPP
