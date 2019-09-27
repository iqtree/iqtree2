#ifndef CONSTRAINTS_IMPL_HPP
#define CONSTRAINTS_IMPL_HPP

#include <terraces/constraints.hpp>

namespace terraces {

namespace utils {
template <class T1, class T2>
struct named_output;
}

std::ostream& operator<<(std::ostream& stream, utils::named_output<constraints, name_map> output);

} // namespace terraces

#endif // CONSTRAINTS_IMPL_HPP
