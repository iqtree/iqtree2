#ifndef TERRACES_CLAMPED_UINT_HPP
#define TERRACES_CLAMPED_UINT_HPP

#include "definitions.hpp"

#include <iosfwd>

namespace terraces {

template <bool except>
class checked_uint {
	index m_value;

public:
	checked_uint(index value = 0);

	checked_uint<except>& operator+=(checked_uint<except> other);
	checked_uint<except>& operator*=(checked_uint<except> other);
	bool is_clamped() const;
	index value() const;
};

template <bool except>
bool operator==(checked_uint<except> a, checked_uint<except> b);

template <bool except>
bool operator!=(checked_uint<except> a, checked_uint<except> b);

template <bool except>
checked_uint<except> operator+(checked_uint<except> a, checked_uint<except> b);

template <bool except>
checked_uint<except> operator*(checked_uint<except> a, checked_uint<except> b);

template <bool except>
std::ostream& operator<<(std::ostream& stream, checked_uint<except> val);

using clamped_uint = checked_uint<false>;
using overflow_except_uint = checked_uint<true>;

} // namespace terraces

#endif // TERRACES_CLAMPED_UINT_HPP
