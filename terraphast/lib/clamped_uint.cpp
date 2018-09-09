#include <terraces/clamped_uint.hpp>

#include <ostream>
#include <terraces/errors.hpp>

#include "intrinsics.hpp"
#include "utils.hpp"

namespace terraces {

namespace {
static constexpr auto max_index = std::numeric_limits<index>::max();
}

template <bool except>
checked_uint<except>::checked_uint(index value) : m_value{value} {}

// clamped version
template <>
checked_uint<false>& checked_uint<false>::operator+=(checked_uint<false> other) {
	if (bits::add_overflow(m_value, other.m_value, m_value)) {
		m_value = max_index;
	}
	return *this;
}

template <>
checked_uint<false>& checked_uint<false>::operator*=(checked_uint<false> other) {
	if (bits::mul_overflow(m_value, other.m_value, m_value)) {
		m_value = max_index;
	}
	return *this;
}

template <>
checked_uint<true>& checked_uint<true>::operator+=(checked_uint<true> other) {
	utils::ensure<tree_count_overflow_error>(
	        !bits::add_overflow(m_value, other.m_value, m_value), "Addition overflowed");
	return *this;
}

template <>
checked_uint<true>& checked_uint<true>::operator*=(checked_uint<true> other) {
	utils::ensure<tree_count_overflow_error>(
	        !bits::mul_overflow(m_value, other.m_value, m_value), "Multiplication overflowed");
	return *this;
}

template <bool except>
bool checked_uint<except>::is_clamped() const {
	return m_value == max_index;
}

template <bool except>
index checked_uint<except>::value() const {
	return m_value;
}

template <bool except>
bool operator==(checked_uint<except> a, checked_uint<except> b) {
	return a.value() == b.value();
}

template <bool except>
bool operator!=(checked_uint<except> a, checked_uint<except> b) {
	return !(a == b);
}

template <bool except>
checked_uint<except> operator+(checked_uint<except> a, checked_uint<except> b) {
	return a += b;
}

template <bool except>
checked_uint<except> operator*(checked_uint<except> a, checked_uint<except> b) {
	return a *= b;
}

template <bool except>
std::ostream& operator<<(std::ostream& stream, checked_uint<except> val) {
	if (val.is_clamped()) {
		stream << ">= ";
	}
	stream << val.value();
	return stream;
}

// explicitly instantiate template class
template class checked_uint<false>;
template class checked_uint<true>;

// explicitly instantiate template functions
template bool operator==(checked_uint<false>, checked_uint<false>);
template bool operator!=(checked_uint<false>, checked_uint<false>);
template checked_uint<false> operator+(checked_uint<false>, checked_uint<false>);
template checked_uint<false> operator*(checked_uint<false>, checked_uint<false>);
template std::ostream& operator<<(std::ostream&, checked_uint<false>);

template bool operator==(checked_uint<true>, checked_uint<true>);
template bool operator!=(checked_uint<true>, checked_uint<true>);
template checked_uint<true> operator+(checked_uint<true>, checked_uint<true>);
template checked_uint<true> operator*(checked_uint<true>, checked_uint<true>);
template std::ostream& operator<<(std::ostream&, checked_uint<true>);

} // namespace terraces
