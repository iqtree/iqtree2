#include <terraces/bigint.hpp>

#include <ostream>

#ifdef USE_GMP
namespace terraces {
big_integer::big_integer(index i) : m_value{i} {}
big_integer& big_integer::operator+=(const big_integer& other) {
	m_value += other.m_value;
	return *this;
}
big_integer& big_integer::operator*=(const big_integer& other) {
	m_value *= other.m_value;
	return *this;
}
bool big_integer::is_clamped() const { return false; }
const mpz_class& big_integer::value() const { return m_value; }

big_integer operator+(const big_integer& a, const big_integer& b) {
	big_integer result = a;
	result += b;
	return result;
}

big_integer operator*(const big_integer& a, const big_integer& b) {
	big_integer result = a;
	result *= b;
	return result;
}

bool operator==(const big_integer& a, const big_integer& b) { return a.value() == b.value(); }

bool operator!=(const big_integer& a, const big_integer& b) { return !(a == b); }

std::ostream& operator<<(std::ostream& stream, const big_integer& val) {
	return stream << val.value();
}

} // namespace terraces
#endif // USE_GMP
