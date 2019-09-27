#ifndef TERRACES_BIGINT_HPP
#define TERRACES_BIGINT_HPP

#include "definitions.hpp"

#ifndef USE_GMP
#include "clamped_uint.hpp"
namespace terraces {
using big_integer = terraces::overflow_except_uint;
}
#else
#include <gmpxx.h>
#include <iosfwd>
namespace terraces {
class big_integer {
	mpz_class m_value;

public:
	big_integer(index i = 0);
	big_integer& operator+=(const big_integer& other);
	big_integer& operator*=(const big_integer& other);
	bool is_clamped() const;
	const mpz_class& value() const;
};
bool operator==(const big_integer& a, const big_integer& b);
bool operator!=(const big_integer& a, const big_integer& b);
big_integer operator+(const big_integer& a, const big_integer& b);
big_integer operator*(const big_integer& a, const big_integer& b);
std::ostream& operator<<(std::ostream& stream, const big_integer& val);
} // namespace terraces
#endif

#endif // TERRACES_BIGINT_HPP
