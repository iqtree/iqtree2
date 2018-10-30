#include <terraces/definitions.hpp>

namespace terraces {
namespace bits {

inline index popcount(index word) { return index(__builtin_popcountll(word)); }

inline index bitscan(index word) { return index(__builtin_ctzll(word)); }

static_assert(sizeof(long long) >= sizeof(index), "intrinsic word sizes incompatible");

inline index rbitscan(index word) {
	return index(std::numeric_limits<long long>::digits - __builtin_clzll(word));
}

inline bool add_overflow(index a, index b, index& result) {
	return __builtin_add_overflow(a, b, &result);
}

inline bool mul_overflow(index a, index b, index& result) {
	return __builtin_mul_overflow(a, b, &result);
}

} // namespace bits
} // namespace terraces
