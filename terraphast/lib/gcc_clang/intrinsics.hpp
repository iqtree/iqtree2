#include <terraces/definitions.hpp>

namespace terraces {
namespace bits {

inline index popcount(index word) { return (index)__builtin_popcountll(word); }

inline index bitscan(index word) { return (index)__builtin_ctzll(word); }

inline index rbitscan(index word) { return (index)(63 - __builtin_clzll(word)); }

inline bool add_overflow(index a, index b, index& result) {
	return __builtin_add_overflow(a, b, &result);
}

inline bool mul_overflow(index a, index b, index& result) {
	return __builtin_mul_overflow(a, b, &result);
}

} // namespace bits
} // namespace terraces
