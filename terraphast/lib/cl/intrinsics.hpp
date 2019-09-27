#include <terraces/definitions.hpp>

#ifdef _WIN64
#pragma intrinsic(_BitScanForward64, _BitScanReverse64, __popcnt64)
#else
#pragma intrinsic(_BitScanForward, _BitScanReverse, __popcnt)
#endif

namespace terraces {
namespace bits {

#ifdef _WIN64
inline index popcount(index word) { return (index)__popcnt64(word); }

inline index bitscan(index word) {
	unsigned long idx;
	_BitScanForward64(&idx, word);
	return (index)idx;
}

inline index rbitscan(index word) {
	unsigned long idx;
	_BitScanReverse64(&idx, word);
	return (index)idx;
}
#else
inline index popcount(index word) { return index(__popcnt(word)); }

inline index bitscan(index word) {
	unsigned long idx;
	_BitScanForward(&idx, word);
	return index(idx);
}

inline index rbitscan(index word) {
	unsigned long idx;
	_BitScanReverse(&idx, word);
	return index(idx);
}
#endif

namespace {
constexpr index max_index = std::numeric_limits<index>::max();
}

inline bool add_overflow(index a, index b, index& result) {
	result = a + b;
	if (max_index - b < a) {
		return true;
	} else {
		return false;
	}
}

inline bool mul_overflow(index a, index b, index& result) {
	result = a * b;
	if (max_index / b < a) {
		return true;
	} else {
		return false;
	}
}

} // namespace bits
} // namespace terraces
