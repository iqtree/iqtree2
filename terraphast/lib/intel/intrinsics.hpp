#include <terraces/definitions.hpp>

namespace terraces {
namespace bits {

#if defined(__LP64__) || defined(_LP64)
inline index popcount(index word) { return index(_popcnt64(word)); }

inline index bitscan(index word) {
    unsigned int idx;
    _BitScanForward64(&idx, word);
    return index(idx);
}

inline index rbitscan(index word) {
    unsigned int idx;
    _BitScanReverse64(&idx, word);
    return index(idx);
}
#else
inline index popcount(index word) { return index(_popcnt(word)); }

inline index bitscan(index word) {
    unsigned int idx;
    _BitScanForward(&idx, word);
    return index(idx);
}

inline index rbitscan(index word) {
    unsigned int idx;
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
