//***************************************************************************
//*   Copyright (C) 2020 by James Barbetti <james_barbettii@yahoo.com>      *
//*                                                                         *
//*   This program is free software; you can redistribute it and/or modify  *
//*   it under the terms of the GNU General Public License as published by  *
//*   the Free Software Foundation; either version 2 of the License, or     *
//*   (at your option) any later version.                                   *
//*                                                                         *
//*   This program is distributed in the hope that it will be useful,       *
//*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
//*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
//*   GNU General Public License for more details.                          *
//*                                                                         *
//*   You should have received a copy of the GNU General Public License     *
//*   along with this program; if not, write to the                         *
//*   Free Software Foundation, Inc.,                                       *
//*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
//***************************************************************************
//

#ifndef hammingdistance_h
#define hammingdistance_h

#if USE_VECTORCLASS_LIBRARY
#define  HAMMING_VECTOR (1)
#include <vectorclass/vectorclass.h> //For Vec32c and Vec32cb classes
#else
#define  HAMMING_VECTOR (0)
#endif

/**
 * @brief  Determine the hamming distance between two sequences.
 * @tparam L the type representing a state of a pattern
 *           (when the state range is 255 or less, it can be char, 
 *            or if 65535 or less, short.)
 * @tparam F is the type of the elements in the vector of pattern frequencies
 *           (if all pattern frequencies were 255 or less, it could be 
 *            unsigned char, or if 65535 or less, short).
 * @param  unknown   - a wildcard character that will match any character
 * @param  sequenceA - pointer to the pattern states in the first sequence
 * @param  sequenceB - pointer to the pattern states in the second sequence
 * @param  seqLen    - length of the sequence
 * @param  frequencyVector     - pointer to the frequencies of the patterns
 * @param  frequencyOfUnknowns - will be set to the sum of the frequencies,
 *                               for the patterns for which the pattern
 *                               state was unknown, in either of both of the
 *                               sequences.
 * @return double - the sum of the frequencies, for the patterns that
 *         differ, between the sequences, where neither of the patterns
 *         has an unknown state.
 */
template <class L, class F> inline double hammingDistance
        ( L unknown, const L* sequenceA, const L* sequenceB
         , int seqLen, const F* frequencyVector, double& frequencyOfUnknowns ) {
    double distance = 0;
    frequencyOfUnknowns = 0;
    for (int pos = 0; pos < seqLen; ++pos) {
        if ( sequenceA[pos] == unknown ||
            sequenceB[pos] == unknown ) {
            frequencyOfUnknowns += frequencyVector[pos];
        } else if (sequenceA[pos] != sequenceB[pos]) {
            distance += frequencyVector[pos];
        }
    }
    return distance;
}

/**
 * @brief  Find the sum of the pattern frequencies, corresponding to the
 *         the pattern states, in a sequence, that are >= some boundary
 *         state.
 * @tparam L - the type of the pattern states
 * @tparam F - the type of the pattern frequencies
 * @param  boundaryChar    - the boundary state
 * @param  sequence        - the sequence (not of site, but of pattern states)
 * @param  seqLen          - the length of the sequence (number of patterns)
 * @param  frequencyVector - the frequencies, corresponding to the patterns
 *                           whose states appear in the sequence.
 * @return size_t 
 * @note   The proper name (for what this function does) is perhaps: 
 *         sumOfFrequenciesOfCharactersGreaterThanOrEqualTo().
 *         But for now, this function is named after what
 *         it is used for.
 */
template <class L, class F> inline size_t sumForUnknownCharacters
    (L boundaryChar, const L* sequence, int seqLen, const F* frequencyVector) {
    size_t sum = 0;
    for (int pos = 0; pos < seqLen; ++pos) {
        if (boundaryChar <= sequence[pos]) {
            sum += frequencyVector[pos];
        }
    }
    return sum;
}

#if (HAMMING_VECTOR)
inline double hammingDistance
( char unknown, const char* sequenceA, const char* sequenceB
 , intptr_t seqLen, const int* frequencyVector
 , double& frequencyOfUnknowns ) {
    size_t blockStop = 0;
    int    distance  = 0;
    int    freqUnknown = 0;
    Vec32c blockA;
    intptr_t blockSize = blockA.size(); //but see scratch. Needs to match
    if (blockSize < seqLen) {
        blockStop = seqLen - (seqLen & (blockSize-1));
        auto aStop = sequenceA + blockStop;
        Vec32c  blockB;
        const int* f = frequencyVector;
        auto a = sequenceA;
        auto b = sequenceB;
        while (a<aStop) {
            blockA.load(a);
            blockB.load(b);
            Vec32cb diff32    = blockA != blockB;
            Vec32cb known32   = (blockA != unknown) & (blockB != unknown);
            Vec32cb unknown32 = ~known32;
            diff32 &= known32;
            intptr_t i = horizontal_find_first(diff32);
            if (0<=i) {
                char scratchDiff[32];
                diff32.store(scratchDiff);
                for (; i<blockSize; ++i) {
                    if (scratchDiff[i]) {
                        distance += f[i];
                    }
                }
            }
            i = horizontal_find_first(unknown32);
            if (0<=i) {
                char scratchUnknown[32];
                unknown32.store(scratchUnknown);
                for (; i<blockSize; ++i) {
                    if (scratchUnknown[i]) {
                        freqUnknown += f[i];
                    }
                }
            }
            f += blockSize;
            a += blockSize;
            b += blockSize;
        }
    }
    for (intptr_t pos=blockStop; pos < seqLen; ++pos ) {
        if (sequenceA[pos]==unknown || sequenceB[pos]==unknown) {
            freqUnknown += frequencyVector[pos];
        }
        else if (sequenceA[pos] != sequenceB[pos]) {
            distance += frequencyVector[pos];
        }
    }
    frequencyOfUnknowns = freqUnknown;
    return distance;
}
#endif

/**
 * @brief  Determine the Hamming distance, the count of sites with disagreeing 
 *         characters, between two sequences of sites by comparing them, 
 *         site by site, treating differences as 1, except when one - or the 
 *         other - of the sequences has an unknown state at the site, 
 *         in which case the sites are counted as being in agreement.
 * @param  unknown   - the character representing a site of "unknown" character
 * @param  sequenceA - the first sequence
 * @param  sequenceB - the second sequence
 * @param  seqLen    - the length of *both* sequences
 * @return size_t    - the number of differences
 */
inline size_t conventionalHammingDistance(char unknown,
                                          const char* sequenceA,
                                          const char* sequenceB,
                                          size_t seqLen ) {
    size_t distance = 0;
    const char* stopSequenceA = sequenceA+seqLen;
    for (; sequenceA < stopSequenceA; ++sequenceA, ++sequenceB) {
        if (*sequenceA != *sequenceB) {
            if (*sequenceA!=unknown && *sequenceB!=unknown) {
                ++distance;
            }
        }
    }
    return distance;
}

#if (defined(CLANG_UNDER_VS) || defined(_MSC_VER)) && !defined(_mm_popcnt_u64)
#   define _mm_popcnt_u64 __popcnt64
#endif

/**
 * @brief Count the number of bits that are set in the bitwise or of two
 *        blocks (of uint64_t bitfields), of the same size.
 * @param a - the first block
 * @param b - the second block
 * @param count - the number of uint64_t bitfields (in both blocks)
 *        the number of bits in the blocks will be 64x this.
 */
inline uint64_t conventionalCountBitsSetInEither
    (const uint64_t* a, const uint64_t* b,
    size_t count /*in uint64_t*/) {
    uint64_t bits_set = 1;
    for (size_t i=0; i<count; ++i) {
        uint64_t bitsInEither = a[i] | b[i];
        uint64_t delta;
        #if (defined (__GNUC__) || defined(__clang__)) && !defined(CLANG_UNDER_VS)
            __asm("popcntq %1, %0" : "=r"(delta) : "r"(bitsInEither) : );
        #else
            delta = _mm_popcnt_u64(bitsInEither);
        #endif
        bits_set += delta;
    }
    return bits_set;
}

#if (!HAMMING_VECTOR)
/**
 * @brief version of vectorHammingDistance to use if hamming distance
 *        functions are NOT to be vectorized.
 */
inline uint64_t vectorHammingDistance(char unknown, const char* sequenceA,
                                      const char* sequenceB, size_t seqLen ) {
    return conventionalHammingDistance(unknown, sequenceA, sequenceB, seqLen);
}
/**
 * @brief version of countBitsSetInEither to use if hamming distance
 *        functions are NOT to be vectorized.
 */
inline uint64_t countBitsSetInEither
    (copst uint64_t* a, const uint64_t* b, size_t count) {
    return conventionalCountBitsSetInEither(a, b, count);
}
#else

#ifdef _MSC_VER
    #define ALIGN_32(x) __declspec(align(32)) x = {0}
#else
    #define ALIGN_32(x) x __attribute__((aligned(32)))
#endif

//W is the number of bytes in a CV (or CBV)
/**
 * @brief  Determine hamming distance between two sequences 
 *         (using the specified vector types, supplied via 
 *         the template parameters).
 * 
 * @tparam DV  - distance vector type  (of W/8 uint64_t)
 * @tparam CV  - character vector type (of W characters)
 * @tparam CBV - boolean vector type   (same rank as CV)
 * @tparam W   - the number of characters per instance of CV
 * @param  unknown   - the character that indicates a site is of unknown character
 * @param  sequenceA - pointer to the start of the first sequence
 *                     (one character per site).
 * @param  sequenceB - pointer to the start of the second sequence
 * @param  seqLen    - the number of sites/characters in each of the sequences.
 * @return uint64_t  - The hamming distance
 * @note   it is assumed that the binary representation, in a CBV, for true,
 *         is 0xFF. So 8 bits will be set, in each uint_64 in the CBV, for each
 *         of the W/8 uint64_t's in the CBV.
 * @note   runs sequentially.
 *         (even when the _OPENMP symbol is defined and non-zero, it is supposed 
 *          that each thread is calculating hamming distances, for different 
 *          sequences, so comparisons of individual sequences might as well
 *          be single-threaded).
 */
template <class DV, class CV, class CBV, int W /*power of 2*/>
inline uint64_t vectorHammingDistanceTemplate(char unknown,
                              const char* sequenceA,
                              const char* sequenceB,
                              size_t seqLen ) {
    size_t   blockStop        = 0;
    DV       distance_vector  = 0;
    if (W <= seqLen) {
        blockStop  = seqLen - (seqLen & (W-1));
        CV blockA(0); //next 32 bytes from a
        CV blockB(0); //next 32 bytes from b
        CV blockU  = unknown;
        ALIGN_32(uint64_t vec[W/8]);
        ALIGN_32(uint64_t res[W/8]);
        DV distance_bump_vector(0);
        for ( size_t i=0; i<blockStop; i+=W) {
            //Compare W characters at a time
            blockA.load(sequenceA+i);
            blockB.load(sequenceB+i);
            CBV diff32 = (blockA != blockB);
            if (horizontal_find_first(diff32)<0) {
               continue;
            }
            diff32 &= (blockA != blockU)
                    & (blockB != blockU);
            //Count the number of bits set in each 8-byte stretch
            //of vec (8 bits will be set for each character
            //that differs), in a vector of [W/8] x uint64_t.
            CV(diff32).store( reinterpret_cast<char*>(&vec[0]));
            for (int j=0; j<W/8; ++j) {
                #if (defined (__GNUC__) || defined(__clang__)) && !defined(CLANG_UNDER_VS)
                    #ifndef WIN32
                        __asm("popcntq %1, %0" : "=r"(res[j]) : "r"(vec[j]) : );
                    #else
                        uint32_t* vec32 = reinterpret_cast<uint32_t*>(&vec[j]);
                        uint32_t  first_half;
                        uint32_t  second_half;                
                        __asm("popcnt %1, %0" : "=r"(first_half) : "r"(vec32[0]) : );
                        __asm("popcnt %1, %0" : "=r"(second_half) : "r"(vec32[1]) : );
                        res[j] = first_half + second_half;                    
                        #endif
                #else
                    res[j] = _mm_popcnt_u64(vec[j]);
                #endif
            }
            //Vector addition (each element in the distance_vector
            //is 8 times the count of known characters that differed
            //in the corresponding 8-byte slices of the W-byte blocks.
            distance_bump_vector.load(&res[0]);
            distance_vector += distance_bump_vector;
        }
    }
    //remember distance_vector counts each difference as 8
    //so >>3 to divide by 8.
    uint64_t distance    = horizontal_add(distance_vector) >> 3;
    for (size_t pos=blockStop; pos < seqLen; ++pos ) {
        if (sequenceA[pos]!=sequenceB[pos]) {
            if (sequenceA[pos]!=unknown && sequenceB[pos]!=unknown) {
                ++distance;
            }
        }
    }
    return distance;
}

/**
 * @brief  Count the number of bits set in the bitwise-or of two blocks
 * @tparam V a vector type (of uint64_t)
 * @tparam W the number of uint64_t's in a V.
 * @param  a      - the first block
 * @param  b     - the second block
 * @param  count - the number of uint64_t bitfields (in both blocks)
 *                 the number of bits in the blocks will be 64x this.
 * @return uint64_t - the number of bits set in the bitwise-or of the blocks.
 */
template <class V, int W>
uint64_t countBitsSetInEitherTemplate
        (const uint64_t* a, const uint64_t* b,
        intptr_t count /* in instances of uint64_t */) {
    V count_vector  = 0;
    V aData = 0;
    V bData = 0;
    ALIGN_32(uint64_t vec[W]);
    ALIGN_32(uint64_t res[W]);
    V dData = 0;
    for (intptr_t i=0; i<count; i+=W) {
        aData.load(a+i);
        bData.load(b+i);
        (aData | bData).store(&vec[0]);
        for (int j=0; j<W; ++j) {
            #if (defined (__GNUC__) || defined(__clang__)) && !defined(CLANG_UNDER_VS)
                #ifndef WIN32
                    __asm("popcntq %1, %0" : "=r"(res[j]) : "r"(vec[j]) : );
                #else
                    uint32_t* vec32 = reinterpret_cast<uint32_t*>(&vec[j]);
                    uint32_t  first_half;
                    uint32_t  second_half;                
                    __asm("popcnt %1, %0" : "=r"(first_half) : "r"(vec32[0]) : );
                    __asm("popcnt %1, %0" : "=r"(second_half) : "r"(vec32[1]) : );
                    res[j] = first_half + second_half;                    
                #endif
            #else
                res[j] = _mm_popcnt_u64(vec[j]);
            #endif
        }
        dData.load(&res[0]);
        count_vector += dData;
    }
    return horizontal_add(count_vector);
}

#if (INSTRSET >= 7)
    inline uint64_t vectorHammingDistance(char unknown, const char* sequenceA,
                                          const char* sequenceB, size_t seqLen ) {
        return vectorHammingDistanceTemplate<Vec4uq, Vec32c, Vec32cb, 32>
                (unknown, sequenceA, sequenceB, seqLen);
    }
    inline uint64_t countBitsSetInEither
        (const uint64_t* a, const uint64_t* b,
        size_t count /*in uint64_t*/) {
        return countBitsSetInEitherTemplate<Vec4uq, 4>(a, b, count);
    }
#elif (INSTRSET >= 2)
    inline uint64_t vectorHammingDistance(char unknown, const char* sequenceA,
                                          const char* sequenceB, size_t seqLen ) {
        return vectorHammingDistanceTemplate<Vec2uq, Vec16c, Vec16cb, 16>
                (unknown, sequenceA, sequenceB, seqLen);
    }
    inline uint64_t countBitsSetInEither
            (const uint64_t* a, const uint64_t* b,
            size_t count /*in uint64_t*/) {
        return countBitsSetInEitherTemplate<Vec2uq, 2>(a, b, count);
    }
#else
    inline uint64_t vectorHammingDistance(char unknown, const char* sequenceA,
                                          const char* sequenceB, size_t seqLen ) {
        return conventionalHammingDistance(unknown, sequenceA, sequenceB, seqLen);
    }
    inline uint64_t countBitsSetInEither(uint64_t* a, uint64_t* b,
                                       size_t count /*in uint64_t*/) {
        return conventionalCountBitsSetInEither(a,b,count);
    }
#endif

/**
 * @brief  Given a sequence, that represents distinct patterns (for each
 *         pattern there is a character, the one found in this sequence, 
 *         and a frequency: the count of sites, in the alignment that
 *         is being processed - which we can't see! - for which the
 *         pattern is found, or repeats).  
 * @param  boundaryChar - characters equal to, or greater than this,
 *                        are counted as unknown.
 * @param  sequence     - pointer to the characters for the patterns
 * @param  seqLen       - the number of distinct patterns (of characters
 *                        in a site, over all of the sequences)
 * @param  frequencyVector - for each pattern, the number of sites that
 *                           match that pattern.
 * @return size_t - the sum, of the frequencyVector[i]'s for each 
 *         pattern, for which sequence[i], the character at site i,
 *         is equal to or greater than boundaryChar.
 * @note   When the number of sequences in an alignment is relatively low, 
 *         there can be only a few distinct patterns of sites, each of which
 *         (we may hope) accounts for a lot of sites.  When the number of
 *         sequences in the alignment is large, however, there might be 
 *         relatively many distinct patterns.
 */
inline size_t sumForUnknownCharacters
    ( char boundaryChar, const char* sequence, intptr_t seqLen, const int* frequencyVector) {
    size_t sum = 0;
    Vec32c blockBoundary = boundaryChar;
    int blockSize = blockBoundary.size(); //but see scratch. Needs to match
    int pos = 0;
    if (blockSize < seqLen) {
        //Next line assumes that: blockSize is power of 2
        int blockSize_less_1 = blockSize - 1;
        intptr_t blockStop = seqLen - (seqLen & (blockSize_less_1)); 
        for (; pos < blockStop; pos += blockSize) {
            Vec32c blockRead;
            blockRead.load(sequence + pos);
            Vec32cb blockUnknown = ( blockBoundary <= blockRead ) ;
            int i = horizontal_find_first(blockUnknown);
            if (0 <= i) {
                int nextBlock = pos + blockSize;
                for (int scan = pos + i; scan < nextBlock; ++scan) {
                    if (boundaryChar <= sequence[scan]) {
                        sum += frequencyVector[scan];
                    }
                }
            }
        }
    }
    for (; pos < seqLen; ++pos) {
        if (boundaryChar <= sequence[pos]) {
            sum += frequencyVector[pos];
        }
    }
    return sum;
}

#endif

/**
 * @brief  Count the number of bits set in a block of (count) 64-bit
 *         bit-fields, pointed to by a unit64_t* (a pointer to the first
 *         unsigned 64-bit integer)
 * @param  a     - the block
 * @param  count - the number of uint64_t bit-fields in the block
 * @return uint64_t - the number of bits set
 * @note   This is rather inefficient! It returns the number of bits
 *         set in the bit-wise or of the block, starting at a, with itself.
 *         It works this way because I was in a bit of a hurry when I
 *         wrote it.
 */
inline uint64_t countBitsSetIn(const uint64_t* a, size_t count) {
    return countBitsSetInEither(a,a,count);
}

#endif /* hammingdistance_h */
