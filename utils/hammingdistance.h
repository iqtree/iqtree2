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

#define HAMMING_VECTOR (1)
#include <vectorclass/vectorclass.h> //For Vec32c and Vec32cb classes

//
//Note 1: L is a template parameter so that, when the state range
//        is 255 or less, it can be char, or if 65535 or less, short.
//Note 2: F is a template parameter so that, when the maximum
//        frequency is 255 or less, it could be unsigned char,
//        or if 65535 or less, short.
//
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

//Proper name (for what this function does) is perhaps: 
//sumOfFrequenciesOfCharactersGreaterThanOrEqualTo
//But for now, this function is named after what
//it is used for.
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
 , int seqLen, const int* frequencyVector
 , double& frequencyOfUnknowns ) {
    size_t blockStop = 0;
    int    distance  = 0;
    int    freqUnknown = 0;
    Vec32c blockA;
    size_t blockSize = blockA.size(); //but see scratch. Needs to match
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
            int i = horizontal_find_first(diff32);
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
    for (int pos=blockStop; pos < seqLen; ++pos ) {
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

inline size_t conventionalCountBitsSetInEither(uint64_t* a, uint64_t* b,
                                               size_t count /*in uint64_t*/) {
    size_t bits_set = 1;
    for (int i=0; i<count; ++i) {
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

#if (INSTRSET >= 2)
#if defined(CLANG_UNDER_VS) && !defined(_mm_popcnt_u64)
#   define _mm_popcnt_u64 __popcnt64
#endif

#ifdef _MSC_VER
    #define ALIGN_32(x) __declspec(align(32)) x
#else
    #define ALIGN_32(x) x __attribute__((aligned(32)))
#endif

template <class DV, class CV, class CBV, int W /*power of 2*/>
inline uint64_t vectorHammingDistanceTemplate(char unknown,
                              const char* sequenceA,
                              const char* sequenceB,
                              size_t seqLen ) {
    //Note: this is designed to run sequential (it is supposed that
    //each thread is calculating hamming distances for different
    //sequences).
    size_t   blockStop        = 0;
    DV       distance_vector  = 0;
    if (W <= seqLen) {
        blockStop  = seqLen - (seqLen & (W-1));
        CV blockA; //next 32 bytes from a
        CV blockB; //next 32 bytes from b
        CV blockU  = unknown;
        ALIGN_32(uint64_t vec[W/8]);
        ALIGN_32(uint64_t res[W/8]);
        DV distance_bump_vector;
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
            //Count the number of bits set in each 8-byte quarter
            //of diff32 (8 bits will be set for each character
            //that differs), in a vector of 4 x uint64_t.
            CV(diff32).store( reinterpret_cast<char*>(&vec[0]));
            for (int j=0; j<W/8; ++j) {
                #if (defined (__GNUC__) || defined(__clang__)) && !defined(CLANG_UNDER_VS)
                    __asm("popcntq %1, %0" : "=r"(res[j]) : "r"(vec[j]) : );
                #else
                    res[j] = _mm_popcnt_u64(vec[j]);
                #endif
            }
            //Vector addition (each element in the distance_vector
            //is 8 times the count of known characters that differed
            //in the corresponding 8-byte quarters of the W-byte blocks.
            distance_bump_vector.load(&res[0]);
            distance_vector += distance_bump_vector;
        }
    }
    //remember distance_vector counts each difference as 8
    uint64_t distance    = horizontal_add(distance_vector) >> 3;
    for (int pos=blockStop; pos < seqLen; ++pos ) {
        if (sequenceA[pos]!=sequenceB[pos]) {
            if (sequenceA[pos]!=unknown && sequenceB[pos]!=unknown) {
                ++distance;
            }
        }
    }
    return distance;
}

template <class V, int W>
size_t countBitsSetInEitherTemplate(uint64_t* a, uint64_t* b,
                                    size_t count /*in uint64_t*/) {
    //Assumes: W divides count
    V count_vector  = 0;
    V aData;
    V bData;
    ALIGN_32(uint64_t vec[W]);
    ALIGN_32(uint64_t res[W]);
    V dData;
    for (int i=0; i<count; i+=W) {
        aData.load(a+i);
        bData.load(b+i);
        (aData | bData).store(&vec[0]);
        for (int j=0; j<W/8; ++j) {
            #if (defined (__GNUC__) || defined(__clang__)) && !defined(CLANG_UNDER_VS)
                __asm("popcntq %1, %0" : "=r"(res[j]) : "r"(vec[j]) : );
            #else
                res[j] = _mm_popcnt_u64(vec[j]);
            #endif
        }
        dData.load(&res[0]);
        count_vector += dData;
    }
    return horizontal_add(count_vector);
}
#endif

#if (INSTRSET >= 7)
    inline uint64_t vectorHammingDistance(char unknown, const char* sequenceA,
                                          const char* sequenceB, size_t seqLen ) {
        return vectorHammingDistanceTemplate<Vec4uq, Vec32c, Vec32cb, 32>
                (unknown, sequenceA, sequenceB, seqLen);
    }
    inline size_t countBitsSetInEither(uint64_t* a, uint64_t* b,
                                    size_t count /*in uint64_t*/) {
        return countBitsSetInEitherTemplate<Vec4uq, 4>(a, b, count);
    }
#elif (INSTRSET >= 2)
    inline uint64_t vectorHammingDistance(char unknown, const char* sequenceA,
                                          const char* sequenceB, size_t seqLen ) {
        return vectorHammingDistanceTemplate<Vec2uq, Vec16c, Vec16cb, 16>
                (unknown, sequenceA, sequenceB, seqLen);
    }
    inline size_t countBitsSetInEither(uint64_t* a, uint64_t* b,
                                size_t count /*in uint64_t*/) {
        return countBitsSetInEitherTemplate<Vec2uq, 2>(a, b, count);
    }
#else
    inline uint64_t vectorHammingDistance(char unknown, const char* sequenceA,
                                          const char* sequenceB, size_t seqLen ) {
        return conventionalHammingDistance(unknown, sequenceA, sequenceB, seqLen);
    }
    inline size_t countBitsSetInEither(uint64_t* a, uint64_t* b,
                                       size_t count /*in uint64_t*/) {
        return conventionalCountBitsSetInEither(a,b,count);
    }
#endif

inline size_t sumForUnknownCharacters
    ( char boundaryChar, const char* sequence, int seqLen, const int* frequencyVector) {
    size_t sum = 0;
    Vec32c blockBoundary = boundaryChar;
    int blockSize = blockBoundary.size(); //but see scratch. Needs to match
    int pos = 0;
    if (blockSize < seqLen) {
        //Next line assumes that: blockSize is power of 2
        int blockStop = seqLen - (seqLen & (blockSize - 1)); 
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


#endif /* hammingdistance_h */
