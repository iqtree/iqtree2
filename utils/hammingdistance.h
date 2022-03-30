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
#define VECTOR_MAD     (0)
#include <vectorclass/vectorclass.h> //For Vec32c and Vec32cb classes

//
//Note 1: L is a template parameter so that, when the state range
//        is 255 or less, it can be char, or if 65535 or less, short.
//Note 2: F is a template parameter so that, when the maximum
//        frequency is 255 or less, it could be unsigned char,
//        or if 65535 or less, short.
//
#if (!HAMMING_VECTOR)
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
#endif

#if (HAMMING_VECTOR)
inline double hammingDistance
( char unknown, const char* sequenceA, const char* sequenceB
 , int seqLen, const int* frequencyVector
 , double& frequencyOfUnknowns ) {
    int blockStop = 0;
    int distance  = 0;
    int freqUnknown = 0;
#if !defined(__ARM_NEON)
    Vec32c blockA;
    int blockSize = blockA.size(); //but see scratch. Needs to match
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
#else
    // NHANLT: FUTURE WORK: implementing method to compute the Hamming distance for Apple Silicon Chipset
    throw "Sorry! Computing of Hamming distance has not yet implemented for Apple Silicon Chipset!";
#endif
    return distance;
}
#endif


#endif /* hammingdistance_h */
