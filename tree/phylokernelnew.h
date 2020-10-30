/*
 * phylokernelnew.h
 * Newly revised kernel based on vectorizing over alignment patterns
 *
 *  Created on: Sept 23, 2016
 *      Author: minh
 */


#if !defined(PHYLOKERNELNEW_H_) || !defined(PHYLOKERNELNEW_STATE_H_)

#define JAMES_VERSION 1

#ifdef KERNEL_FIX_STATES
#   define PHYLOKERNELNEW_STATE_H_
#else
#   define PHYLOKERNELNEW_H_
#endif

#include "phylotree.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/*******************************************************
 *
 * Helper function for vectors and matrix multiplication
 *
 ******************************************************/


#ifdef KERNEL_FIX_STATES

template <class VectorClass> std::string vc_vector_to_string(VectorClass* p, size_t N) {
    std::stringstream s;
    s << "{ ";
    double tmp[VectorClass::size()];
    for (size_t r = 0; r<N; ++r) {
        if (0<r) s << ", ";
        p[r].store(&tmp[0]);
        s << "{";
        for (size_t c = 0 ; c < VectorClass::size(); ++c ) {
            if (0<c) s << ", ";
            s << tmp[c];
        }
        s << "}";
    }
    s << " }";
    return s.str();
}
#endif

/**
    sum of elments of a vector:
    X = A[0] + ... + A[N-1]
    template FMA = true to allow FMA instruction, false otherwise
    @param N number of elements
    @param A vector of size N
    @param[out] X sum of elements of A
*/
#ifndef KERNEL_FIX_STATES
template <class VectorClass, const bool append>
inline void sumVec(const VectorClass *A, VectorClass &X, size_t N)
{
    if (!append) {
        X=0;
    }
    size_t cruft = (N & 3);
    if (cruft<N) {
        VectorClass  V[4]  = { 0, 0, 0, 0 };
        for ( const VectorClass* Astop = A + ( N - cruft ); A<Astop; A+=4 ) {
            V[0] += A[0];
            V[1] += A[1];
            V[2] += A[2];
            V[3] += A[3];
        }
        V[0] += V[2];
        V[1] += V[3];
        V[0] += V[1];
        X    += V[0];
        if (cruft==0) {
            return;
        }
    }
    for ( const VectorClass* Astop = A + cruft; A<Astop; ++A ) {
        X += A[0];
    }
}
#endif

/**
    dotProduct of two vectors A, B
    X = A.B = A[0]*B[0] + ... + A[N-1]*B[N-1]
    template FMA = true to allow FMA instruction, false otherwise
    @param N number of elements
    @param A first vector of size N
    @param B second vector of size N
    @param[out] X dot-product of A and B
*/
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void dotProductVec(const Numeric *A, const VectorClass *B, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductVec(const Numeric *A, const VectorClass *B, VectorClass &X, size_t N)
#endif
{
#if JAMES_VERSION
    if (N==4) {
        X = (A[0]*B[0] + A[1]*B[1]) + (A[2]*B[2] + A[3]*B[3]);
        return;
    }
    size_t cruft = (N & 3);
    if (N==cruft) {
        X = 0;
        //IF cruft==N, N must be 0, 1, 2, or 3
        for (int i=0; i<cruft; ++i) {
            X = mul_add(A[i], B[i], X);
        }
        return;
    }
    X = A[0]*B[0];
    VectorClass v1 = A[1]*B[1];
    VectorClass v2 = A[2]*B[2];
    VectorClass v3 = A[3]*B[3];
    auto Astop = A + (N-cruft);
    A+=4;
    B+=4;
    for (; A<Astop; A+=4, B+=4) {
        X  = mul_add( A[0], B[0], X);
        v1 = mul_add( A[1], B[1], v1);
        v2 = mul_add( A[2], B[2], v2);
        v3 = mul_add( A[3], B[3], v3);
    }
    if (0<cruft) {
        X  = mul_add( A[0], B[0], X);
        if (1<cruft) {
            v1 = mul_add( A[1], B[1], v1);
            if (2<cruft) {
                v2 = mul_add( A[2], B[2], v2);
            }
        }
    }
    X  = ( X + v1 ) + ( v2 + v3 );
    return;
#else
    // quick treatment for small N <= 4
    switch (N) {
    case 1:
        X = A[0]*B[0];
        return;
    case 2:
        X = mul_add(A[1], B[1], A[0]*B[0]);
        return;
    case 3:
        X = mul_add(A[2], B[2], mul_add(A[1], B[1], A[0]*B[0]));
        return;
    case 4:
        X = mul_add(A[1], B[1], A[0]*B[0]) + mul_add(A[3], B[3], A[2]*B[2]);
        return;
    default: break;
    }

    // For N > 4, add the rest with 4-way unrolling
    size_t i, j;
    switch (N % 4) {
    case 0: {
        VectorClass V[4];
        for (j = 0; j < 4; j++)
            V[j] = A[j] * B[j];
        for (i = 4; i < N; i+=4) {
            for (j = 0; j < 4; j++)
                V[j] = mul_add(A[i+j], B[i+j], V[j]);
        }
        X = (V[0]+V[1]) + (V[2]+V[3]);
        break;
    }

    case 1: {
        const size_t Nm1 = N-1;
        VectorClass V[4];
        for (j = 0; j < 4; j++)
            V[j] = A[j] * B[j];
        for (i = 4; i < Nm1; i+=4) {
            for (j = 0; j < 4; j++)
                V[j] = mul_add(A[i+j], B[i+j], V[j]);
        }
        X = mul_add(A[Nm1], B[Nm1], (V[0]+V[1]) + (V[2]+V[3]));
        break;
    }

    case 2: {
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j];
        for (i = 2; i < N; i+=2) {
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j], B[i+j], V[j]);
        }
        X = (V[0]+V[1]);
        break;
    }

    default: {
        // odd number of states
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j];
        for (i = 2; i < N-1; i+=2) {
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j], B[i+j], V[j]);
        }
        X = mul_add(A[N-1], B[N-1], V[0]+V[1]);
        break;
    }
    }
#endif
}

/**
    Dual dotProduct of four vectors A, B, C, D to compute X:
    X = (A.B) * (C.D), where
    A.B = A[0]*B[0] + ... + A[N-1]*B[N-1]
    C.D = C[0]*D[0] + ... + C[N-1]*D[N-1]
    template FMA = true to allow FMA instruction, false otherwise
    @param N number of elements
    @param A first vector of size N
    @param B second vector of size N
    @param C third vector of size N
    @param D fourth vector of size N
    @param[out] X = (A.B) * (C.D)
*/
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void dotProductDualVec(const Numeric *A, const VectorClass *B, const Numeric *C, const VectorClass *D, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductDualVec(const Numeric *A,const  VectorClass *B, const Numeric *C, const VectorClass *D, VectorClass &X, size_t N)
#endif
{
    // quick treatment for small N <= 4
    switch (N) {
    case 1:
        X = (A[0]*B[0])*(C[0]*D[0]);
        return;
    case 2:
        X = mul_add(A[1],B[1],A[0]*B[0]) * mul_add(C[1],D[1],C[0]*D[0]);
        return;
    case 3:
        X = mul_add(A[2], B[2], mul_add(A[1], B[1], A[0]*B[0])) *
            mul_add(C[2], D[2], mul_add(C[1], D[1], C[0]*D[0]));
        return;
    case 4:
        X = (mul_add(A[1], B[1], A[0]*B[0]) + mul_add(A[3], B[3], A[2]*B[2])) *
            (mul_add(C[1], D[1], C[0]*D[0]) + mul_add(C[3], D[3], C[2]*D[2]));
        return;
    default: break;
    }

    // For N > 4, add the rest with 4-way unrolling
    size_t i, j;
    switch (N % 4) {
    case 0: {
        VectorClass AB[4], CD[4];
        for (j = 0; j < 4; j++) {
            AB[j] = A[j] * B[j];
            CD[j] = C[j] * D[j];
        }
        for (i = 4; i < N; i+=4) {

            for (j = 0; j < 4; j++) {
                AB[j] = mul_add(A[i+j],  B[i+j],  AB[j]);
                CD[j] = mul_add(C[i+j], D[i+j], CD[j]);
            }
        }
        X = ((AB[0]+AB[1])+(AB[2]+AB[3])) * ((CD[0]+CD[1])+CD[2]+CD[3]);
        break;
    }

    case 1: {
        const size_t Nm1 = N-1;
        VectorClass AB[4], CD[4];
        for (j = 0; j < 4; j++) {
            AB[j] = A[j] * B[j];
            CD[j] = C[j] * D[j];
        }
        for (i = 4; i < Nm1; i+=4) {

            for (j = 0; j < 4; j++) {
                AB[j] = mul_add(A[i+j],  B[i+j],  AB[j]);
                CD[j] = mul_add(C[i+j], D[i+j], CD[j]);
            }
        }
        X = mul_add(A[Nm1], B[Nm1], (AB[0]+AB[1])+(AB[2]+AB[3])) * mul_add(C[Nm1], D[Nm1], (CD[0]+CD[1])+CD[2]+CD[3]);
        break;
    }

    case 2: {
        VectorClass AB[2], CD[2];
        for (j = 0; j < 2; j++) {
            AB[j] = A[j] * B[j];
            CD[j] = C[j] * D[j];
        }
        for (i = 2; i < N; i+=2) {
            for (j = 0; j < 2; j++) {
                AB[j] = mul_add(A[i+j],  B[i+j],  AB[j]);
                CD[j] = mul_add(C[i+j], D[i+j], CD[j]);
            }
        }
        X = ((AB[0]+AB[1])) * ((CD[0]+CD[1]));
        break;
    }

    default: {
        // odd states
        VectorClass AB[2], CD[2];
        for (j = 0; j < 2; j++) {
            AB[j] = A[j] * B[j];
            CD[j] = C[j] * D[j];
        }
        for (i = 2; i < N-1; i+=2) {
            for (j = 0; j < 2; j++) {
                AB[j] = mul_add(A[i+j],  B[i+j],  AB[j]);
                CD[j] = mul_add(C[i+j], D[i+j], CD[j]);
            }
        }
        AB[0] = mul_add(A[N-1], B[N-1], AB[0]+AB[1]);
        CD[0] = mul_add(C[N-1], D[N-1], CD[0]+CD[1]);
        X = AB[0] * CD[0];
        break;
    }
    }
}

/**
    compute product of a vector A and a matrix M, resulting in a vector X:
    X[i] = A[0]*M[i,0] + ... + A[N-1]*M[i,N-1], for all i = 0,...,N-1
    @param N number of elements
    @param A input vector of size N
    @param M input matrix of size N*N
    @param[out] X output vector of size N
*/
// quick unrolling version of multiplying partial_lh with inv_eigenvector
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void productVecMat(const VectorClass *A, const Numeric *M, VectorClass *X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void productVecMat(const VectorClass *A, const Numeric *M, VectorClass *X, size_t N)
#endif
{
    if (N==4) {
        // manual unrolling
        X[0] = mul_add(A[1],M[1],  A[0]*M[0])  + mul_add(A[3],M[3],  A[2]*M[2]);
        X[1] = mul_add(A[1],M[5],  A[0]*M[4])  + mul_add(A[3],M[7],  A[2]*M[6]);
        X[2] = mul_add(A[1],M[9],  A[0]*M[8])  + mul_add(A[3],M[11], A[2]*M[10]);
        X[3] = mul_add(A[1],M[13], A[0]*M[12]) + mul_add(A[3],M[15], A[2]*M[14]);
        return;
    }
    // quick treatment for small N <= 4
    size_t i;
    switch (N) {
    case 1:
        X[0] = A[0]*M[0];
        return;
    case 2:
        X[0] = mul_add(A[1],M[1],A[0]*M[0]);
        X[1] = mul_add(A[1],M[3],A[0]*M[2]);
        return;
    case 3:
        for (i = 0; i < 3; i++) {
            // manual unrolling
            X[i] = mul_add(A[2],M[2],mul_add(A[1],M[1],A[0]*M[0]));
            M += 3;
        }
        return;
        for (i = 0; i < 4; i++) {
            // manual unrolling
            X[i] = mul_add(A[1],M[1],A[0]*M[0]) + mul_add(A[3],M[3],A[2]*M[2]);
            M += 4;
        }
        return;
    default: break;
    }

    // For N > 4, add the rest with 4-way unrolling
    size_t j, x;

    switch (N % 4) {
    case 0:
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[4];
            for (j = 0; j < 4; j++)
                V[j] = A[j] * M[j];

            for (x = 4; x < N; x+=4) {
                for (j = 0; j < 4; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = (V[0]+V[1])+(V[2]+V[3]);
            M += N;
        }
        break;

    case 1: {
        const size_t Nm1 = N-1;
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[4];
            for (j = 0; j < 4; j++)
                V[j] = A[j] * M[j];

            for (x = 4; x < Nm1; x+=4) {
                for (j = 0; j < 4; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = mul_add(A[Nm1], M[Nm1], (V[0]+V[1])+(V[2]+V[3]));
            M += N;
        }
        break;
    }
    case 2:
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[2];
            for (j = 0; j < 2; j++)
                V[j] = A[j] * M[j];

            for (x = 2; x < N; x+=2) {
                for (j = 0; j < 2; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = (V[0]+V[1]);
            M += N;
        }
        break;
    default:
        // odd number of states
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[2];
            for (j = 0; j < 2; j++)
                V[j] = A[j] * M[j];

            for (x = 2; x < N-1; x+=2) {
                for (j = 0; j < 2; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = mul_add(A[N-1], M[N-1], V[0]+V[1]);
            M += N;
        }
        break;
    }
}


/**
    compute product of a vector A and a matrix M, resulting in a vector X:
    X[i] = A[0]*M[i,0] + ... + A[N-1]*M[i,N-1], for all i = 0,...,N-1
    and also return the maximum of absolute values of X
    @param N number of elements
    @param A input vector of size N
    @param M input matrix of size N*N
    @param[out] X output vector of size N
    @param[out] Xmax max of |X[i]|
*/
// quick unrolling version of multiplying partial_lh with inv_eigenvector
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void productVecMat(const VectorClass *A, const Numeric *M, VectorClass *X, VectorClass &Xmax)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void productVecMat(const VectorClass *A, const Numeric *M, VectorClass *X, VectorClass &Xmax, size_t N)
#endif
{
    if (N==4) {
        // manual unrolling
        X[0] = mul_add(A[1],M[1],A[0]*M[0])   + mul_add(A[3],M[3],A[2]*M[2]);
        Xmax = max(Xmax, abs(X[0]));
        X[1] = mul_add(A[1],M[5],A[0]*M[4])   + mul_add(A[3],M[7],A[2]*M[6]);
        Xmax = max(Xmax, abs(X[1]));
        X[2] = mul_add(A[1],M[9],A[0]*M[8])   + mul_add(A[3],M[11],A[2]*M[10]);
        Xmax = max(Xmax, abs(X[2]));
        X[3] = mul_add(A[1],M[13],A[0]*M[12]) + mul_add(A[3],M[15],A[2]*M[14]);
        Xmax = max(Xmax, abs(X[3]));
        return;
    }
    // quick treatment for small N <= 4
    size_t i;
    switch (N) {
    case 1:
        X[0] = A[0]*M[0];
        Xmax = max(Xmax, abs(X[0]));
        return;
    case 2:
        X[0] = mul_add(A[1],M[1],A[0]*M[0]);
        X[1] = mul_add(A[1],M[3],A[0]*M[2]);
        Xmax = max(Xmax, max(abs(X[0]), abs(X[1])));
        return;
    case 3:
        for (i = 0; i < 3; i++) {
            // manual unrolling
            X[i] = mul_add(A[2],M[2],mul_add(A[1],M[1],A[0]*M[0]));
            Xmax = max(Xmax, abs(X[i]));
            M += 3;
        }
        return;
    default: break;
    }

    // For N > 4, add the rest with 4-way unrolling
    size_t j, x;

    switch (N % 4) {
    case 0:
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[4];
            for (j = 0; j < 4; j++)
                V[j] = A[j] * M[j];

            for (x = 4; x < N; x+=4) {
                for (j = 0; j < 4; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = (V[0]+V[1])+(V[2]+V[3]);
            M += N;
            Xmax = max(Xmax, abs(X[i]));
        }
        break;

    case 1: {
        const size_t Nm1 = N-1;
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[4];
            for (j = 0; j < 4; j++)
                V[j] = A[j] * M[j];

            for (x = 4; x < Nm1; x+=4) {
                for (j = 0; j < 4; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = mul_add(A[Nm1], M[Nm1], (V[0]+V[1])+(V[2]+V[3]));
            M += N;
            Xmax = max(Xmax, abs(X[i]));
        }
        break;
    }
    case 2:
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[2];
            for (j = 0; j < 2; j++)
                V[j] = A[j] * M[j];

            for (x = 2; x < N; x+=2) {
                for (j = 0; j < 2; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = (V[0]+V[1]);
            M += N;
            Xmax = max(Xmax, abs(X[i]));
        }
        break;

    default:
        // odd number of states
        for (i = 0; i < N; i++) {
            // manual unrolling
            VectorClass V[2];
            for (j = 0; j < 2; j++)
                V[j] = A[j] * M[j];

            for (x = 2; x < N-1; x+=2) {
                for (j = 0; j < 2; j++)
                    V[j] = mul_add(A[x+j], M[x+j], V[j]);
            }
            X[i] = mul_add(A[N-1], M[N-1], V[0]+V[1]);
            M += N;
            Xmax = max(Xmax, abs(X[i]));
        }
        break;
    }
}

/**
    compute dot-products of 2 vectors A, B with a single vector D and returns X, Y:
    X +=   A.D = A[0]*D[0] + ... + A[N-1]*D[N-1]
    Y +=   B.D = B[0]*D[0] + ... + B[N-1]*D[N-1]
    @param N number of elements
    @param nstates number of states
    @param A vector of size N
    @param B vector of size N
    @param D vector of size N
    @param[in/out] X += A.D
    @param[in/out] Y += B.D
*/
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void dotProductPairAdd(Numeric *A, Numeric *B, VectorClass *D,
    VectorClass &X, VectorClass &Y)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductPairAdd(Numeric *A, Numeric *B, VectorClass *D,
    VectorClass &X, VectorClass &Y, size_t N)
#endif
{
    if (N == 1) {
        X = mul_add(A[0], D[0], X);
        Y = mul_add(B[0], D[0], Y);
        return;
    }
#if JAMES_VERSION
    VectorClass AD[2], BD[2];
    for (size_t j = 0; j < 2; j++) {
        AD[j] = A[j] * D[j];
        BD[j] = B[j] * D[j];
    }
    size_t iStop = N - (N&1);
    for (size_t i = 2; i < iStop; i+=2) {
        for (size_t j = 0; j < 2; j++) {
            AD[j] = mul_add(A[i+j], D[i+j], AD[j]);
            BD[j] = mul_add(B[i+j], D[i+j], BD[j]);
        }
    }
    if ((N & 1) == 0) {
        X += AD[0] + AD[1];
        Y += BD[0] + BD[1];
    } else {
        X += mul_add(A[N-1], D[N-1], AD[0] + AD[1]);
        Y += mul_add(B[N-1], D[N-1], BD[0] + BD[1]);
    }
#else
    size_t i, j;
    if (N % 2 == 0) {
        VectorClass AD[2], BD[2];
        for (j = 0; j < 2; j++) {
            AD[j] = A[j] * D[j];
            BD[j] = B[j] * D[j];
        }
        for (i = 2; i < N; i+=2) {
            for (j = 0; j < 2; j++) {
                AD[j] = mul_add(A[i+j], D[i+j], AD[j]);
                BD[j] = mul_add(B[i+j], D[i+j], BD[j]);
            }
        }
        X += AD[0] + AD[1];
        Y += BD[0] + BD[1];
    } else {
        // odd states
        VectorClass AD[2], BD[2];
        for (j = 0; j < 2; j++) {
            AD[j] = A[j] * D[j];
            BD[j] = B[j] * D[j];
        }
        for (i = 2; i < N-1; i+=2) {
            for (j = 0; j < 2; j++) {
                AD[j] = mul_add(A[i+j], D[i+j], AD[j]);
                BD[j] = mul_add(B[i+j], D[i+j], BD[j]);
            }
        }
        X += mul_add(A[N-1], D[N-1], AD[0] + AD[1]);
        Y += mul_add(B[N-1], D[N-1], BD[0] + BD[1]);
    }
#endif
}

/**
    compute dot-products of 3 vectors A, B, C with a single vector D and returns X, Y, Z:
    X =   A.D = A[0]*D[0] + ... + A[N-1]*D[N-1]
    Y =   B.D = B[0]*D[0] + ... + B[N-1]*D[N-1]
    Z =   C.D = C[0]*D[0] + ... + C[N-1]*D[N-1]
    @param N number of elements
    @param nstates number of states
    @param A vector of size N
    @param B vector of size N
    @param C vector of size N
    @param D vector of size N
    @param[in/out] X = A.D
    @param[out] Y = B.D
    @param[out] Z = C.D
*/
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t nstates,
          const bool FMA, const bool ADD>
inline void dotProductTriple(const Numeric *A, const Numeric *B,
                             const Numeric *C, const VectorClass *D,
                             VectorClass &X, VectorClass &Y, VectorClass &Z,
                             size_t N)
#else
template <class VectorClass, class Numeric, const bool FMA, const bool ADD>
inline void dotProductTriple(const Numeric *A, const Numeric *B,
                             const Numeric *C, const VectorClass *D,
                             VectorClass &X, VectorClass &Y, VectorClass &Z,
                             size_t N, size_t nstates)
#endif
{
    size_t i, j;
    if ((nstates & 1) == 0) {
        VectorClass AD[2], BD[2], CD[2];
        for (j = 0; j < 2; j++) {
            AD[j] = A[j] * D[j];
            BD[j] = B[j] * D[j];
            CD[j] = C[j] * D[j];
        }
		for (i = 2; i < N; i+=2) {
            for (j = 0; j < 2; j++) {
                AD[j] = mul_add(A[i+j], D[i+j], AD[j]);
                BD[j] = mul_add(B[i+j], D[i+j], BD[j]);
                CD[j] = mul_add(C[i+j], D[i+j], CD[j]);
            }
		}
        if (ADD) {
            X += AD[0] + AD[1];
            Y += BD[0] + BD[1];
            Z += CD[0] + CD[1];
        } else {
            X  = AD[0] + AD[1];
            Y  = BD[0] + BD[1];
            Z  = CD[0] + CD[1];
        }
    } else {
        // odd states
        VectorClass AD[2], BD[2], CD[2];
        for (j = 0; j < 2; j++) {
            AD[j] = A[j] * D[j];
            BD[j] = B[j] * D[j];
            CD[j] = C[j] * D[j];
        }
		for (i = 2; i < N-1; i+=2) {
            for (j = 0; j < 2; j++) {
                AD[j] = mul_add(A[i+j], D[i+j], AD[j]);
                BD[j] = mul_add(B[i+j], D[i+j], BD[j]);
                CD[j] = mul_add(C[i+j], D[i+j], CD[j]);
            }
		}
        if (ADD) {
            X += mul_add(A[N-1], D[N-1], AD[0] + AD[1]);
            Y += mul_add(B[N-1], D[N-1], BD[0] + BD[1]);
            Z += mul_add(C[N-1], D[N-1], CD[0] + CD[1]);
        } else {
            X  = mul_add(A[N-1], D[N-1], AD[0] + AD[1]);
            Y  = mul_add(B[N-1], D[N-1], BD[0] + BD[1]);
            Z  = mul_add(C[N-1], D[N-1], CD[0] + CD[1]);
        }
    }
}


/**
    Given three vectors A, B, C, compute X:
    X = A.B.C = A[0]*B[0]*C[0] + ... + A[N-1]*B[N-1]*C[N-1]
    @param N number of elements
    @param A vector of size N
    @param B vector of size N
    @param C vector of size N
    @param[out] X = A.B.C
*/
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void dotProduct3Vec(const Numeric *A, const VectorClass *B, const VectorClass *C, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProduct3Vec(const Numeric *A, const VectorClass *B, const VectorClass *C, VectorClass &X, size_t N)
#endif
{
    size_t i, j;
    switch (N % 4) {
    case 0: {
        VectorClass V[4];
        for (j = 0; j < 4; j++)
            V[j] = A[j] * B[j] * C[j];
        for (i = 4; i < N; i+=4)
            for (j = 0; j < 4; j++)
                V[j] = mul_add(A[i+j]*B[i+j], C[i+j], V[j]);
        X = (V[0]+V[1])+(V[2]+V[3]);
        break;
    }

    case 1: {
        VectorClass V[4];
        const size_t Nm1 = N-1;
        for (j = 0; j < 4; j++)
            V[j] = A[j] * B[j] * C[j];
        for (i = 4; i < Nm1; i+=4)
            for (j = 0; j < 4; j++)
                V[j] = mul_add(A[i+j]*B[i+j], C[i+j], V[j]);
        X = mul_add(A[Nm1]*B[Nm1], C[Nm1], (V[0]+V[1])+(V[2]+V[3]));
        break;
    }

    case 2: {
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j] * C[j];
        for (i = 2; i < N; i+=2)
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j]*B[i+j], C[i+j], V[j]);
        X = (V[0]+V[1]);
        break;
    }

    default: {
        // odd states
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j] * C[j];
        for (i = 2; i < N-1; i+=2)
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j]*B[i+j], C[i+j], V[j]);
        X = mul_add(A[N-1]*B[N-1], C[N-1], V[0]+V[1]);
        break;
    }
    }
}


/**
    given three vectors A, B, C and a numeric coefficient D, compute X:
    X = exp(A[0]*D)*B[0]*C[0] + ... exp(A[N-1]*D)*B[N-1]*C[N-1]
    @param N number of elements
    @param A vector of size N
    @param B vector of size N
    @param C vector of size N
    @param D coefficient for A
    @param[out] X = exp(A[0]*D)*B[0]*C[0] + ... exp(A[N-1]*D)*B[N-1]*C[N-1]
*/
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void dotProductExp(const VectorClass *A, const VectorClass *B, const VectorClass *C,
                          const Numeric D, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductExp(const VectorClass *A, const VectorClass *B, const VectorClass *C,
                          const Numeric D, VectorClass &X, size_t N)
#endif
{
    size_t i;
    X = exp(A[0]*D)*B[0]*C[0];
    for (i = 1; i < N; i++)
        X = mul_add(exp(A[i]*D), B[i]*C[i], X);
}


/**
    given two vectors A, B and a numeric coefficient D, compute X:
    X = exp(A[0]*D)*B[0] + ... exp(A[N-1]*D)*B[N-1]
    @param N number of elements
    @param A vector of size N
    @param B vector of size N
    @param D coefficient for A
    @param[out] X = exp(A[0]*D)*B[0] + ... exp(A[N-1]*D)*B[N-1]
*/
#ifdef KERNEL_FIX_STATES
template <class VectorClass, class Numeric, const size_t N, const bool FMA>
inline void dotProductExp(const VectorClass *A, const VectorClass *B, Numeric D, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductExp(const VectorClass *A, const VectorClass *B, Numeric D, VectorClass &X, size_t N)
#endif
{
    size_t i;
    X = exp(A[0]*D)*B[0];
    for (i = 1; i < N; i++)
        X = mul_add(exp(A[i]*D), B[i], X);
}

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const size_t nstates>
inline void scaleLikelihood(VectorClass &lh_max, double *invar, double *dad_partial_lh, UBYTE *dad_scale_num,
    size_t ncat_mix)
#else
template <class VectorClass, const bool SAFE_NUMERIC>
inline void scaleLikelihood(VectorClass &lh_max, double *invar, double *dad_partial_lh, UBYTE *dad_scale_num,
    size_t ncat_mix, size_t nstates)
#endif
{
    if (SAFE_NUMERIC) {
        size_t x, i;
        auto underflown = ((lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(invar) == 0.0));
        if (horizontal_or(underflown)) { // at least one site has numerical underflown
            for (x = 0; x < VectorClass::size(); x++)
            if (underflown[x]) {
                // BQM 2016-05-03: only scale for non-constant sites
                // now do the likelihood scaling
                double *partial_lh = &dad_partial_lh[x];
                for (i = 0; i < nstates; i++)
                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                dad_scale_num[x*ncat_mix] += 1;
            }
        }
    } else {
        size_t x, i;
        auto underflown = (lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(invar) == 0.0);
        if (horizontal_or(underflown)) { // at least one site has numerical underflown
            size_t block = ncat_mix * nstates;
            for (x = 0; x < VectorClass::size(); x++)
            if (underflown[x]) {
                double *partial_lh = &dad_partial_lh[x];
                // now do the likelihood scaling
                for (i = 0; i < block; i++) {
                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                }
                dad_scale_num[x] += 1;
            }
        }
    }
}


/*******************************************************
 *
 * Helper function to pre-compute traversal information
 * and buffer to transition matrix
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template<class VectorClass>
void PhyloTree::computePartialInfoWrapper(TraversalInfo &info, double* buffer) {
    computePartialInfo(info, (VectorClass*)buffer);
}

#endif

#ifdef KERNEL_FIX_STATES
template<class VectorClass, const int nstates>
#else
template<class VectorClass>
#endif
void PhyloTree::computePartialInfo(TraversalInfo &info, VectorClass* buffer) {

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif

    size_t c, i, x;
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t block = nstates * ncat_mix;
    size_t tip_block = nstates * model->getNMixtures();
    size_t cat_id[ncat_mix], mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        cat_id[c] = c%ncat;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }
	double *evec = model->getEigenvectors();
	double *eval = model->getEigenvalues();

    PhyloNode* dad          = info.dad;
    PhyloNode* node         = info.dad_branch->getNode();
    double *echild          = info.echildren;
    double *partial_lh_leaf = info.partial_lh_leaves;

    //----------- Non-reversible model --------------

    if (!model->useRevKernel()) {
        size_t nstatesqr = nstates*nstates;
        // non-reversible model
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, child) {
            PhyloNode* childNode = child->getNode();
            // precompute information buffer
            if (child->direction == TOWARD_ROOT) {
                // transpose probability matrix
                double mat[nstatesqr];
                for (c = 0; c < ncat_mix; c++) {
                    double len_child = site_rate->getRate(c%ncat) * child->length;
                    model_factory->computeTransMatrix(len_child, mat, c/denom);
                    double *echild_ptr = &echild[c*nstatesqr];
                    for (i = 0; i < nstates; i++) {
                        for (x = 0; x < nstates; x++)
                            echild_ptr[x] = mat[x*nstates+i];
                        echild_ptr += nstates;
                    }
                }
            } else {
                for (c = 0; c < ncat_mix; c++) {
                    double len_child = site_rate->getRate(c%ncat) * child->length;
                    model_factory->computeTransMatrix(len_child, &echild[c*nstatesqr], c/denom);
                }
            }

            // pre compute information for tip
            if (isRootLeaf(childNode)) {
                for (c = 0; c < ncat_mix; c++) {
                    size_t m = c/denom;
                    model->getStateFrequency(partial_lh_leaf + c*nstates, m);
                }
                partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
            } else if (childNode->isLeaf()) {
                //vector<int>::iterator it;
                if (nstates % VectorClass::size() == 0) {
                    // vectorized version
                    for (int state = 0; state <= aln->STATE_UNKNOWN; state++) {
                        VectorClass *this_tip_partial_lh = (VectorClass*)&tip_partial_lh[state*nstates];
                        double *this_partial_lh_leaf = &partial_lh_leaf[state*block];
                        VectorClass *echild_ptr = (VectorClass*)echild;
                        for (x = 0; x < block; x++) {
                            VectorClass vchild = echild_ptr[0] * this_tip_partial_lh[0];
                            for (i = 1; i < nstates/VectorClass::size(); i++)
                                vchild += echild_ptr[i] * this_tip_partial_lh[i];
                            echild_ptr += nstates/VectorClass::size();
                            this_partial_lh_leaf[x] = horizontal_add(vchild);
                        }
                    }
                } else {
                    // non-vectorized version
                    for (int state = 0; state <= aln->STATE_UNKNOWN; state++) {
                        double* this_tip_partial_lh  = &tip_partial_lh[state*nstates];
                        double* this_partial_lh_leaf = &partial_lh_leaf[state*block];
                        double* echild_ptr = echild;
                        for (x = 0; x < block; x++) {
                            double vchild = 0.0;
                            for (i = 0; i < nstates; i++) {
                                vchild += echild_ptr[i] * this_tip_partial_lh[i];
                            }
                            echild_ptr += nstates;
                            this_partial_lh_leaf[x] = vchild;
                        }
                    }
                }
                partial_lh_leaf += aln->STATE_UNKNOWN * block;
                for (x = 0; x < block; x++) {
                    partial_lh_leaf[x] = 1.0;
                }
                partial_lh_leaf += block;
            }
            echild += block*nstates;
        }
        return;
    } // END non-reversible model

    //----------- Reversible model --------------
    if (nstates % VectorClass::size() == 0) {
        // vectorized version
        VectorClass *expchild = (VectorClass*)buffer;
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, child) {
            PhyloNode*   childNode  = child->getNode();
            VectorClass* echild_ptr = (VectorClass*)echild;
            // precompute information buffer
            for (c = 0; c < ncat_mix; c++) {
                VectorClass len_child = site_rate->getRate(cat_id[c]) * child->getLength(cat_id[c]);
                double *eval_ptr = eval + mix_addr_nstates[c];
                double *evec_ptr = evec + mix_addr[c];
                for (i = 0; i < nstates/VectorClass::size(); i++) {
                    // eval is not aligned!
                    expchild[i] = exp(VectorClass().load_a(&eval_ptr[i*VectorClass::size()]) * len_child);
                }
                for (x = 0; x < nstates; x++) {
                    for (i = 0; i < nstates/VectorClass::size(); i++) {
                        // evec is not be aligned!
                        echild_ptr[i] = (VectorClass().load_a(&evec_ptr[x*nstates+i*VectorClass::size()]) * expchild[i]);
                    }
                    echild_ptr += nstates/VectorClass::size();
                }
            }
            // pre compute information for tip
            if (childNode->isLeaf()) {
                for (int state = 0; state <= aln->STATE_UNKNOWN; state++) {
                    double *this_partial_lh_leaf = partial_lh_leaf + state*block;
                    VectorClass *echild_ptr = (VectorClass*)echild;
                    for (c = 0; c < ncat_mix; c++) {
                        VectorClass *this_tip_partial_lh = (VectorClass*)(tip_partial_lh + state*tip_block + mix_addr_nstates[c]);
                        for (x = 0; x < nstates; x++) {
                            VectorClass vchild = echild_ptr[0] * this_tip_partial_lh[0];
                            for (i = 1; i < nstates/VectorClass::size(); i++) {
                                vchild = mul_add(echild_ptr[i], this_tip_partial_lh[i], vchild);
                            }
                            this_partial_lh_leaf[x] = horizontal_add(vchild);
                            echild_ptr += nstates/VectorClass::size();
                        }
                        this_partial_lh_leaf += nstates;
                    }
                }
                size_t addr = aln->STATE_UNKNOWN * block;
                for (x = 0; x < block; x++) {
                    partial_lh_leaf[addr+x] = 1.0;
                }
                partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
            }
            echild += block*nstates;
        }
    } else {
        // non-vectorized version
        double expchild[nstates];
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, child) {
            // precompute information buffer
            double *echild_ptr = echild;
            for (c = 0; c < ncat_mix; c++) {
                double len_child = site_rate->getRate(cat_id[c]) * child->getLength(cat_id[c]);
                double *eval_ptr = eval + mix_addr_nstates[c];
                double *evec_ptr = evec + mix_addr[c];
                for (i = 0; i < nstates; i++) {
                    expchild[i] = exp(eval_ptr[i]*len_child);
                }
                for (x = 0; x < nstates; x++) {
                    for (i = 0; i < nstates; i++) {
                        echild_ptr[i] = evec_ptr[x*nstates+i] * expchild[i];
                    }
                    echild_ptr += nstates;
                }
            }
            // pre compute information for tip
            PhyloNode* childNode = child->getNode();
            if (childNode->isLeaf()) {
                //vector<int>::iterator it;
                for (int state = 0; state <= aln->STATE_UNKNOWN; state++) {
                    double *this_partial_lh_leaf = partial_lh_leaf + state*block;
                    double *echild_ptr = echild;
                    for (c = 0; c < ncat_mix; c++) {
                        double *this_tip_partial_lh = tip_partial_lh + state*tip_block + mix_addr_nstates[c];
                        for (x = 0; x < nstates; x++) {
                            double vchild = echild_ptr[0] * this_tip_partial_lh[0];
                            for (i = 1; i < nstates; i++) {
                                vchild += echild_ptr[i] * this_tip_partial_lh[i];
                            }
                            this_partial_lh_leaf[x] = vchild;
                            echild_ptr += nstates;
                        }
                        this_partial_lh_leaf += nstates;
                    }
                }
                size_t addr = aln->STATE_UNKNOWN * block;
                for (x = 0; x < block; x++) {
                    partial_lh_leaf[addr+x] = 1.0;
                }
                partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
            }
            echild += block*nstates;
        }
    }
}

#ifdef KERNEL_FIX_STATES
template<class VectorClass, const int nstates>
#else
template<class VectorClass>
#endif
void PhyloTree::computeTraversalInfo(PhyloNode *node, PhyloNode *dad,
                                     LikelihoodBufferSet& buffers, bool compute_partial_lh) {

    if ((tip_partial_lh_computed & 1) == 0)
    {
        computeTipPartialLikelihood();
    }
    traversal_info.clear();
#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    // reserve beginning of buffer_partial_lh for other purpose
    size_t ncat_mix = (model_factory->fused_mix_rate) ? site_rate->getNRate() : site_rate->getNRate()*model->getNMixtures();
    size_t block = aln->num_states * ncat_mix;
    double *buffer = buffers.buffer_partial_lh + block*VectorClass::size()*num_packets + get_safe_upper_limit(block)*(aln->STATE_UNKNOWN+2);

    // more buffer for non-reversible models
    if (!model->useRevKernel()) {
        buffer += get_safe_upper_limit(3*block*nstates);
        buffer += get_safe_upper_limit(block)*(aln->STATE_UNKNOWN+1)*2;
        buffer += block*2*VectorClass::size()*num_packets;
    }

    // sort subtrees for mem save technique
    if (params->lh_mem_save == LM_MEM_SAVE) {
        int node_size = node->computeSize(dad);
        int dad_size  = dad->computeSize(node);
        if (node_size < dad_size) {
            // swap node and dad due to tree size
            std::swap(node, dad);
        }
    }

    PhyloNeighbor* dad_branch  = dad->findNeighbor(node);
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    bool dad_locked  = computeTraversalInfo(dad_branch,  dad,  buffer);
    bool node_locked = computeTraversalInfo(node_branch, node, buffer);
    if (params->lh_mem_save == LM_MEM_SAVE) {
        if (dad_locked) {
            mem_slots.unlock(dad_branch);
        }
        if (node_locked) {
            mem_slots.unlock(node_branch);
        }
    }
    
    if (verbose_mode >= VB_DEBUG && traversal_info.size() > 0) {
        Node *saved = root;
        root = dad;
        hideProgress();
        drawTree(cout);
        showProgress();
        root = saved;
    }

    if (traversal_info.empty()) {
        return;
    }

    if (!model->isSiteSpecificModel()) {
        int num_info = traversal_info.size();
        if (verbose_mode >= VB_DEBUG) {
            hideProgress();
            cout << "traversal order:";
            for (auto it = traversal_info.begin(); it != traversal_info.end(); it++) {
                cout << "  ";
                if (it->dad->isLeaf())
                    cout << it->dad->name;
                else
                    cout << it->dad->id;
                cout << "->";
                if (it->dad_branch->node->isLeaf())
                    cout << it->dad_branch->node->name;
                else
                    cout << it->dad_branch->node->id;
                if (params->lh_mem_save == LM_MEM_SAVE) {
                    if (it->dad_branch->isLikelihoodComputed())
                        cout << " [";
                    else
                        cout << " (";
                    cout << mem_slots.findNei(it->dad_branch) - mem_slots.begin();
                    if (it->dad_branch->isLikelihoodComputed())
                        cout << "]";
                    else
                        cout << ")";
                }
            }
            cout << endl;
            showProgress();
        }
        
        if (!Params::getInstance().buffer_mem_save) {
#ifdef _OPENMP
#pragma omp parallel if (num_info >= 3) num_threads(num_threads)
        {
            VectorClass *buffer_tmp = (VectorClass*)buffer + aln->num_states*omp_get_thread_num();
#pragma omp for schedule(static)
#else
            VectorClass *buffer_tmp = (VectorClass*)buffer;
#endif
            for (int i = 0; i < num_info; i++) {
            #ifdef KERNEL_FIX_STATES
                computePartialInfo<VectorClass, nstates>(traversal_info[i], buffer_tmp);
            #else
                computePartialInfo<VectorClass>(traversal_info[i], buffer_tmp);
            #endif
            }
#ifdef _OPENMP
        }
#endif
        }
    }

    if (compute_partial_lh) {
        vector<size_t> limits;
        size_t orig_nptn = roundUpToMultiple(aln->size(), VectorClass::size());
        size_t nptn      = roundUpToMultiple(orig_nptn+model_factory->unobserved_ptns.size(),VectorClass::size());
        computePatternPacketBounds(VectorClass::size(), num_threads,
                                   num_packets, nptn, limits);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,1) num_threads(num_threads)
        #endif
        for (int packet_id = 0; packet_id < num_packets; ++packet_id) {
            for (auto it = traversal_info.begin(); it != traversal_info.end(); it++) {
                computePartialLikelihood(*it, limits[packet_id], limits[packet_id+1], packet_id, buffers);
            }
        }
        traversal_info.clear();
    }
    return;
}

/*******************************************************
 *
 * NEW! highly-vectorized partial likelihood function
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computePartialLikelihoodSIMD(TraversalInfo &info,
                                             size_t ptn_lower, size_t ptn_upper, int packet_id,
                                             LikelihoodBufferSet& buffers)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computePartialLikelihoodGenericSIMD(TraversalInfo &info,
                                                    size_t ptn_lower, size_t ptn_upper, int packet_id,
                                                    LikelihoodBufferSet& buffers)
#endif
{
    PhyloNeighbor* dad_branch = info.dad_branch;
    PhyloNode*     dad        = info.dad;
    ASSERT(dad);
    PhyloNode*     node       = dad_branch->getNode();
    if (node->isLeaf()) {
        return;
    }

#ifndef KERNEL_FIX_STATES
    size_t nstates       = aln->num_states;
#endif
    const size_t states_square = nstates*nstates;
    size_t orig_nptn     = aln->size();
    size_t max_orig_nptn = roundUpToMultiple(orig_nptn, VectorClass::size());
    size_t nptn          = max_orig_nptn+model_factory->unobserved_ptns.size();

    size_t ncat          = site_rate->getNRate();
    size_t ncat_mix      = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t denom         = (model_factory->fused_mix_rate) ? 1 : ncat;
    
    size_t mix_addr_nstates[ncat_mix];
    size_t mix_addr[ncat_mix];
    for (size_t c = 0; c < ncat_mix; c++) {
        size_t m            = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c]         = mix_addr_nstates[c]*nstates;
    }
    size_t  block        = nstates * ncat_mix;
    size_t  tip_mem_size = max_orig_nptn * nstates;
    size_t  scale_size   = SAFE_NUMERIC ? (ptn_upper-ptn_lower) * ncat_mix : (ptn_upper-ptn_lower);
    size_t  num_leaves   = 0;
    double* eval         = model->getEigenvalues();
    double* evec         = model->getEigenvectors();
    double* inv_evec     = model->getInverseEigenvectors();
    ASSERT( inv_evec!=nullptr && evec!=nullptr );

    // internal node
    PhyloNeighbor* left  = nullptr;
    PhyloNeighbor* right = nullptr; // left & right are two neighbors leading to 2 subtrees
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        // make sure that the partial_lh of children are different!
        ASSERT(dad_branch->partial_lh != nei->partial_lh);
        if (!left) left = nei; else right = nei;
        if (nei->node->isLeaf()) {
            num_leaves++;
        }
    }

    // precomputed buffer to save times
    size_t thread_buf_size        = (2*block+nstates)*VectorClass::size();
    double* buffer_partial_lh_ptr = buffers.buffer_partial_lh
                                  + (getBufferPartialLhSize() - thread_buf_size*num_packets);
    double* echildren             = nullptr;
    double* partial_lh_leaves     = nullptr;

    // pre-compute scaled branch length per category
    double  len_children[ncat*(node->degree()-1)]; // +1 in case num_leaves = 0
    double* len_left  = nullptr;
    double* len_right = nullptr;

    if (SITE_MODEL) {
        double *len_children_ptr = len_children;
        FOR_NEIGHBOR_IT(node, dad, it3) {
            for (size_t c = 0; c < ncat; c++) {
                len_children_ptr[c] = site_rate->getRate(c) * (*it3)->length;
            }
            if (!len_left) {
                len_left  = len_children_ptr;
            } else {
                len_right = len_children_ptr;
            }
            len_children_ptr += ncat;
        }
    } else {
        auto children_size = get_safe_upper_limit(block*nstates*(node->degree()-1));
        auto lh_leaf_size  = get_safe_upper_limit((aln->STATE_UNKNOWN+1)*block*num_leaves);
        if (Params::getInstance().buffer_mem_save) {            
            echildren = aligned_alloc<double>(children_size);
            if (num_leaves > 0) {
                partial_lh_leaves = aligned_alloc<double>(lh_leaf_size);
            } else {
                partial_lh_leaves  = nullptr;
            }
            double *buffer_tmp = aligned_alloc<double>(nstates);
#ifdef KERNEL_FIX_STATES
            computePartialInfo<VectorClass, nstates>(info, (VectorClass*)buffer_tmp);
#else
            computePartialInfo<VectorClass>(info, (VectorClass*)buffer_tmp);
#endif
            aligned_free(buffer_tmp);
        } else {
            echildren         = info.echildren;
            partial_lh_leaves = info.partial_lh_leaves;
        }
    }

    double*    eleft     = echildren;
    double*    eright    = echildren + block * nstates;
    PhyloNode* leftNode  = left->getNode();
    PhyloNode* rightNode = right->getNode();
	if (!leftNode->isLeaf() && rightNode->isLeaf()) {
        std::swap(left,     right);
        std::swap(eleft,    eright);
        std::swap(len_left, len_right);
        std::swap(leftNode, rightNode);
    }
    if (node->degree() > 3) {
        /*--------------------- multifurcating node ------------------*/

        // now for-loop computing partial_lh over all site-patterns
        VectorClass* partial_lh_all = (VectorClass*) &buffer_partial_lh_ptr[thread_buf_size * packet_id];
        double*      vec_tip        = (double*)&partial_lh_all[block];

        for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            for (size_t i = 0; i < block; i++){
                partial_lh_all[i] = 1.0;
            }
            UBYTE *scale_dad = NULL;
            if (SAFE_NUMERIC) {
                scale_dad = dad_branch->scale_num + ptn*ncat_mix;
                memset(scale_dad, 0, sizeof(UBYTE)*ncat_mix*VectorClass::size());
            } else {
                memset(&dad_branch->scale_num[ptn], 0, sizeof(UBYTE)*VectorClass::size());
            }

            // SITE_MODEL variables
            VectorClass*       expchild = partial_lh_all + block;
            const VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
            const VectorClass* evec_ptr = (VectorClass*) &evec[ptn*states_square];
            double*            len_child = len_children;
            VectorClass        vchild;

            // normal model
            double* partial_lh_leaf = partial_lh_leaves;
            double* echild          = echildren;

            FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, child) {
                PhyloNode* childNode = child->getNode();
                if (SITE_MODEL) {
                    UBYTE *scale_child = SAFE_NUMERIC ? child->scale_num + ptn*ncat_mix : NULL;
                    VectorClass *partial_lh = partial_lh_all;
                    if (childNode->isLeaf()) {
                        // external node
                        const VectorClass* tip_partial_lh_child = (VectorClass*) &tip_partial_lh[childNode->id*tip_mem_size + ptn*nstates];
                        for (size_t c = 0; c < ncat; c++) {
                            for (size_t i = 0; i < nstates; i++) {
                                expchild[i] = exp(eval_ptr[i]*len_child[c]) * tip_partial_lh_child[i];
                            }
                            for (size_t x = 0; x < nstates; x++) {
                                const VectorClass* this_evec = &evec_ptr[x*nstates];
#ifdef KERNEL_FIX_STATES
                                dotProductVec<VectorClass, VectorClass, nstates, FMA>(expchild, this_evec, vchild);
#else
                                dotProductVec<VectorClass, VectorClass, FMA>(expchild, this_evec, vchild, nstates);
#endif
                                partial_lh[x] *= vchild;
                            }
                            partial_lh += nstates;
                        }
                    } else {
                        // internal node
                        VectorClass *partial_lh = partial_lh_all;
                        VectorClass *partial_lh_child = (VectorClass*)(child->partial_lh + ptn*block);
                        if (!SAFE_NUMERIC) {
                            for (size_t i = 0; i < VectorClass::size(); i++) {
                                dad_branch->scale_num[ptn+i] += child->scale_num[ptn+i];
                            }
                        }
                        for (size_t c = 0; c < ncat_mix; c++) {
                            if (SAFE_NUMERIC) {
                                for (size_t x = 0; x < VectorClass::size(); x++)
                                    scale_dad[x*ncat_mix+c] += scale_child[x*ncat_mix+c];
                            }
                            // compute real partial likelihood vector
                            for (size_t i = 0; i < nstates; i++) {
                                expchild[i] = exp(eval_ptr[i]*len_child[c]) * partial_lh_child[i];
                            }
                            for (size_t x = 0; x < nstates; x++) {
                                const VectorClass *this_evec = &evec_ptr[x*nstates];
#ifdef KERNEL_FIX_STATES
                                dotProductVec<VectorClass, VectorClass, nstates, FMA>(expchild, this_evec, vchild);
#else
                                dotProductVec<VectorClass, VectorClass, FMA>(expchild, this_evec, vchild, nstates);
#endif
                                partial_lh[x] *= vchild;
                            }
                            partial_lh += nstates;
                            partial_lh_child += nstates;
                        }
                    } // if
                    len_child += ncat;
                } else {
                    // non site specific model
                    auto   stateRow    = this->getConvertedSequenceByNumber(childNode->id);
                    auto   unknown     = aln->STATE_UNKNOWN;
                    UBYTE* scale_child = SAFE_NUMERIC ? child->scale_num + ptn*ncat_mix : NULL;
                    if (childNode->isLeaf()) {
                        // external node
                        // load data for tip
                        for (size_t i = 0; i < VectorClass::size(); i++) {
                            int state;
                            if (ptn+i < orig_nptn) {
                                if (stateRow!=nullptr) {
                                    state = stateRow[ptn+i];
                                } else {
                                    state = (aln->at(ptn+i))[childNode->id];
                                }
                            } else if (ptn+i < max_orig_nptn) {
                                state = unknown;
                            } else if (ptn+i < nptn) {
                                state = model_factory->unobserved_ptns[ptn+i-max_orig_nptn][childNode->id];
                            } else {
                                state = unknown;
                            }
                            double* child_lh     = partial_lh_leaf + block*state;
                            double* this_vec_tip = vec_tip+i;
                            for (size_t c = 0; c < block; c++) {
                                *this_vec_tip = child_lh[c];
                                this_vec_tip += VectorClass::size();
                            }
                        }
                        VectorClass *vtip = (VectorClass*)vec_tip;
                        for (size_t c = 0; c < block; c++) {
                            // compute real partial likelihood vector
                            partial_lh_all[c] *= vtip[c];
                        }
                        partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
                    } else {
                        // internal node
                        VectorClass* partial_lh       = partial_lh_all;
                        VectorClass* partial_lh_child = (VectorClass*)(child->partial_lh + ptn*block);
                        if (!SAFE_NUMERIC) {
                            for (size_t i = 0; i < VectorClass::size(); i++) {
                                dad_branch->scale_num[ptn+i] += child->scale_num[ptn+i];
                            }
                        }
                        double *echild_ptr = echild;
                        for (size_t c = 0; c < ncat_mix; c++) {
                            if (SAFE_NUMERIC) {
                                for (size_t x = 0; x < VectorClass::size(); x++) {
                                    scale_dad[x*ncat_mix+c] += scale_child[x*ncat_mix+c];
                                }
                            }
                            // compute real partial likelihood vector
                            for (size_t x = 0; x < nstates; x++) {
                                VectorClass vchild = echild_ptr[0] * partial_lh_child[0];
    //                            double *echild_ptr = echild + (c*nstatesqr+x*nstates);
                                for (size_t i = 1; i < nstates; i++) {
                                    vchild = mul_add(echild_ptr[i], partial_lh_child[i], vchild);
                                }
                                echild_ptr += nstates;
                                partial_lh[x] *= vchild;
                            }
                            partial_lh += nstates;
                            partial_lh_child += nstates;
                        }
                    } // if
                    echild += block*nstates;
                } // if SITE_MODEL

                /***** now do likelihood rescaling ******/
                if (SAFE_NUMERIC) {
                    VectorClass *partial_lh_tmp = partial_lh_all;
                    for (size_t c = 0; c < ncat_mix; c++) {
                        VectorClass lh_max = 0.0;
                        for (size_t x = 0; x < nstates; x++) {
                            lh_max = max(lh_max,abs(partial_lh_tmp[x]));
                        }
                        // check if one should scale partial likelihoods
                        auto underflown = ((lh_max < SCALING_THRESHOLD)
                                           & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                        if (horizontal_or(underflown)) { // at least one site has numerical underflown
                            for (size_t x = 0; x < VectorClass::size(); x++)
                            if (underflown[x]) {
                                // BQM 2016-05-03: only scale for non-constant sites
                                // now do the likelihood scaling
                                double *partial_lh = (double*)partial_lh_tmp + (x);
                                for (size_t i = 0; i < nstates; i++)
                                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                                dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
                            }
                        }
                        partial_lh_tmp += nstates;
                    }
                } else {
                    // not -safe numeric
                    VectorClass lh_max = 0.0;
                    for (size_t x = 0; x < block; x++)
                        lh_max = max(lh_max,abs(partial_lh_all[x]));
                    auto underflown = (lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                    if (horizontal_or(underflown)) { // at least one site has numerical underflown
                        for (size_t x = 0; x < VectorClass::size(); x++) {
                            if (underflown[x]) {
                                double *partial_lh = (double*)partial_lh_all + (x);
                                // now do the likelihood scaling
                                for (size_t i = 0; i < block; i++) {
                                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                                }
                                // sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                                dad_branch->scale_num[ptn+x] += 1;
                            }
                        }
                    }
                } // if-else
            } // FOR_NEIGHBOR
        
            // compute dot-product with inv_eigenvector
            VectorClass *partial_lh_tmp = partial_lh_all;
            VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;
            double *inv_evec_ptr = SITE_MODEL ? &inv_evec[ptn*states_square] : NULL;
            for (size_t c = 0; c < ncat_mix; c++) {
                if (SITE_MODEL) {
                    // compute dot-product with inv_eigenvector
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, VectorClass, nstates, FMA>(partial_lh_tmp,
                                                                          (VectorClass*)inv_evec_ptr,
                                                                          partial_lh);
#else
                    productVecMat<VectorClass, VectorClass, FMA> (partial_lh_tmp,
                                                                  (VectorClass*)inv_evec_ptr,
                                                                  partial_lh, nstates);
#endif
                } else {
                    inv_evec_ptr = inv_evec + mix_addr[c];
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, double, nstates, FMA>(partial_lh_tmp, inv_evec_ptr,
                                                                     partial_lh, lh_max);
#else
                    productVecMat<VectorClass, double, FMA> (partial_lh_tmp, inv_evec_ptr,
                                                             partial_lh, lh_max, nstates);
#endif
                }
                partial_lh += nstates;
                partial_lh_tmp += nstates;
            }
        } // for ptn

        // end multifurcating treatment
    } else if ( leftNode->isLeaf() && rightNode->isLeaf() ) {
        /*--------------------- TIP-TIP (cherry) case ------------------*/

        const double* partial_lh_left  = SITE_MODEL
            ? &tip_partial_lh[leftNode->id  * tip_mem_size]
            : partial_lh_leaves;
        const double* partial_lh_right = SITE_MODEL
            ? &tip_partial_lh[rightNode->id * tip_mem_size]
            : partial_lh_leaves + (aln->STATE_UNKNOWN+1)*block;

		// scale number must be ZERO
	    memset(dad_branch->scale_num + (SAFE_NUMERIC ? ptn_lower*ncat_mix : ptn_lower), 0, scale_size * sizeof(UBYTE));
        double *vec_left   = buffer_partial_lh_ptr + thread_buf_size * packet_id;

        double *vec_right  =  SITE_MODEL ? &vec_left[nstates*VectorClass::size()] : &vec_left[block*VectorClass::size()];
        VectorClass *partial_lh_tmp = SITE_MODEL ? (VectorClass*)vec_right+nstates : (VectorClass*)vec_right+block;

        auto leftStateRow  = this->getConvertedSequenceByNumber(leftNode->id);
        auto rightStateRow = this->getConvertedSequenceByNumber(rightNode->id);
        auto unknown       = aln->STATE_UNKNOWN;

        for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);

            if (SITE_MODEL) {
                VectorClass* expleft  = (VectorClass*) vec_left;
                VectorClass* expright = (VectorClass*) vec_right;
                VectorClass* vleft    = (VectorClass*) &partial_lh_left[ptn*nstates];
                VectorClass* vright   = (VectorClass*) &partial_lh_right[ptn*nstates];
                const VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                const VectorClass* evec_ptr = (VectorClass*) &evec[ptn*states_square];
                const VectorClass* inv_evec_ptr = (VectorClass*) &inv_evec[ptn*states_square];
                for (size_t c = 0; c < ncat; c++) {
                    for (size_t i = 0; i < nstates; i++) {
                        expleft[i]  = exp(eval_ptr[i]*len_left[c])  * vleft[i];
                        expright[i] = exp(eval_ptr[i]*len_right[c]) * vright[i];
                    }
                    // compute real partial likelihood vector
                    for (size_t x = 0; x < nstates; x++) {
                        const VectorClass* this_evec = evec_ptr + x*nstates;
#ifdef KERNEL_FIX_STATES
                        dotProductDualVec<VectorClass, VectorClass, nstates, FMA>(this_evec, expleft, this_evec, expright, partial_lh_tmp[x]);
#else
                        dotProductDualVec<VectorClass, VectorClass, FMA>(this_evec, expleft, this_evec, expright, partial_lh_tmp[x], nstates);
#endif
                    }
                    // compute dot-product with inv_eigenvector
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, VectorClass, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh);
#else
                    productVecMat<VectorClass, VectorClass, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, nstates);
#endif
                    partial_lh += nstates;
                } // FOR category
            } else {
                VectorClass* vleft  = (VectorClass*)vec_left;
                VectorClass* vright = (VectorClass*)vec_right;
                // load data for tip
                for (size_t x = 0; x < VectorClass::size(); x++) {
                    int leftState;
                    int rightState;
                    if (ptn+x < orig_nptn) {
                        if (leftStateRow!=nullptr) {
                            leftState = leftStateRow[ptn+x];
                        } else {
                            leftState = (aln->at(ptn+x))[leftNode->id];
                        }
                        if (rightStateRow!=nullptr) {
                            rightState =  rightStateRow[ptn+x];
                        } else {
                            rightState = (aln->at(ptn+x))[rightNode->id];
                        }
                    } else if (ptn+x < max_orig_nptn) {
                        leftState  = unknown;
                        rightState = unknown;
                    } else if (ptn+x < nptn) {
                        leftState  = model_factory->unobserved_ptns[ptn+x-max_orig_nptn][leftNode->id];
                        rightState = model_factory->unobserved_ptns[ptn+x-max_orig_nptn][rightNode->id];
                    } else {
                        leftState  = unknown;
                        rightState = unknown;
                    }
                    const double* tip_left       = partial_lh_left  + block*leftState;
                    const double* tip_right      = partial_lh_right + block*rightState;
                    double*       this_vec_left  = vec_left+x;
                    double*       this_vec_right = vec_right+x;
                    for (size_t i = 0; i < block; i++) {
                        *this_vec_left  = tip_left[i];
                        *this_vec_right = tip_right[i];
                        this_vec_left  += VectorClass::size();
                        this_vec_right += VectorClass::size();
                    }
                }

                for (size_t c = 0; c < ncat_mix; c++) {
                    double *inv_evec_ptr = inv_evec + mix_addr[c];
                    // compute real partial likelihood vector
                    for (size_t x = 0; x < nstates; x++) {
                        partial_lh_tmp[x] = vleft[x] * vright[x];
                    }

                    // compute dot-product with inv_eigenvector
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, double, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh);
#else
                    productVecMat<VectorClass, double, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, nstates);
#endif

                    // increase pointer
                    vleft      += nstates;
                    vright     += nstates;
                    partial_lh += nstates;
                } // FOR category
            } // IF SITE_MODEL
        } // FOR LOOP
    } else if (leftNode->isLeaf() && !rightNode->isLeaf()) {

        /*--------------------- TIP-INTERNAL NODE case ------------------*/

        // only take scale_num from the right subtree
        memcpy(
            dad_branch->scale_num + (SAFE_NUMERIC ? ptn_lower*ncat_mix : ptn_lower),
            right->scale_num + (SAFE_NUMERIC ? ptn_lower*ncat_mix : ptn_lower),
            scale_size * sizeof(UBYTE));

        double* partial_lh_left = SITE_MODEL
            ? &tip_partial_lh[leftNode->id * tip_mem_size]
            : partial_lh_leaves;

        double* vec_left = buffer_partial_lh_ptr + thread_buf_size * packet_id;
        VectorClass* partial_lh_tmp = SITE_MODEL
            ? (VectorClass*)vec_left+2*nstates
            : (VectorClass*)vec_left+block;

        auto leftStateRow = this->getConvertedSequenceByNumber(leftNode->id);
        auto unknown      = aln->STATE_UNKNOWN;
        
        for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            VectorClass* partial_lh       = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            const VectorClass* partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;

            if (SITE_MODEL) {
                VectorClass *expleft  = (VectorClass*)vec_left;
                VectorClass *expright = expleft+nstates;
                VectorClass *vleft    = (VectorClass*)&partial_lh_left[ptn*nstates];
                const VectorClass *eval_ptr = (VectorClass*) &eval[ptn*nstates];
                const VectorClass *evec_ptr = (VectorClass*) &evec[ptn*states_square];
                const VectorClass *inv_evec_ptr = (VectorClass*) &inv_evec[ptn*states_square];
                for (size_t c = 0; c < ncat; c++) {
                    for (size_t i = 0; i < nstates; i++) {
                        expleft[i]  = exp(eval_ptr[i]*len_left[c]) * vleft[i];
                        expright[i] = exp(eval_ptr[i]*len_right[c]) * partial_lh_right[i];
                    }
                    // compute real partial likelihood vector
                    for (size_t x = 0; x < nstates; x++) {
                        const VectorClass *this_evec = evec_ptr + x*nstates;
#ifdef KERNEL_FIX_STATES
                        dotProductDualVec<VectorClass, VectorClass, nstates, FMA>(this_evec, expleft, this_evec, expright, partial_lh_tmp[x]);
#else
                        dotProductDualVec<VectorClass, VectorClass, FMA>(this_evec, expleft, this_evec, expright, partial_lh_tmp[x], nstates);
#endif
                    }
                    // compute dot-product with inv_eigenvector
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, VectorClass, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max);
#else
                    productVecMat<VectorClass, VectorClass, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max, nstates);
#endif
                    // check if one should scale partial likelihoods
                    if (SAFE_NUMERIC) {
                        auto underflown = ((lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                        if (horizontal_or(underflown)) { // at least one site has numerical underflown
                            for (size_t x = 0; x < VectorClass::size(); x++)
                            if (underflown[x]) {
                                // BQM 2016-05-03: only scale for non-constant sites
                                // now do the likelihood scaling
                                double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                                for (size_t i = 0; i < nstates; i++) {
                                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                                }
                                dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
                            }
                        }
                    }
                    partial_lh_right += nstates;
                    partial_lh += nstates;
                } // FOR category
            } else {
                //Not SITE_MODEL
                VectorClass *vleft = (VectorClass*)vec_left;
                // load data for tip
                for (size_t x = 0; x < VectorClass::size(); x++) {
                    int state;
                    if (ptn+x < orig_nptn) {
                        if (leftStateRow!=nullptr) {
                            state =  leftStateRow[ptn+x];
                        } else {
                            state = (aln->at(ptn+x))[leftNode->id];
                        }
                    } else if (ptn+x < max_orig_nptn) {
                        state = unknown;
                    } else if (ptn+x < nptn) {
                        state = model_factory->unobserved_ptns[ptn+x-max_orig_nptn][leftNode->id];
                    } else {
                        state = unknown;
                    }
                    double *tip = partial_lh_left + block*state;
                    double *this_vec_left = vec_left+x;
                    for (size_t i = 0; i < block; i++) {
                        *this_vec_left = tip[i];
                        this_vec_left += VectorClass::size();
                    }
                }

                double *eright_ptr = eright;
                for (size_t c = 0; c < ncat_mix; c++) {
                    if (SAFE_NUMERIC) {
                        lh_max = 0.0;
                    }
                    double *inv_evec_ptr = inv_evec + mix_addr[c];
                    // compute real partial likelihood vector
                    for (size_t x = 0; x < nstates; x++) {
                        VectorClass vright;
    #ifdef KERNEL_FIX_STATES
                        dotProductVec<VectorClass, double, nstates, FMA>(eright_ptr, partial_lh_right, vright);
    #else
                        dotProductVec<VectorClass, double, FMA>(eright_ptr, partial_lh_right, vright, nstates);
    #endif
                        eright_ptr += nstates;
                        partial_lh_tmp[x] = vleft[x] * (vright);
                    }

                    // compute dot-product with inv_eigenvector
    #ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, double, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max);
    #else
                    productVecMat<VectorClass, double, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max, nstates);
    #endif
                    // check if one should scale partial likelihoods
                    if (SAFE_NUMERIC) {
                        auto underflown = ((lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                        if (horizontal_or(underflown)) { // at least one site has numerical underflown
                            for (size_t x = 0; x < VectorClass::size(); x++) {
                                if (underflown[x]) {
                                    // BQM 2016-05-03: only scale for non-constant sites
                                    // now do the likelihood scaling
                                    double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                                    for (size_t i = 0; i < nstates; i++)
                                        partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                                    dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
                                }
                            }
                        }
                    }
                    vleft += nstates;
                    partial_lh_right += nstates;
                    partial_lh += nstates;
                } // FOR category
            } // IF SITE_MODEL

            if (!SAFE_NUMERIC) {
                auto underflown = (lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                if (horizontal_or(underflown)) { // at least one site has numerical underflown
                    for (size_t x = 0; x < VectorClass::size(); x++)
                    if (underflown[x]) {
                        double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                        // now do the likelihood scaling
                        for (size_t i = 0; i < block; i++) {
                            size_t j      = i*VectorClass::size();
                            partial_lh[j] = ldexp(partial_lh[j], SCALING_THRESHOLD_EXP);
                        }
//                        sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                        dad_branch->scale_num[ptn+x] += 1;
                    }
                }
            }
        } // big for loop over ptn
    } else {

        /*--------------------- INTERNAL-INTERNAL NODE case ------------------*/

        VectorClass *partial_lh_tmp
            = (VectorClass*)(buffer_partial_lh_ptr + thread_buf_size * packet_id);
		for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            size_t             blockOffset      = ptn*block;
			VectorClass*       partial_lh       = (VectorClass*)(dad_branch->partial_lh + blockOffset);
			const VectorClass* partial_lh_left  = (VectorClass*)(left->partial_lh       + blockOffset);
			const VectorClass* partial_lh_right = (VectorClass*)(right->partial_lh      + blockOffset);
            VectorClass        lh_max           = 0.0;
            UBYTE*             scale_dad;
            const UBYTE*       scale_left;
            const UBYTE*       scale_right;

            if (SAFE_NUMERIC) {
                size_t addr = ptn*ncat_mix;
                scale_dad   = dad_branch->scale_num + addr;
                scale_left  = left->scale_num       + addr;
                scale_right = right->scale_num      + addr;
            } else {
                scale_dad   = dad_branch->scale_num + ptn;
                scale_left  = left->scale_num       + ptn;
                scale_right = right->scale_num      + ptn;
                for (size_t i = 0; i < VectorClass::size(); i++) {
                    scale_dad[i] = scale_left[i] + scale_right[i];
                }
            }

            for (size_t c = 0; c < ncat_mix; c++) {
                if (SAFE_NUMERIC) {
                    lh_max   = 0.0;
                    size_t y = 0;
                    for (size_t x = 0; x < VectorClass::size()
                         ; x++, y += ncat_mix) {
                        scale_dad[y] = scale_left[y] + scale_right[y];
                    }
                }

                if (SITE_MODEL) {
                    // site-specific model
                    VectorClass* expleft  = partial_lh_tmp + nstates;
                    VectorClass* expright = expleft + nstates;
                    VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                    for (size_t i = 0; i < nstates; i++) {
                        expleft[i]  = exp(eval_ptr[i]*len_left[c])  * partial_lh_left[i];
                        expright[i] = exp(eval_ptr[i]*len_right[c]) * partial_lh_right[i];
                    }
                    
                    VectorClass* evec_ptr     = (VectorClass*) &evec[ptn*states_square];
                    VectorClass* inv_evec_ptr = (VectorClass*) &inv_evec[ptn*states_square];
                    for (size_t x = 0; x < nstates; x++) {
                        VectorClass *this_evec = evec_ptr + x*nstates;
#ifdef KERNEL_FIX_STATES
                        dotProductDualVec<VectorClass, VectorClass, nstates, FMA>(this_evec, expleft, this_evec, expright, partial_lh_tmp[x]);
#else
                        dotProductDualVec<VectorClass, VectorClass, FMA>(this_evec, expleft, this_evec, expright, partial_lh_tmp[x], nstates);
#endif
                    }
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, VectorClass, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max);
#else
                    productVecMat<VectorClass, VectorClass, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max, nstates);
#endif
                } else {
                    // normal model
                    double* eleft_ptr    = eleft;
                    double* eright_ptr   = eright;
                    // compute real partial likelihood vector
                    for (size_t x = 0; x < nstates; x++) {
#ifdef KERNEL_FIX_STATES
                        dotProductDualVec<VectorClass, double, nstates, FMA>(eleft_ptr, partial_lh_left, eright_ptr, partial_lh_right, partial_lh_tmp[x]);
#else
                        dotProductDualVec<VectorClass, double, FMA>(eleft_ptr, partial_lh_left, eright_ptr, partial_lh_right, partial_lh_tmp[x], nstates);
#endif
                        eleft_ptr += nstates;
                        eright_ptr += nstates;
                    }
                    
                    // compute dot-product with inv_eigenvector
                    double* inv_evec_ptr = inv_evec + mix_addr[c];
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, double, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max);
#else
                    productVecMat<VectorClass, double, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max, nstates);
#endif
                }
                // check if one should scale partial likelihoods
                if (SAFE_NUMERIC) {
                    auto underflown = ((lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                    if (horizontal_or(underflown)) {
                        for (size_t x = 0; x < VectorClass::size(); x++) {
                            if (underflown[x]) {
                                // BQM 2016-05-03: only scale for non-constant sites
                                // now do the likelihood scaling
                                double *partial_lh = dad_branch->partial_lh
                                    + (ptn*block + c*nstates*VectorClass::size() + x);
                                for (size_t i = 0, j=0; i < nstates; i++, j+=VectorClass::size()) {
                                    partial_lh[j] = ldexp(partial_lh[j], SCALING_THRESHOLD_EXP);
                                }
                                scale_dad[x*ncat_mix] += 1;
                            }
                        }
                    }
                    ++scale_dad;
                    ++scale_left;
                    ++scale_right;
                }
                partial_lh_left  += nstates;
                partial_lh_right += nstates;
                partial_lh       += nstates;
            }

            if (!SAFE_NUMERIC) {
                // check if one should scale partial likelihoods
                auto underflown = (lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                if (horizontal_or(underflown)) { // at least one site has numerical underflown
                    for (size_t x = 0; x < VectorClass::size(); x++) {
                        if (underflown[x]) {
                            double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                            // now do the likelihood scaling
                            for (size_t i = 0; i < block; i++) {
                                partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                            }
                            dad_branch->scale_num[ptn+x] += 1;
                        }
                    }
                }
            }
        } // big for loop over ptn
    }
    if (Params::getInstance().buffer_mem_save) {
        aligned_free(partial_lh_leaves);
        aligned_free(echildren);
        info.partial_lh_leaves = nullptr;
        info.echildren         = nullptr;
    }
}

/*******************************************************
 *
 * NEW! highly-vectorized log-likelihood derivative function
 *
 ******************************************************/


#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodBufferSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad
                                            , size_t ptn_lower, size_t ptn_upper, int packet_id
                                            , LikelihoodBufferSet& buffers)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodBufferGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad
                                                   , size_t ptn_lower, size_t ptn_upper, int packet_id
                                                   , LikelihoodBufferSet& buffers)
#endif
{
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = roundUpToMultiple(orig_nptn,VectorClass::size());
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t block = ncat_mix * nstates;
    size_t tip_block = nstates * model->getNMixtures();
    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (size_t c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    // reserve 3*block for computeLikelihoodDerv
    double *buffer_partial_lh_ptr = buffers.buffer_partial_lh + 3*get_safe_upper_limit(block);
    if (isMixlen()) {
        size_t nmix = getMixlen();
        buffer_partial_lh_ptr += nmix*(nmix+1)*VectorClass::size() + (nmix+3)*nmix*VectorClass::size()*num_packets;
    }

    // first compute partial_lh
    for (auto it = traversal_info.begin(); it != traversal_info.end(); it++) {
        computePartialLikelihood(*it, ptn_lower, ptn_upper, packet_id, buffers);
    }

    if (dad->isLeaf()) {
        // special treatment for TIP-INTERNAL NODE case
        double *tip_partial_lh_node = &tip_partial_lh[dad->id * max_orig_nptn * nstates];
        double *vec_tip = buffer_partial_lh_ptr + tip_block * VectorClass::size() * packet_id;
        auto stateRow = this->getConvertedSequenceByNumber(dad->id);
        auto unknown  = aln->STATE_UNKNOWN;

        size_t offset     = ptn_lower*block;
        size_t offsetStep = block*VectorClass::size();
        for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size(), offset+=offsetStep) {
            VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + offset);
            VectorClass *theta = (VectorClass*)(buffers.theta_all + offset);
            //load tip vector
            if (!SITE_MODEL) {
                for (size_t i = 0; i < VectorClass::size(); i++) {
                    int state;
                    if (ptn+i < orig_nptn) {
                        if (stateRow!=nullptr) {
                            state =  stateRow[ptn+i];
                        } else {
                            state = (aln->at(ptn+i))[dad->id];
                        }
                    } else if (ptn+i < max_orig_nptn) {
                        state = unknown;
                    } else if (ptn+i < nptn) {
                        state = model_factory->unobserved_ptns[ptn+i-max_orig_nptn][dad->id];
                    } else {
                        state = unknown;
                    }
                    double *this_tip_partial_lh = tip_partial_lh + tip_block*state;
                    double *this_vec_tip = vec_tip+i;
                    for (size_t c = 0; c < tip_block; c++) {
                        *this_vec_tip = this_tip_partial_lh[c];
                        this_vec_tip += VectorClass::size();
                    }
                }
            }
            VectorClass *lh_tip;
            if (SITE_MODEL) {
                lh_tip = (VectorClass*)&tip_partial_lh_node[ptn*nstates];
            }
            for (size_t c = 0; c < ncat_mix; c++) {
                if (!SITE_MODEL) {
                    lh_tip = (VectorClass*)(vec_tip + mix_addr_nstates[c]*VectorClass::size());
                }
                for (size_t i = 0; i < nstates; i++) {
                    theta[i] = lh_tip[i] * partial_lh_dad[i];
                }
                partial_lh_dad += nstates;
                theta += nstates;
            }
            if (SAFE_NUMERIC) {
                // numerical scaling per category
                UBYTE *scale_dad;
                UBYTE min_scale;
                for (size_t i = 0; i < VectorClass::size(); i++) {
                    scale_dad = dad_branch->scale_num+(ptn+i)*ncat_mix;
                    min_scale = scale_dad[0];
                    for (size_t c = 1; c < ncat_mix; c++) {
                        min_scale = min(min_scale, scale_dad[c]);
                    }
                    buffers.buffer_scale_all[ptn+i] = min_scale;

                    for (size_t c = 0; c < ncat_mix; c++) {
                        if (scale_dad[c] == min_scale+1) {
                            double *this_theta = &buffers.theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] *= SCALING_THRESHOLD;
                            }
                        } else if (scale_dad[c] > min_scale+1) {
                            double *this_theta = &buffers.theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] = 0.0;
                            }
                        }
                    }
                }
            } else {
                // normal scaling
                for (size_t i = 0; i < VectorClass::size(); i++) {
                    buffers.buffer_scale_all[ptn+i] = dad_branch->scale_num[ptn+i];
                }
            }
            VectorClass *buf = (VectorClass*)(buffers.buffer_scale_all+ptn);
            *buf *= LOG_SCALING_THRESHOLD;

        } // FOR PTN LOOP
//            aligned_free(vec_tip);
    } else {
        //------- both dad and node are internal nodes  --------//
        // now compute theta
        for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            VectorClass* theta = (VectorClass*)(buffers.theta_all + ptn*block);
            VectorClass* partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
            VectorClass* partial_lh_dad  = (VectorClass*)(dad_branch->partial_lh  + ptn*block);
            for (size_t i = 0; i < block; i++) {
                theta[i] = partial_lh_node[i] * partial_lh_dad[i];
            }
            if (SAFE_NUMERIC) {
                // numerical scaling per category
                UBYTE  min_scale;
                UBYTE  sum_scale[ncat_mix];
                size_t ptn_ncat   = ptn*ncat_mix;
                const UBYTE* scale_dad  = dad_branch->scale_num  + ptn_ncat;
                const UBYTE* scale_node = node_branch->scale_num + ptn_ncat;

                for (size_t i = 0; i < VectorClass::size(); i++) {
                    min_scale = sum_scale[0] = scale_dad[0] + scale_node[0];
                    for (size_t c = 1; c < ncat_mix; c++) {
                        sum_scale[c] = scale_dad[c] + scale_node[c];
                        min_scale    = min(min_scale, sum_scale[c]);
                    }
                    buffers.buffer_scale_all[ptn+i] = min_scale;

                    for (size_t c = 0; c < ncat_mix; c++) {
                        if (sum_scale[c] == min_scale+1) {
                            double *this_theta = &buffers.theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] *= SCALING_THRESHOLD;
                            }
                        } else if (sum_scale[c] > min_scale+1) {
                            double *this_theta = &buffers.theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] = 0.0;
                            }
                        }
                    }
                    scale_dad  += ncat_mix;
                    scale_node += ncat_mix;
                }
            } else {
                for (size_t i = 0; i < VectorClass::size(); i++) {
                    buffers.buffer_scale_all[ptn+i] = dad_branch->scale_num[ptn+i] + node_branch->scale_num[ptn+i];
                }
            }
            VectorClass *buf = (VectorClass*)(buffers.buffer_scale_all+ptn);
            *buf *= LOG_SCALING_THRESHOLD;
        } // FOR ptn
    } // internal node
}

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC,
          const int nstates, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervSIMD
    ( PhyloNeighbor* dad_branch, PhyloNode* dad,
      double* df,                double* ddf,
      LikelihoodBufferSet& buffers )
#else
template <class VectorClass, const bool SAFE_NUMERIC,
          const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervGenericSIMD
    ( PhyloNeighbor* dad_branch, PhyloNode* dad,
      double* df,                double* ddf,
      LikelihoodBufferSet& buffers)
#endif
{
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    if (!central_partial_lh) {
        initializeAllPartialLh();
    }
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, buffers, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, buffers, false);
#endif

#ifndef KERNEL_FIX_STATES
    size_t nstates   = aln->num_states;
#endif
    size_t ncat       = site_rate->getNRate();
    size_t ncat_mix   = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block      = ncat_mix * nstates;
//    size_t tip_block = nstates * model->getNMixtures();
    size_t orig_nptn  = aln->size();
    size_t max_orig_nptn = roundUpToMultiple(orig_nptn, VectorClass::size());
    size_t nptn       = max_orig_nptn+model_factory->unobserved_ptns.size();
    ASCType ASC_type  = model_factory->getASC();
    bool ASC_Holder   = (ASC_type == ASC_VARIANT_MISSING || ASC_type == ASC_INFORMATIVE_MISSING);
    bool ASC_Lewis    = (ASC_type == ASC_VARIANT || ASC_type == ASC_INFORMATIVE);

    double *const_df  = nullptr;
    double *const_ddf = nullptr;

    if (ASC_Holder) {
        const_df  = aligned_alloc<double>(get_safe_upper_limit(nptn) - max_orig_nptn);
        const_ddf = aligned_alloc<double>(get_safe_upper_limit(nptn) - max_orig_nptn);
    }
    
    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix], cat_id[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (size_t c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        cat_id[c] = c%ncat;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    double *eval = model->getEigenvalues();
    ASSERT(eval);

    double *buffer_partial_lh_ptr = buffers.buffer_partial_lh;
    vector<size_t> limits;
    computePatternPacketBounds(VectorClass::size(), num_threads,
                               num_packets, nptn, limits);

    ASSERT(buffers.theta_all);

    double *val0 = NULL;
    double *val1 = NULL;
    double *val2 = NULL;
    double cat_rate[ncat];
    double cat_prop[ncat];

    if (SITE_MODEL) {
        for (size_t c = 0; c < ncat; c++) {
            cat_rate[c] = site_rate->getRate(c);
            cat_prop[c] = site_rate->getProp(c);
        }
    } else {
        val0 = buffer_partial_lh_ptr;
        val1 = val0 + get_safe_upper_limit(block);
        val2 = val1 + get_safe_upper_limit(block);
        buffer_partial_lh_ptr += 3*get_safe_upper_limit(block);
        if (nstates % VectorClass::size() == 0) {
            VectorClass *vc_val0 = (VectorClass*)val0;
            VectorClass *vc_val1 = (VectorClass*)val1;
            VectorClass *vc_val2 = (VectorClass*)val2;

            size_t loop_size = nstates/VectorClass::size();
            for (size_t c = 0; c < ncat_mix; c++) {
                size_t       m        = c/denom;
                size_t       mycat    = c%ncat;
                double       len      = dad_branch->getLength(mycat);
                VectorClass* eval_ptr = (VectorClass*)(eval + mix_addr_nstates[c]);
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                double myrate = site_rate->getRate(mycat);
                for (size_t i = 0; i < loop_size; i++) {
                    VectorClass cof   = eval_ptr[i] * myrate;
                    VectorClass val   = exp(cof*len) * prop;
                    VectorClass val1_ = cof*val;
                    vc_val0[i] = val;
                    vc_val1[i] = val1_;
                    vc_val2[i] = cof*val1_;
                }
                vc_val0 += loop_size;
                vc_val1 += loop_size;
                vc_val2 += loop_size;
            }
        } else {
            for (size_t c = 0; c < ncat_mix; c++) {
                size_t  m        = c/denom;
                double* eval_ptr = eval + mix_addr_nstates[c];
                size_t  mycat    = c%ncat;
                double  prop     = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                size_t  addr     = c*nstates;
                double  len      = dad_branch->getLength(mycat);
                for (size_t i = 0; i < nstates; i++) {
                    double cof   = eval_ptr[i]*site_rate->getRate(mycat);
                    double val   = exp(cof*len) * prop;
                    double val1_ = cof*val;
                    val0[addr+i] = val;
                    val1[addr+i] = val1_;
                    val2[addr+i] = cof*val1_;
                }
            }
        }
    }

    double       dad_length = dad_branch->length;
    VectorClass* all_dfvec  = nullptr;
    VectorClass* all_ddfvec = nullptr;

    size_t nmixlen = getMixlen(), nmixlen2 = nmixlen*nmixlen;
    if (isMixlen()) {
        ASSERT(nmixlen == ncat);
        all_dfvec  = (VectorClass*)buffer_partial_lh_ptr;
        all_ddfvec = all_dfvec + nmixlen;
        buffer_partial_lh_ptr += nmixlen*(nmixlen+1)*VectorClass::size();
        for (size_t i = 0; i < nmixlen; i++) all_dfvec[i] = 0.0;
        for (size_t i = 0; i < nmixlen2; i++) all_ddfvec[i] = 0.0;
    }
    
    double all_lh(0.0), all_df(0.0), all_ddf(0.0);
    double all_prob_const(0.0), all_df_const(0.0), all_ddf_const(0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(num_threads) reduction(+:all_lh,all_df,all_ddf,all_prob_const,all_df_const,all_ddf_const)
#endif
    for (int packet_id = 0; packet_id < num_packets; packet_id++) {
        VectorClass my_df(0.0), my_ddf(0.0), vc_prob_const(0.0);
        VectorClass vc_df_const(0.0), vc_ddf_const(0.0);
        size_t ptn_lower = limits[packet_id];
        size_t ptn_upper = limits[packet_id+1];

        if (!buffers.theta_computed) {
            #ifdef KERNEL_FIX_STATES
                computeLikelihoodBufferSIMD<VectorClass, SAFE_NUMERIC, nstates, FMA, SITE_MODEL>
                    ( dad_branch, dad, ptn_lower, ptn_upper, packet_id, buffers );
            #else
                computeLikelihoodBufferGenericSIMD<VectorClass, SAFE_NUMERIC, FMA, SITE_MODEL>
                    ( dad_branch, dad, ptn_lower, ptn_upper, packet_id, buffers );
            #endif
        }

        if (isMixlen()) {
            // mixed branch length model
            VectorClass  lh_ptn;
            VectorClass* df_ptn  = ((VectorClass*)buffer_partial_lh_ptr) + (nmixlen+3)*nmixlen * packet_id; //PROBLEM
            VectorClass* ddf_ptn = df_ptn+nmixlen;

            VectorClass  my_lh(0.0);
            VectorClass* my_df  = df_ptn + nmixlen*2;
            VectorClass* my_ddf = df_ptn + nmixlen*3;
            for (size_t i = 0; i < nmixlen; i++) my_df[i] = 0.0;
            for (size_t i = 0; i < nmixlen2; i++) my_ddf[i] = 0.0;

            for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                for (size_t i = 0; i < nmixlen; i++) {
                    df_ptn[i] = ddf_ptn[i] = 0.0;
                }
                lh_ptn = 0.0;
                const VectorClass* theta    = (VectorClass*)(buffers.theta_all + ptn*block);
                const double* val0_ptr = val0;
                const double* val1_ptr = val1;
                const double* val2_ptr = val2;
                for (size_t c = 0; c < ncat_mix; c++) {
                    size_t i = cat_id[c];
                #ifdef KERNEL_FIX_STATES
                    dotProductTriple<VectorClass, double, nstates, FMA, true>(val0_ptr, val1_ptr, val2_ptr, theta, lh_ptn, df_ptn[i], ddf_ptn[i], nstates);
                #else
                    dotProductTriple<VectorClass, double, FMA, true>(val0_ptr, val1_ptr, val2_ptr, theta, lh_ptn, df_ptn[i], ddf_ptn[i],nstates, nstates);
                #endif
                    val0_ptr += nstates;
                    val1_ptr += nstates;
                    val2_ptr += nstates;
                    theta    += nstates;
                }
                lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
                if (ptn < orig_nptn) {
                    VectorClass freq;
                    freq.load_a(&ptn_freq[ptn]);
                    VectorClass inv_lh_ptn = 1.0 / lh_ptn;

                    // compute gradient (my_df)
                    for (size_t i = 0; i < nmixlen; i++) {
                        df_ptn[i]  *= inv_lh_ptn;
                        ddf_ptn[i] *= inv_lh_ptn;
                        my_df[i]    = mul_add(df_ptn[i], freq, my_df[i]);
                    }
                    // now compute hessian matrix my_ddf
                    for (size_t i = 0; i < nmixlen; i++) {
                        my_ddf[i*nmixlen+i] += nmul_add(df_ptn[i],df_ptn[i], ddf_ptn[i]) * freq;
                        for (size_t c = 0; c < nmixlen; c++) {
                            if (c!=i) {
                                my_ddf[i*nmixlen+c] -= df_ptn[i]*df_ptn[c]*freq;
                            }
                        }
                    }
                    lh_ptn = log(lh_ptn) + VectorClass().load_a(&buffers.buffer_scale_all[ptn]);
                    my_lh  = mul_add(lh_ptn, freq, my_lh);
                } else {
                    ASSERT(0 && "TODO +ASC not supported");
                }
            } // FOR ptn
            all_lh += horizontal_add(my_lh); //handled by reduction clause
        #ifdef _OPENMP
        #pragma omp critical
        #endif
            {
                for (size_t i = 0; i < nmixlen; i++) {
                    all_dfvec[i] += my_df[i];
                }
                for (size_t i = 0; i < nmixlen2; i++) {
                    all_ddfvec[i] += my_ddf[i];
                }
            }
        } else {
            // normal joint branch length model
            for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn;
                //lh_ptn.load_a(&ptn_invar[ptn]);
                VectorClass *theta = (VectorClass*)(buffers.theta_all + ptn*block);
                VectorClass df_ptn, ddf_ptn;

                if (SITE_MODEL) {
                    VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                    lh_ptn = 0.0; df_ptn = 0.0; ddf_ptn = 0.0;
                    for (size_t c = 0; c < ncat; c++) {
                        VectorClass lh_cat(0.0), df_cat(0.0), ddf_cat(0.0);
                        for (size_t i = 0; i < nstates; i++) {
                            VectorClass cof  = eval_ptr[i] * cat_rate[c];
                            VectorClass val  = exp(cof*dad_length)*theta[i];
                            VectorClass val1 = cof*val;
                            lh_cat  += val;
                            df_cat  += val1;
                            ddf_cat  = mul_add(cof, val1, ddf_cat);
                        }
                        lh_ptn   = mul_add(cat_prop[c], lh_cat, lh_ptn);
                        df_ptn   = mul_add(cat_prop[c], df_cat, df_ptn);
                        ddf_ptn  = mul_add(cat_prop[c], ddf_cat, ddf_ptn);
                        theta   += nstates;
                    }
                } else {
            #ifdef KERNEL_FIX_STATES
                    dotProductTriple<VectorClass, double, nstates, FMA, false>(val0, val1, val2, theta, lh_ptn, df_ptn, ddf_ptn, block);
            #else
                    dotProductTriple<VectorClass, double, FMA, false>(val0, val1, val2, theta, lh_ptn, df_ptn, ddf_ptn, block, nstates);
            #endif
                }
                lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);

                if (ptn < orig_nptn) {
                    lh_ptn = 1.0 / lh_ptn;
                    VectorClass df_frac  = df_ptn   * lh_ptn;
                    VectorClass ddf_frac = ddf_ptn  * lh_ptn;
                    VectorClass freq;
                    freq.load_a(&ptn_freq[ptn]);
                    VectorClass tmp1     = df_frac  * freq;
                    VectorClass tmp2     = ddf_frac * freq;
                    my_df  += tmp1;
                    my_ddf += nmul_add(tmp1, df_frac, tmp2);
                } else {
                    // ascertainment bias correction
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                        df_ptn.cutoff(nptn-ptn);
                        ddf_ptn.cutoff(nptn-ptn);
                    }
                    if (horizontal_or(VectorClass().load_a(&buffers.buffer_scale_all[ptn]) != 0.0)) {
                        // some entries are rescaled
                        double* lh_ptn_dbl  = (double*)&lh_ptn;
                        double* df_ptn_dbl  = (double*)&df_ptn;
                        double* ddf_ptn_dbl = (double*)&ddf_ptn;
                        for (size_t i = 0; i < VectorClass::size(); i++) {
                            if (buffers.buffer_scale_all[ptn+i] != 0.0) {
                                lh_ptn_dbl[i]  *= SCALING_THRESHOLD;
                                df_ptn_dbl[i]  *= SCALING_THRESHOLD;
                                ddf_ptn_dbl[i] *= SCALING_THRESHOLD;
                            }
                        }
                    }
                    if (ASC_Holder) {
                        lh_ptn.store_a(&buffers._pattern_lh[ptn]);
                        df_ptn.store_a(&const_df[ptn-max_orig_nptn]);
                        ddf_ptn.store_a(&const_ddf[ptn-max_orig_nptn]);
                    } else {
                        vc_prob_const += lh_ptn;
                        vc_df_const   += df_ptn;
                        vc_ddf_const  += ddf_ptn;
                    }
                }
            } // FOR ptn
            {
                //These adds don't need to be in a critical section, as the
                //shared variables are all mentioned in the reduction clause
                //Note that, since those variables are now double, the
                //horizontal_add operations need to be done here.
                all_df  += horizontal_add(my_df);
                all_ddf += horizontal_add(my_ddf);
                if (ASC_Lewis) {
                    all_prob_const += horizontal_add(vc_prob_const);
                    all_df_const   += horizontal_add(vc_df_const);
                    all_ddf_const  += horizontal_add(vc_ddf_const);
                }
            }
        } // else isMixlen()
    } // FOR packet

    // mark buffer as computed
    buffers.theta_computed = true;

    if (isMixlen()) {
        // mixed branch length model
        for (size_t i = 0; i < nmixlen; i++) {
            df[i] = horizontal_add(all_dfvec[i]);
            ASSERT(std::isfinite(df[i]) && "Numerical underflow for lh-derivative");
        }
        for (size_t i = 0; i < nmixlen2; i++) {
            ddf[i] = horizontal_add(all_ddfvec[i]);
        }
        // NOTE: last entry of df now store log-likelihood!
        df[nmixlen] = all_lh;
        return;
    }

    // normal joint branch length model
    *df  = all_df;
    *ddf = all_ddf;
    if (!std::isfinite(*df)) {
        if (!SAFE_NUMERIC || !warnedAboutNumericalUnderflow) {
            warnedAboutNumericalUnderflow = true;
            hideProgress();
            getModel()->writeInfo(cout);
            getRate()->writeInfo(cout);
            showProgress();
        }
        if (!SAFE_NUMERIC) {
            outError("Numerical underflow (lh-derivative)."
                     " Run again with the safe likelihood kernel via `-safe` option");
        }
    }
    if (ASC_Holder) {
        // Mark Holder's ascertainment bias correction for missing data
        double *const_lh = buffers._pattern_lh + max_orig_nptn;
        size_t step_unobserved_ptns = model_factory->unobserved_ptns.size() / nstates;
        double *const_lh_next  = const_lh + step_unobserved_ptns;
        double *const_df_next  = const_df + step_unobserved_ptns;
        double *const_ddf_next = const_ddf + step_unobserved_ptns;
        for (int step = 1; step < nstates; step++) {
            for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size()) {
                (VectorClass().load_a(&const_lh[ptn]) + VectorClass().load_a(&const_lh_next[ptn])).store_a(&const_lh[ptn]);
                (VectorClass().load_a(&const_df[ptn]) + VectorClass().load_a(&const_df_next[ptn])).store_a(&const_df[ptn]);
                (VectorClass().load_a(&const_ddf[ptn]) + VectorClass().load_a(&const_ddf_next[ptn])).store_a(&const_ddf[ptn]);
            }
            const_lh_next  += step_unobserved_ptns;
            const_df_next  += step_unobserved_ptns;
            const_ddf_next += step_unobserved_ptns;
        }
        // cutoff the last entries if going beyond
        for (size_t ptn = orig_nptn; ptn < max_orig_nptn; ptn++) {
            const_lh[ptn]  = 0.0;
            const_df[ptn]  = 0.0;
            const_ddf[ptn] = 0.0;
        }
        VectorClass sum_df = 0.0, sum_ddf = 0.0;
        
        for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size()) {
            VectorClass prob_const, df_const, ddf_const;
            prob_const.load_a(&const_lh[ptn]);
            df_const.load_a(&const_df[ptn]);
            ddf_const.load_a(&const_ddf[ptn]);
            prob_const = 1.0 - prob_const;
            VectorClass df_frac = df_const / prob_const;
            VectorClass ddf_frac = ddf_const / prob_const;
            sum_df  += VectorClass().load_a(&ptn_freq[ptn]) * df_frac;
            sum_ddf += VectorClass().load_a(&ptn_freq[ptn]) * (ddf_frac + df_frac*df_frac);
        }
        *df  += horizontal_add(sum_df);
        *ddf += horizontal_add(sum_ddf);
        aligned_free(const_ddf);
        aligned_free(const_df);
    } else if (ASC_Lewis) {
        // ascertainment bias correction
        all_prob_const = 1.0 - all_prob_const;
        double df_frac = all_df_const / all_prob_const;
        double ddf_frac = all_ddf_const / all_prob_const;
        size_t nsites = aln->getNSite();
        *df  += nsites * df_frac;
        *ddf += nsites *(ddf_frac + df_frac*df_frac);
    }
    if (!std::isfinite(*df)) {
        if (!warnedAboutNumericalUnderflow) {
            hideProgress();
            cout << "WARNING: Numerical underflow for lh-derivative" << endl;
            showProgress();
            warnedAboutNumericalUnderflow = true;
        }
        *df = *ddf = 0.0;
    }
}

/*******************************************************
 *
 * NEW! highly-vectorized log-likelihood function
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA, const bool SITE_MODEL>
double PhyloTree::computeLikelihoodBranchSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                              LikelihoodBufferSet& buffers)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
double PhyloTree::computeLikelihoodBranchGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                     LikelihoodBufferSet& buffers)
#endif
{
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    
    if (!central_partial_lh) {
        initializeAllPartialLh();
    }
    if (node->isLeaf() ) {
        std::swap(dad,        node);
        std::swap(dad_branch, node_branch);
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, buffers, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, buffers, false);
#endif
    double tree_lh       = 0.0;
#ifndef KERNEL_FIX_STATES
    size_t nstates       = aln->num_states;
#endif
    size_t ncat          = site_rate->getNRate();
    size_t ncat_mix      = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block         = ncat_mix * nstates;
    size_t tip_block     = nstates * model->getNMixtures();
    size_t orig_nptn     = aln->size();
    size_t max_orig_nptn = roundUpToMultiple(orig_nptn, VectorClass::size());
    size_t nptn          = max_orig_nptn + model_factory->unobserved_ptns.size();
    size_t tip_mem_size  = max_orig_nptn * nstates;
    ASCType ASC_type     = model_factory->getASC();
    bool   ASC_Holder    = (ASC_type == ASC_VARIANT_MISSING || ASC_type == ASC_INFORMATIVE_MISSING);
    bool   ASC_Lewis     = (ASC_type == ASC_VARIANT || ASC_type == ASC_INFORMATIVE);

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom         = (model_factory->fused_mix_rate) ? 1 : ncat;

    double *eval         = model->getEigenvalues();
    ASSERT(eval);

    double *val          = nullptr;
    double *buffer_partial_lh_ptr = buffers.buffer_partial_lh;

    double cat_length[ncat];
    double cat_prop[ncat];
    if (SITE_MODEL) {
        for (size_t c = 0; c < ncat; c++) {
            cat_length[c] = site_rate->getRate(c) * dad_branch->length;
            cat_prop[c] = site_rate->getProp(c);
        }
    } else {
        val = buffer_partial_lh_ptr;
        buffer_partial_lh_ptr += get_safe_upper_limit(block);
        if (nstates % VectorClass::size() == 0) {
            size_t loop_size = nstates / VectorClass::size();
            for (size_t c = 0; c < ncat_mix; c++) {
                size_t mycat          = c%ncat;
                size_t m              = c/denom;
                mix_addr_nstates[c]   = m*nstates;
                mix_addr[c]           = mix_addr_nstates[c]*nstates;
                VectorClass* eval_ptr = (VectorClass*)(eval + mix_addr_nstates[c]);
                double len            = site_rate->getRate(mycat)*dad_branch->getLength(mycat);
                double prop           = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                VectorClass *this_val = (VectorClass*)(val + c*nstates);
                for (size_t i = 0; i < loop_size; i++) {
                    this_val[i] = exp(eval_ptr[i]*len) * prop;
                }
            }
        } else {
            for (size_t c = 0; c < ncat_mix; c++) {
                size_t  mycat       = c%ncat;
                size_t  m           = c/denom;
                mix_addr_nstates[c] = m*nstates;
                mix_addr[c]         = mix_addr_nstates[c]*nstates;
                double* eval_ptr    = eval + mix_addr_nstates[c];
                double  len         = site_rate->getRate(mycat)*dad_branch->getLength(mycat);
                double  prop        = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                double* this_val    = val + c*nstates;
                for (size_t i = 0; i < nstates; i++)
                    this_val[i] = exp(eval_ptr[i]*len) * prop;
            }
        }
    }

    double all_tree_lh(0.0);
    double all_prob_const(0.0);

    vector<size_t> limits;
    computePatternPacketBounds(VectorClass::size(), num_threads,
                               num_packets, nptn, limits);
    
    if (dad->isLeaf()  ) {
        // special treatment for TIP-INTERNAL NODE case
        double* partial_lh_node;
        if (SITE_MODEL) {
            partial_lh_node = &tip_partial_lh[dad->id * tip_mem_size];
        }
        else {
            partial_lh_node        = buffer_partial_lh_ptr;
            buffer_partial_lh_ptr += get_safe_upper_limit((aln->STATE_UNKNOWN+1)*block);
//            IntVector states_dad = model->seq_states[dad->id];
//            states_dad.push_back(aln->STATE_UNKNOWN);
            // precompute information from one tip
            if (nstates % VectorClass::size() == 0) {
                // vectorized version
                for (int state = 0; state <= aln->STATE_UNKNOWN; state++) {
                    double*       lh_node    = partial_lh_node + state * block;
                    const double* lh_tip     = tip_partial_lh  + state * tip_block;
                    const double* vc_val_tmp = val;
                    for (size_t c = 0; c < ncat_mix; c++) {
                        const double *this_lh_tip = lh_tip + mix_addr_nstates[c];
                        for (size_t i = 0; i < nstates; i+=VectorClass::size()) {
                            (VectorClass().load_a(&vc_val_tmp[i]) * VectorClass().load_a(&this_lh_tip[i])).store_a(&lh_node[i]);
                        }
                        lh_node += nstates;
                        vc_val_tmp += nstates;
                    }
                }
            } else {
                // non-vectorized version
                for (int state = 0; state <= aln->STATE_UNKNOWN; state++) {
                    double* lh_node = partial_lh_node +state*block;
                    const double* val_tmp = val;
                    const double* this_tip_partial_lh = tip_partial_lh + state*tip_block;
                    for (size_t c = 0; c < ncat_mix; c++) {
                        const double *lh_tip = this_tip_partial_lh + mix_addr_nstates[c];
                        for (size_t i = 0; i < nstates; i++) {
                              lh_node[i] = val_tmp[i] * lh_tip[i];
                        }
                        lh_node += nstates;
                        val_tmp += nstates;
                    }
                }
            }
        }
        
        auto stateRow = this->getConvertedSequenceByNumber(dad->id);
        auto unknown  = aln->STATE_UNKNOWN;
        // now do the real computation
#ifdef _OPENMP
#pragma omp parallel for  schedule(dynamic,1) num_threads(num_threads) reduction(+:all_tree_lh,all_prob_const)
#endif
        for (int packet_id = 0; packet_id < num_packets; packet_id++) {
            VectorClass vc_tree_lh(0.0);
            VectorClass vc_prob_const(0.0);
            size_t ptn_lower = limits[packet_id];
            size_t ptn_upper = limits[packet_id+1];

            // reset memory for _pattern_lh_cat
            memset(buffers._pattern_lh_cat + ptn_lower*ncat_mix, 0, sizeof(double)*(ptn_upper-ptn_lower)*ncat_mix);

            // first compute partial_lh
            for (auto it = traversal_info.begin(); it != traversal_info.end(); ++it) {
                computePartialLikelihood(*it, ptn_lower, ptn_upper, packet_id, buffers);
            }
            double *vec_tip = buffer_partial_lh_ptr + block*VectorClass::size() * packet_id;
            for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass  lh_ptn(0.0);
                VectorClass* lh_cat = (VectorClass*)(buffers._pattern_lh_cat + ptn*ncat_mix);
                const VectorClass* partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
                const VectorClass* lh_node        = SITE_MODEL ? (VectorClass*)&partial_lh_node[ptn*nstates] : (VectorClass*)vec_tip;

                if (SITE_MODEL) {
                    // site-specific model
                    const VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                    for (size_t c = 0; c < ncat; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProductExp<VectorClass, double, nstates, FMA>(eval_ptr, lh_node, partial_lh_dad, cat_length[c], lh_cat[c]);
    #else
                        dotProductExp<VectorClass, double, FMA>(eval_ptr, lh_node, partial_lh_dad, cat_length[c], lh_cat[c], nstates);
    #endif
                        if (SAFE_NUMERIC) {
                            lh_cat[c] *= cat_prop[c];
                        }
                        else {
                            lh_ptn += (lh_cat[c] *= cat_prop[c]);
                        }
                        partial_lh_dad += nstates;
                    }
                } else { // normal model
                    //load tip vector
                    for (size_t i = 0; i < VectorClass::size(); i++) {
                        int state;
                        if (ptn+i < orig_nptn) {
                            if (stateRow!=nullptr) {
                                state =  stateRow[ptn+i];
                            } else {
                                state = (aln->at(ptn+i))[dad->id];
                            }
                        } else if (ptn+i < max_orig_nptn) {
                            state = unknown;
                        } else if (ptn+i < nptn) {
                            state = model_factory->unobserved_ptns[ptn+i-max_orig_nptn][dad->id];
                        } else {
                            state = aln->STATE_UNKNOWN;
                        }
                        const double* lh_tip = partial_lh_node + block*state;
                        double* this_vec_tip = vec_tip+i;
                        for (size_t c = 0; c < block; c++) {
                            *this_vec_tip = lh_tip[c];
                            this_vec_tip += VectorClass::size();
                        }
                    }
                    // compute likelihood per category
                    for (size_t c = 0; c < ncat_mix; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProductVec<VectorClass, VectorClass, nstates, FMA>(lh_node, partial_lh_dad, lh_cat[c]);
    #else
                        dotProductVec<VectorClass, VectorClass, FMA>(lh_node, partial_lh_dad, lh_cat[c], nstates);
    #endif
                        if (!SAFE_NUMERIC) {
                            lh_ptn += lh_cat[c];
                        }
                        lh_node        += nstates;
                        partial_lh_dad += nstates;
                    }
                } // if SITE_MODEL

                // compute scaling factor per pattern
                VectorClass vc_min_scale(0.0);
                double* vc_min_scale_ptr = (double*)&vc_min_scale;
                if (SAFE_NUMERIC) {
                    // numerical scaling per category
                    const UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
                    UBYTE min_scale;
                    for (size_t i = 0; i < VectorClass::size(); i++) {
                        min_scale = scale_dad[0];
                        for (size_t c = 1; c < ncat_mix; c++) {
                            min_scale = min(min_scale, scale_dad[c]);
                        }
                        vc_min_scale_ptr[i] = min_scale;

                        double *this_lh_cat = &buffers._pattern_lh_cat[ptn*ncat_mix + i];
                        for (size_t c = 0; c < ncat_mix; c++) {
                            // rescale lh_cat if necessary
                            if (scale_dad[c] == min_scale+1) {
                                this_lh_cat[c*VectorClass::size()] *= SCALING_THRESHOLD;
                            } else if (scale_dad[c] > min_scale+1) {
                                this_lh_cat[c*VectorClass::size()] = 0.0;
                            }
                        }
                        scale_dad += ncat_mix;
                    }
                    // now take the sum of (rescaled) lh_cat
                    sumVec<VectorClass, true>(lh_cat, lh_ptn, ncat_mix);
                } else {
                    for (size_t i = 0; i < VectorClass::size(); i++) {
                        vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i];
                    }
                }
                vc_min_scale *= LOG_SCALING_THRESHOLD;

                // Sum later to avoid underflow of invariant sites
                lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
                if (ptn < orig_nptn) {
                    lh_ptn = log(lh_ptn) + vc_min_scale;
                    lh_ptn.store_a(&buffers._pattern_lh[ptn]);
                    vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
                } else {
                    // ascertainment bias correction
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                    }
                    // bugfix 2016-01-21, prob_const can be rescaled
                    if (horizontal_or(vc_min_scale != 0.0)) {
                        // some entries are rescaled
                        double *lh_ptn_dbl = (double*)&lh_ptn;
                        for (size_t i = 0; i < VectorClass::size(); i++) {
                            if (vc_min_scale_ptr[i] != 0.0) {
                                lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                            }
                        }
                    }
                    if (ASC_Holder) {
                        lh_ptn.store_a(&buffers._pattern_lh[ptn]);
                    }
                    else {
                        vc_prob_const += lh_ptn;
                    }
                }
            } // FOR PTN
            {
                //These additions are not in a critical section, because the
                //all_tree_lh and all_prob_const variables are mentioned in
                //the reduction clause of the omp parallel for.
                double here_lh = horizontal_add(vc_tree_lh);
                all_tree_lh   += here_lh;
                if (ASC_Lewis) {
                    all_prob_const += horizontal_add(vc_prob_const);
                }
            }
        } // FOR packet
    } else {
    	//-------- both dad and node are internal nodes -----------/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(num_threads) reduction(+:all_tree_lh,all_prob_const)
#endif
        for (int packet_id = 0; packet_id < num_packets; packet_id++) {
            size_t ptn_lower = limits[packet_id];
            size_t ptn_upper = limits[packet_id+1];
            size_t ptn_count = ptn_upper - ptn_lower;

            // reset memory for _pattern_lh_cat
            memset(buffers._pattern_lh_cat + ptn_lower*ncat_mix, 0, sizeof(double)*ptn_count*ncat_mix);

            // first compute partial_lh
            for (auto it = traversal_info.begin(); it != traversal_info.end(); it++) {
                computePartialLikelihood(*it, ptn_lower, ptn_upper, packet_id, buffers);
            }

            VectorClass vc_tree_lh(0.0);
            VectorClass vc_prob_const(0.0);
            for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn(0.0);
                VectorClass *lh_cat = (VectorClass*)(buffers._pattern_lh_cat + ptn*ncat_mix);
                const VectorClass *partial_lh_dad  = (VectorClass*)(dad_branch->partial_lh + ptn*block);
                const VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);

                // compute likelihood per category
                if (SITE_MODEL) {
                    VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                    for (size_t c = 0; c < ncat; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProductExp<VectorClass, double, nstates, FMA>(eval_ptr, partial_lh_node, partial_lh_dad, cat_length[c], lh_cat[c]);
    #else
                        dotProductExp<VectorClass, double, FMA>(eval_ptr, partial_lh_node, partial_lh_dad, cat_length[c], lh_cat[c], nstates);
    #endif
                        if (SAFE_NUMERIC) {
                            lh_cat[c] *= cat_prop[c];
                        } else {
                            lh_ptn += (lh_cat[c] *= cat_prop[c]);
                        }
                        partial_lh_node += nstates;
                        partial_lh_dad += nstates;
                    }
                } else {
                    double *val_tmp = val;
                    for (size_t c = 0; c < ncat_mix; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProduct3Vec<VectorClass, double, nstates, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c]);
    #else
                        dotProduct3Vec<VectorClass, double, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c], nstates);
    #endif
                        if (!SAFE_NUMERIC) {
                            lh_ptn += lh_cat[c];
                        }
                        partial_lh_node += nstates;
                        partial_lh_dad  += nstates;
                        val_tmp         += nstates;
                    }
                } // if SITE MODEL

                // compute the scaling factor per pattern
                VectorClass vc_min_scale(0.0);
                double* vc_min_scale_ptr = (double*)&vc_min_scale;
                if (SAFE_NUMERIC) {
                    const UBYTE *scale_dad  = dad_branch->scale_num  + ptn*ncat_mix;
                    const UBYTE *scale_node = node_branch->scale_num + ptn*ncat_mix;
                    UBYTE sum_scale[ncat_mix];
                    UBYTE min_scale;

                    for (size_t i = 0; i < VectorClass::size(); i++) {
                        min_scale = sum_scale[0] = scale_dad[0] + scale_node[0];
                        for (size_t c = 1; c < ncat_mix; c++) {
                            sum_scale[c] = scale_dad[c] + scale_node[c];
                            min_scale = min(min_scale, sum_scale[c]);
                        }
                        vc_min_scale_ptr[i] = min_scale;
                        double *this_lh_cat = &buffers._pattern_lh_cat[ptn*ncat_mix + i];
                        for (size_t c = 0; c < ncat_mix; c++) {
                            if (sum_scale[c] == min_scale+1) {
                                this_lh_cat[c*VectorClass::size()] *= SCALING_THRESHOLD;
                            } else if (sum_scale[c] > min_scale+1) {
                                // reset if category is scaled a lot
                                this_lh_cat[c*VectorClass::size()] = 0.0;
                            }
                        }
                        scale_dad += ncat_mix;
                        scale_node += ncat_mix;
                    }
                    sumVec<VectorClass, true>(lh_cat, lh_ptn, ncat_mix);
                } else {
                    for (size_t i = 0; i < VectorClass::size(); i++) {
                        vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i] + node_branch->scale_num[ptn+i];
                    }
                } // if SAFE_NUMERIC
                vc_min_scale *= LOG_SCALING_THRESHOLD;

                // Sum later to avoid underflow of invariant sites
                lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);

                if (ptn < orig_nptn) {
                    lh_ptn     = log(lh_ptn) + vc_min_scale;
                    lh_ptn.store_a(&buffers._pattern_lh[ptn]);
                    vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
                } else {
                    // ascertainment bias correction
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                    }
                    // bugfix 2016-01-21, prob_const can be rescaled
                    if (horizontal_or(vc_min_scale != 0.0)) {
                        // some entries are rescaled
                        double *lh_ptn_dbl = (double*)&lh_ptn;
                        for (size_t i = 0; i < VectorClass::size(); i++) {
                            if (vc_min_scale_ptr[i] != 0.0) {
                                lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                            }
                        }
                    }
                    if (ASC_Holder) {
                        lh_ptn.store_a(&buffers._pattern_lh[ptn]);
                    }
                    else {
                        vc_prob_const += lh_ptn;
                    }
                }
            } // FOR LOOP ptn
            {
                //These additions don't need to be in a critical section,
                //because of the reduction(+:all_tree_lh,all_prob_const)
                //clause in the omp parallel for.
                double here_lh = horizontal_add(vc_tree_lh);
                all_tree_lh   += here_lh;
                if (ASC_Lewis) {
                    all_prob_const += horizontal_add(vc_prob_const);
                }
            }
        } // FOR thread
    } // else

    bool justWarned = false;
    tree_lh += all_tree_lh;
    if (!std::isfinite(tree_lh)) {
        if (!warnedAboutNumericalUnderflow) {
            warnedAboutNumericalUnderflow = true;
            justWarned = true;
            hideProgress();
            if (SAFE_NUMERIC) {
                outWarning("Numerical underflow for lh-branch");
            } else {
                outWarning("Numerical underflow for lh-branch, use `-safe` option to avoid this warning");
            }
            showProgress();
        }
        // arbitrarily fix tree_lh if underflown for some sites
        tree_lh = 0.0;
        
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,1) num_threads(num_threads) reduction(+:tree_lh)
        #endif
        for (size_t ptn = 0; ptn < orig_nptn; ++ptn) {
            if (!std::isfinite(buffers._pattern_lh[ptn])) {
                buffers._pattern_lh[ptn] = LOG_SCALING_THRESHOLD*4; // log(2^(-1024))
            }
            tree_lh += buffers._pattern_lh[ptn] * ptn_freq[ptn];
        }
        //Don't we need to set _pattern_lh[...] something?
        if (justWarned) {
            LOG_LINE(VB_MED, "Fixed tree_lh was " << tree_lh);
        }
    }
    // robust phylogeny project, summing log-likelihood over the best sites
    if (params->robust_phy_keep < 1.0) {
        if (ASC_Holder || ASC_Lewis) {
            outError("+ASC not supported with robust phylo");
        }
        size_t sites      = getAlnNSite();
        size_t sites_drop = (int)((1.0-params->robust_phy_keep)*sites);
        size_t site       = 0;
        for (size_t ptn = 0; ptn < orig_nptn; ptn++) {
            for (size_t i = 0; i < ptn_freq[ptn]; i++, site++) {
                _site_lh[site] = buffers._pattern_lh[ptn];
            }
        }
        ASSERT(site == sites);
        nth_element(_site_lh, _site_lh + sites_drop, _site_lh + sites);
        tree_lh = 0.0;
        for (size_t i = sites_drop; i < sites; i++) {
            tree_lh += _site_lh[i];
        }
    } else if (params->robust_median) {
        if (ASC_Holder || ASC_Lewis) {
            outError("+ASC not supported with robust phylo");
        }
        size_t sites = getAlnNSite();
        size_t site = 0;
        for (size_t ptn = 0; ptn < orig_nptn; ptn++) {
            for (int i = 0; i < ptn_freq[ptn]; i++, site++) {
                _site_lh[site] = buffers._pattern_lh[ptn];
            }
        }
        ASSERT(site == sites);
        nth_element(_site_lh, _site_lh + sites/2, _site_lh + sites);
        if ((sites & 1) == 1) {
            // odd number of sites
            tree_lh = _site_lh[sites/2];
        } else {
            // even number of sites
            tree_lh = 0.5*(_site_lh[sites/2] + *max_element(_site_lh, _site_lh + sites/2));
        }
    }
    if (ASC_Holder) {
        // Mark Holder's ascertainment bias correction for missing data
        double *const_lh = buffers._pattern_lh + max_orig_nptn;
        size_t step_unobserved_ptns = model_factory->unobserved_ptns.size() / nstates;
        double *const_lh_next = const_lh + step_unobserved_ptns;
        for (int step = 1; step < nstates; step++, const_lh_next += step_unobserved_ptns) {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size())
                (VectorClass().load_a(&const_lh[ptn]) + VectorClass().load_a(&const_lh_next[ptn])).store_a(&const_lh[ptn]);
        }
        // cutoff the last entries if going beyond
        for (size_t ptn = orig_nptn; ptn < max_orig_nptn; ptn++) {
            const_lh[ptn] = 0.0;
        }
        VectorClass sum_corr = 0.0;
        for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size()) {
            VectorClass prob_variant = log(1.0 - VectorClass().load_a(&const_lh[ptn]));
            (VectorClass().load_a(&buffers._pattern_lh[ptn]) - prob_variant).store_a(&buffers._pattern_lh[ptn]);
            sum_corr += prob_variant*VectorClass().load_a(&ptn_freq[ptn]);
        }
        tree_lh -= horizontal_add(sum_corr);
    } else if (ASC_Lewis) {
    	// ascertainment bias correction
        if (all_prob_const >= 1.0 || all_prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        ASSERT(all_prob_const < 1.0 && all_prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
        //        double inv_const = 1.0 / (1.0-prob_const);
        //        size_t nptn_cat = orig_nptn*ncat;
        //    	for (ptn = 0; ptn < nptn_cat; ptn++)
        //            _pattern_lh_cat[ptn] *= inv_const;
        
        all_prob_const = log(1.0 - all_prob_const);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size()) {
            (VectorClass().load_a(&buffers._pattern_lh[ptn])-all_prob_const).store_a(&buffers._pattern_lh[ptn]);
        }
        tree_lh -= aln->getNSite()*all_prob_const;
        ASSERT(std::isfinite(tree_lh));
    }
    return tree_lh;
}


/*******************************************************
 *
 * NEW! highly-vectorized log-likelihood from buffer
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const int nstates, const bool FMA, const bool SITE_MODEL>
double PhyloTree::computeLikelihoodFromBufferSIMD(LikelihoodBufferSet& buffers)
#else
template <class VectorClass, const bool FMA, const bool SITE_MODEL>
double PhyloTree::computeLikelihoodFromBufferGenericSIMD(LikelihoodBufferSet& buffers)
#endif
{
    ASSERT(buffers.theta_all && buffers.theta_computed);

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = roundUpToMultiple(orig_nptn,VectorClass::size());
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    ASCType ASC_type = model_factory->getASC();
    bool ASC_Holder = (ASC_type == ASC_VARIANT_MISSING || ASC_type == ASC_INFORMATIVE_MISSING);
    bool ASC_Lewis = (ASC_type == ASC_VARIANT || ASC_type == ASC_INFORMATIVE);

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (size_t c = 0; c < ncat_mix; ++c) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    double *eval = model->getEigenvalues();
    ASSERT(eval);

    double *val0 = NULL;
    double cat_length[ncat];
    double cat_prop[ncat];

    if (SITE_MODEL) {
        for (size_t c = 0; c < ncat; ++c) {
            cat_length[c] = site_rate->getRate(c) * current_it->length;
            cat_prop[c] = site_rate->getProp(c);
        }
    } else {
        val0 = buffers.buffer_partial_lh;
        if (nstates % VectorClass::size() == 0) {
            VectorClass *vc_val0 = (VectorClass*)val0;
            size_t loop_size = nstates / VectorClass::size();
            for (size_t c = 0; c < ncat_mix; ++c) {
                size_t m = c/denom;
                VectorClass *eval_ptr = (VectorClass*)(eval + mix_addr_nstates[c]);
                size_t mycat = c%ncat;
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                double len = site_rate->getRate(mycat) * current_it->getLength(mycat);
                for (size_t i = 0; i < loop_size; ++i) {
                    vc_val0[i] = exp(eval_ptr[i] * len) * prop;
                }
                vc_val0 += loop_size;
            }
        } else {
            for (size_t c = 0; c < ncat_mix; ++c) {
                size_t m = c/denom;
                double *eval_ptr = eval + mix_addr_nstates[c];
                size_t mycat = c%ncat;
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                size_t addr = c*nstates;
                for (size_t i = 0; i < nstates; ++i) {
                    double cof = eval_ptr[i]*site_rate->getRate(mycat);
                    double val = exp(cof*current_it->getLength(mycat)) * prop;
                    val0[addr+i] = val;
                }
            }
        }
    }

    double all_tree_lh(0.0), all_prob_const(0.0);

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(num_threads) reduction(+:all_tree_lh,all_prob_const)
    #endif
    for (size_t ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
        VectorClass lh_ptn(0.0);
        VectorClass *theta = (VectorClass*)(buffers.theta_all + ptn*block);
        if (SITE_MODEL) {
            VectorClass *eval_ptr = (VectorClass*)&eval[ptn*nstates];
            for (size_t c = 0; c < ncat; c++) {
                VectorClass lh_cat;
#ifdef KERNEL_FIX_STATES
                dotProductExp<VectorClass, double, nstates, FMA>(eval_ptr, theta, cat_length[c], lh_cat);
#else
                dotProductExp<VectorClass, double, FMA>(eval_ptr, theta, cat_length[c], lh_cat, nstates);
#endif
                lh_ptn = mul_add(lh_cat, cat_prop[c], lh_ptn);
                theta += nstates;
            }
        } else {
            dotProductVec<VectorClass, double, FMA>(val0, theta, lh_ptn, block);
        }

        // Sum later to avoid underflow of invariant sites
        lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);

        VectorClass vc_tree_lh(0.0);
        VectorClass vc_prob_const(0.0);
        if (ptn < orig_nptn) {
            lh_ptn = log(abs(lh_ptn)) + VectorClass().load_a(&buffers.buffer_scale_all[ptn]);
            lh_ptn.store_a(&buffers._pattern_lh[ptn]);
            vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
        } else {
            if (ptn+VectorClass::size() > nptn) {
                // cutoff the last entries if going beyond
                lh_ptn.cutoff(nptn-ptn);
            }
            if (horizontal_or(VectorClass().load_a(&buffers.buffer_scale_all[ptn]) != 0.0)) {
                // some entries are rescaled
                double *lh_ptn_dbl = (double*)&lh_ptn;
                for (size_t i = 0; i < VectorClass::size(); i++) {
                    if (buffers.buffer_scale_all[ptn+i] != 0.0) {
                        lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                    }
                }
            }
            if (ASC_Holder) {
                lh_ptn.store_a(&buffers._pattern_lh[ptn]);
            }
            else {
                vc_prob_const += lh_ptn;
            }
        }
        {
            //Don't need a critical section here, because of the
            //reduction(+:all_tree_lh,all_prob_const) clause on
            //the omp parallel for.
            all_tree_lh += horizontal_add(vc_tree_lh);
            if (ASC_Lewis) {
                all_prob_const += horizontal_add(vc_prob_const);
            }
        }
    }

    double tree_lh = all_tree_lh;
    if (!std::isfinite(tree_lh)) {
        if (!safe_numeric) {
            outError("Numerical underflow (lh-from-buffer)."
                     " Run again with the safe likelihood kernel via `-safe` option");
        } else {
            if (!warnedAboutNumericalUnderflow) {
                hideProgress();
                outWarning("Numerical underflow (from lh-from-buffer).");
                showProgress();
                warnedAboutNumericalUnderflow = true;
            }
        }
    }

    // arbitrarily fix tree_lh if underflown for some sites
    if (!std::isfinite(tree_lh)) {
        tree_lh = 0.0;
        for (size_t ptn = 0; ptn < orig_nptn; ++ptn) {
            if (!std::isfinite(buffers._pattern_lh[ptn])) {
                buffers._pattern_lh[ptn] = LOG_SCALING_THRESHOLD*4; // log(2^(-1024))
            }
            tree_lh += buffers._pattern_lh[ptn] * ptn_freq[ptn];
        }
    }

    if (ASC_Holder) {
        // Mark Holder's ascertainment bias correction for missing data
        double *const_lh = buffers._pattern_lh + max_orig_nptn;
        size_t step_unobserved_ptns = model_factory->unobserved_ptns.size() / nstates;
        double *const_lh_next = const_lh + step_unobserved_ptns;
        for (int step = 1; step < nstates; step++, const_lh_next += step_unobserved_ptns) {
            for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size()) {
                (VectorClass().load_a(&const_lh[ptn]) + VectorClass().load_a(&const_lh_next[ptn])).store_a(&const_lh[ptn]);
            }
        }
        // cutoff the last entries if going beyond
        for (size_t ptn = orig_nptn; ptn < max_orig_nptn; ++ptn) {
            const_lh[ptn] = 0.0;
        }
        VectorClass sum_corr = 0.0;
        for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size()) {
            VectorClass prob_variant = log(1.0 - VectorClass().load_a(&const_lh[ptn]));
            (VectorClass().load_a(&buffers._pattern_lh[ptn]) - prob_variant).store_a(&buffers._pattern_lh[ptn]);
            sum_corr += prob_variant*VectorClass().load_a(&ptn_freq[ptn]);
        }
        tree_lh -= horizontal_add(sum_corr);
    } else if (ASC_Lewis) {
    	// ascertainment bias correction
        if (all_prob_const >= 1.0 || all_prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        ASSERT(all_prob_const < 1.0 && all_prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
        //        double inv_const = 1.0 / (1.0-prob_const);
        //        size_t nptn_cat = orig_nptn*ncat;
        //    	for (ptn = 0; ptn < nptn_cat; ptn++)
        //            _pattern_lh_cat[ptn] *= inv_const;
        
        all_prob_const = log(1.0 - all_prob_const);
        for (size_t ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size()) {
            (VectorClass().load_a(&buffers._pattern_lh[ptn])-all_prob_const).store_a(&buffers._pattern_lh[ptn]);
            //_pattern_lh[ptn] -= prob_const;
        }
        tree_lh -= aln->getNSite()*all_prob_const;
        ASSERT(std::isfinite(tree_lh));
    }
    return tree_lh;
}


#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervMixlenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                double &df, double &ddf,
                                                LikelihoodBufferSet& buffers)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervMixlenGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                       double &df, double &ddf, LikelihoodBufferSet& buffers)
#endif
{
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    if (!central_partial_lh) {
        initializeAllPartialLh();
    }
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, buffers, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, buffers, false);
    size_t nstates = aln->num_states;
#endif

    size_t ncat          = site_rate->getNRate();
    size_t ncat_mix      = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t nmix          = (model_factory->fused_mix_rate) ? 1 : model->getNMixtures();
    size_t block         = ncat_mix * nstates;
    size_t orig_nptn     = aln->size();
    size_t max_orig_nptn = roundUpToMultiple(orig_nptn, VectorClass::size());
    size_t nptn          = max_orig_nptn+model_factory->unobserved_ptns.size();
    ASCType ASC_type     = model_factory->getASC();
    bool ASC_Holder      = (ASC_type == ASC_VARIANT_MISSING || ASC_type == ASC_INFORMATIVE_MISSING);
    bool ASC_Lewis       = (ASC_type == ASC_VARIANT || ASC_type == ASC_INFORMATIVE);
    ASSERT(!ASC_Holder && "Holder's ascertainment bias correction not supported for this mixlen model");

//    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix], cat_id[ncat_mix];
//    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
//    for (c = 0; c < ncat_mix; c++) {
//        size_t m = c/denom;
//        cat_id[c] = c%ncat;
//        mix_addr_nstates[c] = m*nstates;
//        mix_addr[c] = mix_addr_nstates[c]*nstates;
//    }

    double *eval = model->getEigenvalues();
    ASSERT(eval);

    double *buffer_partial_lh_ptr = buffers.buffer_partial_lh;
    vector<size_t> limits;
    computePatternPacketBounds(VectorClass::size(), num_threads,
                               num_packets, nptn, limits);

    ASSERT(buffers.theta_all);

    double *val0 = NULL;
    double *val1 = NULL;
    double *val2 = NULL;
    double cat_rate[ncat];
    double cat_prop[ncat];

    int cur_mixlen = getCurMixture();

    if (SITE_MODEL) {
        for (size_t c = 0; c < ncat; c++) {
            cat_rate[c] = site_rate->getRate(c);
            cat_prop[c] = site_rate->getProp(c);
        }
    } else {
        val0 = buffer_partial_lh_ptr;
        val1 = val0 + get_safe_upper_limit(block);
        val2 = val1 + get_safe_upper_limit(block);
        buffer_partial_lh_ptr += 3*get_safe_upper_limit(block);
        double len = dad_branch->getLength(cur_mixlen);
        for (size_t c = 0; c < nmix; c++) {
            size_t cur_mix = (model_factory->fused_mix_rate) ? cur_mixlen : c;
            double *eval_ptr = eval+cur_mix*nstates;
            double prop = model->getMixtureWeight(cur_mix);
            size_t addr = c*nstates;
            for (size_t i = 0; i < nstates; i++) {
                double cof = eval_ptr[i];
                double val = exp(cof*len) * prop;
                double val1_ = cof*val;
                val0[addr+i] = val;
                val1[addr+i] = val1_;
                val2[addr+i] = cof*val1_;
            }
        }
    }

    double all_df(0.0), all_ddf(0.0), all_prob_const(0.0), all_df_const(0.0), all_ddf_const(0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(num_threads) reduction(+:all_df,all_ddf,all_prob_const,all_df_const,all_ddf_const)
#endif
    for (int packet_id = 0; packet_id < num_packets; packet_id++) {
        VectorClass my_df(0.0), my_ddf(0.0), vc_prob_const(0.0), vc_df_const(0.0), vc_ddf_const(0.0);
        size_t ptn_lower = limits[packet_id];
        size_t ptn_upper = limits[packet_id+1];

        if (!buffers.theta_computed)
        #ifdef KERNEL_FIX_STATES
            computeLikelihoodBufferSIMD<VectorClass, SAFE_NUMERIC, nstates, FMA, SITE_MODEL>(dad_branch, dad, ptn_lower, ptn_upper, packet_id, buffers);
        #else
            computeLikelihoodBufferGenericSIMD<VectorClass, SAFE_NUMERIC, FMA, SITE_MODEL>(dad_branch, dad, ptn_lower, ptn_upper, packet_id, buffers);
        #endif

        // mixed branch length model
        VectorClass lh_ptn;
        VectorClass df_ptn, ddf_ptn;

        for (size_t ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            lh_ptn = df_ptn = ddf_ptn = 0.0;
            VectorClass *theta = ((VectorClass*)(buffers.theta_all + ptn*block)) + cur_mixlen*nstates*nmix;
            double *val0_ptr = val0;
            double *val1_ptr = val1;
            double *val2_ptr = val2;
            for (size_t c = 0; c < nmix; c++) {
            #ifdef KERNEL_FIX_STATES
                dotProductTriple<VectorClass, double, nstates, FMA, true>(val0_ptr, val1_ptr, val2_ptr, theta, lh_ptn, df_ptn, ddf_ptn, nstates);
            #else
                dotProductTriple<VectorClass, double, FMA, true>(val0_ptr, val1_ptr, val2_ptr, theta, lh_ptn, df_ptn, ddf_ptn,nstates, nstates);
            #endif
                val0_ptr += nstates;
                val1_ptr += nstates;
                val2_ptr += nstates;
                theta += nstates;
            }
            lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
            if (ptn < orig_nptn) {
                VectorClass freq;
                freq.load_a(&ptn_freq[ptn]);
                VectorClass inv_lh_ptn = 1.0 / lh_ptn;

                // compute gradient (my_df)
                df_ptn  *= inv_lh_ptn;
                ddf_ptn *= inv_lh_ptn;
                my_df    = mul_add(df_ptn, freq, my_df);
                my_ddf  += nmul_add(df_ptn, df_ptn, ddf_ptn) * freq;
            } else {
                vc_prob_const += lh_ptn;
                vc_df_const   += df_ptn;
                vc_ddf_const  += ddf_ptn;
            }
        } // FOR ptn
        {
            //These additions don't need to be in a critical section, because of the
            //reduction(+:all_df,all_ddf,all_prob_const,all_df_const,all_ddf_const)
            //clause in the omp parallel for.
            all_df  += horizontal_add(my_df);
            all_ddf += horizontal_add(my_ddf);
            if (ASC_Lewis) {
                all_prob_const += horizontal_add(vc_prob_const);
                all_df_const   += horizontal_add(vc_df_const);
                all_ddf_const  += horizontal_add(vc_ddf_const);
            }
        }
    } // FOR packet

    // mark buffer as computed
    buffers.theta_computed = true;

    df  = all_df;
    ddf = all_ddf;

    if (!SAFE_NUMERIC && !std::isfinite(df)) {
        outError("Numerical underflow (lh-derivative-mixlen). Run again with the safe likelihood kernel via `-safe` option");
    }
	if (ASC_Lewis) {
        all_prob_const = 1.0/(1.0 - all_prob_const);
        // ascertainment bias correction
        all_df_const  *= all_prob_const;
        all_ddf_const *= all_prob_const;
        double nsites = aln->getNSite();
        df  += nsites * all_df_const;
        ddf += nsites * (all_ddf_const + all_df_const*all_df_const);
    }
    if (!std::isfinite(df)) {
        cout << "WARNING: Numerical underflow for lh-derivative-mixlen" << endl;
        df = ddf = 0.0;
    }
}

#endif //PHYLOKERNELNEW_H_
