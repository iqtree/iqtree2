/*
 * phylokernelnew.h
 * Newly revised kernel based on vectorizing over alignment patterns
 *
 *  Created on: Sept 23, 2016
 *      Author: minh
 */


#if !defined(PHYLOKERNELNEW_H_) || !defined(PHYLOKERNELNEW_STATE_H_)

#ifdef KERNEL_FIX_STATES
#   define PHYLOKERNELNEW_STATE_H_
#else
#   define PHYLOKERNELNEW_H_
#endif

#include "phylotree.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//#include <thread>

using namespace std;

/*******************************************************
 *
 * Helper function for vectors and matrix multiplication
 *
 ******************************************************/

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
inline void sumVec(VectorClass *A, VectorClass &X, size_t N)
{
    if (N == 1) {
        if (append)
            X += A[0];
        else
            X = A[0];
        return;
    }

    size_t i;
    switch (N % 4) {
    case 0: {
        VectorClass V[4];
        V[0] = A[0];
        V[1] = A[1];
        V[2] = A[2];
        V[3] = A[3];
        for (i = 4; i < N; i+=4) {
            V[0] += A[i];
            V[1] += A[i+1];
            V[2] += A[i+2];
            V[3] += A[i+3];
        }
        if (append)
            X += (V[0] + V[1]) + (V[2] + V[3]);
        else
            X = (V[0] + V[1]) + (V[2] + V[3]);
        break;
    }

    case 2: {
        VectorClass V[2];
        V[0] = A[0];
        V[1] = A[1];
        for (i = 2; i < N; i+=2) {
            V[0] += A[i];
            V[1] += A[i+1];
        }
        if (append)
            X += V[0] + V[1];
        else
            X = V[0] + V[1];
        break;
    }

    default: {
        VectorClass V[2];
        // odd N
        V[0] = A[0];
        V[1] = A[1];
        for (i = 2; i < N-1; i+=2) {
            V[0] += A[i];
            V[1] += A[i+1];
        }
        if (append)
            X += A[N-1] + V[0] + V[1];
        else
            X = A[N-1] + V[0] + V[1];
        break;
    }
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
inline void dotProductVec(Numeric *A, VectorClass *B, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductVec(Numeric *A, VectorClass *B, VectorClass &X, size_t N)
#endif
{
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
inline void dotProductDualVec(Numeric *A, VectorClass *B, Numeric *C, VectorClass *D, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductDualVec(Numeric *A, VectorClass *B, Numeric *C, VectorClass *D, VectorClass &X, size_t N)
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
inline void productVecMat(VectorClass *A, Numeric *M, VectorClass *X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void productVecMat(VectorClass *A, Numeric *M, VectorClass *X, size_t N)
#endif
{
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
    case 4:
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
inline void productVecMat(VectorClass *A, Numeric *M, VectorClass *X, VectorClass &Xmax)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void productVecMat(VectorClass *A, Numeric *M, VectorClass *X, VectorClass &Xmax, size_t N)
#endif
{
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
    case 4:
        for (i = 0; i < 4; i++) {
            // manual unrolling
            X[i] = mul_add(A[1],M[1],A[0]*M[0]) + mul_add(A[3],M[3],A[2]*M[2]);
            Xmax = max(Xmax, abs(X[i]));
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
template <class VectorClass, class Numeric, const size_t nstates, const bool FMA, const bool ADD>
inline void dotProductTriple(Numeric *A, Numeric *B, Numeric *C, VectorClass *D,
    VectorClass &X, VectorClass &Y, VectorClass &Z, size_t N)
#else
template <class VectorClass, class Numeric, const bool FMA, const bool ADD>
inline void dotProductTriple(Numeric *A, Numeric *B, Numeric *C, VectorClass *D,
    VectorClass &X, VectorClass &Y, VectorClass &Z, size_t N, size_t nstates)
#endif
{
    size_t i, j;
    if (nstates % 2 == 0) {
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
inline void dotProduct3Vec(Numeric *A, VectorClass *B, VectorClass *C, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProduct3Vec(Numeric *A, VectorClass *B, VectorClass *C, VectorClass &X, size_t N)
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
inline void dotProductExp(VectorClass *A, VectorClass *B, VectorClass *C, Numeric D, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductExp(VectorClass *A, VectorClass *B, VectorClass *C, Numeric D, VectorClass &X, size_t N)
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
inline void dotProductExp(VectorClass *A, VectorClass *B, Numeric D, VectorClass &X)
#else
template <class VectorClass, class Numeric, const bool FMA>
inline void dotProductExp(VectorClass *A, VectorClass *B, Numeric D, VectorClass &X, size_t N)
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

    PhyloNode *dad = info.dad, *node = (PhyloNode*)info.dad_branch->node;
    double *echild = info.echildren;
    double *partial_lh_leaf = info.partial_lh_leaves;

    //----------- Non-reversible model --------------

    if (!model->isReversible() || params->kernel_nonrev) {
        size_t nstatesqr = nstates*nstates;
        // non-reversible model
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNeighbor *child = (PhyloNeighbor*)*it;
            // precompute information buffer
            if (child->direction == TOWARD_ROOT) {
                // tranpose probability matrix
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
            if (isRootLeaf(child->node)) {
                for (c = 0; c < ncat_mix; c++) {
                    size_t m = c/denom;
                    model->getStateFrequency(partial_lh_leaf + c*nstates, m);
                }
                partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
            } else if (child->node->isLeaf()) {
                vector<int>::iterator it;
                if (nstates % VectorClass::size() == 0) {
                    // vectorized version
                    for (it = aln->seq_states[child->node->id].begin(); it != aln->seq_states[child->node->id].end(); it++) {
                        VectorClass *this_tip_partial_lh = (VectorClass*)&tip_partial_lh[(*it)*nstates];
                        double *this_partial_lh_leaf = &partial_lh_leaf[(*it)*block];
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
                    for (it = aln->seq_states[child->node->id].begin(); it != aln->seq_states[child->node->id].end(); it++) {
                        double *this_tip_partial_lh = &tip_partial_lh[(*it)*nstates];
                        double *this_partial_lh_leaf = &partial_lh_leaf[(*it)*block];
                        double *echild_ptr = echild;
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
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNeighbor *child = (PhyloNeighbor*)*it;
            VectorClass *echild_ptr = (VectorClass*)echild;
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
            if (child->node->isLeaf()) {
                vector<int>::iterator it;

                for (it = aln->seq_states[child->node->id].begin(); it != aln->seq_states[child->node->id].end(); it++) {
                    int state = (*it);
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
//        aligned_free(expchild);
    } else {
        // non-vectorized version
        double expchild[nstates];
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNeighbor *child = (PhyloNeighbor*)*it;
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
            if (child->node->isLeaf()) {
                vector<int>::iterator it;
                for (it = aln->seq_states[child->node->id].begin(); it != aln->seq_states[child->node->id].end(); it++) {
                    int state = (*it);
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

#ifndef KERNEL_FIX_STATES
template<class VectorClass>
inline void computeBounds(int threads, size_t elements, vector<size_t> &limits) {
    limits.reserve(threads+1);
    elements = ((elements+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t rest_elem = elements;
    limits.push_back(0);
    size_t last = 0;
    for (int rest_thread = threads; rest_thread > 1; rest_thread--) {
        size_t block_size = rest_elem/rest_thread;
        if (rest_elem % rest_thread != 0) block_size++;
        // padding to the vector size
        block_size = ((block_size+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();

        last += block_size;
        if (last >= elements)
            break;
        limits.push_back(last);
        rest_elem -= block_size;
    }

    limits.push_back(elements);
    if (limits.size() != threads+1) {
        if (Params::getInstance().num_threads == 0)
            outError("Too many threads may slow down analysis [-nt option]. Reduce threads");
        else
            outError("Too many threads may slow down analysis [-nt option]. Reduce threads or use -nt AUTO to automatically determine it");
    }
}
#endif

#ifdef KERNEL_FIX_STATES
template<class VectorClass, const int nstates>
#else
template<class VectorClass>
#endif
void PhyloTree::computeTraversalInfo(PhyloNode *node, PhyloNode *dad, bool compute_partial_lh) {

    if (!tip_partial_lh_computed)
        computeTipPartialLikelihood();

    traversal_info.clear();
#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    // reserve beginning of buffer_partial_lh for other purpose
    size_t ncat_mix = (model_factory->fused_mix_rate) ? site_rate->getNRate() : site_rate->getNRate()*model->getNMixtures();
    size_t block = aln->num_states * ncat_mix;
    double *buffer = buffer_partial_lh + block*VectorClass::size()*num_threads + get_safe_upper_limit(block)*(aln->STATE_UNKNOWN+2);

    // more buffer for non-reversible models
    if (!model->isReversible() || params->kernel_nonrev) {
        buffer += get_safe_upper_limit(3*block*nstates);
        buffer += get_safe_upper_limit(block)*(aln->STATE_UNKNOWN+1)*2;
        buffer += block*2*VectorClass::size()*num_threads;
    }

    // sort subtrees for mem save technique
    if (params->lh_mem_save == LM_MEM_SAVE) {
//        sortNeighborBySubtreeSize(node, dad);
//        sortNeighborBySubtreeSize(dad, node);
        int node_size = node->computeSize(dad);
        int dad_size = dad->computeSize(node);
//        PhyloNeighbor *dad_branch = (PhyloNeighbor*)dad->findNeighbor(node);
//        PhyloNeighbor *node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
        if (node_size < dad_size) {
            // swap node and dad due to tree size
            PhyloNode *tmp = node;
            node = dad;
            dad = tmp;
        }

    }

    PhyloNeighbor *dad_branch = (PhyloNeighbor*)dad->findNeighbor(node);
    PhyloNeighbor *node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
    bool dad_locked = computeTraversalInfo(dad_branch, dad, buffer);
    bool node_locked = computeTraversalInfo(node_branch, node, buffer);
    if (params->lh_mem_save == LM_MEM_SAVE) {
        if (dad_locked)
            mem_slots.unlock(dad_branch);
        if (node_locked)
            mem_slots.unlock(node_branch);
    }

    if (verbose_mode >= VB_DEBUG && traversal_info.size() > 0) {
        Node *saved = root;
        root = dad;
        drawTree(cout);
        root = saved;
    }

    if (traversal_info.empty())
        return;

    if (!model->isSiteSpecificModel()) {

        int num_info = traversal_info.size();

        if (verbose_mode >= VB_DEBUG) {
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
                    if (it->dad_branch->partial_lh_computed)
                        cout << " [";
                    else
                        cout << " (";
                    cout << mem_slots.findNei(it->dad_branch) - mem_slots.begin();
                    if (it->dad_branch->partial_lh_computed)
                        cout << "]";
                    else
                        cout << ")";
                }
            }
            cout << endl;
        }

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

    if (compute_partial_lh) {
        vector<size_t> limits;
        size_t orig_nptn = ((aln->size()+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
        size_t nptn = ((orig_nptn+model_factory->unobserved_ptns.size()+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
        computeBounds<VectorClass>(num_threads, nptn, limits);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static, 1) num_threads(num_threads)
        #endif
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
                computePartialLikelihood(*it, limits[thread_id], limits[thread_id+1], thread_id);
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
void PhyloTree::computePartialLikelihoodSIMD(TraversalInfo &info, size_t ptn_lower, size_t ptn_upper, int thread_id)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computePartialLikelihoodGenericSIMD(TraversalInfo &info, size_t ptn_lower, size_t ptn_upper, int thread_id)
#endif
{

    PhyloNeighbor *dad_branch = info.dad_branch;
    PhyloNode *dad = info.dad;
    // don't recompute the likelihood
	ASSERT(dad);
//    if (dad_branch->partial_lh_computed & 1)
//        return;
//    dad_branch->partial_lh_computed |= 1;
    PhyloNode *node = (PhyloNode*)(dad_branch->node);


#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    const size_t states_square = nstates*nstates;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
//    size_t max_nptn = ((nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();

//    if (!tip_partial_lh_computed)
//        computeTipPartialLikelihood();

	if (node->isLeaf()) {
//	    dad_branch->lh_scale_factor = 0.0;
		return;
	}
    
    size_t ptn, c;
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }
    size_t i, x;
    size_t block = nstates * ncat_mix;
//    size_t tip_block = nstates * model->getNMixtures();
    size_t tip_mem_size = max_orig_nptn * nstates;
//    size_t scale_size = SAFE_NUMERIC ? max_nptn * ncat_mix : max_nptn;
    size_t scale_size = SAFE_NUMERIC ? (ptn_upper-ptn_lower) * ncat_mix : (ptn_upper-ptn_lower);

	double *evec = model->getEigenvectors();
	double *inv_evec = model->getInverseEigenvectors();
	ASSERT(inv_evec && evec);
	double *eval = model->getEigenvalues();

	// internal node
	PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
	FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNeighbor *nei = (PhyloNeighbor*)(*it);
        // make sure that the partial_lh of children are different!
        ASSERT(dad_branch->partial_lh != nei->partial_lh);
		if (!left) left = nei; else right = nei;
	}

    // precomputed buffer to save times
    size_t thread_buf_size = (2*block+nstates)*VectorClass::size();
    double *buffer_partial_lh_ptr = buffer_partial_lh + (getBufferPartialLhSize() - thread_buf_size*num_threads);
    double *echildren = NULL;
    double *partial_lh_leaves = NULL;

    // pre-compute scaled branch length per category
    double len_children[ncat*(node->degree()-1)]; // +1 in case num_leaves = 0
    double *len_left = NULL, *len_right = NULL;

    if (SITE_MODEL) {
        double *len_children_ptr = len_children;
        FOR_NEIGHBOR_IT(node, dad, it3) {
            for (c = 0; c < ncat; c++) {
                len_children_ptr[c] = site_rate->getRate(c) * (*it3)->length;
            }
            if (!len_left)
                len_left = len_children_ptr;
            else
                len_right = len_children_ptr;
            len_children_ptr += ncat;
        }
    } else {

        echildren = info.echildren;
        partial_lh_leaves = info.partial_lh_leaves;

    }

    double *eleft = echildren, *eright = echildren + block*nstates;

	if (!left->node->isLeaf() && right->node->isLeaf()) {
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
        double *etmp = eleft;
        eleft = eright;
        eright = etmp;
        etmp = len_left;
        len_left = len_right;
        len_right = etmp;
	}

    if (node->degree() > 3) {
        /*--------------------- multifurcating node ------------------*/

        // now for-loop computing partial_lh over all site-patterns
        VectorClass *partial_lh_all = (VectorClass*) &buffer_partial_lh_ptr[thread_buf_size*thread_id];
        double *vec_tip = (double*)&partial_lh_all[block];

        for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            for (i = 0; i < block; i++)
                partial_lh_all[i] = 1.0;
            UBYTE *scale_dad = NULL;
            if (SAFE_NUMERIC) {
                scale_dad = dad_branch->scale_num + ptn*ncat_mix;
                memset(scale_dad, 0, sizeof(UBYTE)*ncat_mix*VectorClass::size());
            } else
                memset(&dad_branch->scale_num[ptn], 0, sizeof(UBYTE)*VectorClass::size());

            // SITE_MODEL variables
            VectorClass *expchild = partial_lh_all + block;
            VectorClass *eval_ptr = (VectorClass*) &eval[ptn*nstates];
            VectorClass *evec_ptr = (VectorClass*) &evec[ptn*states_square];
            double *len_child = len_children;
            VectorClass vchild;

            // normal model
            double *partial_lh_leaf = partial_lh_leaves;
            double *echild = echildren;

            FOR_NEIGHBOR_IT(node, dad, it) {
                if (SITE_MODEL) {
                    PhyloNeighbor *child = (PhyloNeighbor*)*it;
                    UBYTE *scale_child = SAFE_NUMERIC ? child->scale_num + ptn*ncat_mix : NULL;
                    VectorClass *partial_lh = partial_lh_all;
                    if (child->node->isLeaf()) {
                        // external node
                        VectorClass *tip_partial_lh_child = (VectorClass*) &tip_partial_lh[child->node->id*tip_mem_size + ptn*nstates];
                        for (c = 0; c < ncat; c++) {
                            for (i = 0; i < nstates; i++)
                                expchild[i] = exp(eval_ptr[i]*len_child[c]) * tip_partial_lh_child[i];
                            for (x = 0; x < nstates; x++) {
                                VectorClass *this_evec = &evec_ptr[x*nstates];
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
                            for (i = 0; i < VectorClass::size(); i++)
                                dad_branch->scale_num[ptn+i] += child->scale_num[ptn+i];
                        }

                        for (c = 0; c < ncat_mix; c++) {
                            if (SAFE_NUMERIC) {
                                for (x = 0; x < VectorClass::size(); x++)
                                    scale_dad[x*ncat_mix+c] += scale_child[x*ncat_mix+c];
                            }
                            // compute real partial likelihood vector
                            for (i = 0; i < nstates; i++)
                                expchild[i] = exp(eval_ptr[i]*len_child[c]) * partial_lh_child[i];
                            for (x = 0; x < nstates; x++) {
                                VectorClass *this_evec = &evec_ptr[x*nstates];
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
                    PhyloNeighbor *child = (PhyloNeighbor*)*it;
                    UBYTE *scale_child = SAFE_NUMERIC ? child->scale_num + ptn*ncat_mix : NULL;
                    if (child->node->isLeaf()) {
                        // external node
                        // load data for tip
                        for (i = 0; i < VectorClass::size(); i++) {
                            double *child_lh;
                            if (ptn+i < orig_nptn)
                                child_lh = partial_lh_leaf + block*(aln->at(ptn+i))[child->node->id];
                            else if (ptn+i < max_orig_nptn)
                                child_lh = partial_lh_leaf + block*aln->STATE_UNKNOWN;
                            else if (ptn+i < nptn)
                                child_lh = partial_lh_leaf + block*model_factory->unobserved_ptns[ptn+i-max_orig_nptn];
                            else
                                child_lh = partial_lh_leaf + block*aln->STATE_UNKNOWN;
                            double *this_vec_tip = vec_tip+i;
                            for (c = 0; c < block; c++) {
                                *this_vec_tip = child_lh[c];
                                this_vec_tip += VectorClass::size();
                            }
                        }
                        VectorClass *vtip = (VectorClass*)vec_tip;
                        for (c = 0; c < block; c++) {
                            // compute real partial likelihood vector
                            partial_lh_all[c] *= vtip[c];
                        }
                        partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
                    } else {
                        // internal node
                        VectorClass *partial_lh = partial_lh_all;
                        VectorClass *partial_lh_child = (VectorClass*)(child->partial_lh + ptn*block);
                        if (!SAFE_NUMERIC) {
                            for (i = 0; i < VectorClass::size(); i++)
                                dad_branch->scale_num[ptn+i] += child->scale_num[ptn+i];
                        }

                        double *echild_ptr = echild;
                        for (c = 0; c < ncat_mix; c++) {
                            if (SAFE_NUMERIC) {
                                for (x = 0; x < VectorClass::size(); x++)
                                    scale_dad[x*ncat_mix+c] += scale_child[x*ncat_mix+c];
                            }
                            // compute real partial likelihood vector
                            for (x = 0; x < nstates; x++) {
                                VectorClass vchild = echild_ptr[0] * partial_lh_child[0];
    //                            double *echild_ptr = echild + (c*nstatesqr+x*nstates);
                                for (i = 1; i < nstates; i++) {
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
                    for (c = 0; c < ncat_mix; c++) {
                        VectorClass lh_max = 0.0;
                        for (x = 0; x < nstates; x++)
                            lh_max = max(lh_max,abs(partial_lh_tmp[x]));
                        // check if one should scale partial likelihoods
                        auto underflown = ((lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                        if (horizontal_or(underflown)) { // at least one site has numerical underflown
                            for (x = 0; x < VectorClass::size(); x++)
                            if (underflown[x]) {
                                // BQM 2016-05-03: only scale for non-constant sites
                                // now do the likelihood scaling
                                double *partial_lh = (double*)partial_lh_tmp + (x);
                                for (i = 0; i < nstates; i++)
                                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                                dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
                            }
                        }
                        partial_lh_tmp += nstates;
                    }
                } else {
                    // not -safe numeric
                    VectorClass lh_max = 0.0;
                    for (x = 0; x < block; x++)
                        lh_max = max(lh_max,abs(partial_lh_all[x]));
                    auto underflown = (lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                    if (horizontal_or(underflown)) { // at least one site has numerical underflown
                        for (x = 0; x < VectorClass::size(); x++)
                        if (underflown[x]) {
                            double *partial_lh = (double*)partial_lh_all + (x);
                            // now do the likelihood scaling
                            for (i = 0; i < block; i++) {
                                partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                            }
    //                        sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                            dad_branch->scale_num[ptn+x] += 1;
                        }
                    }
                } // if-else

            } // FOR_NEIGHBOR

        
            // compute dot-product with inv_eigenvector
            VectorClass *partial_lh_tmp = partial_lh_all;
            VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;
            double *inv_evec_ptr = SITE_MODEL ? &inv_evec[ptn*states_square] : NULL;
            for (c = 0; c < ncat_mix; c++) {
                if (SITE_MODEL) {
                    // compute dot-product with inv_eigenvector
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, VectorClass, nstates, FMA>(partial_lh_tmp, (VectorClass*)inv_evec_ptr, partial_lh);
#else
                    productVecMat<VectorClass, VectorClass, FMA> (partial_lh_tmp, (VectorClass*)inv_evec_ptr, partial_lh, nstates);
#endif
                } else {
                    inv_evec_ptr = inv_evec + mix_addr[c];
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, double, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max);
#else
                    productVecMat<VectorClass, double, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max, nstates);
#endif
                }
                partial_lh += nstates;
                partial_lh_tmp += nstates;
            }

        } // for ptn

        // end multifurcating treatment
    } else if (left->node->isLeaf() && right->node->isLeaf()) {

        /*--------------------- TIP-TIP (cherry) case ------------------*/

        double *partial_lh_left = SITE_MODEL ? &tip_partial_lh[left->node->id * tip_mem_size] : partial_lh_leaves;
        double *partial_lh_right = SITE_MODEL ? &tip_partial_lh[right->node->id * tip_mem_size] : partial_lh_leaves + (aln->STATE_UNKNOWN+1)*block;

		// scale number must be ZERO
	    memset(dad_branch->scale_num + (SAFE_NUMERIC ? ptn_lower*ncat_mix : ptn_lower), 0, scale_size * sizeof(UBYTE));
        double *vec_left = buffer_partial_lh_ptr + thread_buf_size*thread_id;

        double *vec_right =  SITE_MODEL ? &vec_left[nstates*VectorClass::size()] : &vec_left[block*VectorClass::size()];
        VectorClass *partial_lh_tmp = SITE_MODEL ? (VectorClass*)vec_right+nstates : (VectorClass*)vec_right+block;

		for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);

            if (SITE_MODEL) {
                VectorClass* expleft = (VectorClass*) vec_left;
                VectorClass* expright = (VectorClass*) vec_right;
                VectorClass *vleft = (VectorClass*) &partial_lh_left[ptn*nstates];
                VectorClass *vright = (VectorClass*) &partial_lh_right[ptn*nstates];
                VectorClass *eval_ptr = (VectorClass*) &eval[ptn*nstates];
                VectorClass *evec_ptr = (VectorClass*) &evec[ptn*states_square];
                VectorClass *inv_evec_ptr = (VectorClass*) &inv_evec[ptn*states_square];
                for (c = 0; c < ncat; c++) {
                    for (i = 0; i < nstates; i++) {
                        expleft[i] = exp(eval_ptr[i]*len_left[c]) * vleft[i];
                        expright[i] = exp(eval_ptr[i]*len_right[c]) * vright[i];

                    }
                    // compute real partial likelihood vector
                    for (x = 0; x < nstates; x++) {
                        VectorClass *this_evec = evec_ptr + x*nstates;
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
                VectorClass *vleft = (VectorClass*)vec_left;
                VectorClass *vright = (VectorClass*)vec_right;
                // load data for tip
                for (x = 0; x < VectorClass::size(); x++) {
                    double *tip_left, *tip_right;
                    if (ptn+x < orig_nptn) {
                        tip_left  = partial_lh_left  + block * (aln->at(ptn+x))[left->node->id];
                        tip_right = partial_lh_right + block * (aln->at(ptn+x))[right->node->id];
                    } else if (ptn+x < max_orig_nptn) {
                        tip_left  = partial_lh_left  + block * aln->STATE_UNKNOWN;
                        tip_right = partial_lh_right + block * aln->STATE_UNKNOWN;
                    } else if (ptn+x < nptn) {
                        tip_left  = partial_lh_left  + block * model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                        tip_right = partial_lh_right + block * model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                    } else {
                        tip_left  = partial_lh_left  + block * aln->STATE_UNKNOWN;
                        tip_right = partial_lh_right + block * aln->STATE_UNKNOWN;
                    }
                    double *this_vec_left = vec_left+x;
                    double *this_vec_right = vec_right+x;
                    for (i = 0; i < block; i++) {
                        *this_vec_left = tip_left[i];
                        *this_vec_right = tip_right[i];
                        this_vec_left += VectorClass::size();
                        this_vec_right += VectorClass::size();
                    }
                }


                for (c = 0; c < ncat_mix; c++) {
                    double *inv_evec_ptr = inv_evec + mix_addr[c];
                    // compute real partial likelihood vector
                    for (x = 0; x < nstates; x++) {
                        partial_lh_tmp[x] = vleft[x] * vright[x];
                    }

                    // compute dot-product with inv_eigenvector
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, double, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh);
#else
                    productVecMat<VectorClass, double, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, nstates);
#endif

                    // increase pointer
                    vleft += nstates;
                    vright += nstates;
                    partial_lh += nstates;
                } // FOR category
            } // IF SITE_MODEL
		} // FOR LOOP


	} else if (left->node->isLeaf() && !right->node->isLeaf()) {

        /*--------------------- TIP-INTERNAL NODE case ------------------*/

		// only take scale_num from the right subtree
		memcpy(
            dad_branch->scale_num + (SAFE_NUMERIC ? ptn_lower*ncat_mix : ptn_lower),
            right->scale_num + (SAFE_NUMERIC ? ptn_lower*ncat_mix : ptn_lower),
            scale_size * sizeof(UBYTE));

        double *partial_lh_left = SITE_MODEL ? &tip_partial_lh[left->node->id * tip_mem_size] : partial_lh_leaves;


        double *vec_left = buffer_partial_lh_ptr + thread_buf_size*thread_id;
        VectorClass *partial_lh_tmp = SITE_MODEL ? (VectorClass*)vec_left+2*nstates : (VectorClass*)vec_left+block;

		for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
//            memset(partial_lh, 0, sizeof(VectorClass)*block);
            VectorClass lh_max = 0.0;

            if (SITE_MODEL) {
                VectorClass *expleft = (VectorClass*)vec_left;
                VectorClass *expright = expleft+nstates;
                VectorClass *vleft = (VectorClass*)&partial_lh_left[ptn*nstates];
                VectorClass *eval_ptr = (VectorClass*) &eval[ptn*nstates];
                VectorClass *evec_ptr = (VectorClass*) &evec[ptn*states_square];
                VectorClass *inv_evec_ptr = (VectorClass*) &inv_evec[ptn*states_square];
                for (c = 0; c < ncat; c++) {
                    for (i = 0; i < nstates; i++) {
                        expleft[i] = exp(eval_ptr[i]*len_left[c]) * vleft[i];
                        expright[i] = exp(eval_ptr[i]*len_right[c]) * partial_lh_right[i];
                    }
                    // compute real partial likelihood vector
                    for (x = 0; x < nstates; x++) {
                        VectorClass *this_evec = evec_ptr + x*nstates;
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
                            for (x = 0; x < VectorClass::size(); x++)
                            if (underflown[x]) {
                                // BQM 2016-05-03: only scale for non-constant sites
                                // now do the likelihood scaling
                                double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                                for (i = 0; i < nstates; i++)
                                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                                dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
                            }
                        }
                    }
                    partial_lh_right += nstates;
                    partial_lh += nstates;
                } // FOR category

            } else {
                VectorClass *vleft = (VectorClass*)vec_left;
                // load data for tip
                for (x = 0; x < VectorClass::size(); x++) {
                    double *tip;
                    if (ptn+x < orig_nptn) {
                        tip = partial_lh_left + block*(aln->at(ptn+x))[left->node->id];
                    } else if (ptn+x < max_orig_nptn) {
                        tip = partial_lh_left + block*aln->STATE_UNKNOWN;
                    } else if (ptn+x < nptn) {
                        tip = partial_lh_left + block*model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                    } else {
                        tip = partial_lh_left + block*aln->STATE_UNKNOWN;
                    }
                    double *this_vec_left = vec_left+x;
                    for (i = 0; i < block; i++) {
                        *this_vec_left = tip[i];
                        this_vec_left += VectorClass::size();
                    }
                }

                double *eright_ptr = eright;
                for (c = 0; c < ncat_mix; c++) {
                    if (SAFE_NUMERIC)
                        lh_max = 0.0;
                    double *inv_evec_ptr = inv_evec + mix_addr[c];
                    // compute real partial likelihood vector
                    for (x = 0; x < nstates; x++) {
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
                            for (x = 0; x < VectorClass::size(); x++)
                            if (underflown[x]) {
                                // BQM 2016-05-03: only scale for non-constant sites
                                // now do the likelihood scaling
                                double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                                for (i = 0; i < nstates; i++)
                                    partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                                dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
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
                    for (x = 0; x < VectorClass::size(); x++)
                    if (underflown[x]) {
                        double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                        // now do the likelihood scaling
                        for (i = 0; i < block; i++) {
                            partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                        }
//                        sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                        dad_branch->scale_num[ptn+x] += 1;
                    }
                }
            }

		} // big for loop over ptn

	} else {

        /*--------------------- INTERNAL-INTERNAL NODE case ------------------*/

        VectorClass *partial_lh_tmp = (VectorClass*)(buffer_partial_lh_ptr + thread_buf_size*thread_id);
		for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_left = (VectorClass*)(left->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;
            UBYTE *scale_dad, *scale_left, *scale_right;

            if (SAFE_NUMERIC) {
                size_t addr = ptn*ncat_mix;
                scale_dad = dad_branch->scale_num + addr;
                scale_left = left->scale_num + addr;
                scale_right = right->scale_num + addr;
            } else {
                scale_dad = dad_branch->scale_num + ptn;
                scale_left = left->scale_num + ptn;
                scale_right = right->scale_num + ptn;
                for (i = 0; i < VectorClass::size(); i++)
                    scale_dad[i] = scale_left[i] + scale_right[i];
            }

            double *eleft_ptr = eleft;
            double *eright_ptr = eright;
            VectorClass *expleft, *expright, *eval_ptr, *evec_ptr, *inv_evec_ptr;
            if (SITE_MODEL) {
                expleft = partial_lh_tmp + nstates;
                expright = expleft + nstates;
                eval_ptr = (VectorClass*) &eval[ptn*nstates];
                evec_ptr = (VectorClass*) &evec[ptn*states_square];
                inv_evec_ptr = (VectorClass*) &inv_evec[ptn*states_square];
            }

			for (c = 0; c < ncat_mix; c++) {
                if (SAFE_NUMERIC) {
                    lh_max = 0.0;
                    for (x = 0; x < VectorClass::size(); x++)
                        scale_dad[x*ncat_mix] = scale_left[x*ncat_mix] + scale_right[x*ncat_mix];
                }

                if (SITE_MODEL) {
                    // site-specific model
                    for (i = 0; i < nstates; i++) {
                        expleft[i] = exp(eval_ptr[i]*len_left[c]) * partial_lh_left[i];
                        expright[i] = exp(eval_ptr[i]*len_right[c]) * partial_lh_right[i];
                    }
                    for (x = 0; x < nstates; x++) {
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
                    double *inv_evec_ptr = inv_evec + mix_addr[c];
                    // compute real partial likelihood vector
                    for (x = 0; x < nstates; x++) {
#ifdef KERNEL_FIX_STATES
                        dotProductDualVec<VectorClass, double, nstates, FMA>(eleft_ptr, partial_lh_left, eright_ptr, partial_lh_right, partial_lh_tmp[x]);
#else
                        dotProductDualVec<VectorClass, double, FMA>(eleft_ptr, partial_lh_left, eright_ptr, partial_lh_right, partial_lh_tmp[x], nstates);
#endif
                        eleft_ptr += nstates;
                        eright_ptr += nstates;
                    }
                    
                    // compute dot-product with inv_eigenvector
#ifdef KERNEL_FIX_STATES
                    productVecMat<VectorClass, double, nstates, FMA>(partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max);
#else
                    productVecMat<VectorClass, double, FMA> (partial_lh_tmp, inv_evec_ptr, partial_lh, lh_max, nstates);
#endif
                }

                // check if one should scale partial likelihoods
                if (SAFE_NUMERIC) {
                    auto underflown = ((lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                    if (horizontal_or(underflown))
                        for (x = 0; x < VectorClass::size(); x++)
                        if (underflown[x]) {
                            // BQM 2016-05-03: only scale for non-constant sites
                            // now do the likelihood scaling
                            double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                            for (i = 0; i < nstates; i++)
                                partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                            scale_dad[x*ncat_mix] += 1;
                        }
                    scale_dad++;
                    scale_left++;
                    scale_right++;
                }
                partial_lh_left += nstates;
                partial_lh_right += nstates;
                partial_lh += nstates;
			}

            if (!SAFE_NUMERIC) {
                // check if one should scale partial likelihoods
                auto underflown = (lh_max < SCALING_THRESHOLD) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                if (horizontal_or(underflown)) { // at least one site has numerical underflown
                    for (x = 0; x < VectorClass::size(); x++)
                    if (underflown[x]) {
                        double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                        // now do the likelihood scaling
                        for (i = 0; i < block; i++) {
                            partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                        }
//                        sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                        dad_branch->scale_num[ptn+x] += 1;
                    }
                }
            }

		} // big for loop over ptn

	}
}

/*******************************************************
 *
 * NEW! highly-vectorized log-likelihood derivative function
 *
 ******************************************************/


#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodBufferSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, size_t ptn_lower, size_t ptn_upper, int thread_id)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodBufferGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, size_t ptn_lower, size_t ptn_upper, int thread_id)
#endif
{
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    size_t ptn, i, c;
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
    size_t tip_block = nstates * model->getNMixtures();
    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    // reserve 3*block for computeLikelihoodDerv
    double *buffer_partial_lh_ptr = buffer_partial_lh + 3*get_safe_upper_limit(block);
    if (isMixlen()) {
        size_t nmix = getMixlen();
        buffer_partial_lh_ptr += nmix*(nmix+1)*VectorClass::size() + (nmix+3)*nmix*VectorClass::size()*num_threads;
    }

    // first compute partial_lh
    for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
        computePartialLikelihood(*it, ptn_lower, ptn_upper, thread_id);

    if (dad->isLeaf()) {
        // special treatment for TIP-INTERNAL NODE case
        double *tip_partial_lh_node = &tip_partial_lh[dad->id * max_orig_nptn*nstates];

        double *vec_tip = buffer_partial_lh_ptr + tip_block*VectorClass::size()*thread_id;

        for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
            //load tip vector
            if (!SITE_MODEL)
            for (i = 0; i < VectorClass::size(); i++) {
                double *this_tip_partial_lh;
                if (ptn+i < orig_nptn)
                    this_tip_partial_lh = tip_partial_lh + tip_block*(aln->at(ptn+i))[dad->id];
                else if (ptn+i < max_orig_nptn)
                    this_tip_partial_lh = tip_partial_lh + tip_block*aln->STATE_UNKNOWN;
                else if (ptn+i < nptn)
                    this_tip_partial_lh = tip_partial_lh + tip_block*model_factory->unobserved_ptns[ptn+i-max_orig_nptn];
                else
                    this_tip_partial_lh = tip_partial_lh + tip_block*aln->STATE_UNKNOWN;
                double *this_vec_tip = vec_tip+i;
                for (c = 0; c < tip_block; c++) {
                    *this_vec_tip = this_tip_partial_lh[c];
                    this_vec_tip += VectorClass::size();
                }

            }
            VectorClass *lh_tip;
            if (SITE_MODEL)
                lh_tip = (VectorClass*)&tip_partial_lh_node[ptn*nstates];
            for (c = 0; c < ncat_mix; c++) {
                if (!SITE_MODEL)
                    lh_tip = (VectorClass*)(vec_tip + mix_addr_nstates[c]*VectorClass::size());
                for (i = 0; i < nstates; i++) {
                    theta[i] = lh_tip[i] * partial_lh_dad[i];
                }
                partial_lh_dad += nstates;
                theta += nstates;
            }
            if (SAFE_NUMERIC) {
                // numerical scaling per category
                UBYTE *scale_dad;
                UBYTE min_scale;
                for (i = 0; i < VectorClass::size(); i++) {
                    scale_dad = dad_branch->scale_num+(ptn+i)*ncat_mix;
                    min_scale = scale_dad[0];
                    for (c = 1; c < ncat_mix; c++)
                        min_scale = min(min_scale, scale_dad[c]);

                    buffer_scale_all[ptn+i] = min_scale;

                    for (c = 0; c < ncat_mix; c++) {
                        if (scale_dad[c] == min_scale+1) {
                            double *this_theta = &theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] *= SCALING_THRESHOLD;
                            }
                        } else if (scale_dad[c] > min_scale+1) {
                            double *this_theta = &theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] = 0.0;
                            }
                        }
                    }
                }
            } else {
                // normal scaling
                for (i = 0; i < VectorClass::size(); i++)
                    buffer_scale_all[ptn+i] = dad_branch->scale_num[ptn+i];
            }
            VectorClass *buf = (VectorClass*)(buffer_scale_all+ptn);
            *buf *= LOG_SCALING_THRESHOLD;

        } // FOR PTN LOOP
//            aligned_free(vec_tip);
    } else {
        //------- both dad and node are internal nodes  --------//

        // now compute theta
        for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
            VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
            VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            for (i = 0; i < block; i++)
                theta[i] = partial_lh_node[i] * partial_lh_dad[i];

            if (SAFE_NUMERIC) {
                // numerical scaling per category
                UBYTE min_scale;
                UBYTE sum_scale[ncat_mix];
                size_t ptn_ncat = ptn*ncat_mix; 
                UBYTE *scale_dad = dad_branch->scale_num + ptn_ncat;
                UBYTE *scale_node = node_branch->scale_num + ptn_ncat;

                for (i = 0; i < VectorClass::size(); i++) {
                    min_scale = sum_scale[0] = scale_dad[0] + scale_node[0];
                    for (c = 1; c < ncat_mix; c++) {
                        sum_scale[c] = scale_dad[c] + scale_node[c];
                        min_scale = min(min_scale, sum_scale[c]);
                    }
                    buffer_scale_all[ptn+i] = min_scale;

                    for (c = 0; c < ncat_mix; c++) {
                        if (sum_scale[c] == min_scale+1) {
                            double *this_theta = &theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] *= SCALING_THRESHOLD;
                            }
                        } else if (sum_scale[c] > min_scale+1) {
                            double *this_theta = &theta_all[ptn*block + c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_theta[x*VectorClass::size()] = 0.0;
                            }
                        }
                    }
                    scale_dad += ncat_mix;
                    scale_node += ncat_mix;
                }
            } else {
                for (i = 0; i < VectorClass::size(); i++)
                    buffer_scale_all[ptn+i] = dad_branch->scale_num[ptn+i] + node_branch->scale_num[ptn+i];
            }
            VectorClass *buf = (VectorClass*)(buffer_scale_all+ptn);
            *buf *= LOG_SCALING_THRESHOLD;
        } // FOR ptn
    } // internal node
}

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double *df, double *ddf)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double *df, double *ddf)
#endif
{
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, false);
#endif

//
//    if ((dad_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(dad_branch, dad);
//    if ((node_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(node_branch, node);

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
//    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    bool isASC = model_factory->unobserved_ptns.size() > 0;



    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix], cat_id[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        cat_id[c] = c%ncat;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    double *eval = model->getEigenvalues();
    ASSERT(eval);

    double *buffer_partial_lh_ptr = buffer_partial_lh;
    vector<size_t> limits;
    computeBounds<VectorClass>(num_threads, nptn, limits);

	ASSERT(theta_all);

    double *val0 = NULL;
    double *val1 = NULL;
    double *val2 = NULL;
    double cat_rate[ncat];
    double cat_prop[ncat];


    if (SITE_MODEL) {
        for (c = 0; c < ncat; c++) {
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
            for (c = 0; c < ncat_mix; c++) {
                size_t m = c/denom;
                size_t mycat = c%ncat;
                double len = dad_branch->getLength(mycat);
                VectorClass *eval_ptr = (VectorClass*)(eval + mix_addr_nstates[c]);
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                double myrate = site_rate->getRate(mycat);
                for (i = 0; i < loop_size; i++) {
                    VectorClass cof = eval_ptr[i] * myrate;
                    VectorClass val = exp(cof*len) * prop;
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
            for (c = 0; c < ncat_mix; c++) {
                size_t m = c/denom;
                double *eval_ptr = eval + mix_addr_nstates[c];
                size_t mycat = c%ncat;
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                size_t addr = c*nstates;
                double len = dad_branch->getLength(mycat);
                for (i = 0; i < nstates; i++) {
                    double cof = eval_ptr[i]*site_rate->getRate(mycat);
                    double val = exp(cof*len) * prop;
                    double val1_ = cof*val;
                    val0[addr+i] = val;
                    val1[addr+i] = val1_;
                    val2[addr+i] = cof*val1_;
                }
            }
        }
    }

    double dad_length = dad_branch->length;

    VectorClass all_lh(0.0), all_df(0.0), all_ddf(0.0), all_prob_const(0.0), all_df_const(0.0), all_ddf_const(0.0);
    VectorClass *all_dfvec = NULL;
    VectorClass *all_ddfvec = NULL;

    size_t nmixlen = getMixlen(), nmixlen2 = nmixlen*nmixlen;
    if (isMixlen()) {
        ASSERT(nmixlen == ncat);
        all_dfvec = (VectorClass*)buffer_partial_lh_ptr;
        all_ddfvec = all_dfvec + nmixlen;
        buffer_partial_lh_ptr += nmixlen*(nmixlen+1)*VectorClass::size();
        for (i = 0; i < nmixlen; i++) all_dfvec[i] = 0.0;
        for (i = 0; i < nmixlen2; i++) all_ddfvec[i] = 0.0;
    }

//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) private(ptn, i, c) num_threads(num_threads)
#endif
    for (int thread_id = 0; thread_id < num_threads; thread_id++) {
        VectorClass my_df(0.0), my_ddf(0.0), vc_prob_const(0.0), vc_df_const(0.0), vc_ddf_const(0.0);
        size_t ptn_lower = limits[thread_id];
        size_t ptn_upper = limits[thread_id+1];

        if (!theta_computed)
        #ifdef KERNEL_FIX_STATES
            computeLikelihoodBufferSIMD<VectorClass, SAFE_NUMERIC, nstates, FMA, SITE_MODEL>(dad_branch, dad, ptn_lower, ptn_upper, thread_id);
        #else
            computeLikelihoodBufferGenericSIMD<VectorClass, SAFE_NUMERIC, FMA, SITE_MODEL>(dad_branch, dad, ptn_lower, ptn_upper, thread_id);
        #endif

        if (isMixlen()) {
            // mixed branch length model
            VectorClass lh_ptn;
            VectorClass *df_ptn = ((VectorClass*)buffer_partial_lh_ptr) + (nmixlen+3)*nmixlen*thread_id;
            VectorClass *ddf_ptn = df_ptn+nmixlen;

            VectorClass my_lh(0.0);
            VectorClass *my_df = df_ptn + nmixlen*2;
            VectorClass *my_ddf = df_ptn + nmixlen*3;
            for (i = 0; i < nmixlen; i++) my_df[i] = 0.0;
            for (i = 0; i < nmixlen2; i++) my_ddf[i] = 0.0;

            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                for (i = 0; i < nmixlen; i++)
                    df_ptn[i] = ddf_ptn[i] = 0.0;
//                lh_ptn.load_a(&ptn_invar[ptn]);
                lh_ptn = 0.0;
                VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
                double *val0_ptr = val0;
                double *val1_ptr = val1;
                double *val2_ptr = val2;
                for (c = 0; c < ncat_mix; c++) {
                    i = cat_id[c];
                #ifdef KERNEL_FIX_STATES
                    dotProductTriple<VectorClass, double, nstates, FMA, true>(val0_ptr, val1_ptr, val2_ptr, theta, lh_ptn, df_ptn[i], ddf_ptn[i], nstates);
                #else
                    dotProductTriple<VectorClass, double, FMA, true>(val0_ptr, val1_ptr, val2_ptr, theta, lh_ptn, df_ptn[i], ddf_ptn[i],nstates, nstates);
                #endif
                    val0_ptr += nstates;
                    val1_ptr += nstates;
                    val2_ptr += nstates;
                    theta += nstates;
                }
                lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
//                ASSERT(horizontal_and(lh_ptn > 0.0));

                if (ptn < orig_nptn) {
                    VectorClass freq;
                    freq.load_a(&ptn_freq[ptn]);
                    VectorClass inv_lh_ptn = 1.0 / lh_ptn;

                    // compute gradient (my_df)
                    for (i = 0; i < nmixlen; i++) {
                        df_ptn[i] *= inv_lh_ptn;
                        ddf_ptn[i] *= inv_lh_ptn;
                        my_df[i] = mul_add(df_ptn[i], freq, my_df[i]);
                    }

                    // now compute hessian matrix my_ddf
                    for (i = 0; i < nmixlen; i++) {
                        my_ddf[i*nmixlen+i] += nmul_add(df_ptn[i],df_ptn[i], ddf_ptn[i]) * freq;
                        for (c = 0; c < nmixlen; c++)
                            if (c!=i)
                                my_ddf[i*nmixlen+c] -= df_ptn[i]*df_ptn[c]*freq;
                    }


                    lh_ptn = log(lh_ptn) + VectorClass().load_a(&buffer_scale_all[ptn]);
                    my_lh = mul_add(lh_ptn, freq, my_lh);
                } else {
                    ASSERT(0 && "TODO +ASC not supported");
                }
            } // FOR ptn

        #ifdef _OPENMP
        #pragma omp critical
        #endif
            {
                for (i = 0; i < nmixlen; i++)
                    all_dfvec[i] += my_df[i];
                for (i = 0; i < nmixlen2; i++)
                    all_ddfvec[i] += my_ddf[i];
                all_lh += my_lh;
//                if (isASC) {
//                    all_prob_const += vc_prob_const;
//                    all_df_const += vc_df_const;
//                    all_ddf_const += vc_ddf_const;
//                }
            }

        } else {
            // normal joint branch length model
            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn;
                //lh_ptn.load_a(&ptn_invar[ptn]);
                VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
                VectorClass df_ptn, ddf_ptn;

                if (SITE_MODEL) {
                    VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                    lh_ptn = 0.0; df_ptn = 0.0; ddf_ptn = 0.0;
                    for (c = 0; c < ncat; c++) {
                        VectorClass lh_cat(0.0), df_cat(0.0), ddf_cat(0.0);
                        for (i = 0; i < nstates; i++) {
                            VectorClass cof = eval_ptr[i] * cat_rate[c];
                            VectorClass val = exp(cof*dad_length)*theta[i];
                            VectorClass val1 = cof*val;
                            lh_cat += val;
                            df_cat += val1;
                            ddf_cat = mul_add(cof, val1, ddf_cat);
                        }
                        lh_ptn = mul_add(cat_prop[c], lh_cat, lh_ptn);
                        df_ptn = mul_add(cat_prop[c], df_cat, df_ptn);
                        ddf_ptn = mul_add(cat_prop[c], ddf_cat, ddf_ptn);
                        theta += nstates;

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
                    VectorClass df_frac = df_ptn * lh_ptn;
                    VectorClass ddf_frac = ddf_ptn * lh_ptn;
                    VectorClass freq;
                    freq.load_a(&ptn_freq[ptn]);
                    VectorClass tmp1 = df_frac * freq;
                    VectorClass tmp2 = ddf_frac * freq;
                    my_df += tmp1;
                    my_ddf += nmul_add(tmp1, df_frac, tmp2);
                } else {
                    // ascertainment bias correction
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                        df_ptn.cutoff(nptn-ptn);
                        ddf_ptn.cutoff(nptn-ptn);
                    }
                    if (horizontal_or(VectorClass().load_a(&buffer_scale_all[ptn]) != 0.0)) {
                        // some entries are rescaled
                        double *lh_ptn_dbl = (double*)&lh_ptn;
                        double *df_ptn_dbl = (double*)&df_ptn;
                        double *ddf_ptn_dbl = (double*)&ddf_ptn;
                        for (i = 0; i < VectorClass::size(); i++)
                            if (buffer_scale_all[ptn+i] != 0.0) {
                                lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                                df_ptn_dbl[i] *= SCALING_THRESHOLD;
                                ddf_ptn_dbl[i] *= SCALING_THRESHOLD;
                            }
                    }

                    vc_prob_const += lh_ptn;
                    vc_df_const += df_ptn;
                    vc_ddf_const += ddf_ptn;
                }
            } // FOR ptn
        #ifdef _OPENMP
        #pragma omp critical
        #endif
            {
                all_df += my_df;
                all_ddf += my_ddf;
                if (isASC) {
                    all_prob_const += vc_prob_const;
                    all_df_const += vc_df_const;
                    all_ddf_const += vc_ddf_const;
                }
            }
        } // else isMixlen()
    } // FOR thread

    // mark buffer as computed
    theta_computed = true;

    if (isMixlen()) {
        // mixed branch length model
        for (i = 0; i < nmixlen; i++) {
            df[i] = horizontal_add(all_dfvec[i]);
            ASSERT(std::isfinite(df[i]) && "Numerical underflow for lh-derivative");
        }
        for (i = 0; i < nmixlen2; i++)
            ddf[i] = horizontal_add(all_ddfvec[i]);
        // NOTE: last entry of df now store log-likelihood!
        df[nmixlen] = horizontal_add(all_lh);
        return;
    }

    // normal joint branch length model
    *df = horizontal_add(all_df);
    *ddf = horizontal_add(all_ddf);
    if (!std::isfinite(*df)) {
        getModel()->writeInfo(cout);
        getRate()->writeInfo(cout);
    }

    if (!SAFE_NUMERIC && !std::isfinite(*df))
        outError("Numerical underflow (lh-derivative). Run again with the safe likelihood kernel via `-safe` option");

//    ASSERT(std::isfinite(*df) && "Numerical underflow for lh-derivative");

	if (isASC) {
        double prob_const = 0.0, df_const = 0.0, ddf_const = 0.0;
        prob_const = horizontal_add(all_prob_const);
        df_const = horizontal_add(all_df_const);
        ddf_const = horizontal_add(all_ddf_const);
        // ascertainment bias correction
        prob_const = 1.0 - prob_const;
        double df_frac = df_const / prob_const;
        double ddf_frac = ddf_const / prob_const;
        int nsites = aln->getNSite();
        *df += nsites * df_frac;
        *ddf += nsites *(ddf_frac + df_frac*df_frac);
    }

    if (!std::isfinite(*df)) {
        cout << "WARNING: Numerical underflow for lh-derivative" << endl;
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
double PhyloTree::computeLikelihoodBranchSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
double PhyloTree::computeLikelihoodBranchGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#endif
{
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, false);
#endif
//    if ((dad_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(dad_branch, dad);
//    if ((node_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(node_branch, node);
//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    double tree_lh = 0.0;
#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    size_t tip_mem_size = max_orig_nptn * nstates;
    bool isASC = model_factory->unobserved_ptns.size() > 0;

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;

    double *eval = model->getEigenvalues();
    ASSERT(eval);

//    double *val = aligned_alloc<double>(block);
    double *val = NULL;
    double *buffer_partial_lh_ptr = buffer_partial_lh;


    double cat_length[ncat];
    double cat_prop[ncat];
    if (SITE_MODEL) {
        for (c = 0; c < ncat; c++) {
            cat_length[c] = site_rate->getRate(c) * dad_branch->length;
            cat_prop[c] = site_rate->getProp(c);
        }
    } else {
        val = buffer_partial_lh_ptr;
        buffer_partial_lh_ptr += get_safe_upper_limit(block);
        if (nstates % VectorClass::size() == 0) {
            size_t loop_size = nstates / VectorClass::size();
            for (c = 0; c < ncat_mix; c++) {
                size_t mycat = c%ncat;
                size_t m = c/denom;
                mix_addr_nstates[c] = m*nstates;
                mix_addr[c] = mix_addr_nstates[c]*nstates;
                VectorClass *eval_ptr = (VectorClass*)(eval + mix_addr_nstates[c]);
                double len = site_rate->getRate(mycat)*dad_branch->getLength(mycat);
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                VectorClass *this_val = (VectorClass*)(val + c*nstates);
                for (i = 0; i < loop_size; i++)
                    this_val[i] = exp(eval_ptr[i]*len) * prop;
            }
        } else {
            for (c = 0; c < ncat_mix; c++) {
                size_t mycat = c%ncat;
                size_t m = c/denom;
                mix_addr_nstates[c] = m*nstates;
                mix_addr[c] = mix_addr_nstates[c]*nstates;
                double *eval_ptr = eval + mix_addr_nstates[c];
                double len = site_rate->getRate(mycat)*dad_branch->getLength(mycat);
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                double *this_val = val + c*nstates;
                for (i = 0; i < nstates; i++)
                    this_val[i] = exp(eval_ptr[i]*len) * prop;
            }
        }
    }

    VectorClass all_tree_lh(0.0);
    VectorClass all_prob_const(0.0);

    vector<size_t> limits;
    computeBounds<VectorClass>(num_threads, nptn, limits);

    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
//    	double *partial_lh_node = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block);
        double *partial_lh_node;
        if (SITE_MODEL)
            partial_lh_node = &tip_partial_lh[dad->id * tip_mem_size];
        else {
            partial_lh_node = buffer_partial_lh_ptr;
            buffer_partial_lh_ptr += get_safe_upper_limit((aln->STATE_UNKNOWN+1)*block);
        }

        if (!SITE_MODEL) {
            IntVector states_dad = aln->seq_states[dad->id];
            states_dad.push_back(aln->STATE_UNKNOWN);
            // precompute information from one tip
            if (nstates % VectorClass::size() == 0) {
                // vectorized version
                for (IntVector::iterator it = states_dad.begin(); it != states_dad.end(); it++) {
                    double *lh_node = partial_lh_node + (*it)*block;
                    double *lh_tip = tip_partial_lh + (*it)*tip_block;
                    double *vc_val_tmp = val;
                    for (c = 0; c < ncat_mix; c++) {
                        double *this_lh_tip = lh_tip + mix_addr_nstates[c];
                        for (i = 0; i < nstates; i+=VectorClass::size()) {
                            (VectorClass().load_a(&vc_val_tmp[i]) * VectorClass().load_a(&this_lh_tip[i])).store_a(&lh_node[i]);
                        }
                        lh_node += nstates;
                        vc_val_tmp += nstates;
                    }
                }
            } else {
                // non-vectorized version
                for (IntVector::iterator it = states_dad.begin(); it != states_dad.end(); it++) {
                    double *lh_node = partial_lh_node +(*it)*block;
                    double *val_tmp = val;
                    double *this_tip_partial_lh = tip_partial_lh + (*it)*tip_block;
                    for (c = 0; c < ncat_mix; c++) {
                        double *lh_tip = this_tip_partial_lh + mix_addr_nstates[c];
                        for (i = 0; i < nstates; i++) {
                              lh_node[i] = val_tmp[i] * lh_tip[i];
                        }
                        lh_node += nstates;
                        val_tmp += nstates;
                    }
                }
            }
        }

    	// now do the real computation
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static, 1) num_threads(num_threads)
#endif
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {

            VectorClass vc_tree_lh(0.0), vc_prob_const(0.0);

            size_t ptn_lower = limits[thread_id];
            size_t ptn_upper = limits[thread_id+1];

            // reset memory for _pattern_lh_cat
            memset(_pattern_lh_cat + ptn_lower*ncat_mix, 0, sizeof(double)*(ptn_upper-ptn_lower)*ncat_mix);

            // first compute partial_lh
            for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
                computePartialLikelihood(*it, ptn_lower, ptn_upper, thread_id);

            double *vec_tip = buffer_partial_lh_ptr + block*VectorClass::size()*thread_id;

            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn(0.0);
//                lh_ptn.load_a(&ptn_invar[ptn]);
                VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
                VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
                VectorClass *lh_node = SITE_MODEL ? (VectorClass*)&partial_lh_node[ptn*nstates] : (VectorClass*)vec_tip;

                if (SITE_MODEL) {
                    // site-specific model
                    VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                    for (c = 0; c < ncat; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProductExp<VectorClass, double, nstates, FMA>(eval_ptr, lh_node, partial_lh_dad, cat_length[c], lh_cat[c]);
    #else
                        dotProductExp<VectorClass, double, FMA>(eval_ptr, lh_node, partial_lh_dad, cat_length[c], lh_cat[c], nstates);
    #endif
                        if (SAFE_NUMERIC)
                            lh_cat[c] *= cat_prop[c];
                        else
                            lh_ptn += (lh_cat[c] *= cat_prop[c]);

                        partial_lh_dad += nstates;
                    }
                } else { // normal model
                    //load tip vector
                    for (i = 0; i < VectorClass::size(); i++) {
                        double *lh_tip;
                        if (ptn+i < orig_nptn)
                            lh_tip = partial_lh_node + block*(aln->at(ptn+i))[dad->id];
                        else if (ptn+i < max_orig_nptn)
                            lh_tip = partial_lh_node + block*aln->STATE_UNKNOWN;
                        else if (ptn+i < nptn)
                            lh_tip = partial_lh_node + block*model_factory->unobserved_ptns[ptn+i-max_orig_nptn];
                        else
                            lh_tip = partial_lh_node + block*aln->STATE_UNKNOWN;

                        double *this_vec_tip = vec_tip+i;
                        for (c = 0; c < block; c++) {
                            *this_vec_tip = lh_tip[c];
                            this_vec_tip += VectorClass::size();
                        }

                    }
                    // compute likelihood per category
                    for (c = 0; c < ncat_mix; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProductVec<VectorClass, VectorClass, nstates, FMA>(lh_node, partial_lh_dad, lh_cat[c]);
    #else
                        dotProductVec<VectorClass, VectorClass, FMA>(lh_node, partial_lh_dad, lh_cat[c], nstates);
    #endif
                        if (!SAFE_NUMERIC)
                            lh_ptn += lh_cat[c];
                        lh_node += nstates;
                        partial_lh_dad += nstates;
                    }
                } // if SITE_MODEL

                // compute scaling factor per pattern
                VectorClass vc_min_scale(0.0);
                double* vc_min_scale_ptr = (double*)&vc_min_scale;
                if (SAFE_NUMERIC) {
                    // numerical scaling per category
                    UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
                    UBYTE min_scale;
                    for (i = 0; i < VectorClass::size(); i++) {
    //                    scale_dad = dad_branch->scale_num+(ptn+i)*ncat_mix;
                        min_scale = scale_dad[0];
                        for (c = 1; c < ncat_mix; c++)
                            min_scale = min(min_scale, scale_dad[c]);

                        vc_min_scale_ptr[i] = min_scale;

                        double *this_lh_cat = &_pattern_lh_cat[ptn*ncat_mix + i];
                        for (c = 0; c < ncat_mix; c++) {
                            // rescale lh_cat if neccessary
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
                    for (i = 0; i < VectorClass::size(); i++) {
                        vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i];
                    }
                }
                vc_min_scale *= LOG_SCALING_THRESHOLD;

                // Sum later to avoid underflow of invariant sites
                lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
                if (ptn < orig_nptn) {
                    lh_ptn = log(lh_ptn) + vc_min_scale;
                    lh_ptn.store_a(&_pattern_lh[ptn]);
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
                        for (i = 0; i < VectorClass::size(); i++)
                            if (vc_min_scale_ptr[i] != 0.0)
                                lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                    }
                    vc_prob_const += lh_ptn;
                }
            } // FOR PTN
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                all_tree_lh += vc_tree_lh;
                if (isASC)
                    all_prob_const += vc_prob_const;
            }
        } // FOR thread

    } else {

//        ASSERT(0 && "Don't compute tree log-likelihood from internal branch!");
    	//-------- both dad and node are internal nodes -----------/

#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static, 1) num_threads(num_threads)
#endif
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {

            size_t ptn_lower = limits[thread_id];
            size_t ptn_upper = limits[thread_id+1];

            VectorClass vc_tree_lh(0.0), vc_prob_const(0.0);

            // reset memory for _pattern_lh_cat
            memset(_pattern_lh_cat + ptn_lower*ncat_mix, 0, sizeof(double)*(ptn_upper-ptn_lower)*ncat_mix);

            // first compute partial_lh
            for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
                computePartialLikelihood(*it, ptn_lower, ptn_upper, thread_id);

            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn(0.0);
//                lh_ptn.load_a(&ptn_invar[ptn]);
                VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
                VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
                VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);

                // compute likelihood per category
                if (SITE_MODEL) {
                    VectorClass* eval_ptr = (VectorClass*) &eval[ptn*nstates];
                    for (c = 0; c < ncat; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProductExp<VectorClass, double, nstates, FMA>(eval_ptr, partial_lh_node, partial_lh_dad, cat_length[c], lh_cat[c]);
    #else
                        dotProductExp<VectorClass, double, FMA>(eval_ptr, partial_lh_node, partial_lh_dad, cat_length[c], lh_cat[c], nstates);
    #endif
                        if (SAFE_NUMERIC)
                            lh_cat[c] *= cat_prop[c];
                        else
                            lh_ptn += (lh_cat[c] *= cat_prop[c]);
                        partial_lh_node += nstates;
                        partial_lh_dad += nstates;
                    }
                } else {
                    double *val_tmp = val;
                    for (c = 0; c < ncat_mix; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProduct3Vec<VectorClass, double, nstates, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c]);
    #else
                        dotProduct3Vec<VectorClass, double, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c], nstates);
    #endif
                        if (!SAFE_NUMERIC)
                            lh_ptn += lh_cat[c];
                        partial_lh_node += nstates;
                        partial_lh_dad += nstates;
                        val_tmp += nstates;
                    }
                } // if SITE MODEL


                // compute the scaling factor per pattern
                VectorClass vc_min_scale(0.0);
                double* vc_min_scale_ptr = (double*)&vc_min_scale;
                if (SAFE_NUMERIC) {
                    UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
                    UBYTE *scale_node = node_branch->scale_num + ptn*ncat_mix;
                    UBYTE sum_scale[ncat_mix];
                    UBYTE min_scale;

                    for (i = 0; i < VectorClass::size(); i++) {
                        min_scale = sum_scale[0] = scale_dad[0] + scale_node[0];
                        for (c = 1; c < ncat_mix; c++) {
                            sum_scale[c] = scale_dad[c] + scale_node[c];
                            min_scale = min(min_scale, sum_scale[c]);
                        }
                        vc_min_scale_ptr[i] = min_scale;
                        double *this_lh_cat = &_pattern_lh_cat[ptn*ncat_mix + i];
                        for (c = 0; c < ncat_mix; c++) {
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
                    for (i = 0; i < VectorClass::size(); i++) {
                        vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i] + node_branch->scale_num[ptn+i];
                    }
                } // if SAFE_NUMERIC
                vc_min_scale *= LOG_SCALING_THRESHOLD;

                // Sum later to avoid underflow of invariant sites
                lh_ptn = abs(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);

                if (ptn < orig_nptn) {
                    lh_ptn = log(lh_ptn) + vc_min_scale;
                    lh_ptn.store_a(&_pattern_lh[ptn]);
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
                        for (i = 0; i < VectorClass::size(); i++)
                            if (vc_min_scale_ptr[i] != 0.0)
                                lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                    }
                    vc_prob_const += lh_ptn;
                }
            } // FOR LOOP ptn
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                all_tree_lh += vc_tree_lh;
                if (isASC)
                    all_prob_const += vc_prob_const;
            }
        } // FOR thread
    } // else

    tree_lh += horizontal_add(all_tree_lh);

    if (!SAFE_NUMERIC && !std::isfinite(tree_lh))
        outError("Numerical underflow (lh-branch). Run again with the safe likelihood kernel via `-safe` option");

    if (!std::isfinite(tree_lh))
        outWarning("Numerical underflow for lh-branch");

    // arbitrarily fix tree_lh if underflown for some sites
    if (!std::isfinite(tree_lh)) {
        tree_lh = 0.0;
        for (ptn = 0; ptn < orig_nptn; ptn++) {
          if (!std::isfinite(_pattern_lh[ptn])) {
                _pattern_lh[ptn] = LOG_SCALING_THRESHOLD*4; // log(2^(-1024))
            }
            tree_lh += _pattern_lh[ptn] * ptn_freq[ptn];
        }
    }

    if (isASC) {
    	// ascertainment bias correction
        double prob_const = horizontal_add(all_prob_const);
        if (prob_const >= 1.0 || prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        ASSERT(prob_const < 1.0 && prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
//        double inv_const = 1.0 / (1.0-prob_const);
//        size_t nptn_cat = orig_nptn*ncat;
//    	for (ptn = 0; ptn < nptn_cat; ptn++)
//            _pattern_lh_cat[ptn] *= inv_const;
        
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size())
            (VectorClass().load_a(&_pattern_lh[ptn])-prob_const).store_a(&_pattern_lh[ptn]);
//    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
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
double PhyloTree::computeLikelihoodFromBufferSIMD()
#else
template <class VectorClass, const bool FMA, const bool SITE_MODEL>
double PhyloTree::computeLikelihoodFromBufferGenericSIMD()
#endif
{

	ASSERT(theta_all && theta_computed);

//	double tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
//    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    bool isASC = model_factory->unobserved_ptns.size() > 0;

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
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
        for (c = 0; c < ncat; c++) {
            cat_length[c] = site_rate->getRate(c) * current_it->length;
            cat_prop[c] = site_rate->getProp(c);
        }
    } else {
        val0 = buffer_partial_lh;
        if (nstates % VectorClass::size() == 0) {
            VectorClass *vc_val0 = (VectorClass*)val0;
            size_t loop_size = nstates / VectorClass::size();
            for (c = 0; c < ncat_mix; c++) {
                size_t m = c/denom;
                VectorClass *eval_ptr = (VectorClass*)(eval + mix_addr_nstates[c]);
                size_t mycat = c%ncat;
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                double len = site_rate->getRate(mycat) * current_it->getLength(mycat);
                for (i = 0; i < loop_size; i++) {
                    vc_val0[i] = exp(eval_ptr[i] * len) * prop;
                }
                vc_val0 += loop_size;
            }
        } else {
            for (c = 0; c < ncat_mix; c++) {
                size_t m = c/denom;
                double *eval_ptr = eval + mix_addr_nstates[c];
                size_t mycat = c%ncat;
                double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
                size_t addr = c*nstates;
                for (i = 0; i < nstates; i++) {
                    double cof = eval_ptr[i]*site_rate->getRate(mycat);
                    double val = exp(cof*current_it->getLength(mycat)) * prop;
                    val0[addr+i] = val;
                }
            }
        }
    }

//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;

    VectorClass all_tree_lh(0.0), all_prob_const(0.0);

#ifdef _OPENMP
#pragma omp parallel private(ptn, i, c) num_threads(num_threads)
    {
#endif
        VectorClass vc_tree_lh(0.0), vc_prob_const(0.0);
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
    for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
		VectorClass lh_ptn(0.0);
		VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
        if (SITE_MODEL) {
            VectorClass *eval_ptr = (VectorClass*)&eval[ptn*nstates];
//            lh_ptn.load_a(&ptn_invar[ptn]);
            for (c = 0; c < ncat; c++) {
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

        if (ptn < orig_nptn) {
            lh_ptn = log(abs(lh_ptn)) + VectorClass().load_a(&buffer_scale_all[ptn]);
            lh_ptn.store_a(&_pattern_lh[ptn]);
            vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
        } else {
            // bugfix 2016-01-21, prob_const can be rescaled
//                if (min_scale >= 1)
//                    lh_ptn *= SCALING_THRESHOLD;
//				_pattern_lh[ptn] = lh_ptn;
			// ascertainment bias correction
            if (ptn+VectorClass::size() > nptn) {
                // cutoff the last entries if going beyond
                lh_ptn.cutoff(nptn-ptn);
            }
            if (horizontal_or(VectorClass().load_a(&buffer_scale_all[ptn]) != 0.0)) {
                // some entries are rescaled
                double *lh_ptn_dbl = (double*)&lh_ptn;
                for (i = 0; i < VectorClass::size(); i++)
                    if (buffer_scale_all[ptn+i] != 0.0)
                        lh_ptn_dbl[i] *= SCALING_THRESHOLD;
            }
            vc_prob_const += lh_ptn;
        }
    }
#ifdef _OPENMP
#pragma omp critical
        {
            all_tree_lh += vc_tree_lh;
            if (isASC)
                all_prob_const += vc_prob_const;
        }
    }
#else
    all_tree_lh = vc_tree_lh;
    all_prob_const = vc_prob_const;
#endif

    double tree_lh = horizontal_add(all_tree_lh);

    if (!safe_numeric && !std::isfinite(tree_lh))
        outError("Numerical underflow (lh-from-buffer). Run again with the safe likelihood kernel via `-safe` option");

    ASSERT(std::isfinite(tree_lh) && "Numerical underflow for lh-from-buffer");

    // arbitrarily fix tree_lh if underflown for some sites
    if (!std::isfinite(tree_lh)) {
        tree_lh = 0.0;
        for (ptn = 0; ptn < orig_nptn; ptn++) {
            if (!std::isfinite(_pattern_lh[ptn])) {
                _pattern_lh[ptn] = LOG_SCALING_THRESHOLD*4; // log(2^(-1024))
            }
            tree_lh += _pattern_lh[ptn] * ptn_freq[ptn];
        }
    }

    if (isASC) {
    	// ascertainment bias correction
        double prob_const = horizontal_add(all_prob_const);
        if (prob_const >= 1.0 || prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        ASSERT(prob_const < 1.0 && prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
//        double inv_const = 1.0 / (1.0-prob_const);
//        size_t nptn_cat = orig_nptn*ncat;
//    	for (ptn = 0; ptn < nptn_cat; ptn++)
//            _pattern_lh_cat[ptn] *= inv_const;
        
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size())
            (VectorClass().load_a(&_pattern_lh[ptn])-prob_const).store_a(&_pattern_lh[ptn]);
//    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
		ASSERT(std::isfinite(tree_lh));
    }

    return tree_lh;
}


#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervMixlenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA, const bool SITE_MODEL>
void PhyloTree::computeLikelihoodDervMixlenGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf)
#endif
{
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, false);
#endif

//
//    if ((dad_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(dad_branch, dad);
//    if ((node_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(node_branch, node);

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t nmix = (model_factory->fused_mix_rate) ? 1 : model->getNMixtures();

    size_t block = ncat_mix * nstates;
//    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    bool isASC = model_factory->unobserved_ptns.size() > 0;



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

    double *buffer_partial_lh_ptr = buffer_partial_lh;
    vector<size_t> limits;
    computeBounds<VectorClass>(num_threads, nptn, limits);

	ASSERT(theta_all);

    double *val0 = NULL;
    double *val1 = NULL;
    double *val2 = NULL;
    double cat_rate[ncat];
    double cat_prop[ncat];

    int cur_mixlen = getCurMixture();


    if (SITE_MODEL) {
        for (c = 0; c < ncat; c++) {
            cat_rate[c] = site_rate->getRate(c);
            cat_prop[c] = site_rate->getProp(c);
        }
    } else {
        val0 = buffer_partial_lh_ptr;
        val1 = val0 + get_safe_upper_limit(block);
        val2 = val1 + get_safe_upper_limit(block);
        buffer_partial_lh_ptr += 3*get_safe_upper_limit(block);
        double len = dad_branch->getLength(cur_mixlen);
        for (c = 0; c < nmix; c++) {
            size_t cur_mix = (model_factory->fused_mix_rate) ? cur_mixlen : c;
            double *eval_ptr = eval+cur_mix*nstates;
            double prop = model->getMixtureWeight(cur_mix);
            size_t addr = c*nstates;
            for (i = 0; i < nstates; i++) {
                double cof = eval_ptr[i];
                double val = exp(cof*len) * prop;
                double val1_ = cof*val;
                val0[addr+i] = val;
                val1[addr+i] = val1_;
                val2[addr+i] = cof*val1_;
            }
        }
    }

//    double dad_length = dad_branch->length;

    VectorClass all_df(0.0), all_ddf(0.0), all_prob_const(0.0), all_df_const(0.0), all_ddf_const(0.0);

//    size_t nmixlen = getMixlen(), nmixlen2 = nmixlen*nmixlen;
//    ASSERT(nmixlen == ncat);

//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) private(ptn, i, c) num_threads(num_threads)
#endif
    for (int thread_id = 0; thread_id < num_threads; thread_id++) {
        VectorClass my_df(0.0), my_ddf(0.0), vc_prob_const(0.0), vc_df_const(0.0), vc_ddf_const(0.0);
        size_t ptn_lower = limits[thread_id];
        size_t ptn_upper = limits[thread_id+1];

        if (!theta_computed)
        #ifdef KERNEL_FIX_STATES
            computeLikelihoodBufferSIMD<VectorClass, SAFE_NUMERIC, nstates, FMA, SITE_MODEL>(dad_branch, dad, ptn_lower, ptn_upper, thread_id);
        #else
            computeLikelihoodBufferGenericSIMD<VectorClass, SAFE_NUMERIC, FMA, SITE_MODEL>(dad_branch, dad, ptn_lower, ptn_upper, thread_id);
        #endif

        // mixed branch length model
        VectorClass lh_ptn;
        VectorClass df_ptn, ddf_ptn;

        for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            lh_ptn = df_ptn = ddf_ptn = 0.0;
            VectorClass *theta = ((VectorClass*)(theta_all + ptn*block)) + cur_mixlen*nstates*nmix;
            double *val0_ptr = val0;
            double *val1_ptr = val1;
            double *val2_ptr = val2;
            for (c = 0; c < nmix; c++) {
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
//                ASSERT(horizontal_and(lh_ptn > 0.0));

            if (ptn < orig_nptn) {
                VectorClass freq;
                freq.load_a(&ptn_freq[ptn]);
                VectorClass inv_lh_ptn = 1.0 / lh_ptn;

                // compute gradient (my_df)
                df_ptn *= inv_lh_ptn;
                ddf_ptn *= inv_lh_ptn;
                my_df = mul_add(df_ptn, freq, my_df);
                my_ddf += nmul_add(df_ptn, df_ptn, ddf_ptn) * freq;
            } else {
                vc_prob_const += lh_ptn;
                vc_df_const += df_ptn;
                vc_ddf_const += ddf_ptn;
            }
        } // FOR ptn

    #ifdef _OPENMP
    #pragma omp critical
    #endif
        {
            all_df += my_df;
            all_ddf += my_ddf;
            if (isASC) {
                all_prob_const += vc_prob_const;
                all_df_const += vc_df_const;
                all_ddf_const += vc_ddf_const;
            }
        }

    } // FOR thread

    // mark buffer as computed
    theta_computed = true;

    df = horizontal_add(all_df);
    ddf = horizontal_add(all_ddf);

    if (!SAFE_NUMERIC && !std::isfinite(df))
        outError("Numerical underflow (lh-derivative-mixlen). Run again with the safe likelihood kernel via `-safe` option");

	if (isASC) {
        double prob_const = 0.0, df_const = 0.0, ddf_const = 0.0;
        prob_const = 1.0/(1.0 - horizontal_add(all_prob_const));
        df_const = horizontal_add(all_df_const);
        ddf_const = horizontal_add(all_ddf_const);
        // ascertainment bias correction
        df_const *= prob_const;
        ddf_const *= prob_const;
        double nsites = aln->getNSite();
        df += nsites * df_const;
        ddf += nsites * (ddf_const + df_const*df_const);
    }

    if (!std::isfinite(df)) {
        cout << "WARNING: Numerical underflow for lh-derivative-mixlen" << endl;
        df = ddf = 0.0;
    }
}




#endif //PHYLOKERNELNEW_H_
