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
#ifdef KERNEL_FIX_STATES
template <class VectorClass, const size_t nstates, const bool append>
inline void sumVec(VectorClass *A, VectorClass &X, size_t N)
#else
template <class VectorClass, const bool append>
inline void sumVec(VectorClass *A, VectorClass &X, size_t N)
#endif
{
    size_t i;
    if (N % 4 == 0) {
        VectorClass V[4];
        V[0] = A[0];
        V[1] = A[1];
        V[2] = A[2];
        V[3] = A[3];
        for (i = 4; i < N; i+=4) {
            V[0] += A[i];
            V[1] += A[i+1];
            V[2] += A[i+2];
            V[2] += A[i+3];
        }
        if (append)
            X += (V[0] + V[1]) + (V[2] + V[3]);
        else
            X = (V[0] + V[1]) + (V[2] + V[3]);
    } else if (N % 2 == 0) {
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
    } else {
        // odd N
        VectorClass V[2];
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
    }
}

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
    size_t i, j;
    if (N % 4 == 0) {
        VectorClass V[4];
        for (j = 0; j < 4; j++)
            V[j] = A[j] * B[j];
        for (i = 4; i < N; i+=4) {
            for (j = 0; j < 4; j++)
                V[j] = mul_add(A[i+j], B[i+j], V[j]);
        }
        X = (V[0]+V[1]) + (V[2]+V[3]);
    } else if (N % 2 == 0) {
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j];
        for (i = 2; i < N; i+=2) {
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j], B[i+j], V[j]);
        }
        X = (V[0]+V[1]);
    } else {
        // odd number of states
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j];
        for (i = 2; i < N-1; i+=2) {
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j], B[i+j], V[j]);
        }
        X = mul_add(A[N-1], B[N-1], V[0]+V[1]);
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
    size_t i, j;
    if (N % 4 == 0) {
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
    } else if (N % 2 == 0) {
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
    } else {
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
    size_t i, j, x;

    if (N % 4 == 0) {
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
    } else if (N % 2 == 0){
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
    } else {
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
    size_t i, j, x;

    if (N % 4 == 0) {
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
    } else if (N % 2 == 0){
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
    } else {
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
template <class VectorClass, class Numeric, const size_t nstates, const bool FMA>
inline void dotProductTriple(Numeric *A, Numeric *B, Numeric *C, VectorClass *D,
    VectorClass &X, VectorClass &Y, VectorClass &Z, size_t N)
#else
template <class VectorClass, class Numeric, const bool FMA>
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
        X  = AD[0] + AD[1];
        Y  = BD[0] + BD[1];
        Z  = CD[0] + CD[1];
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
        X  = mul_add(A[N-1], D[N-1], AD[0] + AD[1]);
        Y  = mul_add(B[N-1], D[N-1], BD[0] + BD[1]);
        Z  = mul_add(C[N-1], D[N-1], CD[0] + CD[1]);
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
    if (N % 4 == 0) {
        VectorClass V[4];
        for (j = 0; j < 4; j++)
            V[j] = A[j] * B[j] * C[j];
        for (i = 4; i < N; i+=4)
            for (j = 0; j < 4; j++)
                V[j] = mul_add(A[i+j]*B[i+j], C[i+j], V[j]);
        X = (V[0]+V[1])+(V[2]+V[3]);
    } else if (N % 2 == 0) {
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j] * C[j];
        for (i = 2; i < N; i+=2)
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j]*B[i+j], C[i+j], V[j]);
        X = (V[0]+V[1]);
    } else {
        // odd states
        VectorClass V[2];
        for (j = 0; j < 2; j++)
            V[j] = A[j] * B[j] * C[j];
        for (i = 2; i < N-1; i+=2)
            for (j = 0; j < 2; j++)
                V[j] = mul_add(A[i+j]*B[i+j], C[i+j], V[j]);
        X = mul_add(A[N-1]*B[N-1], C[N-1], V[0]+V[1]);
    }
}

/*******************************************************
 *
 * NEW! highly-vectorized partial likelihood function
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA>
void PhyloTree::computePartialLikelihoodSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA>
void PhyloTree::computePartialLikelihoodGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#endif
{

    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    PhyloNode *node = (PhyloNode*)(dad_branch->node);


#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();

    if (!tip_partial_lh_computed)
        computeTipPartialLikelihood();

	if (node->isLeaf()) {
	    dad_branch->lh_scale_factor = 0.0;
		return;
	}
    
    size_t ptn, c;
    size_t orig_ntn = aln->size();
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
    size_t tip_block = nstates * model->getNMixtures();
    size_t max_nptn = get_safe_upper_limit(nptn);
    size_t scale_size = SAFE_NUMERIC ? max_nptn * ncat_mix : max_nptn;
    
	double *evec = model->getEigenvectors();
	double *inv_evec = model->getInverseEigenvectors();
//    double *inv_evec_trans = aligned_alloc<double>(tip_block * nstates);
    // transpose inv_evec
//    for (c = 0; c < model->getNMixtures(); c++) {
//        double *inv_evec_ptr = inv_evec + nstates*nstates*c;
//        double *inv_evec_trans_ptr = inv_evec_trans + nstates*nstates*c;
//        for (i = 0; i < nstates; i++)
//            for (x = 0; x < nstates; x++)
//                inv_evec_trans_ptr[i*nstates+x] = inv_evec_ptr[x*nstates+i];
//    }
	assert(inv_evec && evec);
	double *eval = model->getEigenvalues();

    dad_branch->lh_scale_factor = 0.0;

	// internal node
	PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
    int num_leaves = 0;
	FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNeighbor *nei = (PhyloNeighbor*)*it;
		if (!left) left = (PhyloNeighbor*)(*it); else right = (PhyloNeighbor*)(*it);
        if ((nei->partial_lh_computed & 1) == 0)
            computePartialLikelihood(nei, node);
//        dad_branch->lh_scale_factor += nei->lh_scale_factor;
        if (nei->node->isLeaf())
            num_leaves ++;
	}

    if (params->lh_mem_save == LM_PER_NODE && !dad_branch->partial_lh) {
        // re-orient partial_lh
        bool done = false;
        FOR_NEIGHBOR_IT(node, dad, it2) {
            PhyloNeighbor *backnei = ((PhyloNeighbor*)(*it2)->node->findNeighbor(node));
            if (backnei->partial_lh) {
                dad_branch->partial_lh = backnei->partial_lh;
                dad_branch->scale_num = backnei->scale_num;
                backnei->partial_lh = NULL;
                backnei->scale_num = NULL;
                backnei->partial_lh_computed &= ~1; // clear bit
                done = true;
                break;
            }
        }
        assert(done && "partial_lh is not re-oriented");
    }

    // precompute buffer to save times
    double *echildren = aligned_alloc<double>(block*nstates*(node->degree()-1));
    double *partial_lh_leaves = NULL;
    if (num_leaves > 0)
        partial_lh_leaves = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block*num_leaves);
    double *echild = echildren;
    double *partial_lh_leaf = partial_lh_leaves;

    if (nstates % VectorClass::size() == 0) {
        // vectorized version
        VectorClass *expchild = (VectorClass*)aligned_alloc<double>(nstates);
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNeighbor *child = (PhyloNeighbor*)*it;
            VectorClass *echild_ptr = (VectorClass*)echild;
            // precompute information buffer
            for (c = 0; c < ncat_mix; c++) {
                VectorClass len_child = site_rate->getRate(c%ncat) * child->length;
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
        aligned_free(expchild);
    } else {
        // non-vectorized version
        double expchild[nstates];
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNeighbor *child = (PhyloNeighbor*)*it;
            // precompute information buffer
            double *echild_ptr = echild;
            for (c = 0; c < ncat_mix; c++) {
                double len_child = site_rate->getRate(c%ncat) * child->length;
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


//    double sum_scale = 0.0;

    double *eleft = echildren, *eright = echildren + block*nstates;
    
	if (!left->node->isLeaf() && right->node->isLeaf()) {
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
        double *etmp = eleft;
        eleft = eright;
        eright = etmp;
	}
    
    if (node->degree() > 3) {
        /*--------------------- multifurcating node ------------------*/

        double *partial_lh_all_dbl = aligned_alloc<double>(block*VectorClass::size()*2);
        double *vec_tip = partial_lh_all_dbl + block*VectorClass::size();

        // now for-loop computing partial_lh over all site-patterns
#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
        for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
            VectorClass *partial_lh_all = (VectorClass*)partial_lh_all_dbl;
            for (i = 0; i < block; i++)
                partial_lh_all[i] = 1.0;
            UBYTE *scale_dad = NULL;
            if (SAFE_NUMERIC) {
                scale_dad = dad_branch->scale_num + ptn*ncat_mix;
                memset(scale_dad, 0, sizeof(UBYTE)*ncat_mix*VectorClass::size());
            } else
                memset(&dad_branch->scale_num[ptn], 0, sizeof(UBYTE)*VectorClass::size());

            double *partial_lh_leaf = partial_lh_leaves;
            double *echild = echildren;

            FOR_NEIGHBOR_IT(node, dad, it) {
                PhyloNeighbor *child = (PhyloNeighbor*)*it;
                UBYTE *scale_child = SAFE_NUMERIC ? child->scale_num + ptn*ncat_mix : NULL;
                if (child->node->isLeaf()) {
                    // external node
                    // load data for tip
                    for (i = 0; i < VectorClass::size(); i++) {
                        double *child_lh;
                        if (ptn+i < orig_ntn)
                            child_lh = partial_lh_leaf + block*(aln->at(ptn+i))[child->node->id];
                        else if (ptn+i < nptn)
                            child_lh = partial_lh_leaf + block*model_factory->unobserved_ptns[ptn+i-orig_ntn];
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
            } // FOR_NEIGHBOR
            
        
            // compute dot-product with inv_eigenvector
            VectorClass *partial_lh_tmp = partial_lh_all;
            VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;
            for (c = 0; c < ncat_mix; c++) {
                if (SAFE_NUMERIC)
                    lh_max = 0.0;
                double *inv_evec_ptr = inv_evec + mix_addr[c];
                for (i = 0; i < nstates; i++) {
					VectorClass res = partial_lh_tmp[0]*inv_evec_ptr[0];
					for (x = 1; x < nstates; x++) {
						res = mul_add(partial_lh_tmp[x], inv_evec_ptr[x], res);
					}
                    inv_evec_ptr += nstates;
                    partial_lh[i] = res;
                    res = abs(res);
                    lh_max = max(lh_max, res);
                }
                // check if one should scale partial likelihoods
                if (SAFE_NUMERIC) {
                    auto underflown = ((lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                    if (horizontal_or(underflown)) { // at least one site has numerical underflown
                        for (x = 0; x < VectorClass::size(); x++)
                        if (underflown[x]) {
                            // BQM 2016-05-03: only scale for non-constant sites
                            // now do the likelihood scaling
                            double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                            for (i = 0; i < nstates; i++)
                                partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
                            dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
                        }
                    }
                }
                partial_lh += nstates;
                partial_lh_tmp += nstates;
            }

            if (!SAFE_NUMERIC) {
                auto underflown = (lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                if (horizontal_or(underflown)) { // at least one site has numerical underflown
                    for (x = 0; x < VectorClass::size(); x++)
                    if (underflown[x]) {
                        double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                        // now do the likelihood scaling
                        for (i = 0; i < block; i++) {
                            partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
                        }
//                        sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                        dad_branch->scale_num[ptn+x] += 1;
                    }
                }
            }

        } // for ptn
//        dad_branch->lh_scale_factor += sum_scale;               
        aligned_free(partial_lh_all_dbl);

        // end multifurcating treatment
    } else if (left->node->isLeaf() && right->node->isLeaf()) {

        /*--------------------- TIP-TIP (cherry) case ------------------*/

        double *partial_lh_left = partial_lh_leaves;
        double *partial_lh_right = partial_lh_leaves + (aln->STATE_UNKNOWN+1)*block;
        // TODO: this is not thread-safe
        double *vec_left = aligned_alloc<double>(block*VectorClass::size()*2 + nstates*VectorClass::size());
        double *vec_right = vec_left + block*VectorClass::size();
        double *partial_lh_dbl = vec_right + block*VectorClass::size();

		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, scale_size * sizeof(UBYTE));
#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass* partial_lh_tmp = (VectorClass*)partial_lh_dbl;
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass *vleft = (VectorClass*)vec_left;
            VectorClass *vright = (VectorClass*)vec_right;

            // load data for tip
            for (x = 0; x < VectorClass::size(); x++) {
                double *tip_left, *tip_right;
                if (ptn+x < orig_ntn) {
                    tip_left  = partial_lh_left  + block * (aln->at(ptn+x))[left->node->id];
                    tip_right = partial_lh_right + block * (aln->at(ptn+x))[right->node->id];
                } else if (ptn+x < nptn) {
                    tip_left  = partial_lh_left  + block * model_factory->unobserved_ptns[ptn+x-orig_ntn];
                    tip_right = partial_lh_right + block * model_factory->unobserved_ptns[ptn+x-orig_ntn];
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
			}
		}

        aligned_free(vec_left);

	} else if (left->node->isLeaf() && !right->node->isLeaf()) {

        /*--------------------- TIP-INTERNAL NODE case ------------------*/

		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, scale_size * sizeof(UBYTE));


        double *partial_lh_left = partial_lh_leaves;

        // TODO: this is not thread-safe
        double *vec_left = aligned_alloc<double>(block*VectorClass::size() + nstates*VectorClass::size());
        double *partial_lh_dbl = vec_left + block*VectorClass::size();


#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass *partial_lh_tmp = (VectorClass*)partial_lh_dbl;
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            memset(partial_lh, 0, sizeof(VectorClass)*block);
            VectorClass *vleft = (VectorClass*)vec_left;
            // load data for tip
            for (x = 0; x < VectorClass::size(); x++) {
                double *tip;
                if (ptn+x < orig_ntn) {
                    tip = partial_lh_left + block*(aln->at(ptn+x))[left->node->id];
                } else if (ptn+x < nptn) {
                    tip = partial_lh_left + block*model_factory->unobserved_ptns[ptn+x-orig_ntn];
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
            VectorClass lh_max = 0.0;
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
                    auto underflown = ((lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                    if (horizontal_or(underflown)) { // at least one site has numerical underflown
                        for (x = 0; x < VectorClass::size(); x++)
                        if (underflown[x]) {
                            // BQM 2016-05-03: only scale for non-constant sites
                            // now do the likelihood scaling
                            double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                            for (i = 0; i < nstates; i++)
                                partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
                            dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
                        }
                    }
                }
                vleft += nstates;
                partial_lh_right += nstates;
                partial_lh += nstates;
			}

            if (!SAFE_NUMERIC) {
                auto underflown = (lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                if (horizontal_or(underflown)) { // at least one site has numerical underflown
                    for (x = 0; x < VectorClass::size(); x++)
                    if (underflown[x]) {
                        double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                        // now do the likelihood scaling
                        for (i = 0; i < block; i++) {
                            partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
                        }
//                        sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                        dad_branch->scale_num[ptn+x] += 1;
                    }
                }
            }

		} // big for loop over ptn
        aligned_free(vec_left);
//		dad_branch->lh_scale_factor += sum_scale;
//		delete [] partial_lh_left;

	} else {

        /*--------------------- INTERNAL-INTERNAL NODE case ------------------*/
        double *partial_lh_dbl = aligned_alloc<double>(nstates*VectorClass::size());

#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass *partial_lh_tmp = (VectorClass*)partial_lh_dbl;
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_left = (VectorClass*)(left->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
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

            VectorClass lh_max = 0.0;
			for (c = 0; c < ncat_mix; c++) {
                if (SAFE_NUMERIC) {
                    lh_max = 0.0;
                    for (x = 0; x < VectorClass::size(); x++)
                        scale_dad[x*ncat_mix] = scale_left[x*ncat_mix] + scale_right[x*ncat_mix];
                }
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
                // check if one should scale partial likelihoods
                if (SAFE_NUMERIC) {
                    auto underflown = ((lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0));
                    if (horizontal_or(underflown))
                        for (x = 0; x < VectorClass::size(); x++)
                        if (underflown[x]) {
                            // BQM 2016-05-03: only scale for non-constant sites
                            // now do the likelihood scaling
                            double *partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
                            for (i = 0; i < nstates; i++)
                                partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
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
                auto underflown = (lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(&ptn_invar[ptn]) == 0.0);
                if (horizontal_or(underflown)) { // at least one site has numerical underflown
                    for (x = 0; x < VectorClass::size(); x++)
                    if (underflown[x]) {
                        double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                        // now do the likelihood scaling
                        for (i = 0; i < block; i++) {
                            partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
                        }
//                        sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn+x];
                        dad_branch->scale_num[ptn+x] += 1;
                    }
                }
            }

		} // big for loop over ptn
//		dad_branch->lh_scale_factor += sum_scale;

        aligned_free(partial_lh_dbl);

	}

    if (partial_lh_leaves)
        aligned_free(partial_lh_leaves);
    aligned_free(echildren);
//    aligned_free(inv_evec_trans);
}

/*******************************************************
 *
 * NEW! highly-vectorized log-likelihood derivative function
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA>
void PhyloTree::computeLikelihoodDervSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA>
void PhyloTree::computeLikelihoodDervGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf)
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
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(node_branch, node);

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
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    size_t maxptn = ((nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    maxptn = max(maxptn, aln->size()+((model_factory->unobserved_ptns.size()+VectorClass::size()-1)/VectorClass::size())*VectorClass::size());

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    double *eval = model->getEigenvalues();
    assert(eval);

	assert(theta_all);
	if (!theta_computed) {
		// precompute theta for fast branch length optimization
	    if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case

            // TODO not thread-safe
            double *vec_tip = aligned_alloc<double>(tip_block*VectorClass::size());

#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static) reduction(+: scale_all)
#endif
	    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
				VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
				VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
                //load tip vector
                for (i = 0; i < VectorClass::size(); i++) {
                    double *this_tip_partial_lh;
                    if (ptn+i < orig_nptn)
                        this_tip_partial_lh = tip_partial_lh + tip_block*(aln->at(ptn+i))[dad->id];
                    else if (ptn+i < nptn)
                        this_tip_partial_lh = tip_partial_lh + tip_block*model_factory->unobserved_ptns[ptn+i-orig_nptn];
                    else
                        this_tip_partial_lh = tip_partial_lh + tip_block*aln->STATE_UNKNOWN;
                    double *this_vec_tip = vec_tip+i;
                    for (c = 0; c < tip_block; c++) {
                        *this_vec_tip = this_tip_partial_lh[c];
                        this_vec_tip += VectorClass::size();
                    }

                }
                for (c = 0; c < ncat_mix; c++) {
                    VectorClass *lh_tip = (VectorClass*)(vec_tip + mix_addr_nstates[c]*VectorClass::size());
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
                    VectorClass *buf = (VectorClass*)(buffer_scale_all+ptn);
                    *buf *= LOG_SCALING_THRESHOLD;
                } else {
                    for (i = 0; i < VectorClass::size(); i++)
                        buffer_scale_all[ptn+i] = dad_branch->scale_num[ptn+i];
                    VectorClass *buf = (VectorClass*)(buffer_scale_all+ptn);
                    *buf *= LOG_SCALING_THRESHOLD;
                }

			}
            aligned_free(vec_tip);
	    } else {
	    	// both dad and node are internal nodes

//	    	size_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static) reduction(+: scale_all)
#endif
	    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
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
                    VectorClass *buf = (VectorClass*)(buffer_scale_all+ptn);
                    *buf *= LOG_SCALING_THRESHOLD;
                } else {
                    for (i = 0; i < VectorClass::size(); i++)
                        buffer_scale_all[ptn+i] = dad_branch->scale_num[ptn+i] + node_branch->scale_num[ptn+i];
                    VectorClass *buf = (VectorClass*)(buffer_scale_all+ptn);
                    *buf *= LOG_SCALING_THRESHOLD;
                }
			}
	    }
        // NO NEED TO copy dummy values!

		theta_computed = true;
	}

    double *val0 = aligned_alloc<double>(block);
    double *val1 = aligned_alloc<double>(block);
    double *val2 = aligned_alloc<double>(block);
	for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        double *eval_ptr = eval + mix_addr_nstates[c];
        size_t mycat = c%ncat;
        double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        size_t addr = c*nstates;
		for (i = 0; i < nstates; i++) {
			double cof = eval_ptr[i]*site_rate->getRate(mycat);
			double val = exp(cof*dad_branch->length) * prop;
			double val1_ = cof*val;
			val0[addr+i] = val;
			val1[addr+i] = val1_;
			val2[addr+i] = cof*val1_;
		}
	}


    VectorClass my_df = 0.0, my_ddf = 0.0, vc_prob_const = 0.0, vc_df_const = 0.0, vc_ddf_const = 0.0;
    double prob_const = 0.0, df_const = 0.0, ddf_const = 0.0;
//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf, prob_const, df_const, ddf_const) private(ptn, i) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
		VectorClass lh_ptn;
        //lh_ptn.load_a(&ptn_invar[ptn]);
		VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
        VectorClass df_ptn, ddf_ptn;

#ifdef KERNEL_FIX_STATES
        dotProductTriple<VectorClass, double, nstates, FMA>(val0, val1, val2, theta, lh_ptn, df_ptn, ddf_ptn, block);
#else
        dotProductTriple<VectorClass, double, FMA>(val0, val1, val2, theta, lh_ptn, df_ptn, ddf_ptn, block, nstates);
#endif
        lh_ptn = abs(lh_ptn + VectorClass().load_a(&ptn_invar[ptn]));
        
        if (ptn < orig_nptn) {
            // TODO ASC overlap
            lh_ptn = 1.0 / lh_ptn;
			VectorClass df_frac = df_ptn * lh_ptn;
			VectorClass ddf_frac = ddf_ptn * lh_ptn;
			VectorClass freq;
            freq.load_a(&ptn_freq[ptn]);
			VectorClass tmp1 = df_frac * freq;
			VectorClass tmp2 = ddf_frac * freq;
			my_df += tmp1;
//			my_ddf += tmp2 - tmp1 * df_frac;
            my_ddf += nmul_add(tmp1, df_frac, tmp2);
		} else {
			// ascertainment bias correction
			vc_prob_const += lh_ptn;
			vc_df_const += df_ptn;
			vc_ddf_const += ddf_ptn;
		}
    }
	df = horizontal_add(my_df);
	ddf = horizontal_add(my_ddf);
    prob_const = horizontal_add(vc_prob_const);
    df_const = horizontal_add(vc_df_const);
    ddf_const = horizontal_add(vc_ddf_const);

    if (!SAFE_NUMERIC && (isnan(df) || isinf(df)))
        outError("Numerical underflow (lh-derivative). Run again with the safe likelihood kernel via `-safe` option");

    assert(!isnan(df) && !isinf(df) && "Numerical underflow for lh-derivative");

	if (orig_nptn < nptn) {
    	// ascertainment bias correction
    	prob_const = 1.0 - prob_const;
    	double df_frac = df_const / prob_const;
    	double ddf_frac = ddf_const / prob_const;
    	int nsites = aln->getNSite();
    	df += nsites * df_frac;
    	ddf += nsites *(ddf_frac + df_frac*df_frac);
    }


    aligned_free(val2);
    aligned_free(val1);
    aligned_free(val0);
}




/*******************************************************
 *
 * NEW! highly-vectorized log-likelihood function
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA>
double PhyloTree::computeLikelihoodBranchSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA>
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
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(node_branch, node);
//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    double tree_lh = 0.0;
    VectorClass vc_tree_lh = 0.0;
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
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;

    double *eval = model->getEigenvalues();
    assert(eval);

    double *val = aligned_alloc<double>(block);
	for (c = 0; c < ncat_mix; c++) {
        size_t mycat = c%ncat;
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
        double *eval_ptr = eval + mix_addr_nstates[c];
		double len = site_rate->getRate(mycat)*dad_branch->length;
		double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        double *this_val = val + c*nstates;
		for (i = 0; i < nstates; i++)
			this_val[i] = exp(eval_ptr[i]*len) * prop;
	}

	double prob_const = 0.0;
    VectorClass vc_prob_const = 0.0;
	memset(_pattern_lh_cat, 0, sizeof(double)*nptn*ncat_mix);

    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
    	double *partial_lh_node = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block);
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

        double *vec_tip = aligned_alloc<double>(block*VectorClass::size());

    	// now do the real computation
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, i, c) schedule(static)
#endif
    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass lh_ptn;
            lh_ptn.load_a(&ptn_invar[ptn]);
            VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
            VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass *lh_node = (VectorClass*)vec_tip;

            //load tip vector
            for (i = 0; i < VectorClass::size(); i++) {
                double *lh_tip;
                if (ptn+i < orig_nptn)
                    lh_tip = partial_lh_node + block*(aln->at(ptn+i))[dad->id];
                else if (ptn+i < nptn)
                    lh_tip = partial_lh_node + block*model_factory->unobserved_ptns[ptn+i-orig_nptn];
                else
                    lh_tip = partial_lh_node + block*aln->STATE_UNKNOWN;

                double *this_vec_tip = vec_tip+i;
                for (c = 0; c < block; c++) {
                    *this_vec_tip = lh_tip[c];
                    this_vec_tip += VectorClass::size();
                }

            }

            VectorClass vc_min_scale(0.0);
            double* vc_min_scale_ptr = (double*)&vc_min_scale;

            if (SAFE_NUMERIC) {
                // numerical scaling per category
                UBYTE *scale_dad;
                UBYTE min_scale;
                for (i = 0; i < VectorClass::size(); i++) {
                    scale_dad = dad_branch->scale_num+(ptn+i)*ncat_mix;
                    min_scale = scale_dad[0];
                    for (c = 1; c < ncat_mix; c++)
                        min_scale = min(min_scale, scale_dad[c]);

                    vc_min_scale_ptr[i] = min_scale;

                    for (c = 0; c < ncat_mix; c++) {
                        if (scale_dad[c] == min_scale+1) {
                            double *this_tip = &vec_tip[c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_tip[x*VectorClass::size()] *= SCALING_THRESHOLD;
                            }
                        } else if (scale_dad[c] > min_scale+1) {
                            double *this_tip = &vec_tip[c*nstates*VectorClass::size() + i];
                            for (size_t x = 0; x < nstates; x++) {
                                this_tip[x*VectorClass::size()] = 0.0;
                            }
                        }
                    }
                }
                vc_min_scale *= LOG_SCALING_THRESHOLD;
            } else {
                for (i = 0; i < VectorClass::size(); i++) {
                    vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i];
                }
                vc_min_scale *= LOG_SCALING_THRESHOLD;
            }

            for (c = 0; c < ncat_mix; c++) {
#ifdef KERNEL_FIX_STATES
                dotProductVec<VectorClass, VectorClass, nstates, FMA>(lh_node, partial_lh_dad, *lh_cat);
#else
                dotProductVec<VectorClass, VectorClass, FMA>(lh_node, partial_lh_dad, *lh_cat, nstates);
#endif
                lh_ptn += *lh_cat;
                lh_node += nstates;
                partial_lh_dad += nstates;
                lh_cat++;
            }

			if (ptn < orig_nptn) {
                lh_ptn = log(abs(lh_ptn)) + vc_min_scale;
				lh_ptn.store_a(&_pattern_lh[ptn]);
				vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
                // TODO: case that ASC overlap
			} else {
                // bugfix 2016-01-21, prob_const can be rescaled
                if (horizontal_or(vc_min_scale != 0.0)) {
                    // some entries are rescaled
                    for (i = 0; i < VectorClass::size(); i++)
                        if (vc_min_scale[i] != 0.0)
                            lh_ptn.insert(i, lh_ptn[i] * SCALING_THRESHOLD);
                }
                vc_prob_const += lh_ptn;
			}
		}
        aligned_free(vec_tip);
		aligned_free(partial_lh_node);
    } else {
    	// both dad and node are internal nodes

#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, i, c) schedule(static)
#endif
    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass lh_ptn;
            lh_ptn.load_a(&ptn_invar[ptn]);
            VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
            VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
            double *val_tmp = val;

            VectorClass vc_min_scale(0.0);
            double* vc_min_scale_ptr = (double*)&vc_min_scale;

            if (SAFE_NUMERIC) {

                for (c = 0; c < ncat_mix; c++) {
#ifdef KERNEL_FIX_STATES
                    dotProduct3Vec<VectorClass, double, nstates, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c]);
#else
                    dotProduct3Vec<VectorClass, double, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c], nstates);
#endif
                    partial_lh_node += nstates;
                    partial_lh_dad += nstates;
                    val_tmp += nstates;
                }

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
                vc_min_scale *= LOG_SCALING_THRESHOLD;
            } else {
                // normal scaling
                for (c = 0; c < ncat_mix; c++) {
#ifdef KERNEL_FIX_STATES
                    dotProduct3Vec<VectorClass, double, nstates, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c]);
#else
                    dotProduct3Vec<VectorClass, double, FMA>(val_tmp, partial_lh_node, partial_lh_dad, lh_cat[c], nstates);
#endif
                    lh_ptn += lh_cat[c];
                    partial_lh_node += nstates;
                    partial_lh_dad += nstates;
                    val_tmp += nstates;
                }
                for (i = 0; i < VectorClass::size(); i++) {
                    vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i] + node_branch->scale_num[ptn+i];
                }
                vc_min_scale *= LOG_SCALING_THRESHOLD;
            }

			if (ptn < orig_nptn) {
                lh_ptn = log(abs(lh_ptn)) + vc_min_scale;
				lh_ptn.store_a(&_pattern_lh[ptn]);
				vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
                // TODO: case that ASC overlap
			} else {
                // bugfix 2016-01-21, prob_const can be rescaled
                if (horizontal_or(vc_min_scale != 0.0)) {
                    // some entries are rescaled
                    for (i = 0; i < VectorClass::size(); i++)
                        if (vc_min_scale[i] != 0.0)
                            lh_ptn.insert(i, lh_ptn[i] * SCALING_THRESHOLD);
                }
                vc_prob_const += lh_ptn;
			}
		}
    }

    tree_lh += horizontal_add(vc_tree_lh);

    if (!SAFE_NUMERIC && (isnan(tree_lh) || isinf(tree_lh)))
        outError("Numerical underflow (lh-branch). Run again with the safe likelihood kernel via `-safe` option");

    assert(!isnan(tree_lh) && !isinf(tree_lh) && "Numerical underflow for lh-branch");

    if (orig_nptn < nptn) {
    	// ascertainment bias correction
        prob_const = horizontal_add(vc_prob_const);
        if (prob_const >= 1.0 || prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        assert(prob_const < 1.0 && prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
//        double inv_const = 1.0 / (1.0-prob_const);
//        size_t nptn_cat = orig_nptn*ncat;
//    	for (ptn = 0; ptn < nptn_cat; ptn++)
//            _pattern_lh_cat[ptn] *= inv_const;
        
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
		assert(!isnan(tree_lh) && !isinf(tree_lh));
    }

    aligned_free(val);
    return tree_lh;
}


/*******************************************************
 *
 * NEW! highly-vectorized log-likelihood from buffer
 *
 ******************************************************/

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA>
double PhyloTree::computeLikelihoodFromBufferSIMD()
#else
template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA>
double PhyloTree::computeLikelihoodFromBufferGenericSIMD()
#endif
{

	assert(theta_all && theta_computed);

//	double tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;
    double tree_lh = 0.0;
    VectorClass vc_tree_lh = 0.0;

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
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    double *eval = model->getEigenvalues();
    assert(eval);

    double *val0 = aligned_alloc<double>(block);
	for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        double *eval_ptr = eval + mix_addr_nstates[c];
        size_t mycat = c%ncat;
        double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        size_t addr = c*nstates;
		for (i = 0; i < nstates; i++) {
			double cof = eval_ptr[i]*site_rate->getRate(mycat);
			double val = exp(cof*current_it->length) * prop;
			val0[addr+i] = val;
		}
	}


    VectorClass vc_prob_const = 0.0;
    double prob_const = 0.0;
//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf, prob_const, df_const, ddf_const) private(ptn, i) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
		VectorClass lh_ptn;
		VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
        dotProductVec<VectorClass, double, FMA>(val0, theta, lh_ptn, block);
        lh_ptn += VectorClass().load_a(&ptn_invar[ptn]);

        if (ptn < orig_nptn) {
            lh_ptn = log(abs(lh_ptn)) + VectorClass().load_a(&buffer_scale_all[ptn]);
            lh_ptn.store_a(&_pattern_lh[ptn]);
            vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
        } else {
            // bugfix 2016-01-21, prob_const can be rescaled
//                if (min_scale >= 1)
//                    lh_ptn *= SCALING_THRESHOLD;
//				_pattern_lh[ptn] = lh_ptn;
            vc_prob_const += lh_ptn;
        }
    }

    tree_lh += horizontal_add(vc_tree_lh);

    if (!SAFE_NUMERIC && (isnan(tree_lh) || isinf(tree_lh)))
        outError("Numerical underflow (lh-from-buffer). Run again with the safe likelihood kernel via `-safe` option");

    assert(!isnan(tree_lh) && !isinf(tree_lh) && "Numerical underflow for lh-from-buffer");

    if (orig_nptn < nptn) {
    	// ascertainment bias correction
        prob_const = horizontal_add(vc_prob_const);
        if (prob_const >= 1.0 || prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        assert(prob_const < 1.0 && prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
//        double inv_const = 1.0 / (1.0-prob_const);
//        size_t nptn_cat = orig_nptn*ncat;
//    	for (ptn = 0; ptn < nptn_cat; ptn++)
//            _pattern_lh_cat[ptn] *= inv_const;
        
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
		assert(!isnan(tree_lh) && !isinf(tree_lh));
    }

    aligned_free(val0);
    return tree_lh;
}


#endif //PHYLOKERNELNEW_H_
