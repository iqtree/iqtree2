#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include <iomanip>
#include <algorithm>
#include <iterator>

#include <boost/tr1/unordered_map.hpp>
#include "axml.h"

#include "aligned_buffer.h"
#include "vec_unit.h"
#include "cycle.h"

namespace sim_reorder {

/* generic function for computing the P matrices, for computing the conditional likelihood at a node p, given child nodes q and r
 we compute P(z1) and P(z2) here */

void makeP(double z1, double z2, double *rptr, double *EI, double *EIGN,
        int numberOfCategories, double *left, double *right, boolean saveMem,
        int maxCat, const int states) {
    int i, j, k,
    /* square of the number of states = P-matrix size */
    statesSquare = states * states;

    /* assign some space for pre-computing and later re-using functions */

    double *lz1 = (double*) malloc(sizeof(double) * states), *lz2 =
            (double*) malloc(sizeof(double) * states), *d1 = (double*) malloc(
                    sizeof(double) * states), *d2 = (double*) malloc(
                            sizeof(double) * states);

    /* multiply branch lengths with eigenvalues */

    for (i = 1; i < states; i++) {
        lz1[i] = EIGN[i] * z1;
        lz2[i] = EIGN[i] * z2;
    }

    /* loop over the number of rate categories, this will be 4 for the GAMMA model and
	 variable for the CAT model */

    for (i = 0; i < numberOfCategories; i++) {
        /* exponentiate the rate multiplied by the branch */

        for (j = 1; j < states; j++) {
            d1[j] = EXP(rptr[i] * lz1[j]);
            d2[j] = EXP(rptr[i] * lz2[j]);
        }

        /* now fill the P matrices for the two branch length values */

        for (j = 0; j < states; j++) {
            /* left and right are pre-allocated arrays */

            left[statesSquare * i + states * j] = 1.0;
            right[statesSquare * i + states * j] = 1.0;

            for (k = 1; k < states; k++) {
                left[statesSquare * i + states * j + k] = d1[k]
                                                             * EI[states * j + k];
                right[statesSquare * i + states * j + k] = d2[k]
                                                              * EI[states * j + k];
            }
        }
    }

    /* if memory saving is enabled and we are using CAT we need to do one additional P matrix
	 calculation for a rate of 1.0 to compute the entries of a column/tree site comprising only gaps */

    if (saveMem) {
        i = maxCat;

        for (j = 1; j < states; j++) {
            d1[j] = EXP (lz1[j]);
            d2[j] = EXP (lz2[j]);
        }

        for (j = 0; j < states; j++) {
            left[statesSquare * i + states * j] = 1.0;
            right[statesSquare * i + states * j] = 1.0;

            for (k = 1; k < states; k++) {
                left[statesSquare * i + states * j + k] = d1[k]
                                                             * EI[states * j + k];
                right[statesSquare * i + states * j + k] = d2[k]
                                                              * EI[states * j + k];
            }
        }
    }

    /* free the temporary buffers */

    free(lz1);
    free(lz2);
    free(d1);
    free(d2);
}

template<typename oiter>
inline void reorder_tip_block(uint8_t *tip_idx, double *tip_vec, size_t block_start,
        const size_t entry_width, oiter ostart, const size_t vw, size_t n) {

//     typedef vector_unit<double,2> vu:
//     typedef vu::vec_t vec_t;
//     
    int i = block_start;
    
    if( i + vw < n ) {
        
//         const int off0 = tip_idx[i] * entry_width;
//         const int off1 = tip_idx[i+1] * entry_width;
//         
//         
//         for (size_t j = 0; j < entry_width; ++j) {
//             //tmp0[j] = tip_vec[off0 + j];
//             //tmp1[j] = tip_vec[off1 + j];
//             
//             *ostart = tip_vec[off0 + j];
//             ++ostart;
//             *ostart = tip_vec[off1 + j];
//             ++ostart;
//             
//         }
//             for (size_t k = 0; k < vw; ++k, ++ostart) {
//                 // do the padding on the fly. most likely very inefficient.
//                 int ti = tip_idx[i+k];
//                 
//                 *ostart = tip_vec[entry_width * ti + j];
//             }
        
        for (size_t j = 0; j < entry_width; ++j) {
            
            for (size_t k = 0; k < vw; ++k, ++ostart) {
                // do the padding on the fly. most likely very inefficient.
                int ti = tip_idx[i+k];
                
                *ostart = tip_vec[entry_width * ti + j];
            }
        }
    } else {
    
        for (size_t j = 0; j < entry_width; ++j) {
            
            for (size_t k = 0; k < vw; ++k, ++ostart) {
                // do the padding on the fly. most likely very inefficient.
                int ti;
                if( i + k < n ) {
                    ti = tip_idx[i+k];
                } else {
                    ti = tip_idx[0];
                }
                *ostart = tip_vec[entry_width * ti + j];
            }
        }
    }
}

template<typename oiter>
void reorder(double *v, size_t width, size_t entry_width, oiter ostart,
        size_t vw) {

    assert( width % vw == 0);

    //    size_t blocks = width % vw;

    for (size_t i = 0; i < width; i += vw) {

        for (size_t j = 0; j < entry_width; ++j) {

            for (size_t k = 0; k < vw; ++k, ++ostart) {
                *ostart = v[(i + k) * entry_width + j];
            }
        }
    }

}

template<typename oiter, typename iiter>
void reorder_back(oiter v, size_t width, size_t states, iiter istart, size_t vw) {

    assert( width % vw == 0);

    //    std::cout << "rob: " << width << " " << states << " " << vw << "\n";
    size_t entry_width = states;

    for (size_t i = 0; i < width; i += vw) {

        for (size_t j = 0; j < entry_width; ++j) {

            for (size_t k = 0; k < vw; ++k, ++istart) {

                v[(i + k) * entry_width + j] = *istart;
            }
        }
    }

}

template<typename oiter>
void reorder_p(oiter o, int* cptr, double *p, size_t states, size_t vw) {
    const size_t s2 = states * states;

    for (size_t i = 0; i < s2; ++i) {
        for (size_t j = 0; j < vw; ++j, ++o) {
            *o = p[cptr[j] * s2 + i];
        }
    }
}

template<typename T>
class delta_equal {
public:
    delta_equal(const T &delta) :
        delta_(delta) {
    }

    inline bool operator()(const T &v1, const T &v2) {
        T diff = v1 - v2;

        return (diff < 0 && diff > -delta_) || (diff >= 0 && diff < delta_);
    }

private:
    const T delta_;
};

template<typename iiter1, typename iiter2>
void print_delta(iiter1 start1, iiter1 end1, iiter2 start2) {

    for (; start1 != end1; ++start1, ++start2) {
        const char *meeeep =
                fabs(*start1 - *start2) > 1e-3 ? "<<<< meeeeeeep" : "";

        std::cout << std::setw(20) << std::left << *start1 << *start2 << meeeep
                << "\n"; // ": " << *start1 - *start2 << "\n";
    }
}

template<typename T>
class gen_inc {

public:

    gen_inc() :
        inc_(0) {
    }

    T operator()() {
        return inc_++;
    }
private:
    T inc_;
};

#if 0
std::tr1::unordered_map<std::pair<int,int>,aligned_buffer<double>*> le_vs;
std::tr1::unordered_map<std::pair<int,int>,aligned_buffer<double>*> ri_vs;

/* The functions here are organized in a similar way as in evaluateGenericSpecial.c
 I provide generic, slow but readable function implementations for computing the
 conditional likelihood arrays at p, given child nodes q and r. Once again we need
 two generic function implementations, one for CAT and one for GAMMA */
template<typename state_t, size_t VW, size_t states>
void newviewCAT_FLEX(int tipCase, state_t *extEV, int *cptr, state_t *x1, state_t *x2, state_t *x3, state_t *tipVector, int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, state_t *left, state_t *right, int *wgt, int *scalerIncrement ) {

    typedef vector_unit<state_t,VW> vu;
    typedef typename vu::vec_t vec_t;

    double *le, *ri, ump_x1, ump_x2, x1px2;

    int i, l, j, scale, addScale = 0;

    const int statesSquare = states * states;

    /* here we switch over the different cases for efficiency, but also because
	 each case accesses different data types.

	 We consider three cases: either q and r are both tips, q or r are tips, and q and r are inner
	 nodes.
     */

    //    double test1[32] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
#if 0
    std::vector<double>test1(32);
    std::generate(test1.begin(), test1.end(), gen_inc<double>());
    double test2[32];

    reorder_back( test2, 8, 4, test1.begin(), 2 );

    std::copy( test2, test2 + 32, std::ostream_iterator<double>(std::cout, " " ) );
    std::cout << "\n";
    getchar();
#endif

    aligned_buffer<double>tip_l( VW * states );
    aligned_buffer<double>tip_r( VW * states );

    aligned_buffer<double>x3_ro( n * states, 0.0 );

    aligned_buffer<double>le_v( VW * states * states );
    aligned_buffer<double>ri_v( VW * states * states );

    std::vector<double> x3_tmp( states * n );

    //std::fill( x3_ro.begin(), x3_ro.end(), 0.0 );

    const size_t nb = n / VW;
    assert( n % VW == 0 );

    switch (tipCase) {

    /* both child nodes of p whether we want to update the conditional likelihood are tips */
    case TIP_TIP:

        /* loop over sites */

        for (size_t ib = 0; ib < nb; ++ib) {
            i = ib * VW;
            /* set a pointer to the P-Matrices for the rate category of this site */
            reorder_p( le_v.begin(), cptr + i, left, states, VW );
            reorder_p( ri_v.begin(), cptr + i, right, states, VW );

            //
            /* pointers to the likelihood entries of the tips q (vl) and r (vr)
			 We will do reading accesses to these values only.
             */
            //double *vl = &(tipVector[states * tipX1[i]]);
            // double *vr = &(tipVector[states * tipX2[i]]);
            //            std::cout << "vl:\n";
            //            std::copy( vl, vl + states, std::ostream_iterator<double>(std::cout, " ") );
            //            std::cout << "\n";
            //            vl = &(tipVector[states * tipX1[i + 1]]);
            //            std::copy( vl, vl + states, std::ostream_iterator<double>(std::cout, " ") );
            //            std::cout << "\n";
            reorder_tip_block(tipX1, tipVector, i, states, tip_l.begin(), VW);
            reorder_tip_block(tipX2, tipVector, i, states, tip_r.begin(), VW);

            //
            //
            //            getchar();
            /* address of the conditional likelihood array entres at site i. This is
			 a writing access to v */
            //double *v = &x3[states * i];
            /* initialize v */
            //            for (l = 0; l < states; l++)
            //                v[l] = 0.0;
            /* loop over states to compute the cond likelihoods at p (v) */

            for (l = 0; l < states; l++) {
                ump_x1 = 0.0;
                ump_x2 = 0.0;

                vec_t ump_x1_v = vu::setzero();
                vec_t ump_x2_v = vu::setzero();

                /* le and ri are the P-matrices */

                for (j = 0; j < states; j++) {
                    //                    ump_x1 += vl[j] * le[l * states + j];
                    //                    ump_x2 += vr[j] * ri[l * states + j];

                    ump_x1_v = vu::add( ump_x1_v,
                            vu::mul( vu::load(tip_l(j * VW)),
                                    vu::load(le_v((l * states + j) * VW))
                            )
                    );

                    ump_x2_v = vu::add( ump_x2_v,
                            vu::mul( vu::load(tip_r(j * VW)),
                                    vu::load(ri_v((l * states + j) * VW))
                            )
                    );

                }

                // x1px2 = ump_x1 * ump_x2;

                const vec_t x1px2_v = vu::mul( ump_x1_v, ump_x2_v );

                /* multiply with matrix of eigenvectors extEV */

                // TODO: try to convince the compiler to keep the v_accs in registers
                for (j = 0; j < states; j++) {

                    vec_t v_acc = vu::load( x3_ro((ib * states + j) * VW));

                    //v[j] += x1px2 * extEV[l * states + j];

                    v_acc = vu::add( v_acc, vu::mul( x1px2_v, vu::set1(extEV[l * states + j])));
                    vu::store( v_acc, x3_ro((ib * states + j) * VW));

                    //std::cout << "store: " << (ib * states + j ) * VW << "\n";
                }
            }
        }

        //        getchar();
#if 0
        reorder_back( x3_tmp.data(), n, states, x3_ro.begin(), VW );
        {
            bool eq = std::equal( x3_tmp.begin(), x3_tmp.end(), x3, delta_equal<double>( 1e-3 ) );

            std::cout << "equal: " << eq << "\n";

            if( !eq ) {
                print_delta( x3_tmp.begin(), x3_tmp.end(), x3 );
                std::cout << "<<<<<<<\n";
                //print_delta( x3_ro.begin(), x3_ro.end(), x3 );
                print_delta( x3_ro.begin(), x3_ro.end(), x3_tmp.begin() );
                getchar();
            }
        }
#endif
        break;
    case TIP_INNER:

        /* same as above, only that now vl is a tip and vr is the conditional probability vector
		 at an inner node. Note that, if we have the case that either q or r is a tip, the
		 nodes will be flipped to ensure that tipX1 always points to the sequence at the tip.
         */

        for (i = 0; i < n; i++) {
            le = &left[cptr[i] * statesSquare];
            ri = &right[cptr[i] * statesSquare];

            /* access tip vector lookup table */
            double *vl = &(tipVector[states * tipX1[i]]);

            /* access conditional likelihoo arrays */
            /* again, vl and vr are reading accesses, while v is a writing access */
            double *vr = &x2[states * i];
            double *v = &x3[states * i];

            /* same as in the loop above */

            for (l = 0; l < states; l++)
                v[l] = 0.0;

            for (l = 0; l < states; l++) {
                ump_x1 = 0.0;
                ump_x2 = 0.0;

                for (j = 0; j < states; j++) {
                    ump_x1 += vl[j] * le[l * states + j];
                    ump_x2 += vr[j] * ri[l * states + j];
                }

                x1px2 = ump_x1 * ump_x2;

                for (j = 0; j < states; j++)
                    v[j] += x1px2 * extEV[l * states + j];
            }

            /* now let's check for numerical scaling.
			 The maths in RAxML are a bit non-standard to avoid/economize on arithmetic operations
			 at the virtual root and for branch length optimization and hence values stored
			 in the conditional likelihood vectors can become negative.
			 Below we check if all absolute values stored at position i of v are smaller
			 than a pre-defined value in axml.h. If they are all smaller we can then safely
			 multiply them by a large, constant number twotothe256 (without numerical overflow)
			 that is also speced in axml.h */

            scale = 1;
            for (l = 0; scale && (l < states); l++)
                scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));

            if (scale) {
                for (l = 0; l < states; l++)
                    v[l] *= twotothe256;

                /* if we have scaled the entries to prevent underflow, we need to keep track of how many scaling
				 multiplications we did per node such as to undo them at the virtual root, e.g., in
				 evaluateGeneric()
				 Note here, that, if we scaled the site we need to increment the scaling counter by the wieght, i.e.,
				 the number of sites this potentially compressed pattern represents ! */

                addScale += wgt[i];
            }
        }
        break;
    case INNER_INNER:

        /* same as above, only that the two child nodes q and r are now inner nodes */

        for (i = 0; i < n; i++) {
            le = &left[cptr[i] * statesSquare];
            ri = &right[cptr[i] * statesSquare];

            /* index conditional likelihood vectors of inner nodes */

            double *vl = &x1[states * i];
            double *vr = &x2[states * i];
            double *v = &x3[states * i];

            for (l = 0; l < states; l++)
                v[l] = 0.0;

            for (l = 0; l < states; l++) {
                ump_x1 = 0.0;
                ump_x2 = 0.0;

                for (j = 0; j < states; j++) {
                    ump_x1 += vl[j] * le[l * states + j];
                    ump_x2 += vr[j] * ri[l * states + j];
                }

                x1px2 = ump_x1 * ump_x2;

                for (j = 0; j < states; j++)
                    v[j] += x1px2 * extEV[l * states + j];
            }

            scale = 1;
            for (l = 0; scale && (l < states); l++)
                scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));

            if (scale) {
                for (l = 0; l < states; l++)
                    v[l] *= twotothe256;

                addScale += wgt[i];
            }
        }
        break;
    default:
        assert(0);
    }

    /* increment the scaling counter by the additional scalings done at node p */

    *scalerIncrement = addScale;
}

#endif

struct arrays {

    bool valid() {
        return !u_x1.empty();
    }

    aligned_buffer<double> u_x1; //(VW * span);
    aligned_buffer<double> u_x2; //(VW * span);

//    aligned_buffer<double> x1_ro; //(span * n, 0.0);
//    aligned_buffer<double> x2_ro; //(span * n, 0.0);
//
//    aligned_buffer<double> x3_ro; //(span * n, 0.0);
    aligned_buffer<double> x3_tmp; //(span * n, 0.0);

    aligned_buffer<double> umpX1;
    aligned_buffer<double> umpX2;

};


template<const size_t VW>
inline int virtual_width( int n ) {
    return (n + 1) / VW * VW;
}

arrays arr;

template<typename state_t, size_t VW, size_t states>
void newviewGAMMA_FLEX(int tipCase, double *x1, double *x2, double *x3,
        double *extEV, double *tipVector, int *ex3, unsigned char *tipX1,
        unsigned char *tipX2, int n, double *left, double *right, int *wgt,
        int *scalerIncrement, const int maxStateValue) {
    typedef vector_unit<state_t, VW> vu;
    typedef typename vu::vec_t vec_t;

    double *v, x1px2, *vl, *vr, al, ar;

    int i, j, l, k, scale, addScale = 0;

    const int statesSquare = states * states, span = states * 4,
            /* this is required for doing some pre-computations that help to save
	 numerical operations. What we are actually computing here are additional lookup tables
	 for each possible state a certain data-type can assume.
	 for DNA with ambuguity coding this is 15, for proteins this is 22 or 23, since there
	 also exist one or two amibguity codes for protein data.
	 Essentially this is very similar to the tip vectors which we also use as lookup tables */
            precomputeLength = maxStateValue * span;

    bool did_scaling = false;
    const bool do_check = false;
    const bool reorder_input = false;


    ticks t1, t2;

    int vn = virtual_width<VW>( n );
    
//     std::cout << "n: " << n << " " << vn << "\n";
    
    //assert( n % VW == 0);

    if (!arr.valid()) {
        arr.u_x1.resize(VW * span);
        arr.u_x2.resize(VW * span);

//        arr.x1_ro.resize(span * n, 0.0);
//        arr.x2_ro.resize(span * n, 0.0);
//
//        arr.x3_ro.resize(span * n, 0.0);
        
        arr.umpX1.resize(precomputeLength);
        arr.umpX2.resize(precomputeLength);

    }
    if( arr.x3_tmp.size() < span * vn ) {
        arr.x3_tmp.resize( span * vn, 0.0 );
    }

    
    switch (tipCase) {
    case TIP_TIP: {
        /* allocate pre-compute memory space */

        //        double *umpX1 = (double*) malloc(sizeof(double) * precomputeLength),
        //                *umpX2 = (double*) malloc(sizeof(double) * precomputeLength);
        /* multiply all possible tip state vectors with the respective P-matrices
         */


//         t1 = getticks();
        for (i = 0; i < maxStateValue; i++) {
            v = &(tipVector[states * i]);

            for (k = 0; k < span; k++) {

//                 arr.umpX1[span * i + k] = 0.0;
//                 arr.umpX2[span * i + k] = 0.0;

// using temp values for accumulation helps in this case (but not for TIP/INNER)
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                
                for (l = 0; l != states; l++) {
                    tmp1 += v[l] * left[k * states + l];
                    tmp2 += v[l] * right[k * states + l];
                }
                arr.umpX1[span * i + k] = tmp1;
                arr.umpX2[span * i + k] = tmp2;

            }
        }

        
        
//         std::vector<uint8_t> last_tip1;
//         std::vector<uint8_t> last_tip2;
//         last_tip1.push_back(uint8_t(-1));
//         last_tip1.push_back(uint8_t(-1));
//         last_tip2.push_back(uint8_t(-1));
//         last_tip2.push_back(uint8_t(-1));
//         
        
        for (i = 0; i < vn; i += VW) {
            /* access the precomputed arrays (pre-computed multiplication of conditional with the tip state)
             */

//             if( !std::equal( last_tip1.begin(), last_tip1.end(), tipX1 + i ) ) {
                reorder_tip_block(tipX1, arr.umpX1.data(), i, span,
                                  arr.u_x1.begin(), VW, n);
//                 std::copy( tipX1 + i, tipX1 + i + VW, last_tip1.begin() );
//             }
            
//             if( !std::equal( last_tip2.begin(), last_tip2.end(), tipX2 + i ) ) {
            
                reorder_tip_block(tipX2, arr.umpX2.data(), i, span,
                                  arr.u_x2.begin(), VW, n);
//                 std::copy( tipX2 + i, tipX2 + i + VW, last_tip2.begin() );
//             }
            /* loop over discrete GAMMA rates */

            for (j = 0; j < 4; j++) {
                /* the rest is the same as for CAT */

                vec_t v_accs[states];

                //                for (k = 0; k < states; k++) {
                //                    //vu::store( vu::setzero(), arr.x3_ro(i * span + (j * states + k) * VW));
                //
                //                    v_accs[k] = vu::setzero();
                //                }
                vec_t v_acc0 = vu::setzero();
                vec_t v_acc1 = vu::setzero();
                vec_t v_acc2 = vu::setzero();
                vec_t v_acc3 = vu::setzero();

                for (k = 0; k < states; k++) {

                    vec_t x1px2_v = vu::mul(
                            vu::load(arr.u_x1((j * states + k) * VW)),
                            vu::load(arr.u_x2((j * states + k) * VW)));

                    //                    for (l = 0; l < states; l++) {
                    //
                    //                      //  double * const __restrict v_addr = arr.x3_ro(i * span + (j * states + l) * VW);
                    //
                    //                       // vec_t v_acc = vu::load( v_addr );
                    //
                    //                        v_accs[l] = vu::add( v_accs[l], vu::mul( x1px2_v, vu::set1(extEV[states * k + l])));
                    //
                    //                        //vu::store( v_acc, v_addr );
                    //                    }

                    // unrolled l-loop
                    v_acc0 = vu::add(v_acc0,
                            vu::mul(x1px2_v, vu::set1(extEV[states * k + 0])));
                    v_acc1 = vu::add(v_acc1,
                            vu::mul(x1px2_v, vu::set1(extEV[states * k + 1])));
                    v_acc2 = vu::add(v_acc2,
                            vu::mul(x1px2_v, vu::set1(extEV[states * k + 2])));
                    v_acc3 = vu::add(v_acc3,
                            vu::mul(x1px2_v, vu::set1(extEV[states * k + 3])));

#if 0
                    {
                        double * vx = x3_ro(i * span + (j * states) * VW);
                        std::copy( vx, vx + span, std::ostream_iterator<double>(std::cout, " " ));
                        std::cout << "<<<\n";

                    }

#endif
                }

                //                for (k = 0; k < states; k++) {
                //                    vu::store( v_accs[k], arr.x3_ro(i * span + (j * states + k) * VW));
                //
                //
                //                }
                vu::store(v_acc0, x3 + (i * span + (j * states + 0) * VW));
                vu::store(v_acc1, x3 + (i * span + (j * states + 1) * VW));
                vu::store(v_acc2, x3 + (i * span + (j * states + 2) * VW));
                vu::store(v_acc3, x3 + (i * span + (j * states + 3) * VW));

            }
        }

        /* free precomputed vectors */
        //
        //        free(umpX1);
        //        free(umpX2);


//         t2 = getticks();


        //        std::copy( x3, x3 + n * span, std::ostream_iterator<double>(std::cout, " "));
        //        std::cout << "<<<<<< I/I\n";

    }
    break;
    case TIP_INNER: {
        /* we do analogous pre-computations as above, with the only difference that we now do them
		 only for one tip vector */

        //        double *umpX1 = (double*) malloc(sizeof(double) * precomputeLength),
        //                *ump_x2 = (double*) malloc(sizeof(double) * states);
        /* precompute P and left tip vector product */

//        if( reorder_input ) {
//            reorder(x2, n, span, arr.x2_ro.begin(), VW);
//        } else {
//            std::copy( x2, x2 + n * span, arr.x2_ro.begin() );
//        }

//         t1 = getticks();

        for (i = 0; i < maxStateValue; i++) {
            v = &(tipVector[states * i]);

            for (k = 0; k < span; k++) {

                arr.umpX1[span * i + k] = 0.0;
//                 double tmp1 = 0.0;
                
                
                for (l = 0; l != states; l++) {
                    arr.umpX1[span * i + k] += v[l] * left[k * states + l];
//                     tmp1 += v[l] * left[k * states + l];
                    
                    
                }
//                 arr.umpX1[span * i + k] = tmp1;
            }
        }



        for (i = 0; i < vn; i+=VW) {
            /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */

            reorder_tip_block(tipX1, arr.umpX1.data(), i, span, arr.u_x1.begin(), VW, n );
            //uX1 = &arr.umpX1[span * tipX1[i]];

            /* loop over discrete GAMMA rates */

            vec_t max_v = vu::setzero();

            for (k = 0; k < 4; k++) {

                //				vr = &(x2[span * i + states * k]);
                //                v = &(x3[span * i + states * k]);
                //
                //                for (l = 0; l < states; l++)
                //                    v[l] = 0;

                vec_t v_acc0 = vu::setzero();
                vec_t v_acc1 = vu::setzero();
                vec_t v_acc2 = vu::setzero();
                vec_t v_acc3 = vu::setzero();

                for (l = 0; l < states; l++) {

                    //                    al = 0.0;
                    //                    ar = 0.0;

                    vec_t ar_v = vu::setzero();

                    for (j = 0; j < states; j++) {

                        ar_v = vu::add(
                                ar_v,
                                vu::mul(
                                        vu::load( x2 + (span * i + (states * k + j)* VW)),
                                        vu::set1( right[k * statesSquare + l * states + j])));

                        //                        al += vl[j] * left[k * statesSquare + l * states + j];
                        //                        ar += vr[j] * right[k * statesSquare + l * states + j];
                    }

                    //                    x1px2 = al * ar;

                    vec_t x1px2_v = vu::mul(
                            vu::load(arr.u_x1((k * states + l) * VW)), ar_v);

                    v_acc0 = vu::add(v_acc0,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 0])));
                    v_acc1 = vu::add(v_acc1,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 1])));
                    v_acc2 = vu::add(v_acc2,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 2])));
                    v_acc3 = vu::add(v_acc3,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 3])));

                    //
                    //                    if( l == 3 ) {
                    //
                    //
                    //                    }

                    //                    for (j = 0; j < states; j++)
                    //                        v[j] += x1px2 * extEV[states * l + j];

                }

                const vec_t max01_v = vu::max(vu::abs(v_acc0), vu::abs(v_acc1));
                const vec_t max23_v = vu::max(vu::abs(v_acc2), vu::abs(v_acc3));
                max_v = vu::max(max_v, vu::max(max01_v, max23_v));

                vu::store(v_acc0, x3 + (i * span + (k * states + 0) * VW));
                vu::store(v_acc1, x3 + (i * span + (k * states + 1) * VW));
                vu::store(v_acc2, x3 + (i * span + (k * states + 2) * VW));
                vu::store(v_acc3, x3 + (i * span + (k * states + 3) * VW));

            }
            /* also do numerical scaling as above. Note that here we need to scale
			 4 * 4 values for DNA or 4 * 20 values for protein data.
			 If they are ALL smaller than our threshold, we scale. Note that,
			 this can cause numerical problems with GAMMA, if the values generated
			 by the four discrete GAMMA rates are too different.

			 For details, see:

			 F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees"

             */

            double max_tmp[VW] __attribute__ ((aligned (32)));

            vu::store(max_v, max_tmp);

            //double max_val = std::max( max_tmp[0], max_tmp[1] );

            bool do_scale = false;

            for (j = 0; j < VW; ++j) {
                if (max_tmp[j] < minlikelihood) {
                    max_tmp[j] = twotothe256;
                    do_scale = true;
                    addScale += wgt[i + j];

                } else {
                    max_tmp[j] = 1.0;
                }
            }

            const vec_t scale_v = vu::load(max_tmp);
            //
            //            for( k = 0; k < span; ++k ) {
            //                max_val = std::max( max_val, fabs(arr.x3_ro[i*span + k]));
            //            }

            //std::cout << "max val: " << max_val << "\n";
            if (do_scale) {

                //                std::cout << "scale sim\n";

                for (k = 0; k < 4; ++k) {
                    double * __restrict base = x3 + (
                            i * span + (k * states + 0) * VW);

                    vec_t v_acc0 = vu::load(base);
                    vec_t v_acc1 = vu::load(base + VW);
                    vec_t v_acc2 = vu::load(base + 2 * VW);
                    vec_t v_acc3 = vu::load(base + 3 * VW);

                    v_acc0 = vu::mul(scale_v, v_acc0);
                    v_acc1 = vu::mul(scale_v, v_acc1);
                    v_acc2 = vu::mul(scale_v, v_acc2);
                    v_acc3 = vu::mul(scale_v, v_acc3);

                    vu::store(v_acc0, base);
                    vu::store(v_acc1, base + VW);
                    vu::store(v_acc2, base + 2 * VW);
                    vu::store(v_acc3, base + 3 * VW);

                }
            }
            did_scaling |= do_scale;

        }

        //        std::copy( x3, x3 + n * span, std::ostream_iterator<double>(std::cout, " "));
        //        std::cout << "<<<<<< T/I\n";
//         t2 = getticks();



    }
    break;
    case INNER_INNER: {
        /* same as above, without pre-computations */

//        if( reorder_input ) {
//            reorder(x1, n, span, arr.x1_ro.begin(), VW);
//            reorder(x2, n, span, arr.x2_ro.begin(), VW);
//        } else {
//            std::copy( x1, x1 + n * span, arr.x1_ro.begin() );
//            std::copy( x2, x2 + n * span, arr.x2_ro.begin() );
//        }
        t1 = getticks();

        for (i = 0; i < vn; i += VW) {
            vec_t max_v = vu::setzero();

            for (k = 0; k < 4; k++) {
                vl = &(x1[span * i + states * k]);
                vr = &(x2[span * i + states * k]);
                //                v = &(x3[span * i + states * k]);
                //
                //                for (l = 0; l < states; l++)
                //                    v[l] = 0;

                vec_t v_acc0 = vu::setzero();
                vec_t v_acc1 = vu::setzero();
                vec_t v_acc2 = vu::setzero();
                vec_t v_acc3 = vu::setzero();

                for (l = 0; l < states; l++) {

                    //                    al = 0.0;
                    //                    ar = 0.0;
                    vec_t al_v = vu::setzero();

                    vec_t ar_v = vu::setzero();

                    for (j = 0; j < states; j++) {
                        al_v = vu::add(
                                al_v,
                                vu::mul(
                                        vu::load(
                                                x1 + (
                                                        span * i
                                                        + (states * k
                                                                + j)
                                                                * VW)),
                                                                vu::set1(
                                                                        left[k * statesSquare
                                                                             + l * states + j])));

                        ar_v = vu::add(
                                ar_v,
                                vu::mul(
                                        vu::load(
                                                x2 + (
                                                        span * i
                                                        + (states * k
                                                                + j)
                                                                * VW)),
                                                                vu::set1(
                                                                        right[k * statesSquare
                                                                              + l * states + j])));

                        //                        al += vl[j] * left[k * statesSquare + l * states + j];
                        //                        ar += vr[j] * right[k * statesSquare + l * states + j];
                    }

                    //                    x1px2 = al * ar;

                    vec_t x1px2_v = vu::mul(al_v, ar_v);

                    v_acc0 = vu::add(v_acc0,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 0])));
                    v_acc1 = vu::add(v_acc1,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 1])));
                    v_acc2 = vu::add(v_acc2,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 2])));
                    v_acc3 = vu::add(v_acc3,
                            vu::mul(x1px2_v, vu::set1(extEV[states * l + 3])));

                    //
                    //                    if( l == 3 ) {
                    //
                    //
                    //                    }

                    //                    for (j = 0; j < states; j++)
                    //                        v[j] += x1px2 * extEV[states * l + j];

                }

                const vec_t max01_v = vu::max(vu::abs(v_acc0), vu::abs(v_acc1));
                const vec_t max23_v = vu::max(vu::abs(v_acc2), vu::abs(v_acc3));
                max_v = vu::max(max_v, vu::max(max01_v, max23_v));

                vu::store(v_acc0, x3 + (i * span + (k * states + 0) * VW));
                vu::store(v_acc1, x3 + (i * span + (k * states + 1) * VW));
                vu::store(v_acc2, x3 + (i * span + (k * states + 2) * VW));
                vu::store(v_acc3, x3 + (i * span + (k * states + 3) * VW));

            }

            //            _mm_prefetch( arr.x3_ro(i * span + (0 * states + 0) * VW), _MM_HINT_T0 );

            //            for( k = 0; k < 4; ++k ) {
            //                vec_t v_acc0 = vu::load( arr.x3_ro(i * span + (k * states + 0) * VW));
            //                vec_t v_acc1 = vu::load( arr.x3_ro(i * span + (k * states + 1) * VW));
            //                vec_t v_acc2 = vu::load( arr.x3_ro(i * span + (k * states + 2) * VW));
            //                vec_t v_acc3 = vu::load( arr.x3_ro(i * span + (k * states + 3) * VW));
            //
            //                const vec_t max01_v = vu::max( vu::abs(v_acc0), vu::abs(v_acc1) );
            //                const vec_t max23_v = vu::max( vu::abs(v_acc2), vu::abs(v_acc3) );
            //                max_v = vu::max( max_v, vu::max( max01_v, max23_v ));
            //            }
#if 1
            //            const vec_t cmp_v = vu::cmp_lt( max_v, vu::set1(minlikelihood));
            //            const vec_t scaletrue_v = vu::bit_and( cmp_v, vu::set1(twotothe256));
            //            const vec_t scalefalse_v = vu::bit_andnot( cmp_v, vu::set1(1));
            //            const vec_t scale_v = vu::bit_or( scaletrue_v, scalefalse_v );

            //const vec_t scale_v = vu::set1(twotothe256);

            double max_tmp[VW] __attribute__ ((aligned (32)));

            vu::store(max_v, max_tmp);

            //double max_val = std::max( max_tmp[0], max_tmp[1] );

            bool do_scale = false;

            for (j = 0; j < VW; ++j) {
                if (max_tmp[j] < minlikelihood) {
                    max_tmp[j] = twotothe256;
                    do_scale = true;
                    addScale += wgt[i + j];

                } else {
                    max_tmp[j] = 1.0;
                }
            }

            const vec_t scale_v = vu::load(max_tmp);
            //
            //            for( k = 0; k < span; ++k ) {
            //                max_val = std::max( max_val, fabs(arr.x3_ro[i*span + k]));
            //            }

            //std::cout << "max val: " << max_val << "\n";
            if (do_scale) {

                //                std::cout << "scale sim\n";

                for (k = 0; k < 4; ++k) {
                    double * __restrict base = x3 + (i * span + (k * states + 0) * VW);
                    vec_t v_acc0 = vu::load(base);
                    vec_t v_acc1 = vu::load(base + VW);
                    vec_t v_acc2 = vu::load(base + 2 * VW);
                    vec_t v_acc3 = vu::load(base + 3 * VW);

                    v_acc0 = vu::mul(scale_v, v_acc0);
                    v_acc1 = vu::mul(scale_v, v_acc1);
                    v_acc2 = vu::mul(scale_v, v_acc2);
                    v_acc3 = vu::mul(scale_v, v_acc3);

                    vu::store(v_acc0, base);
                    vu::store(v_acc1, base + VW);
                    vu::store(v_acc2, base + 2 * VW);
                    vu::store(v_acc3, base + 3 * VW);

                }
            }
            did_scaling |= do_scale;
#endif

        }

        t2 = getticks();



    }
    break;
    default:
        assert(0);
    }

//    if (did_scaling) {
//        std::cout << "sim scaling: " << addScale << "\n";
//    }

   // std::cout << tipCase << " ticks: " << elapsed(t2, t1) << "\n";
//    if( do_check )
//    {
//        reorder_back(arr.x3_tmp.data(), n, span, arr.x3_ro.begin(), VW);
//        bool eq = std::equal(arr.x3_tmp.begin(), arr.x3_tmp.end(), x3,
//                delta_equal<double>(1e-3));
//
//        std::cout << "equal: " << eq << "\n";
//
//        if (!eq) {
//            std::cout << "meeeeeeeep. not equal\n";
//            print_delta(arr.x3_tmp.begin(), arr.x3_tmp.end(), x3);
//            std::cout << "<<<<<<<\n";
//            //print_delta( x3_ro.begin(), x3_ro.end(), x3 );
//            print_delta(arr.x3_ro.begin(), arr.x3_ro.end(),
//                    arr.x3_tmp.begin());
//            throw std::runtime_error("baililng out");
//            getchar();
//
//        }
//    }

    //print_delta(arr.x3_ro.begin(), arr.x3_ro.end(), x3);

    //getchar();
  //  std::copy( arr.x3_ro.begin(), arr.x3_ro.end(), x3 );
    /* as above, increment the global counter that counts scaling multiplications by the scaling multiplications
	 carried out for computing the likelihood array at node p */

    *scalerIncrement = addScale;
}

/* same thing for GAMMA models. The only noteworthy thing here is that we have an additional inner loop over the
   number of discrete gamma rates. The data access pattern is also different since for tip vector accesses through our
   lookup table, we do not distnguish between rates

   Note the different access pattern in TIP_INNER:

   left = &(tipVector[states * tipX1[i]]);
   right = &(x2[span * i + l * states]);

*/

template<typename state_t, size_t VW, size_t states>
void sumGAMMA_FLEX(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
        unsigned char *tipX1, unsigned char *tipX2, int n )
{
    int
    i,
    l,
    k;

    const int
    span = 4 * states;

    double
    *left,
    *right,
    *sum;


    typedef vector_unit<state_t,VW> vu;
    typedef typename vu::vec_t vec_t;


//     assert( n % VW == 0 );
    int vn = virtual_width<VW>(n);
        
    switch(tipCase)
    {
    case TIP_TIP:
        for(i = 0; i < vn; i++)
        {
            left  = &(tipVector[states * tipX1[i]]);
            right = &(tipVector[states * tipX2[i]]);

            for(l = 0; l < 4; l++)
            {
                sum = &sumtable[i * span + l * states];

                for(k = 0; k < states; k++)
                    sum[k] = left[k] * right[k];

            }
        }
        break;
    case TIP_INNER:
        for(i = 0; i < vn; i+=VW)
        {

            reorder_tip_block(tipX1, tipVector, i, states, arr.u_x1.begin(), VW, n );

            //left = &(tipVector[states * tipX1[i]]);

            for(l = 0; l < 4; l++)
            {
                right = &(x2[span * i + (l * states * VW)]);

                sum = &sumtable[i * span + l * states];

                for(k = 0; k < states; k++) {
                   // sum[k] = left[k] * right[k];

                    double tmp[VW] __attribute__ ((aligned (32)));

                    vu::store( vu::mul( vu::load(arr.u_x1(k * VW)), vu::load( right + k * VW )), tmp );
                    sum[k] = tmp[0];
                    sum[k+span] = tmp[1];
                }
            }
        }
        break;
    case INNER_INNER:
        for(i = 0; i < vn; i+=VW)
        {
            for(l = 0; l < 4; l++)
            {
                left  = &(x1[span * i + (l * states * VW)]);
                right = &(x2[span * i + (l * states * VW)]);
                sum   = &(sumtable[i * span + (l * states)]);


                for(k = 0; k < states; k++) {
                    double tmp[VW] __attribute__ ((aligned (32)));

                    vu::store( vu::mul( vu::load(left + k * VW), vu::load( right + k * VW )), tmp );

                    sum[k] = tmp[0];
                    sum[k+span] = tmp[1];
                }
            }
        }
//        for(i = 0; i < n; i++)
//        {
//            for(l = 0; l < 4; l++)
//            {
//                left  = &(x1[span * i + l * states]);
//                right = &(x2[span * i + l * states]);
//                sum   = &(sumtable[i * span + l * states]);
//
//
//                for(k = 0; k < states; k++)
//                    sum[k] = left[k] * right[k];
//            }
//        }
        break;
    default:
        assert(0);
    }
}

}


const static size_t VW = 2;

extern "C" {
void newviewCAT_FLEX_reorder(int tipCase, double *extEV, int *cptr, double *x1,
        double *x2, double *x3, double *tipVector, int *ex3,
        unsigned char *tipX1, unsigned char *tipX2, int n, double *left,
        double *right, int *wgt, int *scalerIncrement, const int states);
void newviewGAMMA_FLEX_reorder(int tipCase, double *x1, double *x2, double *x3,
        double *extEV, double *tipVector, int *ex3, unsigned char *tipX1,
        unsigned char *tipX2, int n, double *left, double *right, int *wgt,
        int *scalerIncrement, const int states, const int maxStateValue);
void sumGAMMA_FLEX_reorder(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
                          unsigned char *tipX1, unsigned char *tipX2, int n, const int states);
void reorder_back( double *x, int n, int span ) {
    sim_reorder::reorder_back( sim_reorder::arr.x3_tmp.begin(), size_t(n), size_t(span), x, VW );
    std::copy( sim_reorder::arr.x3_tmp.begin(), sim_reorder::arr.x3_tmp.begin() + n * span, x );
}


void reorder( double *x, int n, int span ) {
    sim_reorder::reorder( x, size_t(n), size_t(span), sim_reorder::arr.x3_tmp.begin(), VW );
    std::copy( sim_reorder::arr.x3_tmp.begin(), sim_reorder::arr.x3_tmp.begin() + n * span, x );
}

}

void newviewCAT_FLEX_reorder(int tipCase, double *extEV, int *cptr, double *x1,
        double *x2, double *x3, double *tipVector, int *ex3,
        unsigned char *tipX1, unsigned char *tipX2, int n, double *left,
        double *right, int *wgt, int *scalerIncrement, const int states) {

    assert( states == 4);

    //    sim_reorder::newviewCAT_FLEX<double,2,4>(tipCase, extEV,
    //                cptr,
    //                x1, x2, x3, tipVector,
    //                ex3, tipX1, tipX2,
    //                 n, left, right, wgt, scalerIncrement );
    assert(0);
}

void newviewGAMMA_FLEX_reorder(int tipCase, double *x1, double *x2, double *x3,
        double *extEV, double *tipVector, int *ex3, unsigned char *tipX1,
        unsigned char *tipX2, int n, double *left, double *right, int *wgt,
        int *scalerIncrement, const int states, const int maxStateValue) {
    assert( states == 4);
    sim_reorder::newviewGAMMA_FLEX<double,VW,4>( tipCase, x1, x2, x3, extEV, tipVector, ex3, tipX1, tipX2, n, left, right, wgt, scalerIncrement, maxStateValue);
}


void sumGAMMA_FLEX_reorder(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
                          unsigned char *tipX1, unsigned char *tipX2, int n, const int states) {
    assert( states == 4 );
    sim_reorder::sumGAMMA_FLEX<double,VW,4>( tipCase, sumtable, x1, x2, tipVector, tipX1, tipX2, n );

}

