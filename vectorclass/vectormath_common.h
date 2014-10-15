/***************************  vectormath_common.h   ****************************
* Author:        Agner Fog
* Date created:  2014-04-18
* Last modified: 2014-07-23
* Version:       1.14
* Project:       vector classes
* Description:
* Header file containing common code for inline version of mathematical functions.
*
* Theory, methods and inspiration based partially on these sources:
* > Moshier, Stephen Lloyd Baluk: Methods and programs for mathematical functions.
*   Ellis Horwood, 1989.
* > VDT library developed on CERN by Danilo Piparo, Thomas Hauth and
*   Vincenzo Innocente, 2012, https://svnweb.cern.ch/trac/vdt
* > Cephes math library by Stephen L. Moshier 1992,
*   http://www.netlib.org/cephes/
*
* Calculation methods:
* Some functions are using Padé approximations f(x) = P(x)/Q(x)
* Most single precision functions are using Taylor expansions
*
* For detailed instructions, see VectorClass.pdf
*
* (c) Copyright 2014 GNU General Public License http://www.gnu.org/licenses
******************************************************************************/

#ifndef VECTORMATH_COMMON_H
#define VECTORMATH_COMMON_H  1

#ifdef VECTORMATH_LIB_H
#error conflicting header files: vectormath_lib.h for external math functions, other vectormath_xxx.h for inline math functions
#endif

#include <math.h>
#include "vectorclass.h"


/******************************************************************************
               fused multiply-and-add functions
******************************************************************************/

static inline Vec4f mul_add(Vec4f const & a, Vec4f const & b, Vec4f const & c) {
#ifdef __FMA__
    return _mm_fmadd_ps(a, b, c);
#elif defined (__FMA4__)
    return _mm_macc_ps(a, b, c);
#else
    return a * b + c;
#endif
}

static inline Vec4f mul_sub(Vec4f const & a, Vec4f const & b, Vec4f const & c) {
#ifdef __FMA__
    return _mm_fmsub_ps(a, b, c);
#elif defined (__FMA4__)
    return _mm_msub_ps(a, b, c);
#else
    return a * b - c;
#endif
}

// mul_sub_x has extra precision on the intermediate calculations, even if FMA not supported,
// using Veltkamp-Dekker split
static inline Vec4f mul_sub_x(Vec4f const & a, Vec4f const & b, Vec4f const & c) {
#ifdef __FMA__
    return _mm_fmsub_ps(a, b, c);
#elif defined (__FMA4__)
    return _mm_msub_ps(a, b, c);
#else
    // calculate a * b - c with extra precision
    Vec4i upper_mask = -(1 << 12);                         // mask to remove lower 12 bits
    Vec4f a_high = a & Vec4f(reinterpret_f(upper_mask));   // split into high and low parts
    Vec4f b_high = b & Vec4f(reinterpret_f(upper_mask));
    Vec4f a_low  = a - a_high;
    Vec4f b_low  = b - b_high;
    Vec4f r1 = a_high * b_high;                            // this product is exact
    Vec4f r2 = r1 - c;                                     // subtract c from high product
    Vec4f r3 = r2 + (a_high * b_low + b_high * a_low) + a_low * b_low; // add rest of product
    return r3; // + ((r2 - r1) + c);
#endif
}

static inline Vec2d mul_add(Vec2d const & a, Vec2d const & b, Vec2d const & c) {
#ifdef __FMA__
    return _mm_fmadd_pd(a, b, c);
#elif defined (__FMA4__)
    return _mm_macc_pd(a, b, c);
#else
    return a * b + c;
#endif
}

static inline Vec2d mul_sub(Vec2d const & a, Vec2d const & b, Vec2d const & c) {
#ifdef __FMA__
    return _mm_fmsub_pd(a, b, c);
#elif defined (__FMA4__)
    return _mm_msub_pd(a, b, c);
#else
    return a * b - c;
#endif
}

// mul_sub_x has extra precision on the intermediate calculations, even if FMA not supported,
// using Veltkamp-Dekker split
static inline Vec2d mul_sub_x(Vec2d const & a, Vec2d const & b, Vec2d const & c) {
#ifdef __FMA__
    return _mm_fmsub_pd(a, b, c);
#elif defined (__FMA4__)
    return _mm_msub_pd(a, b, c);
#else
    // calculate a * b - c with extra precision
    Vec2q upper_mask = -(1LL << 27);                       // mask to remove lower 27 bits
    Vec2d a_high = a & Vec2d(reinterpret_d(upper_mask));   // split into high and low parts
    Vec2d b_high = b & Vec2d(reinterpret_d(upper_mask));
    Vec2d a_low  = a - a_high;
    Vec2d b_low  = b - b_high;
    Vec2d r1 = a_high * b_high;                            // this product is exact
    Vec2d r2 = r1 - c;                                     // subtract c from high product
    Vec2d r3 = r2 + (a_high * b_low + b_high * a_low) + a_low * b_low; // add rest of product
    return r3; // + ((r2 - r1) + c);
#endif
}

#if MAX_VECTOR_SIZE >= 256

static inline Vec8f mul_add(Vec8f const & a, Vec8f const & b, Vec8f const & c) {
#ifdef __FMA__
    return _mm256_fmadd_ps(a, b, c);
#elif defined (__FMA4__)
    return _mm256_macc_ps(a, b, c);
#else
    return a * b + c;
#endif
    
}

static inline Vec8f mul_sub(Vec8f const & a, Vec8f const & b, Vec8f const & c) {
#ifdef __FMA__
    return _mm256_fmsub_ps(a, b, c);
#elif defined (__FMA4__)
    return _mm256_msub_ps(a, b, c);
#else
    return a * b - c;
#endif    
}

static inline Vec8f mul_sub_x(Vec8f const & a, Vec8f const & b, Vec8f const & c) {
#ifdef __FMA__
    return _mm256_fmsub_ps(a, b, c);
#elif defined (__FMA4__)
    return _mm256_msub_ps(a, b, c);
#else
    // calculate a * b - c with extra precision
    Vec8i upper_mask = -(1 << 12);                       // mask to remove lower 12 bits
    Vec8f a_high = a & Vec8f(reinterpret_f(upper_mask));   // split into high and low parts
    Vec8f b_high = b & Vec8f(reinterpret_f(upper_mask));
    Vec8f a_low  = a - a_high;
    Vec8f b_low  = b - b_high;
    Vec8f r1 = a_high * b_high;                            // this product is exact
    Vec8f r2 = r1 - c;                                     // subtract c from high product
    Vec8f r3 = r2 + (a_high * b_low + b_high * a_low) + a_low * b_low; // add rest of product
    return r3; // + ((r2 - r1) + c);
#endif
}

static inline Vec4d mul_add(Vec4d const & a, Vec4d const & b, Vec4d const & c) {
#ifdef __FMA__
    return _mm256_fmadd_pd(a, b, c);
#elif defined (__FMA4__)
    return _mm256_macc_pd(a, b, c);
#else
    return a * b + c;
#endif
    
}

static inline Vec4d mul_sub(Vec4d const & a, Vec4d const & b, Vec4d const & c) {
#ifdef __FMA__
    return _mm256_fmsub_pd(a, b, c);
#elif defined (__FMA4__)
    return _mm256_msub_pd(a, b, c);
#else
    return a * b - c;
#endif
   
}

static inline Vec4d mul_sub_x(Vec4d const & a, Vec4d const & b, Vec4d const & c) {
#ifdef __FMA__
    return _mm256_fmsub_pd(a, b, c);
#elif defined (__FMA4__)
    return _mm256_msub_pd(a, b, c);
#else
    // calculate a * b - c with extra precision
    Vec4q upper_mask = -(1LL << 27);                       // mask to remove lower 27 bits
    Vec4d a_high = a & Vec4d(reinterpret_d(upper_mask));   // split into high and low parts
    Vec4d b_high = b & Vec4d(reinterpret_d(upper_mask));
    Vec4d a_low  = a - a_high;
    Vec4d b_low  = b - b_high;
    Vec4d r1 = a_high * b_high;                            // this product is exact
    Vec4d r2 = r1 - c;                                     // subtract c from high product
    Vec4d r3 = r2 + (a_high * b_low + b_high * a_low) + a_low * b_low; // add rest of product
    return r3; // + ((r2 - r1) + c);
#endif
}

#endif  // MAX_VECTOR_SIZE >= 256

#if MAX_VECTOR_SIZE >= 512

static inline Vec16f mul_add(Vec16f const & a, Vec16f const & b, Vec16f const & c) {
#if INSTRSET >= 9
    return _mm512_fmadd_ps(a, b, c);
#else
    return Vec16f(mul_add(a.get_low(), b.get_low(), c.get_low()), mul_add(a.get_high(), b.get_high(), c.get_high()));
#endif
}

static inline Vec16f mul_sub(Vec16f const & a, Vec16f const & b, Vec16f const & c) {
#if INSTRSET >= 9
    return _mm512_fmsub_ps(a, b, c);
#else
    return Vec16f(mul_sub(a.get_low(), b.get_low(), c.get_low()), mul_sub(a.get_high(), b.get_high(), c.get_high()));
#endif
}

static inline Vec16f mul_sub_x(Vec16f const & a, Vec16f const & b, Vec16f const & c) {
#if INSTRSET >= 9
    return _mm512_fmsub_ps(a, b, c);
#else
    return Vec16f(mul_sub_x(a.get_low(), b.get_low(), c.get_low()), mul_sub_x(a.get_high(), b.get_high(), c.get_high()));
#endif
}

static inline Vec8d mul_add(Vec8d const & a, Vec8d const & b, Vec8d const & c) {
#if INSTRSET >= 9
    return _mm512_fmadd_pd(a, b, c);
#else
    return Vec8d(mul_add(a.get_low(), b.get_low(), c.get_low()), mul_add(a.get_high(), b.get_high(), c.get_high()));
#endif
}

static inline Vec8d mul_sub(Vec8d const & a, Vec8d const & b, Vec8d const & c) {
#if INSTRSET >= 9
    return _mm512_fmsub_pd(a, b, c);
#else
    return Vec8d(mul_sub(a.get_low(), b.get_low(), c.get_low()), mul_sub(a.get_high(), b.get_high(), c.get_high()));
#endif
}

static inline Vec8d mul_sub_x(Vec8d const & a, Vec8d const & b, Vec8d const & c) {
#if INSTRSET >= 9
    return _mm512_fmsub_pd(a, b, c);
#else
    return Vec8d(mul_sub_x(a.get_low(), b.get_low(), c.get_low()), mul_sub_x(a.get_high(), b.get_high(), c.get_high()));
#endif
}

#endif  // MAX_VECTOR_SIZE >= 512


/******************************************************************************
               define mathematical constants
******************************************************************************/
#define VM_PI       3.14159265358979323846           // pi
#define VM_PI_2     1.57079632679489661923           // pi / 2
#define VM_PI_4     0.785398163397448309616          // pi / 4
#define VM_SQRT2    1.41421356237309504880           // sqrt(2)
#define VM_LOG2E    1.44269504088896340736           // 1/log(2)
#define VM_LOG10E   0.434294481903251827651          // 1/log(10)
#define VM_LN2      0.693147180559945309417          // log(2)
#define VM_LN10     2.30258509299404568402           // log(10)
#define VM_SMALLEST_NORMAL  2.2250738585072014E-308  // smallest normal number, double
#define VM_SMALLEST_NORMALF 1.17549435E-38f          // smallest normal number, float


/******************************************************************************
      templates for producing infinite and nan in desired vector type
******************************************************************************/
template <class VTYPE>
static inline VTYPE infinite_vec();

template <>
inline Vec2d infinite_vec<Vec2d>() {
    return infinite2d();
}

template <>
inline Vec4f infinite_vec<Vec4f>() {
    return infinite4f();
}

#if MAX_VECTOR_SIZE >= 256

template <>
inline Vec4d infinite_vec<Vec4d>() {
    return infinite4d();
}

template <>
inline Vec8f infinite_vec<Vec8f>() {
    return infinite8f();
}

#endif // MAX_VECTOR_SIZE >= 256

#if MAX_VECTOR_SIZE >= 512

template <>
inline Vec8d infinite_vec<Vec8d>() {
    return infinite8d();
}

template <>
inline Vec16f infinite_vec<Vec16f>() {
    return infinite16f();
}

#endif // MAX_VECTOR_SIZE >= 512


// template for producing quiet NAN
template <class VTYPE>
static inline VTYPE nan_vec(int n = 0x100);

template <>
inline Vec2d nan_vec<Vec2d>(int n) {
    return nan2d(n);
}

template <>
inline Vec4f nan_vec<Vec4f>(int n) {
    return nan4f(n);
}

#if MAX_VECTOR_SIZE >= 256

template <>
inline Vec4d nan_vec<Vec4d>(int n) {
    return nan4d(n);
}

template <>
inline Vec8f nan_vec<Vec8f>(int n) {
    return nan8f(n);
}

#endif // MAX_VECTOR_SIZE >= 256

#if MAX_VECTOR_SIZE >= 512

template <>
inline Vec8d nan_vec<Vec8d>(int n) {
    return nan8d(n);
}

template <>
inline Vec16f nan_vec<Vec16f>(int n) {
    return nan16f(n);
}

#endif // MAX_VECTOR_SIZE >= 512

// Define NAN trace values
#define NAN_LOG 0x101  // logarithm for x<0
#define NAN_POW 0x102  // negative number raised to non-integer power
#define NAN_HYP 0x104  // acosh for x<1 and atanh for abs(x)>1


/******************************************************************************
                  templates for polynomials
Using Estrin's scheme to make shorter dependency chains and use FMA, starting
longest dependency chains first.
******************************************************************************/

// template <typedef VECTYPE, typedef CTYPE> 
template <class VTYPE, class CTYPE> 
static inline VTYPE polynomial_2(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2) {
    // calculates polynomial c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE y = x2 * c2 + (x * c1 + c0);
    return y;
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_3(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3) {
    // calculates polynomial c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    return (c2 + c3*x)*x2 + (c1*x + c0);
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_4(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4) {
    // calculates polynomial c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE x4 = x2 * x2;
    return (c2+c3*x)*x2 + ((c0+c1*x) + c4*x4);
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_4n(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3) {
    // calculates polynomial 1*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE x4 = x2 * x2;
    return (c2+c3*x)*x2 + ((c0+c1*x) + x4);
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_5(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5) {
    // calculates polynomial c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE x4 = x2 * x2;
    return (c2+c3*x)*x2 + ((c4+c5*x)*x4 + (c0+c1*x));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_5n(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4) {
    // calculates polynomial 1*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE x4 = x2 * x2;
    return (c2+c3*x)*x2 + ((c4+x)*x4 + (c0+c1*x));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_6(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5, CTYPE c6) {
    // calculates polynomial c6*x^6 + c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE x4 = x2 * x2;
    return  (c4+c5*x+c6*x2)*x4 + ((c2+c3*x)*x2 + (c0+c1*x));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_6n(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5) {
    // calculates polynomial 1*x^6 + c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE x4 = x2 * x2;
    return  (c4+c5*x+x2)*x4 + ((c2+c3*x)*x2 + (c0+c1*x));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_7(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5, CTYPE c6, CTYPE c7) {
    // calculates polynomial c7*x^7 + c6*x^6 + c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x * x;
    VTYPE x4 = x2 * x2;
    return  ((c6+c7*x)*x2 + (c4+c5*x))*x4 + ((c2+c3*x)*x2 + (c0+c1*x));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_8(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5, CTYPE c6, CTYPE c7, CTYPE c8) {
    // calculates polynomial c8*x^8 + c7*x^7 + c6*x^6 + c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x  * x;
    VTYPE x4 = x2 * x2;
    VTYPE x8 = x4 * x4;
    return  ((c6+c7*x)*x2 + (c4+c5*x))*x4 + (c8*x8 + (c2+c3*x)*x2 + (c0+c1*x));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_9(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5, CTYPE c6, CTYPE c7, CTYPE c8, CTYPE c9) {
    // calculates polynomial c9*x^9 + c8*x^8 + c7*x^7 + c6*x^6 + c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x  * x;
    VTYPE x4 = x2 * x2;
    VTYPE x8 = x4 * x4;
    return  (((c6+c7*x)*x2 + (c4+c5*x))*x4 + (c8+c9*x)*x8) + ((c2+c3*x)*x2 + (c0+c1*x));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_10(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5, CTYPE c6, CTYPE c7, CTYPE c8, CTYPE c9, CTYPE c10) {
    // calculates polynomial c10*x^10 + c9*x^9 + c8*x^8 + c7*x^7 + c6*x^6 + c5*x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x  * x;
    VTYPE x4 = x2 * x2;
    VTYPE x8 = x4 * x4;
    return  (((c6+c7*x)*x2 + (c4+c5*x))*x4 + (c8+c9*x+c10*x2)*x8) + ((c2+c3*x)*x2 + (c0+c1*x));
} 

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_13(VTYPE const & x, CTYPE c0, CTYPE c1, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5, CTYPE c6, CTYPE c7, CTYPE c8, CTYPE c9, CTYPE c10, CTYPE c11, CTYPE c12, CTYPE c13) {
    // calculates polynomial c13*x^13 + c12*x^12 + ... + c1*x + c0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x  * x;
    VTYPE x4 = x2 * x2;
    VTYPE x8 = x4 * x4;
    return  ((c8+c9*x) + (c10+c11*x)*x2 + (c12+c13*x)*x4)*x8 + 
        (((c6+c7*x)*x2 + (c4+c5*x))*x4 + ((c2+c3*x)*x2 + (c0+c1*x)));
}

template<class VTYPE, class CTYPE> 
static inline VTYPE polynomial_13m(VTYPE const & x, CTYPE c2, CTYPE c3, CTYPE c4, CTYPE c5, CTYPE c6, CTYPE c7, CTYPE c8, CTYPE c9, CTYPE c10, CTYPE c11, CTYPE c12, CTYPE c13) {
    // calculates polynomial c13*x^13 + c12*x^12 + ... + x + 0
    // VTYPE may be a vector type, CTYPE is a scalar type
    VTYPE x2 = x  * x;
    VTYPE x4 = x2 * x2;
    VTYPE x8 = x4 * x4;
    return  ((c8+c9*x) + (c10+c11*x)*x2 + (c12+c13*x)*x4)*x8 + 
        (((c6+c7*x)*x2 + (c4+c5*x))*x4 + ((c2+c3*x)*x2 + x));
}

#endif
