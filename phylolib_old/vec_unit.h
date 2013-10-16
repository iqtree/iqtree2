/*
 * Copyright (C) 2010 Simon A. Berger
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 */


#ifndef __vec_unit_h
#define __vec_unit_h
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <x86intrin.h>


#ifdef __AVX__
#define HAVE_AVX
//#include <immintrin.h>
#endif

template<class T, size_t W> 
struct vector_unit {
   
};

// vector unit specialization: SSE 8x16bit integer 

template<>
struct vector_unit<short, 8> {

    const static bool do_checks = false;
    
    typedef __m128i vec_t;
    typedef short T;

    
    const static T POS_MAX_VALUE = 0x7fff;
    const static T LARGE_VALUE = 32000;
    const static T SMALL_VALUE = -32000;
    const static T BIAS = 0;
    const static size_t W = 8;
    
    static inline vec_t setzero() {
        return set1(0);
    }
    
    static inline vec_t set1( T val ) {
        return _mm_set1_epi16( val );
    }
    
    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";
        
        _mm_store_si128( (vec_t*)addr, v );
    }
    
    static inline const vec_t load( const T* addr ) {
        return _mm_load_si128( (vec_t*)addr );
    }
    
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        return _mm_and_si128( a, b );
    }
    
    static inline const vec_t bit_or( const vec_t &a, const vec_t &b ) {
        return _mm_or_si128( a, b );
    }
    static inline const vec_t bit_andnot( const vec_t &a, const vec_t &b ) {
        return _mm_andnot_si128( a, b );
    }
    static inline const vec_t bit_invert( const vec_t &a ) {
        return _mm_xor_si128( a, set1(0xffff) );
    }
    
    
//     static inline const vec_t bit_invert( const vec_t &a ) {
//         //return _mm_andnot_pd(a, set1(0xffff));
//     }
    
    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        return _mm_add_epi16( a, b );
    }
    static inline const vec_t adds( const vec_t &a, const vec_t &b ) {
        return _mm_adds_epi16( a, b );
    }
    
    static inline const vec_t sub( const vec_t &a, const vec_t &b ) {
        return _mm_sub_epi16( a, b );
    }
    static inline const vec_t cmp_zero( const vec_t &a ) {
        return _mm_cmpeq_epi16( a, setzero() );
    }
    
    static inline const vec_t cmp_eq( const vec_t &a, const vec_t &b ) {
        return _mm_cmpeq_epi16( a, b );
    }
    
    static inline const vec_t cmp_lt( const vec_t &a, const vec_t &b ) {
     
        return _mm_cmplt_epi16( a, b );
    }
    
    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        return _mm_min_epi16( a, b );
    }
    
    static inline const vec_t max( const vec_t &a, const vec_t &b ) {
        return _mm_max_epi16( a, b );
    }
    
    static inline const vec_t abs_diff( const vec_t &a, const vec_t &b ) {
        // i don't really like this function, as ideally this would just be abs(sub(a,b)), 
        // but there doesn't seem to be a fast way to implement abs on pre SSSSSSE3.
        // The max(sub(a,b),sub(b,a)) work-around seems to be the next-best thing in this special case.
        
        #ifdef __SSSE3__
        return _mm_abs_epi16(sub(a,b));
        #else
        // FIXME: is there a faster method for boring old CPUs?
        return max( sub(a,b), sub(b,a) );
//         #error missing SSSSSSSSSE3
        #endif
    }

    static inline void assert_alignment( T * p ) {
        assert( size_t(p) % 32 == 0 );
    }
    
};



template<>
struct vector_unit<int, 4> {

    const static bool do_checks = false;
    
    typedef __m128i vec_t;
    typedef int T;
    const static T POS_MAX_VALUE = 0x7fffffff;
    const static T LARGE_VALUE = 2100000000;
    const static T SMALL_VALUE = -2100000000;
    const static T BIAS = 0;
    const static size_t W = 4;
    
    static inline vec_t setzero() {
        return set1(0);
    }
    
    static inline vec_t set1( T val ) {
        return _mm_set1_epi32( val );
    }
    
    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";
        
        _mm_store_si128( (vec_t*)addr, v );
    }
    
    static inline const vec_t load( const T* addr ) {
        return _mm_load_si128( (vec_t*)addr );
    }
    
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        return _mm_and_si128( a, b );
    }
    
    static inline const vec_t bit_andnot( const vec_t &a, const vec_t &b ) {
        return _mm_andnot_si128( a, b );
    }
    
//     static inline const vec_t bit_invert( const vec_t &a ) {
//         //return _mm_andnot_pd(a, set1(0xffff));
//     }
    
    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        return _mm_add_epi32( a, b );
    }
    static inline const vec_t adds( const vec_t &a, const vec_t &b ) {
        // there is no saturating add for 32bit int. Just hope that nothing bad will happen.
        // Treat saturation as kind of 'best effort hint'
        return _mm_add_epi32( a, b );
    }
    
    static inline const vec_t sub( const vec_t &a, const vec_t &b ) {
        return _mm_sub_epi32( a, b );
    }
    static inline const vec_t cmp_zero( const vec_t &a ) {
        return _mm_cmpeq_epi32( a, setzero() );
    }
    
    static inline const vec_t cmp_eq( const vec_t &a, const vec_t &b ) {
        return _mm_cmpeq_epi32( a, b );
    }
    
    static inline const vec_t cmp_lt( const vec_t &a, const vec_t &b ) {
     
        return _mm_cmplt_epi32( a, b );
    }
    
    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        // sse 4.1, no shit! what were they smoking...
#ifdef __SSE4_1__
        return _mm_min_epi32( a, b );
#else
        throw std::runtime_error( "missing sse4.1" );
        
// #error missing SSE4.1, find some workaround...
#endif
        
    }
    
    static inline const vec_t max( const vec_t &a, const vec_t &b ) {
#ifdef __SSE4_1__
        return _mm_max_epi32( a, b );
#else      
        throw std::runtime_error( "missing sse4.1" );
// #error missing SSE4.1, find some workaround...
#endif
    }
    
    static inline const vec_t abs_diff( const vec_t &a, const vec_t &b ) {
        // i don't really like this function, as ideally this would just be abs(sub(a,b)), 
        // but there doesn't seem to be a fast way to implement abs on pre SSSSSSE3.
        // The max(sub(a,b),sub(b,a)) work-around seems to be the next-best thing in this special case.
        
        #ifdef __SSSE3__
        return _mm_abs_epi32(sub(a,b));
        #else
        // FIXME: is there a faster method for boring old CPUs?
        return max( sub(a,b), sub(b,a) );
//         #error missing SSSSSSSSSE3
        #endif
    }
    
    static inline void assert_alignment( T * p ) {
        assert( size_t(p) % 32 == 0 );
    }

};
#ifndef WIN32

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu"
#endif

template<>
struct vector_unit<float, 4> {

    const static bool do_checks = false;
    
    typedef __m128 vec_t;
    typedef float T;

    const static int SIGN_MASK_INT = 0x7FFFFFFF;
    
    const static T LARGE_VALUE  = 1e8;
    const static T SMALL_VALUE = -1e8;
    const static T BIAS = 0;
    const static size_t W = 4;
    
    static inline vec_t cast_from_int( const __m128i &iv ) {
    	return _mm_castsi128_ps( iv );
    }

    static inline vec_t setzero() {
        return set1(0);
    }
    
    static inline vec_t set1( T val ) {
        return _mm_set1_ps( val );
    }
    
    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";
        
        _mm_store_ps( (T*)addr, v );
    }
    
    static inline const vec_t load( const T* addr ) {
        return _mm_load_ps( (T*)addr );
    }
    
#if 1 
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        return _mm_and_ps( a, b );
    }
    
    static inline const vec_t bit_andnot( const vec_t &a, const vec_t &b ) {
        return _mm_andnot_ps( a, b );
    }
#endif
//     static inline const vec_t bit_invert( const vec_t &a ) {
//         //return _mm_andnot_pd(a, set1(0xffff));
//     }
    
    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        return _mm_add_ps( a, b );
    }
    static inline const vec_t mul( const vec_t &a, const vec_t &b ) {
		return _mm_mul_ps( a, b );
	}

    static inline const vec_t adds( const vec_t &a, const vec_t &b ) {
        // float add is always saturating, kind of!?
        return _mm_add_ps( a, b );
    }
    static inline const vec_t sub( const vec_t &a, const vec_t &b ) {
        return _mm_sub_ps( a, b );
    }


    static inline const vec_t cmp_zero( const vec_t &a ) {
        return _mm_cmpeq_ps( a, setzero() );
    }
    
    static inline const vec_t cmp_eq( const vec_t &a, const vec_t &b ) {
        // TODO: think about what this really means. Maybe use something epsilon-based as a default for float?
        return _mm_cmpeq_ps( a, b );
    }
    
    static inline const vec_t cmp_lt( const vec_t &a, const vec_t &b ) {
     
        return _mm_cmplt_ps( a, b );
    }
    
    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        return _mm_min_ps( a, b );
    }
    
    static inline const vec_t max( const vec_t &a, const vec_t &b ) {
        return _mm_max_ps( a, b );
    }
    
    static inline const vec_t abs_diff( const vec_t &a, const vec_t &b ) {
        // i don't really like this function, as ideally this would just be abs(sub(a,b)), 
        // but there doesn't seem to be a fast way to implement abs on pre SSSSSSE3.
        // The max(sub(a,b),sub(b,a)) work-around seems to be the next-best thing in this special case.
        
        //unsigned int SIGN_MASK[4] = {0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF};
        //unsigned int SIGN_MASK = 0x7FFFFFFF;
        const float *SIGN_MASK_PTR = (float*)&SIGN_MASK_INT;
        static float SIGN_MASK = *SIGN_MASK_PTR;
        return bit_and(sub(a,b), set1(SIGN_MASK) ); // TODO: could this case any alignment problems?
        //return bit_and(sub(a,b), set1(*((float*)&SIGN_MASK_INT) )); // TODO: could this case any alignment problems?
        
        //return bit_and(sub(a,b), load((float*)SIGN_MASK ));
        
    }
    

};

template<>
struct vector_unit<double, 2> {

    const static bool do_checks = false;

    typedef __m128d vec_t;
    typedef double T;

//    const static uint64_t SIGN_MASK_U64 = 0x7FFFFFFFFFFFFFFF;

    const static T LARGE_VALUE  = 1e8;
    const static T SMALL_VALUE = -1e8;
    const static T BIAS = 0;
    const static size_t W = 2;

    static inline vec_t cast_from_int( const __m128i &iv ) {
        return _mm_castsi128_pd( iv );
    }

    static inline vec_t setzero() {
        return set1(0);
    }

    static inline vec_t set1( T val ) {
        return _mm_set1_pd( val );
    }

    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";

        _mm_store_pd( (T*)addr, v );
    }

    static inline const vec_t load( const T* addr ) {
        return _mm_load_pd( (T*)addr );
    }

#if 1
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        return _mm_and_pd( a, b );
    }

    static inline const vec_t bit_or( const vec_t &a, const vec_t &b ) {
        return _mm_or_pd( a, b );
    }

    static inline const vec_t bit_andnot( const vec_t &a, const vec_t &b ) {
        return _mm_andnot_pd( a, b );
    }
#endif
//     static inline const vec_t bit_invert( const vec_t &a ) {
//         //return _mm_andnot_pd(a, set1(0xffff));
//     }

    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        return _mm_add_pd( a, b );
    }
    static inline const vec_t mul( const vec_t &a, const vec_t &b ) {
        return _mm_mul_pd( a, b );
    }

    static inline const vec_t adds( const vec_t &a, const vec_t &b ) {
        // float add is always saturating, kind of!?
        return _mm_add_pd( a, b );
    }
    static inline const vec_t sub( const vec_t &a, const vec_t &b ) {
        return _mm_sub_pd( a, b );
    }


    static inline const vec_t cmp_zero( const vec_t &a ) {
        return _mm_cmpeq_pd( a, setzero() );
    }

    static inline const vec_t cmp_eq( const vec_t &a, const vec_t &b ) {
        // TODO: think about what this really means. Maybe use something epsilon-based as a default for float?
        return _mm_cmpeq_pd( a, b );
    }

    static inline const vec_t cmp_lt( const vec_t &a, const vec_t &b ) {

        return _mm_cmplt_pd( a, b );
    }

    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        return _mm_min_pd( a, b );
    }

    static inline const vec_t max( const vec_t &a, const vec_t &b ) {
        return _mm_max_pd( a, b );
    }

    static inline const vec_t abs_diff( const vec_t &a, const vec_t &b ) {
        // i don't really like this function, as ideally this would just be abs(sub(a,b)),
        // but there doesn't seem to be a fast way to implement abs on pre SSSSSSE3.
        // The max(sub(a,b),sub(b,a)) work-around seems to be the next-best thing in this special case.

        //unsigned int SIGN_MASK[4] = {0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF};
        //unsigned int SIGN_MASK = 0x7FFFFFFF;

        const uint64_t SIGN_MASK_U64x = 0x7fffffffffffffff;

        const double *SIGN_MASK_PTR = (double*)&SIGN_MASK_U64x;
        double SIGN_MASK = *SIGN_MASK_PTR;
        return bit_and(sub(a,b), set1(SIGN_MASK) ); // TODO: could this case any alignment problems?
        //return bit_and(sub(a,b), set1(*((float*)&SIGN_MASK_INT) )); // TODO: could this case any alignment problems?

        //return bit_and(sub(a,b), load((float*)SIGN_MASK ));


    }

    static inline const vec_t abs( const vec_t &a ) {
//        const uint64_t SIGN_MASK_U64x = 0x7fffffffffffffff;
//
//        const static double *SIGN_MASK_PTR = (double*)&SIGN_MASK_U64x;
//        static double SIGN_MASK = *SIGN_MASK_PTR;
//
//        return bit_and( a, set1(SIGN_MASK));
//
//
        return max( a, sub( setzero(), a ));
    }

    static inline void println( const vec_t & v ) {

        double tmp[2];
        _mm_storeu_pd( tmp, v );
        printf( "%f %f\n", tmp[0], tmp[1] );
    }


};


#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
#ifdef HAVE_AVX
// AVX is the most pointless thing in the world, as far as integers are concerned.
// waiting for AVX5.7

// AVX 16x16bit vector unit
template<>
struct vector_unit<short, 16> {

    const static bool do_checks = false;
    
    typedef __m256i vec_t;
    typedef short T;
    
    const static T SMALL_VALUE = -32000;
    const static T BIAS = 0;
    const static size_t W = 8;
    
    static inline vec_t setzero() {
        return _mm256_setzero_si256();
    }
    
    static inline vec_t set1( T val ) {
        return _mm256_set1_epi16( val );
    }
    
    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";
        
        _mm256_store_si256( (vec_t*)addr, v );
    }
    
    static inline const vec_t load( T* addr ) {
        return _mm256_load_si256( (vec_t*)addr );
    }
    
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        __m128i lowa = _mm256_castsi256_si128(a);
        const __m128i lowb = _mm256_castsi256_si128(b);
        
        lowa = _mm_and_si128( lowa, lowb );
        
        
        __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        higha = _mm_and_si128( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(lowa), higha, 1 );
    }
    
    static inline const vec_t bit_andnot( const vec_t &a, const vec_t &b ) {
        __m128i lowa = _mm256_castsi256_si128(a);
        const __m128i lowb = _mm256_castsi256_si128(b);
        
        lowa = _mm_andnot_si128( lowa, lowb );
        
        
        __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        higha = _mm_andnot_si128( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(lowa), higha, 1 );
    }
    
//     static inline const vec_t bit_invert( const vec_t &a ) {
//         //return _mm_andnot_pd(a, set1(0xffff));
//     }
    
    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        const __m128i lowa = _mm256_extractf128_si256(a,0);
        const __m128i lowb = _mm256_extractf128_si256(b,0);
        
        //lowa = _mm_add_epi16( lowa, lowb );
        
        
        const __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        //higha = _mm_add_epi16( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(_mm_add_epi16( lowa, lowb )),  _mm_add_epi16( higha, highb ), 1 );
    }
    static inline const vec_t sub( const vec_t &a, const vec_t &b ) {
        __m128i lowa = _mm256_castsi256_si128(a);
        const __m128i lowb = _mm256_castsi256_si128(b);
        
        lowa = _mm_sub_epi16( lowa, lowb );
        
        
        __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        higha = _mm_sub_epi16( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(lowa), higha, 1 );
    }
    static inline const vec_t cmp_zero( const vec_t &a ) {
        return cmp_eq( a, setzero() );
    }
    
    static inline const vec_t cmp_eq( const vec_t &a, const vec_t &b ) {
        __m128i lowa = _mm256_castsi256_si128(a);
        const __m128i lowb = _mm256_castsi256_si128(b);
        
        lowa = _mm_cmpeq_epi16( lowa, lowb );
        
        
        __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        higha = _mm_cmpeq_epi16( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(lowa), higha, 1 );
    }
    
    static inline const vec_t cmp_lt( const vec_t &a, const vec_t &b ) {
     
        __m128i lowa = _mm256_castsi256_si128(a);
        const __m128i lowb = _mm256_castsi256_si128(b);
        
        lowa = _mm_cmplt_epi16( lowa, lowb );
        
        
        __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        higha = _mm_cmplt_epi16( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(lowa), higha, 1 );
    }
    
    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        __m128i lowa = _mm256_castsi256_si128(a);
        const __m128i lowb = _mm256_castsi256_si128(b);
        
        lowa = _mm_min_epi16( lowa, lowb );
        
        
        __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        higha = _mm_min_epi16( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(lowa), higha, 1 );
    }
    
    static inline const vec_t max( const vec_t &a, const vec_t &b ) {
        __m128i lowa = _mm256_castsi256_si128(a);
        const __m128i lowb = _mm256_castsi256_si128(b);
        
        lowa = _mm_max_epi16( lowa, lowb );
        
        
        __m128i higha = _mm256_extractf128_si256(a,1);
        const __m128i highb = _mm256_extractf128_si256(b,1);
        
        higha = _mm_max_epi16( higha, highb );
        
        
        return _mm256_insertf128_si256( _mm256_castsi128_si256(lowa), higha, 1 );
    }
};



template<>
struct vector_unit<double, 4> {

    const static bool do_checks = false;

    typedef __m256d vec_t;
    typedef double T;

//    const static uint64_t SIGN_MASK_U64 = 0x7FFFFFFFFFFFFFFF;

    const static T LARGE_VALUE  = 1e8;
    const static T SMALL_VALUE = -1e8;
    const static T BIAS = 0;
    const static size_t W = 2;


    static inline vec_t setzero() {
        return _mm256_setzero_pd();
    }

    static inline vec_t set1( T val ) {
        return _mm256_set1_pd( val );
    }

    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";

        _mm256_store_pd( (T*)addr, v );
    }

    static inline const vec_t load( const T* addr ) {
        return _mm256_load_pd( (T*)addr );
    }

#if 1
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        return _mm256_and_pd( a, b );
    }

    static inline const vec_t bit_or( const vec_t &a, const vec_t &b ) {
        return _mm256_or_pd( a, b );
    }

    static inline const vec_t bit_andnot( const vec_t &a, const vec_t &b ) {
        return _mm256_andnot_pd( a, b );
    }
#endif
//     static inline const vec_t bit_invert( const vec_t &a ) {
//         //return _mm_andnot_pd(a, set1(0xffff));
//     }

    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        return _mm256_add_pd( a, b );
    }
    static inline const vec_t mul( const vec_t &a, const vec_t &b ) {
        return _mm256_mul_pd( a, b );
    }

    static inline const vec_t adds( const vec_t &a, const vec_t &b ) {
        // float add is always saturating, kind of!?
        return _mm256_add_pd( a, b );
    }
    static inline const vec_t sub( const vec_t &a, const vec_t &b ) {
        return _mm256_sub_pd( a, b );
    }


    static inline const vec_t cmp_zero( const vec_t &a ) {
        return _mm256_cmp_pd( a, setzero(), _CMP_EQ_OQ );
    }

    static inline const vec_t cmp_eq( const vec_t &a, const vec_t &b ) {
        // TODO: think about what this really means. Maybe use something epsilon-based as a default for float?
        return _mm256_cmp_pd( a, b, _CMP_EQ_OQ );
    }

    static inline const vec_t cmp_lt( const vec_t &a, const vec_t &b ) {

        return _mm256_cmp_pd( a, b, _CMP_LT_OQ );
    }

    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        return _mm256_min_pd( a, b );
    }

    static inline const vec_t max( const vec_t &a, const vec_t &b ) {
        return _mm256_max_pd( a, b );
    }

    static inline const vec_t abs_diff( const vec_t &a, const vec_t &b ) {
        // i don't really like this function, as ideally this would just be abs(sub(a,b)),
        // but there doesn't seem to be a fast way to implement abs on pre SSSSSSE3.
        // The max(sub(a,b),sub(b,a)) work-around seems to be the next-best thing in this special case.

        //unsigned int SIGN_MASK[4] = {0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF};
        //unsigned int SIGN_MASK = 0x7FFFFFFF;

        const uint64_t SIGN_MASK_U64x = 0x7fffffffffffffff;

        const double *SIGN_MASK_PTR = (double*)&SIGN_MASK_U64x;
        double SIGN_MASK = *SIGN_MASK_PTR;
        return bit_and(sub(a,b), set1(SIGN_MASK) ); // TODO: could this case any alignment problems?
        //return bit_and(sub(a,b), set1(*((float*)&SIGN_MASK_INT) )); // TODO: could this case any alignment problems?

        //return bit_and(sub(a,b), load((float*)SIGN_MASK ));


    }

    static inline const vec_t abs( const vec_t &a ) {
//        const uint64_t SIGN_MASK_U64x = 0x7fffffffffffffff;
//
//        const static double *SIGN_MASK_PTR = (double*)&SIGN_MASK_U64x;
//        static double SIGN_MASK = *SIGN_MASK_PTR;
//
//        return bit_and( a, set1(SIGN_MASK));
//
//
        return max( a, sub( setzero(), a ));
    }

    static inline void println( const vec_t & v ) {

        double tmp[4];
        _mm256_storeu_pd( tmp, v );
        printf( "%f %f %f %f\n", tmp[0], tmp[1], tmp[2], tmp[3] );
    }


};


#endif


// vector unit specialization: SSE 16x8bit integer 

template<>
struct vector_unit<unsigned char, 16> {

    const static bool do_checks = false;
    
    typedef __m128i vec_t;
    typedef unsigned char T;
    
    const static size_t W = 16;
    const static T SMALL_VALUE = 1;
    const static T BIAS = 12;
    
    
    static inline vec_t setzero() {
        return set1(0);
    }
    
    static inline vec_t set1( T val ) {
        return _mm_set1_epi8( val );
    }
    
    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";
        
        _mm_store_si128( (__m128i*)addr, v );
    }
    
    static inline const vec_t load( T* addr ) {
        return (vec_t)_mm_load_si128( (__m128i *) addr );
    }
    
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        return _mm_and_si128( a, b );
    }
    
    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        return _mm_add_epi8( a, b );
    }
    
    static inline const vec_t sub( const vec_t &a, const vec_t &b ) {
        return _mm_sub_epi8( a, b );
    }
    
    static inline const vec_t cmp_zero( const vec_t &a ) {
        return _mm_cmpeq_epi8( a, setzero() );
    }
    
    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        return _mm_min_epu8( a, b );
    }
    static inline const vec_t max( const vec_t &a, const vec_t &b ) {
        return _mm_max_epu8( a, b );
    }
    
    static inline void println( const vec_t &v, T *tmp ) {
        store(v, tmp);
        std::cout << "(";
        for( size_t i = 0; i < W; i++ ) {
         
            std::cout << int(tmp[i]) << ((i < W-1) ? "," : ")");
        }
        
        std::cout << std::endl;
    }
};


template<>
struct vector_unit<short, 1> {

    const static bool do_checks = false;
    
    typedef __m128i vec_t;
    typedef short T;
    
    const static size_t W = 1;
    
    static inline vec_t setzero() {
        return set1(0);
    }
    
    static inline vec_t set1( T val ) {
        return _mm_set1_epi16( val );
    }
    
    static inline void store( const vec_t &v, T *addr ) {

        if( do_checks && addr == 0 ) {
            throw std::runtime_error( "store: addr == 0" );
        }
        //std::cout << "store to: " << addr << "\n";
        
        //_mm_store_pd( (float*)addr, (__m128)v );
        
        *addr = _mm_extract_epi16( v, 0 );
    }
    
    static inline const vec_t load( T* addr ) {
        //return (vec_t)_mm_load_pd( (float *) addr );
        return set1(*addr);
    }
    
    static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
        return _mm_and_si128( a, b );
    }
    
    static inline const vec_t add( const vec_t &a, const vec_t &b ) {
        return _mm_add_epi16( a, b );
    }
    
    static inline const vec_t cmp_zero( const vec_t &a ) {
        return _mm_cmpeq_epi16( a, setzero() );
    }
    
    static inline const vec_t min( const vec_t &a, const vec_t &b ) {
        return _mm_min_epi16( a, b );
    }
};


// vector unit specialization: SSE scalar 16bit integer 

// template<>
// struct vector_unit<short, 1> {
// 
//     const static bool do_checks = false;
//     
//     typedef short vec_t;
//     typedef short T;
//     
//     const static size_t W = 1;
//     
//     static inline vec_t setzero() {
//         return 0;
//     }
//     
//     static inline vec_t set1( T val ) {
//         return val;
//     }
//     
//     static inline void store( const vec_t &v, T *addr ) {
// 
//         if( do_checks && addr == 0 ) {
//             throw std::runtime_error( "store: addr == 0" );
//         }
//         
//         *addr = v;
//         
//     }
//     
//     static inline const vec_t load( T* addr ) {
//         return *addr;
//     }
//     
//     static inline const vec_t bit_and( const vec_t &a, const vec_t &b ) {
//         return  a & b;
//     }
//     
//     static inline const vec_t add( const vec_t &a, const vec_t &b ) {
//         return a + b;
//     }
//     
//     static inline const vec_t cmp_zero( const vec_t &a ) {
//         //return _mm_cmpeq_epi16( a, setzero() );
//         return a == 0 ? 0xffff : 0x0;
//     }
//     
//     static inline const vec_t min( const vec_t &a, const vec_t &b ) {
//         return std::min( a, b );
//     }
// };
#endif
