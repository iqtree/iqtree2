/*
 * Copyright (C) 2011 Simon A. Berger
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
#ifndef __aligned_buffer_h
#define __aligned_buffer_h


#include <cstdlib>
#include <cstddef>
#include <cassert>
#include <stdexcept>

#if 1

#include <vector>

namespace ab_internal_ {
template<typename T, size_t Talign>
class alloc {
    //const static size_t align = 4096;
    
#ifndef WIN32
    struct allocator_posix {
        static inline void *alloc( size_t align, size_t size ) {
            void *ptr;
            int ret = posix_memalign( (void**)&ptr, align, size );
            
            if( ret != 0 ) {
                throw std::runtime_error( "posix_memalign failed" );
            }
            return ptr;
        }

        static inline void free( void *ptr ) {
            std::free( ptr );
        }
    };
    typedef allocator_posix allocator;
#endif
#ifdef WIN32
    struct allocator_ugly {
        static inline void *alloc( size_t align, size_t size ) {
            return _aligned_malloc( size, align );
        }

        static inline void free( void *ptr ) {
            _aligned_free(ptr);
        }
    };
    typedef allocator_ugly allocator;
#endif
public:
    template<typename _Other>
    struct rebind {
        typedef ab_internal_::alloc<_Other,Talign> other;
    };
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;


	alloc( ) {}
	alloc(const alloc<T,Talign>& _Right ) {}

	template<class Other>
	alloc(const alloc<Other,Talign>& _Right ) {}

    pointer allocate( size_type nobj, const void *lh = 0 ) {
        return (pointer) allocator::alloc( Talign, nobj * sizeof(T) );
    }
    
    void deallocate( pointer ptr, size_type nobj ) {
        allocator::free( ptr );
    }
    
    
    void construct( pointer p, const_reference t) { new ((void*) p) T(t); }
    
    
    void destroy( pointer p ){ ((T*)p)->~T(); }
    
    size_type max_size() const {
        return size_t(-1);
    }
};
    
}

template<typename T, size_t alignment = 4096>
struct aligned_buffer : private std::vector<T,ab_internal_::alloc<T,alignment> > {
public:
    typedef typename std::vector<T,ab_internal_::alloc<T,alignment> >::iterator iterator;
    typedef typename std::vector<T,ab_internal_::alloc<T,alignment> >::const_iterator const_iterator;
    
    aligned_buffer() : std::vector<T,ab_internal_::alloc<T,alignment> >() {}
    aligned_buffer( size_t size ) : std::vector<T,ab_internal_::alloc<T,alignment> >(size) {}
    aligned_buffer( size_t size, const T &v ) : std::vector<T,ab_internal_::alloc<T,alignment> >(size, v) {}
    

    using std::vector<T,ab_internal_::alloc<T,alignment> >::begin;
    using std::vector<T,ab_internal_::alloc<T,alignment> >::end;
    using std::vector<T,ab_internal_::alloc<T,alignment> >::size;
    using std::vector<T,ab_internal_::alloc<T,alignment> >::resize;
    using std::vector<T,ab_internal_::alloc<T,alignment> >::reserve;
    using std::vector<T,ab_internal_::alloc<T,alignment> >::push_back;
    using std::vector<T,ab_internal_::alloc<T,alignment> >::data;
    using std::vector<T,ab_internal_::alloc<T,alignment> >::empty;

    
    using std::vector<T,ab_internal_::alloc<T,alignment> >::operator[];
    
    inline T* operator() (ptrdiff_t o) {
        return &(operator[](o));
    }
    
    inline const T* operator() (ptrdiff_t o) const {
        return &(operator[](o));
    }
    inline T* base() {
        return operator()(0);
    }
};


#else
// template<typename T>
// struct aligned_buffer {
//     typedef T* iterator;
//     
//     T* m_ptr;
//     size_t m_size;
//     const static size_t align = 32;
//     
// #ifndef WIN32
// 	struct allocator_posix {
// 		static inline void *alloc( size_t align, size_t size ) {
// 			void *ptr;
// 			int ret = posix_memalign( (void**)&ptr, align, size );
//             
//             if( ret != 0 ) {
//                 throw std::runtime_error( "posix_memalign failed" );
//             }
// 			return ptr;
// 		}
// 
// 		static inline void free( void *ptr ) {
// 			std::free( ptr );
// 		}
// 	};
// 	typedef allocator_posix allocator;
// #endif
// #ifdef WIN32
// 	struct allocator_ugly {
// 		static inline void *alloc( size_t align, size_t size ) {
// 			return _aligned_malloc( size, align );
// 		}
// 
// 		static inline void free( void *ptr ) {
// 			_aligned_free(ptr);
// 		}
// 	};
// 	typedef allocator_ugly allocator;
// #endif
// 
//     aligned_buffer() : m_ptr(0), m_size(0) {}
//     
//     aligned_buffer( size_t size ) : m_ptr(0), m_size(0) {
//          
//         resize( size );
//                 
//     }
//     
//     void resize( size_t ns ) {
//     
//         
//         if( ns != size() ) {
//             allocator::free( m_ptr );
//             
//             m_size = ns;
// 			m_ptr = (T*)allocator::alloc( align, byte_size() );
// 
//             
//         }
//         
//     }
//     
//     size_t size() const {
//         return m_size;
//     }
//     
//     size_t byte_size() {
//         return m_size * sizeof(T);
//     }
//     
//     ~aligned_buffer() {
//         allocator::free( m_ptr );
//     }
//     
//     T *begin() const {
//         return m_ptr;
//     }
//     
//     T* end() const {
//         return begin() + m_size;
//     }
//     
//     inline T* operator() (ptrdiff_t o) const {
//         assert( o < m_size );
//         return begin() + o;
//     }
//     
//     inline T& operator[](ptrdiff_t o) {
//         assert( o < m_size );
//         return *(begin() + o);
//    }
//     
//     aligned_buffer &operator=( const aligned_buffer &other ) {
//         resize(other.size());
//         std::copy( other.begin(), other.end(), begin() );
//         
//         return *this;
//     }
//     
//     aligned_buffer( const aligned_buffer &other ) : m_ptr(0), m_size(0) {
//         resize(other.size());
//         std::copy( other.begin(), other.end(), begin() );
//     }
//         
//         
//     T* base() {
//         return m_ptr;
//     }
//     
// };
#endif
#endif
