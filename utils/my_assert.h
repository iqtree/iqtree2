//
//  assert.h
//  iqtree
//
//  Created by James Barbetti on 17/12/20.
//

#ifndef assert_h
#define assert_h

#include <stdlib.h>
#include <iostream>

// redefine assertion
inline void _my_assert(const char* expression, const char *func, const char* file, int line)
{
    const char *sfile = strrchr(file, '/');
    if (!sfile) {
        sfile = file;
    }
    else {
        ++sfile;
    }
    std::cerr << sfile << ":" << line << ": " << func
        << ": Assertion `" << expression << "' failed." << std::endl;
    abort();
}
 
#ifdef NDEBUG
#define ASSERT(EXPRESSION) ((void)0)
#else
    #if defined(__GNUC__) || defined(__clang__)
        #define ASSERT(EXPRESSION) ((EXPRESSION) ? (void)0 : _my_assert(#EXPRESSION, __PRETTY_FUNCTION__, __FILE__, __LINE__))
    #else
        #define ASSERT(EXPRESSION) ((EXPRESSION) ? (void)0 : _my_assert(#EXPRESSION, __func__, __FILE__, __LINE__))
    #endif
#endif

#endif /* assert_h */
