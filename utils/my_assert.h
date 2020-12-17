//
//  my_assert.h
//  iqtree
//
//_my_assert function, originally created by
//BUI Quang Minh <minh.bui@univie.ac.at> 21-Aug-2016.
//This file: created by James Barbetti on 17-Dec-2020.
//(decentTree needs _my_assert and ASSERT, but
// doesn't need anything from the rest of tools.h)
//

#ifndef assert_h
#define assert_h

#include <stdlib.h>
#include <iostream>

// for MSVC
#ifndef __func__
#define __func__ __FUNCTION__
#endif

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
