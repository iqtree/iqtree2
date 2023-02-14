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
#include <string.h> //for strrchr (needed for GCC 9.2 builds)

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

inline void _my_assert(const std::string &str_expression, const char *func, const char* file, int line) {
    _my_assert(str_expression.c_str(), func, file, line);
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

#if defined(__GNUC__) || defined(__clang__)
    #define FUNCTION_NOT_IMPLEMENTED \
        _my_assert(std::string(__PRETTY_FUNCTION__) + " is not implemented", \
                  __PRETTY_FUNCTION__, __FILE__, __LINE__)
#else
    #define FUNCTION_NOT_IMPLEMENTED \
        _my_assert(std::string(__func__) + " is not implemented", \
                  __func__, __FILE__, __LINE__)
#endif

#endif /* assert_h */
