##################################################################
# IQ-TREE cmake build definition
# Copyright (c) 2012-2015 Bui Quang Minh, Lam-Tung Nguyen
##################################################################

# Windows example usages:
#------------------------
# cmake -G "Visual Studio 12" <source_dir>          (32-bit version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" <source_dir>    (64-bit version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" -T "Intel C++ Compiler XE 15.0" <source_dir>        (64-bit version, compiled with ICC)
# cmake -G "MinGW Makefiles" <source_dir> (TDM-GCC)
# cmake -G "Unix Makefiles" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_MAKE_PROGRAM=mingw32-make <source_dir> (TDM-GCC and clang)

# Linux example usages:
#----------------------
# cmake <source_dir>                           (sequential version)
# cmake <source_dir>      (OpenMP version)
# cmake -DIQTREE_FLAGS="m32-single" <source_dir>      (32-bit sequential version)
# cmake -DIQTREE_FLAGS="m32" <source_dir>  (32-bit OpenMP version)
#
# To compile with CLANG on Linux:
# export CC=/usr/bin/clang
# export CXX=/usr/bin/clang++
# Best practices for setting up CMAKE for diffrent compiler can be found here:
# http://stackoverflow.com/questions/7031126/switching-between-gcc-and-clang-llvm-using-cmake
#
# Mac OSX example usages:
#------------------------
#
# To build OpenMP version one needs to download Clang version 3.7 or later (as of November 2015)
# Then assuming clang3.7 and clang++3.7 are the newly built compilers, then:
# cmake -DCMAKE_C_COMPILER=clang3.7 -DCMAKE_CXX_COMPILER=clang++3.7 <source_dir>   (OpenMP version)
#
# cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DIQTREE_FLAGS=single <source_dir>   (sequential version)
#


# Compile OpenMP version: cmake ....
# Compile 32-bit version: cmake -DIQTREE_FLAGS=m32 ....
# Compile static version: cmake -DIQTREE_FLAGS=static ....
# Compile static OpenMP version: cmake -DIQTREE_FLAGS="static" ....

#NOTE: Static linking with clang windows: make a symlink libgcc_eh.a to libgcc.a (administrator required)
# C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0>mklink libgcc_eh.a libgcc.a

cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

set(GCC_MIN_VERSION "4.8")

project(iqtree)
add_definitions(-DIQ_TREE)

# Find Eigen3 library
if (NOT EIGEN3_INCLUDE_DIR)
    find_package(Eigen3)
    if(NOT EIGEN3_FOUND)
        message(FATAL_ERROR "Eigen3 library not found. Either install it or rerun cmake with -DEIGEN3_INCLUDE_DIR=<installed_eigen3_dir>")
    endif()
endif()
add_definitions(-I${EIGEN3_INCLUDE_DIR})

# Find Boost library
find_package(Boost)
if(Boost_FOUND)
    add_definitions(-I${Boost_INCLUDE_DIRS})
    add_definitions(-DUSE_BOOST)
endif()


# The version number.
set (iqtree_VERSION_MAJOR 1)
set (iqtree_VERSION_MINOR 6)
set (iqtree_VERSION_PATCH "6")

set(BUILD_SHARED_LIBS OFF)

if (CMAKE_C_COMPILER MATCHES "mpic")
    set(IQTREE_FLAGS "${IQTREE_FLAGS} mpi")
endif()

message("IQ-TREE flags : ${IQTREE_FLAGS}")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message("Builde mode   : Release")
endif()

if (CMAKE_GENERATOR MATCHES "Xcode")
    set(CMAKE_XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT "dwarf-with-dsym")
    set(CMAKE_XCODE_ATTRIBUTE_COMPILER_INDEX_STORE_ENABLE "No")
endif()

include_directories("${PROJECT_SOURCE_DIR}")

##################################################################
# Include the Terraphast library
##################################################################
option(USE_TERRAPHAST "Use phylogentic terraces library (terraphast)" ON)

if (USE_TERRAPHAST)
   set(GCC_MIN_VERSION "5.4")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIQTREE_TERRAPHAST")

   option(TERRAPHAST_USE_GMP "" OFF)
   option(TERRAPHAST_BUILD_CLIB "" OFF)
   option(TERRAPHAST_BUILD_APPS "" OFF)
   option(TERRAPHAST_BUILD_TESTS "" OFF)
endif()

##################################################################
# Detect target platforms
##################################################################
if (WIN32)
    message("Target OS     : Windows")
    # build as static binary to run on most machines
    if (IQTREE_FLAGS MATCHES "static")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    add_definitions(-DWIN32)
elseif (APPLE)
    message("Target OS     : Mac OS X")
    # to be compatible back to Mac OS X 10.7
    if (IQTREE_FLAGS MATCHES "oldmac")
        add_definitions("-mmacosx-version-min=10.5")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.5")
    else()
        add_definitions("-mmacosx-version-min=10.7")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.7")
    endif()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif (UNIX)
    message("Target OS     : Unix")
    # build as static binary to run on most machines
    if (NOT IQTREE_FLAGS MATCHES "static")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif()
else()
    # Note that IQ-TREE has NOT been tested on other platforms
    message("Target OS     : Unknown and untested yet")
endif()

##################################################################
# Setup compiler, currently supported GCC, CLANG, MSVC, and ICC
##################################################################

set(GCC "FALSE")   #  GNU compiler
set(CLANG "FALSE") # Clang compiler
set(ICC "FALSE")   # Intel compiler
set(VCC "FALSE")   # MS Visual C Compiler, note that it is different from MSVC variable
# using C++11 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (CMAKE_COMPILER_IS_GNUCXX)
    message("Compiler      : GNU Compiler (gcc)")
    set(GCC "TRUE")
#   set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -fms-extensions -Wno-deprecated")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -ffunction-sections -fdata-sections")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -g -ffunction-sections -fdata-sections")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-dead_strip")
    else()
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,--gc-sections")
    endif()
    # require at least gcc ${GCC_MIN_VERSION}
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_MIN_VERSION)
        message(FATAL_ERROR "GCC version must be at least ${GCC_MIN_VERSION}!")
    endif()
    if (WIN32)
        # disable AVX on Windows due to memory alignment
        set(IQTREE_FLAGS "${IQTREE_FLAGS} novx")
        message("WARNING: AVX is disabled on Windows as GCC does not properly suport memory alignment")
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Compiler      : Clang")
    set(CLANG "TRUE")
#   set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -Wno-nested-anon-types")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffunction-sections -fdata-sections")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ffunction-sections -fdata-sections")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,-dead_strip")
    else()
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Wl,--gc-sections")
    endif()

    # use libc++ per default in MacOS
    if (APPLE)
        SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()

    #remove -rdynamic for Clang under Linux
    if (UNIX AND IQTREE_FLAGS MATCHES "static")
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(VCC "TRUE")
    message("Compiler      : MS Visual C++ Compiler")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message("Compiler      : Intel C++ Compiler (icc)")
    set(ICC "TRUE")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
else()
    message("Compiler      : Unknown and untested yet")
endif()

set(EXE_SUFFIX "")

if (MSVC)
    # MS Visual Studio environment
    message("Exporting MS Visual Studio projects...")
    add_definitions(/MP) # enable multi-processor compilation
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_definitions(/Ot /Oi)
        if (VCC)
            add_definitions(/O2)
        elseif (ICC)
            #add_definitions(/O3)
            add_definitions(/O3)
        endif()
    endif()
endif()

# enable link time optimization
if (IQTREE_FLAGS MATCHES "lto")
    #if (CLANG)
    #    set(COMBINED_FLAGS "${COMBINED_FLAGS} -flto=thin")
    #else()
        set(COMBINED_FLAGS "${COMBINED_FLAGS} -flto")
    #endif()
endif()

##################################################################
# configure MPI compilation
##################################################################

if (IQTREE_FLAGS MATCHES "mpi")
    add_definitions(-D_IQTREE_MPI)
    if (NOT CMAKE_CXX_COMPILER MATCHES "mpi")
        # if not using the MPI compiler wrapper, set own options manually 
        find_package(MPI REQUIRED)
        set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS}")
        set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
        set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
        include_directories(${MPI_C_INCLUDE_PATH})
        include_directories(${MPI_CXX_INCLUDE_PATH})
    endif()
endif()


##################################################################
# Configure PLL build
##################################################################
if (IQTREE_FLAGS MATCHES "pll")
    add_definitions(-DUSING_PLL)
    set(EXE_SUFFIX "${EXE_SUFFIX}-pll")
endif()

##################################################################
# detect 32 or 64 bit binary
##################################################################
set (BINARY32 "FALSE")
if(CMAKE_SIZEOF_VOID_P EQUAL 4 OR IQTREE_FLAGS MATCHES "m32")
    set(BINARY32 "TRUE")
    message("Target binary : 32-bit")
    if (CMAKE_GENERATOR MATCHES "Win64")
        error("Both 32-bit and 64-bit mode cannot be specified")
    endif()
    #SET(EXE_SUFFIX "${EXE_SUFFIX}32")
    if (GCC OR CLANG) 
        set(COMBINED_FLAGS "${COMBINED_FLAGS} -m32")
    endif()
    add_definitions(-DBINARY32)
else()
    message("Target binary : 64-bit")
endif()

if(IQTREE_FLAGS MATCHES "novx")
    add_definitions(-D__NOAVX__)
endif()

##################################################################
# configure OpenMP/PThreads compilation
# change the executable name if compiled for OpenMP parallel version
##################################################################
if (NOT IQTREE_FLAGS MATCHES "single")
    message("OpenMP        : Yes")
    #SET(EXE_SUFFIX "${EXE_SUFFIX}-omp")
    add_definitions(-D_USE_PTHREADS)
    if (MSVC)
        add_definitions(/MT)
    endif()

    if (VCC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
        include_directories("${PROJECT_SOURCE_DIR}/pll") # for PThreads headers
    elseif (ICC)
        if (WIN32)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qopenmp")
            include_directories("${PROJECT_SOURCE_DIR}/pll") # for PThreads headers
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        endif()
    elseif (GCC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -pthread")
    elseif (CLANG)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp -pthread")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp=libomp")
    endif()
else()
    message("OpenMP        : NONE")
endif()


if (IQTREE_FLAGS MATCHES "mpi")
    message("MPI           : Yes")
    SET(EXE_SUFFIX "${EXE_SUFFIX}-mpi")
else()
    message("MPI           : NONE")
endif()

##################################################################
# configure SSE/AVX/FMA instructions
##################################################################

SET(AVX_FLAGS "-D__SSE3 -D__AVX")
if (VCC) 
    set(AVX_FLAGS "${AVX_FLAGS} /arch:AVX")
elseif (CLANG)
    set(AVX_FLAGS "${AVX_FLAGS} -mavx")
elseif (GCC)
    set(AVX_FLAGS "${AVX_FLAGS} -mavx -fabi-version=0")
elseif (ICC)
    if (WIN32)
         set(AVX_FLAGS "${AVX_FLAGS} /arch:avx")
    else()
         set(AVX_FLAGS "${AVX_FLAGS} -mavx")
    endif()
endif()

SET(SSE_FLAGS "-D__SSE3")
if (VCC)
    set(SSE_FLAGS "${SSE_FLAGS} /arch:SSE2 -D__SSE3__")
elseif (GCC OR CLANG)
    set(SSE_FLAGS "${SSE_FLAGS} -msse3")
elseif (ICC)
    if (WIN32)
        set(SSE_FLAGS "${SSE_FLAGS} /arch:sse3")
    else()
        set(SSE_FLAGS "${SSE_FLAGS} -msse3")
    endif()
endif()

SET(FMA_FLAGS "-D__SSE3 -D__AVX")
if (VCC) 
    set(FMA_FLAGS "${FMA_FLAGS} /arch:AVX2")
elseif (CLANG)
    set(FMA_FLAGS "${FMA_FLAGS} -mavx -mfma")
elseif (GCC)
    set(FMA_FLAGS "${FMA_FLAGS} -mavx -fabi-version=0 -mfma")
elseif (ICC) 
    if (WIN32)
         set(FMA_FLAGS "${FMA_FLAGS} /arch:core-avx2")
    else()
         set(FMA_FLAGS "${FMA_FLAGS} -march=core-avx2")
    endif()
endif()

SET(AVX512_FLAGS "-D__SSE3 -D__AVX")
if (VCC)
    message("AVX512 not available in Visual C++")
    #set(AVX512_FLAGS "${AVX512_FLAGS} /arch:AVX512")
elseif (CLANG)
    set(AVX512_FLAGS "${AVX512_FLAGS} -mavx512f -mfma")
elseif (GCC)
    set(AVX512_FLAGS "${AVX512_FLAGS} -mavx512f -mfma")
elseif (ICC) 
    if (WIN32)
         set(AVX512_FLAGS "${AVX512_FLAGS} /QxMIC-AVX512")
    else()
         set(AVX512_FLAGS "${AVX512_FLAGS} -xMIC-AVX512")
    endif()
endif()


# further flag to improve performance

if (IQTREE_FLAGS MATCHES "fma") # AVX+FMA instruction set
    message("Vectorization : AVX+FMA")
    add_definitions(-D__SSE3 -D__AVX) # define both SSE3 and AVX directive
    set(COMBINED_FLAGS "${COMBINED_FLAGS} ${FMA_FLAGS}")
    #SET(EXE_SUFFIX "${EXE_SUFFIX}-fma")
elseif (IQTREE_FLAGS MATCHES "avx") # AVX instruction set
    message("Vectorization : AVX")
    add_definitions(-D__SSE3 -D__AVX) # define both SSE3 and AVX directive
    set(COMBINED_FLAGS "${COMBINED_FLAGS} ${AVX_FLAGS}")
    #SET(EXE_SUFFIX "${EXE_SUFFIX}-avx")
elseif (NOT IQTREE_FLAGS MATCHES "nosse") #SSE intruction set
    if (IQTREE_FLAGS MATCHES "KNL")
        message("Vectorization : SSE3/AVX/AVX2/AVX-512")
        add_definitions(-D__AVX512KNL)
    else()
        message("Vectorization : SSE3/AVX/AVX2")
    endif()
    #add_definitions(-D__SSE3)
    #set(COMBINED_FLAGS "${COMBINED_FLAGS} ${SSE_FLAGS}")
endif()


##################################################################
# Setup compiler flags
##################################################################

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -fno-default-inline -fno-inline -O2 -fno-omit-frame-pointer -g")
set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -O2 -fno-omit-frame-pointer -g")

if(CLANG AND IQTREE_FLAGS MATCHES "static")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -Wl,--allow-multiple-definition")
endif()

if (IQTREE_FLAGS MATCHES "libcxx")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()


if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message("C flags       : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    message("CXX flags     : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("C flags       : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    message("CXX flags     : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Profile")
    message("C flags       : ${CMAKE_C_FLAGS_PROFILE} ")
    message("CXX flags     : ${CMAKE_CXX_FLAGS_PROFILE} ")
endif()

message("LINKER flags  : ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")


if (GCC)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
    set(CMAKE_CXX_FLAGS_MEM "-g -O1")
    set(CMAKE_C_FLAGS_MEM "-g -O1")
elseif (CLANG)
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
    set(CMAKE_CXX_FLAGS_MEM "-g -O1")
    set(CMAKE_C_FLAGS_MEM "-g -O1")
endif()

##################################################################
# check existence of a few basic functions
##################################################################
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists (getrusage HAVE_GETRUSAGE)
check_function_exists (GlobalMemoryStatusEx HAVE_GLOBALMEMORYSTATUSEX)
check_function_exists (strndup HAVE_STRNDUP)
check_function_exists (strtok_r HAVE_STRTOK_R)

find_package(Backtrace)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/iqtree_config.h.in"
  "${PROJECT_BINARY_DIR}/iqtree_config.h"
  )

# add the binary tree to the search path for include files
# so that we will find iqtree_config.h
include_directories("${PROJECT_BINARY_DIR}")

#zlib will be detected for appearance
#include_directories("${PROJECT_BINARY_DIR}/zlib-1.2.7")


if (NOT IQTREE_FLAGS MATCHES "nozlib")
    find_package(ZLIB)
endif()

if(ZLIB_FOUND)
  message ("Using system zlib")
  include_directories(${ZLIB_INCLUDE_DIRS})
else(ZLIB_FOUND)
  message ("Using own zlib-1.2.7")
  include_directories("${PROJECT_BINARY_DIR}/zlib-1.2.7" "${PROJECT_SOURCE_DIR}/zlib-1.2.7")
  add_subdirectory(zlib-1.2.7)
endif(ZLIB_FOUND)

##################################################################
# subdirectories containing necessary libraries for the build
##################################################################
add_subdirectory(pll)
add_subdirectory(ncl)
add_subdirectory(nclextra)
add_subdirectory(utils)
add_subdirectory(pda)
add_subdirectory(lbfgsb)
add_subdirectory(whtest)
add_subdirectory(sprng)
#add_subdirectory(zlib-1.2.7)
add_subdirectory(vectorclass)
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
add_subdirectory(model)
add_subdirectory(gsl)
add_subdirectory(alignment)
add_subdirectory(tree)

if (USE_TERRAPHAST)
   add_subdirectory(terraphast)
   add_subdirectory(terrace)
endif()

##################################################################
# the main executable
##################################################################

add_library(kernelsse tree/phylokernelsse.cpp)

if (NOT BINARY32 AND NOT IQTREE_FLAGS MATCHES "novx")
add_library(kernelavx tree/phylotreeavx.cpp)
add_library(kernelfma tree/phylokernelfma.cpp)
    if (IQTREE_FLAGS MATCHES "KNL")
        add_library(kernelavx512 tree/phylokernelavx512.cpp)
    endif()
endif()

if (IQTREE_FLAGS MATCHES "mpi")
    add_library(mympi utils/TreeCollection.cpp utils/ObjectStream.cpp)
endif()

if (NOT IQTREE_FLAGS MATCHES "single")
    if (BINARY32)
        link_directories(${PROJECT_SOURCE_DIR}/lib32)
    else()
        link_directories(${PROJECT_SOURCE_DIR}/lib)
    endif()
endif()

add_executable(iqtree
main/main.cpp
main/phyloanalysis.cpp
main/phyloanalysis.h
main/phylotesting.cpp
main/phylotesting.h
)

if(Backtrace_FOUND)
  include_directories(${Backtrace_INCLUDE_DIR})
  target_link_libraries(iqtree ${Backtrace_LIBRARY})
endif(Backtrace_FOUND)


if (NOT IQTREE_FLAGS MATCHES "avx" AND NOT IQTREE_FLAGS MATCHES "fma")
    if (NOT IQTREE_FLAGS MATCHES "nosse")
        set_target_properties(iqtree ncl nclextra utils pda lbfgsb whtest sprng vectorclass model alignment tree PROPERTIES COMPILE_FLAGS "${SSE_FLAGS}")
    endif()
    set_target_properties(kernelsse pll PROPERTIES COMPILE_FLAGS "${SSE_FLAGS}")
    if (NOT BINARY32 AND NOT IQTREE_FLAGS MATCHES "novx")
        set_target_properties(kernelavx pllavx PROPERTIES COMPILE_FLAGS "${AVX_FLAGS}")
        set_target_properties(kernelfma PROPERTIES COMPILE_FLAGS "${FMA_FLAGS}")
        if (IQTREE_FLAGS MATCHES "KNL")
            set_target_properties(kernelavx512 PROPERTIES COMPILE_FLAGS "${AVX512_FLAGS}")
        endif()
    endif()
endif()

##################################################################
# setup linking flags
##################################################################

# link special lib for WIN32
if (WIN32)
    set(PLATFORM_LIB "ws2_32")
else()
    set(PLATFORM_LIB "m")
endif()

if (IQTREE_FLAGS MATCHES "libcxx")
    set(STD_LIB "c++abi")
endif()

set(THREAD_LIB "")
if (NOT IQTREE_FLAGS MATCHES "single")
    link_directories(${PROJECT_SOURCE_DIR}/lib)
    if (MSVC)
        if (BINARY32)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LIBPATH:${PROJECT_SOURCE_DIR}/lib32")
            set(THREAD_LIB "pthreadVC2")
        else()
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LIBPATH:${PROJECT_SOURCE_DIR}/lib")
            set(THREAD_LIB "pthreadVC2")
        endif()
    elseif(CLANG AND APPLE)
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${PROJECT_SOURCE_DIR}/libmac")
    elseif(CLANG AND WIN32)
        if (BINARY32)
            target_link_libraries(iqtree ${PROJECT_SOURCE_DIR}/lib32/libiomp5md.dll)
        else()
            target_link_libraries(iqtree ${PROJECT_SOURCE_DIR}/lib/libiomp5md.dll)
        endif()
#        set(THREAD_LIB "ompstatic")
    endif()

    if (CLANG AND BINARY32)
        set (ATOMIC_LIB "atomic")
    endif()

endif()

# basic linking librararies
target_link_libraries(iqtree pll ncl nclextra utils pda lbfgsb whtest sprng vectorclass model gsl alignment tree ${PLATFORM_LIB} ${STD_LIB} ${THREAD_LIB} ${ATOMIC_LIB})

if (USE_TERRAPHAST)
   target_link_libraries(iqtree terrace)
endif()

if (NOT IQTREE_FLAGS MATCHES "nosse")
    target_link_libraries(iqtree kernelsse)
endif()

# MPI libraries
if (IQTREE_FLAGS MATCHES "mpi")
    target_link_libraries(iqtree mympi)
    if (NOT CMAKE_CXX_COMPILER MATCHES "mpi")
        target_link_libraries(iqtree ${MPI_CXX_LIBRARIES})
    endif()
endif()

# SSE, AVX etc. libraries
if (NOT BINARY32 AND NOT IQTREE_FLAGS MATCHES "novx")
    target_link_libraries(iqtree pllavx kernelavx kernelfma)
    if (IQTREE_FLAGS MATCHES "KNL")
        target_link_libraries(iqtree kernelavx512)
    endif()
endif()

# setup the executable name
##################################################################
set_target_properties(iqtree PROPERTIES OUTPUT_NAME "iqtree${EXE_SUFFIX}")

# strip the release build
if (NOT IQTREE_FLAGS MATCHES "nostrip" AND CMAKE_BUILD_TYPE STREQUAL "Release" AND (GCC OR CLANG) AND NOT APPLE) # strip is not necessary for MSVC
    if (WIN32)
        ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND strip $<TARGET_FILE:iqtree>)
    elseif (NOT APPLE)
        ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND strip $<TARGET_FILE:iqtree>)
    endif()
endif()

if (MSVC)
    set (BINARY_DIR "${PROJECT_BINARY_DIR}/Release")
else()
    set (BINARY_DIR "${PROJECT_BINARY_DIR}")
endif()

if (WIN32)
    if (MSVC)
        ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND copy "Release\\iqtree${EXE_SUFFIX}.exe" "Release\\iqtree${EXE_SUFFIX}-click.exe")
    else()
        ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND copy "iqtree${EXE_SUFFIX}.exe" "iqtree${EXE_SUFFIX}-click.exe")
    endif()
endif()

##############################################################
# add the install targets
##############################################################
install (TARGETS iqtree DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/example/models.nex" DESTINATION .)
install (FILES "${PROJECT_SOURCE_DIR}/example/example.phy" DESTINATION .)
install (FILES "${PROJECT_SOURCE_DIR}/example/example.nex" DESTINATION .)
install (FILES "${PROJECT_SOURCE_DIR}/example/example.cf" DESTINATION .)

if (WIN32)
    install (FILES "${BINARY_DIR}/iqtree${EXE_SUFFIX}-click.exe" DESTINATION bin)
    if (NOT IQTREE_FLAGS MATCHES "single" AND MSVC)
        if (BINARY32)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib32/pthreadVC2.dll" DESTINATION bin)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib32/libiomp5md.dll" DESTINATION bin)
        else()
            install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadVC2.dll" DESTINATION bin)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib/libiomp5md.dll" DESTINATION bin)
        endif()
#       install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2.dll" DESTINATION bin)
#       install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2_64.dll" DESTINATION bin)
    endif()

    if (NOT IQTREE_FLAGS MATCHES "single" AND CLANG)
        if (BINARY32)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib32/libiomp5md.dll" DESTINATION bin)
        else()
            install(FILES  "${PROJECT_SOURCE_DIR}/lib/libiomp5md.dll" DESTINATION bin)
        endif()
    endif()

endif()

##############################################################
# build a CPack driven installer package
##############################################################
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${iqtree_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${iqtree_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${iqtree_VERSION_PATCH}")
if(WIN32 OR APPLE)
  set(CPACK_GENERATOR "ZIP")
  set(CPACK_SOURCE_GENERATOR "ZIP")
else()
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR "TGZ")
endif()

#set(CPACK_SOURCE_PACKAGE_FILE_NAME
#  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
set(CPACK_SOURCE_IGNORE_FILES
  "/build.*/;/debug.*/;/examples/;/test_scripts/;/manual/;/.bzr/;~$;/\\\\.svn/;/\\\\.git/;/pllrepo/;${CPACK_SOURCE_IGNORE_FILES}")

set (SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    if (IQTREE_FLAGS MATCHES "oldmac")
        set (SYSTEM_NAME "MacOS10.5")
    else()
        set (SYSTEM_NAME "MacOSX")
    endif()
endif()

if (BINARY32) 
    set (SYSTEM_NAME "${SYSTEM_NAME}32")
endif()

if (IQTREE_FLAGS MATCHES "KNL")
    set (SYSTEM_NAME "${SYSTEM_NAME}KNL")
endif()

set(CPACK_PACKAGE_FILE_NAME 
    "${CMAKE_PROJECT_NAME}${EXE_SUFFIX}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-${SYSTEM_NAME}")

if (NOT APPLE)
    set(CPACK_STRIP_FILES TRUE)
endif()

include (CPack)

#add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
