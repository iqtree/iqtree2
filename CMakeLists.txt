##################################################################
# IQ-TREE cmake build definition
# Copyright (c) 2012-2015 Bui Quang Minh, Lam Tung Nguyen
##################################################################

# Windows example usages:
#------------------------
# cmake -G "Visual Studio 12" <source_dir>			(32-bit version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" <source_dir>	(64-bit version, compiled with MSVC)
# cmake -G "Visual Studio 12 Win64" -T "Intel C++ Compiler XE 15.0" <source_dir>		(64-bit version, compiled with ICC)
# cmake -G "MinGW Makefiles" <source_dir> (TDM-GCC)
# cmake -G "Unix Makefiles" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_MAKE_PROGRAM=mingw32-make <source_dir> (TDM-GCC and clang)

# Linux example usages:
#----------------------
# cmake <source_dir>						   (sequential version)
# cmake -DIQTREE_FLAGS="omp" <source_dir>      (OpenMP version)
# cmake -DIQTREE_FLAGS="m32" <source_dir>      (32-bit sequential version)
# cmake -DIQTREE_FLAGS="m32 omp" <source_dir>  (32-bit OpenMP version)
#
# To compile with CLANG on Linux:
# export CC=/usr/bin/clang
# export CXX=/usr/bin/clang++
# Best practices for setting up CMAKE for diffrent compiler can be found here:
# http://stackoverflow.com/questions/7031126/switching-between-gcc-and-clang-llvm-using-cmake
#
# Mac OSX example usages:
#------------------------
# cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ <source_dir>							(sequential version)
#
# To build OpenMP version one needs to download Clang version 3.7 or later (as of November 2015) 
# Then assuming clang3.7 and clang++3.7 are the newly built compilers, then:
# cmake -DCMAKE_C_COMPILER=clang3.7 -DCMAKE_CXX_COMPILER=clang++3.7 -DIQTREE_FLAGS="omp" <source_dir>		(OpenMP version)
#


# Compile OpenMP version: cmake -DIQTREE_FLAGS=omp ....
# Compile 32-bit version: cmake -DIQTREE_FLAGS=m32 ....
# Compile static version: cmake -DIQTREE_FLAGS=static ....
# Compile static OpenMP version: cmake -DIQTREE_FLAGS="omp static" ....

#NOTE: Static linking with clang windows: make a symlink libgcc_eh.a to libgcc.a (administrator required)
# C:\TDM-GCC-64\lib\gcc\x86_64-w64-mingw32\5.1.0>mklink libgcc_eh.a libgcc.a

cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

##################################################################
# configure MPI compilation
##################################################################
if (IQTREE_FLAGS MATCHES "mpi")
    message("MPI enabled")
	set(IQTREE_FLAGS "${IQTREE_FLAGS} dynamic")
	add_definitions(-D_IQTREE_MPI)
    set(CMAKE_CXX_COMPILER "mpicxx")
	set(CMAKE_C_COMPILER "mpicc")
    if (IQTREE_FLAGS MATCHES "debug")
        add_definitions(-D_MPI_DEBUG)
        message("Debug information enabled")
    endif()
endif()

project(iqtree)
add_definitions(-DIQ_TREE)
# The version number.
set (iqtree_VERSION_MAJOR 1)
set (iqtree_VERSION_MINOR 4)
set (iqtree_VERSION_PATCH "1-rc") 

set(BUILD_SHARED_LIBS OFF)

message("IQ-TREE flags : ${IQTREE_FLAGS}")

if (NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release") 
	message("Builde mode   : Release")
endif()

include_directories("${PROJECT_SOURCE_DIR}")


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
	# to be compatible back to Mac OS X 10.6
	if (IQTREE_FLAGS MATCHES "oldmac") 
		add_definitions("-mmacosx-version-min=10.5") 
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.5")
	else()
		add_definitions("-mmacosx-version-min=10.6") 
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.6")
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

if (CMAKE_COMPILER_IS_GNUCXX) 	
	message("Compiler      : GNU Compiler (gcc)")
	set(GCC "TRUE")
#	set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -fms-extensions -Wno-deprecated")
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -g")
    # require at least gcc 4.6
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.6)
        message(FATAL_ERROR "GCC version must be at least 4.6!")
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	message("Compiler      : Clang")
	set(CLANG "TRUE")
#	set(COMBINED_FLAGS "-Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE -Wno-nested-anon-types")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
	set(CMAKE_C_FLAGS_RELEASE "-O3")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
	set(VCC "TRUE")
	message("Compiler      : MS Visual C++ Compiler")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
	message("Compiler      : Intel C++ Compiler (icc)")
	set(ICC "TRUE")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
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
			add_definitions(/O3)
		endif()
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
	SET(EXE_SUFFIX "${EXE_SUFFIX}32")
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
if (IQTREE_FLAGS MATCHES "omp")
	message("OpenMP      : Yes")
	SET(EXE_SUFFIX "${EXE_SUFFIX}-omp")	
	add_definitions(-D_USE_PTHREADS)
	if (MSVC) 
		add_definitions(/MT)
	endif()
	
	if (VCC) 
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
  		include_directories("${PROJECT_SOURCE_DIR}/pll") # for PThreads headers 
	elseif (ICC)
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qopenmp")
  		if (WIN32)
  			include_directories("${PROJECT_SOURCE_DIR}/pll") # for PThreads headers
  		endif() 
  	elseif (GCC)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  	elseif (CLANG) 
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
  	endif()
else()
	message("OpenMP      : NONE")
endif()


if (IQTREE_FLAGS MATCHES "mpi")
	message("MPI      : Yes")
	SET(EXE_SUFFIX "${EXE_SUFFIX}-mpi")	
else()
	message("MPI      : NONE")
endif()

##################################################################
# configure SSE/AVX/FMA instructions
##################################################################

SET(AVX_FLAGS "-D__AVX")
if (VCC) 
	set(AVX_FLAGS "${AVX_FLAGS} /arch:AVX")
elseif (CLANG)
	set(AVX_FLAGS "${AVX_FLAGS} -mavx")
elseif (GCC)
	set(AVX_FLAGS "${AVX_FLAGS} -mavx -fabi-version=0")
elseif (ICC) 
	if (WIN32)
		 set(AVX_FLAGS "${AVX_FLAGS} /arch:AVX")
	else()
		 set(AVX_FLAGS "${AVX_FLAGS} -mavx")
	endif()
endif()

SET(SSE_FLAGS "")
if (VCC)
	set(SSE_FLAGS "/arch:SSE2 -D__SSE3__")
elseif (GCC OR CLANG)
	set(SSE_FLAGS "-msse3")
elseif (ICC)
	if (WIN32)
		set(SSE_FLAGS "/arch:SSE3")
	else()
		set(SSE_FLAGS "-msse3")
	endif()
endif()

if (IQTREE_FLAGS MATCHES "fma") # AVX+FMA instruction set
 	message("Vectorization : AVX+FMA")
	add_definitions(-D__SSE3 -D__AVX) # define both SSE3 and AVX directive
	if (VCC)
		# Visual C++ has no /mfma flag!, FMA is only included in AVX2 
		set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX2")
	elseif (CLANG)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx -mfma")
	elseif (GCC)
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx -fabi-version=0 -mfma")
	elseif (ICC) 
		if (WIN32)
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} /arch:AVX /Qfma")
		else()
			 set(COMBINED_FLAGS "${COMBINED_FLAGS} -mavx -mfma")
		endif()
	endif()

	SET(EXE_SUFFIX "${EXE_SUFFIX}-fma")

elseif (IQTREE_FLAGS MATCHES "avx") # AVX instruction set
 	message("Vectorization : AVX")
	add_definitions(-D__SSE3 -D__AVX) # define both SSE3 and AVX directive
	set(COMBINED_FLAGS "${COMBINED_FLAGS} ${AVX_FLAGS}")

	SET(EXE_SUFFIX "${EXE_SUFFIX}-avx")
else() #SSE intruction set
	message("Vectorization : SSE3")
	add_definitions(-D__SSE3)

endif()


##################################################################
# Setup compiler flags
##################################################################

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMBINED_FLAGS}")
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -fno-default-inline -fno-inline -O2 -fno-omit-frame-pointer -g")
set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS} -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -O2 -fno-omit-frame-pointer -g")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
	message("C flags    : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}") 
	message("CXX flags  : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}") 
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	message("C flags    : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}") 
	message("CXX flags  : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}") 
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Profile")
	message("C flags    : ${CMAKE_C_FLAGS_PROFILE} ") 
	message("CXX flags  : ${CMAKE_CXX_FLAGS_PROFILE} ") 
endif()


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


##################################################################
# subdirectories containing necessary libraries for the build 
##################################################################
add_subdirectory(pll)
add_subdirectory(ncl)
add_subdirectory(lbfgsb)
add_subdirectory(whtest)
add_subdirectory(sprng)
#add_subdirectory(zlib-1.2.7)
add_subdirectory(vectorclass)
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
find_package(Eigen3)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)
add_subdirectory(model)
add_subdirectory(gsl)

##################################################################
# the main executable
##################################################################

if (NOT BINARY32 AND NOT IQTREE_FLAGS MATCHES "novx")
    add_library(avxkernel phylotreeavx.cpp)
endif()

if (IQTREE_FLAGS MATCHES "mpi")
	add_library(mpi TreeCollection.cpp ObjectStream.cpp)
endif()

add_executable(iqtree
alignment.cpp
alignmentpairwise.cpp
circularnetwork.cpp
eigendecomposition.cpp
greedy.cpp
gss.cpp
#guidedbootstrap.cpp
gurobiwrapper.cpp
gzstream.cpp
hashsplitset.cpp
iqtree.cpp
maalignment.cpp
matree.cpp
mexttree.cpp
mpdablock.cpp
msetsblock.cpp
msplitsblock.cpp
modelsblock.cpp
mtree.cpp
mtreeset.cpp
ncbitree.cpp
ngs.cpp
node.cpp
optimization.cpp
parsmultistate.cpp
pattern.cpp
pda.cpp
pdnetwork.cpp
pdtree.cpp
pdtreeset.cpp
phyloanalysis.cpp
phylonode.cpp
phylosupertree.cpp
phylotree.cpp
phylotreesse.cpp
phylotreepars.cpp
phylokernelsitemodel.cpp
#phylotreeavx.cpp
pruning.cpp
quartet.cpp
split.cpp
splitgraph.cpp
splitset.cpp
stoprule.cpp
superalignment.cpp
superalignmentpairwise.cpp
supernode.cpp
tinatree.cpp
tools.cpp
whtest_wrapper.cpp
lpwrapper.c
pllnni.cpp
phylosupertreeplen.cpp
phylotesting.cpp
ecopd.cpp
ecopdmtreeset.cpp
graph.cpp
candidateset.cpp
checkpoint.cpp
upperbounds.cpp
MPIHelper.cpp
)

if (NOT IQTREE_FLAGS MATCHES "nozlib")
    find_package(ZLIB)
endif()

if(ZLIB_FOUND)
  message ("Using system zlib")
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(iqtree ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  message ("Using own zlib-1.2.7")
  include_directories("${PROJECT_BINARY_DIR}/zlib-1.2.7" "${PROJECT_SOURCE_DIR}/zlib-1.2.7")
  add_subdirectory(zlib-1.2.7)
  target_link_libraries(iqtree zlibstatic)
endif(ZLIB_FOUND)

if (NOT IQTREE_FLAGS MATCHES "avx" AND NOT IQTREE_FLAGS MATCHES "fma")
	set_target_properties(iqtree pll ncl lbfgsb whtest sprng vectorclass model PROPERTIES COMPILE_FLAGS "${SSE_FLAGS}")
	if (NOT BINARY32 AND NOT IQTREE_FLAGS MATCHES "novx")
		set_target_properties(avxkernel pllavx PROPERTIES COMPILE_FLAGS "${AVX_FLAGS}")
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

if(CLANG AND WIN32 AND IQTREE_FLAGS MATCHES "static")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -Wl,--allow-multiple-definition")
endif()

set(THREAD_LIB "")
if (IQTREE_FLAGS MATCHES "omp") 
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
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${PROJECT_SOURCE_DIR}/lib32 libiomp5md.dll")
        else()
            set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${PROJECT_SOURCE_DIR}/lib libiomp5md.dll")
        endif()
#        set(THREAD_LIB "ompstatic")
	endif()
endif()

set (MYMPI_LIB "")
if (IQTREE_FLAGS MATCHES "mpi")
    set (MYMPI_LIB "mpi")
endif()

if (BINARY32 OR IQTREE_FLAGS MATCHES "novx")
    target_link_libraries(iqtree pll ncl lbfgsb whtest sprng vectorclass model ${PLATFORM_LIB} ${STD_LIB} ${THREAD_LIB} ${MYMPI_LIB})	
else()
    target_link_libraries(iqtree pll pllavx ncl lbfgsb whtest sprng vectorclass model avxkernel ${PLATFORM_LIB} ${STD_LIB} ${THREAD_LIB} ${MYMPI_LIB})	
endif()

# setup the executable name 
##################################################################
set_target_properties(iqtree PROPERTIES OUTPUT_NAME "iqtree${EXE_SUFFIX}")

# strip the release build
if (CMAKE_BUILD_TYPE STREQUAL "Release" AND (GCC OR CLANG)) # strip is not necessary for MSVC
	if (WIN32)
		ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND strip $<TARGET_FILE:iqtree>)
	else()
		ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND ${CMAKE_STRIP} $<TARGET_FILE:iqtree>)
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

if (WIN32)
	install (FILES "${BINARY_DIR}/iqtree${EXE_SUFFIX}-click.exe" DESTINATION bin)
	if (EXE_SUFFIX MATCHES "omp" AND MSVC)
        if (BINARY32)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib32/pthreadVC2.dll" DESTINATION bin)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib32/libiomp5md.dll" DESTINATION bin)
        else()
            install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadVC2.dll" DESTINATION bin)
            install(FILES  "${PROJECT_SOURCE_DIR}/lib/libiomp5md.dll" DESTINATION bin)
        endif()
#		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2.dll" DESTINATION bin)
#		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2_64.dll" DESTINATION bin)
	endif()

	if (EXE_SUFFIX MATCHES "omp" AND CLANG)
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
     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
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

set (PROJECT_NAME_SUFFIX "${EXE_SUFFIX}")
#if (NOT IQTREE_FLAGS MATCHES "omp" AND NOT IQTREE_FLAGS MATCHES "avx" AND NOT IQTREE_FLAGS MATCHES "fma")  
#	set (PROJECT_NAME_SUFFIX "${PROJECT_NAME_SUFFIX}-sse") 
#endif()

set(CPACK_PACKAGE_FILE_NAME 
	"${CMAKE_PROJECT_NAME}${PROJECT_NAME_SUFFIX}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-${SYSTEM_NAME}")

set(CPACK_STRIP_FILES TRUE)

include (CPack)

#add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
