#########################################################
# IQ-TREE cmake build definition
#########################################################

cmake_minimum_required(VERSION 2.8)

project(iqtree)
# The version number.
set (iqtree_VERSION_MAJOR 1)
set (iqtree_VERSION_MINOR 0)
set (iqtree_VERSION_PATCH 2)

set(BUILD_SHARED_LIBS OFF)

if (IQTREE_FLAGS) 
	message("IQTREE_FLAGS = ${IQTREE_FLAGS}")
endif()

if (NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Detect platforms, tested for Windows, Mac OS X, and Unix-like
if (WIN32)
	message("Target OS     : Windows")
	# build as static binary to run on most machines
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    add_definitions(-DWIN32)
elseif (APPLE) 
	message("Target OS     : Mac OS X")
	# to be compatible back to Mac OS X 10.6
	add_definitions("-mmacosx-version-min=10.6") 
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=10.6")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
elseif (UNIX) 
	message("Target OS     : Unix")
	# build as static binary to run on most machines
  if (CMAKE_BUILD_TYPE STREQUAL "Release") 	
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  endif() 
else()
	# Note that IQ-TREE has NOT been tested on other platforms
	message("Target OS     : Unknown (untested yet)")
endif()

if (CMAKE_COMPILER_IS_GNUCXX) 
	message("Compiler      : GNU Compiler (gcc)")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	message("Compiler      : Clang")
else()
	message("Compiler      : Unknown (untested yet)")
endif()

# verbose for 32 or 64 bit platform
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
	message("Target binary : possibly 32-bit")
	#SET(EXE_SUFFIX "-32bit")
elseif(CMAKE_CXX_FLAGS MATCHES "m32")
	message("Target binary : 32-bit")
	#SET(EXE_SUFFIX "-32bit")
else()
	message("Target binary : 64-bit")
	#SET(EXE_SUFFIX "")
endif()

# change the executable name if compiled for OpenMP parallel version
if (IQTREE_FLAGS MATCHES "omp")
	message("Parallel      : OpenMP/PThreads")
	SET(EXE_SUFFIX "-omp")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_USE_PTHREADS -pthread")
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
else()
	message("Parallel      : None")
	SET(EXE_SUFFIX "")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release") 
	message("Builde mode   : Release + Strip")
endif()

# check existence of a few basic functions
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists (getrusage HAVE_GETRUSAGE)
check_function_exists (GlobalMemoryStatusEx HAVE_GLOBALMEMORYSTATUSEX)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/iqtree_config.h.in"
  "${PROJECT_BINARY_DIR}/iqtree_config.h"
  )

# add the binary tree to the search path for include files
# so that we will find iqtree_config.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/zlib-1.2.7")

add_definitions(-DIQ_TREE)

if (IQTREE_FLAGS MATCHES "avx")
	set(COMBINED_FLAGS "-msse3 -D__SSE3 -D__AVX -mavx")
	if (CMAKE_COMPILER_IS_GNUCXX) 
		set(COMBINED_FLAGS "${COMBINED_FLAGS} -fabi-version=0")
	endif()
	SET(EXE_SUFFIX "${EXE_SUFFIX}-avx")
else()
	set(COMBINED_FLAGS "-msse3 -D__SSE3")
endif()

if (MSVC) 
	set(COMBINED_FLAGS "${COMBINED_FLAGS} -Wall -pedantic -D_GNU_SOURCE")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fms-extensions -std=c++98 ${COMBINED_FLAGS}")
	set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -fms-extensions ${COMBINED_FLAGS}")
else ()
	set(COMBINED_FLAGS "${COMBINED_FLAGS} -Wall -Wno-unused-function -Wno-sign-compare -pedantic -D_GNU_SOURCE")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fms-extensions -std=c++98 -Wno-deprecated ${COMBINED_FLAGS}")
	set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -fms-extensions ${COMBINED_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g0")
set(CMAKE_C_FLAGS_RELEASE "-O2 -g0")

set(CMAKE_CXX_FLAGS_PROFILE "-fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -fno-default-inline -fno-inline -O0 -fno-omit-frame-pointer -pg")
set(CMAKE_C_FLAGS_PROFILE "-fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -O0 -fno-omit-frame-pointer -pg")

set(CMAKE_CXX_FLAGS_MEM "-g -O1")
set(CMAKE_C_FLAGS_MEM "-g -O1")

if (CMAKE_COMPILER_IS_GNUCXX) 
	set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
	set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline-functions-called-once -fno-default-inline -fno-inline")
else ()
	set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
	set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-inline-functions -fno-inline")
endif()

# subdirectories containing necessary libraries for the build 
add_subdirectory(pll)
add_subdirectory(ncl)
add_subdirectory(whtest)
add_subdirectory(sprng)
add_subdirectory(zlib-1.2.7)
#add_subdirectory(vectorclass)

# the main executable
add_executable(iqtree
alignment.cpp
alignmentpairwise.cpp
circularnetwork.cpp
eigendecomposition.cpp
greedy.cpp
gss.cpp
gtrmodel.cpp
guidedbootstrap.cpp
gurobiwrapper.cpp
gzstream.cpp
hashsplitset.cpp
iqtree.cpp
maalignment.cpp
matree.cpp
mexttree.cpp
modelbin.cpp
modeldna.cpp
modelfactory.cpp
modelnonrev.cpp
modelprotein.cpp
modelset.cpp
modelsubst.cpp
mpdablock.cpp
msetsblock.cpp
msplitsblock.cpp
mtree.cpp
mtreeset.cpp
ncbitree.cpp
ngs.cpp
node.cpp
optimization.cpp
parsmultistate.cpp
partitionmodel.cpp
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
#phylotreeeigen.cpp
pruning.cpp
rategamma.cpp
rategammainvar.cpp
rateheterogeneity.cpp
rateinvar.cpp
ratemeyerdiscrete.cpp
ratemeyerhaeseler.cpp
ratekategory.cpp
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
#modeltest_wrapper.c
nnisearch.cpp
modelcodon.cpp
phylosupertreeplen.cpp
phylotesting.cpp
ecopd.cpp
ecopdmtreeset.cpp
graph.cpp
modelmorphology.cpp
candidateset.cpp
upperbounds.cpp
)


# link special lib for WIN32
if (WIN32)
		target_link_libraries(iqtree pll ncl whtest zlibstatic sprng ws2_32)	
else (WIN32)
	target_link_libraries(iqtree pll ncl whtest zlibstatic sprng m)
endif (WIN32)

# setup the executable name properly
set_target_properties(iqtree PROPERTIES OUTPUT_NAME "iqtree${EXE_SUFFIX}")

# strip the release build
if (CMAKE_BUILD_TYPE STREQUAL "Release") 
	ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND ${CMAKE_STRIP} $<TARGET_FILE:iqtree>)
endif()

if (WIN32)
	ADD_CUSTOM_COMMAND(TARGET iqtree POST_BUILD COMMAND copy "iqtree${EXE_SUFFIX}.exe" "iqtree${EXE_SUFFIX}-click.exe")
endif()

##############################################################
# add the install targets
##############################################################
install (TARGETS iqtree DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/examples/example.phy" DESTINATION examples)
#install (FILES "${PROJECT_SOURCE_DIR}/examples/test.nex" DESTINATION examples)
#install (FILES "${PROJECT_SOURCE_DIR}/examples/test.budget" DESTINATION examples)
#install (FILES "${PROJECT_SOURCE_DIR}/examples/test.area.nex" DESTINATION examples)
#install (FILES "${PROJECT_SOURCE_DIR}/examples/test.area.budget" DESTINATION examples)
if (WIN32)
	install (FILES "${PROJECT_BINARY_DIR}/iqtree${EXE_SUFFIX}-click.exe" DESTINATION bin)
# obsolete: pthread is statically linked, so no dll is required anymore
#	if (EXE_SUFFIX MATCHES "omp")
#		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2.dll" DESTINATION bin)
#		install(FILES  "${PROJECT_SOURCE_DIR}/lib/pthreadGC2_64.dll" DESTINATION bin)
#	endif()
endif()

if (APPLE AND IQTREE_FLAGS MATCHES "omp")
	install (FILES "${PROJECT_SOURCE_DIR}/lib/libiomp5.dylib" DESTINATION bin)
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
if(WIN32)
  set(CPACK_GENERATOR "ZIP")
  set(CPACK_SOURCE_GENERATOR "ZIP")
else(WIN32)
  set(CPACK_GENERATOR "TGZ")
  set(CPACK_SOURCE_GENERATOR "TGZ")
endif(WIN32)
#set(CPACK_SOURCE_PACKAGE_FILE_NAME
#  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
set(CPACK_SOURCE_IGNORE_FILES
  "/build.*/;/debug.*/;/examples/;/manual/;/.bzr/;~$;/\\\\.svn/;/\\\\.git/;${CPACK_SOURCE_IGNORE_FILES}")

set(CPACK_PACKAGE_FILE_NAME 
	"${CMAKE_PROJECT_NAME}${EXE_SUFFIX}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-${CMAKE_SYSTEM_NAME}")

set(CPACK_STRIP_FILES TRUE)

include (CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
