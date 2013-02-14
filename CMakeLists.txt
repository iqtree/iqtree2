cmake_minimum_required(VERSION 2.6)

set(BUILD_SHARED_LIBS OFF)


######## NOTE: run cmake with option "-DIQ_STATIC to build static library 
if (WIN32)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
elseif (NOT APPLE) 
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

IF($ENV{STATIC})
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
	IF(WIN32)
    		SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
	ELSE(WIN32)
    		SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
	ENDIF(WIN32)
ENDIF($ENV{STATIC})

if ($ENV{B32}) 
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
endif ($ENV{B32}) 

if ($ENV{OMP}) 
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
endif ($ENV{OMP}) 

project(iqtree)
# The version number.
set (iqtree_VERSION_MAJOR 0)
set (iqtree_VERSION_MINOR 9)
set (iqtree_VERSION_PATCH 2)


include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists (getrusage HAVE_GETRUSAGE)
#check_function_exists (popen HAVE_POPEN)
#check_function_exists (pclose HAVE_PCLOSE)
check_function_exists (GlobalMemoryStatusEx HAVE_GLOBALMEMORYSTATUSEX)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/iqtree_config.h.in"
  "${PROJECT_BINARY_DIR}/iqtree_config.h"
  )
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

add_definitions(-DIQ_TREE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3 -std=c++98 -Wall -Wno-sign-compare -Wno-deprecated -pedantic -g -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3 -std=c99 -Wall -Wno-sign-compare -pedantic -g -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS")

#if (NOT WIN32)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
#endif(NOT WIN32)
 
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_C_FLAGS_DEBUG "-O0")

if (NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  message("Warning: 32-bit version is not completely working")
else(CMAKE_SIZEOF_VOID_P EQUAL 4)
  message("Info: Compiling on 64-bit system")
endif(CMAKE_SIZEOF_VOID_P EQUAL 4)

add_subdirectory(phylolib)
add_subdirectory(ncl)
add_subdirectory(whtest)
add_subdirectory(sprng)

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
nnisearch.c
)

if (WIN32)
target_link_libraries(iqtree ncl whtest z sprng ws2_32)
else (WIN32)
target_link_libraries(iqtree phylolib ncl whtest z sprng)
endif (WIN32)

##############################################################
# add the install targets
##############################################################
install (TARGETS iqtree DESTINATION bin)
install (FILES "${PROJECT_SOURCE_DIR}/examples/example.phy" DESTINATION examples)
#install (FILES "${PROJECT_BINARY_DIR}/iqtree_config.h"        
#         DESTINATION include)

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
  "/build.*/;/debug.*/;/examples/;/manual/;/.bzr/;~$;/\\\\.svn/;${CPACK_SOURCE_IGNORE_FILES}")

set(CPACK_STRIP_FILES TRUE)

include (CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
