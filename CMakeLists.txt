cmake_minimum_required(VERSION 2.6)

project(iqtree)
# The version number.
set (iqtree_VERSION_MAJOR 0)
set (iqtree_VERSION_MINOR 9)
set (iqtree_VERSION_PATCH 0)


include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (gettimeofday HAVE_GETTIMEOFDAY)

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
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98 -Wall -Wno-deprecated -Wno-sign-compare -Wno-unused")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wall -Wno-sign-compare -Wno-unused")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98 -Wall -Wno-deprecated -Wno-sign-compare -Wno-unused -static-libstdc++")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wall -Wno-sign-compare -Wno-unused -static-libgcc")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0")
set(CMAKE_C_FLAGS_RELEASE "-O3 -g0")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g2")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g2")

if (NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE "Release")
endif()

add_subdirectory(ncl)
add_subdirectory(whtest)
add_subdirectory(modeltest)
add_subdirectory(lp)

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
iqptree.cpp
maalignment.cpp
matree.cpp
mexttree.cpp
modelbin.cpp
modeldna.cpp
modelfactory.cpp
modelnonrev.cpp
modelprotein.cpp
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
split.cpp
splitgraph.cpp
splitset.cpp
stoprule.cpp
substmodel.cpp
superalignment.cpp
superalignmentpairwise.cpp
supernode.cpp
tinatree.cpp
tools.cpp
whtest_wrapper.cpp
lpwrapper.c
modeltest_wrapper.c)

if (WIN32)
target_link_libraries(iqtree ncl whtest modeltest lp z)
else (WIN32)
target_link_libraries(iqtree ncl whtest modeltest lp z dl)
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
  "/build/;/debug/;/examples/;/manual/;/.bzr/;~$;/\\\\.svn/;${CPACK_SOURCE_IGNORE_FILES}")

include (CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
