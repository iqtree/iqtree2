cmake_minimum_required(VERSION 2.6)

project(iqtree)
# The version number.
set (iqtree_VERSION_MAJOR 0)
set (iqtree_VERSION_MINOR 9)


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

#add_definitions(-std=c++0x)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIQ_TREE -O3 -Wno-deprecated -g0 -std=c++0x -Wall -Wno-sign-compare -Wno-unused-but-set-variable")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIQ_TREE -O3 -Wno-deprecated -g0 -std=c++98 -Wall -Wno-sign-compare -Wno-unused")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DIQ_TREE -O3 -g0 -Wall -std=gnu99 -Wno-sign-compare -Wno-unused")

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