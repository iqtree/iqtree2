add_library(utils
eigendecomposition.cpp eigendecomposition.h
gzstream.cpp gzstream.h
optimization.cpp optimization.h
stoprule.cpp stoprule.h
tools.cpp tools.h
pllnni.cpp pllnni.h
checkpoint.cpp checkpoint.h
MPIHelper.cpp MPIHelper.h
starttree.cpp starttree.h
flatmatrix.cpp flatmatrix.h
bionj.cpp
bionj2.cpp bionj2.h
stitchup.cpp
upgma.h nj.h rapidnj.h auctionmatrix.h
clustertree.h distancematrix.h heapsort.h
progress.cpp progress.h
timekeeper.cpp timekeeper.h
timeutil.h hammingdistance.h io.h
operatingsystem.cpp operatingsystem.h io.h my_assert.h)

if(ZLIB_FOUND)
  target_link_libraries(utils ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(utils zlibstatic)
endif(ZLIB_FOUND)

target_link_libraries(utils lbfgsb sprng)

add_executable(decentTree
    decenttree.cpp
    flatmatrix.cpp flatmatrix.h
    starttree.cpp starttree.h
    bionj.cpp
    bionj2.cpp bionj2.h
    upgma.h nj.h rapidnj.h auctionmatrix.h
    clustertree.h distancematrix.h hammingdistance.h
    stitchup.cpp
    
    heapsort.h
    gzstream.cpp gzstream.h
    progress.cpp progress.h
    my_assert.h io.h timeutil.h operatingsystem.cpp operatingsystem.h)

if(ZLIB_FOUND)
  target_link_libraries(decentTree ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(decentTree zlibstatic)
endif(ZLIB_FOUND)

