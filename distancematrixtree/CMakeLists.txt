add_library(distancematrixtree
starttree.cpp starttree.h
flatmatrix.cpp flatmatrix.h
bionj.cpp
bionj2.cpp bionj2.h
stitchup.cpp
upgma.h nj.h rapidnj.h auctionmatrix.h
clustertree.h distancematrix.h)

add_executable(decentTree
    decenttree.cpp
    flatmatrix.cpp flatmatrix.h
    starttree.cpp starttree.h
    bionj.cpp
    bionj2.cpp bionj2.h
    upgma.h nj.h rapidnj.h auctionmatrix.h
    clustertree.h distancematrix.h 
    stitchup.cpp
    ../utils/hammingdistance.h
    ../utils/vectortypes.cpp ../utils/vectortypes.h
    ../utils/heapsort.h
    ../utils/gzstream.cpp ../utils/gzstream.h
    ../utils/progress.cpp ../utils/progress.h
    ../utils/my_assert.h ../utils/safe_io.h ../utils/timeutil.h 
    ../utils/operatingsystem.cpp ../utils/operatingsystem.h)

if(ZLIB_FOUND)
  target_link_libraries(decentTree ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(decentTree zlibstatic)
endif(ZLIB_FOUND)

