add_compile_definitions(USE_PROGRESS_DISPLAY=1)
add_compile_definitions(USE_VECTORCLASS_LIBRARY=1)
add_compile_definitions(USE_GZSTREAM=1)

add_executable(decenttree
    decenttree.cpp
    flatmatrix.cpp flatmatrix.h
    sequence.cpp   sequence.h
    starttree.cpp  starttree.h
    bionj.cpp
    bionj2.cpp     bionj2.h
    upgma.h        nj.h  
    rapidnj.h      auctionmatrix.h
    clustertree.h  distancematrix.h 
    stitchup.cpp  
    
    ../utils/argument.cpp
    ../utils/hammingdistance.h    ../utils/heapsort.h 
    ../utils/timeutil.h           ../utils/safe_io.h
    ../utils/parallel_mergesort.h ../utils/parallel_sort.h
    ../utils/gzstream.cpp         ../utils/gzstream.h
    ../utils/progress.cpp         ../utils/progress.h 
    ../utils/operatingsystem.cpp  ../utils/operatingsystem.h
    ../utils/stringfunctions.cpp  ../utils/stringfunctions.h
    ../utils/vectortypes.cpp      ../utils/vectortypes.h)

add_library(distancematrixtree
starttree.cpp starttree.h 
flatmatrix.cpp flatmatrix.h hashrow.h
clustertree.h upgma.h nj.h rapidnj.h
bionj.cpp bionj2.cpp bionj2.h
)

if(ZLIB_FOUND)
  target_link_libraries(decenttree ${ZLIB_LIBRARIES})
  target_link_libraries(distancematrixtree ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(distancematrixtree zlibstatic)
  target_link_libraries(distancematrixtree zlibstatic)
endif(ZLIB_FOUND)
