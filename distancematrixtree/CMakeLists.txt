add_library(distancematrixtree
starttree.cpp starttree.h 
flatmatrix.cpp flatmatrix.h hashrow.h
clustertree.h upgma.h nj.h rapidnj.h
bionj.cpp bionj2.cpp bionj2.h
)

if(ZLIB_FOUND)
  target_link_libraries(distancematrixtree ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(distancematrixtree zlibstatic)
endif(ZLIB_FOUND)

target_link_libraries(distancematrixtree)
