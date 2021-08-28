add_library(utils
eigendecomposition.cpp eigendecomposition.h
gzstream.cpp gzstream.h
optimization.cpp optimization.h
stoprule.cpp stoprule.h
tools.cpp tools.h
pllnni.cpp pllnni.h
checkpoint.cpp checkpoint.h
MPIHelper.cpp MPIHelper.h
heapsort.h
parallel_sort.h parallel_mergesort.h
progress.cpp progress.h
statefrequency.cpp statefrequency.h
stringfunctions.cpp stringfunctions.h
timekeeper.cpp timekeeper.h
timeutil.h hammingdistance.h safe_io.h
operatingsystem.cpp operatingsystem.h my_assert.h
vectortypes.cpp vectortypes.h)

if(ZLIB_FOUND)
  target_link_libraries(utils ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(utils zlibstatic)
endif(ZLIB_FOUND)

target_link_libraries(utils lbfgsb sprng)
