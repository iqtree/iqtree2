add_library(pda
circularnetwork.cpp circularnetwork.h
greedy.cpp greedy.h
gurobiwrapper.cpp gurobiwrapper.h
hashsplitset.cpp hashsplitset.h
pdnetwork.cpp pdnetwork.h
pruning.cpp pruning.h
split.cpp split.h
splitgraph.cpp splitgraph.h
splitset.cpp splitset.h
lpwrapper.c lpwrapper.h
ecopd.cpp ecopd.h
graph.cpp graph.h
pdtree.cpp pdtree.h
pdtreeset.cpp pdtreeset.h
ecopdmtreeset.cpp ecopdmtreeset.h
)

target_link_libraries(pda nclextra)
