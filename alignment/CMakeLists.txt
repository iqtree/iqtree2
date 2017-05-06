add_library(alignment
pattern.cpp
pattern.h
alignment.cpp
alignment.h
alignmentpairwise.cpp
alignmentpairwise.h
maalignment.cpp
maalignment.h
superalignment.cpp
superalignment.h
superalignmentpairwise.cpp
superalignmentpairwise.h
superalignmentpairwiseplen.cpp
superalignmentpairwiseplen.h
)

target_link_libraries(alignment ncl gsl)
