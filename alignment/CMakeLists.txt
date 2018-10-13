add_library(alignment
statespace.h
statespace.cpp
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
superalignmentunlinked.cpp
superalignmentunlinked.h
)

target_link_libraries(alignment ncl gsl)
