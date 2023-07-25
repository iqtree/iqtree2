add_library(alignment
pattern.cpp
pattern.h
alignment.cpp
alignment.h
alignmentpairwise.cpp
alignmentpairwise.h
alignmentsummary.cpp
alignmentsummary.h
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
sequence.cpp
sequence.h
sequencechunkstr.cpp
sequencechunkstr.h
)

target_link_libraries(alignment simulator ncl gsl model)
