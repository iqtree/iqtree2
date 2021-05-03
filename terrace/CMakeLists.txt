add_library(terrace
terracenode.cpp
terracenode.hpp
terracetree.cpp
terracetree.hpp
presenceabsencematrix.cpp
presenceabsencematrix.hpp
terrace.hpp
terrace.cpp
)

target_link_libraries(terrace tree)
