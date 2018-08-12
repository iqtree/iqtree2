add_library(terrace
terrace.cpp
terrace.h
)

target_link_libraries(terrace tree alignment terraces) # terraces is Terraphast