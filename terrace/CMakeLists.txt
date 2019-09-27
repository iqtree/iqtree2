add_library(terrace
terrace.cpp
terrace.h
)

target_link_libraries(terrace tree alignment terraphast)