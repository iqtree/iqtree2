add_library(main
main.cpp
phyloanalysis.cpp
phyloanalysis.h
phylotesting.cpp
phylotesting.h
treetesting.cpp
treetesting.h
timetree.cpp
timetree.h
)

if (USE_BOOSTER)
    target_link_libraries(main booster)
endif()

if (USE_LSD2)
    target_link_libraries(main lsd2)
endif()

target_link_libraries(main pda whtest vectorclass)

