add_library(main
main.cpp
phyloanalysis.cpp
phyloanalysis.h
phylotesting.cpp
phylotesting.h
treetesting.cpp
treetesting.h
)

if (USE_BOOSTER)
    target_link_libraries(main booster)
endif()

target_link_libraries(main pda whtest vectorclass)

