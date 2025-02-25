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
alisim.cpp
alisim.h
terraceanalysis.cpp
terraceanalysis.h
)

if (USE_BOOSTER)
    target_link_libraries(main booster)
endif()

if (USE_LSD2)
    target_link_libraries(main lsd2)
endif()

if (USE_NN OR USE_OLD_NN)
    target_link_libraries(main nn)
endif()

target_link_libraries(main pda whtest vectorclass terrace)

if (USE_CMAPLE)
	target_link_libraries(main maple)
	
    if (USE_CMAPLE_AA)
        add_library(main-aa
        main.cpp
        phyloanalysis.cpp
        phyloanalysis.h
        phylotesting.cpp
        phylotesting.h
        treetesting.cpp
        treetesting.h
        timetree.cpp
        timetree.h
        alisim.cpp
        alisim.h
        terraceanalysis.cpp
        terraceanalysis.h
        )

        if (USE_BOOSTER)
                target_link_libraries(main-aa booster)
        endif()

        if (USE_LSD2)
                target_link_libraries(main-aa lsd2)
        endif()
        
        target_link_libraries(main-aa pda whtest vectorclass terrace maple-aa)
    endif()
endif()
