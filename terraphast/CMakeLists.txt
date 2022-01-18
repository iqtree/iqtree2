project(terraphast C CXX)

cmake_minimum_required(VERSION 3.0.2)
cmake_policy(SET CMP0054 NEW)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

Option(DEV_ENVIRONMENT "Development environment (extended debugging)" OFF)

#####################################################################
# decide which components need to be built: by default, build everything,
# unless some components were explicitely disabled from the parent project
#####################################################################
Option(TERRAPHAST_USE_GMP "Use GMP" ON)
Option(TERRAPHAST_BUILD_CLIB "Build the C library" ON)
Option(TERRAPHAST_BUILD_APPS "Build the tools" ON)
Option(TERRAPHAST_BUILD_TESTS "Build the tests" ON)
#Option(TERRAPHAST_ARCH_NATIVE "Use -march=native compiler flag" ON)
if (NEON)
	Option(TERRAPHAST_ARCH_NATIVE "Use -march=native compiler flag" OFF)
else()
    Option(TERRAPHAST_ARCH_NATIVE "Use -march=native compiler flag" ON)
endif()

#####################################################################
# build library
#####################################################################
add_library(terraphast
        lib/advanced.cpp
        lib/bigint.cpp
        lib/bipartitions.cpp
        lib/bipartitions.hpp
        lib/bitmatrix.cpp
        lib/bits.hpp
        lib/bitvector.hpp
        lib/clamped_uint.cpp
        lib/constraints.cpp
        lib/constraints_impl.hpp
        lib/errors.cpp
        lib/io_utils.hpp
        lib/multitree.cpp
        lib/multitree.hpp
        lib/multitree_impl.hpp
        lib/multitree_iterator.cpp
        lib/multitree_iterator.hpp
        lib/nodes.cpp
        lib/parser.cpp
        lib/ranked_bitvector.hpp
        lib/rooting.cpp
        lib/simple.cpp
        lib/small_bipartition.hpp
        lib/stack_allocator.hpp
        lib/subtree_extraction.cpp
        lib/subtree_extraction_impl.hpp
        lib/supertree_enumerator.hpp
        lib/supertree_helpers.cpp
        lib/supertree_helpers.hpp
        lib/supertree_variants.hpp
        lib/supertree_variants_debug.hpp
        lib/supertree_variants_multitree.hpp
        lib/trees.cpp
        lib/trees_impl.hpp
        lib/union_find.cpp
        lib/union_find.hpp
        lib/union_find_debug.hpp
        lib/utils.hpp
        lib/validation.cpp
        lib/validation.hpp
        # For QtCreator/CLion/... to show the files
        include/terraces/advanced.hpp
        include/terraces/bigint.hpp
        include/terraces/bitmatrix.hpp
        include/terraces/clamped_uint.hpp
        include/terraces/constraints.hpp
        include/terraces/definitions.hpp
        include/terraces/errors.hpp
        include/terraces/parser.hpp
        include/terraces/rooting.hpp
        include/terraces/simple.hpp
        include/terraces/subtree_extraction.hpp
        include/terraces/trees.hpp
)
target_include_directories(terraphast
        PUBLIC include
        PRIVATE lib
)

if(TERRAPHAST_USE_GMP)
    find_package(GMP)
    if(GMP_FOUND)
        message(STATUS "GMP libraries found")
        target_link_libraries(terraphast gmpxx gmp)
        target_compile_definitions(terraphast PUBLIC USE_GMP)
    else()
        message(FATAL_ERROR "GMP libraries not found! Disable them using -DTERRAPHAST_USE_GMP=OFF")
    endif()
endif()

set(terraces_targets terraphast)

if(TERRAPHAST_BUILD_CLIB)
    add_library(terraces_c
            c_lib/terraces.cpp
            c_include/terraces/terraces.h
    )
    target_include_directories(terraces_c PUBLIC c_include)
    target_link_libraries(terraces_c terraphast)
    if (NOT TERRAPHAST_USE_GMP)
        message(FATAL_ERROR "The C library requires the GMP libraries to build! Enable them using -DTERRAPHAST_USE_GMP=ON")
    endif()

    set(terraces_targets ${terraces_targets} terraces_c)
endif()

#####################################################################
# internal compiler flags
#####################################################################
if(DEV_ENVIRONMENT AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(terraphast PUBLIC _GLIBCXX_DEBUG) # PUBLIC to maintain ABI compatibility
    if(TERRAPHAST_BUILD_CLIB)
        target_compile_definitions(terraces_c PRIVATE _GLIBCXX_DEBUG) # PRIVATE since no stdlib objects are used
    endif()
endif()

#####################################################################
# build tools
#####################################################################
if(TERRAPHAST_BUILD_APPS)
    add_executable(app "app/app.cpp")
    add_executable(validated_run "tools/validated_run.cpp")
    add_executable(verbose_run "tools/verbose_run.cpp")
    add_executable(isomorphic "tools/isomorphic.cpp")
    add_executable(reroot "tools/reroot.cpp")
    add_executable(subtree "tools/subtree.cpp")
    add_executable(tree_gen "tools/tree_gen.cpp")
    add_executable(site_gen "tools/site_gen.cpp")
    add_executable(nwk_to_dot "tools/nwk_to_dot.cpp")
    target_link_libraries(validated_run terraphast)
    target_link_libraries(verbose_run terraphast)
    target_link_libraries(isomorphic terraphast)
    target_link_libraries(reroot terraphast)
    target_link_libraries(subtree terraphast)
    target_link_libraries(tree_gen terraphast)
    target_link_libraries(site_gen terraphast)
    target_link_libraries(nwk_to_dot terraphast)
    target_link_libraries(app terraphast)

    set(terraces_targets ${terraces_targets} app validated_run verbose_run isomorphic reroot subtree tree_gen site_gen nwk_to_dot)
endif()

#####################################################################
# build tests
#####################################################################
if(TERRAPHAST_BUILD_TESTS)
	add_library(Catch INTERFACE)
	target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/catch)

	add_executable(unittests
		test/main.cc
		test/advanced.cpp
		test/bipartitions.cpp
		test/bitmatrix.cpp
		test/bits.cpp
		test/bitvector.cpp
		test/clamped_uint.cpp
		test/constraints.cpp
		test/compile_test.cpp
		test/fast_set.cpp
		test/integration.cpp
		test/multitree_iterator.cpp
		test/parser.cpp
		test/rooting.cpp
		test/small_bipartition.cpp
		test/stack_allocator.cpp
		test/subtree_extraction.cpp
		test/supertree.cpp
		test/trees.cpp
		test/union_find.cpp
		test/util.cpp
		test/validation.cpp
		test/simple.cpp
	)
	target_link_libraries(unittests terraphast Catch)
	if(TERRAPHAST_BUILD_CLIB)
		target_sources(unittests PRIVATE
			test/c_api.cpp
		)
		target_link_libraries(unittests terraces_c)
	endif()
	add_test(NAME unittests COMMAND unittests)
	enable_testing()

	set(terraces_targets ${terraces_targets} unittests)
endif()

set_target_properties(${terraces_targets} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)

#####################################################################
# set platform-specific options, include platform-specific files
#####################################################################
if(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC") OR CLANG_UNDER_VS)
	set(TERRAPHAST_PLATFORM_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cl")
	set(TERRAPHAST_COMPILE_FLAGS -Oi -W4)

	# Most of our files only compile with disabled language extensions for VC++
	# Unfortunately, Catch uses some windows-specific features, so we have to
	# enable these extensions for the Catch main method (more specific: not disable them)
	file(GLOB ALL_SOURCES lib/*.cpp c_lib/*.cpp test/*.cpp tools/*.cpp app/*.cpp)
    if (NOT CLANG_UNDER_VS)
	    set_source_files_properties(${ALL_SOURCES} PROPERTIES COMPILE_FLAGS "-Za")
    endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(TERRAPHAST_PLATFORM_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/lib/intel")
    set(TERRAPHAST_COMPILE_FLAGS -Wall -Wextra -Wconversion -Wsign-conversion)
    if(TERRAPHAST_ARCH_NATIVE)
        set(TERRAPHAST_COMPILE_FLAGS -march=native ${TERRAPHAST_COMPILE_FLAGS})
    endif()
else()
    set(TERRAPHAST_PLATFORM_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/lib/gcc_clang")
    set(TERRAPHAST_COMPILE_FLAGS -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion)
    if(TERRAPHAST_ARCH_NATIVE)
        set(TERRAPHAST_COMPILE_FLAGS -march=native ${TERRAPHAST_COMPILE_FLAGS})
    endif()
endif()

target_include_directories(terraphast PUBLIC "${TERRAPHAST_PLATFORM_INCLUDE}")
if(TERRAPHAST_BUILD_TESTS)
    target_include_directories(unittests PRIVATE "${TERRAPHAST_PLATFORM_INCLUDE}")
endif()
target_compile_options(terraphast PRIVATE "${TERRAPHAST_COMPILE_FLAGS}")
if(TERRAPHAST_BUILD_CLIB)
    target_compile_options(terraces_c PRIVATE "${TERRAPHAST_COMPILE_FLAGS}")
endif()

#####################################################################
# setup clang-tidy
#####################################################################
find_program(CLANG_TIDY_PATH NAMES "clang-tidy")
if(NOT CLANG_TIDY_PATH)
    message(STATUS "clang-tidy not found.")
else()
    message(STATUS "clang-tidy found: ${CLANG_TIDY_PATH}")
    set(CLANG_TIDY_COMMANDLINE "${CLANG_TIDY_PATH}" "-checks=llvm-namespace-comment -fix")
endif()

#set_target_properties(terraphast terraces_c app validated_run verbose_run isomorphic reroot subtree tree_gen site_gen nwk_to_dot unittests PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMANDLINE}")
