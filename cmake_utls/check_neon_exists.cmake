cmake_minimum_required(VERSION 2.8.12)

if (DEFINED NEON)
  return()
endif()

set(NEON 0)

enable_language(C)
enable_language(CXX)

try_compile(test_for_neon_worked ${PROJECT_BINARY_DIR}/neon_test_build ${CMAKE_CURRENT_LIST_DIR}/test_for_neon
  neon_test)

if(test_for_neon_worked)
  set(NEON 1)
endif()
