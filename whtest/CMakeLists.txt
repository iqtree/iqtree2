
include_directories("${PROJECT_BINARY_DIR}")

add_library(whtest 
eigen.c  
eigen_sym.c  
random.c  
weisslambda.c  
weisslambda_sub.c  
whtest.c  
whtest_sub.c
whtest_wrapper.cpp whtest_wrapper.h
)

#install (TARGETS whtest DESTINATION bin)
#install (FILES  eigen.h  eigen_sym.h  random.h  tools.h  weisslambda_sub.h  whtest.h  whtest_sub.h DESTINATION include)
