CC            = gcc
PROFILE_FLAGS = #-pg
DEBUG_FLAGS   = #-Wall -pedantic -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes  -Wdeclaration-after-statement -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef -Wbad-function-cast  
OPT_FLAGS     = -Wall -g -march=native -O2 -fomit-frame-pointer -funroll-loops 
DEFINES      += -D__SIM_SSE3 -D_GNU_SOURCE
# -D__AVX

INCLUDES = -I$(SRC_DIR)   
LIBS = -pthread
CCFLAGS = $(PROFILE_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS)
LFLAGS  = $(PROFILE_FLAGS)

