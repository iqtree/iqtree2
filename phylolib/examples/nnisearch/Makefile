CC           = gcc
CPP          = g++
BUILD_DIR    = ./
SRC_DIR      = ../../
LIB_CFILES   = bipartitionList.c evaluateGenericSpecial.c evaluatePartialGenericSpecial.c makenewzGenericSpecial.c models.c optimizeModel.c restartHashTable.c searchAlgo.c topologies.c trash.c treeIO.c newviewGenericSpecial.c fastDNAparsimony.c randomTree.c recom.c utils.c
LIB_SRC      = $(patsubst %.c, $(SRC_DIR)%.c, $(LIB_CFILES))
LIB_OBJ      = $(patsubst %.c, %.o, $(LIB_CFILES))
CFLAGS       = -I$(SRC_DIR) -std=c99 -c -Wall -g -march=native -O2 -D_GNU_SOURCE -msse3 -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS
CXXFLAGS     = -I$(SRC_DIR) -std=c99 -c -Wall -g -march=native -O2 -D_GNU_SOURCE -msse3 -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS

ifeq ($(mode), debug)
	CFLAGS       = -I$(SRC_DIR) -c -Wall -g -march=native -O0 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS -msse3 -g 
	CXXFLAGS       = -I$(SRC_DIR) -c -Wall -g -march=native -O0 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS -msse3 -g 
endif

all: nnisearch 

nnisearch: $(LIB_OBJ) nnisearch.c
	@echo "===> COMPILE $@"
	$(CC) -I$(SRC_DIR) $+ -g -o $@ -lm

%.o: $(SRC_DIR)%.c
	$(CC) $(CFLAGS) -c $+ -o $@


clean:
	-rm nnisearch *.o

.PHONY: all clean
