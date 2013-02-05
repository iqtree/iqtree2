CC           = gcc
CPP          = g++
BUILD_DIR    = ./
SRC_DIR      = ../../
LIB_CFILES   = bipartitionList.c evaluateGenericSpecial.c evaluatePartialGenericSpecial.c makenewzGenericSpecial.c models.c optimizeModel.c restartHashTable.c searchAlgo.c topologies.c trash.c treeIO.c newviewGenericSpecial.c fastDNAparsimony.c randomTree.c recom.c utils.c
LIB_SRC      = $(patsubst %.c, $(SRC_DIR)%.c, $(LIB_CFILES))
LIB_OBJ      = $(patsubst %.c, %.o, $(LIB_CFILES))
CFLAGS       = -I$(SRC_DIR) -c -Wall -g -march=native -O2 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS -msse3 -g 
CXXFLAGS       = -I$(SRC_DIR) -c -Wall -g -march=native -O2 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS -msse3 -g 

all: msa likelihood ancestral phylip

ancestral: $(LIB_OBJ) ancestral.c
	@echo "===> COMPILE $@"
	$(CC) -I$(SRC_DIR) $+ -o $@ -lm

likelihood: $(LIB_OBJ) likelihood.c
	@echo "===> COMPILE $@"
	$(CC) -I$(SRC_DIR) $+ -o $@ -lm

msa: $(LIB_OBJ) msa.c
	@echo "===> COMPILE $@"
	$(CC) -I$(SRC_DIR) $+ -o $@ -lm -g

phylip: $(LIB_OBJ) phylip_parser.c
	@echo "===> COMPILE $@"
	$(CC) -I$(SRC_DIR) $+ -o $@ -lm -g

%.o: $(SRC_DIR)%.c
	$(CC) $(CFLAGS) -c $+ -o $@


clean:
	-rm msa likelihood ancestral phylip *.o

.PHONY: all clean
