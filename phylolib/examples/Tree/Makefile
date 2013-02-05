CC           = gcc
CPP          = g++ 
BUILD_DIR    = ./
SRC_DIR      = ../../
LIB_CFILES   = bipartitionList.c evaluateGenericSpecial.c evaluatePartialGenericSpecial.c makenewzGenericSpecial.c models.c optimizeModel.c restartHashTable.c searchAlgo.c topologies.c trash.c treeIO.c newviewGenericSpecial.c fastDNAparsimony.c randomTree.c recom.c phylip.c msa_sites.c lexer.c ssort.c xalloc.c utils.c
LIB_SRC      = $(patsubst %.c, $(SRC_DIR)%.c, $(LIB_CFILES))
PHYLIP_DIR   = ../../phylip_parser/
LIB_OBJ      = $(patsubst %.c, %.o, $(LIB_CFILES))
CFLAGS       = -I$(SRC_DIR) -c -Wall -g -march=native -O2 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -DEXPERIMENTAL -D_OPTIMIZED_FUNCTIONS -msse3 -g
CXXFLAGS       = -I$(SRC_DIR) -c -Wall -g -march=native -O2 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS -msse3 -g 

all: read nni read2

nni: $(LIB_OBJ) nni.c
	@echo "===> COMPILE $@"
	$(CC) -g -I$(SRC_DIR) $+ -o $@ -lm

read: $(LIB_OBJ) read.c
	@echo "===> COMPILE $@"
	$(CC) -g -I$(SRC_DIR) $+ -o $@ -lm
	
%.o: $(SRC_DIR)%.c 
	$(CC) $(CFLAGS) -c $+ -o $@

%.o: $(PHYLIP_DIR)%.c
	$(CC) $(CFLAGS) -c $+ -o $@

clean:
	-rm nni read *.o

.PHONY: all clean
