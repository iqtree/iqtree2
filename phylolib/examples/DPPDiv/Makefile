CPP = g++
CXXFLAGS = -DPHYLO_LIBRARY -c -O2
CC  = gcc
DPPDIV_SRC_DIR = dppdiv_src/
PHYLOLIB_SRC_DIR = ../../
OBJ_DIR = obj/
CXXFILES = Alignment.cpp Calibration.cpp dppdiv.cpp MbEigensystem.cpp MbMath.cpp MbRandom.cpp MbTransitionMatrix.cpp Mcmc.cpp Model.cpp Parameter_basefreq.cpp Parameter_cphyperp.cpp Parameter.cpp Parameter_exchangeability.cpp Parameter_expcalib.cpp Parameter_rate.cpp Parameter_shape.cpp Parameter_speciaton.cpp Parameter_tree.cpp Parameter_treescale.cpp
LIB_CFILES   = bipartitionList.c evaluateGenericSpecial.c evaluatePartialGenericSpecial.c makenewzGenericSpecial.c models.c optimizeModel.c restartHashTable.c searchAlgo.c topologies.c trash.c treeIO.c newviewGenericSpecial.c fastDNAparsimony.c randomTree.c recom.c utils.c
PHYLOLIB_OBJ = $(patsubst %.c, $(OBJ_DIR)%.o, $(LIB_CFILES))
DPPDIV_OBJ   = $(patsubst %.cpp, $(OBJ_DIR)%.o, $(CXXFILES))
CFLAGS       = -I$(PHYLOLIB_SRC_DIR) -c -Wall -g -march=native -O2 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS -msse3 -g

all: phylolib dppdiv globals
	$(CPP) $(DPPDIV_OBJ) $(PHYLOLIB_OBJ) $(OBJ_DIR)globals.o -o dppdiv

clean:
	-rm $(OBJ_DIR)*.o
	-rm dppdiv

globals: $(OBJ_DIR)globals.o

$(OBJ_DIR)globals.o: $(DPPDIV_SRC_DIR)globals.c
	$(CC) -I../../ -c $+ -o $(OBJ_DIR)/globals.o

phylolib: $(PHYLOLIB_OBJ)
	@echo "==> Finished compiling PHYLOLIB files"

dppdiv: $(DPPDIV_OBJ)
	@echo "==> Finished compiling DPPDiv files"

$(OBJ_DIR)dppdiv.o: $(DPPDIV_SRC_DIR)dppdiv.cpp
	$(CPP) $(CXXFLAGS) $+ -I$(DPPDIV_SRC_DIR) -I../../ -o $@

$(OBJ_DIR)%.o: $(DPPDIV_SRC_DIR)%.cpp
	$(CPP) $(CXXFLAGS) $+ -I$(DPPDIV_SRC_DIR) -o $@

$(OBJ_DIR)%.o: $(PHYLOLIB_SRC_DIR)%.c
	$(CC) $(CFLAGS) -c $+ -o $@


.PHONY: globals phylolib dppdiv
