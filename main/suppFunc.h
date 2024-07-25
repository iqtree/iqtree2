/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef SUPPFUNC_H
#define SUPPFUNC_H

#include <stdio.h>
#include "tree/phylotree.h"
#include <signal.h>
#include <cstdio>
#include <streambuf>
#include <iostream>
#include <cstdlib>
#include <errno.h>
#include "pda/greedy.h"
#include "pda/pruning.h"
#include "pda/splitgraph.h"
#include "pda/circularnetwork.h"
#include "tree/mtreeset.h"
#include "tree/mexttree.h"
#include "ncl/ncl.h"
#include "nclextra/msetsblock.h"
#include "nclextra/myreader.h"
#include "main/phyloanalysis.h"
#include "main/alisim.h"
#include "tree/matree.h"
#include "obsolete/parsmultistate.h"
#include "alignment/maalignment.h" //added by MA
#include "tree/ncbitree.h"
#include "pda/ecopd.h"
#include "tree/upperbounds.h"
#include "main/terraceanalysis.h"
#include "pda/ecopdmtreeset.h"
#include "pda/gurobiwrapper.h"
#include "utils/timeutil.h"
#include "utils/operatingsystem.h" //for getOSName()
#include <stdlib.h>
#include "vectorclass/instrset.h"

using namespace std;

inline void separator(ostream &out, int type = 0);


void printCopyright(ostream &out);

void printRunMode(ostream &out, RunMode run_mode);

/**
    summarize the running with header
*/
void summarizeHeader(ostream &out, Params &params, bool budget_constraint, InputType analysis_type);

void summarizeFooter(ostream &out, Params &params);

int getMaxNameLen(vector<string> &setName);

void printPDUser(ostream &out, Params &params, PDRelatedMeasures &pd_more);

void summarizeTree(Params &params, PDTree &tree, vector<PDTaxaSet> &taxa_set,
                   PDRelatedMeasures &pd_more);

void printTaxaSet(Params &params, vector<PDTaxaSet> &taxa_set, RunMode cur_mode);

/**
    run PD algorithm on trees
*/
void runPDTree(Params &params);

void checkSplitDistance(ostream &out, PDNetwork &sg);

/**
    check if the set are nested and there are no multiple optimal sets.
    If yes, return the ranking as could be produced by a greedy algorithm
*/
bool makeRanking(vector<SplitSet> &pd_set, IntVector &indices, IntVector &ranking);

void printNexusSets(const char *filename, PDNetwork &sg, vector<SplitSet> &pd_set);

void computeTaxaFrequency(SplitSet &taxa_set, DoubleVector &freq);
/**
    summarize the running results
*/
void summarizeSplit(Params &params, PDNetwork &sg, vector<SplitSet> &pd_set, PDRelatedMeasures &pd_more, bool full_report);

void printGainMatrix(char *filename, mmatrix(double) &delta_gain, int start_k);

/**
    run PD algorithm on split networks
*/
void runPDSplit(Params &params);

void printSplitSet(SplitGraph &sg, SplitIntMap &hash_ss);

void readTaxaOrder(char *taxa_order_file, StrVector &taxa_order);

void calcTreeCluster(Params &params);

void printTaxa(Params &params);

void printAreaList(Params &params);

void scaleBranchLength(Params &params);

void calcDistribution(Params &params);

void printRFDist(string filename, double *rfdist, int n, int m, int rf_dist_mode, bool print_msg = true);

void computeRFDistExtended(const char *trees1, const char *trees2, const char *filename);

void computeRFDistSamePair(const char *trees1, const char *trees2, const char *filename);

void computeRFDist(Params &params);


void testInputFile(Params &params);

/**MINH ANH: for some statistics about the branches on the input tree*/
void branchStats(Params &params);

/**MINH ANH: for comparison between the input tree and each tree in a given set of trees*/
void compare(Params &params);

/**MINH ANH: to compute 'guided bootstrap' alignment*/
void guidedBootstrap(Params &params);

/**MINH ANH: to compute the probability of an alignment given the multinomial distribution of patterns frequencies derived from a reference alignment*/
void computeMulProb(Params &params);

void processNCBITree(Params &params);

/* write simultaneously to cout/cerr and a file */
class outstreambuf : public streambuf {
public:
    outstreambuf* open( const char* name, ios::openmode mode = ios::out);
    bool is_open();
    outstreambuf* close();
    ~outstreambuf() { close(); }
    streambuf *get_fout_buf() {
        return fout_buf;
    }
    streambuf *get_cout_buf() {
        return cout_buf;
    }
    ofstream *get_fout() {
        return &fout;
    }
    
protected:
    ofstream fout;
    streambuf *cout_buf;
    streambuf *fout_buf;
    virtual int     overflow( int c = EOF);
    virtual int     sync();
};

class errstreambuf : public streambuf {
public:
    void init(streambuf *fout_buf) {
        this->fout_buf = fout_buf;
        cerr_buf = cerr.rdbuf();
        cerr.rdbuf(this);
        new_line = true;
    }
    
    void reset() {
        cerr.rdbuf(cerr_buf);
    }
    
    ~errstreambuf() {
        cerr.rdbuf(cerr_buf);
    }
    
protected:
    streambuf *cerr_buf;
    streambuf *fout_buf;
    bool new_line;
    
    virtual int overflow( int c = EOF) {
        if (new_line)
            cerr_buf->sputn("ERROR: ", 7);
        if (cerr_buf->sputc(c) == EOF) {
            new_line = false;
            if (c == '\n') new_line = true;
            return EOF;
        }
        if ((Params::getInstance().suppress_output_flags & OUT_LOG)) {
            new_line = false;
            if (c == '\n') new_line = true;
            return c;
        }
        if (new_line)
            fout_buf->sputn("ERROR: ", 7);
        new_line = false;
        if (c == '\n') new_line = true;
        if (fout_buf->sputc(c) == EOF) return EOF;
        return c;
    }
    
    virtual int sync() {
        cerr_buf->pubsync();
        if (Params::getInstance().suppress_output_flags & OUT_LOG)
            return 0;
        return fout_buf->pubsync();
    }
};

class muststreambuf : public streambuf {
public:
    void init(streambuf *cout_buf, streambuf *fout_buf) {
        this->fout_buf = fout_buf;
        this->cout_buf = cout_buf;
    }
    
protected:
    streambuf *cout_buf;
    streambuf *fout_buf;
    
    virtual int overflow( int c = EOF) {
        if (cout_buf->sputc(c) == EOF) {
            return EOF;
        }
        if (fout_buf->sputc(c) == EOF) return EOF;
        return c;
    }
    
    virtual int sync() {
        cout_buf->pubsync();
        return fout_buf->pubsync();
    }
};

inline outstreambuf _out_buf;
inline errstreambuf _err_buf;
inline muststreambuf _must_buf;
inline ostream cmust(&_must_buf);
inline string _log_file;
inline int _exit_wait_optn = FALSE;

extern "C" void startLogFile(bool append_log);

extern "C" void endLogFile();

void funcExit(void);

extern "C" void funcAbort(int signal_number);

extern "C" void getintargv(int *argc, char **argv[]);

/*********************************************************************************************************************************
    Olga: ECOpd - phylogenetic diversity with ecological constraint: choosing a viable subset of species which maximizes PD/SD
*********************************************************************************************************************************/

void processECOpd(Params &params);

void collapseLowBranchSupport(char *user_file, char *split_threshold_str);

#endif /* SUPPFUNC_H */
