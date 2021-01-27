/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
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


#include "whtest_wrapper.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "whtest/whtest.h"

#ifdef __cplusplus
}
#endif

void addArg(int &argc, char **argv, const char *arg) {
	size_t len = strlen(arg);
	argv[argc] = new char[len+1];
#ifndef _MSC_VER
	strcpy(argv[argc], arg);
#else
	strcpy_s(argv[argc], len + 1, arg);
#endif
	argc++;
}

#ifdef _MSC_VER
#define format_string sprintf_s
#else
#define format_string sprintf
#endif

int WHTest_old(Params &params, PhyloTree &tree) {
	int argc = 0;
	char *argv[10];
	char tmp[100];
	addArg(argc, argv, "WHTest");
	addArg(argc, argv, params.aln_file);
	addArg(argc, argv, "-a");
	format_string(tmp, "%f", tree.getModelFactory()->site_rate->getGammaShape());
	addArg(argc, argv, tmp);
	return WHTest_run(argc, argv);
}

int WHTest(Params &params, IQTree &tree) {

	int retval;
	int nseq  = tree.aln->getNSeq32();
    int nsite = tree.aln->getNSite32(); 


	WHT_setAlignmentSize(nseq, nsite);
	WHT_allocateMemory();
	for (int i = 0; i < nseq; i++) {
		for (int j = 0; j < nsite; j++) {
			WHT_setSequenceSite(i, j, (*tree.aln)[tree.aln->getPatternID(j)][i]);
		}
	}			
	for (int i = 0; i < nseq; i++) {
		WHT_setSequenceName(i, tree.aln->getSeqName(i).c_str());
	}
	double gamma_shape = tree.getModelFactory()->site_rate->getGammaShape();
	if (gamma_shape == 0) {
		gamma_shape = 100.0;
	}
	//WHT_setParams(params.whtest_simulations, gamma_shape, params.out_prefix, tree.dist_matrix);
    WHT_setParams(static_cast<int>(params.whtest_simulations), gamma_shape, params.out_prefix.c_str(), NULL);
	retval = WHTest_run(0, NULL);
	WHT_getResults(&params.whtest_delta, &params.whtest_delta_quantile, &params.whtest_p_value);
	return retval;
}
