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

void addArg(int &argc, char **argv, char *arg) {
	argv[argc] = new char[strlen(arg)+1];
	strcpy(argv[argc], arg);
	argc++;
}

int WHTest_old(Params &params, PhyloTree &tree) {
	int argc = 0;
	char *argv[10];
	char tmp[100];
	addArg(argc, argv, (char*)"WHTest");
	addArg(argc, argv, params.aln_file);
	addArg(argc, argv, (char*)"-a");
	sprintf(tmp, "%f", tree.getModelFactory()->site_rate->getGammaShape());
	addArg(argc, argv, tmp);
	return WHTest_run(argc, argv);
}

int WHTest(Params &params, IQTree &tree) {

	int retval;
	size_t nseq = tree.aln->getNSeq();
    size_t nsite = tree.aln->getNSite(); 


	WHT_setAlignmentSize(nseq, nsite);
	WHT_allocateMemory();
	for (size_t i = 0; i < nseq; i++)
		for (size_t j = 0; j < nsite; j++)
			WHT_setSequenceSite(i, j, (*tree.aln)[tree.aln->getPatternID(j)][i]);
			
	for (size_t i = 0; i < nseq; i++)
		WHT_setSequenceName(i, tree.aln->getSeqName(i).c_str());
	double gamma_shape = tree.getModelFactory()->site_rate->getGammaShape();
	if (gamma_shape == 0) gamma_shape = 100.0;
	//WHT_setParams(params.whtest_simulations, gamma_shape, params.out_prefix, tree.dist_matrix);
	WHT_setParams(params.whtest_simulations, gamma_shape, params.out_prefix, NULL);
	retval = WHTest_run(0, NULL);
	WHT_getResults(&params.whtest_delta, &params.whtest_delta_quantile, &params.whtest_p_value);
	return retval;
}
