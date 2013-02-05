/* 
 * DPPDiv version 1.0b source code (git: 9c0ac3d2258f89827cfe9ba2b5038f0f656b82c1)
 * Copyright 2009-2011
 * Tracy Heath(1,2,3) (NSF postdoctoral fellowship in biological informatics DBI-0805631)
 * Mark Holder(1)
 * John Huelsenbeck(2)
 *
 * (1) Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS 66045
 * (2) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (3) email: tracyh@berkeley.edu
 *
 * DPPDiv is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * Some of this code is from publicly available source by John Huelsenbeck
 */

#include "Parameter.h"
#include "Parameter_cphyperp.h"
#include "Parameter_rate.h"
#include "Parameter_tree.h"
#include "MbRandom.h"
#include "Model.h"
#include "util.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;


Cphyperp::Cphyperp(MbRandom *rp, Model *mp, double ga, double gb, int nn, int pm, bool fixcp) : Parameter(rp, mp) {
	
	numNodes = nn - 1; 
	double expA = calculateFromPriorMean(pm, numNodes);
	if(ga > 0.0 && gb > 0.0){
		gammaAlpha = ga;
		gammaBeta = gb;
		expA = gammaAlpha / gammaBeta;
	}
	else{
		gammaAlpha = 2.0;
		gammaBeta = gammaAlpha / expA;
	}
	if(fixcp)
		currentCP = expA;
	else
		currentCP = ranPtr->gammaRv(gammaAlpha, gammaBeta);
	
	cout << "Expected cp = " << expA << endl;
	cout << "Initialized = " << currentCP << endl;
	name = "CP";
}

Cphyperp::~Cphyperp(void) {
	
}

Cphyperp& Cphyperp::operator=(const Cphyperp &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void Cphyperp::clone(const Cphyperp &c) {
	
	gammaAlpha = c.gammaAlpha;
	gammaBeta = c.gammaBeta;
	currentCP = c.currentCP;
	numNodes = c.numNodes;
}

void Cphyperp::print(std::ostream & o) const {
	
	o << "Concentration parameter: ";
	o << fixed << setprecision(4) << currentCP << " ";
	o << endl;
}

double Cphyperp::update(double &oldLnL) {
	
	NodeRate *nr = modelPtr->getActiveNodeRate();
	int k = nr->getNumRateGroups();
	int n = numNodes;
	double oldAlpha = nr->getConcenParam();
	
	MbVector<double> z(2);
	MbVector<double> f(2);
	z[0] = oldAlpha + 1.0;
	z[1] = (double)n;
	ranPtr->dirichletRv(z, f);
	double eta = f[0];
	
	double u = ranPtr->uniformRv();
	double x = ( gammaAlpha + (double)k - 1.0 ) / ( (double)n * (gammaBeta - log(eta)) );
	double newAlpha;
	if((u / (1.0 - u)) < x)
		newAlpha = ranPtr->gammaRv(gammaAlpha + k, gammaBeta - log(eta));
	else
		newAlpha = ranPtr->gammaRv(gammaAlpha + k - 1.0, gammaBeta - log(eta));

	currentCP = newAlpha;
	nr->setConcenParam(newAlpha);
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	Tree *t = modelPtr->getActiveTree();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	return 0.0;
}

double Cphyperp::lnPrior(void) {
	
	return 0.0;
}

string Cphyperp::writeParam(void){
	
	stringstream ss;
	ss << "Concentration parameter: " << currentCP << endl;
	string outp = ss.str();
	return outp;
}
