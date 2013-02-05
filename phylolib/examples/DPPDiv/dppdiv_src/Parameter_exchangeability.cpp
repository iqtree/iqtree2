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

#include "MbRandom.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_exchangeability.h"
#include "Parameter_tree.h"
#include <iostream>
#include <iomanip>

using namespace std;

Exchangeability::Exchangeability(MbRandom *rp, Model *mp) : Parameter(rp, mp) {

	rates = MbVector<double>(6);
	alpha = MbVector<double>(6);
	for (int i=0; i<6; i++)
		alpha[i] = 1.0;
	ranPtr->dirichletRv(alpha, rates);
	alpha0 = 800.0;
	name = "RM";
}

Exchangeability::~Exchangeability(void) {

}

Exchangeability& Exchangeability::operator=(const Exchangeability &b) {

	if (this != &b)
		clone(b);
	return *this;
}

void Exchangeability::clone(const Exchangeability &b) {

	for (int i=0; i<6; i++)
		rates[i] = b.rates[i];
}

void Exchangeability::print(std::ostream & o) const {

	o << "Substitution Rates: ";
	for (int i=0; i<6; i++)
		o << fixed << setprecision(4) << rates[i] << " ";
	o << endl;
}

double Exchangeability::update(double &oldLnL) {

	MbVector<double> aForward(6);
	MbVector<double> aReverse(6);
	MbVector<double> oldRates(6);
	MbVector<double> newRates(6);
	
	for (int i=0; i<6; i++){
		oldRates[i] = rates[i];
		aForward[i] = rates[i] * alpha0;
	}
		
	ranPtr->dirichletRv(aForward, newRates);
	
	double sum = 0.0;
	for(int i=0; i<6; i++){
		if(newRates[i] < 0.000001)
			newRates[i] = 0.000001;
		sum += newRates[i];
	}
	for(int i=0; i<6; i++)
		newRates[i] /= sum;	
	
	for (int i=0; i<6; i++)
		rates[i] = newRates[i];
	
	for (int i=0; i<6; i++)
		aReverse[i] = newRates[i] * alpha0;
		
	double lnProposalRatio = ranPtr->lnDirichletPdf(aReverse, oldRates) - ranPtr->lnDirichletPdf(aForward, newRates);
	Tree *t = modelPtr->getActiveTree();
	t->flipAllCls();
	t->flipAllTis();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->upDateRateMatrix();
	modelPtr->setTiProb();

	return lnProposalRatio;
}

double Exchangeability::lnPrior(void) {

	return ranPtr->lnDirichletPdf(alpha, rates);
}

string Exchangeability::writeParam(void){
	
	stringstream ss;
	ss << "Substitution Rates: ";
	for (int i=0; i<6; i++)
		ss << fixed << setprecision(4) << rates[i] << " ";
	ss << endl;
	string outp = ss.str();
	return outp;
}

