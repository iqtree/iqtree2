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
 * Birth-death model from Gernhard (2008), similar to implementation in BEAST
*/

#include "Parameter.h"
#include "Parameter_speciaton.h"
#include "Parameter_treescale.h"
#include "Parameter_tree.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

Speciation::Speciation(MbRandom *rp, Model *mp, double bdr, double bda) : Parameter(rp, mp) {
		
	maxdivV = 30000.0;
	name = "SP";
	relativeDeath = bda;
	netDiversificaton = bdr;
	treeTimePrior = modelPtr->getTreeTimePriorNum();
	if(relativeDeath < 0.0 && netDiversificaton < 0.0){
		relativeDeath = ranPtr->uniformRv();
		netDiversificaton = ranPtr->uniformRv() * maxdivV;
	}
	if(relativeDeath >= 1.0){
		cerr << "ERROR: The relative death rate (-bda) is >= 1 " << endl;
		exit(1);
	}
	if(netDiversificaton <= 0.0){
		cerr << "ERROR: The net diversification rate (-bdr) is <= 0 " << endl;
		exit(1);
	}
	if(treeTimePrior == 1){
		relativeDeath = 0.0;
		netDiversificaton = 0.0;
	}	
	else if(treeTimePrior == 2)
		relativeDeath = 0.0;
	else if(treeTimePrior == 4)
		cout << "Speciaton parameters are fixed to: b/d = " << relativeDeath << " , b-d = " << netDiversificaton << endl;
	if(netDiversificaton >= maxdivV)
		netDiversificaton = maxdivV * maxdivV / netDiversificaton;
}

Speciation::~Speciation(void) {
	
}

Speciation& Speciation::operator=(const Speciation &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void Speciation::clone(const Speciation &c) {
	
	relativeDeath = c.relativeDeath;
	netDiversificaton = c.netDiversificaton;
	name = "SP";
}

void Speciation::print(std::ostream & o) const {
	if(treeTimePrior > 1){
		o << "Speciaton parameters: d/b = ";
		o << fixed << setprecision(4) << relativeDeath << " , b-d = ";
		o << fixed << setprecision(4) << netDiversificaton << " ";
		o << endl;
	}
}

double Speciation::update(double &oldLnL) {
	
	Tree *t = modelPtr->getActiveTree();
	double oldtreeprob = getLnTreeProb(t); 
	double lnProposalRatio = 0.0;
	
	
	if(treeTimePrior == 2){
		lnProposalRatio += updateNetDivRate();
		relativeDeath = 0.0;
	}
	else if(treeTimePrior == 3){
		if(ranPtr->uniformRv() < 0.5)
			lnProposalRatio += updateRelDeathRt();
		else
			lnProposalRatio += updateNetDivRate();
	}
	
	double newtreeprob = getLnTreeProb(t); 
	double lnPriorRatio = (newtreeprob - oldtreeprob);
	double lnR = lnPriorRatio + lnProposalRatio;
	
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	return lnR;
}

double Speciation::updateRelDeathRt(void) {
	
	double rdwindow = 0.2;
	double oldRD = relativeDeath;
	double u;
	double newRD;
	u = ranPtr->uniformRv(-0.5,0.5) * (rdwindow);
	newRD = oldRD + u;
	bool validV = false;
	do{
		if(newRD < 0.0)
			newRD = 0.0 - newRD;
		else if(newRD > 0.9999)
			newRD = (2 * 0.9999) - newRD;
		else
			validV = true;
	}while(!validV);
	relativeDeath = newRD;
	return 0.0;
}

double Speciation::updateNetDivRate(void) {
	
	double lpr = 0.0;
	double oldND = netDiversificaton;
	double newND;	
	double tuning = log(2.0);
	double rv = ranPtr->uniformRv();
	double c = tuning * (rv - 0.5);
	newND = oldND * exp(c);
	double minV = 0.0001;
	double maxV = maxdivV;
	bool validV = false;
	do{
		if(newND < minV)
			newND = minV * minV / newND;
		else if(newND > maxV)
			newND = maxV * maxV / newND;
		else
			validV = true;
	} while(!validV);
	netDiversificaton = newND;
	lpr = c; 
	return lpr;
}

double Speciation::lnPrior(void) {
	
	return 0.0;
}

string Speciation::writeParam(void){
	
	stringstream ss;
	ss << "Speciation parameters: m/l = " << fixed << setprecision(4) << relativeDeath << " , l-m = " 
	   << fixed << setprecision(4) << netDiversificaton << endl;
	return ss.str();
}

double Speciation::getLnTreeProb(Tree *t) {
	
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		int ntax = t->getNumTaxa();
		double c1 = (ntax - 1) * log(netDiversificaton); 
		double nps = t->getTreeCBDNodePriorProb(netDiversificaton, relativeDeath);
		return c1 + nps;
	}
	else{
		int ntax = t->getNumTaxa();
		double c1 = (ntax - 1) * log(netDiversificaton) + ntax * log(1 - relativeDeath);
		double nps = t->getTreeCBDNodePriorProb(netDiversificaton, relativeDeath);
		return c1 + nps;
	}
	return 0.0;
}

