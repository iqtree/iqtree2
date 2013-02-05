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
#include "Parameter_expcalib.h"
#include "Parameter_speciaton.h"
#include "Parameter_treescale.h"
#include "Parameter_tree.h"
#include "MbMath.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

Treescale::Treescale(MbRandom *rp, Model *mp, double sv, double yb, double ob, int dt, bool calib, bool exhpc) : Parameter(rp, mp) {
	oldBound = ob;
	yngBound = yb;
	scaleVal = sv;
	name = "RH";
	treeTimePrior = modelPtr->getTreeTimePriorNum();
	retune = false;
	tsPriorD = dt;
	isUserCalib = calib;
	exponHPCalib = exhpc;
	expoRate = 1.0;
	
	tuning = sv * 0.1;
	if(tsPriorD == 1){
		exponHPCalib = false;
		if(yngBound > 0.0 && oldBound > 0.0){
			isBounded = true;
			if(yngBound != oldBound){
				tuning = ((oldBound + yngBound) * 0.5) * 0.1;
				double windowSize = 2 * tuning;
				double calibSize = oldBound - yngBound;
				while(windowSize > calibSize * 0.85){
					tuning = tuning * 0.9;
					windowSize = 2 * tuning;
				}
				cout << "Calib window = " << calibSize << "  Tuning window = " << windowSize << endl;
			}
		}
	}
	else if(tsPriorD == 2){
		isBounded = true;

		if(isUserCalib)
			expoRate = modelPtr->getRootNExpRate();
		else
			expoRate = 1.0 / (yngBound * 2.0);
		
		oldBound = 1000000.0;
		tuning = ((yngBound + (sv * 1.2)) * 0.5) * 0.1;
		double windowSize = 2 * tuning;
		cout << "Exponential rate = " << expoRate << "  E[TS] = " << (1 / expoRate) + yngBound
			<< "  offset = " << yngBound << "  Tuning window = " << windowSize << endl;
	}
}

Treescale::~Treescale(void) {
	
}

Treescale& Treescale::operator=(const Treescale &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void Treescale::clone(const Treescale &c) {

	scaleVal = c.scaleVal;
	oldBound = c.oldBound;
	yngBound = c.yngBound;
	isBounded = c.isBounded;
	tuning = c.tuning;
	name = "RH";
}

void Treescale::print(std::ostream & o) const {

	o << "Root height parameter: ";
	o << fixed << setprecision(4) << scaleVal << " ";
	o << endl;
}

double Treescale::update(double &oldLnL) {
		
	Tree *t = modelPtr->getActiveTree();
	Node *rt = t->getRoot();
	
	double oldtreeprob = getLnTreeProb(t);
	
	if(retune && tuning > scaleVal)
		tuning = scaleVal * 0.5;
		
	double limO = scaleVal + tuning;
	double limY = scaleVal - tuning;
	double lftHt = rt->getLft()->getNodeDepth() * scaleVal;
	double rhtHt = rt->getRht()->getNodeDepth() * scaleVal;
	double lowBound = lftHt;
	if(lowBound < rhtHt)
		lowBound = rhtHt;
	double hiBound = limO;
	
	if(isBounded){
		if(lowBound < yngBound)
			lowBound = yngBound;
		if(hiBound > oldBound)
			hiBound = oldBound;
	}
	
	double oldRH, newRH;
	oldRH = scaleVal;
		
	double u = ranPtr->uniformRv(-0.5,0.5) * (limO - limY);
	newRH = oldRH + u;
	while(newRH < lowBound || newRH > hiBound){
		if(newRH < lowBound)
			newRH = (2 * lowBound) - newRH;
		if(newRH > hiBound)
			newRH = (2 * hiBound) - newRH;
	}
		
	double scaleRatio = oldRH / newRH;
	int numNodes = t->getNumNodes();
	for(int i=0; i<numNodes; i++){
		Node *p = t->getNodeByIndex(i);
		if(p->getIsLeaf() == false && p != rt){
			double oldP = p->getNodeDepth();
			double newP = oldP * scaleRatio;
			p->setNodeDepth(newP);
		}
	}
	
	scaleVal = newRH;
	t->setTreeScale(scaleVal);
	t->setAllNodeBranchTimes();
		
	double lnPriorRatio = 0.0; 
	double newtreeprob = getLnTreeProb(t);
	lnPriorRatio += (newtreeprob - oldtreeprob);
	if(tsPriorD == 2){
		if(exponHPCalib)
			expoRate = t->getRootCalibExpRate();
		lnPriorRatio += lnExponentialTSPriorRatio(newRH, oldRH);
	}
	double lnProposalRatio = 0.0; 
	
	double jacobian = 0.0;
	if(treeTimePrior < 2)
		jacobian = (log(oldRH) - log(newRH)) * (t->getNumTaxa() - 2);
	
	t->flipAllCls();
	t->flipAllTis();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	
	return lnPriorRatio + lnProposalRatio + jacobian;
}

double Treescale::lnPrior(void) {
	
	return 0.0;
}

double Treescale::lnExponentialTSPriorRatio(double newTS, double oldTS) {
	
	double offSt = yngBound;
	double numr = 0.0;
	double dnom = 0.0;
	numr = ranPtr->lnExponentialPdf(expoRate, newTS - offSt);
	dnom = ranPtr->lnExponentialPdf(expoRate, oldTS - offSt);
	return numr - dnom;
}


string Treescale::writeParam(void){
	
	stringstream ss;
	ss << "Root height parameter: " << scaleVal << " [+/- " << tuning << "]" << endl;
	return ss.str();
}

double Treescale::getLnTreeProb(Tree *t) {
	
	if(treeTimePrior == 1)
		return 0.0;
	else		
		return t->getTreeCBDNodePriorProb();
	return 0.0;
}






