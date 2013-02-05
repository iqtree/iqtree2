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
#include "Parameter_tree.h"
#include "MbRandom.h"
#include "Model.h"
#include "util.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

void LambdaTable::updateNodesAtTable(Tree *t){
	
	for (set<int>::const_iterator d=diners.begin(); d != diners.end(); d++){
		int cnID = (*d);
		Node *c = t->getNodeByIndex(cnID);
		c->setNodeExpCalRate(lambda);
	}
}

double LambdaTable::getPriorPrForDiners(Tree *t){
	
	double prob = 0.0;
	double tScl = t->getTreeScale();
	for (set<int>::const_iterator d=diners.begin(); d != diners.end(); d++){
		int cnID = (*d);
		Node *c = t->getNodeByIndex(cnID);
		double nMin = c->getNodeYngTime();
		double nDep = c->getNodeDepth();
		double nDelt = (nDep * tScl) - nMin;
		prob += randP->lnExponentialPdf(lambda, nDelt);
	}
	return prob;
}

void LambdaTable::printLambdaTableInfo(){
	
	cout << "Lambda Cat (" << lambda << ") --> { ";
	for (set<int>::const_iterator d=diners.begin(); d != diners.end(); d++){
		int cnID = (*d);
		cout << cnID << " ";
	}
	cout << "}" << endl;
}

ExpCalib::ExpCalib(MbRandom *rp, Model *mp, bool dphplc, int dphpng, double ts, bool ghp) : Parameter(rp, mp){
	
	betaAlph1 = 1.0;
	betaAlph2 = 20.0;
	epsilonValue = 0.09; 
	majorityExpParm = 0.02 * ts; 
	outlieExpParm = 0.12 * ts;
	dpmLambdaExpParm = 0.07 * ts; 
	name = "EC";
	dpmLHP = dphplc;
	dpmCP = 0.001;
	gammaHPVals = ghp;
	prNGrpsDPM = dphpng;
	
	if(gammaHPVals){
		gHPAlph = 2.0;
		gHPBetaM = gHPAlph * majorityExpParm;
		gHPBetaO = gHPAlph * outlieExpParm;
		gHPBetaDP = gHPAlph * dpmLambdaExpParm;
		curMajorityLambda = ranPtr->gammaRv(gHPAlph, gHPBetaM);
		curOutlieLambda = ranPtr->gammaRv(gHPAlph, gHPBetaO);
	}
	else{
		curMajorityLambda = ranPtr->exponentialRv(majorityExpParm);
		curOutlieLambda = ranPtr->exponentialRv(outlieExpParm);
	}
	
	if(dpmLHP)
		cout << "DPM fossil calibration lambda parameter = " << dpmLambdaExpParm << endl;
	else
		cout << "Contamination model fossil calibration lambda parameters: lambda1 = " << majorityExpParm << ", lambda2 = " << outlieExpParm << endl;
}

ExpCalib::~ExpCalib(void){
	
}

ExpCalib& ExpCalib::operator=(const ExpCalib &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void ExpCalib::clone(const ExpCalib &c) {
	
	betaAlph1 = c.betaAlph1;
	betaAlph2 = c.betaAlph2;
	epsilonValue = c.epsilonValue;
	majorityExpParm = c.majorityExpParm;
	outlieExpParm = c.outlieExpParm;
	curMajorityLambda = c.curMajorityLambda;
	curOutlieLambda = c.curOutlieLambda;
	calibNodeList = c.calibNodeList;
	dpmLambdaHyp = c.dpmLambdaHyp;
	dpmLHP = c.dpmLHP;
	dpmCP = c.dpmCP;
	dpmLambdaExpParm = c.dpmLambdaExpParm;
	name = "EC";
}

void ExpCalib::print(std::ostream & o) const {
	
	if(gammaHPVals){
		if(dpmLHP){
			o << "Exponential calibration parameters (DP model):  = [expected = ";
			o << dpmLambdaExpParm  << ", gHPAlph = " << gHPAlph << ", gHPBetaDP = " << gHPBetaDP << "]";		
		}
		else{
			o << "Exponential calibration parameters (contamination model):  = [gHPAlph = ";
			o << gHPAlph  << ", expected M = " << majorityExpParm << ", gHPBetaM = " << gHPBetaM;
			o << ", expected O = " << outlieExpParm << ", gHPBetaO = " << gHPBetaO << "]";
		}
	}
	else{
		if(dpmLHP){
			o << "Exponential calibration parameters (DP model):  = [dpmLambdaExpParm = ";
			o << dpmLambdaExpParm  << "]";		
		}
		else{
			o << "Exponential calibration parameters (contamination model):  = [l1 = ";
			o << curMajorityLambda << ", l2 = " << curOutlieLambda;
			o << ", eps = " << epsilonValue << "]";
		}
	}
	o << endl;
}

double ExpCalib::update(double &oldLnL) {

	if(dpmLHP)
		updateDPMHyperPrior();
	else
		updateContamination();
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	Tree *t = modelPtr->getActiveTree();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	return 0.0;
}

void ExpCalib::updateContamination() {
	
	double tScl = modelPtr->getActiveTree()->getTreeScale();
	nodeDeltas.clear();
	for(vector<Node *>::iterator v = calibNodeList.begin(); v != calibNodeList.end(); v++){
		double nMin = (*v)->getNodeYngTime();
		double nDep = (*v)->getNodeDepth();
		double nDelt = (nDep * tScl) - nMin;
		nodeDeltas.push_back(nDelt);
	}
	
	updateMajorityLambda();
	updateOutlierLambda();
	updateNodeTaintedClassAssignment(tScl);
	
}

void ExpCalib::updateDPMHyperPrior() {
	
	Tree *t = modelPtr->getActiveTree();
	double tuning = log(2.0);
	for(vector<LambdaTable *>::iterator lt=dpmLambdaHyp.begin(); lt != dpmLambdaHyp.end(); lt++){
		double oldTL = (*lt)->getTableLambda();
		double oldPr = (*lt)->getPriorPrForDiners(t);
		double newTL = oldTL * exp(tuning * (ranPtr->uniformRv() - 0.5));
		double minV = 0.0001;
		double maxV = 1000;
		bool validV = false;
		do{
			if(newTL < minV)
				newTL = minV * minV / newTL;
			else if(newTL > maxV)
				newTL = maxV * maxV / newTL;
			else
				validV = true;
		} while(!validV);
		(*lt)->setTableLambda(newTL);
		double newPr = (*lt)->getPriorPrForDiners(t);
		
		double lpr = (newPr - oldPr) + (log(newTL) - log(oldTL));
		if(gammaHPVals)
			lpr += (ranPtr->lnGammaPdf(gHPAlph, gHPBetaDP, newPr) - ranPtr->lnGammaPdf(gHPAlph, gHPBetaDP, oldPr));
		else
			lpr += (ranPtr->lnExponentialPdf(dpmLambdaExpParm, newPr) - ranPtr->lnExponentialPdf(dpmLambdaExpParm, oldPr));
		
		double r = modelPtr->safeExponentiation(lpr);
		if(ranPtr->uniformRv() >= r){
			(*lt)->setTableLambda(oldTL);
		}
	}
	
	const int nAux = 3;
	const double lnConOvNAux = log(dpmCP / nAux);
	LambdaTable **auxilLTs = new LambdaTable*[nAux];
	for(vector<Node *>::iterator v = calibNodeList.begin(); v != calibNodeList.end(); v++){
		int calID = (*v)->getIdx();
		LambdaTable *origLT = findLambTabWCalID(calID);
		origLT->removeDiner(calID);
		if(origLT->getNumDiners() == 0)
			removeDPMLambdaTable(origLT);
		
		vector<double> lnProb;
		lnProb.reserve(dpmLambdaHyp.size() + nAux);
		for(vector<LambdaTable *>::iterator lt=dpmLambdaHyp.begin(); lt != dpmLambdaHyp.end(); lt++){
			const int nDiners = (*lt)->getNumDiners();
			(*lt)->addDiner(calID);
			double lnLTProb = getProbsAcrossAllTables(); 
			lnProb.push_back( log(nDiners) + lnLTProb );
			(*lt)->removeDiner(calID);
		}
		

		for(int i=0; i<nAux; i++){
			double auxLV;
			if(gammaHPVals)
				auxLV = ranPtr->gammaRv(gHPAlph, gHPBetaDP);
			else
				auxLV = ranPtr->exponentialRv(dpmLambdaExpParm);
			auxilLTs[i] = new LambdaTable(ranPtr, auxLV);
		}

		for(int i=0; i<nAux; i++){
			LambdaTable *tempLT = auxilLTs[i];
			dpmLambdaHyp.push_back(tempLT);
			tempLT->addDiner(calID);
			double lnLTProb = getProbsAcrossAllTables(); 
			lnProb.push_back( lnConOvNAux + lnLTProb );
			tempLT->removeDiner(calID);
			removeDPMLambdaTable(tempLT);
		}
		
		normalizeVector(lnProb);
		unsigned whichTable = ranPtr->categoricalRv(&lnProb[0], lnProb.size());
		LambdaTable *newLT = NULL;
		if(whichTable < dpmLambdaHyp.size()){
			newLT = dpmLambdaHyp[whichTable];
			newLT->addDiner(calID);
		}
		else{
			newLT = auxilLTs[whichTable-dpmLambdaHyp.size()];
			dpmLambdaHyp.push_back(newLT);
			newLT->addDiner(calID);
		}
		for(int i=0; i<nAux; i++){
			if(auxilLTs[i] != newLT)
				delete auxilLTs[i];
		}
	}
	
	delete [] auxilLTs;
	assignDPMLambdasToCals();
}



void ExpCalib::updateMajorityLambda() {
		
	double lpr = 0.0;
	double lnPriorRat = 0.0;
	double oldMajL = curMajorityLambda;
	double newMajL;	
	double tuning = log(2.0);
	newMajL = oldMajL * exp(tuning * (ranPtr->uniformRv() - 0.5));
	double minV = 0.0000001;
	double maxV = 1000000;
	bool validV = false;
	do{
		if(newMajL < minV)
			newMajL = minV * minV / newMajL;
		else if(newMajL > maxV)
			newMajL = maxV * maxV / newMajL;
		else
			validV = true;
	} while(!validV);
	lpr = log(newMajL) - log(oldMajL);

	double prNum, prDen;
	if(gammaHPVals){
		prNum = ranPtr->lnGammaPdf(gHPAlph, gHPBetaM, newMajL);
		prDen = ranPtr->lnGammaPdf(gHPAlph, gHPBetaM, oldMajL);
	}
	else{
		prNum = ranPtr->lnExponentialPdf(majorityExpParm, newMajL);
		prDen = ranPtr->lnExponentialPdf(majorityExpParm, oldMajL);
	}
	double priorValue = (prNum - prDen);
	
	for(int i=0; i<nodeDeltas.size(); i++){
		double outLV = epsilonValue * ranPtr->exponentialPdf(curOutlieLambda, nodeDeltas[i]);
		double nv = (1 - epsilonValue) * ranPtr->exponentialPdf(newMajL, nodeDeltas[i]);
		double dv = (1 - epsilonValue) * ranPtr->exponentialPdf(oldMajL, nodeDeltas[i]);
		priorValue += (log(outLV + nv) - log(outLV + dv));
	}

	lnPriorRat = priorValue;
	double lnR = lnPriorRat + lpr;
	double r = modelPtr->safeExponentiation(lnR);
	if ( ranPtr->uniformRv() < r )
		curMajorityLambda = newMajL;
	else
		curMajorityLambda = oldMajL;
	
}

void ExpCalib::updateOutlierLambda() {
		
	double lpr = 0.0;
	double lnPriorRat = 0.0;
	double oldOutL = curOutlieLambda;
	double newOutL;	
	double tuning = log(2.0);
	newOutL = oldOutL * exp(tuning * (ranPtr->uniformRv() - 0.5));
	double minV = 0.0000001;
	double maxV = 1000000;
	bool validV = false;
	do{
		if(newOutL < minV)
			newOutL = minV * minV / newOutL;
		else if(newOutL > maxV)
			newOutL = maxV * maxV / newOutL;
		else
			validV = true;
	} while(!validV);
	lpr = log(newOutL) - log(oldOutL);
	
	double prNum, prDen;
	if(gammaHPVals){
		prNum = ranPtr->lnGammaPdf(gHPAlph, gHPBetaO, newOutL);
		prDen = ranPtr->lnGammaPdf(gHPAlph, gHPBetaO, oldOutL);
	}
	else{
		prNum = ranPtr->lnExponentialPdf(outlieExpParm, newOutL);
		prDen = ranPtr->lnExponentialPdf(outlieExpParm, oldOutL);
	}
	
	double priorValue = (prNum - prDen);
	
	for(int i=0; i<nodeDeltas.size(); i++){
		double majLV = (1 - epsilonValue) * ranPtr->exponentialPdf(curMajorityLambda, nodeDeltas[i]);
		double nv = epsilonValue * ranPtr->exponentialPdf(newOutL, nodeDeltas[i]);
		double dv = epsilonValue * ranPtr->exponentialPdf(oldOutL, nodeDeltas[i]);
		priorValue += (log(nv + majLV) - log(dv + majLV));
	}
	
	lnPriorRat = priorValue;
	double lnR = lnPriorRat + lpr;
	double r = modelPtr->safeExponentiation(lnR);
	if ( ranPtr->uniformRv() < r )
		curOutlieLambda = newOutL;
	else
		curOutlieLambda = oldOutL;
	
}

void ExpCalib::updateNodeTaintedClassAssignment(double tScl){
	
	
	for(int i=0; i<calibNodeList.size(); i++){
		vector<double> lnProb;
		Node *p = calibNodeList[i];
		double nMin = p->getNodeYngTime();
		double nDep = p->getNodeDepth();
		double nDelta = (nDep * tScl) - nMin;
		double taintPr = log(epsilonValue) + ranPtr->lnExponentialPdf(curOutlieLambda, nDelta);
		double cleanPr = log(1 - epsilonValue) + ranPtr->lnExponentialPdf(curMajorityLambda, nDelta);
		lnProb.push_back(taintPr);
		lnProb.push_back(cleanPr);
		normalizeVector(lnProb);
		unsigned whichClass = ranPtr->categoricalRv(&lnProb[0], lnProb.size());
		if(whichClass == 0){ 
			p->setIsContaminatedFossil(true);
			p->setNodeExpCalRate(curOutlieLambda);
		}
		else{
			p->setIsContaminatedFossil(false);
			p->setNodeExpCalRate(curMajorityLambda);
		}
		lnProb.clear();
	}
	
}


void ExpCalib::updateEpsilonValue() {
		

	double lpr = 0.0;
	double lnPriorRat = 0.0;
	double oldEp = epsilonValue;
	double newEp = ranPtr->uniformRv();
	double prNum = ranPtr->betaPdf(betaAlph1, betaAlph2, newEp);
	double prDen = ranPtr->betaPdf(betaAlph1, betaAlph2, oldEp);

	for(int i=0; i<nodeDeltas.size(); i++){
		double outPDF = ranPtr->exponentialPdf(curOutlieLambda, nodeDeltas[i]);
		double majPDF = ranPtr->exponentialPdf(curMajorityLambda, nodeDeltas[i]);
		double nv = (newEp * outPDF) + ((1 - newEp) * majPDF);
		double dv = (oldEp * outPDF) + ((1 - oldEp) * majPDF);
		prNum += log(nv);
		prDen += log(dv);
	}
	
	lnPriorRat = prNum - prDen;
	double lnR = lnPriorRat + lpr;
	double r = modelPtr->safeExponentiation(lnR);
	if ( ranPtr->uniformRv() < r )
		epsilonValue = newEp;
	else
		epsilonValue = oldEp;
}


double ExpCalib::getLambdaForNode() {
	if(ranPtr->uniformRv() < epsilonValue)
		return curOutlieLambda;
	else
		return curMajorityLambda;
}

string ExpCalib::writeParam(void){
	
	stringstream ss;
	ss << "Exponential calibration parameters:  = [l1 = ";
	ss << curMajorityLambda << ", l2 = " << curOutlieLambda;
	ss << ", eps = " << epsilonValue << "]\n";
	string outp = ss.str();
	return outp;
}

void ExpCalib::getAllExpHPCalibratedNodes(){
	
	Tree *t = modelPtr->getActiveTree();
	calibNodeList = t->getListOfCalibratedNodes();
	if(dpmLHP)
		initializeDPMLambdasToCals();

}

bool ExpCalib::getIsLambdaContaminationClass(double l){
	
	if(l == curMajorityLambda)
		return false;
	else
		return true;
}

void ExpCalib::initializeDPMLambdasToCals(){
	
	int numCalNodes = calibNodeList.size();
	if(numCalNodes > 2){
		if(prNGrpsDPM > numCalNodes){
			cerr << "ERROR: the prior mean number of calibration clusters is > the number of calibrated nodes!" << endl;
			exit(1);
		}		
		dpmCP = calculateFromPriorMean(prNGrpsDPM, numCalNodes);
		int seated = 0;
		for(vector<Node *>::iterator v = calibNodeList.begin(); v != calibNodeList.end(); v++){
			double prNewTable = dpmCP / (seated + dpmCP);
			int ndID = (*v)->getIdx();
			if(ranPtr->uniformRv() < prNewTable){
				double lmda;
				if(gammaHPVals)
					lmda = ranPtr->gammaRv(gHPAlph, gHPBetaDP);
				else
					lmda = ranPtr->exponentialRv(dpmLambdaExpParm);
				LambdaTable *lt = new LambdaTable(ranPtr, lmda);
				dpmLambdaHyp.push_back(lt);
				lt->addDiner(ndID);
			}
			else{
				double u = ranPtr->uniformRv();
				double sum = 0.0;
				for(vector<LambdaTable *>::iterator p=dpmLambdaHyp.begin(); p != dpmLambdaHyp.end(); p++){
					sum += (double)((*p)->getNumDiners()) / seated;
					if(u < sum){
						(*p)->addDiner(ndID);
						break;
					}
				}
			}
			seated++;
		}
		assignDPMLambdasToCals();
	}
	else
		dpmLHP = false;

}



void ExpCalib::assignDPMLambdasToCals(){
	
	Tree *t = modelPtr->getActiveTree();
	for(vector<LambdaTable *>::iterator p=dpmLambdaHyp.begin(); p != dpmLambdaHyp.end(); p++){
		(*p)->updateNodesAtTable(t);
	}
	
}

LambdaTable* ExpCalib::findLambTabWCalID(int ix){
	
	LambdaTable *g = NULL;
	for(vector<LambdaTable *>::iterator lt=dpmLambdaHyp.begin(); lt != dpmLambdaHyp.end(); lt++){
		if((*lt)->getIsDinerSeated(ix)){
			g = (*lt);
			break;
		}
	}
	return g;
}


void ExpCalib::removeDPMLambdaTable(LambdaTable *g){
	
	for(vector<LambdaTable *>::iterator lt=dpmLambdaHyp.begin(); lt != dpmLambdaHyp.end(); lt++){
		if((*lt) == g){
			dpmLambdaHyp.erase(lt);
			break;
		}
	}
}


double ExpCalib::getProbsAcrossAllTables(){
	
	Tree *t = modelPtr->getActiveTree();
	double lnProbTables = 0.0;
	for(vector<LambdaTable *>::iterator lt=dpmLambdaHyp.begin(); lt != dpmLambdaHyp.end(); lt++){
		lnProbTables += (*lt)->getPriorPrForDiners(t);
	}
	return lnProbTables;
}




