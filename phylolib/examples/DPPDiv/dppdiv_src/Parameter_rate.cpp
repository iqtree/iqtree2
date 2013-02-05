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
#include "Parameter_rate.h"
#include "Parameter_tree.h"
#include "util.h"
#include <iostream>
#include <string>

using namespace std;

void RateGroup::print(std::ostream &ss) const {
	
	ss << "Rate Group Elements: (" << rate << ") -> ";
	for (set<int>::const_iterator p=rateElements.begin(); p != rateElements.end(); p++)
		ss << *p << " ";
	ss << '\n';
}

void RateGroup::updateRelevantNodesinTre(Tree *t)  {
	
	for (set<int>::const_iterator p=rateElements.begin(); p != rateElements.end(); p++)
		t->updateToRootClsTis(*p);
}

void RateGroup::coutElementsInGroup()  {
	
	for (set<int>::const_iterator p=rateElements.begin(); p != rateElements.end(); p++){
		cout << (*p) << " ";
	}
}

void RateGroup::setRateValuesForEachElem(Tree *t)  {
	
	for (set<int>::const_iterator p=rateElements.begin(); p != rateElements.end(); p++){
		Node *b = t->getNodeByIndex((*p));
		b->setRtGrpVal(rate);
		b->setRtGrpIdx(indx);
	}
}

NodeRate::NodeRate(MbRandom *rp, Model *mp, int nn, double a, double b, 
				   double cp, double fxrt, int rm) : Parameter(rp, mp) {
	alpha             = a;
	beta              = b;
	concentrationParm = cp;
	numNodes = nn;
	name = "DP";
	numAuxTabs = 5;
	
	subRateModelType = rm;
	
	int ntax = (nn + 1) / 2;
	rootID = ntax;
	
	// initialize from prior
	if(subRateModelType == 2){
		double fixrate = fxrt;
		if(fixrate < 0.0){
			fixrate = ranPtr->gammaRv(alpha, beta);
			cout << "Strict clock, initial substitution rate = " << fixrate << endl;
		}
		else
			cout << "Strict clock, rates fixed to = " << fixrate << endl;
		RateGroup *rg = new RateGroup(fixrate);
		rateGroups.push_back( rg );
		for (int i=0; i<numNodes; i++){
			if(i != rootID)
				rg->addRateElement(i);
		}
	}
	else if(subRateModelType == 3){
		for (int i=0; i<numNodes; i++){
			if(i != rootID){  
				int nodesin = 0;
				if(i < rootID) nodesin = i;
				else nodesin = i - 1;
				RateGroup *rg = new RateGroup(ranPtr->gammaRv(alpha, beta));
				rateGroups.push_back( rg );
				rg->addRateElement(i);
				
			}
		}
	}
	else{
		for (int i=0; i<numNodes; i++){
			if(i != rootID){  
				int nodesin = 0;
				if(i < rootID) nodesin = i;
				else nodesin = i - 1;
				double probNewTable = concentrationParm / (nodesin + concentrationParm);
				if ( ranPtr->uniformRv() < probNewTable ){
					RateGroup *rg = new RateGroup(ranPtr->gammaRv(alpha, beta));
					rateGroups.push_back( rg );
					rg->addRateElement(i);
				}
				else{
					double u = ranPtr->uniformRv();
					double sum = 0.0;
					for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
						sum += (double)((*p)->getNumRateElements()) / nodesin;
						if ( u < sum ){
							(*p)->addRateElement(i);
							break;
						}
					}
				}
			}
		}
	}
	labelTables();
	print(std::cout);
}

NodeRate::~NodeRate(void) {

	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++)
		delete (*p);
}

NodeRate& NodeRate::operator=(const NodeRate &r) {

	if (this != &r)
		clone(r);
	return *this;
}

void NodeRate::clone(const NodeRate &r) {

}

double NodeRate::update(double &oldLnL) {
	
	if(subRateModelType == 2)
		return updateCLOCK(oldLnL);
	else if(subRateModelType == 3)
		return updateUnCorrGamma(oldLnL);
	else
		return updateDPM(oldLnL);
}

double NodeRate::updateDPM(double &oldLnL) {

	Tree *t = modelPtr->getActiveTree();
	double oldLike = oldLnL;
	
	t->upDateAllCls();
	t->upDateAllTis();

	const double tuning = log(2.0);
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		double oldR = (*p)->getRate();
		double newR = oldR * exp(tuning*(ranPtr->uniformRv()-0.5));
		(*p)->setRate(newR);
		(*p)->updateRelevantNodesinTre(t);
		modelPtr->setTiProb();
		double newLnL = modelPtr->lnLikelihood();
		double lnR = (newLnL-oldLike) + 
		             (ranPtr->lnGammaPdf(alpha, beta, newR)-ranPtr->lnGammaPdf(alpha, beta, oldR)) + 
					 (log(newR)-log(oldR));
		double r = modelPtr->safeExponentiation(lnR);
		if ( ranPtr->uniformRv() < r ){
			oldLike = newLnL;
		}
		else{
			(*p)->setRate(oldR);
			
			(*p)->updateRelevantNodesinTre(t);
		}
	}

	const int numAuxiliary = numAuxTabs;
	const double lnConcOverNumAux = log(concentrationParm/numAuxiliary);
	RateGroup **auxiliaryRateGroups = new RateGroup*[numAuxiliary];
	
	for (int i=0; i<numNodes; i++){
		if(i != rootID) {
			RateGroup *origGroup = findRateGroupWithElementIndexed( i );
			origGroup->removeRateElement( i );
			if (origGroup->getNumRateElements() == 0)
				removeRateGroup(origGroup);
			
			vector<double> lnProb;
			lnProb.reserve(rateGroups.size() + numAuxiliary);
			for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
				const int numSeatedElements = (*p)->getNumRateElements();
				(*p)->addRateElement(i);
				t->updateToRootClsTis(i);
				modelPtr->setTiProb();
				double rglnl = modelPtr->lnLikelihood();
				lnProb.push_back( log(numSeatedElements) + rglnl );
				(*p)->removeRateElement(i);
			}

			for (int j=0; j<numAuxiliary; j++)
				auxiliaryRateGroups[j] = new RateGroup(ranPtr->gammaRv(alpha, beta));
			for (int j=0; j<numAuxiliary; j++){
				RateGroup *tempGrp = NULL;  
				tempGrp = auxiliaryRateGroups[j]; 
				rateGroups.push_back( tempGrp ); 
				auxiliaryRateGroups[j]->addRateElement(i);
				t->updateToRootClsTis(i);
				modelPtr->setTiProb(); 
				double rglnl = modelPtr->lnLikelihood();
				lnProb.push_back( lnConcOverNumAux + rglnl );
				auxiliaryRateGroups[j]->removeRateElement(i);
				removeRateGroup(tempGrp); 
			}
						
			normalizeVector(lnProb);
			unsigned whichTable = ranPtr->categoricalRv(&lnProb[0], lnProb.size());
			RateGroup *newGroup = NULL;
			if ( whichTable < rateGroups.size() ){
				newGroup = rateGroups[whichTable];
				newGroup->addRateElement(i);
			}
			else{
				newGroup = auxiliaryRateGroups[whichTable-rateGroups.size()];
				rateGroups.push_back( newGroup );
				newGroup->addRateElement(i);  
			}
				
			for (int j=0; j<numAuxiliary; j++){
				if (auxiliaryRateGroups[j] != newGroup)
					delete auxiliaryRateGroups[j];
			}
			t->updateToRootClsTis(i);
		}
	}
	
	delete [] auxiliaryRateGroups;
	t->flipAllCls();
	t->flipAllTis();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	oldLnL = modelPtr->lnLikelihood(); 
	labelTables();
	setRatesForNodes(t);
	modelPtr->setLnLGood(true);
	return 0.0;
}

double NodeRate::updateCLOCK(double &oldLnL) {
	
	Tree *t = modelPtr->getActiveTree();
	double oldLike = oldLnL;
	
	t->upDateAllCls();
	t->upDateAllTis();
	
	const double tuning = log(2.0);
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		double oldR = (*p)->getRate();
		double newR = oldR * exp(tuning*(ranPtr->uniformRv()-0.5));
		(*p)->setRate(newR);
		(*p)->updateRelevantNodesinTre(t);
		modelPtr->setTiProb();
		double newLnL = modelPtr->lnLikelihood();
		double lnR = (newLnL-oldLike) + 
		(ranPtr->lnGammaPdf(alpha, beta, newR)-ranPtr->lnGammaPdf(alpha, beta, oldR)) + 
		(log(newR)-log(oldR));
		double r = modelPtr->safeExponentiation(lnR);
		if ( ranPtr->uniformRv() < r )
			oldLike = newLnL;
		else{
			(*p)->setRate(oldR);
			(*p)->updateRelevantNodesinTre(t);
		}
	}
	t->flipAllCls();
	t->flipAllTis();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	oldLnL = modelPtr->lnLikelihood(); 
	labelTables();
	setRatesForNodes(t);
	modelPtr->setLnLGood(true);
	
	return 0.0;
}

double NodeRate::updateUnCorrGamma(double &oldLnL) {
	
	Tree *t = modelPtr->getActiveTree();
	double oldLike = oldLnL;
	
	t->upDateAllCls();
	t->upDateAllTis();
	
	const double tuning = log(2.0);
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		double oldR = (*p)->getRate();
		double newR = oldR * exp(tuning*(ranPtr->uniformRv()-0.5));
		(*p)->setRate(newR);
		(*p)->updateRelevantNodesinTre(t);
		modelPtr->setTiProb();
		double newLnL = modelPtr->lnLikelihood();
		double lnR = (newLnL-oldLike) + 
		(ranPtr->lnGammaPdf(alpha, beta, newR)-ranPtr->lnGammaPdf(alpha, beta, oldR)) + 
		(log(newR)-log(oldR));
		double r = modelPtr->safeExponentiation(lnR);
		if ( ranPtr->uniformRv() < r )
			oldLike = newLnL;
		else{
			(*p)->setRate(oldR);			
			(*p)->updateRelevantNodesinTre(t);
		}
	}

	t->flipAllCls();
	t->flipAllTis();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	oldLnL = modelPtr->lnLikelihood(); 
	labelTables();
	setRatesForNodes(t);
	modelPtr->setLnLGood(true);
	return 0.0;
}

void NodeRate::labelTables(void) {

	int n = 0;
	for (vector<RateGroup *>::const_iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		(*p)->setTableIndex(n);
		n++;
	}
}


double NodeRate::lnPrior(void) {

	return 0.0;
}

void NodeRate::print(std::ostream & o) const {

	for (vector<RateGroup *>::const_iterator p=rateGroups.begin(); p != rateGroups.end(); p++)
		(*p)->print(o);
}

int NodeRate::getTableNumForNodeIndexed(int idx) {

	int num = -1;
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		if ( (*p)->isElementPresent(idx) == true ){
			num = (*p)->getTableIndex();
			break;
		}
	}
	if(num == -1) {
		cerr << "ERROR: index " << idx << " not in rateGroups" << endl;
		exit(1);
	}
	return num;
}

double NodeRate::getRateForNodeIndexed(int idx) {
	
	double theRate = 0.0;
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		if ( (*p)->isElementPresent(idx) == true ){
			theRate = (*p)->getRate();
			break;
		}
	}
	if(theRate == 0.0) {
		cerr << "ERROR: index " << idx << " not in rateGroups" << endl;
		exit(1);
	}
	return theRate;
}

RateGroup* NodeRate::findRateGroupWithElementIndexed(int idx) {

	RateGroup *theRateGroup = NULL; 
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		if ( (*p)->isElementPresent(idx) == true ){
			theRateGroup = (*p);
			break;
		}
	}
	return theRateGroup;
}

void NodeRate::removeRateGroup(RateGroup *g) {

	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++){
		if ( (*p) == g ){
			rateGroups.erase( p );
			break;
		}
	}
}

double NodeRate::getAverageRate(){
	
	double sumRt = 0.0;
	for(int i=0; i<numNodes; i++) { 
		if(i != rootID)
			sumRt += getRateForNodeIndexed(i);
	}
	return sumRt / (numNodes - 1);
}

string NodeRate::writeParam(void){
	
	stringstream ss;
	ss << "Total Rate Groups = " << getNumRateGroups() << " (" << concentrationParm << ")\n";
	ss << "Average Rate = " << getAverageRate() << "\n";
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++)
		(*p)->print(ss);
	string outp = ss.str();
	return outp;
}

void NodeRate::setRatesForNodes(Tree *t){
	
	for (vector<RateGroup *>::iterator p=rateGroups.begin(); p != rateGroups.end(); p++)
		(*p)->setRateValuesForEachElem(t);
}





