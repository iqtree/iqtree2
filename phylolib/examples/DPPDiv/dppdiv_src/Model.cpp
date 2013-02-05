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

#include "Alignment.h"
#include "MbRandom.h"
#include "MbTransitionMatrix.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_basefreq.h"
#include "Parameter_exchangeability.h"
#include "Parameter_expcalib.h"
#include "Parameter_rate.h"
#include "Parameter_shape.h"
#include "Parameter_speciaton.h"
#include "Parameter_tree.h"
#include "Parameter_cphyperp.h"
#include "Parameter_treescale.h"
#include "Calibration.h"
#include "util.h"
#include <string>
#include <vector>
#include <fstream>

using namespace std;

Model::Model(MbRandom *rp, Alignment *ap, string ts, double pm, double ra, double rb, 
			 double hal, double hbe, bool ubl, bool alnm, int offmv, bool rndNo, 
			 string clfn, int nodpr, double bdr, double bda, double fxclkrt, bool roofix,
			 bool sfb, bool ehpc, bool dphpc, int dphpng, bool gamhp, int rmod, bool fxmod) {

	ranPtr       = rp;
	alignmentPtr = ap;
	priorMeanN   = pm;
	ranPtr->getSeed(startS1, startS2);
	runUnderPrior = false;
	treeTimePrior = nodpr;
	myCurLnL = 0.0;
	lnLGood = false;
	calibfilen = clfn;
	fixRootHeight = roofix;
	zeroNodeTimeMove = false;
	exponCalibHyperParm = ehpc;
	exponDPMCalibHyperParm = dphpc;
	turnedOffMove = offmv;
	fixedClockRate = fxclkrt;
	fixSomeModParams = fxmod;
	double initRootH = 1.0;
	rHtY = 0.0;
	rHtO = 0.0;
	int tsPrDist = 1;
	rootNExpRate = -1.0;
	bool rtCalib = false;
	if(calibfilen.empty() == false){
		initRootH = readCalibFile();
		Calibration *rCal = getRootCalibration();
		if(rCal != NULL && fixRootHeight == false){
			rHtY = rCal->getYngTime();
			rHtO = rCal->getOldTime();
			tsPrDist = rCal->getPriorDistributionType();
			rtCalib = true;
			rootNExpRate =  rCal->getCalExponRate();
		}
		else if(rCal == NULL && fixRootHeight == false){
			rHtY = rHtY / 1.1;
			rHtO = -1.0;
			tsPrDist = 2;
		}
	}
	
	cout << "\nStarting with seeds: { " << startS1 << " , " << startS2 << " } \n\n";
	
	numGammaCats = 4;
	numPatterns  = alignmentPtr->getNumChar();
	
	cpfix = false;
	if(turnedOffMove == 5)
		cpfix = true;
	else if(turnedOffMove == 6)
		cpfix = true;
	if(rmod > 1)
		cpfix = true;

	int nn = 2*alignmentPtr->getNumTaxa()-1;
	if(pm > nn - 1){
		cerr << "ERROR: the prior on the mean number of tables cannot exceed the number of nodes in the tree!" << endl;
		exit(1);
	}
	Cphyperp *conp = new Cphyperp(ranPtr, this, hal, hbe, nn, pm, cpfix);
	ExpCalib *excal = new ExpCalib(ranPtr, this, dphpc, dphpng, initRootH, gamhp);
	NodeRate *nr = new NodeRate(ranPtr, this, nn, ra, rb, conp->getCurrentCP(), fxclkrt, rmod);
	for (int i=0; i<2; i++){ 
		parms[i].push_back( new Basefreq(ranPtr, this, 4, fxmod) );					// base frequency parameter
		parms[i].push_back( new Exchangeability(ranPtr, this) );				// rate parameters of the GTR model
		parms[i].push_back( new Shape(ranPtr, this, numGammaCats, 2.0, fxmod) );		// gamma shape parameter for rate variation across sites
		parms[i].push_back( new Tree(ranPtr, this, alignmentPtr, ts, ubl, alnm, rndNo, calibrs, initRootH, sfb, ehpc, excal) );    // rooted phylogenetic tree
		parms[i].push_back( nr );												// restaurant containing node rates
		parms[i].push_back( conp );												// hyper prior on DPP concentration parameter
		parms[i].push_back( new Treescale(ranPtr, this, initRootH, rHtY, rHtO, tsPrDist, rtCalib, ehpc) ); // the tree scale prior
		parms[i].push_back( new Speciation(ranPtr, this, bdr, bda) );			// hyper prior on diversification for cBDP speciation
		parms[i].push_back( excal );											// hyper prior exponential node calibration parameters
	}
	numParms = parms[0].size();
	activeParm = 0;
	for (int i=0; i<numParms; i++)
		*parms[0][i] = *parms[1][i];

	for (int i=0; i<numParms; i++)
		parms[0][i]->print(std::cout);
	
		
	setUpdateProbabilities(true);	
	if(ehpc)
		excal->getAllExpHPCalibratedNodes();
	initializeConditionalLikelihoods();
	initializeTransitionProbabilityMatrices();
	tiCalculator = new MbTransitionMatrix( getActiveExchangeability()->getRate(), getActiveBasefreq()->getFreq(), true );
	
	setTiProb();
	myCurLnL = lnLikelihood();
	cout << "lnL = " << myCurLnL << endl;
}

Model::~Model(void) {

	delete [] cls;
	delete tiCalculator;
	for (int i=0; i<2; i++){
		delete [] clPtr[i];
		delete [] tis[i][0];
		delete [] tis[i];
	}
}

Basefreq* Model::getActiveBasefreq(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Basefreq *derivedPtr = dynamic_cast<Basefreq *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Tree* Model::getActiveTree(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Tree *derivedPtr = dynamic_cast<Tree *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Treescale* Model::getActiveTreeScale(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Treescale *derivedPtr = dynamic_cast<Treescale *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}


Exchangeability* Model::getActiveExchangeability(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Exchangeability *derivedPtr = dynamic_cast<Exchangeability *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Shape* Model::getActiveShape(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Shape *derivedPtr = dynamic_cast<Shape *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

NodeRate* Model::getActiveNodeRate(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		NodeRate *derivedPtr = dynamic_cast<NodeRate *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Speciation* Model::getActiveSpeciation(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Speciation *derivedPtr = dynamic_cast<Speciation *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Cphyperp* Model::getActiveCphyperp(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Cphyperp *derivedPtr = dynamic_cast<Cphyperp *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

ExpCalib* Model::getActiveExpCalib(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		ExpCalib *derivedPtr = dynamic_cast<ExpCalib *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}


void Model::initializeConditionalLikelihoods(void) {

	int nNodes = 2*alignmentPtr->getNumTaxa()-1;
	int nChar  = alignmentPtr->getNumChar();
	int sizeOneNode = nChar * numGammaCats * 4;
	int sizeOneSpace = nNodes * sizeOneNode;
	cls = new double[2 * sizeOneSpace];
	for (int i=0; i<2*sizeOneSpace; i++)
		cls[i] = 0.0;
	for (int i=0; i<2; i++)
		{
		clPtr[i] = new double*[nNodes];
		for (int j=0; j<nNodes; j++)
			clPtr[i][j] = &cls[ i * sizeOneSpace + j * sizeOneNode ];
		}
		
	for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
		{
		double *cl0 = clPtr[0][i];
		double *cl1 = clPtr[1][i];
		for (int j=0; j<alignmentPtr->getNumChar(); j++)
			{
			int nucCode = alignmentPtr->getNucleotide(i, j);
			int possibleNucs[4];
			alignmentPtr->getPossibleNucs(nucCode, possibleNucs);
			for (int k=0; k<numGammaCats; k++)
				{
				for (int s=0; s<4; s++)
					{
					cl0[s] = (double)possibleNucs[s];
					cl1[s] = (double)possibleNucs[s];
					}
				cl0 += 4;
				cl1 += 4;
				}
			}
		}
}

void Model::initializeTransitionProbabilityMatrices(void) {

	int nNodes = 2*alignmentPtr->getNumTaxa()-1;
	for (int i=0; i<2; i++)
		{
		tis[i] = new MbMatrix<double>*[nNodes];
		tis[i][0] = new MbMatrix<double>[numGammaCats*nNodes];
		for (int j=1; j<nNodes; j++)
			tis[i][j] = tis[i][j-1] + numGammaCats;
		}
	for (int i=0; i<2; i++)
		for (int j=0; j<nNodes; j++)
			for (int k=0; k<numGammaCats; k++)
				tis[i][j][k] = MbMatrix<double>(4,4);
}

double Model::lnLikelihood(void) {

	if(runUnderPrior){
		myCurLnL = 0.0;
		return 0.0;
	}
	Tree *t = getActiveTree();
	MbMatrix<double> *tL = new MbMatrix<double>[numGammaCats];
	MbMatrix<double> *tR = new MbMatrix<double>[numGammaCats];
	
	for (int n=0; n<t->getNumNodes(); n++) {
		Node *p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getIsClDirty() == true) {
			double *clL = clPtr[p->getLft()->getActiveCl()][p->getLft()->getIdx()];
			double *clR = clPtr[p->getRht()->getActiveCl()][p->getRht()->getIdx()];
			double *clP = clPtr[p->getActiveCl()          ][p->getIdx()          ];
			for (int k=0; k<numGammaCats; k++) {
				tL[k] = tis[p->getLft()->getActiveTi()][p->getLft()->getIdx()][k];
				tR[k] = tis[p->getRht()->getActiveTi()][p->getRht()->getIdx()][k];
			}
			for (int c=0; c<numPatterns; c++) {
				for (int k=0; k<numGammaCats; k++) {
#					if 0
					for (int i=0; i<4; i++) {
						double sumL = 0.0, sumR = 0.0;
						for (int j=0; j<4; j++) {
							sumL += tL[k][i][j] * clL[j];
							sumR += tR[k][i][j] * clR[j];
						}
						clP[i] = sumL * sumR;
					}
#					else
					double sumL = 0.0, sumR = 0.0;
					sumL += tL[k][0][0] * clL[0] + tL[k][0][1] * clL[1] + tL[k][0][2] * clL[2] +tL[k][0][3] * clL[3];
					sumR += tR[k][0][0] * clR[0] + tR[k][0][1] * clR[1] + tR[k][0][2] * clR[2] +tR[k][0][3] * clR[3];
					clP[0] = sumL * sumR;
					
					sumL = 0.0, sumR = 0.0;
					sumL += tL[k][1][0] * clL[0] + tL[k][1][1] * clL[1] + tL[k][1][2] * clL[2] +tL[k][1][3] * clL[3];
					sumR += tR[k][1][0] * clR[0] + tR[k][1][1] * clR[1] + tR[k][1][2] * clR[2] +tR[k][1][3] * clR[3];
					clP[1] = sumL * sumR;
					
					sumL = 0.0, sumR = 0.0;
					sumL += tL[k][2][0] * clL[0] + tL[k][2][1] * clL[1] + tL[k][2][2] * clL[2] +tL[k][2][3] * clL[3];
					sumR += tR[k][2][0] * clR[0] + tR[k][2][1] * clR[1] + tR[k][2][2] * clR[2] +tR[k][2][3] * clR[3];
					clP[2] = sumL * sumR;
					
					sumL = 0.0, sumR = 0.0;
					sumL += tL[k][3][0] * clL[0] + tL[k][3][1] * clL[1] + tL[k][3][2] * clL[2] +tL[k][3][3] * clL[3];
					sumR += tR[k][3][0] * clR[0] + tR[k][3][1] * clR[1] + tR[k][3][2] * clR[2] +tR[k][3][3] * clR[3];
					clP[3] = sumL * sumR;
#					endif
					
					clL += 4; 
					clR += 4;
					clP += 4;
				}
			}
			p->setIsClDirty(false);
		}
	}
		
	Node *r = t->getRoot();
	MbVector<double> f = getActiveBasefreq()->getFreq();
	double *clP = clPtr[r->getActiveCl()][r->getIdx()];
	double catProb = 1.0 / numGammaCats;
	double lnL = 0.0;
	for (int c=0; c<numPatterns; c++){
		double siteProb = 0.0;
#		if 0
		for (int k=0; k<4; k++){
			for (int i=0; i<4; i++)
				siteProb += clP[i] * f[i] ;
			clP += 4;
		}
#		else
		siteProb += clP[0] * f[0] + clP[1] * f[1] + clP[2] * f[2] + clP[3] * f[3] ;
		clP += 4;
		siteProb += clP[0] * f[0] + clP[1] * f[1] + clP[2] * f[2] + clP[3] * f[3] ;
		clP += 4;
		siteProb += clP[0] * f[0] + clP[1] * f[1] + clP[2] * f[2] + clP[3] * f[3] ;
		clP += 4;
		siteProb += clP[0] * f[0] + clP[1] * f[1] + clP[2] * f[2] + clP[3] * f[3] ;
		clP += 4;
#		endif
		siteProb *= catProb;
		lnL += alignmentPtr->getNumSitesOfPattern(c) * log(siteProb);
	}

	delete [] tL;
	delete [] tR;
	myCurLnL = lnL;
	return lnL;
}

Parameter* Model::pickParmToUpdate(void) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	Parameter *parm = NULL;
	for (unsigned i=0; i<updateProb.size(); i++)
		{
		sum += updateProb[i];
		if ( u < sum )
			{
			parm = parms[activeParm][i];
			break;
			}
		}
	return parm;
}

void Model::printTis(std::ostream & o) const {

	int nNodes = 2*alignmentPtr->getNumTaxa()-1;
	for (int j=0; j<nNodes; j++)
		{
		o << "Node " << j << endl;
		for (int a=0; a<4; a++)
			{
			for (int i=0; i<2; i++)
				{
				for (int k=0; k<numGammaCats; k++)
					{
					for (int b=0; b<4; b++)
						{
						o << fixed << setprecision(10) << tis[i][j][k][a][b] << " ";
						}
					}
				}
			o << '\n';
			}
		}
	o.flush();	
}

double Model::safeExponentiation(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}

void Model::setTiProb(void) {

	Tree *t     = getActiveTree();
	Shape *s    = getActiveShape();
	NodeRate *r = getActiveNodeRate();
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node *p = t->getDownPassNode(n);
		if (p->getAnc() != NULL && p->getIsTiDirty() == true) 
			{
			setTiProb(p, s, r);
			p->setIsTiDirty(false);
			}
		}
}

void Model::setTiProb(Node *p, Shape *s, NodeRate *r) {

	int activeTi = p->getActiveTi();
	int idx      = p->getIdx();
	double branchProportion =  p->getAnc()->getNodeDepth() - p->getNodeDepth();
	double rP = r->getRateForNodeIndexed(p->getIdx());
	double v = branchProportion * rP;
	
	if(rP == 0.0){
		cerr << "ERROR: Problem rP = 0" << endl;
		exit(1);
	}
	
#	if 0
	for (int k=0; k<4; k++){
		double rt = s->getRate(k);
		tis[activeTi][idx][k] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][k] );
	}
#	else
	double rt = s->getRate(0);
	tis[activeTi][idx][0] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][0] );
	
	rt = s->getRate(1);
	tis[activeTi][idx][1] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][1] );
	
	rt = s->getRate(2);
	tis[activeTi][idx][2] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][2] );
	
	rt = s->getRate(3);
	tis[activeTi][idx][3] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][3] );
#	endif
}

void Model::setNodeRateGrpIndxs(void) {
	
	Tree *t     = getActiveTree();
	NodeRate *r = getActiveNodeRate();
	int rtID = t->getRoot()->getIdx();  
	for(int n=0; n<t->getNumNodes(); n++){
		if(n != rtID){
			Node *p = t->getNodeByIndex(n);
			int tn = r->getTableNumForNodeIndexed(n);
			double srt = r->getRateForNodeIndexed(n);
			p->setRtGrpIdx(tn);
			p->setRtGrpVal(srt);
		}
	}
}

void Model::updateAccepted(void) {

	int from = activeParm, to;
	if (from == 0)
		to = 1;
	else
		to = 0;
	for (int i=0; i<numParms; i++)
		*parms[to][i] = *parms[from][i];
}

void Model::updateRejected(void) {

	int to = activeParm, from;
	if (to == 0)
		from = 1;
	else
		from = 0;
	for (int i=0; i<numParms; i++)
		*parms[to][i] = *parms[from][i];

}

void Model::upDateRateMatrix(void) {

	tiCalculator->updateQ( getActiveExchangeability()->getRate(), getActiveBasefreq()->getFreq() );
}

double Model::getMyCurrLnL(void) {
	
	if(lnLGood){
		lnLGood = false;
		return myCurLnL;
	}
	else
		return lnLikelihood();
}


double Model::readCalibFile(void) {
		
	cout << "\nCalibrations:" << endl;
	bool rootIs = false;
	Calibration *rooCal;
	string ln = getLineFromFile(calibfilen, 1);
	int nlins = atoi(ln.c_str());
	int nnodes = alignmentPtr->getNumTaxa() - 1;
	string *calList = new string[nlins];
	for(int i=0; i<nlins; i++){
		calList[i] = getLineFromFile(calibfilen, i+2);
		Calibration *cal = new Calibration(calList[i]);
		calibrs.push_back(cal);
		if(cal->getIsRootCalib()){
			rooCal = cal;
			rootIs = true;
		}
	}
	delete [] calList;

	double initTScale = 1.0;
	double yb = 0.0;
	double ob = 0.0;
	if(rootIs){
		if(rooCal->getPriorDistributionType() == 1){
			yb = rooCal->getYngTime();
			ob = rooCal->getOldTime();
			if(yb == ob){
				initTScale = yb;
				fixRootHeight = true;
			}
			else{
				initTScale = yb + (ranPtr->uniformRv() * (ob - yb));
				fixRootHeight = false;
			}
		}
		else if(rooCal->getPriorDistributionType() == 2){
			fixRootHeight = false;
			yb = rooCal->getYngTime();
			double expMean = yb * 0.2;
			initTScale = yb + ranPtr->exponentialRv(1 / expMean);
		}
	}
	else{
		yb = 0.0;
		fixRootHeight = false;
		for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
			double tmpv;
			if((*v)->getPriorDistributionType() == 1)
				tmpv = (*v)->getOldTime();
			else if((*v)->getPriorDistributionType() == 2)
				tmpv = (*v)->getYngTime() * 1.1;
			if(tmpv > yb)
				yb = tmpv;
		}
		ob = yb + (yb * 2);
		double tsc = yb + (ranPtr->uniformRv() * (ob - yb));
		initTScale = tsc;
		rHtY = yb;
		rHtO = ob;
	}
	
	if(nlins == nnodes){
		bool fixall = true;
		for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
			double tmpo = (*v)->getOldTime();
			double tmpy = (*v)->getYngTime();
			if(tmpo != tmpy){
				fixall = false;
				break;
			}
		}
		zeroNodeTimeMove = fixall;
	}
	else if(nlins == nnodes-1 && rootIs == false){
		bool fixall = true;
		for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
			double tmpo = (*v)->getOldTime();
			double tmpy = (*v)->getYngTime();
			if(tmpo != tmpy){
				fixall = false;
				break;
			}
		}
		zeroNodeTimeMove = fixall;
	}
	
	cout << "\nInitial root height : " << initTScale <<  " [" << yb << ", " << ob << "]" << endl;
	return initTScale;
}

Calibration* Model::getRootCalibration(void) {

	for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
		if((*v)->getIsRootCalib())
			return (*v);
	}
	return NULL;
}

void Model::setUpdateProbabilities(bool initial) {

	double bfp, srp, shp, ntp, dpp, cpa, tsp, spp, ehp;
	if(initial){
		bfp = 0.3;
		srp = 0.3;
		shp = 0.3;
		ntp = 0.4;
		dpp = 0.5;
		cpa = 0.3;
		tsp = 0.3;
		spp = 0.4;
		ehp = 0.0;
	}
	else{
		bfp = 0.2;
		srp = 0.2;
		shp = 0.2;
		ntp = 0.4;
		dpp = 0.5;
		cpa = 0.3;
		tsp = 0.4;
		spp = 0.4;
		ehp = 0.0;
	}
	if(turnedOffMove == 1)
		bfp = 0.0;
	else if(turnedOffMove == 2)
		srp = 0.0;
	else if(turnedOffMove == 3)
		shp = 0.0;
	else if(turnedOffMove == 4){
		ntp = 0.0;
		tsp = 0.0;
		spp = 0.0;
		ehp = 0.0;
	}
	else if(turnedOffMove == 5){
		dpp = 0.0;
		cpa = 0.0;
		setNodeRateGrpIndxs(); 
		*parms[0][3] = *parms[1][3]; 
	}
	else if(turnedOffMove == 6 || cpfix == true) 
		cpa = 0.0;
	else if(turnedOffMove == 8)
		spp = 0.0;
	if(fixRootHeight)
		tsp = 0.0;
	if(treeTimePrior == 1 || treeTimePrior == 4)
		spp = 0.0;
	if(zeroNodeTimeMove == 1){
		ntp = 0.0;
		cout << "All internal node times are fixed" << endl;
	}
	if(fixedClockRate < 0.0)
		cpa = 0.0;
	if(exponCalibHyperParm){
		if(initial)
			ehp = 0.3;
		else
			ehp = 0.4;
	}
	if(fixSomeModParams){
		bfp = 0.0;
		shp = 0.0;
	}
	updateProb.clear();
	updateProb.push_back(bfp); // 1 basefreq
	updateProb.push_back(srp); // 2 sub rates
	updateProb.push_back(shp); // 3 gamma shape
	updateProb.push_back(ntp); // 4 node times
	updateProb.push_back(dpp); // 5 dpp rates
	updateProb.push_back(cpa); // 6 concentration parameter
	updateProb.push_back(tsp); // 7 tree scale parameter
	updateProb.push_back(spp); // 8 speciation parameters
	updateProb.push_back(ehp); // 9 exponential calibration hyper priors
	double sum = 0.0;
	for (unsigned i=0; i<updateProb.size(); i++)
		sum += updateProb[i];
	for (unsigned i=0; i<updateProb.size(); i++)
		updateProb[i] /= sum;
}

