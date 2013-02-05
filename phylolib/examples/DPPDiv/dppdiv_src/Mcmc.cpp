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
#include "Mcmc.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_basefreq.h"
#include "Parameter_exchangeability.h"
#include "Parameter_expcalib.h"
#include "Parameter_rate.h"
#include "Parameter_tree.h"
#include "Parameter_shape.h"
#include "Parameter_speciaton.h"
#include "util.h"
#include <iomanip>
#include <iostream>
#include <ctime>

using namespace std;

Mcmc::Mcmc(MbRandom *rp, Model *mp, int nc, int pf, int sf, string ofp, bool wdf, bool modUpP) {

	ranPtr          = rp;
	modelPtr        = mp;
	numCycles       = nc;
	printFrequency  = pf;
	sampleFrequency = sf;
	fileNamePref    = ofp;
	writeInfoFile   = wdf;
	printratef		= false;
	modUpdateProbs  = modUpP;
	runChain();
}

void Mcmc::runChain(void) {
	
	string pFile = fileNamePref + ".p"; // parameter file name
	string treFile = fileNamePref + ".t"; // tree file name
	string figTFile = fileNamePref + ".ant.tre"; // write to a file with the nodes colored by their rate classes
	string dFile = fileNamePref + ".info.out";
	string ndFile = fileNamePref + ".nodes.out"; // info about nodes
	string mtxFile = fileNamePref + ".rates.out"; // info about nodes
	ofstream pOut(pFile.c_str(), ios::out);
	ofstream tOut(treFile.c_str(), ios::out);
	ofstream fTOut(figTFile.c_str(), ios::out);
	ofstream nOut(ndFile.c_str(), ios::out);
	ofstream mxOut;
	ofstream dOut;
	
	writeCalibrationTree();
	if(writeInfoFile)
		dOut.open(dFile.c_str(), ios::out);
	if(printratef)
		mxOut.open(mtxFile.c_str(), ios::out);
		
	double oldLnLikelihood = modelPtr->lnLikelihood();
	
	// verbose logging
	if(writeInfoFile){
		dOut << "Running MCMC with:\n";
		dOut << "   Starting seeds = { " << modelPtr->getStartingSeed1() << " , " << modelPtr->getStartingSeed2() << " } \n";
		dOut << "   # Gens = " << numCycles << "\n";
		dOut << "   Prior mean # groups = " << modelPtr->getPriorMeanV() << "\n";
		dOut << "   lnL = " << oldLnLikelihood << "\n";
		printAllModelParams(dOut);
	}
	
	int modifyUProbsGen = (int)numCycles * 0.5;
	for (int n=1; n<=numCycles; n++){
		if(modUpdateProbs && n == modifyUProbsGen)
			modelPtr->setUpdateProbabilities(false);

		modelPtr->switchActiveParm();
		Parameter *parm = modelPtr->pickParmToUpdate();
		
		double prevlnl = oldLnLikelihood;
		double lnPriorProposalRatio = parm->update(oldLnLikelihood);
		
		double newLnLikelihood = modelPtr->getMyCurrLnL(); 
		double lnLikelihoodRatio = newLnLikelihood - oldLnLikelihood;
		
		double lnR = lnLikelihoodRatio + lnPriorProposalRatio;
		double r = safeExponentiation(lnR);
		
		bool isAccepted = false;
		if ( ranPtr->uniformRv() < r )
			isAccepted = true;
		
		if ( n % printFrequency == 0 || n == 1){
			cout << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
			if(writeInfoFile){
				dOut << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
				dOut << n << " -- " << parm->writeParam();
			}
		}
		
		if (isAccepted == true){
			oldLnLikelihood = newLnLikelihood;
			modelPtr->updateAccepted();
		}
		else{
			modelPtr->updateRejected();
			Tree *t = modelPtr->getActiveTree(); 
			t->flipAllCls();
			t->flipAllTis();
			t->upDateAllCls();
			t->upDateAllTis();
			modelPtr->upDateRateMatrix();
			modelPtr->setTiProb();
		}
		
		if(n < 10){ 
			Tree *t = modelPtr->getActiveTree(); 
			t->setNodeRateValues();
		}
		
		if ( n % sampleFrequency == 0 || n == 1)
			sampleChain(n, pOut, tOut, fTOut, nOut, oldLnLikelihood);
			
	}
	cout << "   Markov chain completed." << endl;
	pOut.close();
	tOut.close();
	fTOut.close();
	dOut.close();
	nOut.close();
	mxOut.close();
}

double Mcmc::safeExponentiation(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}

void Mcmc::sampleChain(int gen, ofstream &paraOut, ofstream &treeOut, ofstream &figTOut, 
					   ofstream &nodeOut, double lnl) {

	Basefreq *f = modelPtr->getActiveBasefreq();
	Exchangeability *e = modelPtr->getActiveExchangeability();
	NodeRate *nr = modelPtr->getActiveNodeRate();
	Tree *t = modelPtr->getActiveTree();
	Shape *sh = modelPtr->getActiveShape();
	Speciation *sp = modelPtr->getActiveSpeciation();
	ExpCalib *hpex;
	bool expHPCal = modelPtr->getExponCalibHyperParm();
	bool dpmHPCal = modelPtr->getExponDPMCalibHyperParm();
	if(expHPCal)
		hpex = modelPtr->getActiveExpCalib();
	
	if(gen == 1){
		paraOut << "Gen\tlnLikelihood\tf(A)\tf(C)\tf(G)\tf(T)";
		paraOut << "\tr(AC)\tr(AG)\tr(AT)\tr(CG)\tr(CT)\tr(GT)\tshape\tave rate\tnum rate groups\tconc param\n";
		treeOut << "#NEXUS\nbegin trees;\n";
		figTOut << "#NEXUS\nbegin trees;\n";
		nodeOut << "Gen\tlnL";
		nodeOut << "\tNetDiv(b-d)\tRelativeDeath(d/b)\tPr(speciation)\tave.subrate\tnum.DPMgroups\tDPM.conc";
		if(expHPCal){
			if(dpmHPCal)
				nodeOut << "\texpHP.dpmConP\texpHP.dpmNumLs";
			else
				nodeOut << "\texpHP.lambda1\texpHP.lambda2\texpHP.epsilon";
		}
			
		nodeOut << t->getNodeInfoNames();
		if(expHPCal)
			nodeOut << t->getCalNodeInfoNames();
		nodeOut << "\n";
	}

	paraOut << gen << "\t" << lnl;
	for(int i=0; i<f->getNumStates(); i++)
		paraOut << "\t" << f->getFreq(i);
	for(int i=0; i<6; i++)
		paraOut << "\t" << e->getRate(i);
	paraOut << "\t" << sh->getAlphaSh();
	paraOut << "\t" << nr->getAverageRate();
	paraOut << "\t" << nr->getNumRateGroups();
	paraOut << "\t" << nr->getConcenParam();
	paraOut << "\n";
	
	treeOut << "  tree t" << gen << " = ";
	treeOut << t->getTreeDescription() << "\n";
	
	figTOut << "  tree t" << gen << " = ";
	figTOut << t->getFigTreeDescription() << "\n";
	
	nodeOut << gen << "\t" << lnl;
	nodeOut << "\t" << sp->getNetDiversification();
	nodeOut << "\t" << sp->getRelativeDeath();
	nodeOut << "\t" << t->getTreeSpeciationProbability();
	nodeOut << "\t" << nr->getAverageRate();
	nodeOut << "\t" << nr->getNumRateGroups();
	nodeOut << "\t" << nr->getConcenParam();
	if(expHPCal){
		if(dpmHPCal){
			nodeOut << "\t" << hpex->getDPMExpHPConcentParam();
			nodeOut << "\t" << hpex->getNumLambdaTables();
		}
		else{
			nodeOut << "\t" << hpex->getCurMajorityLambda();
			nodeOut << "\t" << hpex->getCurOutlieLambda();
			nodeOut << "\t" << hpex->getEpsilonValue();
		}
	}
	nodeOut << t->getNodeInfoList();
	if(expHPCal)
		nodeOut << t->getCalNodeInfoList();
	nodeOut << "\n";
	
	if(gen == numCycles){
		treeOut << "end;\n";
		figTOut << "end;\n";
		figTOut << "\nbegin figtree;\n";
		figTOut << "    set appearance.branchColorAttribute=\"rate_cat\";\n";
		figTOut << "    set appearance.branchLineWidth=2.0;\n";
		figTOut << "    set scaleBar.isShown=false;\n";
		figTOut << "end;\n";
	}
}

void Mcmc::printAllModelParams(ofstream &dOut){
	
	dOut << "\n--------------------------------------------------\n";
	dOut << "Initial: \n";
	dOut << modelPtr->getActiveBasefreq()->writeParam();
	dOut << modelPtr->getActiveExchangeability()->writeParam();
	dOut << modelPtr->getActiveShape()->writeParam();
	dOut << modelPtr->getActiveNodeRate()->writeParam();
	dOut << modelPtr->getActiveTree()->writeParam();
	dOut << "--------------------------------------------------\n\n";
}

void Mcmc::writeCalibrationTree(){
	
	modelPtr->setNodeRateGrpIndxs();
	Tree *t = modelPtr->getActiveTree();
	if(t->getIsCalibratedTree()){
		string ctfn = fileNamePref + ".CALIB.tre";
		ofstream cto(ctfn.c_str(), ios::out);
		cto << "#NEXUS\nbegin trees;\n";
		cto << "    tree calib_init = ";
		cto << t->getCalibInitialTree() << "\nend;\n";
		cto << "\nbegin figtree;\n";
		cto << "    set appearance.branchColorAttribute=\"User selection\";\n";
		cto << "    set appearance.branchLineWidth=3.0;\n";
		cto << "    set nodeBars.isShown=true;\n";
		cto << "    set nodeBars.barWidth=27.0;\n";
		cto << "    set nodeLabels.isShown=true;\n";
		cto << "    set nodeLabels.displayAttribute=\"Node ages\";\n";
		cto << "    set nodeLabels.fontSize=18;\n";
		cto << "    set nodeLabels.fontStyle=2;\n";
		cto << "    set rectilinearLayout.rootLength=0;\n";
		cto << "    set scaleAxis.isShown=true;\n";
		cto << "    set scaleAxis.reverseAxis=true;\n";
		cto << "    set tipLabels.fontSize=20;\n";
		cto << "    set scaleBar.isShown=false;\n";
		cto << "end;\n";
	}
}

void Mcmc::sampleRtsFChain(int gen, std::ofstream &rOut){
	
	NodeRate *nr = modelPtr->getActiveNodeRate();
	Tree *t = modelPtr->getActiveTree();
	if(gen == 1){
		rOut << "Gen\tDPM.conc";
		rOut << t->getDownPNodeInfoNames();
		rOut << "\n";
	}
	
	rOut << gen << "\t" << nr->getConcenParam();
	rOut << t->getDownPNodeInfoList();
	rOut << "\n";
}

