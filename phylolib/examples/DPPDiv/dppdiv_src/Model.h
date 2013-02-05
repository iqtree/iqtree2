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

#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>
#include "MbMatrix.h"

class Calibration;
class Alignment;
class Basefreq;
class Exchangeability;
class MbRandom;
class MbTransitionMatrix;
class Node;
class NodeRate;
class Parameter;
class Shape;
class Speciation;
class Tree;
class Treescale;
class Cphyperp;
class ExpCalib;
class Model {

	enum TreeDirection 
		{
		DOWN_TREE_DIR = 0,
		LEFT_TREE_DIR = 1,
		RIGHT_TREE_DIR = 2,
		DIRTY_FLAG = 3
		};
	public:
										Model(MbRandom *rp, Alignment *ap, std::string ts, double pm, 
											  double ra, double rb, double hal, double hbe, bool ubl, 
											  bool alnm, int offmv, bool rndNo, std::string clfn, int nodpr, 
											  double bdr, double bda, double fxclkrt, bool roofix,
											  bool sfb, bool ehpc, bool dphpc, int dphpng, bool gamhp, int rmod,
											  bool fxmod); 
										~Model(void);
		double							lnLikelihood(void);
		double							getPriorMeanV(void) { return priorMeanN; }
		Basefreq*						getActiveBasefreq(void);
		Tree*							getActiveTree(void);
		Treescale*						getActiveTreeScale(void);
		Exchangeability*				getActiveExchangeability(void);
		Shape*							getActiveShape(void);
		NodeRate*						getActiveNodeRate(void);
		Speciation*						getActiveSpeciation(void);
		Cphyperp*						getActiveCphyperp(void);
		ExpCalib*						getActiveExpCalib(void);
		Parameter*						pickParmToUpdate(void);
		void							printTis(std::ostream &) const;
		void							setTiProb(void);
		void							setTiProb(Node *p, Shape *s, NodeRate *r);
		void							setNodeRateGrpIndxs(void);
		void							upDateRateMatrix(void);
		void							updateAccepted(void);
		void							updateRejected(void);
		double							safeExponentiation(double lnX);
		void							switchActiveParm(void) { (activeParm == 0 ? activeParm = 1 : activeParm = 0); }
		seedType						getStartingSeed1() { return startS1; }
		seedType						getStartingSeed2() { return startS2; }
		void							setRunUnderPrior(bool b) { runUnderPrior = b; }
		void							setLnLGood(bool b) { lnLGood = b; }
		double							getMyCurrLnL(void);
		void							setMyCurrLnl(double d) { myCurLnL = d; }
		int								getTreeTimePriorNum() { return treeTimePrior; }
		double							getRootNExpRate() { return rootNExpRate; }
		bool							getExponCalibHyperParm() { return exponCalibHyperParm; }
		bool							getExponDPMCalibHyperParm() { return exponDPMCalibHyperParm; }
		void							setUpdateProbabilities(bool initial);
		
	private:
		void							initializeConditionalLikelihoods(void);
		void							initializeTransitionProbabilityMatrices(void);
		double							readCalibFile();
		Calibration*					getRootCalibration();
		
		MbRandom						*ranPtr;
		Alignment						*alignmentPtr;
		int								numGammaCats;
		std::vector<Parameter *>		parms[2];
		double							*cls;
		double							**clPtr[2];
		MbTransitionMatrix				*tiCalculator;
		int								activeParm;
		std::vector<double>				updateProb;
		int								numParms;
		int								numPatterns;
		MbMatrix<double>				**tis[2];
		double							priorMeanN;
		seedType						startS1, startS2;
		bool							runUnderPrior;
		double							fixedClockRate;
		
		typedef double ** TransitionMat;
		typedef TransitionMat * TransitionMatArray;
		TransitionMatArray * tiMatArrays;
		TreeDirection * tiMatFlags;
		
		
		double							myCurLnL;
		bool							lnLGood;
		int								turnedOffMove;
		
		std::string						calibfilen;
		std::vector<Calibration*>		calibrs;
		bool							fixRootHeight;
		int								treeTimePrior;
		bool							zeroNodeTimeMove;
		double							rHtY;
		double							rHtO;
		double							rootNExpRate;
		bool							exponCalibHyperParm;
		bool							exponDPMCalibHyperParm;
		bool							fixSomeModParams;
		bool							cpfix;
};

#endif
