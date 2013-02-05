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

#ifndef PARAMETER_EXPCALIB_H
#define PARAMETER_EXPCALIB_H

#include <vector>
#include <set>
#include <string>
#include <ostream>

class MbRandom;
class Model;
class Node;
class Tree;
class LambdaTable {
	
	public:
							LambdaTable(MbRandom *rp, double r) :randP(rp), lambda(r), indx(-1) {}
		void				addDiner(int idx) { diners.insert(idx); }
		void				removeDiner(int idx) { diners.erase(idx); }
		double				getTableLambda(void) const { return lambda; }
		bool				getIsDinerSeated(int idx) const { return diners.find(idx) != diners.end(); }
		void				setTableLambda(double x) { lambda = x; }
		int					getNumDiners(void) const { return (int)diners.size(); }
		void				print(std::ostream &ss) const;
		std::set<int>		getDinerList() { return diners; }
		
		int					getTableIndx(void) { return indx; }
		void				setTableIndx(int i) { indx = i; }
		
		void				updateNodesAtTable(Tree *t);
		double				getPriorPrForDiners(Tree *t);
		void				printLambdaTableInfo();
				
	private:
		double				lambda;
		std::set<int>		diners;
		int					indx;
		MbRandom			*randP;
};


class ExpCalib : public Parameter {

	public:
										ExpCalib(MbRandom *rp, Model *mp, bool dphplc, int dphpng,
												 double ts, bool ghp);
										~ExpCalib(void);
		ExpCalib						&operator=(const ExpCalib &c);
		void							clone(const ExpCalib &c);
		double							update(double &oldLnL);
		void							print(std::ostream & o) const;
		double							lnPrior(void){ return 0.0; }
		std::string						writeParam();
		
		double							getEpsilonValue() { return epsilonValue; }
		double							getCurMajorityLambda() { return curMajorityLambda; }
		double							getCurOutlieLambda() { return curOutlieLambda; }
		double							getLambdaForNode();
		
		void							getAllExpHPCalibratedNodes();
		bool							getIsLambdaContaminationClass(double l);
		double							getDPMExpHPConcentParam() { return dpmCP; }
		int								getNumLambdaTables() { return dpmLambdaHyp.size(); }
	
	private:
		double							epsilonValue;
		double							betaAlph1, betaAlph2;
		double							curMajorityLambda, curOutlieLambda;
		double							majorityExpParm, outlieExpParm;
		std::vector<Node *>				calibNodeList;
		std::vector<double>				nodeDeltas;
		std::vector<LambdaTable *>		dpmLambdaHyp;
		bool							dpmLHP;
		double							dpmCP;
		double							dpmLambdaExpParm;
		int								prNGrpsDPM;
		bool							gammaHPVals;
		double							gHPAlph, gHPBetaM, gHPBetaO, gHPBetaDP;
		
		void							updateEpsilonValue();
		void							updateMajorityLambda();
		void							updateOutlierLambda();
		void							updateNodeTaintedClassAssignment(double tScl);
		void							updateContamination();
		void							updateDPMHyperPrior();
		void							initializeDPMLambdasToCals();
		void							assignDPMLambdasToCals();
		LambdaTable*					findLambTabWCalID(int ix);
		void							removeDPMLambdaTable(LambdaTable *g);
		double							getProbsAcrossAllTables();
		
	
};	


#endif