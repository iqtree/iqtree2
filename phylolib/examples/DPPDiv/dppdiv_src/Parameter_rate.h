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

#ifndef PARAMETER_RATE_H
#define PARAMETER_RATE_H

#include <vector>
#include <set>
#include <string>
#include <ostream>
#include "MbVector.h"

class MbRandom;
class Tree;

/*----------------------------------------------
 class RateGroup
 ------------------------------------------------*/
class RateGroup {

	public:
							RateGroup(double r) :rate(r), indx(-1) {}
		void				addRateElement(int idx) { rateElements.insert(idx); }
		void				removeRateElement(int idx) { rateElements.erase(idx); }
		double				getRate(void) const { return rate; }
		bool				isElementPresent(int idx) const { return rateElements.find(idx) != rateElements.end(); }
		void				setRate(double x) { rate = x; }
		int					getNumRateElements(void) const { return rateElements.size(); }
		void				print(std::ostream &ss) const;
		
		int					getTableIndex(void) { return indx; }
		void				setTableIndex(int i) { indx = i; }
		
		void				updateRelevantNodesinTre(Tree *t);
		void				coutElementsInGroup();
		void				setRateValuesForEachElem(Tree *t);
							
	private:
		double				rate;
		std::set<int>		rateElements;
		int					indx;
};

class Model;

/*----------------------------------------------
 class NodeRate
 ------------------------------------------------*/
class NodeRate : public Parameter {

	public:
									NodeRate(MbRandom *rp, Model *mp, int nn, double a, double b, 
											 double cp, double fxrt, int rm);
									~NodeRate(void); 
		NodeRate					&operator=(const NodeRate &r);
		void						clone(const NodeRate &r);
		RateGroup*					findRateGroupWithElementIndexed(int idx);
		double						update(double &oldLnL);
		double						updateDPM(double &oldLnL);
		double						updateCLOCK(double &oldLnL);
		double						updateUnCorrGamma(double &oldLnL);
		double						lnPrior(void);
		void						print(std::ostream &) const;
		int							getTableNumForNodeIndexed(int idx);
		double						getRateForNodeIndexed(int idx);
		void						removeRateGroup(RateGroup *g);
		double						getAverageRate(void);
		std::string					writeParam(void);
		int							getNumRateGroups(void) { return rateGroups.size(); }
		double						getConcenParam(){ return concentrationParm; }
		void						setConcenParam(double a) { concentrationParm = a; }
		int							getNumberofNodes(){ return numNodes; }
		void						setNumAuxTabs(int n) { numAuxTabs = n; }
												
	private:
		void						labelTables(void);
		void						setRatesForNodes(Tree *t);

		double						alpha;
		double						beta;
		int							numNodes;
		double						concentrationParm;
		std::vector<RateGroup *>	rateGroups;
		int							rootID;
		int							numAuxTabs;
		int							subRateModelType;

};

#endif
