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

#ifndef PARAMETER_TREESCALE_H
#define PARAMETER_TREESCALE_H

class MbRandom;
class Model;
class ExpCalib;
class Tree;
class Treescale : public Parameter {
	
	public:
							Treescale(MbRandom *rp, Model *mp, double sv, double yb, double ob, int dt, bool calib, bool exhpc);
							~Treescale(void);
		Treescale			&operator=(const Treescale &c);
		void				clone(const Treescale &c);
		double				update(double &oldLnL);
		void				print(std::ostream & o) const;
		double				lnPrior(void);
		double				lnExponentialTSPriorRatio(double newTS, double oldTS);
		std::string			writeParam(void);
		double				getScaleValue() { return scaleVal; }
		void				setScaleValue(double s) { scaleVal = s; }
		double				getLnTreeProb(Tree *t);
		double				calcLnPriorProbPlusJaco(Tree *t, double otp);
		
	private:
		double				scaleVal;
		double				oldBound;
		double				yngBound;
		double				tuning;
		bool				isBounded;
		int					treeTimePrior;
		bool				retune;
		int					tsPriorD;
		bool				isUserCalib;
		double				expoRate;
		bool				exponHPCalib;
};

#endif