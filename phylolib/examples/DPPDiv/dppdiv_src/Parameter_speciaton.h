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

#ifndef PARAMETER_SPECIATION_H
#define PARAMETER_SPECIATION_H

class MbRandom;
class Model;
class Tree;
class Speciation : public Parameter {
	
	public:
							Speciation(MbRandom *rp, Model *mp, double bdr, double bda);
							~Speciation(void);
		Speciation			&operator=(const Speciation &c);
		void				clone(const Speciation &c);
		double				update(double &oldLnL);
		void				print(std::ostream & o) const;
		double				lnPrior(void);
		std::string			writeParam(void);
		double				getRelativeDeath() { return relativeDeath; }
		void				setRelativeDeath(double v) { relativeDeath = v; }
		double				getNetDiversification() { return netDiversificaton; }
		void				setNetDiversification(double v) { netDiversificaton = v; }
		double				getLnTreeProb(Tree *t);
		
	private:
		double				relativeDeath;		// mu / lambda
		double				netDiversificaton; // lambda - mu
		int					treeTimePrior;
		double				maxdivV;
		
		double				updateRelDeathRt();
		double				updateNetDivRate();
};


#endif