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

#ifndef MCMC_H
#define MCMC_H

#include <string>
#include <fstream>

class MbRandom;
class Model;
class Mcmc {

	public:
						Mcmc(MbRandom *rp, Model *mp, int nc, int pf, int sf, 
							 std::string ofp, bool wdf, bool modUpP);
							
	private:
		void			runChain(void);
		double			safeExponentiation(double lnX);
		void			sampleChain(int gen, std::ofstream &paraOut, std::ofstream &treeOut, 
									std::ofstream &figTOut, std::ofstream &nodeOut, double lnl);
		void			sampleRtsFChain(int gen, std::ofstream &rOut);
		void			printAllModelParams(std::ofstream &dOut);
		void			writeCalibrationTree();
		int				numCycles;
		int				printFrequency;
		int				sampleFrequency;
		MbRandom		*ranPtr;
		Model			*modelPtr;
		std::string		fileNamePref;
		bool			writeInfoFile;
		bool			printratef;
		bool			modUpdateProbs;
};

#endif
