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

#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <string>


class Calibration {
	private:
		std::string		txn1, txn2;
		double			youngtime, oldtime;
		int				nodeIDX;
		bool			isRootCal;
		int				prDistType; // 1 = Uniform, 2 = offset Exp
		double			exponRate, exponMean;
		
	public:
						Calibration(std::string calstr);
		std::string		getTxN1() { return txn1; }
		std::string		getTxN2() { return txn2; }
		double			getYngTime() { return youngtime; }
		double			getOldTime() { return oldtime; }
		void			setNodeIndex(int i) { nodeIDX = i; }
		int				getNodeIndex() { return nodeIDX; }
		void			setIsRootCalib(bool b) { isRootCal = b; }
		bool			getIsRootCalib() { return isRootCal; }
		int				getPriorDistributionType() { return prDistType; }
		double			getCalExponRate() { return exponRate; }
		double			getCalExponMean() { return exponMean; }
};


#endif
