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

#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>

class MbRandom;
class Model;
class Parameter {

	public:
								Parameter(MbRandom *rp, Model *mp); 
		virtual					~Parameter(){}

		Parameter				&operator=(Parameter &p);
		std::string				getName(void) { return name; }
		void					setName(std::string s) { name = s; }
		virtual double			update(double &oldLnL)=0;
		virtual double			lnPrior(void)=0;
		virtual void			print(std::ostream &) const = 0;
		virtual std::string		writeParam(void)=0;
						
	protected:
		std::string				name;
		MbRandom				*ranPtr;
		Model					*modelPtr;
};

#endif
