//	Copyright (C) 1999-2003 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library) version 2.0.
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef NCL_NXSDISTANCEDATUM_H
#define NCL_NXSDISTANCEDATUM_H

/*----------------------------------------------------------------------------------------------------------------------
|	This class stores pairwise distance values. It has no public access functions, reflecting the fact that it is 
|	manipulated strictly by its only friend class, the NxsDistancesBlock class.
*/
class NxsDistanceDatum
	{
	friend class NxsDistancesBlock;

	public:

					NxsDistanceDatum();
		virtual		~NxsDistanceDatum();

	private:

		double		value;		/* the pairwise distance value stored */
		bool		missing;	/* true if there is missing data for this pair */
	};

typedef NxsDistanceDatum DistanceDatum;

#endif
