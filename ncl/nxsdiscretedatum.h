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

#ifndef NCL_NXSDISCRETEDATUM_H
#define NCL_NXSDISCRETEDATUM_H

/*----------------------------------------------------------------------------------------------------------------------
|	Class for holding discrete states in a matrix. Note that there is no way to access the variables of this class 
|	since they are all private and there are no public access functions. This class is designed to be manipulated by 
|	the class NxsDiscreteMatrix, which is the only class that has been designated a friend of NxsDiscreteDatum. 
|	The variable `states' is NULL if there is missing data, and non-NULL for any other state. If `states' is non-NULL,
|	the first cell is used to store the number of states. This will be 0 if the state is the gap state, 1 if the state
|	is unambiguous and nonpolymorphic (and not the gap state of course), and 2 or higher if there is either 
|	polymorphism or uncertainty. If polymorphism or uncertainty apply, it becomes necessary to store information about 
|	which of these two situations holds. Thus, the last cell in the array is set to either 1 (polymorphism) or 0 
|	(uncertainty). While a little complicated, this scheme has the following benefits:
|~
|	o if the state is missing, the only memory allocated is for a pointer (`states')
|	o if the state is unambiguous and not polymorphic, no storage is used for keeping track of whether polymorphism or 
|	  uncertainty holds
|	o it allows for a virtually unlimited number of states, which is important if it is to be general enough to store 
|	  microsatellite data for a NxsAllelesBlock object, for example.
|~
|	Supposing the gap symbol is '-', the missing data symbol is '?', and the symbols list is "ACGT", the following 
|	table shows the status of the states variable under several different possible data matrix entries:
|>
|	Matrix entry        states array
|	--------------------------------
|	     ?              NULL
|	     -              [0]
|	     G              [1][2]
|	(AG) polymorphic    [2][0][2][1]
|	{AG} ambiguous      [2][0][2][0]
|	--------------------------------
|>	
*/
class NxsDiscreteDatum
	{
	friend class NxsDiscreteMatrix;

	public:

					NxsDiscreteDatum();
		virtual		~NxsDiscreteDatum();

		void		CopyFrom(const NxsDiscreteDatum &other);

	private:

		unsigned	*states;	/* holds information about state for a single taxon-character combination */
	};

typedef NxsDiscreteDatum DiscreteDatum;

#endif
