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

#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `states' to NULL.
*/
NxsDiscreteDatum::NxsDiscreteDatum()
	{
	states = NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes memory associated with `states' (if any was allocated).
*/
NxsDiscreteDatum::~NxsDiscreteDatum()
	{
	if (states != NULL)
		delete [] states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes this NxsDiscreteDatum object an exact copy of `other'. Useful for dealing with matchchar symbols in a matrix.
*/
void NxsDiscreteDatum::CopyFrom(
  const NxsDiscreteDatum &other)	/* the source NxsDiscreteDatum object  */
	{
	if (states != NULL)
		{
		delete [] states;
		states = NULL;
		}

	if (other.states == NULL)
		return;

	unsigned sz = other.states[0];
	if (sz == 0)
		{
		// First element of other.states is zero, indicating that the gap state is present
		//
		states = new unsigned[1];
		states[0] = 0;
		}

	else if (sz == 1)
		{
		// First element of other.states is one, indicating that there is just one state for this
		// taxon-character combination. With just one state, no need to worry about either 
		// ambiguity or polymorphism.
		//
		states = new unsigned[2];
		states[0] = 1;
		states[1] = other.states[1];
		}

	else
		{
		// First element of other.states is greater than 1, indicating that ambiguity or 
		// polymorphism is present.
		//
		states = new unsigned[sz + 2];
		states[0] = sz;
		for (unsigned i = 1; i <= sz; i++)
			states[i] = other.states[i];

		// Copy the polymorphism indicator element.
		//
		states[sz + 1] = other.states[sz + 1];
		}
	}
