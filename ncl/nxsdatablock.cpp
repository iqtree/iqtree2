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
|	Sets `id' to "DATA" and `newtaxa' to true, and calls the base class (NxsCharactersBlock) constructor.
*/
NxsDataBlock::NxsDataBlock(
  NxsTaxaBlock *tb,			/* the taxa block object for storing taxon labels */
  NxsAssumptionsBlock *ab)	/* the assumptions block object for storing exsets */
  : NxsCharactersBlock(tb, ab)
	{
	id = "DATA";
	newtaxa = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Reset function of the parent class (NxsCharactersBlock) and resets `newtaxa' to true in preparation for 
|	reading another DATA block.
*/
void NxsDataBlock::Reset()
	{
	NxsCharactersBlock::Reset();
	newtaxa = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts this NxsDataBlock object into a NxsCharactersBlock object, storing the result in the supplied 
|	NxsCharactersBlock object. This NxsDataBlock object will subsequently say it is empty when asked, and 
|	the `formerly_datablock' data member of `charactersblock' will be set to true.
*/
void NxsDataBlock::TransferTo(
  NxsCharactersBlock &charactersblock)	/* the NxsCharactersBlock object that will receive all the data from this object */
	{
	charactersblock.Reset();
	charactersblock.Consume((NxsCharactersBlock &)(*this));
	}
