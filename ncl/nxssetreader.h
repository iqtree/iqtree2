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

#ifndef NCL_NXSSETREADER_H
#define NCL_NXSSETREADER_H

/*----------------------------------------------------------------------------------------------------------------------
|	A class for reading NEXUS set objects and storing them in a set of int values. The NxsUnsignedSet `nxsset' will be 
|	cleared, and `nxsset' will be built up as the set is read, with each element in the list storing a 
|	member of the set (ranges are stored as individual elements). This class handles set descriptions of the following 
|	form:
|>
|	4-7 15 20-.\3;
|>
|	The above set includes every number from 4 to 7 (inclusive), 15 and every third number from 20 to max, where `max' 
|	would ordinarily be set to either the last character (if `settype' is `NxsSetReaderEnum::charset') or the last 
|	taxon (if `settype' is `NxsSetReaderEnum::taxset'). If `max' equaled 30, the example above would be stored as
|	follows (remember that internally the numbers are stored with offset 0, even though in the NEXUS data file the
|	numbers always start at 1.
|>
|	3, 4, 5, 6, 14, 19, 22, 25, 28
|>
|	The following example of how NxsSetReader is used comes from the NxsCharactersBlock::HandleEliminate function:
|>
|	NxsSetReader(token, ncharTotal, eliminated, *this, NxsSetReader::charset).Run();
|>
|	This reads in a set of eliminated characters from a NEXUS data file, storing the resulting set in the data member
|	`eliminated'. In this case `max' is set to `ncharTotal' (the total number of characters), and the block reference
|	is set to the NxsCharactersBlock object, which provides a 
*/
class NxsSetReader
	{
	public:
		
		enum NxsSetReaderEnum	/* For use with the variable `settype' */
			{
			generic = 1,		/* means expect a generic set (say, characters weights) */
			charset,			/* means expect a character set */
			taxset				/* means expect a taxon set */
			};

						NxsSetReader(NxsToken &t, unsigned maxValue, NxsUnsignedSet &iset, NxsBlock &nxsblk, unsigned type);

		bool			Run();

	protected:

		bool			AddRange(unsigned first, unsigned last, unsigned modulus = 0);

	private:

		unsigned		GetTokenValue();

		NxsBlock		&block;		/* reference to the block object used for looking up labels */
		NxsToken		&token;		/* reference to the token being used to parse the NEXUS data file */
		NxsUnsignedSet	&nxsset;	/* reference to the NxsUnsignedSet set being read */
		unsigned		max;		/* maximum number of elements in the set */
		unsigned		settype;	/* the type of set being read (see the NxsSetReaderEnum enumeration) */
	};

typedef NxsSetReader SetReader;

#endif
