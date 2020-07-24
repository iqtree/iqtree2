//	Copyright (C) 1999-2003 Paul O. Lewis and Mark T. Holder
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

#ifndef NCL_NXSINDENT_H
#define NCL_NXSINDENT_H

#include <ostream> //for std::ostream

/*----------------------------------------------------------------------------------------------------------------------
|	Manipulator for use in indenting text `leftMarg' characters.
*/
class Indent
	{
	public:
					Indent(unsigned i);

		unsigned	leftMarg;	/* the amount by which to indent */
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `leftMarg' to `i'.
*/
inline Indent::Indent(
  unsigned i)	/* the amount (in characters) by which to indent */
	{
	leftMarg = i;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Output operator for the Indent manipulator.
*/
inline std::ostream &operator <<(
  std::ostream &o,		/* the ostream object */
  const Indent &i)	/* the Indent object to be sent to `o' */
	{
#if defined (HAVE_PRAGMA_UNUSED)
#	pragma unused(i)
#endif
	return o;
	}

#endif
