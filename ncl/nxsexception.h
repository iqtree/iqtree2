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

#ifndef NCL_NXSEXCEPTION_H
#define NCL_NXSEXCEPTION_H

#include "nxsstring.h" //for NxsString
#include "nxsdefs.h"   //for file_pos

class NxsToken;

/*----------------------------------------------------------------------------------------------------------------------
|	Exception class that conveys a message specific to the problem encountered.
*/
class NxsException
	{
	public:
		NxsString	msg;	/* NxsString to hold message */
		file_pos	pos;	/* current file position */
		long		line;	/* current line in file */
		long		col;	/* column of current line */

		explicit NxsException(const NxsString &s, file_pos fp = 0, long fl = 0L, long fc = 0L);
		NxsException(const NxsString &s, const NxsToken &t);
	};

typedef NxsException XNexus;

#endif
