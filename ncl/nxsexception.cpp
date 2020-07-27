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
|	Copies 's' to msg and sets line, col and pos to the current line, column and position in the file where parsing
|	stopped.
*/
NxsException::NxsException(
  const NxsString &s,	/* the message for the user */
  file_pos fp,	/* the current file position */
  long fl,		/* the current file line */
  long fc)		/* the current file column */
	{
	pos		= fp;
	line	= fl;
	col		= fc;
	msg		= s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a NxsException object with the specified message, getting file position information from the NxsToken.
*/
NxsException::NxsException(
  const NxsString &s,		/* message that describes the error */
  const NxsToken &t)		/* NxsToken that was supplied the last token (the token that caused the error) */
	{
	msg		= s; 
	pos		= t.GetFilePosition();
	line	= t.GetFileLine();
	col		= t.GetFileColumn();
  	}
