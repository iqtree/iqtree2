//	Copyright (C) 1999-2002 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library).
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
|	Sets the base class data member `id' to the name of the block (i.e. "EMPTY") in NEXUS data files.
*/
NxsEmptyBlock::NxsEmptyBlock()
	{
	id = "EMPTY";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing needs to be done.
*/
NxsEmptyBlock::~NxsEmptyBlock()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The code here is identical to the base class version (simply returns 0), so the code here should either be modified
|	or this derived version eliminated altogether. Under what circumstances would you need to modify the default code, 
|	you ask? This function should be modified to something meaningful if this derived class needs to construct and run
|	a NxsSetReader object to read a set involving characters. The NxsSetReader object may need to use this function to
|	look up a character label encountered in the set. A class that overrides this method should return the character 
|	index in the range [1..`nchar']; i.e., add one to the 0-offset index.
*/
unsigned NxsEmptyBlock::CharLabelToNumber(
  NxsString s)	/* the character label to be translated to character number */
	{
	return 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when the END or ENDBLOCK command needs to be parsed from within the EMPTY block. Basically just checks to 
|	make sure the next token in the data file is a semicolon.
*/
void NxsEmptyBlock::HandleEndblock(
  NxsToken &token)	/* the token used to read from `in' */
	{
	// Get the semicolon following END or ENDBLOCK token
	//
	token.GetNextToken();

	if(!token.Equals(";"))
		{
		errormsg = "Expecting ';' to terminate the END or ENDBLOCK command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the ability to read everything following the block name (which is read by the NxsReader 
|	object) to the END or ENDBLOCK statement. Characters are read from the input stream `in'. Overrides the pure 
|	virtual function in the base class.
*/
void NxsEmptyBlock::Read(
  NxsToken &token)	/* the token used to read from `in'*/
	{
	isEmpty = false;

	// This should be the semicolon after the block name
	//
	token.GetNextToken();

	if (!token.Equals(";"))
		{
		errormsg = "Expecting ';' after ";
		errormsg += id;
		errormsg += " block name, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	for(;;)
		{
		token.GetNextToken();

		if (token.Equals("END"))
			{
			HandleEndblock(token);
			break;
			}

		else if(token.Equals("ENDBLOCK"))
			{
			HandleEndblock(token);
			break;
			}

		else
			{
			SkippingCommand(token.GetToken());

			do
				{
				token.GetNextToken();
				}
			while (!token.AtEOF() && !token.Equals(";"));

			if (token.AtEOF())
				{
				errormsg = "Unexpected end of file encountered";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `isEmpty' to true in preparation for reading a new EMPTY block. Overrides the pure virtual function in the 
|	base class.
*/
void NxsEmptyBlock::Reset()
	{
	isEmpty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function outputs a brief report of the contents of this EMPTY block. Overrides the pure virtual function in 
|	the base class.
*/
void NxsEmptyBlock::Report(
  ostream &out)	/* the output stream to which to write the report */
	{
	out << endl;
	out << id << " block contains...";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when an unknown command named `commandName' is about to be skipped. This version of the 
|	function (which is identical to the base class version) does nothing (i.e., no warning is issued that a command 
|	was unrecognized). Modify this virtual function to provide such warnings to the user (or eliminate it altogether 
|	since the base class version already does what this does). 
*/
void NxsEmptyBlock::SkippingCommand(
  NxsString commandName)	/* the name of the command being skipped */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The code here is identical to the base class version (simply returns 0), so the code here should either be modified
|	or this derived version eliminated altogether. Under what circumstances would you need to modify the default code,
|	you ask? This function should be modified to something meaningful if this derived class needs to construct and run 
|	a NxsSetReader object to read a set involving taxa. The NxsSetReader object may need to use this function to look 
|	up a taxon label encountered in the set. A class that overrides this method should return the taxon index in the 
|	range [1..ntax]; i.e., add one to the 0-offset index.
*/
unsigned NxsEmptyBlock::TaxonLabelToNumber(
  NxsString s)	/* the taxon label to be translated to a taxon number */
	{
	return 0;
	}

