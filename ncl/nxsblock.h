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
#ifndef NCL_NXSBLOCK_H
#define NCL_NXSBLOCK_H

#include "nxsstring.h"	//for NxsString
#include "nxstoken.h"   //for NxsToken
class NxsReader;

/*----------------------------------------------------------------------------------------------------------------------
|	This is the base class from which all block classes are derived. A NxsBlock-derived class encapsulates a Nexus block
|	(e.g. DATA block, TREES block, etc.). The abstract virtual function Read must be overridden for each derived class 
|	to provide the ability to read everything following the block name (which is read by the NxsReader object) to the 
|	end or endblock statement. Derived classes must provide their own data storage and access functions. The abstract
|	virtual function Report must be overridden to provide some feedback to user on contents of block. The abstract
|	virtual function Reset must be overridden to empty the block of all its contents, restoring it to its 
|	just-constructed state.
*/
class NxsBlock
	{
	friend class NxsReader;

	public:
							NxsBlock();
		virtual				~NxsBlock();

		void				SetNexus(NxsReader *nxsptr);

		NxsString			GetID();
		bool				IsEmpty();

		void				Enable();
		void				Disable();
		bool				IsEnabled();
		bool				IsUserSupplied();

		virtual unsigned	CharLabelToNumber(NxsString s);
		virtual unsigned	TaxonLabelToNumber(NxsString s);

		virtual void		SkippingCommand(NxsString commandName);

		virtual void		Report(std::ostream &out);
		virtual void		Reset();

		NxsString			errormsg;			/* workspace for creating error messages */

	protected:
		bool				isEmpty;			/* true if this object is currently storing data */
		bool				isEnabled;			/* true if this block is currently ebabled */
		bool				isUserSupplied;		/* true if this object has been read from a file; false otherwise */
		NxsReader			*nexus;				/* pointer to the Nexus file reader object */
		NxsBlock			*next;				/* pointer to next block in list */
		NxsString			id;					/* holds name of block (e.g., "DATA", "TREES", etc.) */
				
		virtual void		Read(NxsToken &token);
	};

#endif


