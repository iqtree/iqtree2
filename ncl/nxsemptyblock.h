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
#ifndef NCL_NXSEMPTYBLOCK_H
#define NCL_NXSEMPTYBLOCK_H

/*----------------------------------------------------------------------------------------------------------------------
|	This is a template that can be used to create a class representing a NEXUS block. Here are the steps to follow if
|	you wish to create a new block specifically for use with your particular application. Suppose your application is
|	called Phylome and you want to create a private block called a PHYLOME block that can appear in NEXUS data files
|	and contains commands for your program.
|~
|	o Copy the files nxsemptyblock.h and nxsemptyblock.cpp and rename them (e.g. nxsphylomeblock.h and 
|	  nxsphylomeblock.cpp)
|	o In nxsphylomeblock.h and nxsphylomeblock.cpp, replace all instances of EMPTY (case-sensitive, whole word search)
|	  with PHYLOME
|	o In nxsphylomeblock.h, replace both instances of NCL_NXSEMPTYBLOCK_H at the top of the file with
|	  NCL_NXSPHYLOMEBLOCK_H
|	o In nxsphylomeblock.h and nxsphylomeblock.cpp, replace all instances of NxsEmptyBlock (case-sensitive, whole word
|	  search) with NxsPhylomeBlock
|	o Modify the Read function in nxsphylomeblock.cpp to interpret what comes after the BEGIN PHYLOME command in the
|	  NEXUS data file
|	o Modify the CharLabelToNumber and TaxonLabelToNumber if you need to read in sets of characters or taxa, 
|	  respectively. These functions provide a way for NxsSetReader objects to translate character or taxon labels to
|	  the corresponding numbers. If you do not need these capabilities, then it is safe to just delete these functions
|	  from nxsphylomeblock.h and nxsphylomeblock.cpp because they are no different that the base class versions
|	o Modify the SkippingCommand function if you want to notify users when commands within the PHYLOME block are not 
|	  recognized and are being skipped
|	o In nxsphylomeblock.h, replace this comment with something meaningful for your class. Start off with something
|	  like "This class handles reading and storage for the NEXUS block PHYLOME. It overrides the member functions 
|	  Read and Reset, which are abstract virtual functions in the base class NxsBlock"
|~
|	Adding a new data member? Don't forget to:
|~
|	o Describe it in the class declaration using a C-style comment. 
|	o Initialize it (unless it is self-initializing) in the constructor and reinitialize it in the Reset function.
|	o Describe the initial state in the constructor documentation. 
|	o Delete memory allocated to it in both the destructor and Reset function. 
|	o Report it in some way in the Report function. 
|~
*/
class NxsEmptyBlock
  : public NxsBlock
	{
	public:

						NxsEmptyBlock();
		virtual			~NxsEmptyBlock();

		virtual void	Report(ostream &out);

	protected:

		void			SkippingCommand(NxsString commandName);
		unsigned		TaxonLabelToNumber(NxsString s);
		unsigned		CharLabelToNumber(NxsString s);
		void			HandleEndblock(NxsToken &token);
		virtual void	Read(NxsToken &token);
		virtual void	Reset();
	};

#endif
