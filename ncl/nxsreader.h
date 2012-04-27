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

#ifndef NCL_NXSREADER_H
#define NCL_NXSREADER_H

/*----------------------------------------------------------------------------------------------------------------------
|	This is the class that orchestrates the reading of a NEXUS data file. An object of this class should be created, 
|	and objects of any block classes that are expected to be needed should be added to `blockList' using the Add 
|	member function. The Execute member function is then called, which reads the data file until encountering a block 
|	name, at which point the correct block is looked up in `blockList' and that object's Read method called. 
*/
class NxsReader
	{
	public:
		enum	NxsTolerateFlags	/* Flags used with data member tolerate used to allow some flexibility with respect to the NEXUS format */
			{
			allowMissingInEquate	= 0x0001,	/* if set, equate symbols are allowed for missing data symbol */
			allowPunctuationInNames	= 0x0002	/* if set, some punctuation is allowed within tokens representing labels for taxa, characters, and sets */
			};

						NxsReader();
		virtual			~NxsReader();

		bool			BlockListEmpty();
		unsigned		PositionInBlockList(NxsBlock *b);
		void			Add(NxsBlock *newBlock);
		void			Detach(NxsBlock *newBlock);
		void			Reassign(NxsBlock *oldb, NxsBlock *newb);
		void			Execute(NxsToken& token, bool notifyStartStop = true);

		virtual void	DebugReportBlock(NxsBlock &nexusBlock);

		const char			*NCLNameAndVersion();
		const char			*NCLCopyrightNotice();
		const char			*NCLHomePageURL();

		virtual void	ExecuteStarting();
		virtual void	ExecuteStopping();

		virtual bool	EnteringBlock(NxsString blockName);
		virtual void	ExitingBlock(NxsString blockName);

		virtual void	OutputComment(const NxsString &comment);

		virtual void	NexusError(NxsString msg, file_pos pos, long line, long col);

		virtual void	SkippingDisabledBlock(NxsString blockName);
		virtual void	SkippingBlock(NxsString blockName);

	protected:

		NxsBlock		*blockList;	/* pointer to first block in list of blocks */
		NxsBlock		*currBlock;	/* pointer to current block in list of blocks */
	};

typedef NxsBlock NexusBlock;
typedef NxsReader Nexus;

#endif

