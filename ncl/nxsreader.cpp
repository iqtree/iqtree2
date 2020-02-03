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
|	Initializes both `blockList' and `currBlock' to NULL.
*/
NxsReader::NxsReader()
	{
	blockList	= NULL;
	currBlock	= NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
NxsReader::~NxsReader()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `newBlock' to the end of the list of NxsBlock objects growing from `blockList'. If `blockList' points to NULL,
|	this function sets `blockList' to point to `newBlock'. Calls SetNexus method of `newBlock' to inform `newBlock' of
|	the NxsReader object that now owns it. This is useful when the `newBlock' object needs to communicate with the 
|	outside world through the NxsReader object, such as when it issues progress reports as it is reading the contents
|	of its block.
*/
void NxsReader::Add(
  NxsBlock *newBlock)	/* a pointer to an existing block object */
	{
	assert(newBlock != NULL);

	newBlock->SetNexus(this);

	if (!blockList)
		blockList = newBlock;
	else
		{
		// Add new block to end of list
		//
		NxsBlock *curr;
		for (curr = blockList; curr && curr->next;)
			curr = curr->next;
		assert(curr && !curr->next);
		curr->next = newBlock;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns position (first block has position 0) of block `b' in `blockList'. Returns UINT_MAX if `b' cannot be found
|	in `blockList'.
*/
unsigned NxsReader::PositionInBlockList(
  NxsBlock *b)	/* a pointer to an existing block object */
	{
	unsigned pos = 0;
	NxsBlock *curr = blockList;

	for (;;)
		{
		if (curr == NULL || curr == b)
			break;
		pos++;
		curr = curr->next;
		}

	if (curr == NULL)
		pos = UINT_MAX;

	return pos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reassign should be called if a block (`oldb') is about to be deleted (perhaps to make way for new data). Create 
|	the new block (`newb') before deleting `oldb', then call Reassign to replace `oldb' in `blockList' with `newb'. 
|	Assumes `oldb' exists and is in `blockList'.
*/
void NxsReader::Reassign(
  NxsBlock *oldb,	/* a pointer to the block object soon to be deleted */
  NxsBlock *newb)	/* a pointer to oldb's replacement */
	{
	NxsBlock *prev = NULL;
	NxsBlock *curr = blockList;
	newb->SetNexus(this);

	for (;;)
		{
		if (curr == NULL || curr == oldb)
			break;
		prev = curr;
		curr = curr->next;
		}

	assert(curr != NULL);

	newb->next = curr->next;
	if (prev == NULL) 
		blockList = newb;
	else 
		prev->next = newb;
	curr->next = NULL;
	curr->SetNexus(NULL);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `blockList' data member still equals NULL, returns true; otherwise, returns false. `blockList' will not be equal
|	to NULL if the Add function has been called to add a block object to the list.
*/
bool NxsReader::BlockListEmpty()
	{
	return (blockList == NULL ? true : false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function was created for purposes of debugging a new NxsBlock. This version does nothing; to create an active
|	DebugReportBlock function, override this version in the derived class and call the Report function of `nexusBlock'.
|	This function is called whenever the main NxsReader Execute function encounters the [&spillall] command comment 
|	between blocks in the data file. The Execute function goes through all blocks and passes them, in turn, to this 
|	DebugReportBlock function so that their contents are displayed. Placing the [&spillall] command comment between
|	different versions of a block allows multiple blocks of the same type to be tested using one long data file. Say 
|	you are interested in testing whether the normal, transpose, and interleave format of a matrix can all be read 
|	correctly. If you put three versions of the block in the data file one after the other, the second one will wipe out
|	the first, and the third one will wipe out the second, unless you have a way to report on each one before the next 
|	one is read. This function provides that ability.
*/
void NxsReader::DebugReportBlock(
  NxsBlock &nexusBlock)	/* the block that should be reported */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(nexusBlock)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Detaches `oldBlock' from the list of NxsBlock objects growing from `blockList'. If `blockList' itself points to 
|	`oldBlock', this function sets `blockList' to point to `oldBlock->next'. Note: the object pointed to by `oldBlock' 
|	is not deleted, it is simply detached from the linked list. No harm is done in Detaching a block pointer that has 
|	already been detached previously; if `oldBlock' is not found in the block list, Detach simply returns quietly. If 
|	`oldBlock' is found, its SetNexus object is called to set the NxsReader pointer to NULL, indicating that it is no 
|	longer owned by (i.e., attached to) a NxsReader object.
*/
void NxsReader::Detach(
  NxsBlock *oldBlock)	/* a pointer to an existing block object */
	{
	assert(oldBlock != NULL);

	// Return quietly if there are not blocks attached
	//
	if (blockList == NULL)
		return;

	if (blockList == oldBlock) 
		{
		blockList = oldBlock->next;
		oldBlock->SetNexus(NULL);
		}
	else 
		{
		// Bug fix MTH 6/17/2002: old version detached intervening blocks as well
		//
		NxsBlock *curr = blockList;
		for (; curr->next != NULL && curr->next != oldBlock;)
			curr = curr->next;

		// Line below can be uncommented to find cases where Detach function is 
		// called for pointers that are not in the linked list. If line below is
		// uncommented, the part of the descriptive comment that precedes this
		// function about "...simply returns quietly" will be incorrect (at least
		// in the Debugging version of the program where asserts are active).
		//
		//assert(curr->next == oldBlock);

		if (curr->next == oldBlock) 
			{
			curr->next = oldBlock->next;
			oldBlock->SetNexus(NULL);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the NxsReader object when a block named `blockName' is entered. Allows derived class overriding this
|	function to notify user of progress in parsing the NEXUS file. Also gives program the opportunity to ask user if it
|	is ok to purge data currently contained in this block. If user is asked whether existing data should be deleted, and
|	the answer comes back no, then then the overrided function should return false, otherwise it should return true.
|	This (base class) version always returns true.
*/
bool NxsReader::EnteringBlock(
  NxsString blockName)	/* the name of the block just entered */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the NxsReader object when a block named `blockName' is being exited. Allows derived class overriding this
|	function to notify user of progress in parsing the NEXUS file.
*/
void NxsReader::ExitingBlock(
  NxsString blockName)	/* the name of the block being exited */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads the NxsReader data file from the input stream provided by `token'. This function is responsible for reading 
|	through the name of a each block. Once it has read a block name, it searches `blockList' for a block object to 
|	handle reading the remainder of the block's contents. The block object is responsible for reading the END or 
|	ENDBLOCK command as well as the trailing semicolon. This function also handles reading comments that are outside 
|	of blocks, as well as the initial "#NEXUS" keyword. The `notifyStartStop' argument is provided in case you do not 
|	wish the ExecuteStart and ExecuteStop functions to be called. These functions are primarily used for creating and 
|	destroying a dialog box to show progress, and nested Execute calls can thus cause problems (e.g., a dialog box is 
|	destroyed when the inner Execute calls ExecuteStop and the outer Execute still expects the dialog box to be 
|	available). Specifying `notifyStartStop' false for all the nested Execute calls thus allows the outermost Execute 
|	call to control creation and destruction of the dialog box.
*/
void NxsReader::Execute(
  NxsToken	&token,				/* the token object used to grab NxsReader tokens */
  bool		notifyStartStop)	/* if true, ExecuteStarting and ExecuteStopping will be called */
	{
	char id_str[256];
	currBlock = NULL;

	bool disabledBlock = false;
	NxsString errormsg;

	try
		{
		token.GetNextToken();
		}
	catch (NxsException x)
		{
		NexusError(token.errormsg, 0, 0, 0);
		return;
		}

	if (!token.Equals("#NEXUS"))
		{
		errormsg = "Expecting #NEXUS to be the first token in the file, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		return;
		}

	if (notifyStartStop)
		ExecuteStarting();

	for (;;)
		{
		token.SetLabileFlagBit(NxsToken::saveCommandComments);
		token.GetNextToken();

		if (token.AtEOF())
			break;

		if (token.Equals("BEGIN"))
			{
			disabledBlock = false;
			token.GetNextToken();

			for (currBlock = blockList; currBlock != NULL; currBlock = currBlock->next)
				{
				if (token.Equals(currBlock->GetID()))
					{
					if (currBlock->IsEnabled()) 
						{
						strcpy(id_str, currBlock->GetID().c_str());
						bool ok_to_read = EnteringBlock(id_str);
						if (!ok_to_read) 
							currBlock = NULL;
						else
							{
							currBlock->Reset();

							// We need to back up currBlock, because the Read statement might trigger
							// a recursive call to Execute (if the block contains instructions to execute 
							// another file, then the same NxsReader object may be used and any member fields (e.g. currBlock)
							//  could be trashed.
							//
							NxsBlock *tempBlock = currBlock;	

							try 
								{
								currBlock->Read(token);
								currBlock = tempBlock;
								}

							catch (NxsException x) 
								{
								currBlock = tempBlock;
								if (currBlock->errormsg.length() > 0)
									NexusError(currBlock->errormsg, x.pos, x.line, x.col);
								else
									NexusError(x.msg, x.pos, x.line, x.col);
								currBlock = NULL;
								return;
								}	// catch (NxsException x) 
							ExitingBlock(id_str /*currBlock->GetID()*/);
							}	// else
						}	// if (currBlock->IsEnabled()) 

					else
						{
						disabledBlock = true;
						SkippingDisabledBlock(token.GetToken());
						}
					break;
					}	// if (token.Equals(currBlock->GetID()))
				}	// for (currBlock = blockList; currBlock != NULL; currBlock = currBlock->next)

			if (currBlock == NULL)
				{
				token.BlanksToUnderscores();
				NxsString currBlockName = token.GetToken();

				if (!disabledBlock) 
					SkippingBlock(currBlockName);

				for (;;)
					{
                    token.SetLabileFlagBit(token.hyphenNotPunctuation);
					token.GetNextToken();

					if (token.Equals("END") || token.Equals("ENDBLOCK")) 
						{
						token.GetNextToken();

						if (!token.Equals(";")) 
							{
							errormsg = "Expecting ';' after END or ENDBLOCK command, but found ";
							errormsg += token.GetToken();
							errormsg += " instead";
							NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
							return;
							}
						break;
						}

					if (token.AtEOF()) 
						{
						errormsg = "Encountered end of file before END or ENDBLOCK in block ";
						errormsg += currBlockName;
						NexusError(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						return;
						}
					}	// for (;;)
				}	// if (currBlock == NULL)
			currBlock = NULL;
			}	// if (token.Equals("BEGIN"))

		else if (token.Equals("&SHOWALL"))
			{
			for (NxsBlock*  showBlock = blockList; showBlock != NULL; showBlock = showBlock->next)
				{
				DebugReportBlock(*showBlock);
				}
			}

		else if (token.Equals("&LEAVE"))
			{
			break;
			}

		} // for (;;)

	if (notifyStartStop)
		ExecuteStopping();

	currBlock = NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing the copyright notice for the NxsReader Class Library, useful for reporting the use of 
|	this library by programs that interact with the user.
*/
const char *NxsReader::NCLCopyrightNotice()
	{
	return NCL_COPYRIGHT;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing the URL for the NxsReader Class Library internet home page.
*/
const char *NxsReader::NCLHomePageURL()
	{
	return NCL_HOMEPAGEURL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing the name and current version of the NxsReader Class Library, useful for reporting the 
|	use of this library by programs that interact with the user.
*/
const char *NxsReader::NCLNameAndVersion()
	{
	return NCL_NAME_AND_VERSION;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called just after Execute member function reads the opening "#NEXUS" token in a NEXUS data file. Override this 
|	virtual base class function if your application needs to do anything at this point in the execution of a NEXUS data
|	file (e.g. good opportunity to pop up a dialog box showing progress). Be sure to call the Execute function with the
|	`notifyStartStop' argument set to true, otherwise ExecuteStarting will not be called.
|	
*/
void NxsReader::ExecuteStarting()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when Execute member function encounters the end of the NEXUS data file, or the special comment [&LEAVE] is
|	found between NEXUS blocks. Override this virtual base class function if your application needs to do anything at 
|	this point in the execution of a NEXUS data file (e.g. good opportunity to hide or destroy a dialog box showing 
|	progress). Be sure to call the Execute function with the `notifyStartStop' argument set to true, otherwise 
|	ExecuteStopping will not be called.
*/
void NxsReader::ExecuteStopping()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when an error is encountered in a NEXUS file. Allows program to give user details of the error as well as 
|	the precise location of the error.
*/
void NxsReader::NexusError(
  NxsString	msg,	/* the error message to be displayed */
  file_pos	pos,	/* the current file position */
  long	line,	/* the current file line */
  long	col)	/* the current column within the current file line */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(msg, pos, line, col)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function may be used to report progess while reading through a file. For example, the NxsAllelesBlock class 
|	uses this function to report the name of the population it is currently reading so the user doesn't think the 
|	program has hung on large data sets.
*/
void NxsReader::OutputComment(
  const NxsString &comment)	/* a comment to be shown on the output */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(comment)
#	endif
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when an unknown block named `blockName' is about to be skipped. Override this pure virtual
|	function to provide an indication of progress as the NEXUS file is being read.
*/
void NxsReader::SkippingBlock(
  NxsString blockName)	/* the name of the block being skipped */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when a disabled block named `blockName' is encountered in a NEXUS data file being executed.
|	Override this pure virtual function to handle this event in an appropriate manner. For example, the program may 
|	wish to inform the user that a data block was encountered in what is supposed to be a tree file.
*/
void NxsReader::SkippingDisabledBlock(
  NxsString blockName)	/* the name of the disabled block being skipped */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(blockName)
#	endif
	}

