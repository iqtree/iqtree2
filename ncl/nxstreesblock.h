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

#ifndef NCL_NXSTREESBLOCK_H
#define NCL_NXSTREESBLOCK_H

/*----------------------------------------------------------------------------------------------------------------------
|	This class handles reading and storage for the NEXUS block TREES. It overrides the member functions Read and Reset,
|	which are abstract virtual functions in the base class NxsBlock. The translation table (if one is supplied) is 
|	stored in the `translateList'. The tree names are stored in `treeName' and the tree descriptions in 
|	`treeDescription'. Information about rooting of trees is stored in `rooted'. Note that no checking is done to 
|	ensure that the tree descriptions are valid. The validity of the tree descriptions could be checked after the TREES
|	block has been read (but before the next block in the file has been read) by overriding the NxsReader::ExitingBlock
|	member function, but no functionality for this is provided by the NCL. Below is a table showing the correspondence
|	between the elements of a TREES block and the variables and member functions that can be used to access each piece 
|	of information stored. 
|>
|	NEXUS command     Data members    Member functions
|	-----------------------------------------------------
|	TRANSLATE         translateList
|	
|	TREE              treeName        GetTreeName
|	                                  GetTreeDescription
|	                                  GetNumTrees
|	                                  GetNumDefaultTree
|	                                  IsDefaultTree
|	
|	                  rooted          IsRootedTree
|	-----------------------------------------------------
|>
*/
class NxsTreesBlock 
  : public NxsBlock
	{
 	public:
							NxsTreesBlock(NxsTaxaBlock *tb);
		virtual				~NxsTreesBlock();

				void		ReplaceTaxaBlockPtr(NxsTaxaBlock *tb);
				unsigned	GetNumDefaultTree();
				unsigned	GetNumTrees();
				NxsString	GetTreeName(unsigned i);
				NxsString	GetTreeDescription(unsigned i);
				NxsString	GetTranslatedTreeDescription(unsigned i);
				bool		IsDefaultTree(unsigned i);
				bool		IsRootedTree(unsigned i);
		virtual void		Report(std::ostream &out);
		virtual void		BriefReport(NxsString &s);
		virtual void		Reset();

	protected :

		NxsStringMap		translateList;		/* storage for translation table (if any) */
		NxsStringVector		treeName;			/* storage for tree names */
		NxsStringVector		treeDescription;	/* storage for tree descriptions */
		NxsBoolVector		rooted;				/* stores information about rooting for each tree */
		NxsTaxaBlock		*taxa;				/* pointer to existing NxsTaxaBlock object */
		unsigned			ntrees;				/* number of trees stored */
		unsigned			defaultTree;		/* 0-offset index of default tree specified by user, or 0 if user failed to specify a default tree using an asterisk in the NEXUS data file */

		virtual	void		Read(NxsToken &token);
		void				HandleTreeDescription(NxsToken &token, bool utree);
	};

typedef NxsTreesBlock TreesBlock;

#endif
