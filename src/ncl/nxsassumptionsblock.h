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

#ifndef NCL_ASSUMPTIONSBLOCK_H
#define NCL_ASSUMPTIONSBLOCK_H

/*----------------------------------------------------------------------------------------------------------------------
|	This class handles reading and storage for the NxsReader block ASSUMPTIONS. It overrides the member functions Read 
|	and Reset, which are abstract virtual functions in the base class NxsBlock. Adding a new data member? Don't forget
|	to:
|~
|	o Describe it in the class declaration using a C-style comment.
|	o Initialize it (unless it is self-initializing) in the constructor and re-initialize it in the Reset function.
|	o Describe the initial state in the constructor documentation.
|	o Delete memory allocated to it in both the destructor and Reset function.
|	o Report it in some way in the Report function.
|~
*/
class NxsAssumptionsBlock
  : public NxsBlock
	{
	public:
							NxsAssumptionsBlock(NxsTaxaBlock *t);
		virtual				~NxsAssumptionsBlock();

		void				ReplaceTaxaBlockPtr(NxsTaxaBlock *tb);
		void				SetCallback(NxsCharactersBlock *p);

		int					GetNumCharSets();
		void				GetCharSetNames(NxsStringVector &names);
		NxsUnsignedSet		&GetCharSet(NxsString nm);
		NxsString			GetDefCharSetName();

		int					GetNumTaxSets();
		void				GetTaxSetNames(NxsStringVector &names);
		NxsUnsignedSet		&GetTaxSet(NxsString nm);
		NxsString			GetDefTaxSetName();

		int					GetNumExSets();
		void				GetExSetNames(NxsStringVector &names);
		NxsUnsignedSet		&GetExSet(NxsString nm);
		NxsString			GetDefExSetName();
		void				ApplyExSet(NxsString nm);

		virtual void		Report(std::ostream& out);
		virtual void		Reset();

	private:
		NxsTaxaBlock		*taxa;				/* pointer to the NxsTaxaBlock object */
		NxsCharactersBlock	*charBlockPtr;		/* pointer to the NxsCharactersBlock-derived object to be notified in the event of exset changes */

	protected:
		NxsUnsignedSetMap	charsets;			/* the variable storing charsets */
		NxsUnsignedSetMap	taxsets;			/* the variable storing taxsets */
		NxsUnsignedSetMap	exsets;				/* the variable storing exsets */

		NxsString			def_charset;		/* the default charset */
		NxsString			def_taxset;			/* the default taxset */
		NxsString			def_exset;			/* the default exset */

	protected:
		void				HandleCharset(NxsToken& token);
		void				HandleEndblock(NxsToken& token);
		void				HandleExset(NxsToken& token);
		void				HandleTaxset(NxsToken& token);
		virtual void		Read(NxsToken& token);
		virtual unsigned	TaxonLabelToNumber(NxsString s);
	};

typedef NxsAssumptionsBlock AssumptionsBlock;	// for backward compatibility

#endif

