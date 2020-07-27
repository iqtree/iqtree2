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
#ifndef NCL_NXSTAXABLOCK_H
#define NCL_NXSTAXABLOCK_H

#include <iostream>    //for std::ostream
#include "nxsstring.h"
#include "nxsdefs.h"   //for NxsBoolVector
#include "nxstoken.h"  //for NxsToken
#include "nxsblock.h"  //for NxsBlock
/*----------------------------------------------------------------------------------------------------------------------
|	This class handles reading and storage for the NxsReader block TAXA. It overrides the member functions Read and 
|	Reset, which are abstract virtual functions in the base class NxsBlock. The taxon names are stored in an vector of
|	strings (taxonLabels) that is accessible through the member functions GetTaxonLabel(int), AddTaxonLabel(NxsString), 
|	ChangeTaxonLabel(int, NxsString), and GetNumTaxonLabels().
*/
class NxsTaxaBlock
  : public NxsBlock
	{
	friend class NxsDataBlock;
	friend class NxsAllelesBlock;
	friend class NxsCharactersBlock;
	friend class NxsDistancesBlock;

	public:
							NxsTaxaBlock();
		virtual				~NxsTaxaBlock();

		virtual unsigned	AddTaxonLabel(NxsString s);
		void  				ChangeTaxonLabel(unsigned i, NxsString s);
		unsigned			FindTaxon(NxsString label);
		bool  				IsAlreadyDefined(NxsString label);
		unsigned			GetMaxTaxonLabelLength();
		unsigned			GetNumTaxonLabels();
		NxsString 			GetTaxonLabel(unsigned i);
		bool 				NeedsQuotes(unsigned i);
		virtual void		Report(std::ostream &out);
		virtual void 		Reset();

		class NxsX_NoSuchTaxon {};	/* thrown if FindTaxon cannot locate a supplied taxon label in the taxonLabels vector */

	protected:
		unsigned		ntax;			/* number of taxa */
		NxsStringVector	taxonLabels;	/* storage for list of taxon labels */
		NxsBoolVector 	needsQuotes;	/* needsQuotes[i] true if label i needs to be quoted when output */

		virtual void 	Read(NxsToken &token);

	private:
		void 			SetNtax(unsigned n);
	};

// The following typedef maintains compatibility with existing code.
// The TaxaBlock class name is deprecated; please use NxsTaxaBlock instead.
//
typedef NxsTaxaBlock TaxaBlock;

#endif
