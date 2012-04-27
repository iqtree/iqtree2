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

#ifndef NCL_NXSDISTANCESBLOCK_H
#define NCL_NXSDISTANCESBLOCK_H

class NxsDistanceDatum;

/*----------------------------------------------------------------------------------------------------------------------
|	This class handles reading and storage for the NEXUS block DISTANCES. It overrides the member functions Read and 
|	Reset, which are abstract virtual functions in the base class NxsBlock. Below is a table showing the correspondence 
|	between the elements of a DISTANCES block and the variables and member functions that can be used to access each 
|	piece of information stored.
|>
|	NEXUS command   Command attribute  Data Members        Member Functions
|	------------------------------------------------------------------------
|	DIMENSIONS      NEWTAXA            newtaxa
|	
|	                NTAX               ntax                GetNtax 
|	
|	                NCHAR              nchar               GetNchar
|	
|	FORMAT          TRIANGLE           triangle            GetTriangle
|	                                                       IsUpperTriangular
|	                                                       IsLowerTriangular
|	                                                       IsRectangular
|	
|	                [NO]DIAGONAL       diagonal            IsDiagonal
|	
|	                [NO]LABELS         labels              IsLabels
|	
|	                MISSING            missing             GetMissingSymbol
|	
|	                INTERLEAVE         interleave          IsInterleave
|	
|	                TAXLABELS          (stored in the      (access through
|					                   NxsTaxaBlock        data member taxa)
|									   object)  
|	
|	MATRIX                             matrix              GetDistance
|	                                                       IsMissing
|	                                                       SetMissing
|	                                                       SetDistance
|	------------------------------------------------------------------------
|>
*/
class NxsDistancesBlock
  : public NxsBlock
	{
	public:
							NxsDistancesBlock(NxsTaxaBlock *t);
		virtual				~NxsDistancesBlock();

		double				GetDistance(unsigned i, unsigned j);
		char				GetMissingSymbol();
		unsigned			GetNchar();
		unsigned			GetNtax();
		unsigned			GetTriangle();
		bool				IsRectangular();
		bool				IsDiagonal();
		bool				IsInterleave();
		bool				IsLabels();
		bool				IsLowerTriangular();
		bool				IsMissing(unsigned i, unsigned j);
		bool				IsUpperTriangular();
		virtual void		Report(std::ostream &out);
		virtual void		Reset();
		void				SetDistance(unsigned i, unsigned j, double d);
		void				SetMissing(unsigned i, unsigned j);
		void				SetNchar(unsigned i);

		enum NxsDistancesBlockEnum		/* used by data member triangle to determine which triangle(s) of the distance matrix is/are occupied */
			{
			upper			= 1,		/* matrix is upper-triangular */
			lower			= 2,		/* matrix is lower-triangular */
			both			= 3			/* matrix is rectangular */
			};

	protected:

		void				HandleDimensionsCommand(NxsToken &token);
		void				HandleFormatCommand(NxsToken &token);
		void				HandleMatrixCommand(NxsToken &token);
		bool				HandleNextPass(NxsToken &token, unsigned &offset);
		void				HandleTaxlabelsCommand(NxsToken &token);
		virtual void		Read(NxsToken &token);

	private:

		NxsTaxaBlock		*taxa;		/* pointer to NxsTaxaBlock object that stores the taxon labels */

		bool				newtaxa;	/* true if new taxa were named in this DISTANCES block */
		unsigned			ntax;		/* number of taxa (determines dimensions of the matrix) */
		unsigned			nchar;		/* the number of characters used in generating the pairwise distances */

		bool				diagonal;	/* true if diagonal elements provided when reading in DISTANCES block */
		bool				interleave;	/* true if interleave format used when reading in DISTANCES block */
		bool				labels;		/* true if taxon labels were provided when reading in DISTANCES block */

		int					triangle;	/* indicates whether matrix is upper triangular, lower triangular, or rectangular, taking on one of the elements of the NxsDistancesBlockEnum enumeration */

		char				missing;	/* the symbol used to represent missing data (e.g. '?') */

		NxsDistanceDatum	**matrix;	/* the structure used for storing the pairwise distance matrix */
		unsigned			*taxonPos;	/* array holding 0-offset index into the NxsTaxaBlock list of taxon labels (used to ensure that order of taxa is same for each interleaved block) */
	};

typedef NxsDistancesBlock	DistancesBlock;
#define IsBoth				IsRectangular

#endif

