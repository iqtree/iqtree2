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

#ifndef NCL_NXSCHARACTERSBLOCK_H
#define NCL_NXSCHARACTERSBLOCK_H

#include <iostream> //for std::ostream

class NxsTaxaBlock;
class NxsAssumptionsBlock;

/*----------------------------------------------------------------------------------------------------------------------
|	This class handles reading and storage for the NEXUS block CHARACTERS. It overrides the member functions Read and 
|	Reset, which are abstract virtual functions in the base class NxsBlock. The issue of bookkeeping demands a careful
|	explanation. Users are allowed to control the number of characters analyzed either by "eliminating" or "excluding"
|	characters. Characters can be eliminated (by using the ELIMINATE command) at the time of execution of the data 
|	file, but not thereafter. Characters can, however, be excluded at any time after the data are read. No storage is 
|	provided for eliminated characters, whereas excluded characters must be stored because at any time they could be 
|	restored to active status. Because one can depend on eliminated characters continuing to be eliminated, it would 
|	be inefficient to constantly have to check whether a character has been eliminated. Hence, the characters are 
|	renumbered so that one can efficiently traverse the entire range of non-eliminated characters. The original range 
|	of characters will be hereafter denoted [0..`ncharTotal'), whereas the new, reduced range will be denoted 
|	[0..`nchar'). The two ranges exactly coincide if `ncharTotal' = `nchar' (i.e., no ELIMINATE command was specified
|	in the CHARACTERS block. The possibility for eliminating and excluding characters creates a very confusing situation
|	that is exacerbated by the fact that character indices used in the code begin at 0 whereas character numbers in the
|	data file begin at 1. The convention used hereafter will be to specify "character number k" when discussing 
|	1-offset character numbers in the data file and either "character index k" or simply "character k" when discussing 
|	0-offset character indices.
|	
|	There are several functions (and data structures) that provide services related to keeping track of the 
|	correspondence between character indices in the stored data matrix compared to character numbers in the original 
|	data file. The array `charPos' can be used to find the index of one of the original characters in the matrix. 
|	The function GetCharPos provides public access to the protected `charPos' array. For example, if character 9 
|	(= character number 10) was the only one eliminated, GetCharPos(9) would return UINT_MAX indicating that that 
|	character 9 does not now exist. GetCharPos(10) returns 9 indicating that character 10 in the data file corresponds 
|	to character 9 in the stored data matrix. All public functions in which a character number must be supplied (such 
|	as GetInternalRepresentation) assume that the character number is the current position of the character in the data
|	matrix. This allows one to quickly traverse the data matrix without having to constantly check whether or not a 
|	character was eliminated. Note that GetNChar returns `nchar', not `ncharTotal', and this function should be used 
|	to obtain the end point for a traversal of characters in the matrix. Other functions requiring a (current) character
|	index are: 
|>
|	GetInternalRepresentation
|	GetNumStates
|	GetNumStates
|	GetObsNumStates
|	GetOrigCharIndex
|	GetOrigCharNumber
|	GetState
|	HandleNextState
|	HandleTokenState
|	IsGapState
|	IsMissingState
|	IsPolymorphic
|	ShowStateLabels
|>
|	The function IsEliminated is exceptional in requiring (by necessity) the original character index. The function 
|	GetOrigCharIndex returns the original character index for any current character index. This is useful only when 
|	outputting information that will be seen by the user, and in this case, it is really the character number that 
|	should be output. To get the original character number, either add 1 to GetOrigCharIndex or call GetOrigCharNumber
|	function (which simply returns GetOrigCharIndex + 1).
|	
|	A character may be excluded by calling the function ExcludeCharacter and providing the current character index or 
|	by calling the function ApplyExset and supplying an exclusion set comprising original character indices. These 
|	functions manipulate a bool array, `activeChar', which can be queried using one of two functions: IsActiveChar
|	or IsExcluded. The array `activeChar' is `nchar' elements long, so IsActiveChar and IsExcluded both accept only 
|	current character indices. Thus, a normal loop through all characters in the data matrix should look something 
|	like this:
|>
|	for(unsigned j = 0; j < nchar; j++)
|		{
|		if (IsExcluded(j))
|			continue;
|		.
|		.
|		.
|		}
|>
|	A corresponding set of data structures and functions exists to provide the same services for taxa. Thus, `ntax'
|	 holds the current number of taxa, whereas `ntaxTotal' holds the number of taxa specified in the TAXA block. 
|	If data is provided in the MATRIX command for all taxa listed in the TAXA block, ntax will be equal to `ntaxTotal'.
|	If data is not provided for some of the taxa, the ones left out are treated just like eliminated characters. The 
|	function GetTaxonPos can be used to query the `taxonPos' array, which behaves like the `charPos' array does for 
|	characters: UINT_MAX for element `i' means that the taxon whose original index was `i' has been eliminated and no 
|	data is stored for it in the matrix. Otherwise, GetTaxonPos(i) returns the current index corresponding to the taxon 
|	with an original index of `i'. The function GetNTax returns `ntax', whereas GetNTaxTotal must be used to gain 
|	access to `ntaxTotal' (but this is seldom necessary). The functions GetOrigTaxonIndex and GetOrigTaxonNumber behave 
|	like their character counterparts, GetOrigCharIndex and GetOrigCharNumber. Like characters, taxa can be temporarily
|	inactivated so that they do not participate in any analyses.until they are reactivated by the user. Inactivation 
|	of a taxon is refered to as deleting the taxon, whereas restoring a taxon means reactivating it. Thus, the 
|	ApplyDelset, DeleteTaxon, and RestoreTaxon functions correspond to the ApplyExset, ExcludeCharacter, and 
|	IncludeCharacter functions for characters. To query whether a taxon is currently deleted, use either 
|	IsActiveTaxon or IsDeleted. A normal loop across all active taxa can be constructed as follows:
|>
|	for (unsigned i = 0; i < ntax; i++)
|		{
|		if (IsDeleted(i))
|			continue;
|		.
|		.
|		.
|		}
|>
|	Below is a table showing the correspondence between the elements of a CHARACTERS block in a NEXUS file and the
|	variables and member functions of the NxsCharactersBlock class that can be used to access each piece of information
|	stored. Items in parenthesis should be viewed as "see also" items.
|>
|	NEXUS         Command        Data           Member
|   Command       Atribute       Member         Functions
|	---------------------------------------------------------------------
|	DIMENSIONS    NEWTAXA        newtaxa
|	
|	              NTAX           ntax           GetNTax
|                                (ntaxTotal)    (GetNumMatrixRows)
|	
|	              NCHAR          nchar          GetNChar
|	                             (ncharTotal)   (GetNumMatrixCols)
|	
|	FORMAT        DATATYPE       datatype       GetDataType
|	
|	              RESPECTCASE    respectingCase IsRespectCase
|	
|	              MISSING        missing        GetMissingSymbol
|	
|	              GAP            gap            GetGapSymbol
|	
|	              SYMBOLS        symbols        GetSymbols
|	
|	              EQUATE         equates        GetEquateKey
|	                                            GetEquateValue
|	                                            GetNumEquates
|	
|	              MATCHCHAR      matchchar      GetMatchcharSymbol
|	
|	              (NO)LABELS     labels         IsLabels
|	
|	              TRANSPOSE      transposing    IsTranspose
|	
|	              INTERLEAVE     interleaving   IsInterleave
|	
|	              ITEMS          (Note: only STATES implemented)
|	
|	              STATESFORMAT   (Note: only STATESPRESENT implemented)
|	
|	              (NO)TOKENS     tokens         IsTokens
|	
|	ELIMINATE                    eliminated     IsEliminated
|	                                            GetNumEliminated
|	
|	MATRIX                       matrix         GetState
|	                                            GetInternalRepresentation
|	                                            GetNumStates
|	                                            GetNumMatrixRows
|	                                            GetNumMatrixCols
|	                                            IsPolymorphic
|>
*/
class NxsString;      //James B. (Added forward declare needed in VS builds, 23-Jul-2020)
class NxsCharactersBlock
  : public NxsBlock
	{
	friend class NxsAssumptionsBlock;

	public:

		enum DataTypesEnum		/* values used to represent different basic types of data stored in a CHARACTERS block, and used with the data member `datatype' */
			{
			standard = 1,		/* indicates `matrix' holds characters with arbitrarily-assigned, discrete states, such as discrete morphological data */
			dna,				/* indicates `matrix' holds DNA sequences (states A, C, G, T) */
			rna,				/* indicates `matrix' holds RNA sequences (states A, C, G, U) */
			nucleotide,			/* indicates `matrix' holds nucleotide sequences */
			protein,			/* indicates `matrix' holds amino acid sequences */
			continuous			/* indicates `matrix' holds continuous data */
			};

								NxsCharactersBlock(NxsTaxaBlock *tb, NxsAssumptionsBlock *ab);
		virtual					~NxsCharactersBlock();

		unsigned				ApplyDelset(NxsUnsignedSet &delset);
		unsigned				ApplyExset(NxsUnsignedSet &exset);
		unsigned				ApplyIncludeset(NxsUnsignedSet &inset);
		unsigned				ApplyRestoreset(NxsUnsignedSet &restoreset);
		unsigned				GetCharPos(unsigned origCharIndex);
		unsigned				GetTaxPos(unsigned origTaxonIndex);
		unsigned				GetDataType();
		int						GetInternalRepresentation(unsigned i, unsigned j, unsigned k = 0);
		unsigned				GetNTax();
		unsigned				GetNChar();
		unsigned				GetNCharTotal();
		unsigned				GetNTaxTotal();
		unsigned				GetNumActiveChar();
		unsigned				GetNumActiveTaxa();
		unsigned				GetNumEliminated();
		unsigned				GetNumEquates();
		unsigned				GetNumMatrixCols();
		unsigned				GetNumMatrixRows();
		unsigned				GetNumStates(unsigned i, unsigned j);
		unsigned				GetOrigCharIndex(unsigned j);
		unsigned				GetOrigCharNumber(unsigned j);
		unsigned				GetOrigTaxonIndex(unsigned j);
		unsigned				GetOrigTaxonNumber(unsigned j);
		char					GetGapSymbol();
		char					GetMatchcharSymbol();
		char					GetMissingSymbol();
        NxsDiscreteMatrix       *GetMatrix() { return matrix; }
		bool					IsGapState(unsigned i, unsigned j);
		bool					IsInterleave();
		bool					IsLabels();
		bool					IsMissingState(unsigned i, unsigned j);
		bool					IsPolymorphic(unsigned i, unsigned j);
		bool					IsRespectCase();
		bool					IsTokens();
		bool					IsTranspose();
		bool					IsEliminated(unsigned origCharIndex);
		void					Consume(NxsCharactersBlock &other);
		void					ExcludeCharacter(unsigned i);
		void					IncludeCharacter(unsigned i);
		bool					IsActiveChar(unsigned j);
		bool					IsExcluded(unsigned j);
		void					DeleteTaxon(unsigned i);
		void					RestoreTaxon(unsigned i);
		bool					IsActiveTaxon(unsigned i);
		bool					IsDeleted(unsigned i);
		void					ShowStateLabels(std::ostream &out, unsigned i, unsigned c, unsigned first_taxon = -1);
		unsigned				GetStateSymbolIndex(unsigned i, unsigned j, unsigned k = 0);	// added by mth for standard data types
		char					GetState(unsigned i, unsigned j, unsigned k = 0);
		char					*GetSymbols();
		bool					*GetActiveTaxonArray();
		bool					*GetActiveCharArray();
		NxsString				GetCharLabel(unsigned i);
		NxsString				GetStateLabel(unsigned i, unsigned j);
		NxsString				GetTaxonLabel(unsigned i);
		virtual unsigned		CharLabelToNumber(NxsString s);
		virtual unsigned		TaxonLabelToNumber(NxsString s);
		virtual unsigned		GetMaxObsNumStates();
		virtual unsigned		GetObsNumStates(unsigned j);
		virtual void			DebugShowMatrix(std::ostream &out, bool use_matchchar, const char *marginText = 0);
		virtual void			Report(std::ostream &out);
		virtual void			Reset();

		NxsTaxaBlock			*taxa;				/* pointer to the TAXA block in which taxon labels are stored */

	protected:

		void					BuildCharPosArray(bool check_eliminated = false);
		bool					IsInSymbols(char ch);
		void					HandleCharlabels(NxsToken &token);
		void					HandleCharstatelabels(NxsToken &token);
		void					HandleDimensions(NxsToken &token, NxsString newtaxaLabel, NxsString ntaxLabel, NxsString ncharLabel);
		void					HandleEliminate(NxsToken &token);
		void					HandleEndblock(NxsToken &token, NxsString charToken);
		virtual void			HandleFormat(NxsToken &token);
		virtual void			HandleMatrix(NxsToken &token);
		virtual bool			HandleNextState(NxsToken &token, unsigned i, unsigned c);
		virtual void			HandleStdMatrix(NxsToken &token);
		virtual unsigned		HandleTokenState(NxsToken &token, unsigned c);
		virtual void			HandleTransposedMatrix(NxsToken &token);
		virtual void			Read(NxsToken &token);
		unsigned				PositionInSymbols(char ch);
		void					HandleStatelabels(NxsToken &token);
		void					HandleTaxlabels(NxsToken &token);
		void					ResetSymbols();
		void					ShowStates(ostream &out, unsigned i, unsigned j);
		void					WriteStates(NxsDiscreteDatum &d, char *s, unsigned slen);

		NxsAssumptionsBlock		*assumptionsBlock;	/* pointer to the ASSUMPTIONS block in which exsets, taxsets and charsets are stored */

		unsigned				ntax;				/* number of rows in matrix (same as `ntaxTotal' unless fewer taxa appeared in CHARACTERS MATRIX command than were specified in the TAXA block, in which case `ntaxTotal' > `ntax') */
		unsigned				ntaxTotal;			/* number of taxa (same as `ntax' unless fewer taxa appeared in CHARACTERS MATRIX command than were specified in the TAXA block, in which case `ntaxTotal' > `ntax') */
		unsigned				nchar;				/* number of columns in matrix (same as `ncharTotal' unless some characters were eliminated, in which case `ncharTotal' > `nchar') */
		unsigned				ncharTotal;			/* total number of characters (same as `nchar' unless some characters were eliminated, in which case `ncharTotal' > `nchar') */

		bool					newtaxa;			/* true if NEWTAXA keyword encountered in DIMENSIONS command */
		bool					newchar;			/* true unless CHARLABELS or CHARSTATELABELS command read */

		bool					formerly_datablock;	/* true if this object was originally read in as a DATA block rather than as a CHARACTERS block, false otherwise */
		bool					respectingCase;		/* if true, RESPECTCASE keyword specified in FORMAT command */
		bool					transposing;		/* indicates matrix will be in transposed format */
		bool					interleaving;		/* indicates matrix will be in interleaved format */
		bool					tokens;				/* if false, data matrix entries must be single symbols; if true, multicharacter entries are allows */
		bool					labels;				/* indicates whether or not labels will appear on left side of matrix */

		char					missing;			/* missing data symbol */
		char					gap;				/* gap symbol for use with molecular data */
		char					matchchar;			/* match symbol to use in matrix */

		char					*symbols;			/* list of valid character state symbols */

		NxsStringMap			equates;			/* list of associations defined by EQUATE attribute of FORMAT command */

		NxsDiscreteMatrix		*matrix;			/* storage for discrete data */
		unsigned				*charPos;			/* maps character numbers in the data file to column numbers in matrix (necessary if some characters have been eliminated) */
		unsigned				*taxonPos;			/* maps taxon numbers in the data file to row numbers in matrix (necessary if fewer taxa appear in CHARACTERS block MATRIX command than are specified in the TAXA block) */
		NxsUnsignedSet			eliminated;			/* array of (0-offset) character numbers that have been eliminated (will remain empty if no ELIMINATE command encountered) */

		bool					*activeChar;		/* `activeChar[i]' true if character `i' not excluded; `i' is in range [0..`nchar') */
		bool					*activeTaxon;		/* `activeTaxon[i]' true if taxon `i' not deleted; `i' is in range [0..`ntax') */

		NxsStringVector			charLabels;			/* storage for character labels (if provided) */
		NxsStringVectorMap		charStates;			/* storage for character state labels (if provided) */

	private:

		DataTypesEnum			datatype;			/* flag variable (see datatypes enum) */
	};

typedef NxsCharactersBlock CharactersBlock;

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes taxon whose 0-offset current index is `i'. If taxon has already been deleted, this function has no effect.
*/
inline void NxsCharactersBlock::DeleteTaxon(
  unsigned i)	/* index of taxon to delete in range [0..`ntax') */
	{
	activeTaxon[i] = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Excludes character whose 0-offset current index is `i'. If character has already been excluded, this function has 
|	no effect.
*/
inline void NxsCharactersBlock::ExcludeCharacter(
  unsigned i)	/* index of character to exclude in range [0..`nchar') */
	{
	activeChar[i] = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `activeChar' data member (pointer to first element of the `activeChar' array). Access to this protected 
|	data member is necessary in certain circumstances, such as when a NxsCharactersBlock object is stored in another 
|	class, and that other class needs direct access to the `activeChar' array even though it is not derived from 
|	NxsCharactersBlock.
*/
inline bool *NxsCharactersBlock::GetActiveCharArray()
	{
		return activeChar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `activeTaxon' data member (pointer to first element of the `activeTaxon' array). Access to this protected 
|	data member is necessary in certain circumstances, such as when a NxsCharactersBlock object is stored in another 
|	class, and that other class needs direct access to the `activeTaxon' array even though it is not derived from 
|	NxsCharactersBlock.
*/
inline bool *NxsCharactersBlock::GetActiveTaxonArray()
	{
	return activeTaxon;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns label for character `i', if a label has been specified. If no label was specified, returns string 
|	containing a single blank (i.e., " ").
*/
inline NxsString NxsCharactersBlock::GetCharLabel(
  unsigned i)	/* the character in range [0..`nchar') */
	{
	NxsString s = " ";
	if (static_cast<unsigned>(i) < charLabels.size())
		s = charLabels[i];
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current index of character in matrix. This may differ from the original index if some characters were 
|	removed using an ELIMINATE command. For example, character number 9 in the original data matrix may now be at 
|	position 8 if the original character 8 was eliminated. The parameter `origCharIndex' is assumed to range from 
|	0 to `ncharTotal' - 1.
*/
inline unsigned NxsCharactersBlock::GetCharPos(
  unsigned origCharIndex)	/* original index of character in range [0..`ncharTotal' - 1) */
	{
	assert(charPos);
	assert(origCharIndex >= 0);
	assert(origCharIndex < ncharTotal);

	return charPos[origCharIndex];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the gap symbol currently in effect. If no gap symbol specified, returns '\0'.
*/
inline char NxsCharactersBlock::GetGapSymbol()
	{
	return gap;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current index of taxon in matrix. This may differ from the original index if some taxa were listed in the 
|	TAXA block but not in the DATA or CHARACTERS block. The parameter `origTaxonIndex' is assumed to range from 0 to 
|	`ntaxTotal' - 1.
*/
inline unsigned NxsCharactersBlock::GetTaxPos(
  unsigned origTaxonIndex)	/* original index of taxon */
	{
	assert(taxonPos);
	assert(origTaxonIndex >= 0);
	assert(origTaxonIndex < ntaxTotal);

	return taxonPos[origTaxonIndex];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `datatype'.
*/
inline unsigned NxsCharactersBlock::GetDataType()
	{
	return (int)datatype;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `matchchar' symbol currently in effect. If no `matchchar' symbol specified, returns '\0'.
*/
inline char NxsCharactersBlock::GetMatchcharSymbol()
	{
	return matchchar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns internal representation of the state for taxon `i', character `j'. In the normal situation, `k' is 0 meaning
|	there is only one state with no uncertainty or polymorphism. If there are multiple states, specify a number in the 
|	range [0..n) where n is the number of states returned by the GetNumStates function. Use the IsPolymorphic 
|	function to determine whether the multiple states correspond to uncertainty in state assignment or polymorphism in 
|	the taxon. The value returned from this function is one of the following:
|~
|	o -3 means gap state (see note below)
|	o -2 means missing state (see note below)
|	o an integer 0 or greater is internal representation of a state
|~
|	Note: gap and missing states are actually represented internally in a different way; for a description of the actual
|	internal representation of states, see the documentation for NxsDiscreteDatum.
*/
inline int NxsCharactersBlock::GetInternalRepresentation(
  unsigned i,	/* the taxon in range [0..`ntax') */
  unsigned j,	/* the character in range [0..`nchar') */
  unsigned k)	/* the 0-offset index of state to return */
	{
	if (IsGapState(i, j))
		return -3;
	else if (IsMissingState(i, j))
		return -2;
	else
		return matrix->GetState(i, j, k);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the missing data symbol currently in effect. If no missing data symbol specified, returns '\0'.
*/
inline char NxsCharactersBlock::GetMissingSymbol()
	{
	return missing;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `nchar'.
*/
inline unsigned NxsCharactersBlock::GetNChar()
	{
	return nchar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ncharTotal'.
*/
inline unsigned NxsCharactersBlock::GetNCharTotal()
	{
	return ncharTotal;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ntax'.
*/
inline unsigned NxsCharactersBlock::GetNTax()
	{
	return ntax;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ntaxTotal'.
*/
inline unsigned NxsCharactersBlock::GetNTaxTotal()
	{
	return ntaxTotal;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of characters eliminated with the ELIMINATE command.
*/
inline unsigned NxsCharactersBlock::GetNumEliminated()
	{
	return (ncharTotal - nchar);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of stored equate associations.
*/
inline unsigned NxsCharactersBlock::GetNumEquates()
	{
	return static_cast<unsigned>(equates.size());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of actual columns in `matrix'. This number is equal to `nchar', but can be smaller than 
|	`ncharTotal' since the user could have eliminated some of the characters.
*/
inline unsigned NxsCharactersBlock::GetNumMatrixCols()
	{
	return nchar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of actual rows in `matrix'. This number is equal to `ntax', but can be smaller than `ntaxTotal'
|	since the user did not have to provide data for all taxa specified in the TAXA block.
*/
inline unsigned NxsCharactersBlock::GetNumMatrixRows()
	{
	return ntax;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of states for taxon `i', character `j'.
*/
inline unsigned NxsCharactersBlock::GetNumStates(
  unsigned i,	/* the taxon in range [0..`ntax') */
  unsigned j)	/* the character in range [0..`nchar') */
	{
	return matrix->GetNumStates(i, j);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of states for character `j' over all taxa. Note: this function is rather slow, as it must walk 
|	through each row, adding the states encountered to a set, then finally returning the size of the set. Thus, if this 
|	function is called often, it would be advisable to initialize an array using this function, then refer to the array 
|	subsequently.
*/
inline unsigned NxsCharactersBlock::GetObsNumStates(
  unsigned j)	/* the character in range [0..`nchar') */
	{
	return matrix->GetObsNumStates(j);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the original character number (used in the NEXUS data file) in the range [1..`ncharTotal']. Will be equal 
|	to `j' + 1 unless some characters were eliminated.
*/
inline unsigned NxsCharactersBlock::GetOrigCharNumber(
  unsigned j)	/* the character in range [0..`nchar') */
	{
	return (1 + GetOrigCharIndex(j));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the original taxon number (used in the NEXUS data file) in the range [1..`ntaxTotal']. Will be equal to 
|	`i' + 1 unless data was not provided for some taxa listed in a preceding TAXA block.
*/
inline unsigned NxsCharactersBlock::GetOrigTaxonNumber(
  unsigned i)	/* the character in range [0..`ntax') */
	{
	return (1 + GetOrigTaxonIndex(i));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns symbol from symbols list representing the state for taxon `i' and character `j'. The normal situation in 
|	which there is only one state with no uncertainty or polymorphism is represented by `k' = 0. If there are multiple 
|	states, specify a number in the range [0..n) where n is the number of states returned by the GetNumStates function.
|	Use the IsPolymorphic function to determine whether the multiple states correspond to uncertainty in state 
|	assignment or polymorphism in the taxon. Assumes `symbols' is non-NULL.
*/
inline char NxsCharactersBlock::GetState(
  unsigned i,	/* the taxon in range [0..`ntax') */
  unsigned j,	/* the character in range [0..`nchar') */
  unsigned k)	/* the 0-offset index of the state to return */
	{
	assert(symbols);
	char state_char = '\0';

	//unsigned symbolsLen = strlen(symbols);
	unsigned p = matrix->GetState(i, j, k);
	assert(p < strlen(symbols));
	state_char = *(symbols + p);

	return state_char;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns data member `symbols'.  Warning: returned value may be NULL.
*/
inline char *NxsCharactersBlock::GetSymbols()
	{
	return symbols;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns label for taxon number `i' (`i' ranges from 0 to `ntax' - 1).
*/
inline NxsString NxsCharactersBlock::GetTaxonLabel(
  unsigned i)	/* the taxon's position in the taxa block */
	{
	NxsString s = taxa->GetTaxonLabel(i);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Includes character whose 0-offset current index is `i'. If character is already active, this function has no effect.
*/
inline void NxsCharactersBlock::IncludeCharacter(
  unsigned i)	/* index of character to include in range [0..`nchar') */
	{
	activeChar[i] = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if character `j' is active. If character `j' has been excluded, returns false. Assumes `j' is in the 
|	range [0..`nchar').
*/
inline bool NxsCharactersBlock::IsActiveChar(
  unsigned j)	/* the character in question, in the range [0..`nchar') */
	{
	assert(j >= 0);
	assert(j < nchar);

	return activeChar[j];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxon `i' is active. If taxon `i' has been deleted, returns false. Assumes `i' is in the range 
|	[0..`ntax').
*/
inline bool NxsCharactersBlock::IsActiveTaxon(
  unsigned i)	/* the taxon in question, in the range [0..`ntax') */
	{
	assert(i >= 0);
	assert(i < ntax);

	return activeTaxon[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxon number i has been deleted, false otherwise.
*/
inline bool NxsCharactersBlock::IsDeleted(
  unsigned i)	/* the taxon in question, in the range [0..`ntax') */
	{
	return !IsActiveTaxon(i);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if character `j' has been excluded. If character `j' is active, returns false. Assumes `j' is in the 
|	range [0..`nchar').
*/
inline bool NxsCharactersBlock::IsExcluded(
  unsigned j)	/* the character in question, in the range [0..`nchar') */
	{
	return !IsActiveChar(j);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the state at taxon `i', character `j' is the gap state, false otherwise. Assumes `matrix' is 
|	non-NULL.
*/
inline bool NxsCharactersBlock::IsGapState(
  unsigned i,	/* the taxon, in range [0..`ntax') */
  unsigned j)	/* the character, in range [0..`nchar') */
	{
	assert(matrix);
	return matrix->IsGap(i, j);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if INTERLEAVE was specified in the FORMAT command, false otherwise.
*/
inline bool NxsCharactersBlock::IsInterleave()
	{
	return interleaving;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if LABELS was specified in the FORMAT command, false otherwise.
*/
inline bool NxsCharactersBlock::IsLabels()
	{
	return labels;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the state at taxon `i', character `j' is the missing state, false otherwise. Assumes `matrix' is 
|	non-NULL.
*/
inline bool NxsCharactersBlock::IsMissingState(
  unsigned i,	/* the taxon, in range [0..`ntax') */
  unsigned j)	/* the character, in range [0..`nchar') */
	{
	assert(matrix);
	return matrix->IsMissing(i, j);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxon `i' is polymorphic for character `j', false otherwise. Assumes `matrix' is non-NULL. Note 
|	that return value will be false if there is only one state (i.e., one cannot tell whether there is uncertainty 
|	using this function).
*/
inline bool NxsCharactersBlock::IsPolymorphic(
  unsigned i,	/* the taxon in range [0..`ntax') */
  unsigned j)	/* the character in range [0..`nchar') */
	{
	assert(matrix);
	return matrix->IsPolymorphic(i, j);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if RESPECTCASE was specified in the FORMAT command, false otherwise.
*/
inline bool NxsCharactersBlock::IsRespectCase()
	{
	return respectingCase;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if TOKENS was specified in the FORMAT command, false otherwise.
*/
inline bool NxsCharactersBlock::IsTokens()
	{
	return tokens;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if TRANSPOSE was specified in the FORMAT command, false otherwise.
*/
inline bool NxsCharactersBlock::IsTranspose()
	{
	return transposing;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores taxon whose 0-offset current index is `i'. If taxon is already active, this function has no effect.
*/
inline void NxsCharactersBlock::RestoreTaxon(
  unsigned i)	/* index of taxon to restore in range [0..`ntax') */
	{
	activeTaxon[i] = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Shows the states for taxon `i', character `j', on the stream `out'. Uses `symbols' array to translate the states 
|	from the way they are stored (as integers) to the symbol used in the original data matrix. Assumes `i' is in the 
|	range [0..`ntax') and `j' is in the range [0..`nchar'). Also assumes `matrix' is non-NULL.
*/
inline void NxsCharactersBlock::ShowStates(
  ostream &out,	/* the stream on which to show the state(s) */
  unsigned i,	/* the (0-offset) index of the taxon in question */
  unsigned j)	/* the (0-offset) index of the character in question */
	{
	assert(i >= 0);
	assert(i < ntax);
	assert(j >= 0);
	assert(j < nchar);
	assert(matrix);

	char s[NCL_MAX_STATES + 3];
	WriteStates(matrix->GetDiscreteDatum(i, j), s, NCL_MAX_STATES + 3);

	out << s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts a taxon label to a number corresponding to the taxon's position within the list maintained by the 
|	NxsTaxaBlock object. This method overrides the virtual function of the same name in the NxsBlock base class. If 
|	`s' is not a valid taxon label, returns the value 0.
*/
inline unsigned NxsCharactersBlock::TaxonLabelToNumber(
  NxsString s)	/* the taxon label to convert */
	{
	unsigned i;
	try
		{
		i = 1 + taxa->FindTaxon(s);
		}
	catch(NxsTaxaBlock::NxsX_NoSuchTaxon)
		{
		i = 0;
		}

	return i;
	}





#endif
