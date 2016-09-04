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
|	Sets `id' to "DISTANCES", `taxa' to `t', `triangle' to `NxsDistancesBlockEnum::lower', `missing' to '?', `matrix' 
|	and `taxonPos' to NULL, `labels' and `diagonal' to true, `newtaxa' and `interleave' to false, and `ntax' and `nchar'
|	to 0. Assumes `t' is non-NULL.
*/
NxsDistancesBlock::NxsDistancesBlock(
  NxsTaxaBlock *t)	/* the NxsTaxaBlock that will keep track of taxon labels */
  : NxsBlock()
	{
	assert(t != NULL);
	taxa		= t;
	id = "DISTANCES";
	ntax		= 0;
	nchar		= 0;
	diagonal	= true;
	labels		= true;
	newtaxa		= false;
	interleave	= false;
	triangle	= NxsDistancesBlockEnum(lower);
	missing		= '?';
	matrix		= NULL;
	taxonPos	= NULL;
}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes `matrix' and `taxonPos' arrays.
*/
NxsDistancesBlock::~NxsDistancesBlock()
	{
	if (matrix != NULL)
		delete matrix;
	if (taxonPos != NULL)
		delete [] taxonPos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when DIMENSIONS command needs to be parsed from within the DISTANCES block. Deals with everything after the 
|	token DIMENSIONS up to and including the semicolon that terminates the DIMENSIONS command.
*/
void NxsDistancesBlock::HandleDimensionsCommand(
  NxsToken &token)	/* the token used to read from `in' */
	{
	for (;;)
		{
		token.GetNextToken();

		// Token should either be ';' or the name of a subcommand
		//
		if (token.Equals(";"))
			break;

		else if (token.Equals("NEWTAXA"))
			{
			ntax = 0;
			newtaxa = 1;
			}

		else if (token.Equals("NTAX"))
			{
			if (!newtaxa)
				{
				errormsg = "Must specify NEWTAXA before NTAX if new taxa are being defined";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			// This should be the equals sign
			//
			token.GetNextToken();

			if (!token.Equals("="))
				{
				errormsg = "Expecting '=' but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			// This should be the number of taxa
			//
			token.GetNextToken();
			ntax = atoi(token.GetToken().c_str());
			}

		else if (token.Equals("NCHAR"))
			{
			// This should be the equals sign
			//
			token.GetNextToken();

			if (!token.Equals("="))
				{
				errormsg = "Expecting '=' but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			// This should be the number of characters
			//
			token.GetNextToken();
			nchar = atoi(token.GetToken().c_str());
			}
		}

	if (ntax == 0)
		ntax = taxa->GetNumTaxonLabels();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when FORMAT command needs to be parsed from within the DISTANCES block. Deals with everything after the 
|	token FORMAT up to and including the semicolon that terminates the FORMAT command.
*/
void NxsDistancesBlock::HandleFormatCommand(
  NxsToken &token)	/* the token used to read from `in' */
	{
	for (;;)
		{
		// This should either be ';' or the name of a subcommand
		//
		token.GetNextToken();

		if (token.Equals(";"))
			break;

		else if (token.Equals("TRIANGLE"))
			{
			// This should be the equals sign
			//
			token.GetNextToken();

			if (!token.Equals("="))
				{
				errormsg = "Expecting '=' but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			// This should be LOWER, UPPER, or BOTH
			//
			token.GetNextToken();

			if (token.Equals("LOWER"))
				triangle = NxsDistancesBlockEnum(lower);
			else if (token.Equals("UPPER"))
				triangle = NxsDistancesBlockEnum(upper);
			else if (token.Equals("BOTH"))
				triangle = NxsDistancesBlockEnum(both);
			else
				{
				errormsg = "Expecting UPPER, LOWER, or BOTH but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			}

		else if (token.Equals("DIAGONAL"))
			{
			diagonal = 1;
			}

		else if (token.Equals("NODIAGONAL"))
			{
			diagonal = 0;
			}

		else if (token.Equals("LABELS"))
			{
			labels = 1;
			}

		else if (token.Equals("NOLABELS"))
			{
			labels = 0;
			}

		else if (token.Equals("INTERLEAVE"))
			{
			interleave = 1;
			}

		else if (token.Equals("NOINTERLEAVE"))
			{
			interleave = 0;
			}

		else if (token.Equals("MISSING"))
			{
			// This should be the equals sign
			//
			token.GetNextToken();

			if (!token.Equals("="))
				{
				errormsg = "Expecting '=' but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			// This should be the missing data symbol
			//
			token.GetNextToken();

			if (token.GetTokenLength() != 1)
				{
				errormsg = "Missing data symbol specified (";
				errormsg += token.GetToken();
				errormsg += ") is invalid (must be a single character)";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			missing = token.GetToken()[0];
			}

		else
			{
			errormsg = "Token specified (";
			errormsg += token.GetToken();
			errormsg += ") is an invalid subcommand for the FORMAT command";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called from within HandleMatrix, this function is used to deal with interleaved matrices. It is called once for 
|	each pass through the taxa. The local variable `jmax' records the number of columns read in the current interleaved 
|	page and is used to determine the offset used for j in subsequent pages.
*/
bool NxsDistancesBlock::HandleNextPass(
  NxsToken &token,	/* the token we are using for reading the data file */
  unsigned &offset)	/* the offset */
	{
	unsigned i, j, k, jmax = 0; 
	bool done = false;

	unsigned i_first = 0;
	if (triangle == NxsDistancesBlockEnum(lower))
		i_first = offset;

	unsigned i_last = ntax;

	for (i = i_first; i < i_last; i++)
		{
		// Deal with taxon label if provided. Here are the four situations we need to deal with:
		//   newtaxa  (offset > 0)  handled by
		//      0           0         case 1
		//      0           1         case 1
		//      1           0         case 2
		//      1           1         case 1
		//
		if (labels && (!newtaxa || offset > 0))
			{
			// Case 1: Expecting taxon labels, and also expecting them to already be in taxa
			//
			do
				{
				token.SetLabileFlagBit(NxsToken::newlineIsToken);
				token.GetNextToken();
				}
			while(token.AtEOL());

			try
				{
				// Look up position of taxon in NxsTaxaBlock list
				//
				k = taxa->FindTaxon(token.GetToken());

				// Array taxonPos is initialized to UINT_MAX and filled in as taxa are encountered
				//
				if (taxonPos[i] == UINT_MAX)
					{
					taxonPos[i] = k;
					}
				else if (taxonPos[i] != k)
					{
					errormsg = "Taxon labeled ";
					errormsg += token.GetToken();
					errormsg += " is out of order compared to previous interleave pages";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
				}

			catch (NxsTaxaBlock::NxsX_NoSuchTaxon)
				{
				errormsg = "Could not find ";
				errormsg += token.GetToken();
				errormsg += " among taxa previously defined";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			}

		else if (labels && newtaxa)
			{
			// Case 2: Expecting taxon labels, and also expecting taxa block to be empty
			//
			do
				{
				token.SetLabileFlagBit(NxsToken::newlineIsToken);
				token.GetNextToken();
				}
			while(token.AtEOL());

			taxa->AddTaxonLabel(token.GetToken());
			taxonPos[i] = i;
			}

		// Now deal with the row of distance values
		//
		unsigned true_j = 0;
		for (j = 0; j < ntax; j++)
			{
			if (i == ntax - 1 && j == ntax - 1)
				{
				done = true;
				}

			if ((i == ntax - 1) && (true_j == ntax - 1))
				{
				done = true;
				break;
				}

			if (i == ntax-1 && !diagonal && triangle == NxsDistancesBlockEnum(upper))
				{
				done = true;
				break;
				}

			if (!diagonal && triangle == NxsDistancesBlockEnum(lower) && j == ntax - offset - 1)
				{
				done = true;
				break;
				}

			token.SetLabileFlagBit(NxsToken::newlineIsToken);
			token.GetNextToken();

			if (token.AtEOL())
				{
				if (j > jmax)
					{
					jmax = j;
					if (!diagonal && triangle == NxsDistancesBlockEnum(upper) && i >= offset)
						jmax++;
					if (interleave && triangle == NxsDistancesBlockEnum(upper))
						i_last = jmax + offset;
					}
				break;
				}

			true_j = j + offset;
			if (triangle == NxsDistancesBlockEnum(upper) && i > offset)
				true_j += (i - offset);
			if (!diagonal && triangle == NxsDistancesBlockEnum(upper) && i >= offset)
				true_j++;

			if (true_j == ntax)
				{
				errormsg = "Too many distances specified in row just read in";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			string t = token.GetToken();
			if (token.GetTokenLength() == 1 && t[0] == missing)
				SetMissing(i, true_j);
			else
				SetDistance(i, true_j, atof(t.c_str()));
			}
		}

	offset += jmax;

	return done;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when MATRIX command needs to be parsed from within the DISTANCES block. Deals with everything after the 
|	token MATRIX up to and including the semicolon that terminates the MATRIX command.
*/
void NxsDistancesBlock::HandleMatrixCommand(
  NxsToken &token)	/* the token used to read from `in' */
	{
	unsigned i;
	unsigned prev_ntax = ntax;

	if (ntax == 0)
		ntax = taxa->GetNumTaxonLabels();

	if (ntax == 0)
		{
		errormsg = "MATRIX command cannot be read if NTAX is zero";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	if (triangle == NxsDistancesBlockEnum(both) && !diagonal)
		{
		errormsg = "Cannot specify NODIAGONAL and TRIANGLE=BOTH at the same time";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	if (newtaxa)
		taxa->Reset();

	// Allocate taxonPos array, deleting it first if previously allocated
	//
	if (taxonPos != NULL)
		{
		delete [] taxonPos;
		}

	taxonPos = new unsigned[ntax];
	
	for (i = 0; i < ntax; i++)
		taxonPos[i] = UINT_MAX;

	// Allocate matrix array, deleting it first if previously allocated
	//
	if (matrix != NULL)
		{
		assert(prev_ntax > 0);
		for (i = 0; i < prev_ntax; i++)
			delete [] matrix[i];
		delete [] matrix;
		}

	matrix = new NxsDistanceDatum*[ntax];
	for (i = 0; i < ntax; i++)
		matrix[i] = new NxsDistanceDatum[ntax];

	unsigned offset = 0;
	bool done = false;
	while (!done)
		{
		done = HandleNextPass(token, offset);
		}

	// Token should be equal to the terminating semicolon
	//
	token.GetNextToken();

	if (!token.Equals(";"))
		{
		errormsg = "Expecting ';' to terminate MATRIX command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when TAXLABELS command needs to be parsed from within the DISTANCES block. Deals with everything after the 
|	token TAXLABELS up to and including the semicolon that terminates the TAXLABELS command.
*/
void NxsDistancesBlock::HandleTaxlabelsCommand(
  NxsToken &token)	/* the token used to read from `in' */
	{
	if (!newtaxa)
		{
		errormsg = "NEWTAXA must have been specified in DIMENSIONS command to use the TAXLABELS command in a ";
		errormsg += id;
		errormsg += " block";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	if (ntax == 0)
		{
		errormsg = "NTAX must be specified before TAXLABELS command";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	for (unsigned i = 0; i < ntax; i++)
		{
        token.SetLabileFlagBit(NxsToken::hyphenNotPunctuation + NxsToken::preserveUnderscores);
		token.GetNextToken();
		taxa->AddTaxonLabel(token.GetToken());
		}

	// This should be terminating semicolon
	//
	token.GetNextToken(); 

	if (!token.Equals(";"))
		{
		errormsg = "Expecting ';' to terminate TAXLABELS command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	// Some may object to setting newtaxa to false here, because then the
	// fact that new taxa were specified in this DISTANCES block rather than in
	// a preceding TAXA block is lost.  This will only be important if we wish to
	// recreate the original data file, which I don't anticipate anyone doing with
	// this code (too difficult to remember all comments, the order of blocks in
	// the file, etc.)
	//
	newtaxa = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the ability to read everything following the block name (which is read by the NEXUS object)
|	to the end or endblock statement. Characters are read from the input stream in. Overrides the abstract virtual 
|	function in the base class.
*/
void NxsDistancesBlock::Read(
  NxsToken &token)	/* the token used to read from `in' */
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

	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("DIMENSIONS"))
			{
			HandleDimensionsCommand(token);
			}

		else if (token.Equals("FORMAT"))
			{
			HandleFormatCommand(token);
			}

		else if (token.Equals("TAXLABELS"))
			{
			HandleTaxlabelsCommand(token);
			}

		else if (token.Equals("MATRIX"))
			{
			HandleMatrixCommand(token);
			}

		else if (token.Equals("END"))
			{
			// Get the semicolon following END
			//
			token.GetNextToken();

			if (!token.Equals(";"))
				{
				errormsg = "Expecting ';' to terminate the END command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			break;
			}

		else if (token.Equals("ENDBLOCK"))
			{
			// Get the semicolon following ENDBLOCK
			//
			token.GetNextToken();

			if (!token.Equals(";"))
				{
				errormsg = "Expecting ';' to terminate the ENDBLOCK command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

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
|	This function outputs a brief report of the contents of this taxa block. Overrides the abstract virtual function in 
|	the base class.
*/
void NxsDistancesBlock::Report(
  ostream &out)	/* the output stream to which to write the report */
	{
	unsigned ntaxTotal = ntax;

	if (ntaxTotal == 0)
		ntaxTotal = taxa->GetNumTaxonLabels();

	out << endl;
	out << id << " block contains ";
	if (ntaxTotal == 0)
		{
		out << "no taxa" << endl;
		}
	else if (ntaxTotal == 1)
		out << "one taxon" << endl;
	else
		out << ntaxTotal << " taxa" << endl;

	if (IsLowerTriangular())
		out << "  Matrix is lower-triangular" << endl;
	else if (IsUpperTriangular())
		out << "  Matrix is upper-triangular" << endl;
	else
		out << "  Matrix is rectangular" << endl;

	if (IsInterleave())
		out << "  Matrix is interleaved" << endl;
	else 
		out << "  Matrix is non-interleaved" << endl;

	if (IsLabels())
		out << "  Taxon labels provided" << endl;
	else
		out << "  No taxon labels provided" << endl;

	if (IsDiagonal())
		out << "  Diagonal elements specified" << endl;
	else 
		out << "  Diagonal elements not specified" << endl;

	out << "  Missing data symbol is " << missing << endl;

	if (ntax == 0)
		return;

	out.setf(ios::fixed, ios::floatfield);
	out.setf(ios::showpoint);
	for (unsigned i = 0; i < ntax; i++)
		{
		if (labels)
			out << setw(20) << taxa->GetTaxonLabel(i);
		else
			out << "\t\t";

		for (unsigned j = 0; j < ntax; j++)
			{
			if (triangle == NxsDistancesBlockEnum(upper) && j < i)
				{
				out << setw(12) << " ";
				}
			else if (triangle == NxsDistancesBlockEnum(lower) && j > i)
				continue;
			else if (!diagonal && i == j)
				{
				out << setw(12) << " ";
				}
			else if (IsMissing(i, j))
				out << setw(12) << missing;
			else
				out << setw(12) << setprecision(5) << GetDistance(i, j);
			}

		out << endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Flushes taxonLabels and sets ntax to 0 in preparation for reading a new TAXA block.
*/
void NxsDistancesBlock::Reset()
	{
	// Reset base class data members that could have changed
	//
	errormsg.clear();
	isEnabled      = true;
	isEmpty        = true;
	isUserSupplied = false;

	if (matrix != NULL)
		{
		for (unsigned i = 0; i < ntax; i++)
			delete [] matrix[i];
		delete [] matrix;
		matrix = NULL;
		}

	if (taxonPos != NULL)
		delete [] taxonPos;
	taxonPos = NULL;

	ntax        = 0;
	nchar       = 0;
	diagonal    = true;
	labels      = true;
	newtaxa     = false;
	interleave  = false;
	missing     = '?';
	triangle    = NxsDistancesBlockEnum(lower);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of ntax.
*/
unsigned NxsDistancesBlock::GetNtax()
	{
	return ntax;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of nchar.
*/
unsigned NxsDistancesBlock::GetNchar()
	{
	return nchar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the (`i', `j')th element of `matrix'. Assumes `i' and `j' are both in the range [0..`ntax') 
|	and the distance stored at `matrix[i][j]' is not missing. Also assumes `matrix' is not NULL.
*/
double NxsDistancesBlock::GetDistance(
  unsigned i,	/* the row */
  unsigned j)	/* the column */
	{
	assert(i >= 0);
	assert(i < ntax);
	assert(j >= 0);
	assert(j < ntax);
	assert(matrix != NULL);

	return matrix[i][j].value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `missing'.
*/
char NxsDistancesBlock::GetMissingSymbol()
	{
	return missing;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `triangle'.
*/
unsigned NxsDistancesBlock::GetTriangle()
	{
	return triangle;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the value of `triangle' is NxsDistancesBlockEnum(both), false otherwise.
*/
bool NxsDistancesBlock::IsRectangular()
	{
	return (triangle == NxsDistancesBlockEnum(both) ? true : false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the value of triangle is NxsDistancesBlockEnum(upper), false otherwise.
*/
bool NxsDistancesBlock::IsUpperTriangular()
	{
	return (triangle == NxsDistancesBlockEnum(upper) ? true : false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the value of triangle is NxsDistancesBlockEnum(lower), false otherwise.
*/
bool NxsDistancesBlock::IsLowerTriangular()
	{
	return (triangle == NxsDistancesBlockEnum(lower) ? true : false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of diagonal.
*/
bool NxsDistancesBlock::IsDiagonal()
	{
	return diagonal;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of interleave.
*/
bool NxsDistancesBlock::IsInterleave()
	{
	return interleave;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of labels.
*/
bool NxsDistancesBlock::IsLabels()
	{
	return labels;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the (`i',`j')th distance is missing. Assumes `i' and `j' are both in the range [0..`ntax') and 
|	`matrix' is not NULL.
*/
bool NxsDistancesBlock::IsMissing(
  unsigned i,	/* the row */
  unsigned j)	/* the column */
	{
	assert(i >= 0);
	assert(i < ntax);
	assert(j >= 0);
	assert(j < ntax);
	assert(matrix != NULL);

	return (bool)(matrix[i][j].missing);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the (`i',`j')th matrix element to `d' and `missing' to false . Assumes `i' and `j' are both in 
|	the range [0..`ntax') and `matrix' is not NULL.
*/
void NxsDistancesBlock::SetDistance(
  unsigned i,	/* the row */
  unsigned j,	/* the column */
  double d)		/* the distance value */
	{
	assert(i >= 0);
	assert(i < ntax);
	assert(j >= 0);
	assert(j < ntax);
	assert(matrix != NULL);

	matrix[i][j].value = d;
	matrix[i][j].missing = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the (`i', `j')th `matrix' element to missing. Assumes `i' and `j' are both in the range 
|	[0..`ntax') and `matrix' is not NULL.
*/
void NxsDistancesBlock::SetMissing(
  unsigned i,	/* the row */
  unsigned j)	/* the column */
	{
	assert(i >= 0);
	assert(i < ntax);
	assert(j >= 0);
	assert(j < ntax);
	assert(matrix != NULL);

	matrix[i][j].missing = 1;
	matrix[i][j].value = 0.0;
	}

 /*----------------------------------------------------------------------------------------------------------------------
|	Sets `nchar' to `n'.
*/
void NxsDistancesBlock::SetNchar(
  unsigned n)	/* the number of characters */
	{
	nchar = n;
	}
