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
|	Sets id = "ASSUMPTIONS", charBlockPtr = NULL, and taxa = t. Assumes taxa is non-NULL.
*/
NxsAssumptionsBlock::NxsAssumptionsBlock(
  NxsTaxaBlock *t)	/* pointer to the taxa block */
	{
	assert(t);
	taxa			= t;
	charBlockPtr	= NULL;
	id				= "ASSUMPTIONS";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing needs to be done in the destructor.
*/
NxsAssumptionsBlock::~NxsAssumptionsBlock()
{
}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes data member taxa point to 'tb'. Assumes tb is non-NULL.
*/
void NxsAssumptionsBlock::ReplaceTaxaBlockPtr(
  NxsTaxaBlock *tb)	/* pointer to new NxsTaxaBlock object */
	{
	assert(tb);
	taxa = tb;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of character sets stored.
*/
int NxsAssumptionsBlock::GetNumCharSets()
	{
	return (int)charsets.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Erases 'names' vector, then fills 'names' with the names of all stored character sets.
*/
void NxsAssumptionsBlock::GetCharSetNames(
  NxsStringVector &names)	/* the vector in which to store the names */
	{
	names.erase(names.begin(), names.end());
	NxsUnsignedSetMap::const_iterator i;
	for (i = charsets.begin(); i != charsets.end(); i++)
	names.push_back((*i).first);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns reference to character set having name 'nm'.
*/
NxsUnsignedSet &NxsAssumptionsBlock::GetCharSet(
  NxsString nm)	/* the name of the character set to return */
	{
	return charsets[nm];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns name of default character set. If returned string has zero length, then no default character set was defined
|	in the data set.
*/
NxsString NxsAssumptionsBlock::GetDefCharSetName()
	{
	return def_charset;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of taxon sets stored.
*/
int NxsAssumptionsBlock::GetNumTaxSets()
	{
	return (int)taxsets.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Erases 'names' vector, then fills 'names' with the names of all stored taxon sets.
*/
void NxsAssumptionsBlock::GetTaxSetNames(
  NxsStringVector &names)	/* the vector in which to store the names */
	{
	names.erase(names.begin(), names.end());
	NxsUnsignedSetMap::const_iterator i;
	for (i = taxsets.begin(); i != taxsets.end(); i++)
		names.push_back((*i).first);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns reference to taxon set having name 'nm'.
*/
NxsUnsignedSet &NxsAssumptionsBlock::GetTaxSet(
  NxsString nm)	/* the name of the taxon set to return */
	{
	return taxsets[nm];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns name of default taxon set. If returned string has zero length, then no default taxon set was defined in the
|	data set.
*/
NxsString NxsAssumptionsBlock::GetDefTaxSetName()
	{
	return def_taxset;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of exclusion sets stored.
*/
int NxsAssumptionsBlock::GetNumExSets()
	{
	return (int)exsets.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Erases names, then fills names with the names of all stored exclusion sets.
*/
void NxsAssumptionsBlock::GetExSetNames(
  NxsStringVector &names)	/* the vector in which to store the names */
	{
	names.erase(names.begin(), names.end());
	NxsUnsignedSetMap::const_iterator i;
	for (i = exsets.begin(); i != exsets.end(); i++)
		names.push_back((*i).first);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns reference to exclusion set having name 'nm'.
*/
NxsUnsignedSet &NxsAssumptionsBlock::GetExSet(
  NxsString nm)	/* the name of the exclusion set to return */
	{
	return exsets[nm];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns name of default exclusion set. If returned string has zero length, then no default exclusion set was defined
|	in the data set.
*/
NxsString NxsAssumptionsBlock::GetDefExSetName()
	{
	return def_exset;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Applies exclusion set having name 'nm' by calling the ApplyExset method of the NxsCharactersBlock or 
|	NxsCharactersBlock-derived object stored in the charBlockPtr pointer (which will be whichever block last called the 
|	NxsAssumptionsBlock::SetCallback method).
*/
void NxsAssumptionsBlock::ApplyExSet(
  NxsString nm)	/* the name of the exclusion set to apply */
	{
	assert(charBlockPtr);
	charBlockPtr->ApplyExset(exsets[nm]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads and stores information contained in the command CHARSET within an ASSUMPTIONS block.
*/
void NxsAssumptionsBlock::HandleCharset(
  NxsToken &token)	/* the token used to read from in */
	{
	bool asterisked = false;

	// Next token should be either an asterisk or the name of a charset
	//
	token.GetNextToken();

	if (token.Equals("*"))
		{
		asterisked = true;
		token.GetNextToken();
		}

	// Token now stored should be the name of a charset
	//
	NxsString charset_name = token.GetToken();

	// Now grab the equals sign
	//
	token.GetNextToken();
	if (!token.Equals("="))
		{
		errormsg = "Expecting '=' in CHARSET definition but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	assert(charBlockPtr);
	NxsCharactersBlock &charBlock = *charBlockPtr;
	NxsUnsignedSet s;
	int totalChars = charBlock.GetNCharTotal();
	NxsSetReader(token, totalChars, s, charBlock, NxsSetReader::charset).Run();

	charsets[charset_name] = s;

	if (asterisked)
		def_charset = charset_name;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when the END or ENDBLOCK command needs to be parsed from within the ASSUMPTIONS block. Basically just checks
|	to make sure the next token in the data file is a semicolon.
*/
void NxsAssumptionsBlock::HandleEndblock(
  NxsToken &token)	/* the token used to read from in */
	{
	// Get the semicolon following END or ENDBLOCK token
	//
	token.GetNextToken();

	if (!token.Equals(";"))
		{
		errormsg = "Expecting ';' to terminate the END or ENDBLOCK command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads and stores information contained in the command EXSET within an ASSUMPTIONS block. If EXSET keyword is 
|	followed by an asterisk, last read NxsCharactersBlock or NxsCharactersBlock-derived object is notified of the 
|	characters to be excluded (its ApplyExset function is called).
*/
void NxsAssumptionsBlock::HandleExset(
  NxsToken &token)	/* the token used to read from in */
	{
	bool asterisked = false;

	// Next token should be either an asterisk or the name of an exset
	//
	token.GetNextToken();

	if (token.Equals("*"))
		{
		asterisked = true;
		token.GetNextToken();
		}

	// Token now stored should be the name of an exset
	//
	NxsString exset_name = token.GetToken();

	// Now grab the equals sign
	//
	token.GetNextToken();
	if (!token.Equals("="))
		{
		errormsg = "Expecting '=' in EXSET definition but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	assert(charBlockPtr);
	NxsCharactersBlock &charBlock = *charBlockPtr;
	NxsUnsignedSet s;
	int totalChars = charBlock.GetNCharTotal();
	NxsSetReader(token, totalChars, s, charBlock, NxsSetReader::charset).Run();

	exsets[exset_name] = s;

	if (asterisked)
		{
		def_exset = exset_name;
		charBlock.ApplyExset(s);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads and stores information contained in the command TAXSET within an ASSUMPTIONS block.
*/
void NxsAssumptionsBlock::HandleTaxset(
  NxsToken &token)	/* the token used to read from in */
	{
	bool asterisked = false;

	// Next token should be either an asterisk or the name of a taxset
	//
	token.GetNextToken();

	if (token.Equals("*"))
		{
		asterisked = true;
		token.GetNextToken();
		}

	// Token now stored should be the name of a taxset
	//
	NxsString taxset_name = token.GetToken();

	// Now grab the equals sign
	//
	token.GetNextToken();
	if (!token.Equals("="))
		{
		errormsg = "Expecting '=' in TAXSET definition but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	NxsUnsignedSet s;
	int totalTaxa = taxa->GetNumTaxonLabels();
	NxsSetReader(token, totalTaxa, s, *this, NxsSetReader::taxset).Run();

	taxsets[taxset_name] = s;

	if (asterisked)
		def_taxset = taxset_name;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the ability to read everything following the block name (which is read by the NxsReader 
|	object) to the end or ENDBLOCK statement. Characters are read from the input stream in. Overrides the pure virtual
|	function in the base class.
*/
void NxsAssumptionsBlock::Read(
  NxsToken &token)	/* the token used to read from in */
	{
	isEmpty = false;
	isUserSupplied = true;

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

	for(;;)
		{
		token.GetNextToken();

		if (token.Equals("EXSET"))
			{
			HandleExset(token);
			}
		else if (token.Equals("TAXSET"))
			{
			HandleTaxset(token);
			}
		else if (token.Equals("CHARSET"))
			{
			HandleCharset(token);
			}
		else if (token.Equals("END"))
			{
			HandleEndblock(token);
			break;
			}
		else if (token.Equals("ENDBLOCK"))
			{
			HandleEndblock(token);
			break;
			}
		else
			{
			SkippingCommand(token.GetToken());
			do
				{
				token.GetNextToken();
				} while(!token.AtEOF() && !token.Equals(";"));

			if (token.AtEOF())
				{
				errormsg = "Unexpected end of file encountered";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			}
		}	// for(;;)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Prepares for reading a new ASSUMPTIONS block. Overrides the pure virtual function in the base class.
*/
void NxsAssumptionsBlock::Reset()
	{
	exsets.erase(exsets.begin(), exsets.end());
	taxsets.erase(taxsets.begin(), taxsets.end());
	charsets.erase(charsets.begin(), charsets.end());
	def_taxset.clear();
	def_charset.clear();
	def_exset.clear();
	errormsg.clear();
	isEnabled		= true;
	isEmpty			= true;
	isUserSupplied	= false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function outputs a brief report of the contents of this ASSUMPTIONS block. Overrides the pure virtual function
|	in the base class.
*/
void NxsAssumptionsBlock::Report(
  ostream &out)	/* the output stream to which to write the report */
	{
	out << endl;
	out << id << " block contains the following:" << endl;

	if (charsets.empty())
		out << "  No character sets were defined" << endl;
	else
		{
		NxsUnsignedSetMap::const_iterator charsets_iter = charsets.begin();
		if (charsets.size() == 1)
			{
			out << "  1 character set defined:" << endl;
			out << "    " << (*charsets_iter).first << endl;
			}
		else
			{
			out << "  " << charsets.size() << " character sets defined:" << endl;
			for (; charsets_iter != charsets.end(); charsets_iter++)
				{
				NxsString nm = (*charsets_iter).first;
				out << "    " << nm;
				if (nm == def_charset)
					out << " (default)";
				out << endl;
				}
			}
		}	// if (charsets.empty()) ... else

	if (taxsets.empty())
		out << "  No taxon sets were defined" << endl;
	else
		{
		NxsUnsignedSetMap::const_iterator taxsets_iter = taxsets.begin();
		if (taxsets.size() == 1)
			{
			out << "  1 taxon set defined:" << endl;
			out << "    " << (*taxsets_iter).first << endl;
			}
		else
			{
			out << "  " << taxsets.size() << " taxon sets defined:" << endl;
			for (; taxsets_iter != taxsets.end(); taxsets_iter++)
				{
				NxsString nm = (*taxsets_iter).first;
				out << "    " << nm;
				if (nm == def_taxset)
					out << " (default)";
				out << endl;
				}
			}
		}	// if (taxsets.empty()) ... else

	if (exsets.empty())
		out << "  No exclusion sets were defined" << endl;
	else
		{
		NxsUnsignedSetMap::const_iterator exsets_iter = exsets.begin();
		if (exsets.size() == 1)
			{
			out << "  1 exclusion set defined:" << endl;
			out << "    " << (*exsets_iter).first << endl;
			}
		else
			{
			out << "  " << exsets.size() << " exclusion sets defined:" << endl;
			for (; exsets_iter != exsets.end(); exsets_iter++)
				{
				NxsString nm = (*exsets_iter).first;
				out << "    " << nm;
				if (nm == def_exset)
				out << " (default)";
				out << endl;
				}
			}
		}

	out << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	A CHARACTERS, DATA, or ALLELES block can call this function to specify that it is to receive notification when the 
|	current taxon or character set changes (e.g., an "EXSET *" command is read or a program requests that one of the 
|	predefined taxon sets, character sets, or exsets be applied). Normally, a NxsCharactersBlock-derived object calls 
|	this function upon entering its MATRIX command, since when that happens it becomes the primary data-containing block.
*/
void NxsAssumptionsBlock::SetCallback(
  NxsCharactersBlock* p)	/* the object to be called in the event of a change in character status */
	{
	charBlockPtr = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts a taxon label to a number corresponding to the taxon's position within the list maintained by the 
|	NxsTaxaBlock object. This method overrides the virtual function of the same name in the NxsBlock base class. If s 
|	is not a valid taxon label, returns the value 0.
*/
unsigned NxsAssumptionsBlock::TaxonLabelToNumber(
  NxsString s)	/* the taxon label to convert */
	{
	int i;
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
