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
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `id' to "TREES", `ntrees' to 0, `defaultTree' to 0, and `taxa' to `tb'. Assumes `tb' is non-NULL.
*/
NxsTreesBlock::NxsTreesBlock(
  NxsTaxaBlock *tb)	/* the NxsTaxaBlock object to be queried for taxon names appearing in tree descriptions */
  : NxsBlock(), taxa(tb)
	{
	assert(tb != NULL);

	id			= "TREES";
	taxa		= tb;
	ntrees		= 0;
	defaultTree	= 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears `translateList', `rooted', `treeName' and `treeDescription'.
*/
NxsTreesBlock::~NxsTreesBlock()
	{
	translateList.clear();
	rooted.clear();
	treeName.clear();
	treeDescription.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes data member `taxa' point to `tb' rather than the NxsTaxaBlock object it was previously pointing to. Assumes 
|	`tb' is non-NULL.
*/
void NxsTreesBlock::ReplaceTaxaBlockPtr(
  NxsTaxaBlock *tb)		/* pointer to new NxsTaxaBlock object (does not attempt to delete the object previously pointed to) */
	{
	assert(tb != NULL);

	taxa = tb;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Takes control from the Read member function when a TREE or UTREE command is encountered. If a TREE command is found,
|	the HandleTreeDescription member function is called with `utree' equal to false. If a UTREE command is found, 
|	`utree' equals true.
*/
void NxsTreesBlock::HandleTreeDescription(
  NxsToken &token,	/* the token used to read from `in' */
  bool utree)		/* true if handling UTREE command, false if handling TREE command */
	{
	// Start off assuming that there will be no command comments contradicting
	// the rooted/unrooted status impied by the use of the TREE/UTREE command
	//
	bool tree_is_unrooted = utree;

	// This should be either an asterisk or a tree name
	//
	token.GetNextToken();

	if (token.Equals("*"))
		{
		// ntrees is not incremented until entire tree command has been read
		//
		defaultTree = ntrees; 

		// This should be tree name
		//
		token.GetNextToken();
		}

	// Save the tree name as the key
	//
	NxsString skey = token.GetToken();

	// This should be an equals sign
	//
	token.GetNextToken();

	if (!token.Equals("="))
		{
		errormsg = "Expecting '=' after tree name in TREE command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	// This should be either a tree description or a command comment specifying
	// whether this tree is to be rooted ([&R]) or unrooted ([&U]).
	//
	token.SetLabileFlagBit(NxsToken::saveCommandComments);
	token.SetLabileFlagBit(NxsToken::parentheticalToken);
	token.GetNextToken();

	NxsString s = token.GetToken();
	NxsString cmdName = (utree ? "UTREE" : "TREE");
	if (s.size() < 2)
		{
		errormsg = "Expecting command comment or tree description in ";
		errormsg += cmdName;
		errormsg += " command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	if (s[0] == '&') 
		{
		// Command comment found
		//
		if (s[1] == 'R' || s[1] == 'r')
			{
			tree_is_unrooted = false;
			}
		else if (s[1] == 'U' || s[1] == 'u')
			{
			tree_is_unrooted = true;
			}
		else
			{
			errormsg = "[";
			errormsg += token.GetToken();
			errormsg += "] is not a valid command comment in a ";
			errormsg += cmdName;
			errormsg += " command";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

		// This should be only the tree description
		//
		token.SetLabileFlagBit(NxsToken::parentheticalToken);
		token.GetNextToken();
		}

	NxsString sval = token.GetToken();

	// This should be a semicolon
	//
	token.GetNextToken();

	if (!token.Equals(";"))
		{
		errormsg = "Expecting ';' to terminate the ";
		errormsg += cmdName;
		errormsg += " command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	ntrees++;
	treeName.push_back(skey);
	treeDescription.push_back(sval);

	if (tree_is_unrooted)
		rooted.push_back(false);
	else
		rooted.push_back(true);

	assert(rooted.size() == (unsigned)ntrees);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the ability to read everything following the block name (which is read by the NxsReader 
|	object) to the END or ENDBLOCK command. Characters are read from the input stream `in'. Overrides the abstract 
|	virtual function in the base class.
*/
void NxsTreesBlock::Read(
  NxsToken &token)	/* the token used to read from `in' */
	{
	isEmpty = false;
	isUserSupplied = true;

	// This should be the semicolon after the block name
	//
	token.GetNextToken();

	if (!token.Equals(";"))
		{
		errormsg.PrintF("Expecting ';' after TREES block name, but found %s instead.", token.GetTokenAsCStr());
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("TRANSLATE")) 
			{
			// Note that numEntries will be 0 if no taxa block
			// has been created
			//
			unsigned numEntries = taxa->GetNumTaxonLabels();
			bool building_taxa_block = (numEntries == 0);

			for (unsigned k = 0; ; k++) 
				{
				if (numEntries > 0 && k == numEntries)
					break;

				// Create the Association

				// Get the key
				//
				token.GetNextToken();
				NxsString skey = token.GetToken();

				// Get the value
				//
				token.GetNextToken();
				NxsString sval = token.GetToken();

				// Add the taxon label to the TAXA block (if we are building one)
				// and check to make sure taxon label is one that is in the taxa
				// block (if we are not building up the taxa block as we go)
				//
				if (building_taxa_block)
					taxa->AddTaxonLabel(sval);
				else
					{
					try
						{
						taxa->FindTaxon(sval);
						}
					catch(NxsTaxaBlock::NxsX_NoSuchTaxon)
						{
						errormsg.clear();
						errormsg.PrintF("The taxon \"%s\" was found in the TRANSLATE command of the TREES block, but was not found in the TAXA block", sval.c_str());
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}
					}

				// Add the Association object to the translate list
				//
				translateList[skey] = sval;

				// This should be a comma, unless we are at the last pair, in
				// which case it should be a semicolon. If it is a semicolon,
				// and a TAXA block exists, we should have read in exactly
				// numEntries taxon labels at this point
				//
				token.GetNextToken();

				if (token.Equals(";")) 
					{
					if (numEntries > 0 && k != numEntries - 1)
						{
						errormsg.clear();
						errormsg.PrintF("There were %d entries in TRANSLATE command but only %d taxa in the TAXA block.", k + 1, numEntries);
						errormsg += "\nThe number of TRANSLATE entries should equal the number of taxa.";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}
					break;
					}

				else if (!token.Equals(","))
					{
					errormsg.clear();
					errormsg.PrintF("Expecting ',' to terminate each number/name pair in TRANSLATE command, but found %s instead.", token.GetTokenAsCStr());
					errormsg += "\nPerhaps there were fewer taxa in the tree file than previously defined.";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
				}	// for (unsigned k = 0; ; k++) 
			}	// if (token.Equals("TRANSLATE")) 

		else if (token.Equals("TREE")) 
			{
			HandleTreeDescription(token, false);
			}	

		else if (token.Equals("UTREE")) 
			{
			HandleTreeDescription(token, true);
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
			}	// else if (token.Equals("END")) 

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
			}	// else if (token.Equals("ENDBLOCK")) 

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
		}	// for (;;)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Flushes `treeName', `treeDescription', `translateList' and `rooted', and sets `ntrees' and `defaultTree' both to 0
|	in preparation for reading a new TREES block.
*/
void NxsTreesBlock::Reset()
	{
	// Reset base class data members that could have changed
	//
	errormsg.clear();
	isEnabled      = true;
	isEmpty        = true;
	isUserSupplied = false;

	ntrees			= 0;
	defaultTree		= 0;

	treeName.clear();
	treeDescription.clear();
	translateList.clear();
	rooted.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the 0-offset index of the default tree, which will be 0 if there is only one tree stored or no trees 
|	stored. If more than one tree is stored, the default tree will be the one specifically indicated by the user (using
|	an asterisk in the data file), or 0 if the user failed to specify.
*/
unsigned NxsTreesBlock::GetNumDefaultTree()
	{
	return defaultTree;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of trees stored in this NxsTreesBlock object.
*/
unsigned NxsTreesBlock::GetNumTrees()
	{
	return ntrees;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the name of the tree stored at position `i' in `treeName'. Assumes that `i' will be in the range 
|	[0..`ntrees').
*/
NxsString NxsTreesBlock::GetTreeName(
  unsigned i)	/* the index of the tree for which the name is to be returned */
	{
	assert(i >= 0);
	assert(i < ntrees);

	return treeName[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the description of the tree stored at position `i' in `treeName'. Assumes that `i' will be in the range 
|	[0..ntrees). Node numbers will be translated to names in the resulting tree description. Use GetTreeDescription if 
|	translation is not desired.
*/
NxsString NxsTreesBlock::GetTranslatedTreeDescription(
  unsigned i)	/* the index of the tree for which the description is to be returned */
	{
	assert(i >= 0);
	assert(i < ntrees);

	//bool adding_labels = (taxa->GetNumTaxonLabels() == 0);

	// s is the original tree definition string 
	//
	NxsString s = treeDescription[i];
	unsigned slen = s.size();
	assert(slen > 1);

	// x is the new tree definition string that will be built
	// using s as the template
	//
	NxsString x;
	x += s[0];

	//bool inside_tip_label = false;
	for (unsigned k = 1; k < slen; k++)
		{
		char prev = s[k - 1];
		char curr = s[k];

		if (isdigit(curr) && (prev == '(' || prev == ','))
			{
			// Discovered a number where a taxon label should be in the tree description
			// Read entire number and then look it up in the translateList
			//
			NxsString ns;
			ns += curr;
			for (;;)
				{
				curr = s[k+1];
				prev = s[k++];
				if (isdigit(curr))
					ns += curr;
				else
					{
					--k;
					break;
					}
				}
			NxsString nss = translateList[ns];
			x += nss;
			}
		else
			x += curr;
		}

	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the description of the tree stored at position `i' in `treeDescription'. Assumes that `i' will be in the 
|	range [0..`ntrees').
*/
NxsString NxsTreesBlock::GetTreeDescription(
  unsigned i)	/* the index of the tree for which the description is to be returned */
	{
	assert(i >= 0);
	assert(i < ntrees);

	return treeDescription[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the `i'th tree (0-offset) is the default tree, false otherwise. Assumes that `i' will be in the 
|	range [0..ntrees).
*/
bool NxsTreesBlock::IsDefaultTree(
  unsigned i)	/* the index of the tree in question */
	{
	assert(i >= 0);
	assert(i < ntrees);

	if (i == GetNumDefaultTree())
		return true;
	else
		return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the `i'th tree (0-offset) is rooted, false otherwise. Assumes that `i' will be in the 
|	range [0..ntrees).
*/
bool NxsTreesBlock::IsRootedTree(
  unsigned i)	/* the index of the tree in question */
	{
	assert(i >= 0);
	assert(i < ntrees);

	return rooted[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function outputs a brief report of the contents of this block. Overrides the abstract virtual function in the 
|	base class.
*/
void NxsTreesBlock::Report(
  ostream &out)	/* the output stream to which to write the report */
	{
	out << endl;
	out << id << " block contains ";
	if (ntrees == 0)
		{
		out << "no trees" << endl;
		}
	else if (ntrees == 1)
		out << "one tree" << endl;
	else
		out << ntrees << " trees" << endl;

	if (ntrees == 0)
		return;

	for (unsigned k = 0; k < ntrees; k++)
		{
		out << '\t' << (k+1) << '\t' << treeName[k];
		out << "\t(";
		if (rooted[k])
			out << "rooted";
		else
			out << "unrooted";
		if (defaultTree == k)
			out << ",default tree)" << endl;
		else
			out << ')' << endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Outputs a brief description of this block's contents to the referenced NxsString. An example of the output of this 
|	command is shown below:
|>
|	TREES block contains 102 trees
|>
*/
void NxsTreesBlock::BriefReport(
  NxsString &s)	/* reference to the string in which to store the contents of the brief report */
	{
	s = "\n\n";
	s += id;
	s += " block contains ";
	if (ntrees == 0)
		s += "no trees\n";
	else if (ntrees == 1)
		s += "one tree\n";
	else
		{
		s += ntrees;
		s += " trees\n";
		}
	}
