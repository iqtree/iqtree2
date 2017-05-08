/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "msplitsblock.h"
#include "pda/split.h"
#include "pda/splitgraph.h"

MSplitsBlock::MSplitsBlock(SplitGraph *asgraph)
		: NxsBlock()
{
	nsplits = 0;
	ntaxa = 0;
	id = "SPLITS";
	sgraph = asgraph;
}


void MSplitsBlock::Report(ostream &out)
{
	sgraph->report(out);
}

void MSplitsBlock::Reset()
{
	errormsg.clear();
	isEmpty			= true;
	isEnabled		= true;
	isUserSupplied	= false;

	ntaxa			= 0;
	nsplits = 0;
	sgraph->clear();
}

void MSplitsBlock::AddSplit(NxsToken &token)
{
	token.SetLabileFlagBit(token.hyphenNotPunctuation);
	// this should be the weight of split
	token.GetNextToken();

	NxsString str = token.GetToken();

	// get the next token
	token.GetNextToken();

	// if continuing number....
	if (token.GetToken() == "+")
	{
		str += token.GetToken();
		token.GetNextToken();
		str += token.GetToken();
		token.GetNextToken();
		//cout << str << " ";
	}

	double weight = atof(str.c_str());

	//cout << token.GetToken() << " split weight = " << weight << endl;

	vector<int> taxa_list;

	//@pol should check to make sure this is not punctuation
	while (!token.AtEOF() && !token.Equals(","))
	{
		int index = atoi(token.GetToken().c_str());
		if (index < 1 || index > ntaxa)
		{
			errormsg = "Taxon index should be greater than zero and less than ";
			errormsg += (ntaxa+1);
			errormsg += "(";
			errormsg += token.GetToken();
			errormsg += " was specified)";
			throw NxsException(errormsg, token);
		}
		
		taxa_list.push_back(index-1);

		token.GetNextToken();
	}

	if (token.AtEOF())
	{
		errormsg = "Unexpected end of file encountered";
		throw NxsException(errormsg, token);
	}

	Split *split = new Split(ntaxa, weight, taxa_list);	
	sgraph->push_back(split);
}

void MSplitsBlock::Read(NxsToken &token)
{

	//int nominal_ntax	= 0;
	//int nominal_nsplits = 0;

	// This should be the semicolon after the block name
	//
	token.GetNextToken();

	if (!token.Equals(";"))
	{
		errormsg = "Expecting ';' after TAXA block name, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token);
	}

	for (;;)
	{
		token.GetNextToken();

		if (token.Equals("DIMENSIONS"))
		{
			// This should be the NTAX keyword
			//
			token.GetNextToken();

			if (!token.Equals("NTAX"))
			{
				errormsg = "Expecting NTAX keyword, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}

			// This should be the equals sign
			//
			token.GetNextToken();

			if (!token.Equals("="))
			{
				errormsg = "Expecting '=', but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}

			// This should be the number of taxa
			//
			token.GetNextToken();

			ntaxa = atoi(token.GetToken().c_str());
			if (ntaxa <= 0)
			{
				errormsg = "NTAX should be greater than zero (";
				errormsg += token.GetToken();
				errormsg += " was specified)";
				throw NxsException(errormsg, token);
			}

			// This should be the NSPLITS keyword
			//
			token.GetNextToken();

			if (!token.Equals("NSPLITS"))
			{
				errormsg = "Expecting NSPLITS keyword, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}

			// This should be the equals sign
			//
			token.GetNextToken();

			if (!token.Equals("="))
			{
				errormsg = "Expecting '=', but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}

			// This should be the number of taxa
			//
			token.GetNextToken();

			nsplits = atoi(token.GetToken().c_str());
			if (nsplits <= 0)
			{
				errormsg = "NSPLITS should be greater than zero (";
				errormsg += token.GetToken();
				errormsg += " was specified)";
				throw NxsException(errormsg, token);
			}


			// This should be the terminating semicolon
			//
			token.GetNextToken();

			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate DIMENSIONS command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}
		}	// if (token.Equals("DIMENSIONS"))

		else if (token.Equals("CYCLE"))
		{
			if (ntaxa <= 0)
			{
				errormsg = "DIMENSIONS must be specified before CYCLE command";
				throw NxsException(errormsg, token);
			}
			token.GetNextToken();
			while (!token.AtEOF() && !token.Equals(";")) {
				int tax = atoi(token.GetToken().c_str());
				if (tax <= 0 || tax > ntaxa)
				{
					errormsg = "taxon index in CYCLE should be between 1 and";
					errormsg += ntaxa;
					errormsg += " (";
					errormsg += token.GetToken();
					errormsg += " was specified)";
					throw NxsException(errormsg, token);
				}
				cycle.push_back(tax-1);
				token.GetNextToken();
			}
			if (cycle.size() != ntaxa) {
				errormsg = "Not all taxa in CYCLE are included";
				throw NxsException(errormsg, token);
			}
			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate CYCLE command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}
		}
		else if (token.Equals("MATRIX"))
		{
			if (nsplits <= 0)
			{
				errormsg = "NSPLITS must be specified before MATRIX command";
				throw NxsException(errormsg, token);
			}

			for (int i = 0; i < nsplits; i++)
			{
				AddSplit(token);
			}

			// This should be terminating semicolon
			//
			token.GetNextToken();

			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate MATRIX command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}
		}	// if (token.Equals("MATRIX"))

		else if (token.Equals("END") || token.Equals("ENDBLOCK"))
		{
			// Get the semicolon following END
			//
			token.GetNextToken();

			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate the ENDBLOCK command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token);
			}
			break;
		}	// if (token.Equals("END") || token.Equals("ENDBLOCK"))

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
				throw NxsException(errormsg, token);
			}
		}	// token not END, ENDBLOCK, MATRIX, or DIMENSIONS
	}	// GetNextToken loop

}


MSplitsBlock::~MSplitsBlock()
{}


