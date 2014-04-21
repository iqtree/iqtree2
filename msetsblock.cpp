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
#include "msetsblock.h"

MSetsBlock::MSetsBlock()
 : NxsBlock()
{
	id = "SETS";
}


MSetsBlock::~MSetsBlock()
{
	for (TaxaSetNameVector::reverse_iterator it = sets.rbegin(); it != sets.rend(); it++) {
		//cout << (*it)->name << endl;
		delete *it;
	}
	sets.clear();
}


void MSetsBlock::Report(ostream &out)
{
	int nsets = getNSets();
	out << "Number of sets: " << nsets << endl;
	for (TaxaSetNameVector::iterator i = sets.begin(); i != sets.end(); i++) {
		out << "Set " << (*i)->name << " contains: ";
		for (vector<string>::iterator it = (*i)->taxlist.begin(); it != (*i)->taxlist.end(); it++)
			out << (*it) << "  ";
		out << endl;
	}
}

void MSetsBlock::Reset()
{
	for (TaxaSetNameVector::reverse_iterator it = sets.rbegin(); it != sets.rend(); it++)
		delete *it;
	sets.clear();
}

void MSetsBlock::Read(NxsToken &token)
{
	// This should be the semicolon after the block name
	//

	token.GetNextToken();

	if (!token.Equals(";"))
	{
		errormsg = "Expecting ';' after SETS block name, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	for (;;)
	{
		token.GetNextToken();

		if (token.Equals("TAXSET"))
		{
			// This should be the NTAX keyword
			//
			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextToken();

			//sets.resize(sets.size()+1);
			TaxaSetName *myset = new TaxaSetName;
			sets.push_back(myset);

			myset->name = token.GetToken();

			token.GetNextToken();

			if (!token.Equals("="))
			{
				errormsg = "Expecting '=', but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextToken();
			do {
				myset->taxlist.push_back(token.GetToken());
				token.SetLabileFlagBit(NxsToken::preserveUnderscores);
				token.GetNextToken();
			} while (!token.AtEOF() && !token.Equals(";"));

			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate PARAMETERS command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

		} // if (token.Equals("TAXSET"))

		else if (token.Equals("CHARSET"))
		{
			// This should be the NTAX keyword
			//
			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextToken();


			//sets.resize(sets.size()+1);
			CharSet *myset = new CharSet;
			charsets.push_back(myset);
			myset->aln_file = "";
			myset->model_name = "";
			myset->position_spec = "";
			myset->sequence_type = "";
			myset->char_partition = "";

			myset->name = token.GetToken();

			token.GetNextToken();
			if (!token.Equals("="))
			{
				errormsg = "Expecting '=', but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextContiguousToken(';');
			myset->position_spec = token.GetToken();

			// separate position_spec into alignment name if ':' exists
			size_t pos = myset->position_spec.find(':');
			if (pos != string::npos) {
				myset->aln_file = myset->position_spec.substr(0, pos);
				myset->position_spec = myset->position_spec.substr(pos+1);
			}

			token.GetNextToken();
			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate PARAMETERS command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

		} // if (token.Equals("CHARSET"))
		else if (token.Equals("CHARPARTITION"))
		{
			// This should be the NTAX keyword
			//
			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextToken();
			string partition_name = token.GetToken();
			token.GetNextToken();
			if (!token.Equals("="))
			{
				errormsg = "Expecting '=', but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}
			token.SetLabileFlagBit(NxsToken::preserveUnderscores);
			token.GetNextToken();
			do {
				string model_name = "";
				while (!token.AtEOF() && !token.Equals(":")) {
					model_name += token.GetToken();
					token.SetLabileFlagBit(NxsToken::preserveUnderscores);
					token.GetNextToken();
				}

				if (!token.Equals(":"))
				{
					errormsg = "Expecting ':' or ',' but found ";
					errormsg += token.GetToken();
					errormsg += " instead";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
				string charset_name;
				token.SetLabileFlagBit(NxsToken::preserveUnderscores);
				token.GetNextToken();
				charset_name = token.GetToken();
				CharSet *myset = findCharSet(charset_name);
				if (!myset)
				{
					errormsg = "CharSet ";
					errormsg += token.GetToken();
					errormsg += " not found";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
				myset->model_name = model_name;
				myset->char_partition = partition_name;
				token.GetNextToken();
				if (!token.Equals(",") && !token.Equals(";"))
				{
					errormsg = "Expecting ',' or ';', but found ";
					errormsg += token.GetToken();
					errormsg += " instead";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
				if (token.Equals(";"))
					break;
				else
					token.GetNextToken();
			} while (!token.AtEOF() && !token.Equals(";"));
			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate CHARPARTITION command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}
		}
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
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
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
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}
		}	// token not END, ENDBLOCK, COST
	}	// GetNextToken loop

}

CharSet *MSetsBlock::findCharSet(string name) {
	for (vector<CharSet*>::iterator it = charsets.begin(); it != charsets.end(); it++)
		if ((*it)->name == name) return (*it);
	return NULL;
}

int MSetsBlock::findArea(string &name) {
	for (int i = 0; i < sets.size(); i++)
		if (sets[i]->name == name) return i;
	return -1;
}
