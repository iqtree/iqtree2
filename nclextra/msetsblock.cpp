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
#include "utils/stringfunctions.h" //for convert_double
#include "utils/tools.h"

CharSet::CharSet(): tree_len(0.0), is_ASC(false) {
}

CharSet::CharSet(const std::string& name) : CharSet() {
	model_name = name;
	trimString(model_name);
	//            std::transform(info.model_name.begin(), info.model_name.end(), 
	//                           info.model_name.begin(), ::toupper);
	
	is_ASC = model_name.substr(0,4) == "ASC_";
	if (is_ASC) {
		model_name.erase(0, 4);
	}
	StateFreqType freq = StateFreqType::FREQ_UNKNOWN;
	if (model_name.find_first_of("*+{") == string::npos ) {
		if (*model_name.rbegin() == 'F' && model_name != "DAYHOFF") {
			freq = StateFreqType::FREQ_EMPIRICAL;
			model_name.erase(model_name.length()-1);
		} else if (*model_name.rbegin() == 'X' && model_name != "LG4X") {
			freq = StateFreqType::FREQ_ESTIMATE;
			model_name.erase(model_name.length()-1);
		}
	}
}

void CharSet::setSequenceTypeAndModelNameFromString(std::string input) {
	if (input.empty()) {
		outError("Please give model names in partition file!");
	}
	if (input == "BIN") {
		sequence_type = "BIN";
		model_name    = "GTR2";
	} else if (input == "DNA") {
		sequence_type = "DNA";
		model_name    = "GTR";
	} else if (input == "MULTI") {
		sequence_type = "MORPH";
		model_name    = "MK";
	} else if (startsWith(input,"CODON")) {
		sequence_type = input;
		model_name    = "GY";
	} else {
		sequence_type = "AA";
		if (*input.begin() == '[') {
			if (*input.rbegin() != ']') {
				outError("User-defined protein model should be [myProtenSubstitutionModelFileName]");
			}
			model_name = input.substr(1, input.length()-2);
		}
	}
}

void CharSet::adjustModelName(const std::string& rate_type) {
	if (freq == StateFreqType::FREQ_EMPIRICAL)
		model_name += "+F";
	else if (freq == StateFreqType::FREQ_ESTIMATE)
		model_name += "+FO";
	if (is_ASC) {
		model_name += "+ASC";
	}
	model_name += rate_type;
}

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

    for (vector<CharSet* >::reverse_iterator it2 = charsets.rbegin(); it2 != charsets.rend(); it2++)
        delete *it2;
        
    charsets.clear();
}


void MSetsBlock::Report(ostream &out)
{
    int nsets = static_cast<int>(getNSets());
    out << "Number of sets: " << nsets << endl;
    for (TaxaSetNameVector::iterator i = sets.begin(); i != sets.end(); i++) {
        out << "Set " << (*i)->name << " contains: ";
        for (auto it = (*i)->taxlist.begin(); it != (*i)->taxlist.end(); it++) {
            out << (*it) << "  ";
        }
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
                string taxname = token.GetToken();
                renameString(taxname);
				myset->taxlist.push_back(taxname);
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
            token.SetLabileFlagBit(NxsToken::preserveUnderscores + NxsToken::hyphenNotPunctuation);
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
            trimString(myset->position_spec);
			pos = myset->position_spec.find(',');
            if (pos != string::npos && isalpha(myset->position_spec[0])) {
                myset->sequence_type = myset->position_spec.substr(0, pos);
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
            trimString(myset->aln_file);
            trimString(myset->char_partition);
            trimString(myset->model_name);
            trimString(myset->position_spec);
            trimString(myset->sequence_type);
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
                if (token.Equals("{")) {
                    token.GetNextToken();
                    myset->tree_len = convert_double(token.GetToken().c_str());
                    token.GetNextToken();
                    if (!token.Equals("}")) {
                        errormsg = "Expecting '}', but found ";
                        errormsg += token.GetToken();
                        errormsg += " instead";
                        throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
                    }
                    token.GetNextToken();
                }
                
                if (!token.Equals(",") && !token.Equals(";"))
				{
					errormsg = "Expecting ',' or ';', but found ";
					errormsg += token.GetToken();
					errormsg += " instead";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
                token.SetLabileFlagBit(NxsToken::preserveUnderscores);
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
            errormsg = "Unknown command ";
            errormsg += token.GetToken();
            throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
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

void MSetsBlock::SkippingCommand(NxsString commandName) {
    cout << "WARNING: Skipping unknown command " << commandName << endl;
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
