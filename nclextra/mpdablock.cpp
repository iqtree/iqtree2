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
#include "mpdablock.h"
#include "pda/split.h"
#include "pda/splitgraph.h"

MPdaBlock::MPdaBlock(SplitGraph *asgraph)
 : NxsBlock()
{
	budget = -1;
	min_budget = -1;
	sub_size = 0;
	cost_constrained = false;
	id = "PDA";
	sgraph = asgraph;
}


void MPdaBlock::Report(ostream &out)
{
	out << "Budget = " << budget << endl;
	out << "Taxa Costs = ";
	for (DoubleVector::iterator it = costs.begin(); it != costs.end(); it++)
		out << *it << " ";
	out << endl;
}

void MPdaBlock::Reset()
{
	errormsg.clear();
	isEmpty			= true;
	isEnabled		= true;
	isUserSupplied	= false;

}

void MPdaBlock::Read(NxsToken &token)
{

	int ntax = sgraph->getNTaxa();
	if (ntax <= 0) {
		errormsg = "PDA Block should be preceded by Splits Block";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		
	}

	//int nominal_ntax	= 0;
	//int nominal_nsplits = 0;

	// This should be the semicolon after the block name
	//
	token.GetNextToken();

	if (!token.Equals(";"))
	{
		errormsg = "Expecting ';' after PDA block name, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	for (;;)
	{
		token.GetNextToken();

		if (token.Equals("PARAMETERS"))
		{
			// This should be the NTAX keyword
			//
			token.GetNextToken();

			do {

				if (token.Equals("BUDGET")) {
					// This should be the equals sign
					//
					token.GetNextToken();
		
					if (!token.Equals("="))
					{
						errormsg = "Expecting '=', but found ";
						errormsg += token.GetToken();
						errormsg += " instead";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
		
					// This should be the integer budget
					//
					token.GetNextToken();
		
					budget = convert_double(token.GetToken().c_str());
					if (budget <= 0)
					{
						errormsg = "BUDGET should be greater than zero (";
						errormsg += token.GetToken();
						errormsg += " was specified)";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
				} else if (token.Equals("MIN BUDGET")) {
					// This should be the equals sign
					//
					token.GetNextToken();
		
					if (!token.Equals("="))
					{
						errormsg = "Expecting '=', but found ";
						errormsg += token.GetToken();
						errormsg += " instead";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
		
					// This should be the integer budget
					//
					token.GetNextToken();
		
					min_budget = convert_double(token.GetToken().c_str());
					if (budget < 0)
					{
						errormsg = "MIN_BUDGET should be greater than or equal to zero (";
						errormsg += token.GetToken();
						errormsg += " was specified)";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
		
				} else if (token.Equals("BUDGET CONSTRAINED")) {
					cost_constrained = true;
				} else if (token.Equals("K")) {
					// This should be the equals sign
					//
					token.GetNextToken();
		
					if (!token.Equals("="))
					{
						errormsg = "Expecting '=', but found ";
						errormsg += token.GetToken();
						errormsg += " instead";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
		
					// This should be the integer budget
					//
					token.GetNextToken();
		
					sub_size = atoi(token.GetToken().c_str());
					if (sub_size <= 1)
					{
						errormsg = "K should be greater than 1 (";
						errormsg += token.GetToken();
						errormsg += " was specified)";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
				} else
				{
					errormsg = "Invalid PARAMETERS command: ";
					errormsg += token.GetToken();
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

				token.GetNextToken();

			} while (!token.AtEOF() && !token.Equals(";"));

			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate PARAMETERS command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

		} // if (token.Equals("PARAMETERS"))

		else if (token.Equals("TAXCOSTS")) {

			costs.resize(ntax, -1);
			// This should be taxon name
			//
			token.GetNextToken();

			do {				
				int tax_id = -1;

				try {
					tax_id = sgraph->getTaxa()->FindTaxon(token.GetToken());
				} catch (NxsTaxaBlock::NxsX_NoSuchTaxon) {
					tax_id = -1;
				}

				if (tax_id < 0)
				{
					errormsg = "Taxon is not found (";
					errormsg += token.GetToken();
					errormsg += " was specified)";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
		
				// This should be the cost of taxon
				//
				token.GetNextToken();

				int taxcost = convert_double(token.GetToken().c_str());
				if (taxcost < 0)
				{
					errormsg = "Taxon cost should be greater than or equal to zero (";
					errormsg += token.GetToken();
					errormsg += " was specified)";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
				costs[tax_id] = taxcost;

				token.GetNextToken();
			} while (!token.AtEOF() && !token.Equals(";"));

			// This should be the terminating semicolon
			//
			//token.GetNextToken();

			if (!token.Equals(";"))
			{
				errormsg = "Expecting ';' to terminate TAXCOSTS command, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

			for (int i = 0; i < ntax; i++)
				if (costs[i] < 0) {
					costs[i] = 0;
					cout << "WARNING: taxon " << sgraph->getTaxa()->GetTaxonLabel(i)
						<< "has no cost! set to 0." << endl;
				}
		}	// if (token.Equals("TAXCOSTS"))


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

void MPdaBlock::readBudgetFile(Params &params) {
	ifstream in;
	in.exceptions(ios::failbit | ios::badbit);
	cout << "Reading budget information file " << params.budget_file << "..." << endl;
	NxsString taxname;
	int ntaxa = sgraph->getNTaxa() - params.is_rooted;
	int i;

	try {
		costs.resize(ntaxa, -1);
		in.open(params.budget_file);
		in >> budget;
		if (budget < 0) 
			throw "Negative total budget.";
		for (i = 0; i < ntaxa && !in.eof(); i++) {
			double taxcost;
			int tax_id = -1;
			taxname = "";
			in.exceptions(ios::badbit);
			in >> taxname;
			in.exceptions(ios::failbit | ios::badbit);
			if (taxname == "") break;
			in >> taxcost;
			if (taxcost < 0) 
				throw "Negative taxa preservation cost.";
			tax_id = sgraph->getTaxa()->FindTaxon(taxname);
			costs[tax_id] = taxcost;
		}
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (NxsTaxaBlock::NxsX_NoSuchTaxon) {
		outError(ERR_NO_TAXON, taxname);
	} catch (const char *str) {
		outError(str);
	} catch (...) {
		// anything else
		outError(ERR_READ_ANY);
	}

	for (i = 0; i < ntaxa; i++)
		if (costs[i] < 0) {
			costs[i] = 0;
			cout << "WARNING: taxon " << sgraph->getTaxa()->GetTaxonLabel(i)
				<< "has no cost! set to 0." << endl;
		}
	cost_constrained = true;
}

void MPdaBlock::readBudgetAreaFile(Params &params) {
	ifstream in;
	in.exceptions(ios::failbit | ios::badbit);
	cout << "Reading budget for areas information file " << params.budget_file << "..." << endl;
	string areaname;
	int nareas = sgraph->getNAreas();
	int i;

	try {
		costs.resize(nareas, -1);
		in.open(params.budget_file);
		in >> budget;
		if (budget < 0) 
			throw "Negative total budget.";
		for (i = 0; i < nareas && !in.eof(); i++) {
			double areacost;
			int area_id = -1;
			areaname = "";
			in.exceptions(ios::badbit);
			in >> areaname;
			in.exceptions(ios::failbit | ios::badbit);
			if (areaname == "") break;
			in >> areacost;
			if (areacost < 0) 
				throw "Negative taxa preservation cost.";
			area_id = sgraph->getSetsBlock()->findArea(areaname);
			if (area_id < 0)
				outError(ERR_NO_AREA, areaname);
			costs[area_id] = areacost;
		}
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (const char *str) {
		outError(str);
	} catch (...) {
		// anything else
		outError(ERR_READ_ANY);
	}

	for (i = 0; i < nareas; i++)
		if (costs[i] < 0) {
			costs[i] = 0;
			cout << "WARNING: area " << sgraph->getSetsBlock()->getSet(i)->name
				<< "has no cost! set to 0." << endl;
		}
	cost_constrained = true;
}


void MPdaBlock::SkippingCommand(NxsString commandName) {
	cout << "   Skipping unknown command (" << commandName << ")..." << endl;
}


MPdaBlock::~MPdaBlock()
{
}
