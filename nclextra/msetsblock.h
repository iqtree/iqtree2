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
#ifndef MSETSBLOCK_H
#define MSETSBLOCK_H

#include "ncl/ncl.h"

/**
	a taxa set with name
*/
class TaxaSetName {
public:
	/**
		set name
	*/
	string name;
	
	/**
		string vector of taxa names
	*/
	vector<string> taxlist;
};

typedef vector<TaxaSetName*> TaxaSetNameVector;

/**
 * a charset
 */
class CharSet {
public:
    
    /** constructor */
    CharSet() {
        tree_len = 0.0;
    }
    
	/** charset name */
	string name;

	/** positions specification, e.g., 1-500\3 501-502 */
	string position_spec;

	/** name of model associated with charset, e.g., GTR+G */
	string model_name;

	/** alignment name */
	string aln_file;

	/** sequence type */
	string sequence_type;

	/** name of CharPartition where this charset is included*/
	string char_partition;
    
    /** total tree length for this charset */
    double tree_len;
};


/**
Sets Block of Nexus file parser

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class MSetsBlock : public NxsBlock
{
public:

	/**
		constructor, assigning an associated splits graph
	*/
    MSetsBlock();

	/**
		destructor
	*/
    virtual ~MSetsBlock();

	/**
		print info to an output stream
		@param out output stream, cout for output to screen
	*/
	virtual void Report(ostream &out);

	/**
		reset the block
	*/
	virtual void Reset();

    /**
     called when some commands are skipped
     @param commandName command name
     */
    virtual void SkippingCommand(NxsString commandName);

	/**
		@return the number of sets
	*/
	size_t getNSets() const { return sets.size(); }

	/**
		@param id set id
		@return reference to the corresponding set
	*/
	inline TaxaSetName *getSet(int id) { return sets[id]; }

	/**
		@return vector of all taxa set
	*/
	inline TaxaSetNameVector *getSets() { return &sets; }

	/**
		@param name an area name
		@return ID of the area with that name, -1 if not found
	*/
	int findArea(string &name);


	/** list of charsets with (possible) models */
	vector<CharSet* > charsets;

	/**
	 * return CharSet with a name from charsets, NULL if not found
	 */
	CharSet *findCharSet(string name);

protected:

	/**
		main method to read block from file
		@param token a token reader
	*/
	virtual void Read(NxsToken &token);


	/**
		list of taxa set names
	*/
	TaxaSetNameVector sets;

};

#endif
