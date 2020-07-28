//
// C++ Interface: pattern
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PATTERN_H
#define PATTERN_H

#include "phylo-yaml/statespace.h"

using namespace std;
using namespace PML;

const int PAT_CONST       = 1; // const site pattern, e.g. AAAAAA, CC-C-CCCC
const int PAT_INVARIANT   = 2; // invariant site pattern, including const patterns and e.g., GS--G-GGG (S = G/C)
const int PAT_INFORMATIVE = 4; // parsimony informative sites
const int PAT_VARIANT     = 8; // variant site pattern

/**
	Site-patterns in a multiple sequence alignment
	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class Pattern : public vector<StateType>
{
public:
	/** 
		constructor
	*/
    Pattern();

    /**
     constructor
     */
    Pattern(int nseq, int freq = 1);

    Pattern(const Pattern &pat);

    /**
		@param num_states number of states of the model
		@return the number of ambiguous character incl. gaps 
	*/
	int computeAmbiguousChar(int num_states);

	/**
		@param num_states number of states of the model
		@return the number of gaps 
	*/
	int computeGapChar(int num_states, int STATE_UNKNOWN) const;

//    Pattern &operator= (Pattern pat);

	/** 
		destructor
	*/
    virtual ~Pattern();

    inline bool isConst() const {
        return (flag & PAT_CONST) != 0;
    }

    inline bool isInvariant() const {
        return (flag & PAT_INVARIANT) != 0;
    }

    inline bool isInformative() const {
        return (flag & PAT_INFORMATIVE) != 0;
    }

	/**
		frequency appearance of the pattern
	*/
	int frequency;

	/**
		true if this is a constant pattern
		2015-03-04: is_const will also be true for pattern like "AA-A--AAA"
	*/
//	bool is_const;
    
    /** true if pattern is informative, false otherwise */
//    bool is_informative;

    int flag;

	/** 2015-03-04: if is_const is true, this will store the const character for the pattern */
	char const_char;

    /** number of different character states */
    int num_chars;
};

#endif
