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

#include <iostream>
#include <string>

using namespace std;

/**
	Site-patterns in a multiple sequence alignment
	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class Pattern : public string
{
public:
	/** 
		constructor
	*/
    Pattern();

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
	int computeGapChar(int num_states, int STATE_UNKNOWN);

//    Pattern &operator= (Pattern pat);

	/** 
		destructor
	*/
    virtual ~Pattern();

	/**
		frequency appearance of the pattern
	*/
	int frequency;

	/**
		true if this is a constant pattern
		2015-03-04: is_const will also be true for pattern like "AA-A--AAA"
	*/
	bool is_const;
    
    /** true if pattern is informative, false otherwise */
    bool is_informative;

	/** 2015-03-04: if is_const is true, this will store the const character for the pattern */
	char const_char;

    /** number of different character states */
    int num_chars;
};

#endif
