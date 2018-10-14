
//
// C++ Interface: StateSpace
//
// Description:
//
//
// Author: BUI Quang Minh (c) 2018
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef STATESPACE_H
#define STATESPACE_H

#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>
#include "utils/tools.h"
#include "yaml-cpp/yaml.h"

/**
 StateType as 32-bit unsigned int
 */
typedef uint32_t StateType;

enum SeqType {
    SEQ_DNA, SEQ_PROTEIN, SEQ_BINARY, SEQ_MORPH, SEQ_MULTISTATE, SEQ_CODON, SEQ_POMO, SEQ_UNKNOWN
};

// IMPORTANT: refactor STATE_UNKNOWN
//const char STATE_UNKNOWN = 126;

// TODO DS: This seems like a significant restriction.
/* PoMo: STATE_INVALID is not handled in PoMo.  Set STATE_INVALID to
 127 to remove warning about comparison to char in alignment.cpp.
 This is important if the maximum N will be increased above 21
 because then the state space is larger than 127 and we have to
 think about something else. */
/* const unsigned char STATE_INVALID = 255; */
const unsigned char STATE_INVALID = 127;

#ifdef USE_HASH_MAP
typedef unordered_map<string, int> StringIntMap;
typedef unordered_map<string, StateType> StringStateMap;
typedef unordered_map<string, double> StringDoubleHashMap;
typedef unordered_map<uint32_t, uint32_t> IntIntMap;
#else
typedef map<string, int> StringIntMap;
typedef map<string, StateType> StringStateMap;
typedef map<string, double> StringDoubleHashMap;
typedef map<uint32_t, uint32_t> IntIntMap;
#endif


/**
 general class defining state space
 */
class StateSpace {
public:
    /** constructor */
    StateSpace() {}

    /** destructor */
    ~StateSpace() {}

    /**
     initialise from a state definition string
     */
    void parseStateSpace(YAML::Node datatype);

    void initStateSpace(SeqType seqtype);

protected:

    /** state space name */
    string space_name;

    /** map from raw state string to state ID */
    StringStateMap states;
    
    /** map from ambiguous states to vector of state ID */
    unordered_map<string, vector<StateType> >equate;
    
    /** vector of the same size as states to translate to another state space */
    StrVector translate;
};

#endif
