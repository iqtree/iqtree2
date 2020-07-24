
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

namespace PML {

/**
 StateType as 32-bit unsigned int
 */
typedef uint32_t StateType;

typedef vector<StateType> StateVector;

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
typedef unordered_map<StateType, string> StateStringMap;
typedef unordered_map<string, double> StringDoubleHashMap;
typedef unordered_map<uint32_t, uint32_t> IntIntMap;
#else
typedef map<string, int> StringIntMap;
typedef map<string, StateType> StringStateMap;
typedef map<StateType, string> StateStringMap;
typedef map<string, double> StringDoubleHashMap;
typedef map<uint32_t, uint32_t> IntIntMap;
#endif


/**
 general class defining state space
 */
class StateSpace {
public:
    /** constructor */
    StateSpace();

    /** destructor */
    ~StateSpace();

    /** convert a raw string to single state ID */
    StateType toState(string str);
    
    /**
    convert the entire string into vector of states
    @param[in] str input string
    @param[out] str_states output vector of StateType
    */
    void toState(string &str, StateVector &str_states);
    
    /** convert a state back to raw string */
    string toString(StateType state);

    /**
    check if a state is unknown (missing or gap)
    */
    bool isUnknown(StateType state);

    /** get number of states */
    inline int getNStates() { return num_states; }

    /** get all number of states incl. missing/gap/ambiguous states */
    inline size_t getNAllStates() { return states.size(); }

    /**
     initialise from a state definition string
     @param datatype a YAML::Node structure
     */
    void parseStateSpace(YAML::Node datatype);

    /**
     initialise state space from a SeqType
     @param seqtype sequence type
    */
    void initStateSpace(SeqType seqtype);

    /**
    reset state space
    */
    void resetStateSpace();

    /** number of state */
    int num_states;

protected:

    /** state space name */
    string space_name;

    /** number of state */
    int num_all_states;

    /** map from raw state string to state ID */
    StringStateMap states;

    /** map from state ID to raw state string */
    StateStringMap raw_states;

    /** map from ambiguous states to vector of state ID */
    unordered_map<StateType, StateVector>equate;
    
    /** vector of the same size as states to translate to another state space */
    StrVector translate;

private:

    /** minimum length of state string */
    int min_state_len;

    /** maximum length of state string */
    int max_state_len;

};

} // namespace PML

#endif
