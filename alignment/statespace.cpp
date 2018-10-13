//
// C++ Implementation: StateSpace
//
// Description:
//
//
// Author: BUI Quang Minh(C) 2018
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "statespace.h"
#include "yaml-cpp/yaml.h"

const char* builtin_state_spaces = R"(
### DNA data definition ###
- datatype: DNA
state: &Nucleotide [ A, C, G, T ]  # anchor to Nucleotide
  missing: [ N, "?" ]
  gap: "-"
  equate:
    U: T      # T and U are the same
    R: [A, G] # R is interpreted as A or G
    Y: [C, T]
    W: [A, T]
    S: [G, C]
    M: [A, C]
    K: [G, T]
    B: [C, G, T]
    H: [A, C, T]
    D: [A, G, T]
    V: [A, G, C]

### Amino-acid data definition ###
- datatype: AA
  state: [ A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V ]
  missing: [ X, "?", "*" ]
  gap: "-"
  equate:
    B: [ N, D ]
    Z: [ Q, E ]
    J: [ I, L ]

### Binary (0/1) data ###
- datatype: BIN
  state: [ 0, 1 ]
  missing: "?"
  gap: "-"

### Morphological data ###
- datatype: MORPH
  state: [ 0..9, A..Z ]
  missing: "?"
  gap: "-"

### Codon data with standard genetic code ###

- datatype: CODON
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  drop: [ TAA, TAG, TGA ]  # remove stop codons
  missing: [ "NNN", "???" ]
  gap: "---"
  translate: [ ]
  translateTo: AminoAcid

)";

/**
 parse a string with range (e.g. 1..5) to a vector of string
 */
void parseStates(string state, StrVector &states) {
    states.clear();
    size_t pos;
    if ((pos = state.find("..")) == string::npos) {
        states.push_back(state);
        return;
    }
    string first = state.substr(0, pos);
    string last = state.substr(pos+2);
    trimString(first);
    trimString(last);
    if (first.length() == 1 && last.length() == 1 && first[0] < last[0]) {
        for (char ch = first[0]; ch <= last[0]; ch++)
            states.push_back(string(1,ch));
    } else {
        states.push_back(state);
    }
}

void StateSpace::initStateSpace(string name) {
    YAML::Node spaceDef = YAML::Load(builtin_state_spaces);
    ASSERT(spaceDef.IsSequence());
    for (auto it = spaceDef.begin(); it != spaceDef.end(); it++) {
        auto datatype = (*it);
        YAML::Node node;
        if (!(node = datatype["datatype"]))
            continue;
        if (node.as<string>() != name)
            continue;
        cout << datatype << endl;
        // definition found
        // parse state: symbols
        if (!(node = datatype["state"]))
            outError("No 'state: ' for datatype " + name);
        if (!node.IsSequence()) {
            cerr << node;
            outError(" is not a list");
        }
        StateType stateID = 0;
        YAML::const_iterator state;
        for (state = node.begin(); state != node.end(); state++) {
            StrVector str;
            parseStates(state->as<string>(), str);
            for (auto sit = str.begin(); sit != str.end(); sit++, stateID++) {
                states[*sit] = stateID;
            }
        }
        cout << states.size() << " states" << endl;

        // parse equate: symbols
        if ((node = datatype["equate"])) {
            if (!node.IsMap()) {
                cerr << node;
                outError(" is not a map");
            }
            for (auto nit = node.begin(); nit != node.end(); nit++) {
                string key = nit->first.as<string>();
                auto value = nit->second;
                if (value.IsScalar()) {
                    if (states.find(value.as<string>()) == states.end()) {
                        cerr << *nit;
                        outError(" undefined");
                    }
                    equate[key] = { states[value.as<string>()] };
                    continue;
                }
                if (!value.IsSequence()) {
                    cerr << value;
                    outError("is not a list");
                }
                for (state = value.begin(); state != value.end(); state++) {
                    if (states.find(state->as<string>()) == states.end()) {
                        cerr << *state;
                        outError(" undefined");
                    }
                    if (equate.find(key) == equate.end())
                        equate[key] = { states[state->as<string>()] };
                    else
                        equate[key].push_back(states[state->as<string>()]);
                }
            } // for Node
            cout << equate.size() << " ambiguous states" << endl;
        } // equate

        // parse missing: symbols
        if ((node = datatype["missing"])) {
            if (node.IsScalar())
                equate[node.as<string>()] = { (StateType)states.size() };
            else if (node.IsSequence()) {
                for (state = node.begin(); state != node.end(); state++)
                    equate[state->as<string>()] = { (StateType)states.size() };
            } else {
                cerr << node;
                outError(" invalid");
            }
        }

        // parse gap: symbols
        if ((node = datatype["gap"])) {
            if (node.IsScalar())
                equate[node.as<string>()] = { (StateType)states.size() };
            else if (node.IsSequence()) {
                for (state = node.begin(); state != node.end(); state++)
                    equate[state->as<string>()] = { (StateType)states.size() };
            } else {
                cerr << node;
                outError(" invalid");
            }
        }
        
        cout << equate.size() << " ambiguous states" << endl;

    } // for spaceDef
}

void StateSpace::initStateSpace(SeqType seqtype) {
    try {
        switch (seqtype) {
            case SEQ_DNA: initStateSpace("CODON");
                break;
            default:
                break;
        }
    } catch (YAML::Exception &e) {
        outError(e.what());
    }
}
