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

namespace PML {

const char* const ERR_NOT_A_LIST = "list '[...]' expected";
const char* const ERR_NOT_A_MAP = "'key: value' pairs expected";
const char* const ERR_UNDEFINED_STATE = "undefined state";
const char* const ERR_STRING_LIST = "string or list [...] expected";
const char* const ERR_TRANSLATE_LENGTH = "translate length different from #states";

const char* const KEY_DATATYPE = "datatype";
const char* const KEY_STATE = "state";
const char* const KEY_MISSING = "missing";
const char* const KEY_GAP = "gap";
const char* const KEY_EQUATE = "equate";
const char* const KEY_TRANSLATE = "translate";

const char* builtin_state_spaces = R"(
### DNA data definition ###
- datatype: DNA
  state: &Nucleotide [ A, C, G, T ]  # anchor to Nucleotide
  missing: &NTmissing [ N, "?" ]
  gap: &NTgap "-"
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

### RY data definition ###
- datatype: RY
  state: [ R, Y ] # R=AG, Y=CT
  missing: [ N, "?", W, S, M, K, B, H, D, V ]
  gap: "-"
  equate:
    A: R
    C: Y
    G: R
    T: Y
    U: Y

### Morphological data ###
- datatype: MORPH
  state: [ 0..9, A..Z ]
  missing: "?"
  gap: "-"

### Codon data with standard genetic code ###
- datatype: CODON
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: [ K, N, K, N, T, T, T, T, R, S, R, S, I, I, M, I,
               Q, H, Q, H, P, P, P, P, R, R, R, R, L, L, L, L,
               E, D, E, D, A, A, A, A, G, G, G, G, V, V, V, V,
               X, Y, X, Y, S, S, S, S, X, C, W, C, L, F, L, F ]

### Codon data with Vertebrate Mitochondrial code ###
- datatype: CODON2
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Vertebrate Mitochondrial code ###
- datatype: CODON2
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Yeast Mitochondrial code ###
- datatype: CODON3
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Mold, Protozoan code ###
- datatype: CODON4
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Invertebrate Mitochondrial code ###
- datatype: CODON5
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Ciliate, Dasycladacean and Hexamita Nuclear code ###
- datatype: CODON6
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF

### Codon data with Echinoderm and Flatworm Mitochondrial code ###
- datatype: CODON9
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Euplotid Nuclear code ###
- datatype: CODON10
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF

### Codon data with Bacterial, Archaeal and Plant Plastid code ###
- datatype: CODON11
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF

### Codon data with Alternative Yeast Nuclear code ###
- datatype: CODON12
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF

### Codon data with Ascidian Mitochondrial code ###
- datatype: CODON13
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Alternative Flatworm Mitochondrial code ###
- datatype: CODON14
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF

### Codon data with Blepharisma Nuclear code ###
- datatype: CODON15
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF

### Codon data with Chlorophycean Mitochondrial code ###
- datatype: CODON16
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF

### Codon data with Trematode Mitochondrial code ###
- datatype: CODON21
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Scenedesmus obliquus mitochondrial code ###
- datatype: CODON22
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF

### Codon data with Thraustochytrium Mitochondrial code ###
- datatype: CODON23
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF

### Codon data with Pterobranchia mitochondrial code ###
- datatype: CODON24
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF

### Codon data with Candidate Division SR1 and Gracilibacteria code ###
- datatype: CODON25
  state: [ *Nucleotide, *Nucleotide, *Nucleotide ]  # reference to Nucleotide
  missing: [ *NTmissing, *NTmissing, *NTmissing ]
  gap: [ *NTgap, *NTgap, *NTgap ]
  translate: KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF


)";

StateSpace::StateSpace() {
    num_states = 0;
    num_all_states = 0;
}

StateSpace::~StateSpace() {

}

bool StateSpace::isUnknown(StateType state) {
    return (state == num_states);
}

StateType StateSpace::toState(string str) {
    StringStateMap::iterator it;
    it = states.find(str);
    if (it == states.end())
        throw str + " is not a valid state symbol";
    return it->second;
}
    
void StateSpace::toState(string &str, StateVector &str_states) {
    size_t pos;
    for (pos = 0; pos < str.length();) {
        bool found = false;
        for (int len = min_state_len; len <= max_state_len; len++) {
            auto it = states.find(str.substr(pos, len));
            if (it == states.end())
                continue;
            found = true;
            str_states.push_back(it->second);
            pos += len;
            break;
        }
        if (!found)
            throw str.substr(pos, max_state_len) + " is not a valid state symbol";
    }
}

string StateSpace::toString(StateType state) {
    auto it = raw_states.find(state);
    ASSERT(it != raw_states.end());
    return it->second;
}

/**
 parse a string with range (e.g. 1..5) to a vector of string
 */
void parseRange(string str, StrVector &list) {
    size_t pos;
    if ((pos = str.find("..")) == string::npos) {
        list.push_back(str);
        return;
    }
    string first = str.substr(0, pos);
    string last = str.substr(pos+2);
    trimString(first);
    trimString(last);
    if (first.length() == 1 && last.length() == 1 && first[0] < last[0]) {
        for (char ch = first[0]; ch <= last[0]; ch++)
            list.push_back(string(1,ch));
    } else {
        list.push_back(str);
    }
}

/**
 parse a list into a vector of string
 */
void parseList(YAML::const_iterator first, YAML::const_iterator last, StrVector &list) {
    StrVector this_list;
    if (first->IsScalar())
        parseRange(first->Scalar(), this_list);
    else if (first->IsSequence()) {
        for (auto it = first->begin(); it != first->end(); it++) {
            parseRange(it->Scalar(), this_list);
        }
    } else {
        throw YAML::Exception(first->Mark(), ERR_STRING_LIST);
    }
    StrVector last_list;
    first++;
    if (first != last)
        parseList(first, last, last_list);
    else
        last_list = { "" };
    for (auto sit = this_list.begin(); sit != this_list.end(); sit++)
        for (auto sit2 = last_list.begin(); sit2 != last_list.end(); sit2++ )
            list.push_back(*sit + *sit2);
}

/**
 parse a YAML::Node into a list of strings
 @param extend_length TRUE to make vector of characters if list has length 1
 */
void parseStringList(YAML::Node node, StrVector &list, bool extend_length = false) {
    if (node.IsScalar()) {
        // scalar assumed to be string
        parseRange(node.Scalar(), list);
    } else if (node.IsSequence()) {
        YAML::const_iterator it;
        // check if a sequence of scalars
        bool all_scalars = true;
        for (it = node.begin(); it != node.end(); it++)
            if (!it->IsScalar()) {
                all_scalars = false;
                break;
            }
        
        if (all_scalars) {
            for (it = node.begin(); it != node.end(); it++)
                parseRange(it->Scalar(), list);
        } else {
            // now it can be a sequence of sequences, merge them together
            parseList(node.begin(), node.end(), list);
        }
    } else {
        throw YAML::Exception(node.Mark(), ERR_STRING_LIST);
    }
    
    if (list.size() == 1 && extend_length) {
        // single list, convert to vector of characters
        for (auto i = list[0].begin()+1; i != list[0].end(); i++)
            list.push_back(string(1,*i));
        list[0] = list[0].substr(0,1);
    }
}

void StateSpace::resetStateSpace() {
    space_name = "";
    num_states = 0;
    num_all_states = 0;
    states.clear();
    raw_states.clear();
    equate.clear();
    translate.clear();
    min_state_len = max_state_len = 0;
}

void StateSpace::parseStateSpace(YAML::Node datatype) {
    if (!datatype.IsMap())
        throw YAML::Exception(datatype.Mark(), ERR_NOT_A_MAP);
    if (!datatype[KEY_DATATYPE])
        throw YAML::Exception(datatype.Mark(), "'datatype: XXX' declaration not found");
    resetStateSpace();
    space_name = datatype[KEY_DATATYPE].Scalar();
    // definition found
    // parse state: symbols
    if (!datatype[KEY_STATE])
        throw YAML::Exception(datatype.Mark(), "datatype does not have 'state: [...]'");
    StrVector allstates;
    parseStringList(datatype[KEY_STATE], allstates);
    if (allstates.size() < 2)
        throw YAML::Exception(datatype[KEY_STATE].Mark(), "state space must have at least 2 states");
    StateType stateID = 0;
    for (auto sit = allstates.begin(); sit != allstates.end(); sit++, stateID++) {
        states[*sit] = stateID;
        raw_states[stateID] = *sit;
    }
    num_states = stateID;
    
    if (verbose_mode >= VB_MED)
        cout << states.size() << " " << KEY_STATE << endl;

    // parse missing: symbols
    if (datatype[KEY_MISSING]) {
        StrVector list;
        parseStringList(datatype[KEY_MISSING], list);
        for (auto i = list.begin(); i != list.end(); i++) {
            states[*i] = stateID;
            raw_states[stateID] = *i;
        }
    }
    
    // parse gap: symbols
    if (datatype[KEY_GAP]) {
        StrVector list;
        parseStringList(datatype[KEY_GAP], list);
        for (auto i = list.begin(); i != list.end(); i++) {
            states[*i] = stateID;
            raw_states[stateID] = *i;
        }
    }
    
    stateID++;
    
    // parse equate: symbols
    YAML::Node node_equate;
    if ((node_equate = datatype[KEY_EQUATE])) {
        if (!node_equate.IsMap())
            throw YAML::Exception(node_equate.Mark(), ERR_NOT_A_MAP);
        for (auto nit = node_equate.begin(); nit != node_equate.end(); nit++) {
            string key = nit->first.Scalar();
            states[key] = stateID;
            auto value = nit->second;
            StrVector values;
            parseStringList(value, values);
            for (auto i = values.begin(); i != values.end(); i++) {
                if (states.find(*i) == states.end())
                    throw YAML::Exception(value.Mark(), ERR_UNDEFINED_STATE);
                if (equate.find(stateID) == equate.end())
                    equate[stateID] = { states[*i] };
                else
                    equate[stateID].push_back(states[*i]);
            }
            if (equate[stateID].size() == 1) {
                // map to just one state, so it's not an ambiguous state
                states[key] = equate[stateID][0];
                equate.erase(stateID);
            } else {
                // increase number of states
                raw_states[stateID] = key;
                stateID++;
            }
        } // for Node
        if (verbose_mode >= VB_MED)
            cout << equate.size() << " ambiguous states" << endl;
    } // equate

    // parse translate
    if (datatype[KEY_TRANSLATE]) {
        parseStringList(datatype[KEY_TRANSLATE], translate, true);
        if (translate.size() != num_states)
            throw YAML::Exception(datatype[KEY_TRANSLATE].Mark(), ERR_TRANSLATE_LENGTH);
    }
    
    num_all_states = stateID;
    min_state_len = max_state_len = states.begin()->first.length();
    for (auto i = states.begin(); i != states.end(); i++) {
        if (min_state_len > i->first.length())
            min_state_len = i->first.length();
        if (max_state_len < i->first.length())
            max_state_len = i->first.length();
    }
}

void StateSpace::initStateSpace(SeqType seqtype) {
    
    string name;
    switch (seqtype) {
    case SEQ_DNA: name = "DNA"; break;
    case SEQ_CODON: name = "CODON"; break;
    case SEQ_MORPH: name = "MORPH"; break;
    case SEQ_BINARY: name = "BIN"; break;
    case SEQ_PROTEIN: name = "AA"; break;
    case SEQ_MULTISTATE: name = "MULTI"; break;
    case SEQ_POMO: outError("Unhandled POMO state space"); break;
    case SEQ_UNKNOWN: ASSERT(0);
    }
    
    try {
        YAML::Node spaceDef = YAML::Load(builtin_state_spaces);
        if (!spaceDef.IsSequence())
            throw YAML::Exception(spaceDef.Mark(), ERR_NOT_A_LIST);
        for (auto it = spaceDef.begin(); it != spaceDef.end(); it++)
        {
            auto datatype = *it;
            if (!(datatype[KEY_DATATYPE]))
                continue;
            if (datatype[KEY_DATATYPE].Scalar() == name) {
                parseStateSpace(datatype);
                break;
            }
        }
    } catch (YAML::Exception &e) {
        outError(e.what());
    }
}

} // namespace PML
