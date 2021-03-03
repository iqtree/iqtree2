//
//  statefrequency.h
//  This file, created by James Barbetti on 4/2/21.
//  But, most of the content was taken from utils/tools.h
//  (functions and types from other places are tagged to
//  indicate where they are from).

#ifndef statefrequency_h
#define statefrequency_h

#include <string>

//This enumerated type came from phyloYAML/statespace.h
enum SeqType {
    SEQ_DNA, SEQ_PROTEIN, SEQ_BINARY, SEQ_MORPH, SEQ_MULTISTATE, SEQ_CODON, SEQ_POMO, SEQ_UNKNOWN
};

//This function came from main/phylotesting.cpp
std::string getSeqTypeName(SeqType seq_type);

//This function came from alignment/alignment.cpp
//(where it was Alignment::getSeqType)
SeqType getSeqType(const char *sequence_type);

/**
        State frequency type
 */
enum StateFreqType {
    FREQ_UNKNOWN, FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, FREQ_ESTIMATE,
    FREQ_CODON_1x4, FREQ_CODON_3x4, FREQ_CODON_3x4C, // special frequency for codon model
    FREQ_MIXTURE, // mixture-frequency model
    // FREQ_DNA_RY has pi_A+pi_G = 0.5 = pi_C+pi_T. Similarly WS pairs (AT)(CG),
    // MK pairs (AC)(GT) in same way.
    FREQ_DNA_RY, FREQ_DNA_WS, FREQ_DNA_MK,
    // in following, digits indicate which frequencies must equal each other
    // (in ACGT order), e.g. 2131 means pi_C=pi_T (pi_A, pi_G unconstrained)
    FREQ_DNA_1112, FREQ_DNA_1121, FREQ_DNA_1211, FREQ_DNA_2111,
    FREQ_DNA_1122, FREQ_DNA_1212, FREQ_DNA_1221,
    FREQ_DNA_1123, FREQ_DNA_1213, FREQ_DNA_1231,
    FREQ_DNA_2113, FREQ_DNA_2131, FREQ_DNA_2311,
};

/*
 * Given a model name, look in it for "+F..." and
 * determine the StateFreqType. Returns FREQ_UNKNOWN if
 * unable to find a good +F... specifier
 */
StateFreqType parseStateFreqFromPlusF(std::string model_name);

bool parseStateFrequencyTypeName(const std::string& name,
                                 StateFreqType& freq);

bool parseStateFrequencyTypeName(const char* name,
                                 StateFreqType& freq);

/*
 * Given a string of 4 digits, return a StateFreqType according to
 * equality constraints expressed by those digits.
 * E.g. "1233" constrains pi_G=pi_T (ACGT order, 3rd and 4th equal)
 * which results in FREQ_DNA_2311. "5288" would give the same result.
 */
StateFreqType parseStateFreqDigits(std::string digits);

/*
 * All params in range [0,1]
 * returns true if base frequencies have changed as a result of this call
 */
bool freqsFromParams(double *freq_vec, double *params, StateFreqType freq_type);

/*
 * For given freq_type, derives frequency parameters from freq_vec
 * All parameters are in range [0,1] (assuming freq_vec is valid)
 */
void paramsFromFreqs(double *params, double *freq_vec, StateFreqType freq_type);

/*
 * Given a DNA freq_type and a base frequency vector, alter the
 * base freq vector to conform with the constraints of freq_type
 */
void forceFreqsConform(double *base_freq, StateFreqType freq_type);

/*
 * For given freq_type, how many parameters are needed to
 * determine frequenc vector?
 * BQM 2017-04-28: works for DNA and other data types
 */
 int nFreqParams(StateFreqType freq_type);

/*
 * For freq_type, and given every base must have frequency >= min_freq, set upper
 * and lower bounds for parameters.
 */
 void setBoundsForFreqType(double *lower_bound,
                           double *upper_bound,
                           bool *bound_check,
                           double min_freq,
                           StateFreqType freq_type);


// This function came from model/modelmarkov.h
std::string freqTypeString(StateFreqType freq_type, SeqType seq_type, bool full_str);


#endif /* statefrequency_h */
