//
//  statefrequency.cpp
//  alignment
//
//  Created by James Barbetti on 4/2/21.
//

#include "statefrequency.h"
#include "tools.h" //for outError function, and ASSERT macro


//This function came from main/phylotesting.cpp
std::string getSeqTypeName(SeqType seq_type) {
    switch (seq_type) {
        case SEQ_BINARY:     return "binary";
        case SEQ_DNA:        return "DNA";
        case SEQ_PROTEIN:    return "protein";
        case SEQ_CODON:      return "codon";
        case SEQ_MORPH:      return "morphological";
        case SEQ_POMO:       return "PoMo";
        case SEQ_UNKNOWN:    return "unknown";
        case SEQ_MULTISTATE: return "MultiState";
        default:             return "unknown";
    }
}

//This function came from alignment/alignment.cpp
//(where it was Alignment::getSeqType)
SeqType getSeqType(const char *sequence_type) {
    SeqType user_seq_type = SEQ_UNKNOWN;
    if (strcmp(sequence_type, "BIN") == 0) {
        user_seq_type = SEQ_BINARY;
    } else if (strcmp(sequence_type, "NT") == 0 || strcmp(sequence_type, "DNA") == 0) {
        user_seq_type = SEQ_DNA;
    } else if (strcmp(sequence_type, "AA") == 0 || strcmp(sequence_type, "PROT") == 0) {
        user_seq_type = SEQ_PROTEIN;
    } else if (strncmp(sequence_type, "NT2AA", 5) == 0) {
        user_seq_type = SEQ_PROTEIN;
    } else if (strcmp(sequence_type, "NUM") == 0 || strcmp(sequence_type, "MORPH") == 0) {
        user_seq_type = SEQ_MORPH;
    } else if (strcmp(sequence_type, "TINA") == 0 || strcmp(sequence_type, "MULTI") == 0) {
        user_seq_type = SEQ_MULTISTATE;
    } else if (strncmp(sequence_type, "CODON", 5) == 0) {
        user_seq_type = SEQ_CODON;
    }
    return user_seq_type;
}



/*
 * Given a model name, look in it for "+F..." and
 * determine the StateFreqType. Returns FREQ_EMPIRICAL if
 * unable to find a good +F... specifier
 */
StateFreqType parseStateFreqFromPlusF(std::string model_name) {
//    StateFreqType freq_type = FREQ_EMPIRICAL;

    // BQM 2017-05-02: change back to FREQ_UNKNOWN to resemble old behavior
    StateFreqType freq_type = FREQ_UNKNOWN;
    size_t plusFPos;
    if (model_name.find("+F1X4") != std::string::npos)
        freq_type = FREQ_CODON_1x4;
    else if (model_name.find("+F3X4C") != std::string::npos)
        freq_type = FREQ_CODON_3x4C;
    else if (model_name.find("+F3X4") != std::string::npos)
        freq_type = FREQ_CODON_3x4;
    else if (model_name.find("+FQ") != std::string::npos)
        freq_type = FREQ_EQUAL;
    else if (model_name.find("+FO") != std::string::npos)
        freq_type = FREQ_ESTIMATE;
    else if (model_name.find("+FU") != std::string::npos)
        freq_type = FREQ_USER_DEFINED;
    else if (model_name.find("+FRY") != std::string::npos)
        freq_type = FREQ_DNA_RY;
    else if (model_name.find("+FWS") != std::string::npos)
        freq_type = FREQ_DNA_WS;
    else if (model_name.find("+FMK") != std::string::npos)
        freq_type = FREQ_DNA_MK;
    else if ((plusFPos = model_name.find("+F")) != std::string::npos) {
        freq_type = FREQ_EMPIRICAL;
        // Now look for +F#### where #s are digits
        if (model_name.length() > plusFPos+2 && isdigit(model_name[plusFPos+2]))
        try {
            // throws if string is not 4 digits
            freq_type = parseStateFreqDigits(model_name.substr(plusFPos+2,4));
        } catch (const char *str) {
            // +F exists, but can't parse it as anything else
            outError(str);
        }
    }
    return(freq_type);
}

bool parseStateFrequencyTypeName(const std::string& name,
                                 StateFreqType& freq) {
    if (name=="q" || name=="EQ") {
        freq = FREQ_EQUAL;
    }
    else if (name=="c" || name=="EM") {
        freq = FREQ_EMPIRICAL;
    }
    else if (name=="o" || name=="ES") {
        freq = FREQ_ESTIMATE;
    }
    else if (name=="u" || name=="UD") {
        freq = FREQ_USER_DEFINED;
    }
    else if (name=="ry" || name=="RY") {
        freq = FREQ_DNA_RY;
    }
    else if (name=="ws" || name=="WS") {
        freq = FREQ_DNA_WS;
    }
    else if (name=="mk" || name=="MK") {
        freq = FREQ_DNA_MK;
    }
    else {
        return false;
    }
    return true;
}

bool parseStateFrequencyTypeName(const char* name,
                                 StateFreqType& freq) {
    std::string str(name);
    return parseStateFrequencyTypeName(str, freq);
}

/*
 * Given a string of 4 digits, return a StateFreqType according to
 * equality constraints expressed by those digits.
 * E.g. "1233" constrains pi_G=pi_T (ACGT order, 3rd and 4th equal)
 * which results in FREQ_DNA_2311. "5288" would give the same result.
 */

StateFreqType parseStateFreqDigits(std::string digits) {
    bool good = true;
    if (digits.length()!=4) {
        good = false;
    } else {
        // Convert digits to canonical form, first occuring digit becomes 1 etc.
        int digit_order[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        int first_found = 0;
        for (int i=0; i<4; ++i) {
            int digit = digits[i]-'0';
            if (digit<0 || digit>9) {
    good = false; // found a non-digit
    break;
            }
            if (digit_order[digit]==-1) {
    // haven't seen this digit before
    digit_order[digit]=++first_found;
            }
            // rewrite digit in canonical form
            digits[i] = '0'+digit_order[digit];
        }
        // e.g. if digits was "5288", digit_order will end up as {-1,-1,2,-1,-1,1,-1,-1,3,-1}
    }
    if (!good) throw "Use -f <c | o | u | q | ry | ws | mk | <digit><digit><digit><digit>>";
    StateFreqType freq_type = FREQ_UNKNOWN;
    // Now just exhaustively list all canonical digit possibilities
    if (digits.compare("1111")==0) {
        freq_type = FREQ_EQUAL;
    } else if (digits.compare("1112")==0) {
        freq_type = FREQ_DNA_1112;
    } else if (digits.compare("1121")==0) {
        freq_type = FREQ_DNA_1121;
    } else if (digits.compare("1211")==0) {
        freq_type = FREQ_DNA_1211;
    } else if (digits.compare("1222")==0) {
        freq_type = FREQ_DNA_2111;
    } else if (digits.compare("1122")==0) {
        freq_type = FREQ_DNA_1122;
    } else if (digits.compare("1212")==0) {
        freq_type = FREQ_DNA_1212;
    } else if (digits.compare("1221")==0) {
        freq_type = FREQ_DNA_1221;
    } else if (digits.compare("1123")==0) {
        freq_type = FREQ_DNA_1123;
    } else if (digits.compare("1213")==0) {
        freq_type = FREQ_DNA_1213;
    } else if (digits.compare("1231")==0) {
        freq_type = FREQ_DNA_1231;
    } else if (digits.compare("1223")==0) {
        freq_type = FREQ_DNA_2113;
    } else if (digits.compare("1232")==0) {
        freq_type = FREQ_DNA_2131;
    } else if (digits.compare("1233")==0) {
        freq_type = FREQ_DNA_2311;
    } else if (digits.compare("1234")==0) {
        freq_type = FREQ_ESTIMATE;
    } else
        throw ("Unrecognized canonical digits - Can't happen"); // paranoia is good.
    return freq_type;
}


/*
 * All params in range [0,1]
 * returns true if base frequencies have changed as a result of this call
 */

bool freqsFromParams(double *freq_vec, double *params, StateFreqType freq_type) {

    // BQM 2017-05-02: Note that only freq for A, C, G are free parameters and stored
    // in params, whereas freq_T is not free and should be handled properly

    double pA, pC, pG, pT; // base freqs
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        return false;
    case FREQ_ESTIMATE:
    // Minh: in code review, please pay extra attention to ensure my treadment of FREQ_ESTIMATE is equivalent to old treatment.
    // BQM: DONE!
        pA=params[0];
        pC=params[1];
        pG=params[2];
        //pT=1-pA-pC-pG;
        pT=freq_vec[3];
        break;
    case FREQ_DNA_RY:
        pA = params[0]/2;
        pG = 0.5-pA;
        pC = params[1]/2;
        pT = 0.5-pC;
        break;
    case FREQ_DNA_WS:
        pA = params[0]/2;
        pT = 0.5-pA;
        pC = params[1]/2;
        pG = 0.5-pC;
        break;
    case FREQ_DNA_MK:
        pA = params[0]/2;
        pC = 0.5-pA;
        pG = params[1]/2;
        pT = 0.5-pG;
        break;
    case FREQ_DNA_1112:
        pA = pC = pG = params[0]/3;
        pT = 1-3*pA;
        break;
    case FREQ_DNA_1121:
        pA = pC = pT = params[0]/3;
        pG = 1-3*pA;
        break;
    case FREQ_DNA_1211:
        pA = pG = pT = params[0]/3;
        pC = 1-3*pA;
        break;
    case FREQ_DNA_2111:
        pC = pG = pT = params[0]/3;
        pA = 1-3*pC;
        break;
    case FREQ_DNA_1122:
        pA = params[0]/2;
        pC = pA;
        pG = 0.5-pA;
        pT = pG;
        break;
    case FREQ_DNA_1212:
        pA = params[0]/2;
        pG = pA;
        pC = 0.5-pA;
        pT = pC;
        break;
    case FREQ_DNA_1221:
        pA = params[0]/2;
        pT = pA;
        pC = 0.5-pA;
        pG = pC;
        break;
    case FREQ_DNA_1123:
        pA = params[0]/2;
        pC = pA;
        pG = params[1]*(1-2*pA);
        pT = 1-pG-2*pA;
        break;
    case FREQ_DNA_1213:
        pA = params[0]/2;
        pG = pA;
        pC = params[1]*(1-2*pA);
        pT = 1-pC-2*pA;
        break;
    case FREQ_DNA_1231:
        pA = params[0]/2;
        pT = pA;
        pC = params[1]*(1-2*pA);
        pG = 1-pC-2*pA;
        break;
    case FREQ_DNA_2113:
        pC = params[0]/2;
        pG = pC;
        pA = params[1]*(1-2*pC);
        pT = 1-pA-2*pC;
        break;
    case FREQ_DNA_2131:
        pC = params[0]/2;
        pT = pC;
        pA = params[1]*(1-2*pC);
        pG = 1-pA-2*pC;
        break;
    case FREQ_DNA_2311:
        pG = params[0]/2;
        pT = pG;
        pA = params[1]*(1-2*pG);
        pC = 1-pA-2*pG;
        break;
    default:
        throw("Unrecognized freq_type in freqsFromParams - can't happen");
    }

    // To MDW, 2017-05-02: please make sure that frequencies are positive!
    // Otherwise, numerical issues will occur.

    bool changed = freq_vec[0]!=pA || freq_vec[1]!=pC || freq_vec[2]!=pG || freq_vec[3]!=pT;
    if (changed) {
        freq_vec[0]=pA;
        freq_vec[1]=pC;
        freq_vec[2]=pG;
        freq_vec[3]=pT;
    }
    return(changed);
}

/*
 * For given freq_type, derives frequency parameters from freq_vec
 * All parameters are in range [0,1] (assuming freq_vec is valid)
 */

void paramsFromFreqs(double *params, double *freq_vec, StateFreqType freq_type) {
    double pA = freq_vec[0]; // These just improve code readability
    double pC = freq_vec[1];
    double pG = freq_vec[2];
//    double pT = freq_vec[3]; // pT is not used below
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        break; // freq_vec never changes
    case FREQ_ESTIMATE:
        params[0]=pA;
        params[1]=pC;
        params[2]=pG;
        break;
    case FREQ_DNA_RY:
        params[0]=2*pA;
        params[1]=2*pC;
        break;
    case FREQ_DNA_WS:
        params[0]=2*pA;
        params[1]=2*pC;
        break;
    case FREQ_DNA_MK:
        params[0]=2*pA;
        params[1]=2*pG;
        break;
    case FREQ_DNA_1112:
        params[0]=3*pA;
        break;
    case FREQ_DNA_1121:
        params[0]=3*pA;
        break;
    case FREQ_DNA_1211:
        params[0]=3*pA;
        break;
    case FREQ_DNA_2111:
        params[0]=3*pC;
        break;
    case FREQ_DNA_1122:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1212:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1221:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1123:
        params[0]=2*pA;
        params[1]=pG/(1-params[0]);
        break;
    case FREQ_DNA_1213:
        params[0]=2*pA;
        params[1]=pC/(1-params[0]);
        break;
    case FREQ_DNA_1231:
        params[0]=2*pA;
        params[1]=pC/(1-params[0]);
        break;
    case FREQ_DNA_2113:
        params[0]=2*pC;
        params[1]=pA/(1-params[0]);
        break;
    case FREQ_DNA_2131:
        params[0]=2*pC;
        params[1]=pA/(1-params[0]);
        break;
    case FREQ_DNA_2311:
        params[0]=2*pG;
        params[1]=pA/(1-params[0]);
        break;
    default:
        throw("Unrecognized freq_type in paramsFromFreqs - can't happen");
    }
}

/*
 * Given a DNA freq_type and a base frequency vector, alter the
 * base freq vector to conform with the constraints of freq_type
 */
void forceFreqsConform(double *base_freq, StateFreqType freq_type) {
    double pA = base_freq[0]; // These just improve code readability
    double pC = base_freq[1];
    double pG = base_freq[2];
    double pT = base_freq[3];
    double scale;
    switch (freq_type) {
    case FREQ_EQUAL:
        // this was already handled, thus not necessary to check here
//        base_freq[0] = base_freq[1] = base_freq[2] = base_freq[3] = 0.25;
//        break;
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
    case FREQ_ESTIMATE:
        break; // any base_freq is legal
    case FREQ_DNA_RY:
        scale = 0.5/(pA+pG);
        base_freq[0] = pA*scale;
        base_freq[2] = pG*scale;
        scale = 0.5/(pC+pT);
        base_freq[1] = pC*scale;
        base_freq[3] = pT*scale;
        break;
    case FREQ_DNA_WS:
        scale = 0.5/(pA+pT);
        base_freq[0] = pA*scale;
        base_freq[3] = pT*scale;
        scale = 0.5/(pC+pG);
        base_freq[1] = pC*scale;
        base_freq[2] = pG*scale;
        break;
    case FREQ_DNA_MK:
        scale = 0.5/(pA+pC);
        base_freq[0] = pA*scale;
        base_freq[1] = pC*scale;
        scale = 0.5/(pG+pT);
        base_freq[2] = pG*scale;
        base_freq[3] = pT*scale;
        break;
    case FREQ_DNA_1112:
        base_freq[0]=base_freq[1]=base_freq[2]=(pA+pC+pG)/3;
        break;
    case FREQ_DNA_1121:
        base_freq[0]=base_freq[1]=base_freq[3]=(pA+pC+pT)/3;
        break;
    case FREQ_DNA_1211:
        base_freq[0]=base_freq[2]=base_freq[3]=(pA+pG+pT)/3;
        break;
    case FREQ_DNA_2111:
        base_freq[1]=base_freq[2]=base_freq[3]=(pC+pG+pT)/3;
        break;
    case FREQ_DNA_1122:
        base_freq[0]=base_freq[1]=(pA+pC)/2;
        base_freq[2]=base_freq[3]=(pG+pT)/2;
        break;
    case FREQ_DNA_1212:
        base_freq[0]=base_freq[2]=(pA+pG)/2;
        base_freq[1]=base_freq[3]=(pC+pT)/2;
        break;
    case FREQ_DNA_1221:
        base_freq[0]=base_freq[3]=(pA+pT)/2;
        base_freq[1]=base_freq[2]=(pC+pG)/2;
        break;
    case FREQ_DNA_1123:
        base_freq[0]=base_freq[1]=(pA+pC)/2;
        break;
    case FREQ_DNA_1213:
        base_freq[0]=base_freq[2]=(pA+pG)/2;
        break;
    case FREQ_DNA_1231:
        base_freq[0]=base_freq[3]=(pA+pT)/2;
        break;
    case FREQ_DNA_2113:
        base_freq[1]=base_freq[2]=(pC+pG)/2;
        break;
    case FREQ_DNA_2131:
        base_freq[1]=base_freq[3]=(pC+pT)/2;
        break;
    case FREQ_DNA_2311:
        base_freq[2]=base_freq[3]=(pG+pT)/2;
        break;
    default:
        throw("Unrecognized freq_type in forceFreqsConform - can't happen");
    }
    ASSERT(base_freq[0]>=0 && base_freq[1]>=0 && base_freq[2]>=0 && base_freq[3]>=0 && fabs(base_freq[0]+base_freq[1]+base_freq[2]+base_freq[3]-1)<1e-7);
}

/*
 * For given freq_type, how many parameters are needed to
 * determine frequenc vector?
 * Currently, this is for DNA StateFreqTypes only.
 */

int nFreqParams(StateFreqType freq_type) {
    switch (freq_type) {
    case FREQ_DNA_1112:
    case FREQ_DNA_1121:
    case FREQ_DNA_1211:
    case FREQ_DNA_2111:
    case FREQ_DNA_1122:
    case FREQ_DNA_1212:
    case FREQ_DNA_1221:
        return(1);
    case FREQ_DNA_RY:
    case FREQ_DNA_WS:
    case FREQ_DNA_MK:
    case FREQ_DNA_1123:
    case FREQ_DNA_1213:
    case FREQ_DNA_1231:
    case FREQ_DNA_2113:
    case FREQ_DNA_2131:
    case FREQ_DNA_2311:
        return(2);
    default:
        return 0; // BQM: don't care about other cases
    }
}

/*
 * For freq_type, and given every base must have frequency >= min_freq, set upper
 * and lower bounds for parameters.
 */
 void setBoundsForFreqType(double *lower_bound,
                           double *upper_bound,
                           bool *bound_check,
                           double min_freq,
                           StateFreqType freq_type) {
    // Sanity check: if min_freq==0, lower_bound=0 and upper_bound=1
    // (except FREQ_ESTIMATE, which follows legacy code way of doing things.)
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        break; // There are no frequency determining parameters
    case FREQ_DNA_1112:
    case FREQ_DNA_1121:
    case FREQ_DNA_1211:
    case FREQ_DNA_2111:
        // one frequency determining parameter
        lower_bound[0] = 3*min_freq;
        upper_bound[0] = 1-min_freq;
        bound_check[0] = true;
        break;
    case FREQ_DNA_1122:
    case FREQ_DNA_1212:
    case FREQ_DNA_1221:
        // one frequency determining parameter
        lower_bound[0] = 2*min_freq;
        upper_bound[0] = 1-2*min_freq;
        bound_check[0] = true;
        break;
    case FREQ_DNA_RY:
    case FREQ_DNA_WS:
    case FREQ_DNA_MK:
        // two frequency determining parameters
        lower_bound[0] = lower_bound[1] = 2*min_freq;
        upper_bound[0] = upper_bound[1] = 1-2*min_freq;
        bound_check[0] = bound_check[1] = true;
    break;
    case FREQ_DNA_1123:
    case FREQ_DNA_1213:
    case FREQ_DNA_1231:
    case FREQ_DNA_2113:
    case FREQ_DNA_2131:
    case FREQ_DNA_2311:
        // two frequency determining parameters
        lower_bound[0] = 2*min_freq;
        upper_bound[0] = 1-2*min_freq;
    lower_bound[1] = min_freq/(1-2*min_freq);
        upper_bound[1] = (1-3*min_freq)/(1-2*min_freq);
        bound_check[0] = bound_check[1] = true;
    break;
        /* NOTE:
     * upper_bound[1] and lower_bound[1] are not perfect. Some in-bounds parameters
         * will give base freqs for '2' or '3' base below minimum. This is
         * the best that can be done without passing min_freq to freqsFromParams
         */
    case FREQ_ESTIMATE:
        lower_bound[0] = lower_bound[1] = lower_bound[2] = min_freq;
        upper_bound[0] = upper_bound[1] = upper_bound[2] = 1;
        bound_check[0] = bound_check[1] = bound_check[2] = false;
        break;
    default:
        throw("Unrecognized freq_type in setBoundsForFreqType - can't happen");
    }
}

/*
 * Note: This function came from model/modelmarkov.cpp
 * For freq_type, return a "+F" string specifying that freq_type.
 * Note not all freq_types accomodated.
 * Inverse of this occurs in ModelFactory::ModelFactory,
 * where +F... suffixes on model names get parsed.
 */
string freqTypeString(StateFreqType freq_type,
                      SeqType seq_type, bool full_str) {
    switch(freq_type) {
    case FREQ_UNKNOWN:    return("");
    case FREQ_USER_DEFINED:
        if (seq_type == SEQ_PROTEIN)
            return "";
        else
            return "+FU";
    case FREQ_EQUAL:
        if (seq_type == SEQ_DNA && !full_str)
            return "";
        else
            return "+FQ";
    case FREQ_EMPIRICAL:  return "+F";
    case FREQ_ESTIMATE:
        return "+FO";
    case FREQ_CODON_1x4:  return("+F1X4");
    case FREQ_CODON_3x4:  return("+F3X4");
    case FREQ_CODON_3x4C: return("+F3X4C");
    case FREQ_MIXTURE:  return(""); // no idea what to do here - MDW
    case FREQ_DNA_RY:   return("+FRY");
    case FREQ_DNA_WS:   return("+FWS");
    case FREQ_DNA_MK:   return("+FMK");
    case FREQ_DNA_1112: return("+F1112");
    case FREQ_DNA_1121: return("+F1121");
    case FREQ_DNA_1211: return("+F1211");
    case FREQ_DNA_2111: return("+F2111");
    case FREQ_DNA_1122: return("+F1122");
    case FREQ_DNA_1212: return("+F1212");
    case FREQ_DNA_1221: return("+F1221");
    case FREQ_DNA_1123: return("+F1123");
    case FREQ_DNA_1213: return("+F1213");
    case FREQ_DNA_1231: return("+F1231");
    case FREQ_DNA_2113: return("+F2113");
    case FREQ_DNA_2131: return("+F2131");
    case FREQ_DNA_2311: return("+F2311");
    default: throw("Unrecognized freq_type in freqTypeString - can't happen");
    }
}
