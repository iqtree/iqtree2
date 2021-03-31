//
//  modelinfo.cpp
//  Created by James Barbetti on 28-Jan-21
//

#include "modelinfo.h"
#include "modelsubst.h"      //for OPEN_BRACKET and CLOSE_BRACKET
#include <utils/my_assert.h> //for ASSERT macro
#include <utils/stringfunctions.h> //for convert_int
#include <utils/tools.h>     //for outError

ModelInfoFromName::ModelInfoFromName(std::string name): model_name(name) {}
ModelInfoFromName::ModelInfoFromName(const char* name): model_name(name) {}

namespace {
    bool startsWith(std::string s, const char* front) {
        auto frontLen = strlen(front);
        return (s.substr(0, frontLen) == front);
    }
    bool contains(std::string s, const char* pattern) {
        return s.find(pattern) != std::string::npos;
    }
    std::string::size_type findSubStr(const std::string &name,
                                      const std::string sub1,
                                      const std::string sub2) {
        std::string::size_type pos1, pos2;
        for (pos1 = 0; pos1 != std::string::npos; pos1++) {
            pos1 = name.find(sub1, pos1);
            if (pos1 == std::string::npos)
                break;
            if (pos1+2 >= name.length() || !isalpha(name[pos1+2])) {
                break;
            }
        }
        for (pos2 = 0; pos2 != std::string::npos; pos2++) {
            pos2 = name.find(sub2, pos2);
            if (pos2 == std::string::npos)
                break;
            if (pos2+2 >= name.length() ||!isalpha(name[pos2+2]))
                break;
        }
        if (pos1 != std::string::npos && pos2 != std::string::npos) {
            return std::min(pos1, pos2);
        } else if (pos1 != std::string::npos)
            return pos1;
        else
            return pos2;
    }
    std::string::size_type findSubStr(const std::string &name,
                                      const char * sub1,
                                      const char* sub2) {
        std::string s1(sub1);
        std::string s2(sub2);
        return findSubStr(name, s1, s2);
    }

    std::string::size_type posPOMO(const std::string &model_name) {
        return findSubStr(model_name, "+P", "*P");
    }
};

std::string ModelInfoFromName::getFreeRateParameters
    (int& num_rate_cats, bool& fused_mix_rate) const {
     string::size_type posR = model_name.find("+R"); // FreeRate model
     string::size_type posR2 = model_name.find("*R"); // FreeRate model

    std::string freerate_params;
    if (posR != string::npos && posR2 != string::npos) {
        auto posFirst = (posR < posR2) ? posR : posR2;
        cout << "NOTE: both +R and *R were specified, continue with "
             << model_name.substr(posFirst,2) << endl;
    }
    if (posR2 != string::npos && posR2 < posR) {
        posR = posR2;
        fused_mix_rate = true;
    }
    // FreeRate model
    int end_pos = 0;
    if (model_name.length() > posR+2 && isdigit(model_name[posR+2])) {
        num_rate_cats = convert_int(model_name.substr(posR+2).c_str(), end_pos);
        if (num_rate_cats < 1) {
            outError("Wrong number of rate categories");
        }
    }
    if (model_name.length() > posR+2+end_pos &&
        model_name[posR+2+end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posR);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        auto param_start = posR+3+end_pos;
        auto param_len   = close_bracket-posR-3-end_pos;
        freerate_params  = model_name.substr(param_start, param_len);
    } else if (model_name.length() > posR+2+end_pos &&
               model_name[posR+2+end_pos] != '+') {
        outError("Wrong model name ", model_name);
    }
    return freerate_params;
}

std::string ModelInfoFromName::getFrequencyMixtureParams
    (std::string& freq_str) const {
    // first handle mixture frequency
    std::string::size_type posfreq = model_name.find("+FMIX");
    size_t close_bracket;
    freq_str = model_name;

    if (posfreq != string::npos) {
        string fmix_str;
        size_t last_pos = freq_str.find_first_of("+*", posfreq+1);

        if (last_pos == string::npos) {
            fmix_str = freq_str.substr(posfreq);
            freq_str = freq_str.substr(0, posfreq);
        } else {
            fmix_str = freq_str.substr(posfreq, last_pos-posfreq);
            freq_str = freq_str.substr(0, posfreq)
                     + freq_str.substr(last_pos);
        }
        if (fmix_str[5] != OPEN_BRACKET)
            outError("Mixture-frequency must start with +FMIX{");
        close_bracket = fmix_str.find(CLOSE_BRACKET);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", fmix_str);
        }
        if (close_bracket != fmix_str.length()-1) {
            outError("Wrong close bracket position ", fmix_str);
        }
        return fmix_str.substr(6, close_bracket-6);
    }
    return "";
}

void ModelInfoFromName::getFrequencyOptions(std::string& freq_str,
                                            StateFreqType& freq_type,
                                            std::string& freq_params,
                                            bool& optimize_mixmodel_weight) const {
    // then normal frequency
    freq_str = model_name;
    std::string::size_type posfreq = freq_str.find("+FO");
    if (posfreq== string::npos) {
        posfreq = freq_str.find("+Fo");
        if (posfreq == string::npos) {
            posfreq = freq_str.find("+F");
        }
    }
    if (posfreq != string::npos) {
        string fstr;
        size_t last_pos = freq_str.find_first_of("+*", posfreq+1);
        if (last_pos == string::npos) {
            fstr = freq_str.substr(posfreq);
            freq_str = freq_str.substr(0, posfreq);
        } else {
            fstr = freq_str.substr(posfreq, last_pos-posfreq);
            freq_str = freq_str.substr(0, posfreq)
                     + freq_str.substr(last_pos);
        }
        if (fstr.length() > 2 && fstr[2] == OPEN_BRACKET) {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with user-defined frequency"
                         " is not allowed");
            }
            auto close_bracket = fstr.find(CLOSE_BRACKET);
            if (close_bracket == string::npos)
                outError("Close bracket not found in ", fstr);
            if (close_bracket != fstr.length()-1)
                outError("Wrong close bracket position ", fstr);
            freq_type = FREQ_USER_DEFINED;
            freq_params = fstr.substr(3, close_bracket-3);
        } else if (fstr == "+FC" || fstr == "+Fc" || fstr == "+F") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "empirical," + freq_params;
                optimize_mixmodel_weight = true;
            } else
                freq_type = FREQ_EMPIRICAL;
        } else if (fstr == "+FU" || fstr == "+Fu") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with user-defined frequency"
                         " is not allowed");
            else
                freq_type = FREQ_USER_DEFINED;
        } else if (fstr == "+FQ" || fstr == "+Fq") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with equal frequency"
                         " is not allowed");
            else
            freq_type = FREQ_EQUAL;
        } else if (fstr == "+FO" || fstr == "+Fo") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "optimize," + freq_params;
                optimize_mixmodel_weight = true;
            } else
                freq_type = FREQ_ESTIMATE;
    } else if (fstr == "+F1x4" || fstr == "+F1X4") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_CODON_1x4;
        } else if (fstr == "+F3x4" || fstr == "+F3X4") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_CODON_3x4;
        } else if (fstr == "+F3x4C" || fstr == "+F3x4c" ||
                   fstr == "+F3X4C" || fstr == "+F3X4c") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_CODON_3x4C;
        } else if (fstr == "+FRY") {
        // MDW to Minh: I don't know how these should interact with FREQ_MIXTURE,
        // so as nearly everything else treats it as an error, I do too.
        // BQM answer: that's fine
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_DNA_RY;
        } else if (fstr == "+FWS") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_DNA_WS;
        } else if (fstr == "+FMK") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_DNA_MK;
        } else {
            // might be "+F####" where # are digits
            try {
                freq_type = parseStateFreqDigits(fstr.substr(2));
                // throws an error if not in +F#### format
            } catch (...) {
                outError("Unknown state frequency type ",fstr);
            }
        }
        //model_str = model_str.substr(0, posfreq);
    }
}

void ModelInfoFromName::getGammaParameters(int& num_rate_cats,
                                           double& gamma_shape) const {
    string::size_type posG = model_name.find("+G");
    string::size_type posG2 = model_name.find("*G");
    if (posG != string::npos && posG2 != string::npos) {
        auto posFirst = (posG < posG2) ? posG : posG2;
        cout << "NOTE: both +G and *G were specified, continue with "
             << model_name.substr(posFirst,2) << endl;
    }
    if (posG2 != string::npos && posG2 < posG) {
        posG = posG2;
    }
    int end_pos = 0;
    if (model_name.length() > posG+2 &&
        isdigit(model_name[posG+2])) {
        auto rest     = model_name.substr(posG+2);
        num_rate_cats = convert_int(rest.c_str(), end_pos);
        if (num_rate_cats < 1) {
            outError("Wrong number of rate categories");
        }
    }
    if (model_name.length() > posG+2+end_pos &&
        model_name[posG+2+end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posG);
        if (close_bracket == std::string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        gamma_shape = convert_double(model_name.substr(posG+3+end_pos, close_bracket-posG-3-end_pos).c_str());
        //if (gamma_shape < MIN_GAMMA_SHAPE ||
        //    gamma_shape > MAX_GAMMA_SHAPE) {
        //    stringstream str;
        //    str << "Gamma shape parameter " << gamma_shape << "out of range ["
        //        << MIN_GAMMA_SHAPE << ',' << MAX_GAMMA_SHAPE << "]" << endl;
        //    outError(str.str());
        //}
    } else if (model_name.length() > posG+2+end_pos &&
               model_name[posG+2+end_pos] != '+') {
        outError("Wrong model name ", model_name);
    }
}


std::string ModelInfoFromName::getHeterotachyParameters
    (bool is_mixture_model, int& num_rate_cats,
     bool& fused_mix_rate) const {
    string::size_type posH  = model_name.find("+H"); // heterotachy model
    string::size_type posH2 = model_name.find("*H"); // heterotachy model

    if (posH != string::npos && posH2 != string::npos) {
        auto posFirst = (posH < posH2) ? posH : posH2;
        cout << "NOTE: both +H and *H were specified, continue with "
             << model_name.substr(posFirst,2) << endl;
    }
    if (posH2 != string::npos && posH2 < posH) {
        posH = posH2;
        fused_mix_rate = true;
    }
    std::string heterotachy_params;
    // Heterotachy model
    int end_pos = 0;
    if (model_name.length() > posH+2 && isdigit(model_name[posH+2])) {
        auto rest = model_name.substr(posH+2);
        num_rate_cats = convert_int(rest.c_str(), end_pos);
        if (num_rate_cats < 1) {
            outError("Wrong number of rate categories");
        }
    } else if (!is_mixture_model || !fused_mix_rate) {
        outError("Please specify number of heterotachy classes (e.g., +H2)");
    }
    if (model_name.length() > posH+2+end_pos &&
        model_name[posH+2+end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posH);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        auto hetero_start  = posH+3+end_pos;
        auto hetero_len    = close_bracket-posH-3-end_pos;
        heterotachy_params = model_name.substr(hetero_start, hetero_len);
    } else if (model_name.length() > posH+2+end_pos &&
               model_name[posH+2+end_pos] != '+') {
        outError("Wrong model name ", model_name);
    }
    return heterotachy_params;
}

double ModelInfoFromName::getProportionOfInvariantSites() const {
    string::size_type posI = model_name.find("+I");
    if (posI == string::npos) {
        outError("Cannot determine proportion of invariant sites"
                 " for model ", model_name);
    }
    // invariant site model
    if (model_name.length() > posI+2 &&
        model_name[posI+2] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posI);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        std::string num = model_name.substr(posI+3, close_bracket-posI-3);
        double p_invar_sites = convert_double(num.c_str());
        if (p_invar_sites < 0 || p_invar_sites >= 1) {
            outError("p_invar must be in [0,1)");
        }
        return p_invar_sites;
    } else if (model_name.length() > posI+2 &&
               model_name[posI+2] != '+' &&
               model_name[posI+2] != '*') {
        outError("Wrong model name ", model_name);
    }
    return 0;
}


bool ModelInfoFromName::hasAscertainmentBiasCorrection() const {
    return model_name.find("+ASC") != std::string::npos;
}

bool ModelInfoFromName::hasRateHeterotachy() const {
    return findSubStr(model_name, "+H", "*H") != std::string::npos;
}

bool ModelInfoFromName::isFreeRate() const {
    return findSubStr(model_name, "+R", "*R") != std::string::npos;
}

bool ModelInfoFromName::isFrequencyMixture() const {
    return model_name.find("+FMIX") != string::npos;
}

bool ModelInfoFromName::isGammaModel() const {
    std::string::size_type posG  = model_name.find("+G");
    std::string::size_type posG2 = model_name.find("*G");
    if (posG != std::string::npos &&
        posG2 != std::string::npos) {
        stringstream s;
        auto posFirst = (posG < posG2) ? posG : posG2;
        s << "NOTE: both +G and *G were specified, continue with "
          << model_name.substr(posFirst,2);
        outWarning(s.str());
        if (posG2 < posG) {
            //Todo: port: posG = posG2;
            //fused_mix_rate = true;
        }
    }
    return posG  != std::string::npos ||
           posG2 != std::string::npos;
}

bool ModelInfoFromName::isInvariantModel() const {
    return model_name.find("+I") != std::string::npos;
}

bool ModelInfoFromName::isMixtureModel() const {
    return startsWith(model_name, "MIX");
}

bool ModelInfoFromName::isModelFinder() const {
    return model_name.empty() || startsWith(model_name, "TEST")
    || startsWith(model_name, "MF");
}

bool ModelInfoFromName::isModelFinderOnly() const {
    return contains(model_name,"ONLY")
    || (startsWith(model_name,"MF") &&
        !startsWith(model_name,"MFP"));
}

bool ModelInfoFromName::isPolymorphismAware() const {
    return posPOMO(model_name) != std::string::npos;
}

bool ModelInfoFromName::isWeissAndVonHaeselerTest() const {
    return model_name == "WHTEST";
}

ASCType ModelInfoFromName::extractASCType
    (std::string& leftover_name) const {
    auto posasc = model_name.find("+ASC_INF");
    if (posasc != std::string::npos) {
        leftover_name = model_name.substr(0, posasc)
                      + model_name.substr(posasc+8);
        return ASC_INFORMATIVE;
    }
    posasc = model_name.find("+ASC_MIS");
    if (posasc != std::string::npos) {
        leftover_name = model_name.substr(0, posasc)
                      + model_name.substr(posasc+8);
        return ASC_VARIANT_MISSING;
    }
    posasc = model_name.find("+ASC");
    ASSERT( posasc != std::string::npos);
    leftover_name = model_name.substr(0, posasc)
                  + model_name.substr(posasc+4);
    return ASC_VARIANT;
}

std::string ModelInfoFromName::extractMixtureModelList
    (std::string& leftover_name) const {
    ASSERT(startsWith(model_name, "MIX"));
    if (model_name[3] != OPEN_BRACKET) {
        outError("Mixture model name must start with 'MIX{'");
    }
    if (model_name.rfind(CLOSE_BRACKET) != model_name.length()-1) {
        outError("Close bracket not found at the end of ", model_name);
    }
    leftover_name = "MIX";
    //length is at least 5, since string starts MIX{ and ends },
    //so the next line is safe:
    return model_name.substr(4, model_name.length()-5);
}

std::string ModelInfoFromName::extractPolymorphicHeterozygosity
            (std::string& leftover_name) const {
    auto p_pos = findSubStr(model_name, "+P", "*P");
    ASSERT( p_pos != std::string::npos );
    leftover_name = model_name;
    std::string pomo_heterozygosity;
    if (model_name[p_pos+2] == '{') {
        std::string::size_type close_bracket = model_name.find("}", p_pos);
        if (close_bracket == std::string::npos) {
            std::cout << "Model string: " << model_name << std::endl;
            outError("No closing bracket in PoMo parameters.");
        }
        else {
            auto het_start = p_pos+3;
            auto het_len   = close_bracket-p_pos-3;
            pomo_heterozygosity = model_name.substr(het_start,het_len);
            leftover_name = model_name.substr(0, p_pos)
                          + model_name.substr(close_bracket+1);
        }
    }
    else {
        leftover_name = model_name.substr(0, p_pos)
                      + model_name.substr(p_pos + 2);
    }
    return pomo_heterozygosity;
}

void ModelInfoFromName::updateName(const std::string& name) {
    model_name = name;
}


YAMLFileParameter::YAMLFileParameter() : value(0.0) {
}

class YAMLFileLoader {
public:
    
    const char* file_path;
    const std::string model_name;
   
    YAMLFileLoader(const char* path): file_path(path) {
    }
        
    std::string stringScalar(const YAML::Node& node,
                                    const char* key) {
        auto cite = node[key];
        return cite ? cite.Scalar() : "";
    }
    
    bool booleanScalar(const YAML::Node& node,
                              const char* key) {
        std::string s = string_to_lower(stringScalar(node, key));
        return s == "true" || s == "yes" || s == "t" || s == "y" || s == "1";
    }
    
    void complainIfNot(bool check_me,
                              std::string error_message) {
        if (!check_me) {
            outError(error_message);
        }
    }
    
    double toDouble(const YAML::Node& i, double default_val) {
        if (!i.IsScalar()) {
            return default_val;
        }
        std::string double_string = i.Scalar();
        return convert_double_nothrow(double_string.c_str(), default_val);
    }
    
    ParameterRange parseRange(const YAML::Node& node, const char* key) {
        ParameterRange range;
        auto r = node[key];
        if (!r) {
            return range;
        }
        //Todo: what if range is a string?
        if (r.IsSequence()) {
            int ix = 0;
            for (auto i : r) {
                if (ix==0) {
                    range.first  = toDouble(i, 0);
                } else if (ix==1) {
                    range.second = toDouble(i, 0);
                } else {
                    //Throw: Range ought to be a two-value sequence
                    //       Three values... one too many.
                    //(for now, 3rd and later items in sequence
                    // are just ignored).
                }
                ++ix;
            }
            range.is_set = ( 1 < ix );
        }
        return range;
    }
    
    void parseYAMLModelParameters(const YAML::Node& params,
                                  ModelInfoFromYAMLFile& info) {
        //
        //Assumes: params is a sequence of parameter declarations
        //
        for (auto param: params) {
            YAMLFileParameter p;
            p.name      = stringScalar(param, "name");
            p.type_name = string_to_lower(stringScalar(param, "type"));
            double dv   = 0.0; //default initial value
            if (p.type_name=="frequency") {
                dv = 0.25;  //Todo: should be 1.0 divided by number of states
                            //determined from the data type (info.data_type_name ?)
                            //Or 1 divided by the number of parameters.
            } else if (p.type_name=="rate") {
                dv = 1.0;
            } else if (p.type_name=="weight") {
                dv = 0.5;    //Todo: Should be 1.0 divided by # of parameters
            }
            p.range  = parseRange  (param, "range");
            std::string value_string = stringScalar(param, "initValue");
            p.value  = convert_double_nothrow(value_string.c_str(), dv);
            std::cout << "Parsed parameter " << p.name
                      << " of type " << p.type_name
                      << ", with range " << p.range.first
                      << " to " << p.range.second
                      << ", and initial value " << p.value << std::endl;
            info.parameters.emplace_back(p);
        }
    }
    
    template <class S>
    void dumpRateMatrixTo(const ModelInfoFromYAMLFile& info, S &out) {
        std::stringstream dump;
        for (auto r : info.rate_matrix_expressions) {
            for (auto c: r) {
                dump << c << " : ";
            }
            dump << "\n";
        }
        out << "Rate matrix for " << info.model_name << " is...\n" << dump.str();
    }
    
    /* Example of a rate matrix
     
     rateMatrix:
     - [      ,     r(1)*r(2), r(2)*f(3), r(3)*f(4) ]
     - [ r(1)*f(1),          , r(4)*f(3), r(5)*f(4) ]
     - [ r(2)*f(1), r(4)*f(2),          , f(4)      ]
     - [ r(3)*f(1), r(5)*f(2), f(3)     ,           ]

     */
    
    void parseRateMatrix(const YAML::Node& rate_matrix,
                         ModelInfoFromYAMLFile& info) {
        //Assumes rate_matrix is a sequence (of rows)
        for (auto row : rate_matrix)
        {
            ++info.rate_matrix_rank;
            std::stringstream s;
            s << "Row " << info.rate_matrix_rank << " of rate matrix "
              << " for model " << info.model_name
              << " in " << info.model_file_path;
            std::string context = s.str();
            complainIfNot(row.IsSequence(),
                          context + " is not a sequence" );
            StrVector expr_row;
            for (auto col : row) {
                if (col.IsNull()) {
                    expr_row.emplace_back("");
                    //Gets whatever hasn't been assigned in this row
                }
                else if (!col.IsScalar()) {
                    std::stringstream s2;
                    s2 << "Column " << (expr_row.size()+1)
                       << " of " << context << " is not a scalar";
                    outError(s2.str());
                } else {
                    expr_row.emplace_back(col.Scalar());
                }
            }
            info.rate_matrix_expressions.emplace_back(expr_row);
        }
        dumpRateMatrixTo(info, cout);
    }
        
    void parseYAMLSubstitutionModel(const YAML::Node& substitution_model,
                                    const std::string& name_of_model,
                                    ModelInfoFromYAMLFile& info) {
        info.model_file_path = file_path;
        info.model_name      = name_of_model;
        info.citation        = stringScalar(substitution_model,  "citation");
        info.reversible      = booleanScalar(substitution_model, "reversible");
        //
        //Todo: read off the datatype (if there is one).
        //
        info.data_type_name  = stringScalar(substitution_model,  "datatype");
        
        //
        //Todo: extract other information from the subsstitution model.
        //      Such as parameters and rate matrices and so forth
        //
        auto params = substitution_model["parameters"];
        if (params) {
            complainIfNot(params.IsSequence(),
                          "Parameters of model " + model_name +
                          " in file " + file_path + " not a sequence");
            parseYAMLModelParameters(params, info);
        }
        
        //
        //Todo: It should be okay for a model *not* to specify a rateMatrix,
        //      so long as it subclasses another model that *does* specify one.
        //
        auto rateMatrix = substitution_model["rateMatrix"];
        complainIfNot(rateMatrix, "Model " + model_name +
                      " in file " + file_path +
                      " does not specify a rateMatrix" );
        parseRateMatrix(rateMatrix, info);
        
        //
        //Todo: How do you specify paramters
        //
        auto stateFrequency = substitution_model["stateFrequency"];
        if (stateFrequency) {
            //Check that dimension of the specified parameter is the
            //same as the rank of the rate matrix (it must be!).
        } else {
            //If we have parameters with a type of frequency, we're all good.
        }
    }
};

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile(): rate_matrix_rank(0) {
    
}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile(const std::string& path)
    : model_file_path(path) {
}

void ModelInfoFromYAMLFile::updateName(const std::string& name) {
    model_name = name;
}

void ModelListFromYAMLFile::loadFromFile (const char* file_path) {
    YAML::Node yaml_model_list = YAML::LoadFile(file_path);
    YAMLFileLoader loader(file_path);
    try {
        if (!yaml_model_list.IsSequence()) {
            throw YAML::Exception(yaml_model_list.Mark(), "list '[...]' expected");
        }
        for (auto node : yaml_model_list) {
            if (!(node["substitutionmodel"])) {
                continue;
            }
            std::string yaml_model_name = node["substitutionmodel"].Scalar();
            std::cout << "Parsing YAML model " << yaml_model_name << std::endl;
            ModelInfoFromYAMLFile &y = models_found[yaml_model_name] = ModelInfoFromYAMLFile();
            loader.parseYAMLSubstitutionModel(node, yaml_model_name, y);
        }
    }
    catch (YAML::Exception &e) {
        outError(e.what());
    }
}

bool ModelListFromYAMLFile::isModelNameRecognized (const char* model_name) {
    return models_found.find(std::string(model_name)) != models_found.end();
}

ModelMarkov* ModelListFromYAMLFile::getModelByName
    (const char* model_name,   PhyloTree *tree,
     const char* model_params, StateFreqType freq_type,
     const char* freq_params,  PhyloTree* report_to_tree) {
    FUNCTION_NOT_IMPLEMENTED;
    return nullptr;
}
