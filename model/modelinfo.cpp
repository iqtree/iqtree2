//
//  modelinfo.cpp
//  Created by James Barbetti on 28-Jan-21
//

#include "modelinfo.h"
#include "modelsubst.h"      //for OPEN_BRACKET and CLOSE_BRACKET
#include "modelfileloader.h"
#include "modeldna.h"        //for ModelDNA
#include "modelexpression.h" //for InterpretedExpression

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
            if (pos1 == std::string::npos) {
                break;
            }
            if (pos1+2 >= name.length() || !isalpha(name[pos1+2])) {
                break;
            }
        }
        for (pos2 = 0; pos2 != std::string::npos; pos2++) {
            pos2 = name.find(sub2, pos2);
            if (pos2 == std::string::npos) {
                break;
            }
            if (pos2+2 >= name.length() ||!isalpha(name[pos2+2])) {
                break;
            }
        }
        if (pos1 != std::string::npos && pos2 != std::string::npos) {
            return std::min(pos1, pos2);
        } else if (pos1 != std::string::npos) {
            return pos1;
        } else {
            return pos2;
        }
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
        std::cout << "NOTE: both +R and *R were specified, continue with "
                  << model_name.substr(posFirst,2) << std::endl;
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
        if (fmix_str[5] != OPEN_BRACKET) {
            outError("Mixture-frequency must start with +FMIX{");
        }
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
            if (close_bracket == string::npos) {
                outError("Close bracket not found in ", fstr);
            }
            if (close_bracket != fstr.length()-1) {
                outError("Wrong close bracket position ", fstr);
            }
            freq_type = FREQ_USER_DEFINED;
            freq_params = fstr.substr(3, close_bracket-3);
        } else if (fstr == "+FC" || fstr == "+Fc" || fstr == "+F") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "empirical," + freq_params;
                optimize_mixmodel_weight = true;
            } else {
                freq_type = FREQ_EMPIRICAL;
            }
        } else if (fstr == "+FU" || fstr == "+Fu") {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with user-defined frequency"
                         " is not allowed");
            }
            else {
                freq_type = FREQ_USER_DEFINED;
            }
        } else if (fstr == "+FQ" || fstr == "+Fq") {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with equal frequency"
                         " is not allowed");
            }
            else {
                freq_type = FREQ_EQUAL;
            }
        } else if (fstr == "+FO" || fstr == "+Fo") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "optimize," + freq_params;
                optimize_mixmodel_weight = true;
            } else {
                freq_type = FREQ_ESTIMATE;
            }
        } else if (fstr == "+F1x4" || fstr == "+F1X4") {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with " + fstr + " is not allowed");
            }
            else {
                freq_type = FREQ_CODON_1x4;
            }
        } else if (fstr == "+F3x4" || fstr == "+F3X4") {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with " + fstr + " is not allowed");
            }
            else {
                freq_type = FREQ_CODON_3x4;
            }
        } else if (fstr == "+F3x4C" || fstr == "+F3x4c" ||
                   fstr == "+F3X4C" || fstr == "+F3X4c") {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with " + fstr + " is not allowed");
            }
            else {
                freq_type = FREQ_CODON_3x4C;
            }
        } else if (fstr == "+FRY") {
            // MDW to Minh: I don't know how these should interact with FREQ_MIXTURE,
            // so as nearly everything else treats it as an error, I do too.
            // BQM answer: that's fine
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with " + fstr + " is not allowed");
            }
            else {
                freq_type = FREQ_DNA_RY;
            }
        } else if (fstr == "+FWS") {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with " + fstr + " is not allowed");
            }
            else {
                freq_type = FREQ_DNA_WS;
            }
        } else if (fstr == "+FMK") {
            if (freq_type == FREQ_MIXTURE) {
                outError("Mixture frequency with " + fstr + " is not allowed");
            }
            else {
                freq_type = FREQ_DNA_MK;
            }
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
        std::cout << "NOTE: both +G and *G were specified, continue with "
                  << model_name.substr(posFirst,2) << std::endl;
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
        //        << MIN_GAMMA_SHAPE << ',' << MAX_GAMMA_SHAPE << "]\n";
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
        std::cout << "NOTE: both +H and *H were specified, continue with "
                  << model_name.substr(posFirst,2) << std::endl;
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
    return model_name.empty() ||
           startsWith(model_name, "TEST") ||
           startsWith(model_name, "MF");
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

YAMLFileParameter::YAMLFileParameter() : is_subscripted(false), value(0.0) {
}

std::string YAMLFileParameter::getSubscriptedVariableName(int subscript) const {
    std::stringstream subscripted_name;
    subscripted_name << name << "(" << subscript << ")";
    return subscripted_name.str();
}

ModelVariable::ModelVariable(): value(0) {
}

ModelVariable::ModelVariable(ModelParameterType t,
                             const ModelParameterRange& r,
                             double v)
    : range(r), type(t), value(v) {
}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile()
    : rate_matrix_rank(0), frequency_type(FREQ_UNKNOWN) {
}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile(const std::string& path)
    :  model_file_path(path), rate_matrix_rank(0), frequency_type(FREQ_UNKNOWN) {
}

void ModelInfoFromYAMLFile::updateName(const std::string& name) {
    model_name = name;
}

void ModelInfoFromYAMLFile::addParameter(const YAMLFileParameter& p) {
    parameters.emplace_back(p);
    if (p.is_subscripted) {
        for (int i=p.minimum_subscript; i<=p.maximum_subscript; ++i) {
            std::string var_name = p.getSubscriptedVariableName(i);
            variables[var_name] = ModelVariable(p.type, p.range, p.value);
        }
    } else {
        variables[p.name] = ModelVariable(p.type, p.range, p.value);
    }
}

bool ModelInfoFromYAMLFile::isFrequencyParameter(const std::string& param_name) const {
    for ( auto p : parameters ) {
        if (string_to_lower(p.name) == string_to_lower(param_name)) {
            return p.type == ModelParameterType::FREQUENCY;
        }
    }
    return false;
}


void ModelInfoFromYAMLFile::setBounds(int param_count, double *lower_bound,
                                      double *upper_bound, bool *bound_check) const {
    int i = 1; //Rate parameter numbering starts at 1, see ModelMarkov
    for ( auto p : parameters ) {
        if (p.type == ModelParameterType::RATE) {
            for (int sub = p.minimum_subscript; sub <= p.maximum_subscript; ++sub) {
                ASSERT( i<= param_count );
                lower_bound[i] = p.range.first;
                upper_bound[i] = p.range.second;
                bound_check[i] = false;
                ++i;
            }
        }
    }
}

void ModelInfoFromYAMLFile::updateVariables(const double* updated_values,
                                            int param_count) {
    int i = 1; //Rate parameter numbering starts at 1, see ModelMarkov
    for ( auto p : parameters ) {
        if (p.type == ModelParameterType::RATE) {
            for (int sub = p.minimum_subscript;
                 sub <= p.maximum_subscript; ++sub) {
                ASSERT( i<= param_count );
                std::string var_name  = p.getSubscriptedVariableName(sub);
                double      var_value = updated_values[i];
                this->variables[var_name].value = var_value;
                ++i;
            }
        }
    }
}



void ModelListFromYAMLFile::loadFromFile (const char* file_path) {
    YAML::Node yaml_model_list = YAML::LoadFile(file_path);
    ModelFileLoader loader(file_path);
    try {
        if (!yaml_model_list.IsSequence()) {
            throw YAML::Exception(yaml_model_list.Mark(),
                                  "list '[...]' expected");
        }
        for (auto node : yaml_model_list) {
            if (!(node["substitutionmodel"])) {
                continue;
            }
            std::string yaml_model_name = node["substitutionmodel"].Scalar();
            std::cout << "Parsing YAML model " << yaml_model_name << std::endl;
            ModelInfoFromYAMLFile &y = models_found[yaml_model_name]
                                     = ModelInfoFromYAMLFile();
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

class YAMLModelDNA: public ModelDNA {
protected:
    ModelInfoFromYAMLFile model_info;
public:
    typedef ModelDNA super;
    YAMLModelDNA(const char *model_name, string model_params,
                 StateFreqType freq, string freq_params,
                 PhyloTree *tree, PhyloTree* report_to_tree,
                 const ModelInfoFromYAMLFile& info)
        : super(model_name, model_params, freq,
                freq_params, tree, report_to_tree),
          model_info(info) {
    }
    virtual void setBounds(double *lower_bound, double *upper_bound,
                            bool *bound_check) {
        //ASSERT(is_reversible && "setBounds should only be called"
        //       "on subclass of ModelMarkov");
        int ndim = getNDim();
        for (int i = 1; i <= ndim; ++i) {
            //std::cout << variables[i] << std::endl;
            lower_bound[i] = MIN_RATE;
            upper_bound[i] = MAX_RATE;
            bound_check[i] = false;
        }
        model_info.setBounds(ndim, lower_bound,
                             upper_bound, bound_check);
    }
    virtual bool getVariables(double *variables) {
        bool changed = false;
        if (num_params > 0) {
            int num_all = static_cast<int>(param_spec.length());
            //if (verbose_mode >= VB_MAX) {
                for (int i = 1; i <= num_params; i++) {
                    std::cout << "  estimated rates[" << i << "] changing from "
                              << rates[i] << " to " << variables[i] << std::endl;
                }
            //}
            for (int i = 0; i < num_all; i++) {
                changed |= (rates[i] != variables[i]);
                rates[i] = variables[i];
            }
        }
        if (freq_type == FREQ_ESTIMATE) {
            int ndim = getNDim();
            auto read_freq = variables+(ndim-num_states+2);
            for (int i=0; i<num_states-1; ++i) {
                if (state_freq[i]!=read_freq[i]) {
                    std::cout << "  estimated freqs[" << i << "] changing from "
                              << state_freq[i] << " to " << read_freq[i] << std::endl;
                    changed = true;
                }
                state_freq[i] = read_freq[i];
                //if (verbose_mode >= VB_MAX) {
                //}
            }
        } else {
            changed |= freqsFromParams(state_freq,variables+num_params+1,freq_type);
        }
        //if (verbose_mode >= VB_MAX) {
        std::cout << std::endl;
        //}
        if (changed) {
            model_info.updateVariables(variables, getNDim());
        }
        return changed;
    }
    virtual void setVariables(double *variables) {
        if (num_params > 0) {
            for (int i=0; i<num_params; ++i) {
                variables[i] = rates[i];
            }
        }
        if (freq_type == FREQ_ESTIMATE) {
            // 2015-09-07: relax the sum of state_freq to be 1,
            // this will be done at the end of optimization
            int ndim = getNDim();
            memcpy(variables+(ndim-num_states+2), state_freq,
                   (num_states-1)*sizeof(double));
        } else {
            paramsFromFreqs(variables+num_params+1,
                            state_freq, freq_type);
        }
    }
};


ModelMarkov* ModelListFromYAMLFile::getModelByName(const char* model_name,   PhyloTree *tree,
                                                   const char* model_params, StateFreqType freq_type,
                                                   const char* freq_params,  PhyloTree* report_to_tree) {
    ModelInfoFromYAMLFile& model_info = models_found[model_name];
    std::cout << "Model Params " << model_params
              << " Freq Params " << freq_params << std::endl;
    
    ModelMarkov* model = nullptr;
    string dummy_rate_params;
    string dummy_freq_params;
    
    if (freq_type == FREQ_UNKNOWN) {
        freq_type = model_info.frequency_type;
    }
    //Todo: other data types, based on ModelBIN,
    //      ModelCodon, ModelPoMo, ModelProtein
    //if (model_info.data_type_name=="DNA") {
    model = new YAMLModelDNA("", dummy_rate_params, freq_type,
                             dummy_freq_params, tree,
                             report_to_tree, model_info);
    
    ASSERT( model_info.rate_matrix_rank == model->num_states);
    
    DoubleVector rates;
    const char* separator = "";
    for (int row=0; row<model_info.rate_matrix_rank; ++row) {
        separator = " ";
        for (int col=0; col<model_info.rate_matrix_rank; ++col) {
            if (col!=row) {
                std::string expr_string = model_info.rate_matrix_expressions[row][col];
                typedef ModelExpression::InterpretedExpression Interpreter;
                Interpreter interpreter(model_info, expr_string);
                double entry = interpreter.evaluate();
                rates.push_back(entry);
                std::cout << separator << entry;
            }
        }
    }
    std::cout << std::endl;
    model->setRateMatrix(rates.data());
    
    //model_parameters = new double [num_params];
    //memset(model_parameters, 0, sizeof(double)*num_params);
    //this->setRates();

    return model;
}
