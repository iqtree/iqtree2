//
//  modelinfo.cpp
//  Created by James Barbetti on 28-Jan-21
//

#include "modelinfo.h"
#include "modelsubst.h"      //for OPEN_BRACKET and CLOSE_BRACKET
#include <utils/tools.h>     //for outError
#include <utils/my_assert.h> //for ASSERT macro

ModelInfoFromName::ModelInfoFromName(std::string name): model_name(name) {}
ModelInfoFromName::ModelInfoFromName(const char* name): model_name(name) {}
ModelInfoFromName::~ModelInfoFromName() {}

namespace {
    bool startsWith(std::string s, const char* front) {
        auto frontLen = strlen(front);
        return (s.substr(0, frontLen) == front);
    }
    bool contains(std::string s, const char* pattern) {
        return s.find(pattern) != std::string::npos;
    }
    std::string::size_type findSubStr(const std::string &name, const std::string sub1,
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
    std::string::size_type findSubStr(const std::string &name, const char * sub1,
                                  const char* sub2) {
        std::string s1(sub1);
        std::string s2(sub2);
        return findSubStr(name, s1, s2);
    }

    std::string::size_type posPOMO(const std::string &model_name) {
        return findSubStr(model_name, "+P", "*P");
    }
};

std::string ModelInfoFromName::getFreeRateParameters(int& num_rate_cats, bool& fused_mix_rate) const {
     string::size_type posR = model_name.find("+R"); // FreeRate model
     string::size_type posR2 = model_name.find("*R"); // FreeRate model

    std::string freerate_params;
    if (posR != string::npos && posR2 != string::npos) {
        cout << "NOTE: both +R and *R were specified, continue with "
        << ((posR < posR2)? model_name.substr(posR,2) : model_name.substr(posR2,2)) << endl;
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
    if (model_name.length() > posR+2+end_pos && model_name[posR+2+end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posR);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        freerate_params = model_name.substr(posR+3+end_pos, close_bracket-posR-3-end_pos).c_str();
    } else if (model_name.length() > posR+2+end_pos && model_name[posR+2+end_pos] != '+') {
        outError("Wrong model name ", model_name);
    }
    return freerate_params;
}

void ModelInfoFromName::getGammaParameters(int& num_rate_cats,
                                           double& gamma_shape) const {
    string::size_type posG = model_name.find("+G");
    string::size_type posG2 = model_name.find("*G");
    if (posG != string::npos && posG2 != string::npos) {
        cout << "NOTE: both +G and *G were specified, continue with "
            << ((posG < posG2)? model_name.substr(posG,2) : model_name.substr(posG2,2)) << endl;
    }
    if (posG2 != string::npos && posG2 < posG) {
        posG = posG2;
    }
    int end_pos = 0;
    if (model_name.length() > posG+2 && isdigit(model_name[posG+2])) {
        num_rate_cats = convert_int(model_name.substr(posG+2).c_str(), end_pos);
        if (num_rate_cats < 1) {
            outError("Wrong number of rate categories");
        }
    }
    if (model_name.length() > posG+2+end_pos && model_name[posG+2+end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posG);
        if (close_bracket == std::string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        gamma_shape = convert_double(model_name.substr(posG+3+end_pos, close_bracket-posG-3-end_pos).c_str());
//            if (gamma_shape < MIN_GAMMA_SHAPE || gamma_shape > MAX_GAMMA_SHAPE) {
//                stringstream str;
//                str << "Gamma shape parameter " << gamma_shape << "out of range ["
//                        << MIN_GAMMA_SHAPE << ',' << MAX_GAMMA_SHAPE << "]" << endl;
//                outError(str.str());
//            }
    } else if (model_name.length() > posG+2+end_pos && model_name[posG+2+end_pos] != '+') {
        outError("Wrong model name ", model_name);
    }
}


std::string ModelInfoFromName::getHeterotachyParameters(bool is_mixture_model,
                                                        int& num_rate_cats,
                                                        bool& fused_mix_rate) const {
    string::size_type posH  = model_name.find("+H"); // heterotachy model
    string::size_type posH2 = model_name.find("*H"); // heterotachy model

    if (posH != string::npos && posH2 != string::npos) {
        cout << "NOTE: both +H and *H were specified, continue with "
            << ((posH < posH2)? model_name.substr(posH,2) : model_name.substr(posH2,2)) << endl;
    }
    if (posH2 != string::npos && posH2 < posH) {
        posH = posH2;
        fused_mix_rate = true;
    }

    std::string heterotachy_params;
    // Heterotachy model
    int end_pos = 0;
    if (model_name.length() > posH+2 && isdigit(model_name[posH+2])) {
        num_rate_cats = convert_int(model_name.substr(posH+2).c_str(), end_pos);
            if (num_rate_cats < 1) outError("Wrong number of rate categories");
    } else {
        if (!is_mixture_model || !fused_mix_rate) {
            outError("Please specify number of heterotachy classes (e.g., +H2)");
        }
    }
    if (model_name.length() > posH+2+end_pos && model_name[posH+2+end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posH);
        if (close_bracket == string::npos)
            outError("Close bracket not found in ", model_name);
        heterotachy_params = model_name.substr(posH+3+end_pos, close_bracket-posH-3-end_pos).c_str();
    } else if (model_name.length() > posH+2+end_pos && model_name[posH+2+end_pos] != '+')
        outError("Wrong model name ", model_name);
    return heterotachy_params;
}

double ModelInfoFromName::getProportionOfInvariantSites() const {
    string::size_type posI = model_name.find("+I");
    if (posI == string::npos) {
        outError("Cannot determine proportion of invariant sites for model ", model_name);
    }
    // invariant site model
    if (model_name.length() > posI+2 && model_name[posI+2] == OPEN_BRACKET) {
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
    } else if (model_name.length() > posI+2 && model_name[posI+2] != '+'
               && model_name[posI+2] != '*') {
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
    if (posG != std::string::npos && posG2 != std::string::npos) {
        stringstream s;
        s << "NOTE: both +G and *G were specified, continue with "
            << ((posG < posG2) ? model_name.substr(posG,2) : model_name.substr(posG2,2));
        outWarning(s.str());
        if (posG2 < posG) {
            //Todo: port: posG = posG2;
            //fused_mix_rate = true;
        }
    }
    return posG != std::string::npos || posG2 != std::string::npos;
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
    || (startsWith(model_name,"MF") && !startsWith(model_name,"MFP"));
}

bool ModelInfoFromName::isPolymorphismAware() const {
    return posPOMO(model_name) != std::string::npos;
}

bool ModelInfoFromName::isWeissAndVonHaeselerTest() const {
    return model_name == "WHTEST";
}

ASCType ModelInfoFromName::extractASCType(std::string& leftover_name) const {
    auto posasc = model_name.find("+ASC_INF");
    if (posasc != std::string::npos) {
        leftover_name = model_name.substr(0, posasc) + model_name.substr(posasc+8);
        return ASC_INFORMATIVE;
    }
    posasc = model_name.find("+ASC_MIS");
    if (posasc != std::string::npos) {
        leftover_name = model_name.substr(0, posasc) + model_name.substr(posasc+8);
        return ASC_VARIANT_MISSING;
    }
    posasc = model_name.find("+ASC");
    ASSERT( posasc != std::string::npos);
    leftover_name = model_name.substr(0, posasc) + model_name.substr(posasc+4);
    return ASC_VARIANT;
}

std::string ModelInfoFromName::extractMixtureModelList(std::string& leftover_name) const {
    ASSERT(startsWith(model_name, "MIX"));
    if (model_name[3] != OPEN_BRACKET)
        outError("Mixture model name must start with 'MIX{'");
    if (model_name.rfind(CLOSE_BRACKET) != model_name.length()-1)
        outError("Close bracket not found at the end of ", model_name);
    leftover_name = "MIX";
    //length is at least 5, since string starts MIX{ and ends },
    //so the next line is safe:
    return model_name.substr(4, model_name.length()-5);
}

std::string ModelInfoFromName::extractPolymorphicHeterozygosity(std::string& leftover_name) const {
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
            pomo_heterozygosity = model_name.substr(p_pos+3,close_bracket-p_pos-3);
            leftover_name = model_name.substr(0, p_pos) + model_name.substr(close_bracket+1);
        }
    }
    else {
        leftover_name = model_name.substr(0, p_pos) + model_name.substr(p_pos + 2);
    }
    return pomo_heterozygosity;
}

void ModelInfoFromName::updateName(const std::string& name) {
    model_name = name;
}

