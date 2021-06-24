//
// modelinfofromname.cpp
// 
/***************************************************************************
 *   Created by James Barbetti on 28-Jan-2021                              *
 *   james_barbetti@yahoo.com                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "modelinfo.h"
#include "modelsubst.h"            //for OPEN_BRACKET and CLOSE_BRACKET
#include <utils/stringfunctions.h> //for convert_double

namespace {
    std::string::size_type findSubStr(const std::string& name,
        const std::string sub1,
        const std::string sub2) {
        std::string::size_type pos1, pos2;
        for (pos1 = 0; pos1 != std::string::npos; pos1++) {
            pos1 = name.find(sub1, pos1);
            if (pos1 == std::string::npos) {
                break;
            }
            if (pos1 + 2 >= name.length() || !isalpha(name[pos1 + 2])) {
                break;
            }
        }
        for (pos2 = 0; pos2 != std::string::npos; pos2++) {
            pos2 = name.find(sub2, pos2);
            if (pos2 == std::string::npos) {
                break;
            }
            if (pos2 + 2 >= name.length() || !isalpha(name[pos2 + 2])) {
                break;
            }
        }
        if (pos1 != std::string::npos && pos2 != std::string::npos) {
            return (pos1 < pos2) ? pos1 : pos2;
        }
        else if (pos1 != std::string::npos) {
            return pos1;
        }
        else {
            return pos2;
        }
    }

    std::string::size_type findSubStr(const std::string& name,
        const char* sub1,
        const char* sub2) {
        std::string s1(sub1);
        std::string s2(sub2);
        return findSubStr(name, s1, s2);
    }

    std::string::size_type posPOMO(const std::string& model_name) {
        return findSubStr(model_name, "+P", "*P");
    }
};

ModelInfoFromName::ModelInfoFromName(std::string name) : model_name(name) {}
ModelInfoFromName::ModelInfoFromName(const char* name) : model_name(name) {}

std::string ModelInfoFromName::getFreeRateParameters
(int& num_rate_cats, bool& fused_mix_rate) const {
    string::size_type posR = model_name.find("+R"); // FreeRate model
    string::size_type posR2 = model_name.find("*R"); // FreeRate model

    std::string freerate_params;
    if (posR != string::npos && posR2 != string::npos) {
        auto posFirst = (posR < posR2) ? posR : posR2;
        std::cout << "NOTE: both +R and *R were specified, continue with "
            << model_name.substr(posFirst, 2) << std::endl;
    }
    if (posR2 != string::npos && posR2 < posR) {
        posR = posR2;
        fused_mix_rate = true;
    }
    // FreeRate model
    int end_pos = 0;
    if (model_name.length() > posR + 2 && isdigit(model_name[posR + 2])) {
        num_rate_cats = convert_int(model_name.substr(posR + 2).c_str(), end_pos);
        if (num_rate_cats < 1) {
            outError("Wrong number of rate categories");
        }
    }
    if (model_name.length() > posR + 2 + end_pos &&
        model_name[posR + 2 + end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posR);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        auto param_start = posR + 3 + end_pos;
        auto param_len = close_bracket - posR - 3 - end_pos;
        freerate_params = model_name.substr(param_start, param_len);
    }
    else if (model_name.length() > posR + 2 + end_pos &&
        model_name[posR + 2 + end_pos] != '+') {
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
        size_t last_pos = freq_str.find_first_of("+*", posfreq + 1);

        if (last_pos == string::npos) {
            fmix_str = freq_str.substr(posfreq);
            freq_str = freq_str.substr(0, posfreq);
        }
        else {
            fmix_str = freq_str.substr(posfreq, last_pos - posfreq);
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
        if (close_bracket != fmix_str.length() - 1) {
            outError("Wrong close bracket position ", fmix_str);
        }
        return fmix_str.substr(6, close_bracket - 6);
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
    if (posfreq == string::npos) {
        posfreq = freq_str.find("+Fo");
        if (posfreq == string::npos) {
            posfreq = freq_str.find("+F");
        }
    }
    if (posfreq != string::npos) {
        string fstr;
        size_t last_pos = freq_str.find_first_of("+*", posfreq + 1);
        if (last_pos == string::npos) {
            fstr = freq_str.substr(posfreq);
            freq_str = freq_str.substr(0, posfreq);
        }
        else {
            fstr = freq_str.substr(posfreq, last_pos - posfreq);
            freq_str = freq_str.substr(0, posfreq)
                + freq_str.substr(last_pos);
        }
        if (fstr.length() > 2 && fstr[2] == OPEN_BRACKET) {
            extractUserDefinedFrequency(fstr, freq_type, freq_params);
        }
        struct Rule {
            const char*   param;
            const char*   mixture_prefix;
            StateFreqType freq_type;
            const char*   mixture_error;
        } rules[] = {
            { "+FC",    "empirical,", StateFreqType::FREQ_EMPIRICAL,    "" }, 
            { "+F",     "empirical,", StateFreqType::FREQ_EMPIRICAL,    "" },
            { "+FU",    "",           StateFreqType::FREQ_USER_DEFINED, "user-defined" },
            { "+FQ",    "",           StateFreqType::FREQ_EQUAL,        "equal"},
            { "+FO",    "optimize,",  StateFreqType::FREQ_ESTIMATE,     ""},
            { "+F1X4",  "",           StateFreqType::FREQ_CODON_1x4,    "1X4"},
            { "+F3X4",  "",           StateFreqType::FREQ_CODON_3x4,    "3X4"},
            { "+F3X4C", "",           StateFreqType::FREQ_CODON_3x4C,   "3X4C"},
            { "+FRY",   "",           StateFreqType::FREQ_DNA_RY,       "RY"},
            { "+FWS",   "",           StateFreqType::FREQ_DNA_WS,       "WS"},
            { "+FMK",   "",           StateFreqType::FREQ_DNA_MK,       "MK"},
        };
        std::string upper_fstr = string_to_upper(fstr);
        for (const Rule& rule : rules) {
            if (upper_fstr==rule.param) {
                if (freq_type == StateFreqType::FREQ_MIXTURE) {
                    if (strlen(rule.mixture_error)!=0) {
                        std::stringstream complaint;
                        complaint << "Mixture frequency with " 
                                  << rule.mixture_error 
                                  << " frequency is not allowed";
                        outError(complaint.str());
                    } else {
                        freq_params = rule.mixture_prefix + freq_params;
                        optimize_mixmodel_weight = true;
                    }
                }
                else {
                    freq_type = rule.freq_type;
                }
                return;
            }
        }
        // might be "+F####" where # are digits
        try {
            freq_type = parseStateFreqDigits(fstr.substr(2));
            // throws an error if not in +F#### format
        }
        catch (...) {
            outError("Unknown state frequency type ", fstr);
        }
        //model_str = model_str.substr(0, posfreq);
    }
}

void ModelInfoFromName::extractUserDefinedFrequency(std::string fstr,
                                                    StateFreqType& freq_type,
                                                    std::string& freq_params) const {
    if (freq_type == StateFreqType::FREQ_MIXTURE) {
        outError("Mixture frequency with user-defined frequency"
            " is not allowed");
    }
    auto close_bracket = fstr.find(CLOSE_BRACKET);
    if (close_bracket == string::npos) {
        outError("Close bracket not found in ", fstr);
    }
    if (close_bracket != fstr.length() - 1) {
        outError("Wrong close bracket position ", fstr);
    }
    freq_type = StateFreqType::FREQ_USER_DEFINED;
    freq_params = fstr.substr(3, close_bracket - 3);
}

void ModelInfoFromName::getGammaParameters(int& num_rate_cats,
    double& gamma_shape) const {
    string::size_type posG = model_name.find("+G");
    string::size_type posG2 = model_name.find("*G");
    if (posG != string::npos && posG2 != string::npos) {
        auto posFirst = (posG < posG2) ? posG : posG2;
        std::cout << "NOTE: both +G and *G were specified, continue with "
            << model_name.substr(posFirst, 2) << std::endl;
    }
    if (posG2 != string::npos && posG2 < posG) {
        posG = posG2;
    }
    int end_pos = 0;
    if (model_name.length() > posG + 2 &&
        isdigit(model_name[posG + 2])) {
        auto rest = model_name.substr(posG + 2);
        num_rate_cats = convert_int(rest.c_str(), end_pos);
        if (num_rate_cats < 1) {
            outError("Wrong number of rate categories");
        }
    }
    if (model_name.length() > posG + 2 + end_pos &&
        model_name[posG + 2 + end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posG);
        if (close_bracket == std::string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        auto gamma_str = model_name.substr(posG + 3 + end_pos,
            close_bracket - posG - 3 - end_pos);
        gamma_shape = convert_double(gamma_str.c_str());
        //if (gamma_shape < MIN_GAMMA_SHAPE ||
        //    gamma_shape > MAX_GAMMA_SHAPE) {
        //    stringstream str;
        //    str << "Gamma shape parameter " << gamma_shape << "out of range ["
        //        << MIN_GAMMA_SHAPE << ',' << MAX_GAMMA_SHAPE << "]\n";
        //    outError(str.str());
        //}
    }
    else if (model_name.length() > posG + 2 + end_pos &&
        model_name[posG + 2 + end_pos] != '+') {
        outError("Wrong model name ", model_name);
    }
}

int ModelInfoFromName::getNumberOfCategories(std::string model_name, 
                                             string::size_type posH, 
                                             bool is_mixture_model,
                                             bool fused_mix_rate,
                                             int &end_pos) const {
    int number_of_categories = 0;
    if (model_name.length() > posH + 2 && isdigit(model_name[posH + 2])) {
        auto rest = model_name.substr(posH + 2);
        number_of_categories = convert_int(rest.c_str(), end_pos);
        if (number_of_categories < 1) {
            std::stringstream complaint;
            complaint << "Wrong number of rate categories"
                      << " (" << number_of_categories << ")";
            outError(complaint.str());
        }
    }
    else if (!is_mixture_model || !fused_mix_rate) {
        outError("Please specify number of heterotachy classes (e.g., +H2)");
    }
    return number_of_categories;
}

std::string ModelInfoFromName::getOtherHeterotachyParameters
                ( std::string model_name, string::size_type posH, 
                  int end_pos) const {
    std::string heterotachy_params;
    if (model_name.length() > posH + 2 + end_pos &&
        model_name[posH + 2 + end_pos] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posH);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        auto hetero_start = posH + 3 + end_pos;
        auto hetero_len = close_bracket - posH - 3 - end_pos;
        heterotachy_params = model_name.substr(hetero_start, hetero_len);
    }
    else if (model_name.length() > posH + 2 + end_pos &&
        model_name[posH + 2 + end_pos] != '+') {
        outError("Wrong model name ", model_name);
    }
    return heterotachy_params;
}

std::string ModelInfoFromName::getHeterotachyParameters
(bool is_mixture_model, int& num_rate_cats,
    bool& fused_mix_rate) const {
    string::size_type posH  = model_name.find("+H"); // heterotachy model
    string::size_type posH2 = model_name.find("*H"); // heterotachy model

    if (posH != string::npos && posH2 != string::npos) {
        auto posFirst = (posH < posH2) ? posH : posH2;
        std::cout << "NOTE: both +H and *H were specified, continue with "
                  << model_name.substr(posFirst, 2) << std::endl;
    }
    if (posH2 != string::npos && posH2 < posH) {
        posH = posH2;
        fused_mix_rate = true;
    }


    int end_pos = 0;
    // Heterotachy model
    num_rate_cats = getNumberOfCategories(model_name, posH, is_mixture_model,
                                          fused_mix_rate, end_pos);
    return getOtherHeterotachyParameters(model_name, posH, end_pos);

}

double ModelInfoFromName::getProportionOfInvariantSites() const {
    string::size_type posI = model_name.find("+I");
    if (posI == string::npos) {
        outError("Cannot determine proportion of invariant sites"
            " for model ", model_name);
    }
    // invariant site model
    if (model_name.length() > posI + 2 &&
        model_name[posI + 2] == OPEN_BRACKET) {
        auto close_bracket = model_name.find(CLOSE_BRACKET, posI);
        if (close_bracket == string::npos) {
            outError("Close bracket not found in ", model_name);
        }
        std::string num = model_name.substr(posI + 3, close_bracket - posI - 3);
        double p_invar_sites = convert_double(num.c_str());
        if (p_invar_sites < 0 || p_invar_sites >= 1) {
            outError("p_invar must be in [0,1)");
        }
        return p_invar_sites;
    }
    else if (model_name.length() > posI + 2 &&
        model_name[posI + 2] != '+' &&
        model_name[posI + 2] != '*') {
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
    std::string::size_type posG = model_name.find("+G");
    std::string::size_type posG2 = model_name.find("*G");
    if (posG != std::string::npos &&
        posG2 != std::string::npos) {
        stringstream s;
        auto posFirst = (posG < posG2) ? posG : posG2;
        s << "NOTE: both +G and *G were specified, continue with "
            << model_name.substr(posFirst, 2);
        outWarning(s.str());
        if (posG2 < posG) {
            //Todo: port: posG = posG2;
            //fused_mix_rate = true;
        }
    }
    return posG != std::string::npos ||
        posG2 != std::string::npos;
}

bool ModelInfoFromName::isInvariantModel() const {
    return model_name.find("+I") != std::string::npos;
}

bool ModelInfoFromName::isKategoryModel() const {
    return model_name.find("+K") != string::npos;
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
    return contains(model_name, "ONLY")
        || (startsWith(model_name, "MF") &&
            !startsWith(model_name, "MFP"));
}

bool ModelInfoFromName::isPolymorphismAware() const {
    return posPOMO(model_name) != std::string::npos;
}

bool ModelInfoFromName::isWeissAndVonHaeselerTest() const {
    return model_name == "WHTEST";
}

int  ModelInfoFromName::getKategoryRateCount(int rate_count, int min_count) const {
    auto posX = model_name.find("+K");
    if (model_name.length() > posX+2 && isdigit(model_name[posX+2])) {
        rate_count = convert_int(model_name.substr(posX+2).c_str());
        if (rate_count < min_count) {
            outError("Wrong number of rate categories");
        }
    }
    return rate_count;
}

ASCType ModelInfoFromName::extractASCType
(std::string& leftover_name) const {
    auto posasc = model_name.find("+ASC_INF");
    if (posasc != std::string::npos) {
        leftover_name = model_name.substr(0, posasc)
            + model_name.substr(posasc + 8);
        return ASCType::ASC_INFORMATIVE;
    }
    posasc = model_name.find("+ASC_MIS");
    if (posasc != std::string::npos) {
        leftover_name = model_name.substr(0, posasc)
            + model_name.substr(posasc + 8);
        return ASCType::ASC_VARIANT_MISSING;
    }
    posasc = model_name.find("+ASC");
    ASSERT(posasc != std::string::npos);
    leftover_name = model_name.substr(0, posasc)
        + model_name.substr(posasc + 4);
    return ASCType::ASC_VARIANT;
}

std::string ModelInfoFromName::extractMixtureModelList
(std::string& leftover_name) const {
    ASSERT(startsWith(model_name, "MIX"));
    if (model_name[3] != OPEN_BRACKET) {
        outError("Mixture model name must start with 'MIX{'");
    }
    if (model_name.rfind(CLOSE_BRACKET) != model_name.length() - 1) {
        outError("Close bracket not found at the end of ", model_name);
    }
    leftover_name = "MIX";
    //length is at least 5, since string starts MIX{ and ends },
    //so the next line is safe:
    return model_name.substr(4, model_name.length() - 5);
}

std::string ModelInfoFromName::extractPolymorphicHeterozygosity
(std::string& leftover_name) const {
    auto p_pos = findSubStr(model_name, "+P", "*P");
    ASSERT(p_pos != std::string::npos);
    leftover_name = model_name;
    std::string pomo_heterozygosity;
    if (model_name[p_pos + 2] == '{') {
        std::string::size_type close_bracket = model_name.find("}", p_pos);
        if (close_bracket == std::string::npos) {
            std::cout << "Model string: " << model_name << std::endl;
            outError("No closing bracket in PoMo parameters.");
        }
        else {
            auto het_start = p_pos + 3;
            auto het_len = close_bracket - p_pos - 3;
            pomo_heterozygosity = model_name.substr(het_start, het_len);
            leftover_name = model_name.substr(0, p_pos)
                + model_name.substr(close_bracket + 1);
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
