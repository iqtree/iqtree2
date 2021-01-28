//
//  modelinfo.cpp
//  Created by James Barbetti on 28-Jan-21
//

#include "modelinfo.h"
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

bool ModelInfoFromName::hasRateHeterotachy() const {
    return findSubStr(model_name, "+H", "*H") != std::string::npos;
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

