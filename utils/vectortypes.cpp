#include "vectortypes.h"
#include <string>
#include <sstream>     //for std::stringstream

StrVector::StrVector(size_t size): super(size) {
}

bool StrVector::contains(const char* find_me)        const {
    for (std::string s : *this) {
        if (s==find_me) {
            return true;
        }
    }
    return false;
}

bool StrVector::contains(const std::string& find_me) const {
    for (std::string s : *this) {
        if (s==find_me) {
            return true;
        }
    }
    return false;
}

std::string StrVector::join(const char* sep) const {
    std::stringstream result;
    const char* next_sep = "";
    for (const std::string& s : *this) {
        result << next_sep << s;
        next_sep = sep;
    }
    return result.str();
}

std::string StrVector::join(const std::string& sep) const {
    return join(sep.c_str());
}
