#include "vectortypes.h"
#include <string>

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
