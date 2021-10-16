#include "nametoidmap.h"

NameToIDMap::NameToIDMap(const StrVector& src) {
    int count = static_cast<int>(src.size());
    for (int id = 0 ; id < count; ++id) {
        operator[](src[id]) = id;
    }
}

bool NameToIDMap::contains(const std::string& name) const {
    return find(name)!=end();
}

bool NameToIDMap::contains(const char* name) const {
    return find(name)!=end();
}


