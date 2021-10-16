#pragma once
#ifndef NAME_TO_ID_MAP
#define NAME_TO_ID_MAP
#include "vectortypes.h"
#include <map>

//Mapping names to ids
class NameToIDMap: public std::map<std::string, int> {
public:
    NameToIDMap()  = default;
    ~NameToIDMap() = default;
    NameToIDMap(const NameToIDMap& rhs) = default;
    NameToIDMap& operator=(const NameToIDMap& rhs) = default;
    explicit NameToIDMap(const StrVector& src);
    bool contains(const std::string& name) const;
    bool contains(const char*        name) const;
};

#endif //NAME_TO_ID_MAP