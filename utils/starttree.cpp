//  Copyright (C) 2020, James Barbetti.
//
//  LICENSE:
//* This program is free software; you can redistribute it and/or modify
//* it under the terms of the GNU General Public License as published by
//* the Free Software Foundation; either version 2 of the License, or
//* (at your option) any later version.
//*
//* This program is distributed in the hope that it will be useful,
//* but WITHOUT ANY WARRANTY; without even the implied warranty of
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//* GNU General Public License for more details.
//*
//* You should have received a copy of the GNU General Public License
//* along with this program; if not, write to the
//* Free Software Foundation, Inc.,
//* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//

#include "starttree.h"
#include <iostream>

namespace StartTree {

extern void addBioNJ2009TreeBuilders(Factory& f);
extern void addBioNJ2020TreeBuilders(Factory& f);

Factory::Factory() {
}

Factory::~Factory() {
    for (auto it=mapOfTreeBuilders.begin(); it!=mapOfTreeBuilders.end(); ++it) {
        delete it->second;
    }
    mapOfTreeBuilders.clear();
}

size_t Factory::getBuilderCount() {
    return mapOfTreeBuilders.size();
}

Factory& Factory::getInstance() {
    static Factory instance;
    if (instance.getBuilderCount()==0) {
        addBioNJ2009TreeBuilders(instance);
        addBioNJ2020TreeBuilders(instance);
    }
    return instance;
}

void Factory::addBuilder(const std::string& name, BuilderInterface* builder) {
    mapOfTreeBuilders [ name ] = builder;
}

BuilderInterface* Factory::getBuilder(const std::string& name) {
    auto found = mapOfTreeBuilders.find(name);
    return ( found == mapOfTreeBuilders.end())
            ? nullptr
            : found->second;
}

BuilderInterface* Factory::getBuilder(const char* name) {
    std::string s(name);
    return getBuilder(s);
}

void Factory::advertiseTreeBuilder(BuilderInterface* builder) {
    std::string name = builder->getName();
    addBuilder( name, builder );
}

BuilderInterface* Factory::getTreeBuilderByName(const std::string& name) {
    return getInstance().getBuilder(name);
}

};
