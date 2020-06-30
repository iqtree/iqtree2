//
//  starttree.cpp
//  alignment
//
//  Created by James Barbetti on 26/6/20.
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
