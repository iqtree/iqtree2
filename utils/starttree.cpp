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
#include <iostream> //for std::cout
#include <sstream>  //for std::stringstream

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
        BuilderInterface *bench = new BenchmarkingTreeBuilder(instance, "BENCHMARK", "Benchmark");
        instance.addBuilder(bench->getName(), bench);
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

std::string Factory::getListOfTreeBuilders() const {
    std::stringstream list;
    for (auto it=mapOfTreeBuilders.begin(); it!=mapOfTreeBuilders.end(); ++it) {
        auto builder = (*it).second;
        list << builder->getName() << ": " << builder->getDescription() << "\n";
    }
    return list.str();
}


void Factory::advertiseTreeBuilder(BuilderInterface* builder) {
    std::string name = builder->getName();
    addBuilder( name, builder );
}

void Factory::setNameOfDefaultTreeBuilder(const char* name) {
    nameOfDefaultTreeBuilder = name;
}

const std::string& Factory::getNameOfDefaultTreeBuilder() {
    return getInstance().nameOfDefaultTreeBuilder;
}

BuilderInterface* Factory::getTreeBuilderByName(const std::string& name) {
    return getInstance().getBuilder(name);
}

BenchmarkingTreeBuilder::BenchmarkingTreeBuilder(Factory& f, const char* nameToUse, const char *descriptionToGive)
    : name(nameToUse), description(descriptionToGive)
    , isOutputToBeZipped(false) {
    for (auto it=f.mapOfTreeBuilders.begin(); it!=f.mapOfTreeBuilders.end(); ++it) {
        if (!it->second->getName().empty()) {
            builders.push_back(it->second);
        }
    }
}

const std::string& BenchmarkingTreeBuilder::getName() {
    return name;
}
const std::string& BenchmarkingTreeBuilder::getDescription() {
    return description;
}

bool BenchmarkingTreeBuilder::constructTree
    ( const std::string &distanceMatrixFilePath
     , const std::string & newickTreeFilePath) {
        bool result = (!builders.empty());
        for (auto it=builders.begin(); it!=builders.end(); ++it) {
            (*it)->setZippedOutput(isOutputToBeZipped);
            result &= (*it)->constructTree(distanceMatrixFilePath, newickTreeFilePath);
        }
        return result;
    }

void BenchmarkingTreeBuilder::setZippedOutput(bool zipIt) {
    isOutputToBeZipped = zipIt;
}

bool BenchmarkingTreeBuilder::constructTreeInMemory
    ( const std::vector<std::string> &sequenceNames
    , double *distanceMatrix
    , const std::string & newickTreeFilePath) {
        bool ok = false;
        #ifdef _OPENMP
            int maxThreads = omp_get_max_threads();
        #endif
        for (auto it=builders.begin(); it!=builders.end(); ++it) {
            double startTime = getRealTime();
            #ifdef _OPENMP
                omp_set_num_threads(1);
            #endif
            (*it)->beSilent();
            bool succeeded = (*it)->constructTreeInMemory(sequenceNames, distanceMatrix, newickTreeFilePath);
            double elapsed = getRealTime() - startTime;
            if (succeeded) {
                ok = true;
                std::cout.precision(6);
                std::cout << (*it)->getName() << " \t" << elapsed;
                #ifdef _OPENMP
                for (int t=2; t<=maxThreads; ++t) {
                    omp_set_num_threads(t);
                    startTime = getRealTime();
                    ok &= (*it)->constructTreeInMemory(sequenceNames, distanceMatrix, newickTreeFilePath);
                    elapsed = getRealTime() - startTime;
                    std::cout << "\t" << (elapsed);
                }
                #endif
                std::cout << std::endl;
            }
        }
        return true;
    }
};
