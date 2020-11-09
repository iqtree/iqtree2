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
    , isOutputToBeZipped(false), silent(false) {
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

bool BenchmarkingTreeBuilder::isBenchmark() const {
    return true;
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

void BenchmarkingTreeBuilder::beSilent() {
    silent = true;
}

namespace {
    std::string formatPositiveNumber(double n, size_t w) {
        std::stringstream s;
        s.precision( (w>=2) ? (w-2) : 0);
        s << n;
        std::string t = s.str();
        if (1.0 <= n ) {
            if ( t.length() < w ) {
                return std::string(w-t.length(), ' ') + t;
            } else {
                return t;
            }
        }
        if (t.length() == w || w < 2 ) {
            return t;
        }
        if (t.length() < w  ) {
            return t + std::string(w-t.length(), '0');
        }
        bool carry = (t[w] >= '5');
        if (carry) {
            int j = w - 1;
            for (; t[j] == '9'; --j ) {
                t[j] = '0';
            }
            ++t[j];
        }
        return t.substr(0,w);
    }
}

bool BenchmarkingTreeBuilder::constructTreeInMemory
    ( const std::vector<std::string> &sequenceNames
    , const double *distanceMatrix
    , const std::string & newickTreeFilePath) {
        bool ok = false;
        #ifdef _OPENMP
            int maxThreads = omp_get_max_threads();
        #endif
        size_t maxNameLen = 0;
        for (auto it=builders.begin(); it!=builders.end(); ++it) {
            auto builder_name = (*it)->getName();
            if ( maxNameLen < builder_name.length() ) {
                maxNameLen = builder_name.length();
            }
        }
        std::string padding(maxNameLen, '\x20');
        std::cout.width(0);
        for (auto it=builders.begin(); it!=builders.end(); ++it) {
            auto   builder_name = (*it)->getName();
            (*it)->beSilent();
            #ifdef _OPENMP
                omp_set_num_threads(1);
            #endif
            double startTime = getRealTime();
            bool   succeeded = (*it)->constructTreeInMemory(sequenceNames, distanceMatrix, newickTreeFilePath);
            double elapsed   = getRealTime() - startTime;
            if (succeeded) {
                ok = true;
                std::cout << (builder_name + padding).substr(0, maxNameLen);
                std::cout << "\t" << formatPositiveNumber(elapsed,7);
                #ifdef _OPENMP
                for (int t=2; t<=maxThreads; ++t) {
                    std::cout.flush();
                    omp_set_num_threads(t);
                    startTime = getRealTime();
                    ok &= (*it)->constructTreeInMemory(sequenceNames, distanceMatrix, newickTreeFilePath);
                    elapsed = getRealTime() - startTime;
                    std::cout << "\t" << formatPositiveNumber(elapsed,7);
                }
                #endif
                std::cout << std::endl;
            }
        }
        return true;
    }
};

