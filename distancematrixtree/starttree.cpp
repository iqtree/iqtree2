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
#include <utils/stringfunctions.h> //for contains

namespace StartTree {

#ifndef USE_BIONJ_2009
#define USE_BIONJ_2009 (1)
#endif
#if (USE_BIONJ_2009)
    extern void addBioNJ2009TreeBuilders(Registry& f);
#endif
extern void addBioNJ2020TreeBuilders(Registry& f);

//extern void addStitchupTreeBuilders (Registry& f);

Registry::Registry() {
}

Registry::~Registry() {
    for (auto it=mapOfTreeBuilders.begin();
         it!=mapOfTreeBuilders.end(); ++it) {
        delete it->second;
    }
    mapOfTreeBuilders.clear();
}

size_t Registry::getBuilderCount() {
    return mapOfTreeBuilders.size();
}

Registry& Registry::getInstance() {
    static Registry instance;
    if (instance.getBuilderCount()==0) {
        #if (USE_BIONJ_2009)
            addBioNJ2009TreeBuilders(instance);
        #endif
        addBioNJ2020TreeBuilders(instance);
        //addStitchupTreeBuilders(instance);
        BuilderInterface *bench = new BenchmarkingTreeBuilder(instance, "BENCHMARK", "Benchmark");
        instance.addBuilder(bench->getName(), bench);
    }
    return instance;
}

void Registry::addBuilder(const std::string& name, BuilderInterface* builder) {
    auto found = mapOfTreeBuilders.find(name);
    if ( found != mapOfTreeBuilders.end()) {
        if (found->second != builder ) {
            //If replacing an existing distance matrix phylogenetic inferefence
            //implementation (owned by the registry), then the existing one
            //needs to be deleted (otherwise, any memory allocated for it
            //will be leaked).
            delete found->second;
            found->second = nullptr;
        }
    }
    mapOfTreeBuilders [ name ] = builder;
}

BuilderInterface* Registry::getBuilder(const std::string& name) const {
    auto found = mapOfTreeBuilders.find(name);
    return ( found == mapOfTreeBuilders.end())
            ? nullptr
            : found->second;
}

BuilderInterface* Registry::getBuilder(const char* name) const {
    std::string s(name);
    return getBuilder(s);
}

std::string Registry::getListOfTreeBuilders() const {
    std::stringstream list;
    for (auto it=mapOfTreeBuilders.begin();
         it!=mapOfTreeBuilders.end(); ++it) {
        auto builder = (*it).second;
        list << builder->getName() << ": "
             << builder->getDescription() << "\n";
    }
    return list.str();
}

StrVector Registry::getVectorOfTreeBuilderNames(bool withDescriptions) const {
    StrVector vector;
    for (auto it=mapOfTreeBuilders.begin();
         it!=mapOfTreeBuilders.end(); ++it) {
        auto builder     = (*it).second;
        std::string name = builder->getName();
        if (withDescriptions) {
            name += ": ";
            name += builder->getDescription();
        }
        vector.emplace_back(name);
    }
    return vector;
}

void Registry::advertiseTreeBuilder(BuilderInterface* builder) {
    std::string name = builder->getName();
    addBuilder( name, builder );
}

void Registry::setNameOfDefaultTreeBuilder(const char* name) {
    nameOfDefaultTreeBuilder = name;
}

const std::string& Registry::getNameOfDefaultTreeBuilder() {
    return getInstance().nameOfDefaultTreeBuilder;
}

BuilderInterface* Registry::getTreeBuilderByName(const std::string& name) {
    return getInstance().getBuilder(name);
}

BuilderInterface* Registry::getDefaultTreeBuilder() const {
    return getBuilder(nameOfDefaultTreeBuilder);
}

BuilderInterface* Registry::getTreeBuilderByName(const char* name) {
    if (name!=nullptr && strlen(name)!=0) {
        return getInstance().getBuilder(name);
    } 
    else {
        return getInstance().getDefaultTreeBuilder();
    }
}

BenchmarkingTreeBuilder::BenchmarkingTreeBuilder(Registry& f, const char* nameToUse,
                                                 const char *descriptionToGive)
    : name(nameToUse), description(descriptionToGive)
    , isOutputToBeZipped(false), silent(false), precision(4)
    , subtreeOnly(false) {
    for (auto it=f.mapOfTreeBuilders.begin(); it!=f.mapOfTreeBuilders.end(); ++it) {
        if (!it->second->getName().empty()) {
            builders.push_back(it->second);
        }
    }
}

const std::string& BenchmarkingTreeBuilder::getName() const {
    return name;
}

const std::string& BenchmarkingTreeBuilder::getDescription() {
    return description;
}

bool BenchmarkingTreeBuilder::isBenchmark() const {
    return true;
}

bool BenchmarkingTreeBuilder::constructTree
    ( const std::string& distanceMatrixFilePath
    , const std::string& newickTreeFilePath) {
    bool result = (!builders.empty());
    for (auto it=builders.begin(); it!=builders.end(); ++it) {
        (*it)->setIsRooted(isRooted);
        (*it)->setZippedOutput(isOutputToBeZipped);
        (*it)->setSubtreeOnly(subtreeOnly);
        result &= (*it)->constructTree(distanceMatrixFilePath,
                                       newickTreeFilePath);
    }
    return result;
}

bool BenchmarkingTreeBuilder::setZippedOutput(bool zipIt) {
    isOutputToBeZipped = zipIt;
    return true;
}

void BenchmarkingTreeBuilder::beSilent() {
    silent = true;
}

bool BenchmarkingTreeBuilder::setAppendFile(bool /*appendIt*/) {
    return false;
}

bool BenchmarkingTreeBuilder::setSubtreeOnly(bool /*appendIt*/) {
    return false;
}

bool BenchmarkingTreeBuilder::setPrecision(int precisionToUse) {
    precision = precisionToUse;
    return true;
}

bool BenchmarkingTreeBuilder::setIsRooted(bool rootIt) {
    isRooted = rootIt;
    return true;
}

namespace {
    std::string formatPositiveNumber(double n, intptr_t w) {
        std::stringstream s;
        s.precision( (w>=2) ? (w-2) : 0);
        s << n;
        std::string t = s.str();

        if (contains(t,"e")) {
            return t;
        }

        if (1.0 <= n ) {
            if ( t.length() < (size_t)w ) {
                return std::string(w-t.length(), ' ') + t;
            } else {
                return t;
            }
        }
        if (t.length() == (size_t)w || w < 2 ) {
            return t;
        }
        if (t.length() < (size_t)w  ) {
            return t + std::string(w-t.length(), '0');
        }
        bool carry = (t[w] >= '5');
        if (carry) {
            intptr_t j = w - 1;
            for (; t[j] == '9'; --j ) {
                t[j] = '0';
            }
            ++t[j];
        }
        return t.substr(0,w);
    }
}

bool BenchmarkingTreeBuilder::constructTreeStringInMemory
    ( const StrVector& sequenceNames
    , const double*    distanceMatrix
    , std::string&     output_string) {
    return false;
}

bool BenchmarkingTreeBuilder::constructTreeInMemory
( const StrVector &sequenceNames
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
        (*it)->setPrecision(precision);
        #ifdef _OPENMP
            omp_set_num_threads(1);
        #endif
        double startTime = getRealTime();
        bool   succeeded = (*it)->constructTreeInMemory(sequenceNames,
                                                        distanceMatrix,
                                                        newickTreeFilePath);
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
                ok &= (*it)->constructTreeInMemory(sequenceNames,
                                                   distanceMatrix,
                                                   newickTreeFilePath);
                elapsed = getRealTime() - startTime;
                std::cout << "\t" << formatPositiveNumber(elapsed,7);
            }
            #endif
            std::cout << std::endl;
        }
    }
    return ok;
}

bool BenchmarkingTreeBuilder::constructTreeAndAppendToStream
    ( const StrVector& sequenceNames
    , const double*    distanceMatrix
    , std::iostream&   newickTreeFilePath) {
        return false;
}

};

