//
//  DecentTree.cpp
//
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

#include <string>      //for std::string
#include <iostream>    //for std::cout
#include "progress.h"  //for progress_display::setProgressDisplay()
#include "starttree.h" //for StartTree::Factory
#include "operatingsystem.h" //for getOSName
#include "distancematrix.h"  //for loadDistanceMatrixInto

#define PROBLEM(x) if (1) problems = problems + x + ".\n"; else 0

namespace {
    bool endsWith(const std::string s, const char* suffix) {
        auto suffixLen = strlen(suffix);
        if (s.length() < suffixLen) {
            return false;
        }
        return s.substr(s.length()-suffixLen, suffixLen) == suffix;
    }
};

class FlatMatrix {
private:
    std::vector<std::string> sequenceNames;
    size_t                   rowCount;
    double*                  distanceMatrix;
public:
    FlatMatrix(): rowCount(0), distanceMatrix(nullptr) {
    }
    virtual ~FlatMatrix() {
        delete [] distanceMatrix;
        distanceMatrix = nullptr;
    }
    const std::vector<std::string> getSequenceNames() const {
        return sequenceNames;
    }
    void setSize(size_t rows) {
        delete [] distanceMatrix;
        rowCount = rows;
        distanceMatrix = new double [ rowCount * rowCount ];
    }
    size_t getSize() {
        return rowCount;
    }
    const double* getDistanceMatrix() const {
        return distanceMatrix;
    }
    double cell(size_t r, size_t c) const {
        return distanceMatrix[r * rowCount + c];
    }
    double& cell(size_t r, size_t c) {
        return distanceMatrix[r * rowCount + c];
    }
    void addCluster(const std::string& clusterName) {
        sequenceNames.emplace_back(clusterName);
    }
};

void showBanner() {
    std::cout << "\nDecentTree for " << getOSName() << "\n";
    std::cout << "Based on algorithms (UPGMA, NJ, BIONJ) proposed by Sokal & Michener [1958], Saitou & Nei [1987], Gascuel [2009]\n";
    std::cout << "Incorporating (in NJ-R and BIONJ-R) techniques proposed by Simonson, Mailund, and Pedersen [2011]\n";
    std::cout << "Developed by Olivier Gascuel [2009], Hoa Sien Cuong [2009], James Barbetti [2020]\n";
    std::cout << "(To suppress this banner pass -no-banner)\n";
}

void showUsage() {
    std::cout << "\nUsage: DecentTree -in [mldist] -out [newick] -t [algorithm] -nt [threadcount] (-gz) (-no-banner)\n";
    std::cout << "[mldist] is the path of a distance matrix file (which may be in .gz format)\n";
    std::cout << "[newick] is the path to write the newick tree file to (if it ends in .gz it will be compressed)\n";
    std::cout << "[algorithm] is one of the following, supported, distance matrix algorithms:\n";
    std::cout << "[threadcount] is the number of threads, which should be between 1 and the number of CPUs.\n";
    std::cout << StartTree::Factory::getInstance().getListOfTreeBuilders();
}

int main(int argc, char* argv[]) {
    std::string problems;
    progress_display::setProgressDisplay(true); //Displaying progress bars
    std::string algorithmName  = StartTree::Factory::getNameOfDefaultTreeBuilder();
    std::string inputFilePath;
    std::string outputFilePath;
    bool isOutputZipped     = false;
    bool isOutputSuppressed = false;
    bool isBannerSuppressed = false;
    int  threadCount        = 0;
    bool beSilent           = false;
    for (int argNum=1; argNum<argc; ++argNum) {
        std::string arg     = argv[argNum];
        std::string nextArg = (argNum+1<argc) ? argv[argNum+1] : "";
        if (arg=="-in") {
            inputFilePath = nextArg;
            ++argNum;
        }
        else if (arg=="-t") {
            if (START_TREE_RECOGNIZED(nextArg)) {
                algorithmName = nextArg;
            } else {
                PROBLEM("Algorithm name " + nextArg + " not recognized");
                PROBLEM("Recognized distance matrix algorithms are:");
                PROBLEM(StartTree::Factory::getInstance().getListOfTreeBuilders());
            }
            ++argNum;
        }
        else if (arg=="-out") {
            outputFilePath = nextArg;
            ++argNum;
        }
        else if (arg=="-no-out") {
            isOutputSuppressed = true;
        }
        else if (arg=="-gz") {
            isOutputZipped = true;
        }
        else if (arg=="-no-banner") {
            isBannerSuppressed = true;
        }
        else if (arg=="-nt") {
            if ( nextArg.empty() || nextArg[0]<'0' || '9'<nextArg[0] ) {
                PROBLEM("-nt argument should be followed by a numeric thread count");
                break;
            }
            threadCount = atol(nextArg.c_str());
            ++argNum;
        }
        else if (arg=="-q") {
            isBannerSuppressed = true;
            beSilent = true;
        }
        else {
            PROBLEM("Unrecognized command-line argument, " + arg);
            break;
        }
    }
    if (argc==1) {
        if (!isBannerSuppressed) {
            showBanner();
        }
        showUsage();
        return 0;
    }
    if (inputFilePath.empty()) {
        PROBLEM("Input (mldist) file should be specified via -in [filepath.mldist]");
    }
    if (outputFilePath.empty() && !isOutputSuppressed) {
        PROBLEM("Ouptut (newick format) filepath should be specified via -out [filepath.newick]");
    }
    else if (inputFilePath==outputFilePath) {
        PROBLEM("Input file and output file paths are the same (" + inputFilePath + ")");
    }
    if (!problems.empty()) {
        std::cerr << problems;
        return 1;
    }
    if (!isBannerSuppressed) {
        showBanner();
    }
    if (threadCount!=0) {
#ifdef _OPENMP
        int maxThreadCount = omp_get_max_threads();
        if (maxThreadCount < threadCount ) {
            std::cerr << "Warning: Requested number of threads, " << threadCount
                << " is greater than the maximum, " << maxThreadCount << "." << std::endl;
            std::cerr << "Warning: " << maxThreadCount << " threads will be used." << std::endl;
            threadCount = maxThreadCount;
        }
        omp_set_num_threads(threadCount);
#else
        std::cerr << "Warning: -nt argument, requesting " << threadCount
            << " thread can not be honoured (Open MP is not enabled in this build)." << std::endl;
        std::cerr << "Warning: Distance matrix processing will be single-threaded." << std::endl;
#endif
    }
    StartTree::BuilderInterface* algorithm = StartTree::Factory::getTreeBuilderByName(algorithmName);
    if (algorithm==nullptr) {
        std::cerr << "Tree builder algorithm was unexpectedly null (internal logic error)." << std::endl;
        return 1;
    }
    algorithm->setZippedOutput(isOutputZipped || endsWith(outputFilePath,".gz"));
    if (beSilent) {
        algorithm->beSilent();
    }
    if (algorithm->isBenchmark()) {
        FlatMatrix m;
        loadDistanceMatrixInto(inputFilePath, false, m);
        algorithm->constructTreeInMemory(m.getSequenceNames(),
                                         m.getDistanceMatrix(),
                                         outputFilePath);
    }  else {
        bool succeeded = algorithm->constructTree(inputFilePath, outputFilePath);
        if (!succeeded) {
            std::cerr << "Tree construction failed." << std::endl;
            return 1;
        }
    }
    return 0;
}
