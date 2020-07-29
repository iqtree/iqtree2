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

int main(int argc, char* argv[]) {
    std::string problems;
    progress_display::setProgressDisplay(true); //Displaying progress bars
    std::string algorithmName  = StartTree::Factory::getNameOfDefaultTreeBuilder();
    std::string inputFilePath;
    std::string outputFilePath;
    bool isOutputZipped = false;
    for (int argNum=1; argNum<argc; ++argNum) {
        std::string arg = argv[argNum];
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
            }
            ++argNum;
        }
        else if (arg=="-out") {
            outputFilePath = nextArg;
            ++argNum;
        }
        else if (arg=="-gz") {
            isOutputZipped = true;
        }
        else {
            PROBLEM("Unrecognized command-line argument, " + arg);
            break;
        }
    }
    if (inputFilePath.empty()) {
        PROBLEM("Input (mldist) file should be specified via -in [filepath.mldist]");
    }
    if (outputFilePath.empty()) {
        PROBLEM("Ouptut (newick format) filepath should be specified via -out [filepath.newick]");
    }
    else if (inputFilePath==outputFilePath) {
        PROBLEM("Input file and output file paths are the same (" + inputFilePath + ")");
    }
    if (!problems.empty()) {
        std::cerr << problems;
        return 1;
    }
    StartTree::BuilderInterface* algorithm = StartTree::Factory::getTreeBuilderByName(algorithmName);
    algorithm->setZippedOutput(isOutputZipped || endsWith(outputFilePath,".gz"));
    if (!algorithm->constructTree(inputFilePath, outputFilePath)) {
        std::cerr << "Tree construction failed.";
        return 1;
    }
    return 0;
}
