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

#include <string>            //for std::string
#include <iostream>          //for std::cout
#include <math.h>            //for log
#include <utils/progress.h>        //for progress_display::setProgressDisplay()
#include <utils/operatingsystem.h> //for getOSName
#include <utils/hammingdistance.h> //for hammingDistance
#include <utils/argument.h>        //for Argument, ArgumentMap, et al.
#include "starttree.h"       //for StartTree::Factory
#include "flatmatrix.h"      //for FlatMatrix
#include "distancematrix.h"  //for loadDistanceMatrixInto
#include "sequence.h"        //for 
#if USE_GZSTREAM
#include <utils/gzstream.h>
#endif

#define PROBLEM(x) if (1) problems << x << ".\n"; else 0

namespace {
    bool endsWith(const std::string& s, const char* suffix) {
        auto suffixLen = strlen(suffix);
        if (s.length() < suffixLen) {
            return false;
        }
        return s.substr(s.length()-suffixLen, suffixLen) == suffix;
    }

    std::string string_to_lower(const std::string& input_string) {
        std::string answer = input_string;
        std::transform(answer.begin(), answer.end(), answer.begin(),
                    []( char c){ return std::tolower(c); });
        return answer;
    }

    bool correcting_distances      = true;
    bool is_DNA                    = true;
    bool numbered_names            = false;
    bool filter_problem_sequences  = false;
    char unknown_char              = 'N';
    int  precision                 = 8;
    int  compression_level         = 9;
    std::string msaOutputPath; //write .msa formatted version of .fasta input here
    std::string alphabet;      //defaults to ACGT
    std::string unknown_chars; //defaults to .~_-?N
    std::string format             = "square.interleaved";
    bool        interleaved_format = true;

    std::string stripName;              //characters to strip from names
    std::string nameReplace("_");       //characters to replace stripepd chars with, in names
    std::string truncateName;           //truncate names when you see one of these characters
                                        //e.g. to make IQTree happy, truncate at space " ".
};

void showBanner() {
    std::cout << "\nDecentTree for " << getOSName() << "\n";
    std::cout << "Based on algorithms (UPGMA, NJ, BIONJ) proposed by Sokal & Michener [1958],\n";
    std::cout << "Saitou & Nei [1987], Gascuel [1997] and [2009]\n";
    std::cout << "Incorporating (in NJ-R and BIONJ-R) techniques proposed by Simonson, Mailund, and Pedersen [2011]\n";
    std::cout << "Developed by Olivier Gascuel [2009], Hoa Sien Cuong [2009], James Barbetti [2020]\n";
    std::cout << "(To suppress this banner pass -no-banner)\n";
}

void showUsage() {
    std::cout << "Usage: decenttree (-fasta [fastapath] (-strip-name [stripped]) \n";
    std::cout << "       (-name-replace [reps]) (-truncate-name-at [chars])\n";
    std::cout << "       (-uncorrected) (-no-matrix) (-dist-out [distout]\n";
    std::cout << "       (-alphabet [states]) (-unknown [chars]) (-not-dna))\n";
    std::cout << "       -in [mldist] (-c [level]) (-f [prec]) -out [newick] -t [algorithm]\n";
    std::cout << "       (-nt [threads]) (-gz) (-no-banner) (-q)\n";
    std::cout << "Arguments in parentheses () are optional.\n";
    std::cout << "[fastapath]  is the path of a .fasta format file specifying genetic sequences\n";
    std::cout << "             (which may be in .gz format)\n";
    std::cout << "             (by default, the character indicating an unknown state is 'N')\n";
    std::cout << "[stripped]   is a list of characters to replace in taxon names, e.g. \" /\"\n";
    std::cout << "[rep]        is a list of characters to replace them with e.g. \"_\"\n";
    std::cout << "             (may be shorter than [strippped]; if so first character is the default.\n";
    std::cout << "[distout]    is the path, of a file, into which the distance matrix is to be written\n";
    std::cout << "             (possibly in a .gz format)\n";
    std::cout << "[states]     are the characters for each site\n";
    std::cout << "[chars]      are the characters that indicate a site has an unknown character.\n";
    std::cout << "[mldist]     is the path of a distance matrix file (which may be in .gz format)\n";
    std::cout << "[newick]     is the path to write the newick tree file to (if it ends in .gz it will be compressed)\n";
    std::cout << "[threads]    is the number of threads, which should be between 1 and the number of CPUs.\n";
    std::cout << "-q           asks for quiet (less progress reporting).\n";
    std::cout << "-uncorrected turns off Jukes-Cantor distance correction (only relevant if -fasta supplied).\n";
    std::cout << "-not-dna     indicates number of states to be determined from input (if -fasta supplied).\n";
    std::cout << "-num         indicates sequence names will be replaced with A1, A2, ... in outputs.\n";
    std::cout << "[level]      is a compression level between 1 and 9 (default 5)\n";
    std::cout << "[prec]       is a precision (default 6)\n";
    std::cout << "[algorithm]  is one of the following, supported, distance matrix algorithms:\n";

    std::cout << StartTree::Factory::getInstance().getListOfTreeBuilders();
}

bool loadSequenceDistancesIntoMatrix(Sequences& sequences,
                                     const std::vector<char>&   is_site_variant,
                                     bool report_progress, FlatMatrix& m) {
    SequenceLoader loader(unknown_char, is_DNA, sequences, correcting_distances,
                          precision, compression_level, format, 
                          is_site_variant, report_progress);
    bool success = loader.loadSequenceDistances(m);
    return success;
}

bool writeMSAOutputFile(Sequences& sequences, const std::string& msaPath)
{
    //This is writing a format that fastme likes for its inputs.
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    bool success = false;
    try {
        out.open(msaPath.c_str());
        size_t len = sequences.size()==0 ? 0 : sequences[0].sequenceLength();
        out << sequences.size() << " " << len << "\n";
        
        //Todo: The following is for the interleaved sequence format
        //      that fastme likes.
        if (interleaved_format) {
            for (size_t pos=0; pos<len; pos+=60) {
                for (size_t i=0; i<sequences.size(); ++i) {
                    const std::string& seq_data = sequences[i].sequenceData();
                    if (pos==0) {
                        out << sequences.getFormattedName(i);
                    } else {
                        out << "      ";
                    }
                    size_t posEnd = pos + 60;
                    if (len<posEnd) posEnd = len;
                    for (size_t block=pos; block<posEnd; block+=10) {
                        size_t blockEnd = block + 10;
                        if (posEnd<blockEnd) {
                            posEnd = blockEnd;
                        }
                        out << " " << seq_data.substr(block, blockEnd-block);
                    }
                    out << "\n";
                }
                out << "\n";
            }
        } else {
            //Flat
            for (size_t i=0; i<sequences.size(); ++i) {
                out << sequences.getFormattedName(i) << " "
                    << sequences.at(i).sequenceData() << "\n";
            }
        }
        success = true;
    }
    catch (std::ios::failure& f) {
        std::cerr << "I/O error trying to write"
            << " MSA format file: " << msaPath << std::endl;
    }
    catch (...) {
        std::cerr << "Unexpected error trying to write"
            << " MSA format file: " << msaPath << std::endl;
    }
    out.close();
    return success;
}

void removeProblematicSequences(Sequences& sequences,
                                FlatMatrix& m) {
    size_t count_problem_sequences = sequences.countOfProblematicSequences();
    if (count_problem_sequences==0) {
        return;
    }
    Sequences  old_sequences(numbered_names);
    FlatMatrix old_matrix;
    
    std::cout << "Removing " << count_problem_sequences
        << " problematic sequences." << std::endl;
    
    double startTime = getRealTime();
    std::swap(sequences, old_sequences);
    std::swap(old_matrix, m);
    m.setSize(old_sequences.size() - count_problem_sequences);
    size_t rNew = 0;
    for (size_t r=0; r<old_sequences.size(); ++r) {
        if (!old_sequences[r].isProblematic()) {
            sequences.emplace_back(old_sequences[r].getName());
            m.sequenceName(rNew) = old_sequences[r].getName();
            size_t cNew = 0;
            for (size_t c=0; c<old_sequences.size(); ++c) {
                if (!old_sequences[c].isProblematic()) {
                    m.cell(rNew, cNew) = old_matrix.cell(r, c);
                    ++cNew;
                }
            }
            ++rNew;
        }
    }
    std::cout << "Removed sequences"
        << " in " << (getRealTime() - startTime)
        << " wall-clock seconds." << std::endl;
}

bool prepInput(const std::string& fastaFilePath,
               const std::string& phylipFilePath,
               const std::string& matrixInputFilePath,
               bool  reportProgress,
               const std::string& distanceOutputFilePath,
               Sequences& sequences, bool loadMatrix,
               FlatMatrix& m) {
    if (!fastaFilePath.empty() || !phylipFilePath.empty() ) {
        std::vector<char> is_site_variant;
        if (!sequences.loadAlignment(fastaFilePath, phylipFilePath,
                                     alphabet, unknown_char,
                                     reportProgress, is_site_variant)) {
            return false;
        }
        fixUpSequenceNames(truncateName, stripName, nameReplace, sequences);
        if (loadMatrix) {
            if (!loadSequenceDistancesIntoMatrix
                            (sequences, is_site_variant,
                            reportProgress, m)) {
                return false;
            }
        }
        if (filter_problem_sequences) {
            removeProblematicSequences(sequences, m);
        }
        if (!msaOutputPath.empty()) {
            writeMSAOutputFile(sequences, msaOutputPath);
        }
        if (!distanceOutputFilePath.empty()) {
            if (loadMatrix) {
                useNumberedNamesIfAskedTo(numbered_names, m);
                return m.writeToDistanceFile(format, precision,
                                 compression_level,
                                 distanceOutputFilePath );
            }
            else {
                SequenceLoader loader(unknown_char, is_DNA, sequences, 
                                      correcting_distances, precision, 
                                      compression_level, format,
                                      is_site_variant, reportProgress);
                bool success = loader.writeDistanceMatrixToFile(numbered_names, distanceOutputFilePath);
                return success;
            }
        }
    } else if (!matrixInputFilePath.empty()) {
        if (loadDistanceMatrixInto(matrixInputFilePath,
                                    reportProgress, m)) {
            fixUpSequenceNames(truncateName, stripName, nameReplace, m);
            if (!distanceOutputFilePath.empty()) {
                return m.writeToDistanceFile(format, precision,
                                             compression_level,
                                             distanceOutputFilePath );
            }
            return true;
        } 
        else {
            return false;
        }
    }
    else {
        return false;
    }
    return true;
}

template <class T> void range_restrict(const T lo, const T hi, T& restrict_me) {
    if (restrict_me < lo) {
        restrict_me = lo;
    } else if (hi < restrict_me) {
        restrict_me = hi;
    }
}

class DecentTreeOptions {
public:
    std::stringstream problems;
    std::string algorithmName;
    std::string fastaFilePath;          //fasta alignment
    std::string phylipFilePath;         //phylip alignment 
    std::string inputFilePath;          //phylip distance matrix formats are supported
    std::string outputFilePath;         //newick tree format
    std::string distanceOutputFilePath; //phylip distance matrix format

    bool isOutputZipped            = false;
    bool isOutputSuppressed        = false;
    bool isOutputToStandardOutput  = false; //caller asked for newick tree to go to std::cout
    bool isBannerSuppressed        = false;
    int  threadCount               = 0;
    bool beSilent                  = false;
    bool isMatrixToBeLoaded        = true;  //set to false if caller passes -no-matrix
    bool isTreeConstructionSkipped = false;
    bool beVerbose                 = false;

    ArgumentMap arg_map;

    DecentTreeOptions()
        : algorithmName(StartTree::Factory::getNameOfDefaultTreeBuilder()) {
        isOutputZipped            = false;
        isOutputSuppressed        = false;
        isOutputToStandardOutput  = false; //caller asked for newick tree to go to std::cout
        isBannerSuppressed        = false;
        threadCount               = 0;
        beSilent                  = false;
        isMatrixToBeLoaded        = true;  //set to false if caller passes -no-matrix
        isTreeConstructionSkipped = false;
        initializeArgumentMap();
    }

    void initializeArgumentMap() {
        arg_map << new StringArgument("-fasta", "fasta file path",             fastaFilePath);
        arg_map << new StringArgument("-phylip", "phylip alignment file path", phylipFilePath);
        arg_map << new StringArgument("-in",    "distance matrix file path",   inputFilePath);
        arg_map << new StringArgument("-dist",  "distance matrix file path",   inputFilePath);
        arg_map << new IntArgument   ("-c",     "compression level between 1 and 9", 
                                    compression_level);
        arg_map << new IntArgument   ("-f",     "precision level between 4 and 15",
                                    precision);
        arg_map << new StringArgument("-out-format", "output format" 
                                    "(e.g. square, upper, or lower)", format);
        arg_map << new StringArgument("-msa-out",  "msa format file path", msaOutputPath);
        arg_map << new StringArgument("-dist-out", "distance matrix file path", distanceOutputFilePath);
        arg_map << new SwitchArgument("-no-matrix", isMatrixToBeLoaded, false);
        arg_map << new StringArgument("-strip-name", "list of characters to strip from name",
                                      stripName);
        arg_map << new StringArgument("-truncate-name-at", "list of truncation characters",
                                      truncateName);
        arg_map << new StringArgument("-name-replace", "list of characters to replace"
                                    " those stripped from names", nameReplace);
        arg_map << new StringArgument("-out", "output file path", outputFilePath);
        arg_map << new SwitchArgument("-no-out",      isOutputSuppressed,       true);
        arg_map << new SwitchArgument("-std-out",     isOutputToStandardOutput, true);
        arg_map << new SwitchArgument("-gz",          isOutputZipped,           true);
        arg_map << new SwitchArgument("-no-banner",   isBannerSuppressed,       true);
        arg_map << new SwitchArgument("-uncorrected", correcting_distances,     false);
        arg_map << new IntArgument   ("-nt", "thread count", threadCount);
        arg_map << new SwitchArgument("-q",           beSilent,                 true);
        arg_map << new SwitchArgument("-filter",      filter_problem_sequences, true);
        arg_map << new StringArgument("-alphabet", "list of characters", alphabet);
        arg_map << new StringArgument("-unknown",  "list of characters", unknown_chars);
        arg_map << new SwitchArgument("-num",         numbered_names,           true);
        arg_map << new SwitchArgument("-not-dna",     is_DNA,                   false);
        arg_map << new SwitchArgument("-v",           beVerbose,                true);
    }

    void processCommandLineOptions(int argc, char* argv[]) {
        for (int argNum=1; argNum<argc; ++argNum) {
            std::string arg     = argv[argNum];
            std::string nextArg = (argNum+1<argc) ? argv[argNum+1] : "";
            
            Argument* argument = arg_map.findByName(arg);
            if (argument!=nullptr) {
                argument->accept(arg, nextArg, argv, argc, argNum, problems);
            }
            else if (arg=="-t") {
                if (START_TREE_RECOGNIZED(nextArg)) {
                    algorithmName = nextArg;
                } else if (string_to_lower(nextArg)=="none") {
                    isTreeConstructionSkipped = true;
                } else {
                    PROBLEM("Algorithm name " + nextArg + " not recognized");
                    PROBLEM("Recognized distance matrix algorithms are:");
                    PROBLEM(StartTree::Factory::getInstance().getListOfTreeBuilders());
                }
                ++argNum;
            }
            else {
                PROBLEM("Unrecognized command-line argument, " + arg);
                break;
            }
        }
        range_restrict(0, 9,  compression_level );
        range_restrict(1, 15, precision );
        format = string_to_lower(format);
        if (isOutputToStandardOutput) {
            outputFilePath = "STDOUT";
        }
        isBannerSuppressed |= beSilent;
        if (alphabet.empty() && is_DNA) {
            alphabet = "ACGT";
        }
        if (unknown_chars.empty()) {
            unknown_chars = ".~_-?N";
            unknown_char  = 'N';
        } else {
            unknown_char = unknown_chars[unknown_chars.length()-1];
        }
    }

    bool checkCommandLineOptions() {
        if (inputFilePath.empty() && fastaFilePath.empty() && phylipFilePath.empty()) {
            PROBLEM("Input (mldist) file should be specified via -in [filepath.mldist]");
            PROBLEM("or alignment (fasta) file may be specified via -fasta [filepath.fasta]");
            PROBLEM("or alignment (phylip) file, via -phylip [filpeath.phy]");
        }
        if (!fastaFilePath.empty() && !phylipFilePath.empty()) {
            PROBLEM("Cannot specify both a fast file path (with -fasta)");
            PROBLEM("and a phylip file path (with -phylip)");
        }
        if (outputFilePath.empty() && !isOutputSuppressed && !isOutputToStandardOutput) {
            PROBLEM("Ouptut (newick format) filepath should be specified via -out [filepath.newick]");
            PROBLEM("or output can be sent to standard output, or suppressed, via -std-out or -no-out");
        }
        else if (!inputFilePath.empty() && inputFilePath==outputFilePath) {
            PROBLEM("Input file and output file paths are the same (" + inputFilePath + ")");
        }
        if (fastaFilePath.empty() && phylipFilePath.empty() && !isMatrixToBeLoaded) {
            PROBLEM("If distance matrix is not be loaded, an alignment file must be specified");
        }
        if (!problems.str().empty()) {
            std::cerr << problems.str();
            return false;
        }
        return true;
    }
    void configureThreading() {
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
    }
};

int obeyCommandLineOptions(DecentTreeOptions& options);

int main(int argc, char* argv[]) {
    #if USE_PROGRESS_DISPLAY
    progress_display::setProgressDisplay(true); //Displaying progress bars
    #endif
    DecentTreeOptions options;

    options.processCommandLineOptions(argc, argv);

    if (argc==1) {
        showBanner();
        showUsage();
        return 0;
    }
    if (!options.checkCommandLineOptions()) {
        return 1;
    }
    if (!options.isBannerSuppressed) {
        showBanner();
    }
    options.configureThreading();
    return obeyCommandLineOptions(options);
}

int obeyCommandLineOptions(DecentTreeOptions& options) {
    StartTree::BuilderInterface* algorithm = 
        StartTree::Factory::getTreeBuilderByName(options.algorithmName);
    if (!options.isTreeConstructionSkipped) {
        if (algorithm==nullptr) {
            std::cerr << "Tree builder algorithm was unexpectedly null"
                      << " (internal logic error)." << std::endl;
            return 1;
        }
        algorithm->setZippedOutput(options.isOutputZipped || 
                                   endsWith(options.outputFilePath,".gz"));
        if (options.beSilent) {
            algorithm->beSilent();
        }
        algorithm->setPrecision(precision);
    }
    Sequences  sequences(numbered_names);
    FlatMatrix m;
    bool succeeded = prepInput(options.fastaFilePath, options.phylipFilePath,
                               options.inputFilePath, !algorithm->isBenchmark(),
                               options.distanceOutputFilePath,
                               sequences, options.isMatrixToBeLoaded, m);
    if (options.isTreeConstructionSkipped) {
        succeeded = true;
    }
    else if (succeeded && options.isMatrixToBeLoaded) {
        succeeded = algorithm->constructTreeInMemory(m.getSequenceNames(),
                                                     m.getDistanceMatrix(),
                                                     options.outputFilePath);    
    }
    else if (!options.inputFilePath.empty()) {
        succeeded = algorithm->constructTree(options.inputFilePath, 
                                             options.outputFilePath);
    }
    else if (!options.isMatrixToBeLoaded && 
             !options.fastaFilePath.empty() && 
             !options.phylipFilePath.empty() &&
             !options.outputFilePath.empty() ) {
        succeeded = algorithm->constructTree(options.distanceOutputFilePath, 
                                             options.outputFilePath);
    }
    else {
        std::cerr << "Distance matrix calculation failed." << std::endl;
        return 1;
    }
    if (!succeeded) {
        std::cerr << "Tree construction failed." << std::endl;
        return 1;
    }
    return 0;
}
