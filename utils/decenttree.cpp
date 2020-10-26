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
#include "progress.h"        //for progress_display::setProgressDisplay()
#include "starttree.h"       //for StartTree::Factory
#include "operatingsystem.h" //for getOSName
#include "distancematrix.h"  //for loadDistanceMatrixInto
#include "hammingdistance.h" //for hammingDistance
#include "gzstream.h"

#define PROBLEM(x) if (1) problems << x << ".\n"; else 0

const char unknown_char = 'N';

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
    std::cout << "\nUsage: DecentTree (-fasta [fastapath]) -in [mldist] -out [newick] -t [algorithm] (-nt [threadcount]) (-gz) (-no-banner) (-q)\n";
    std::cout << "Arguments in parentheses () are optional.\n";
    std::cout << "[fastapath] is the path of a .fasta format file specifying genetic sequences (which may be in .gz format)\n";
    std::cout << "[mldist] is the path of a distance matrix file (which may be in .gz format)\n";
    std::cout << "[newick] is the path to write the newick tree file to (if it ends in .gz it will be compressed)\n";
    std::cout << "[threadcount] is the number of threads, which should be between 1 and the number of CPUs.\n";
    std::cout << "-q asks for quiet (less progress reporting).\n";
    std::cout << "[algorithm] is one of the following, supported, distance matrix algorithms:\n";
    std::cout << StartTree::Factory::getInstance().getListOfTreeBuilders();
}

bool processSequenceLine(std::string &sequence, std::string &line, size_t line_num) {
    //Note: this is based on processSeq from IQTree's alignment/ alignment.cpp
    //(except it returns false rather than throwing exceptions), and writes
    //errors to std::cerr.
    for (auto it = line.begin(); it != line.end(); it++) {
        if ((*it) <= ' ') continue;
        if (isalnum(*it) || (*it) == '-' || (*it) == '?'|| (*it) == '.' || (*it) == '*' || (*it) == '~')
            sequence.append(1, toupper(*it));
        else if (*it == '(' || *it == '{') {
            while (*it != ')' && *it != '}' && it != line.end()) {
                it++;
            }
            if (it == line.end()) {
                std::cerr << "Line " << line_num << ": No matching close-bracket ) or } found";
                return false;
            }
            sequence.append(1, unknown_char);
            #if (0)
                std::cerr << "NOTE: Line " << line_num
                    << ": " << line.substr(start_it-line.begin(), (it-start_it)+1)
                    << " is treated as unknown character" << std::endl;
            #endif
        } else {
            std::cerr << "Line " << line_num << ": Unrecognized character "  + std::string(1,*it);
            return false;
        }
    }
    return true;
}

bool checkLastTwoSequenceLengths(const std::vector<std::string>& sequences) {
    if (2<=sequences.size()) {
        const std::string& last        = sequences.back();
        auto               last_length = last.length();
        const std::string& penultimate = sequences[sequences.size()-2];
        if (last_length != penultimate.length()) {
            //std::cout << last << std::endl << std::endl;
            //std::cout << penultimate << std::endl << std::endl;
            std::cerr << "Sequence " << (sequences.size())
                << " had length ("          << last_length          << ")"
                << " different from that (" << penultimate.length() << ")"
                << " of the previous sequence." << std::endl;
            return false;
        }
    }
    return true;
}

double correctDistance(double obs_dist, double num_states) {
    double z = num_states / (num_states-1);
    double x = 1.0 - (z * obs_dist);
    if (x <= 0) {
        return 10; //Todo: parameter should control this
    }
    return -log(x) / z;
}

bool loadSequenceDistancesIntoMatrix(const std::vector<std::string>& seq_names,
                                     const std::vector<std::string>& sequences,
                                     const std::vector<char>&   is_site_variant,
                                     bool report_progress, FlatMatrix& m) {
    size_t rank      = seq_names.size();
    size_t rawSeqLen = sequences.front().length();
    m.setSize(rank);
    size_t seqLen    = 0; //# of characters that actually vary between two sequences
    for (auto it=is_site_variant.begin(); it!=is_site_variant.end(); ++it) {
        seqLen += *it;
    }
    for (size_t row=0; row<rank; ++row) {
        m.addCluster(seq_names[row]);
    }
    size_t     unkLen        = ((seqLen+255)/256)*4;
    char*      buffer        = new char  [ seqLen * rank];
    char**     sequence_data = new char* [ rank ];
    uint64_t*  unk_buffer    = new uint64_t  [ unkLen * rank];
    uint64_t** unknown_data  = new uint64_t* [ rank ];
    memset(unk_buffer, 0, unkLen * rank);
    {
        const char* task = report_progress ? "Extracting variant sites": "";
        progress_display extract_progress(rank, task, "extracted from", "sequence");
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t row=0; row<rank; ++row) {
            sequence_data[row] = buffer + seqLen * row;
            const char* site = sequences[row].data();
            char* write = sequence_data[row] ;
            if (seqLen<rawSeqLen) {
                for (size_t col=0; col<rawSeqLen; ++col) {
                    *write = site[col];
                    write += is_site_variant[col];
                }
            } else {
                memcpy(write, site, rawSeqLen);
            }
            
            //calculate bit array that indicates (with 1s) which
            //sites in the sequence were unknown
            unknown_data[row]     = unk_buffer + unkLen * row;
            const char* read_site = sequence_data[row];
            uint64_t*   write_unk = unknown_data[row];
            size_t           bits = 0;
            size_t            unk = 0;
            for (int col=0; col<seqLen; ++col) {
                unk <<= 1;
                if (read_site[col] == unknown_char ) {
                    ++unk;
                }
                if (++bits == 64) {
                    bits = 0;
                    *write_unk = unk;
                    ++write_unk;
                    unk = 0;
                }
            }
            if (bits!=0) {
                *write_unk = unk;
            }
            if ((row%100)==0) {
                extract_progress += 100;
            }
        }
        extract_progress.done();
    }
    {
        const char* task = report_progress ? "Calculating distances": "";
        progress_display progress( rank*(rank-1)/2, task );
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (size_t row=0; row<rank; ++row) {
            for (size_t col=row+1; col<rank; ++col) {
                uint64_t char_distance = vectorHammingDistance
                                         (unknown_char, sequence_data[row],
                                          sequence_data[col], seqLen);
                uint64_t count_unknown = countBitsSetInEither
                                         (unknown_data[row], unknown_data[col],
                                          unkLen);
                double distance  = 0;
                if (count_unknown<rawSeqLen) {
                    distance  = correctDistance(char_distance, rawSeqLen - count_unknown);
                }
                m.cell(row, col) = distance;
                m.cell(col, row) = distance;
            }
            progress += (rank-row);
        }
        progress.done();
    }
    delete [] unknown_data;
    delete [] unk_buffer;
    delete [] sequence_data;
    delete [] buffer;
    return true;
}

struct SiteInfo {
public:
    char minState;
    char maxState;
    size_t unknownCount;
    SiteInfo(): unknownCount(0) {
    }
    inline void handle(size_t sequenceIndex, char state) {
        if (state==unknown_char) {
            ++unknownCount;
        } else if (unknownCount == sequenceIndex) {
            minState = maxState = state;
        } else if (state<minState) {
            minState = state;
        } else if (maxState<state) {
            maxState = state;
        }
    }
};

bool loadAlignmentIntoDistanceMatrix(const std::string& alignmentFilePath,
                                     bool report_progress,
                                     FlatMatrix &m) {
    pigzstream in(report_progress ? "fasta" : "");
    in.open(alignmentFilePath.c_str(), std::ios::binary | std::ios::in);
    if (!in.is_open()) {
        std::cerr << "Unable to open alignment file " << alignmentFilePath << std::endl;
        return false;
    }
    size_t line_num = 0;
    std::vector<std::string> seq_names;
    std::vector<std::string> sequences;
    for (; !in.eof(); line_num++) {
        std::string line;
        safeGetLine(in, line);
        if (line == "") {
            continue;
        }
        if (line[0] == '>') { // next sequence
            auto pos = line.find_first_of("\n\r");
            seq_names.emplace_back(line.substr(1, pos-1));
            std::string& str = seq_names.back();
            str.erase(0, str.find_first_not_of(" \n\r\t"));
            str.erase(str.find_last_not_of(" \n\r\t")+1);
            if (!checkLastTwoSequenceLengths(sequences)) {
                return false;
            }
            sequences.push_back("");
            continue;
        }
        // read sequence contents
        else if (sequences.empty()) {
            std::cerr << "First line must begin with '>' to define sequence name" << std::endl;
            return false;
        }
        else if (!processSequenceLine(sequences.back(), line, line_num)) {
            return false;
        }
    }
    if (!checkLastTwoSequenceLengths(sequences)) {
        return false;
    }
    in.close();
    
    std::vector<char>   is_site_variant;
    std::vector<size_t> sequence_odd_site_count;
    {
        size_t seqLen   = sequences.front().length();
        std::vector<SiteInfo> sites;
        sites.resize(seqLen);
        SiteInfo* siteData = sites.data();
        
        size_t seqCount = sequences.size();
        for (size_t s=0; s<seqCount; ++s) {
            const char* sequence = sequences[s].data();
            for (size_t i=0; i<seqLen; ++i) {
                siteData[i].handle(s, sequence[i]);
            }
        }
        
        is_site_variant.resize(seqLen, 0);
        sequence_odd_site_count.resize(seqCount, 0);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i=0; i<seqLen; ++i) {
            SiteInfo* info = siteData + i;
            if (info->unknownCount==seqCount) {
                continue;
            }
            if (info->minState==info->maxState) {
                continue;
            }
            is_site_variant[i] = 1;
        }
    }
    
    return loadSequenceDistancesIntoMatrix(seq_names, sequences,
                                           is_site_variant,
                                           report_progress, m);
}

int main(int argc, char* argv[]) {
    std::stringstream problems;
    progress_display::setProgressDisplay(true); //Displaying progress bars
    std::string algorithmName  = StartTree::Factory::getNameOfDefaultTreeBuilder();
    std::string alignmentFilePath;
    std::string inputFilePath;
    std::string outputFilePath;
    bool isOutputZipped           = false;
    bool isOutputSuppressed       = false;
    bool isOutputToStandardOutput = false; //caller asked for newick tree to go to std::cout
    bool isBannerSuppressed       = false;
    int  threadCount              = 0;
    bool beSilent                 = false;
    for (int argNum=1; argNum<argc; ++argNum) {
        std::string arg     = argv[argNum];
        std::string nextArg = (argNum+1<argc) ? argv[argNum+1] : "";
        if (arg=="-fasta") {
            if (nextArg.empty()) {
                PROBLEM("-fasta should be followed by a file path");
            }
            alignmentFilePath = nextArg;
            ++argNum;
        }
        else if (arg=="-in" || arg=="-dist") {
            if (nextArg.empty()) {
                PROBLEM("should be followed by a file path");
            }
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
        else if (arg=="-std-out") {
            //Write output to standard output
            outputFilePath = "STDOUT";
            isOutputToStandardOutput = true;
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
    if (inputFilePath.empty() && alignmentFilePath.empty()) {
        PROBLEM("Input (mldist) file should be specified via -in [filepath.mldist]");
        PROBLEM("Or alignment (fasta) file may be specified via -fasta [filepath.fasta]");
    }
    if (outputFilePath.empty() && !isOutputSuppressed && !isOutputToStandardOutput) {
        PROBLEM("Ouptut (newick format) filepath should be specified via -out [filepath.newick]");
        PROBLEM("Or output can be sent to standard output, or suppressed, via -std-out or -no-out");
    }
    else if (!inputFilePath.empty() && inputFilePath==outputFilePath) {
        PROBLEM("Input file and output file paths are the same (" + inputFilePath + ")");
    }
    if (!problems.str().empty()) {
        std::cerr << problems.str();
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
    FlatMatrix m;
    bool succeeded = false;
    if (algorithm->isBenchmark()) {
        if (!alignmentFilePath.empty()) {
            succeeded = loadAlignmentIntoDistanceMatrix(alignmentFilePath, false, m);
        } else {
            succeeded = loadDistanceMatrixInto(inputFilePath, false, m);
        }
        if (succeeded) {
            succeeded = algorithm->constructTreeInMemory(m.getSequenceNames(),
                                                         m.getDistanceMatrix(),
                                                         outputFilePath);
        }
    }  else {
        if (!alignmentFilePath.empty()) {
            if (loadAlignmentIntoDistanceMatrix(alignmentFilePath, true, m)) {
                succeeded = algorithm->constructTreeInMemory
                    ( m.getSequenceNames(), m.getDistanceMatrix(),
                      outputFilePath );
            }
            exit(succeeded ? 0 : 1);
        }
        if (!succeeded) {
            succeeded = algorithm->constructTree(inputFilePath, outputFilePath);
        }
        if (!succeeded) {
            std::cerr << "Tree construction failed." << std::endl;
            return 1;
        }
    }
    return 0;
}
