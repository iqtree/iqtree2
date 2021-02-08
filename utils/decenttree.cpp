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
#include "flatmatrix.h"      //for FlatMatrix
#include "distancematrix.h"  //for loadDistanceMatrixInto
#include "hammingdistance.h" //for hammingDistance
#if USE_GZSTREAM
#include "gzstream.h"
#endif

#define PROBLEM(x) if (1) problems << x << ".\n"; else 0

namespace {
    bool endsWith(const std::string s, const char* suffix) {
        auto suffixLen = strlen(suffix);
        if (s.length() < suffixLen) {
            return false;
        }
        return s.substr(s.length()-suffixLen, suffixLen) == suffix;
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
    std::cout << "Usage: decenttree (-fasta [fastapath] (-uncorrected)\n";
    std::cout << "       (-alphabet [states]) (-unknown [chars]) (-not-dna))\n";
    std::cout << "       -in [mldist] (-c [level]) (-f [prec]) -out [newick] -t [algorithm]\n";
    std::cout << "       (-nt [threads]) (-gz) (-no-banner) (-q)\n";
    std::cout << "Arguments in parentheses () are optional.\n";
    std::cout << "[fastapath]  is the path of a .fasta format file specifying genetic sequences\n";
    std::cout << "             (which may be in .gz format)\n";
    std::cout << "             (by default, the character indicating an unknown state is 'N')\n";
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

bool processSequenceLine(const std::vector<int> &in_alphabet,
                         std::string &sequence,
                         std::string &line, size_t line_num) {
    //Note: this is based on processSeq from IQTree's alignment/ alignment.cpp
    //(except it returns false rather than throwing exceptions), and writes
    //errors to std::cerr.
    for (auto it = line.begin(); it != line.end(); it++) {
        if ((*it) <= ' ') continue;
        if (isalnum(*it) || (*it) == '-' || (*it) == '?'|| (*it) == '.' 
            || (*it) == '*' || (*it) == '~') {
            auto c = toupper(*it);
            if (!in_alphabet[c]) {
                c = unknown_char;
            }
            sequence.append(1, c);
        }
        else if (*it == '(' || *it == '{') {
            while (*it != ')' && *it != '}' && it != line.end()) {
                it++;
            }
            if (it == line.end()) {
                std::cerr << "Line " << line_num 
                    << ": No matching close-bracket ) or } found";
                return false;
            }
            sequence.append(1, unknown_char);
            #if (0)
                std::cerr << "NOTE: Line " << line_num
                    << ": " << line.substr(start_it-line.begin(), (it-start_it)+1)
                    << " is treated as unknown character" << std::endl;
            #endif
        } else {
            std::cerr << "Line " << line_num 
                << ": Unrecognized character "  + std::string(1,*it);
            return false;
        }
    }
    return true;
}

inline double correctDistance(double char_dist, double chars_compared,
                              double num_states) {
    double obs_dist = (0<chars_compared)
                    ? (char_dist/chars_compared) : 0.0;
    double z = num_states / (num_states-1);
    double x = 1.0 - (z * obs_dist);
    if (x <= 0) {
        return 10; //Todo: parameter should control this
    }
    return -log(x) / z;
}

inline double uncorrectedDistance(double char_dist,
                                  double chars_compared) {
    return (0<chars_compared)
        ? (char_dist/chars_compared) : 0.0;
}

struct Sequence {
protected:
    std::string name;
    std::string sequence_data;
    bool        is_problematic;

public:
    explicit Sequence(const std::string& seq_name)
        : name(seq_name), is_problematic(false) {}
    size_t sequenceLength()           const { return sequence_data.size(); }
    const char* data()                const { return sequence_data.data(); }
    const std::string& sequenceData() const { return sequence_data; }
          std::string& sequenceData()       { return sequence_data; }
    const std::string& getName()      const { return name; }
    bool  isProblematic()             const { return is_problematic; }
    void  markAsProblematic()               { is_problematic = true; }
};

struct Sequences: public std::vector<Sequence> {
    bool checkLastTwoSequenceLengths() const {
        if (2<=size()) {
            const std::string& last        = back().sequenceData();
            auto               last_length = last.length();
            const std::string& penultimate = at(size()-2).sequenceData();
            if (last_length != penultimate.length()) {
                //std::cout << last << std::endl << std::endl;
                //std::cout << penultimate << std::endl << std::endl;
                std::cerr << "Sequence " << (size())
                    << " had length ("          << last_length          << ")"
                    << " different from that (" << penultimate.length() << ")"
                    << " of the previous sequence." << std::endl;
                return false;
            }
        }
        return true;
    }
    size_t countOfProblematicSequences() {
        size_t count = 0;
        for (size_t i=0; i<size(); ++i) {
            if (at(i).isProblematic()) {
                ++count;
            }
        }
        return count;
    }
    std::string getFormattedName(size_t i) {
        if (numbered_names) {
            std::stringstream number_name;
            number_name << "A" << (i+1); //"A1", "A2", ...
            return number_name.str();
        } else {
            return at(i).getName();
        }
    }
};

bool loadSequenceDistancesIntoMatrix(Sequences& sequences,
                                     const std::vector<char>&   is_site_variant,
                                     bool report_progress, FlatMatrix& m) {
    intptr_t rank      = sequences.size();
    size_t   rawSeqLen = sequences.front().sequenceLength();
    size_t   seqLen    = 0; //# of characters that actually vary between two sequences
    m.setSize(rank);
    for (auto it=is_site_variant.begin(); it!=is_site_variant.end(); ++it) {
        seqLen += *it;
    }
    for (intptr_t row=0; row<rank; ++row) {
        m.addCluster(sequences[row].getName());
    }
    size_t     unkLen        = ((seqLen+255)/256)*4;
    char*      buffer        = new char      [ seqLen * rank ];
    char**     sequence_data = new char*     [ rank ];
    uint64_t*  unk_buffer    = new uint64_t  [ unkLen * rank ];
    uint64_t** unknown_data  = new uint64_t* [ rank ];
    memset(unk_buffer, 0, unkLen * rank);
    
    {
        #if USE_PROGRESS_DISPLAY
        const char* task = report_progress ? "Extracting variant sites": "";
        progress_display extract_progress(rank, task, "extracted from", "sequence");
        #else
        double extract_progress = 0.0;
        #endif
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t row=0; row<rank; ++row) {
            sequence_data[row] = buffer + seqLen * row;
            const char* site   = sequences[row].data();
            char*       write  = sequence_data[row] ;
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
                extract_progress += 100.0;
            }
        }
        #if USE_PROGRESS_DISPLAY
        extract_progress.done();
        #endif
    }
    
    //Determine the number of states (needed for correcting distances)
    double num_states = 0.0;
    if (is_DNA) {
        num_states = 4;
    }
    else
    {
        std::vector<size_t> char_counts;
        char_counts.resize(256, 0);
        auto char_count_array = char_counts.data();
        const unsigned char* start_buffer = reinterpret_cast<unsigned char*>(buffer);
        const unsigned char* end_buffer   = start_buffer + rank * seqLen;
        for (const unsigned char* scan=start_buffer; scan<end_buffer; ++scan) {
            ++char_count_array[*scan];
        }
        for (int i=0; i<256; ++i) {
            num_states += (char_counts[i]==0) ? 0 : 1;
        }
        if (0<char_counts[unknown_char]) {
            --num_states;
        }
    }
    
    {
        #if USE_PROGRESS_DISPLAY
        const char* task = report_progress ? "Calculating distances": "";
        progress_display progress( rank*(rank-1)/2, task );
        #else
        double progress = 0.0;
        #endif
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t row=0; row<rank; ++row) {
            for (intptr_t col=row+1; col<rank; ++col) {
                uint64_t char_distance = vectorHammingDistance
                                         (unknown_char, sequence_data[row],
                                          sequence_data[col], seqLen);
                uint64_t count_unknown = countBitsSetInEither
                                         (unknown_data[row], unknown_data[col],
                                          unkLen);
                double   distance      = 0;
                intptr_t adjSeqLen     = rawSeqLen - count_unknown;
                if (0<adjSeqLen) {
                    if (correcting_distances) {
                        distance = correctDistance((double)char_distance, (double)adjSeqLen, (double)num_states);
                    } else {
                        distance = uncorrectedDistance((double)char_distance, (double)adjSeqLen);
                    }
                    if (distance < 0) {
                        distance = 0;
                    }
                } else {
                    bool eitherMarked = sequences[row].isProblematic()
                                     || sequences[col].isProblematic();
                    if (!eitherMarked)
                    {
                        //Cannot calculate distance between these two sequences.
                        //Mark one of them (the one with more unknowns) as problematic.
                        uint64_t unkRow = countBitsSetIn(unknown_data[row], unkLen);
                        uint64_t unkCol = countBitsSetIn(unknown_data[col], unkLen);
                        intptr_t zap    = (unkCol < unkRow) ? row : col;
                        sequences[zap].markAsProblematic();
                    }
                }
                m.cell(row, col) = distance;
                m.cell(col, row) = distance;
            }
            progress += (rank-row);
        }
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
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

bool loadAlignment(const std::string& alignmentFilePath,
                   bool report_progress, Sequences &sequences,
                   std::vector<char>& is_site_variant)
{
    #if USE_GZSTREAM
    pigzstream    in(report_progress ? "fasta" : "");
    #else
    std::ifstream in;
    #endif
    in.open(alignmentFilePath.c_str(), std::ios::binary | std::ios::in);
    if (!in.is_open()) {
        std::cerr << "Unable to open alignment file " 
            << alignmentFilePath << std::endl;
        return false;
    }
    size_t line_num = 0;
    std::vector<int> in_alphabet;
    in_alphabet.resize(256, 0);
    for (auto alpha=alphabet.begin(); alpha!=alphabet.end(); ++alpha) {
        in_alphabet[*alpha] = 1;
    }
    for (; !in.eof(); line_num++) {
        std::string line;
        safeGetLine(in, line);
        if (line == "") {
            continue;
        }
        if (line[0] == '>') { // next sequence
            auto pos = line.find_first_of("\n\r");
            std::string str = line.substr(1, pos-1);
            str.erase(0, str.find_first_not_of(" \n\r\t"));
            str.erase(str.find_last_not_of(" \n\r\t")+1);
            if (!sequences.checkLastTwoSequenceLengths()) {
                return false;
            }
            sequences.emplace_back(str);
            continue;
        }
        // read sequence contents
        else if (sequences.empty()) {
            std::cerr << "First line must begin with '>'"
                << " to define sequence name" << std::endl;
            return false;
        }
        else if (!processSequenceLine(in_alphabet, sequences.back().sequenceData(),
                                      line, line_num)) {
            return false;
        }
    }
    if (!sequences.checkLastTwoSequenceLengths()) {
        return false;
    }
    in.close();
    
    std::vector<size_t> sequence_odd_site_count;
    {
        size_t seqLen   = sequences.front().sequenceLength();
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
    return true;
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
    catch (std::ios::failure f) {
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

bool writeDistanceMatrixToFile(FlatMatrix& m,
                               const std::string& distFilePath) {
    if (numbered_names) {
        auto name_count = m.getSequenceNames().size();
        for (size_t i=0; i<name_count; ++i) {
            std::stringstream name;
            name << "A" << (i+1); //"A1", "A2", ...
            m.sequenceName(i) = name.str();
        }
    }
    return m.writeToDistanceFile(format, precision,
                                 compression_level,
                                 distFilePath );
}

void removeProblematicSequences(Sequences& sequences,
                                FlatMatrix& m) {
    size_t count_problem_sequences = sequences.countOfProblematicSequences();
    if (count_problem_sequences==0) {
        return;
    }
    Sequences  old_sequences;
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

bool prepInput(const std::string& alignmentInputFilePath,
               const std::string& matrixInputFilePath,
               bool  reportProgress,
               const std::string& distanceOutputFilePath,
               Sequences& sequences, FlatMatrix& m) {
    if (!alignmentInputFilePath.empty()) {
        Sequences sequences;
        std::vector<char> is_site_variant;
        if (!loadAlignment(alignmentInputFilePath, reportProgress,
            sequences, is_site_variant)) {
            return false;
        }
        if (!loadSequenceDistancesIntoMatrix
                         (sequences, is_site_variant,
                          reportProgress, m)) {
            return false;
        }
        if (filter_problem_sequences) {
            removeProblematicSequences(sequences, m);
        }
        if (!msaOutputPath.empty()) {
            writeMSAOutputFile(sequences, msaOutputPath);
        }
        if (!distanceOutputFilePath.empty()) {
            writeDistanceMatrixToFile(m, distanceOutputFilePath);
        }
    } else if (!matrixInputFilePath.empty()) {
        if (!loadDistanceMatrixInto(matrixInputFilePath,
                                    reportProgress, m)) {
            return false;
        }
    }
    else {
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::stringstream problems;
    #if USE_PROGRESS_DISPLAY
    progress_display::setProgressDisplay(true); //Displaying progress bars
    #endif
    std::string algorithmName  = StartTree::Factory::getNameOfDefaultTreeBuilder();
    std::string alignmentFilePath; //only .fasta format is supported
    std::string inputFilePath;     //phylip distance matrix formats are supported
    std::string outputFilePath;    //newick tree format
    std::string distanceOutputFilePath; //phylip distance matrix format
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
                PROBLEM(arg + " should be followed by a file path");
            }
            inputFilePath = nextArg;
            ++argNum;
        }
        else if (arg=="-c") {
            if (nextArg.empty()) {
                PROBLEM(arg + " should be followed by compression level between 1 and 9");
            }
            compression_level = atoi(arg.c_str());
            if (compression_level<0) compression_level = 0;
            if (9<compression_level) compression_level = 9;
            ++argNum;
        }
        else if (arg=="-f") {
            if (nextArg.empty()) {
                PROBLEM(arg + " should be followed by precision level between 4 and 15");
            }
            precision = atoi(arg.c_str());
            if (15<precision) precision=15;
            if (precision<1)  precision=1;
            ++argNum;
        }
        else if (arg=="-msa-out") {
            if (nextArg.empty()) {
                PROBLEM(arg + " should be followed by a file path");
            }
            msaOutputPath = nextArg;
            ++argNum;
        }
        else if (arg=="-dist-out") {
            if (nextArg.empty()) {
                PROBLEM(arg + " should be followed by a file path");
            }
            distanceOutputFilePath = nextArg;
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
        else if (arg=="-uncorrected") {
            correcting_distances = false;
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
        else if (arg=="-filter") {
            filter_problem_sequences = true;
        }
        else if (arg=="-alphabet") {
            if (nextArg.empty()) {
                PROBLEM(arg + " should be followed by a list of characters");
            }
            alphabet = nextArg;
            ++argNum;
        }
        else if (arg=="-unknown") {
            if (nextArg.empty()) {
                PROBLEM(arg + " should be followed by a list of characters");
            }
            unknown_chars = nextArg;
            ++argNum;
        }
        else if (arg=="-num") {
            numbered_names = true;
        }
        else if (arg=="-not-dna") {
            is_DNA = false;
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
    if (alphabet.empty() && is_DNA) {
        alphabet = "ACGT";
    }
    if (unknown_chars.empty()) {
        unknown_chars = ".~_-?N";
        unknown_char  = 'N';
    } else {
        unknown_char = unknown_chars[unknown_chars.length()-1];
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
        std::cerr << "Tree builder algorithm was unexpectedly null"
            << " (internal logic error)." << std::endl;
        return 1;
    }
    algorithm->setZippedOutput(isOutputZipped || endsWith(outputFilePath,".gz"));
    if (beSilent) {
        algorithm->beSilent();
    }
    algorithm->setPrecision(precision);
    Sequences  sequences;
    FlatMatrix m;
    bool succeeded = prepInput(alignmentFilePath, inputFilePath,
                               !algorithm->isBenchmark(),
                               distanceOutputFilePath,
                               sequences, m);
    if (succeeded) {
        succeeded = algorithm->constructTreeInMemory(m.getSequenceNames(),
                                                     m.getDistanceMatrix(),
                                                     outputFilePath);
    }
    else if (!inputFilePath.empty()) {
        succeeded = algorithm->constructTree(inputFilePath, outputFilePath);
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
