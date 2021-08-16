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

bool isNucleotideOrMissing(const char c) {
    return isalnum(c) ||  c == '-' || c == '?'|| 
           c == '.' || c == '*' || c == '~';
}

bool isOpeningBracket(const char c) {
    return c == '(' || c == '{';
}

bool processSequenceLine(const std::vector<int> &in_alphabet,
                         std::string &sequence,
                         std::string &line, size_t line_num) {
    //Note: this is based on processSeq from IQTree's alignment/ alignment.cpp
    //(except it returns false rather than throwing exceptions), and writes
    //errors to std::cerr.
    for (auto it = line.begin(); it != line.end(); ++it) {
        auto c = toupper(*it);
        if (c  <= ' ') {
            continue;
        }
        if (isNucleotideOrMissing(c)) {
            if (!in_alphabet[c]) {
                c = unknown_char;
            }
            sequence.append(1, c);
        }
        else if (isOpeningBracket(c)) {
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

void useNumberedNamesIfAskedTo(FlatMatrix& m) {
    if (numbered_names) {
        auto name_count = m.getSequenceNames().size();
        for (size_t i=0; i<name_count; ++i) {
            std::stringstream name;
            name << "A" << (i+1); //"A1", "A2", ...
            m.setSequenceName(i, name.str());
        }
    }
}

struct Sequence {
protected:
    std::string name;
    std::string sequence_data;
    bool        is_problematic;

public:
    explicit Sequence(const std::string& seq_name)
        : name(seq_name), is_problematic(false) {}
    size_t sequenceLength()           const    { return sequence_data.size(); }
    const char* data()                const    { return sequence_data.data(); }
    const std::string& sequenceData() const    { return sequence_data; }
          std::string& sequenceData()          { return sequence_data; }
    const std::string& getName()      const    { return name; }
    void  setName(const std::string& new_name) { name = new_name; }
    void  setName(const char* new_name)        { name = new_name; }
    bool  isProblematic()             const    { return is_problematic; }
    void  markAsProblematic()                  { is_problematic = true; }
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
    intptr_t getSize() const {
        return size();
    }
    const std::string& getSequenceName(size_t i) const {
        return at(i).getName();
    }
    void setSequenceName(size_t i, const std::string& new_name) {
        at(i).setName(new_name);
    }
};

class SequenceLoader {
protected:
    Sequences&               sequences;
    const std::vector<char>& is_site_variant;
    bool                     report_progress;
    intptr_t rank;
    size_t   rawSeqLen;
    size_t   seqLen    = 0; //# of characters that actually vary between two sequences

    //Serialized data drawn from the sequences (N=rank, P=number of variable sites)
    size_t     unkLen        ; //Call this U
    char*      buffer        ; //All of the sequence data, back to back (NP bytes)
    char**     sequence_data ; //Array of N pointers into buffer (one per sequence)
    uint64_t*  unk_buffer    ; //An array of bits, indicating which sites are unknown (NU bytes)
    uint64_t** unknown_data  ; //Array of N pointers into unk_buffer (one per sequence)
    double     num_states    ; //Number of states

    void setUpSerializedData() {
        if (unknown_data!=nullptr) {
            return; //It's already been set up.
        }
        unkLen        = ((seqLen+255)/256)*4;
        buffer        = new char      [ seqLen * rank ];
        sequence_data = new char*     [ rank ];
        unk_buffer    = new uint64_t  [ unkLen * rank ];
        unknown_data  = new uint64_t* [ rank ];
        memset(unk_buffer, 0, unkLen * rank);
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

    void getNumberOfStates() {
        //Determine the number of states (needed for correcting distances)
        num_states = 0.0;
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
    }

public:
    SequenceLoader(Sequences& sequences_to_load,
                          const std::vector<char>& site_variant,
                          bool report_progress_while_loading)
        : sequences(sequences_to_load), is_site_variant(site_variant)
        , report_progress(report_progress_while_loading), unkLen(0)
        , buffer(nullptr), sequence_data(nullptr)
        , unk_buffer(nullptr), unknown_data(nullptr) {
        rank      = sequences.size();
        rawSeqLen = sequences.front().sequenceLength();
        seqLen    = 0;
        for (auto it=is_site_variant.begin(); it!=is_site_variant.end(); ++it) {
            seqLen += *it;
        }
    }

    ~SequenceLoader() {
        delete [] unknown_data;
        delete [] unk_buffer;
        delete [] sequence_data;
        delete [] buffer;
    }

    inline double getDistanceBetweenSequences(intptr_t row, intptr_t col) const {
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
                distance = correctDistance(static_cast<double>(char_distance), 
                                           static_cast<double>(adjSeqLen), 
                                           static_cast<double>(num_states));
            } else {
                distance = uncorrectedDistance(static_cast<double>(char_distance), 
                                               static_cast<double>(adjSeqLen));
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
        return distance;    
    }

    bool loadSequenceDistances(FlatMatrix& m) {
        m.setSize(rank);
        for (intptr_t row=0; row<rank; ++row) {
            m.addCluster(sequences[row].getName());
        }
        setUpSerializedData();
        getNumberOfStates();
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
                double distance = getDistanceBetweenSequences(row, col);
                m.cell(row, col) = distance;
                m.cell(col, row) = distance;
            }
            progress += (rank-row);
        }
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
        return true;
    }

    bool writeDistanceMatrixToFile(const std::string& filePath) {
        setUpSerializedData();
        getNumberOfStates();
        bool   isTriangle = format.find("lower") != std::string::npos ||
                            format.find("upper") != std::string::npos;
        double halfIfTriangle = isTriangle ? 0.5 : 1.0;
        double calculations   = static_cast<double>(rank) 
                              * static_cast<double>(rank) * halfIfTriangle;

        #if USE_PROGRESS_DISPLAY
        const char* task = report_progress ? "Writing distance matrix file": "";
        progress_display progress(calculations, task );
        #else
        progress_display progress = 0.0;
        #endif

        class FakeMatrix : public FlatMatrix {        
            //This pretends to be a flat matrix.
        protected:
            SequenceLoader&   owner;
            progress_display& show_progress;
        public:
            typedef FlatMatrix super;
            FakeMatrix(SequenceLoader& my_owner, progress_display& progress_bar)
                : super(), owner(my_owner), show_progress(progress_bar) {}
            virtual void setSize(intptr_t rows) { rowCount=rows;}
            virtual void appendRowDistancesToLine(intptr_t nseqs,    intptr_t seq1, 
                                                  intptr_t rowStart, intptr_t rowStop,
                                                  std::stringstream& line) const {
                std::vector<double> distance_row_vector;
                distance_row_vector.resize(rowStop-rowStart, 0);
                double*  distance_row = distance_row_vector.data();
                intptr_t column_count = rowStop - rowStart;
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic)
                #endif
                for (intptr_t column = 0; column < column_count; ++column) {
                    intptr_t seq2      = column + rowStart;
                    distance_row[seq2] = owner.getDistanceBetweenSequences(seq1, seq2);
                }
                for (intptr_t column = 0; column < column_count; ++column) {
                    if (distance_row[column] <= 0) {
                        line << " 0";
                    } else {
                        line << " " << distance_row[column];
                    }
                }
                show_progress += static_cast<double>(column_count);
            }
        } m(*this, progress);
        m.setSize(rank);
        for (intptr_t row=0; row<rank; ++row) {
            m.addCluster(sequences[row].getName());
        }
        useNumberedNamesIfAskedTo(m);

        m.writeToDistanceFile(format, precision,
                              compression_level,
                              filePath);
        return true;
    }
};

bool loadSequenceDistancesIntoMatrix(Sequences& sequences,
                                     const std::vector<char>&   is_site_variant,
                                     bool report_progress, FlatMatrix& m) {
    SequenceLoader loader(sequences, is_site_variant, report_progress);
    bool success = loader.loadSequenceDistances(m);
    return success;
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

bool loadSequences(const std::string& alignmentFilePath,
                   bool report_progress, Sequences &sequences) {
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
    for (; !in.eof(); ++line_num) {
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
    return true;
}

bool loadAlignment(const std::string& alignmentFilePath,
                   bool report_progress, Sequences &sequences,
                   std::vector<char>& is_site_variant)
{
    if (!loadSequences(alignmentFilePath, report_progress, 
                       sequences)) {
        return false;
    }
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

void setUpReplacementArray(const std::string& chars_to_strip, 
                           const std::string& replacement_chars,
                           char *in_char_to_out_char /*array of 256*/) {
    unsigned char ch=0;
    for (unsigned int c=0; c<256; ++c, ++ch) {
        in_char_to_out_char[c] = ch;
    }
    size_t strip_count = chars_to_strip.length();
    size_t rep_count   = replacement_chars.length();
    for (size_t i=0; i<strip_count; ++i) {
        char   ch_in  = chars_to_strip[i];
        size_t ix_in  = (unsigned char)(ch_in);
        size_t j      = ( (i<rep_count) ? i : 0);
        char   ch_out = replacement_chars[j];
        in_char_to_out_char[ix_in] = ch_out;
    }
}

template <class S>
void truncateSequenceNames(const char* truncation_chars,
                           S&          sequences) {
    intptr_t seq_count = sequences.getSize();
    IntVector trunc(256, 0);
    for (; *truncation_chars; ++truncation_chars) {
        int ch = (static_cast<int>(*truncation_chars) & 255);
        trunc[ch] = 1;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t seq=0; seq<seq_count; ++seq) {
        std::string old_name = sequences.getSequenceName(seq);
        auto it = old_name.begin();
        for (; it != old_name.end(); ++it) {
            int ch = (static_cast<int>(*it) & 255);
            if (trunc[ch] !=0 ) {
                break;
            }
        }
        if (it != old_name.end()) {
            std::string new_name = old_name.substr(0, it-old_name.begin());
            sequences.setSequenceName(seq, new_name);
        }
    }
    //What about duplicated sequence names?
    //Might as well warn the caller that it's happened
}

template <class S>
void fixUpSequenceNames(const std::string& truncation_chars,
                        const std::string& chars_to_strip, 
                        const std::string& replacement_chars,
                        S&                 sequences) {
    if (!truncation_chars.empty()) {
        truncateSequenceNames(truncation_chars.c_str(), sequences);
    }
    if (chars_to_strip.empty() || replacement_chars.empty()) {
        return;
    }
    char in_char_to_out_char[256];
    setUpReplacementArray(chars_to_strip, replacement_chars, in_char_to_out_char);

    intptr_t seq_count = sequences.getSize();
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t seq=0; seq<seq_count; ++seq) {
        std::string old_name = sequences.getSequenceName(seq);
        std::string new_name = old_name;
        for (size_t i=0; i<new_name.length(); ++i) {
            char ch_in   = old_name[i];
            size_t ix_in = (unsigned char)(ch_in);
            new_name[i]  = in_char_to_out_char[ix_in];
        }
        if (new_name!=old_name) {
            sequences.setSequenceName(seq, new_name);
        }
    }
}

bool prepInput(const std::string& alignmentInputFilePath,
               const std::string& matrixInputFilePath,
               bool  reportProgress,
               const std::string& distanceOutputFilePath,
               Sequences& sequences, bool loadMatrix,
               FlatMatrix& m) {
    if (!alignmentInputFilePath.empty()) {
        std::vector<char> is_site_variant;
        if (!loadAlignment(alignmentInputFilePath, reportProgress,
            sequences, is_site_variant)) {
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
                useNumberedNamesIfAskedTo(m);
                return m.writeToDistanceFile(format, precision,
                                 compression_level,
                                 distanceOutputFilePath );
            }
            else {
                SequenceLoader loader(sequences, is_site_variant, reportProgress);
                bool success = loader.writeDistanceMatrixToFile(distanceOutputFilePath);
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

class Argument {
public:
    std::string name;
    explicit Argument(const char* arg_name): name(arg_name) {}
    virtual ~Argument() = default;
    virtual void accept(const std::string& arg, const std::string& nextArg, 
                        char* argv[], int argc, int &argNum,
                        std::stringstream& problems) = 0;
};

class StringArgument: public Argument {
public:
    typedef Argument super;
    std::string  description;
    std::string& mapped_to;
    StringArgument(const char* argument_name, 
                   const char* argument_description, 
                   std::string& variable) 
        : super(argument_name), description(argument_description),
          mapped_to(variable) { }
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, 
                int &argNum, std::stringstream& problems) override {
        ++argNum;
        if (argNum==argc) {
            problems << name << " should be followed by " << description << "\n";
        }
        mapped_to = nextArg;
    }
};

class IntArgument: public Argument {
    std::string description;
    int& int_var;
public:
    typedef Argument super;
    IntArgument(const char* arg_name, const char* desc, int& var) 
        : super(arg_name), description(desc), int_var(var) { }
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, int &argNum, 
                std::stringstream& problems) override {
        ++argNum;
        if (argNum==argc) {
            problems << name << " should be followed by " << description << ".\n";
        }
        int_var = atoi(nextArg.c_str());
    }
};

class SwitchArgument: public Argument {
private:
    bool& switch_var;
    bool  switch_setting;
public:
    typedef Argument super;
    SwitchArgument(const char* arg_name, bool& var, bool setting)
        : super(arg_name), switch_var(var), switch_setting(setting) { }
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, int &argNum,
                std::stringstream& problems) override {
        switch_var = switch_setting;
    }
};

class ArgumentMap: std::map<std::string, Argument*> {
public:
    ArgumentMap& operator << (Argument* arg) {
        insert(value_type(arg->name, arg));
        return *this;
    }
    ~ArgumentMap() {
        for (auto it=begin(); it!=end(); ++it) {
            delete it->second;
        }
        clear();
    }
    Argument* findByName(const std::string& name) {
        auto it = find(name);
        if (it == end()) {
            return nullptr;
        }
        return it->second;
    }
};

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
    std::string alignmentFilePath;      //only .fasta format is supported
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
        arg_map << new StringArgument("-fasta", "fasta file path", alignmentFilePath);
        arg_map << new StringArgument("-in",    "distance matrix file path", inputFilePath);
        arg_map << new StringArgument("-dist",  "distance matrix file path", inputFilePath);
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
        if (alignmentFilePath.empty() && !isMatrixToBeLoaded) {
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
        algorithm->setZippedOutput(options.isOutputZipped || endsWith(options.outputFilePath,".gz"));
        if (options.beSilent) {
            algorithm->beSilent();
        }
        algorithm->setPrecision(precision);
    }
    Sequences  sequences;
    FlatMatrix m;
    bool succeeded = prepInput(options.alignmentFilePath, options.inputFilePath,
                               !algorithm->isBenchmark(),
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
             !options.alignmentFilePath.empty() && 
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
