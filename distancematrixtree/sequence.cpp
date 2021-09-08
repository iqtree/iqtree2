//
//  sequence.cpp
//
//  Copyright (C) 2021, James Barbetti.
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


#include "sequence.h"

#include <iostream>
#include <sstream>      //for std::stringstream
#include <math.h>       //for log()

#include <utils/progress.h>        //for progress_display
#include <utils/hammingdistance.h> //for vectorHamingDistance()
#include <utils/gzstream.h>        //for pigzstream class
#include <utils/safe_io.h>         //for safeGetLine()

bool isNucleotideOrMissing(const char c) {
    return isalnum(c) ||  c == '-' || c == '?'|| 
           c == '.' || c == '*' || c == '~';
}

bool isOpeningBracket(const char c) {
    return c == '(' || c == '{';
}

bool processSequenceLine(const std::vector<int> &in_alphabet,
                         char unknown_char, std::string &sequence,
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
            while (it != line.end() && *it != ')' && *it != '}') {
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

double correctedDistance(double char_dist,  double chars_compared,
                         double num_states, double max_distance) {
    double obs_dist = ( 0 < chars_compared )
                    ? ( char_dist / chars_compared ) : 0.0;
    double z = num_states / ( num_states - 1.0 );
    double x = 1.0 - ( z * obs_dist );
    double d = (x<=0) ? max_distance : -log(x)/z;
    if (max_distance<=d) {
        d = max_distance;
    }
    return d;
}

double uncorrectedDistance(double char_dist,
                           double chars_compared,
                           double max_distance) {
    return (0<chars_compared)
        ? (char_dist/chars_compared) : 0.0;
}

void useNumberedNamesIfAskedTo(bool numbered_names, FlatMatrix& m) {
    if (numbered_names) {
        auto name_count = m.getSequenceNames().size();
        for (size_t i=0; i<name_count; ++i) {
            std::stringstream name;
            name << "A" << (i+1); //"A1", "A2", ...
            m.setSequenceName(i, name.str());
        }
    }
}

Sequence::Sequence(const std::string& seq_name)
    : name(seq_name), is_problematic(false) {}
size_t Sequence::sequenceLength()           const    { return sequence_data.size(); }
const char* Sequence::data()                const    { return sequence_data.data(); }
const std::string& Sequence::sequenceData() const    { return sequence_data; }
        std::string& Sequence::sequenceData()          { return sequence_data; }
const std::string& Sequence::getName()      const    { return name; }
void  Sequence::setName(const std::string& new_name) { name = new_name; }
void  Sequence::setName(const char* new_name)        { name = new_name; }
bool  Sequence::isProblematic()             const    { return is_problematic; }
void  Sequence::markAsProblematic()                  { is_problematic = true; }

Sequences::Sequences(bool number_names): numbered_names(number_names) {
}
bool Sequences::checkLastTwoSequenceLengths() const {
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
size_t Sequences::countOfProblematicSequences() {
    size_t count = 0;
    for (size_t i=0; i<size(); ++i) {
        if (at(i).isProblematic()) {
            ++count;
        }
    }
    return count;
}
std::string Sequences::getFormattedName(size_t i) {
    if (numbered_names) {
        std::stringstream number_name;
        number_name << "A" << (i+1); //"A1", "A2", ...
        return number_name.str();
    } else {
        return at(i).getName();
    }
}
intptr_t Sequences::getSize() const {
    return size();
}
const std::string& Sequences::getSequenceName(size_t i) const {
    return at(i).getName();
}
void Sequences::setSequenceName(size_t i, const std::string& new_name) {
    at(i).setName(new_name);
}
bool Sequences::loadSequencesFromFasta(const std::string& fastaFilePath,
                                       const std::string& alphabet, bool unknown_char,
                                       bool report_progress) {
    #if USE_GZSTREAM
    pigzstream    in(report_progress ? "fasta" : "");
    #else
    std::ifstream in;
    #endif
    in.open(fastaFilePath.c_str(), std::ios::binary | std::ios::in);
    if (!in.is_open()) {
        std::cerr << "Unable to open alignment file " 
            << fastaFilePath << std::endl;
        return false;
    }
    size_t line_num = 1;
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
            if (!checkLastTwoSequenceLengths()) {
                return false;
            }
            emplace_back(str);
            continue;
        }
        // read sequence contents
        else if (empty()) {
            //Skip over it.
        }
        else if (!processSequenceLine(in_alphabet, unknown_char,
                                      back().sequenceData(),
                                      line, line_num)) {
            return false;
        }
    }
    if (!checkLastTwoSequenceLengths()) {
        return false;
    }
    in.close();
    return true;
}

bool Sequences::loadSequencesFromPhylip(const std::string& phylipFilePath,
                                        const std::string& alphabet, bool unknown_char,
                                        bool report_progress) {
    #if USE_GZSTREAM
    pigzstream    in(report_progress ? "phylip" : "");
    #else
    std::ifstream in;
    #endif
    in.open(phylipFilePath.c_str(), std::ios::binary | std::ios::in);
    if (!in.is_open()) {
        std::cerr << "Unable to open alignment file " 
            << phylipFilePath << std::endl;
        return false;
    }
    size_t line_num = 1;
    std::vector<int> in_alphabet;
    in_alphabet.resize(256, 0);
    for (auto alpha=alphabet.begin(); alpha!=alphabet.end(); ++alpha) {
        in_alphabet[*alpha] = 1;
    }
    size_t num_sequences       = 0;
    size_t sequence_length     = 0;
    bool   have_read_names     = 0;
    size_t name_length         = 0; //Number of characters to use for sequence name
    size_t sequence_num        = 0; //Ordinal sequence # to read next

    for (; !in.eof(); ++line_num) {
        std::string line;
        safeGetLine(in, line);
        if (line_num == 1) {
            //Read the heder line
            std::stringstream linestream(line);
            linestream >> num_sequences;
            linestream >> sequence_length;
            if (num_sequences < 1 || sequence_length < 1 ) {
                in.close();
                std::cerr << "Number of sequences " << num_sequences
                          << " or Sequence length " << sequence_length
                          << " was invalid.";
                return false;
            }
            continue;
        }
        if (line == "") {
            if (sequence_num!=0 && sequence_num!=num_sequences) {
                in.close();
                std::cerr << "Too few sequences (" << sequence_num << ") specified "
                          << " before blank line, at line " << line_num
                          << ". Expected " << num_sequences << ".";
                return false;
            }
            have_read_names     = true;
            sequence_num        = 0;
            continue;            
        }
        if (!have_read_names) {
            processPhylipSequenceName(line_num, sequence_num, 
                                      line, name_length);
        }        
        sequence_num %= num_sequences;
        std::string& seq_string = at(sequence_num).sequenceData();
        if (!processSequenceLine(in_alphabet, unknown_char,
                                 seq_string, line, line_num) ||
            !validateInterleaving(phylipFilePath, line_num, sequence_num) ) {
            in.close();
            return false;
        }
        ++sequence_num;
    }
    in.close();
    return validateLoadFromPhylip(phylipFilePath, num_sequences, 
                                  sequence_length);
}

bool Sequences::validateInterleaving(const std::string& phylipFilePath,
                                     size_t line_num, size_t sequence_num) {                                        
    if (0<sequence_num) {
        //Should we be checking that interleaving is consistent?
        //Or... is this being too fussy?
        if (at(sequence_num-1).sequenceLength() !=
            at(sequence_num).sequenceLength() ) {
            std::cerr << "Inconsistent interleaving at line " << line_num
                        << " of phylip multi-sequence alignment " << phylipFilePath << "."
                        << "\nSequence " << (sequence_num) << " length "
                        << " was " << at(sequence_num-1).sequenceLength() 
                        << " but sequence " << (sequence_num+1) << " length "
                        << " was " << at(sequence_num).sequenceLength() << ".";
            return false;
        }
    }
    return true;
}

bool Sequences::validateLoadFromPhylip(const std::string& phylipFilePath,
                                       size_t num_sequences, size_t sequence_length) {
    size_t sequence_num;
    if (size() != num_sequences) {
        std::cerr << "Only read " << size() << " sequences from " << phylipFilePath << "."
                  << "\nExpected to read " << num_sequences << ".";
        return false;
    }
    sequence_num = 1;
    for (Sequence& seq : *this) {
        if (seq.sequenceLength() != sequence_length) {
            std::cerr << "In " << phylipFilePath << ", "
                      << " sequence " << sequence_num << " had length "
                      << " " << seq.sequenceLength() << ","
                      << " but expected length " << sequence_length;
            return false;
        }
        ++sequence_num;
    }
    return true;
}

bool Sequences::processPhylipSequenceName(int line_num, int sequence_num, 
                                          std::string& line, size_t& name_length) {
    auto line_length = line.length();
    if (sequence_num==0) {
        //Scan for first white space
        for (name_length = 0; name_length<line_length; ++name_length) {
            auto ch = line[name_length];
            if (ch==' ' || ch=='\t' || ch=='\r' || ch=='\n') {
                break;
            }
        }
        //Scan for first non-white space after that
        for (;name_length<line_length; ++name_length) {
            auto ch = line[name_length];
            if (ch!=' ' && ch!='\t') {
                break;
            }
        }
    } else {
        if (line_length<name_length) {
            std::cerr << "Sequence at line " << line_num << " did not have"
                      << " a name with the expected length (" << name_length << ").\n";
            return false;
        }
    }

    std::string name(line.substr(0, name_length));
    name.erase(name.find_last_not_of(" \n\r\t")+1);
    name.erase(0, name.find_first_not_of(" \n\r\t"));
    emplace_back(name);
    line = line.substr(name_length, line_length-name_length);
    return true;
}

bool Sequences::loadAlignment(const std::string& fastaFilePath,
                              const std::string& phylipFilePath,
                              const std::string& alphabet, char unknown_char,
                              bool report_progress, 
                              std::vector<char>& is_site_variant) {
    //Assumes: either fastFilePath or 
    if (!fastaFilePath.empty()) {
        if (!loadSequencesFromFasta(fastaFilePath, alphabet, 
                                    unknown_char, report_progress)) {
            return false;
        }
    } else if (!phylipFilePath.empty()) {
        if (!loadSequencesFromPhylip(phylipFilePath, alphabet, 
                                     unknown_char, report_progress)) {
            return false;
        }
    } else {
        std::cerr << "Alignment file format not recognized.\n";
        return false;
    }
    if (size()<2) {
        std::cerr << "Cannot calculate distance matrix for a matrix"
                  << " of only " << size() << " sequences.";
        return false;
    }
    std::vector<size_t> sequence_odd_site_count;
    {
        size_t seqLen   = front().sequenceLength();
        std::vector<SiteInfo> sites;
        sites.resize(seqLen);
        SiteInfo* siteData = sites.data();
        
        size_t seqCount = size();
        for (size_t s=0; s<seqCount; ++s) {
            const char* sequence = at(s).data();
            for (size_t i=0; i<seqLen; ++i) {
                siteData[i].handle(unknown_char, s, sequence[i]);
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

void SequenceLoader::setUpSerializedData() {
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
        uint64_t          unk = 0;
        for (int col=0; col<seqLen; ++col) {
            unk <<= 1;
            if (read_site[col] == unknown_char ) {
                ++unk;
            }
            if ((col&63)==63) {
                *write_unk = unk;
                ++write_unk;
                unk = 0;
            }
        }
        if (unk!=0) {
            *write_unk = unk;
        }
        #if USE_PROGRESS_DISPLAY
        if ((row%100)==0) {
            extract_progress += 100.0;
        }
        #endif
    }
    #if USE_PROGRESS_DISPLAY
    extract_progress.done();
    #endif
}

void SequenceLoader::getNumberOfStates() {
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

SequenceLoader::SequenceLoader(char unknown, bool isDNA, 
                               double maximum_distance,
                               Sequences& sequences_to_load, 
                               bool use_corected_distances,
                               int  precision_to_use, int compression, 
                               const std::string& output_format_to_use,
                               const std::vector<char>& site_variant,
                               bool report_progress_while_loading)
    : unknown_char(unknown), is_DNA(isDNA)
    , max_distance(maximum_distance)
    , correcting_distances(use_corected_distances)
    , output_format(output_format_to_use)
    , precision(precision_to_use)
    , compression_level(compression)
    , sequences(sequences_to_load)
    , is_site_variant(site_variant)
    , report_progress(report_progress_while_loading), unkLen(0)
    , buffer(nullptr), sequence_data(nullptr)
    , unk_buffer(nullptr), unknown_data(nullptr) {
    rank      = sequences.size();
    rawSeqLen = sequences.front().sequenceLength();
    seqLen    = 0;
    for (auto it=is_site_variant.begin(); it!=is_site_variant.end(); ++it) {
        seqLen += *it;
    }
    #if (0)
    std::cout << "Number of invariant sites " << (rawSeqLen-seqLen) << "\n";
    std::cout << "Number of variable sites " << (seqLen) << "\n";
    std::cout << "Total number of sites " << rawSeqLen << "\n";
    #endif
}

SequenceLoader::~SequenceLoader() {
    delete [] unknown_data;
    delete [] unk_buffer;
    delete [] sequence_data;
    delete [] buffer;
}

double SequenceLoader::getDistanceBetweenSequences(intptr_t row, intptr_t col) const {
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
            distance = correctedDistance(static_cast<double>(char_distance), 
                                         static_cast<double>(adjSeqLen), 
                                         static_cast<double>(num_states),
                                         max_distance);
        } else {
            distance = uncorrectedDistance(static_cast<double>(char_distance), 
                                           static_cast<double>(adjSeqLen),
                                           max_distance);
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

bool SequenceLoader::loadSequenceDistances(FlatMatrix& m) {
    m.setSize(rank);
    for (intptr_t row=0; row<rank; ++row) {
        m.addCluster(sequences[row].getName());
    }
    setUpSerializedData();
    getNumberOfStates();
    #if USE_PROGRESS_DISPLAY
    const char* task = report_progress ? "Calculating distances": "";
    progress_display progress( rank*(rank-1)/2, task );
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
        #if USE_PROGRESS_DISPLAY
        progress += (rank-row);
        #endif
    }
    #if USE_PROGRESS_DISPLAY
    progress.done();
    #endif
    return true;
}

bool SequenceLoader::writeDistanceMatrixToFile(bool numbered_names,
                                               const std::string& filePath) {
    setUpSerializedData();
    getNumberOfStates();

    #if USE_PROGRESS_DISPLAY
    bool   isTriangle = output_format.find("lower") != std::string::npos ||
                        output_format.find("upper") != std::string::npos;
    double halfIfTriangle = isTriangle ? 0.5 : 1.0;
    double calculations   = static_cast<double>(rank) 
                          * static_cast<double>(rank) * halfIfTriangle;
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
    useNumberedNamesIfAskedTo(numbered_names, m);

    m.writeToDistanceFile(output_format, precision,
                          compression_level, false,
                          filePath);
    return true;
}

SiteInfo::SiteInfo(): minState('\0'), maxState('\0'), unknownCount(0) {
}

void SiteInfo::handle(char unknown_char, size_t sequenceIndex, char state) {
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
