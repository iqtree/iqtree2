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
#include <utils/stringfunctions.h> //for contains()

/**
 * @brief  Determines whether a character represents
 *         a nucleotide OR a position where one is 
 *         missing or unknown.
 * @param  c the character
 * @return true if it is
 * @return false if it snot
 */
bool isNucleotideOrMissing(const char c) {
    return isalnum(c) ||  c == '-' || c == '?'|| 
           c == '.' || c == '*' || c == '~';
}

/**
 * @brief  returns true if a character's an opening bracket
 *         (curly or smooth, but not square)
 * @param  c - the character
 * @return true if it is
 * @return false if not
 */
bool isOpeningBracket(const char c) {
    return c == '(' || c == '{';
}

/**
 * @brief Process a line containing characters representing characters, at
 *        sites, in a sequence (from an input file), and write characters 
 *        to a sequence (suitable for using elsewhere in the program).
 * @param in_alphabet  a vector of (at least) 256 ints
 *                     indicating whether a character is
 *                     considered to be in an alphabet
 *                     (nonzero if yes, zero if no).
 * @param unknown_char the character used to represent
 *                     an unknown nucleotide
 * @param sequence     the sequence being appended, with
 *                     the in-alphabet characters (plus
 *                     unknown_char).
 * @param line         the text of the line 
 * @param line_num     the line number in the file
 *                     (used in error messages)
 * @return true on success
 * @return false on failure
 * @note   ambiguous characters, represented with a list of possible 
 *         characters, concatenated together and enclosed in brackets 
 *         {} or () are all treated as entirely unkown.
 *         (This isn't good, but the sequence-reading functionality in 
 *          decentTree was added, primarily so that I could derive  
 *          distance matrices from interleaved alignment files, to 
 *          generate inputs that didn't have zero distances between 
 *          sequences, so I could test decentTree against programs that
 *          don't tolerate off-diagonal zeroes in the distance matrix,
 *          and so that I could fairly test it against other programs that
 *          expect interleaved phylip format alignment input files).
 */
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

/**
 * @brief  Determine a corrected distance between two sequences (given the
 *         parameters supplied - see below).
 * @param  char_dist      number of known characters that differ
 * @param  chars_compared number of characters that were compared (excludes 
 *                        sites that had an unknown character in either
 *                        of the sequences being distanced, in the alignment).
 * @param  num_states     the number of legal known states
 *                        (e.g. for DNA there are 4)
 * @param  max_distance   the maximum distance that is allowed
 *                        (if the corrected distance is more, it
 *                         will be set to max_distance)
 * @return double - the Jukes/Cantor corrected distance
 * @note   Generally, distance matrix algorithms work better if they are
 *         provided with uncorrected distances (Perhaps because the 
 *         "metric" is closer to a "euclidean" metric, and 3-taxon distance
 *         triplets are less likely to violate the triangle inequality?)
 * @note   it is assumed that num_states is a whole number greater than 1,
 *         but this is NOT checked.
 * @note   it is assumed char_dist, chars_compared, are non-negative,
 *         but this is NOT checked.
 * @note   if chars_compared is 0 the distance returned will be
 *         min(max_distance,0).
 */
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

/**
 * @brief  calculate the uncorrected distance between two sequences
 * @param  char_dist      number of known characters that differ
 * @param  chars_compared number of characters that were compared (excludes 
 *                        sites that had an unknown character in either
 *                        of the sequences being distanced, in the alignment).
 * @param  max_distance   the maximum distance that is allowed
 *                        (if the corrected distance is more, it
 *                         will be set to max_distance)
 * @return double - the uncorrected distance
 * @note   it is assumed that char_dist and chars_compared are non-negative,
 *         but this is NOT checked.
 * @note   if chars_compared is 0 the distance returned will be
 *         min(max_distance,0).
 */
double uncorrectedDistance(double char_dist,
                           double chars_compared,
                           double max_distance) {
    return (0<chars_compared)
        ? (char_dist/chars_compared) : 0.0;
}

/**
 * @brief Replace actual sequence names with numbered names, if asked to do so
 * @param numbered_names - true to do it, false not to
 * @param m - a reference to the FlatMatrix instance, in which the names
 *            might be replaced with numbered names.
 * @note  This function exists because it made it easier to test 
 *        decentTree against other programs that were fussier (some, much
 *        fussier) about the length of, or the legal characters in, sequence
 *        names.
 * @note  The names will be A1, A2, A3, ... through A[n] is the number 
 *        of sequences.
 */
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

/**
 * @brief  Check whether the current length, of the last two sequences
 *         is the same.
 * @return true  - if so
 * @return false - if not (note; an error message is also written to 
 *         std::cerr, if they are not the same length)
 * @note   This is used for both non-interleaved and interleaved alignment
 *         formats.  Even in interleaved formats, the number of characters
 *         added in each "column" of interleaving, should be the same for
 *         each of the sequences.  And this can be called.
 */
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

/**
 * @brief  Return the number of sequences that are "problematic".
 *         A sequence (with sequence index y) is problematic, if there is any 
 *         other sequence (index x), such that sequences x and y share no
 *         known characters, and (take a deep breath!) x has fewer unknown 
 *         characters than y, OR (x and y have the same number of unknown 
 *         characters, AND x<y).
 * @return size_t - the number of sequences marked as problematic
 * @note   The function that identifies problematic sequences is
 *         getDistanceBetweenSequences().
 */
size_t Sequences::countOfProblematicSequences() {
    size_t count = 0;
    for (size_t i=0; i<size(); ++i) {
        if (at(i).isProblematic()) {
            ++count;
        }
    }
    return count;
}

/**
 * @brief  read a name (or perhaps a numbered name) from an alignment,
 *         for a sequence (identified by sequence index)
 * @param  i - the sequence number
 * @return std::string - the name of sequence (or "A[i+1]"
 */
std::string Sequences::getFormattedName(size_t i) const {
    if (numbered_names) {
        std::stringstream number_name;
        number_name << "A" << (i+1); //"A1", "A2", ...
        return number_name.str();
    } else {
        return at(i).getName();
    }
}

/**
 * @brief  read a name (or perhaps a numbered name) from an alignment,
 *         for a sequence (identified by sequence index), and pad it
 *         (out to pad_length characters) if it is shorter than pad_length.
 * @param  i - the sequence number
 * @param  pad_length - the padding length
 * @return std::string - the name of sequence (or "A[i+1]"
 */
std::string Sequences::getFormattedName(size_t i, size_t pad_length) const {
    std::string result;
    if (numbered_names) {
        std::stringstream number_name;
        number_name << "A" << (i+1); //"A1", "A2", ...
        result = number_name.str();
    } else {
        result = at(i).getName();
    }
    if (pad_length<=result.size()) {
        return result;
    }
    result += std::string(pad_length-result.size(),' ');
    return result;    
}

size_t Sequences::getLengthOfLongestFormattedName() const {
    size_t longest_length = 0;
    for (int i=0; i<size(); ++i) {
        std::string s = getFormattedName(i);
        longest_length = std::max(longest_length, s.size());
    }
    return longest_length;
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

/**
 * @brief  Load a sequence alignment from a fasta format file
 * @param  fastaFilePath   - the path to the file
 * @param  alphabet        - the (nucleotide?) alphabet to use
 * @param  unknown_char    - the character that represents sites of 
 *                           unknown state.
 * @param  report_progress - true if progress is to be reported
 * @return true  - on success
 * @return false - on failure
 * @note   if the USE_GZSTREAM symbol is defined, and non-zero, the
 *         fasta file may be gzip compressed.  Otherwise it has to
 *         be uncompressed.
 */
bool Sequences::loadSequencesFromFasta(const std::string& fastaFilePath,
                                       const std::string& alphabet, 
                                       char  unknown_char,
                                       bool  report_progress) {
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
    in_alphabet.resize(256, alphabet.empty() ? 1 : 0);
    if (!alphabet.empty()) {
        for (auto alpha=alphabet.begin(); alpha!=alphabet.end(); ++alpha) {
            in_alphabet[*alpha] = 1;
        }
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

/**
 * @brief  Read the first line of a Phylip Alignment file (which is the number of sequences,
 *         white space, and the length of all the sequences)
 * 
 * @tparam S                the stream type
 * @param  in               the stream itself
 * @param  num_sequences    reference - will have the number of sequences
 *                          in the alignment assigned to it (on success)
 * @param  sequence_length  reference - will have the length of the sequences
 *                          (every sequence must have the same length)
 *                          assigned to it (on success)
 * @return true             on success
 * @return false            on failure (an error message will have been written
 *                          to standard output).
 * 
 * @note   it is assumed that failbit is not already set on the stream when this 
 *         function is called.
 */
template <class S> bool readFirstLineOfPhylipAlignmentFile(S& in, size_t& num_sequences, size_t& sequence_length) {
    if (in.eof()) {
        std::cerr << "Sequence file was empty.";
        return false;
    }
    std::string line;
    safeGetLine(in, line);
    //Read the header line
    std::stringstream linestream(line);
    linestream >> num_sequences;
    if (in.fail()) {
        std::cerr << "Could not read number of sequences.";
        return false;
    }
    linestream >> sequence_length;
    if (in.fail()) {
        std::cerr << "Could not read sequence length.";
        return false;
    }
    if (num_sequences < 1 || sequence_length < 1 ) {
        std::cerr << "Number of sequences " << num_sequences
                    << " or Sequence length " << sequence_length
                    << " was invalid.";
        return false;
    }
    return true;
}

/**
 * @brief  Load a sequence alignment from a phylip format file
 * @param  phylipFilePath   - the path to the file
 * @param  alphabet        - the (nucleotide?) alphabet to use
 * @param  unknown_char    - the character that represents sites of 
 *                           unknown state.
 * @param  report_progress - true if progress is to be reported
 * @return true  - on success
 * @return false - on failure
 * @note   if the USE_GZSTREAM symbol is defined, and non-zero, the
 *         phylip file may be gzip compressed.  Otherwise it has to
 *         be uncompressed.
 * @note   both un-interleaved and interleaved Phylip formats are
 *         supported.
 */
bool Sequences::loadSequencesFromPhylip(const std::string& phylipFilePath,
                                        const std::string& alphabet, 
                                        char  unknown_char,
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
    in_alphabet.resize(256, alphabet.empty() ? 1 : 0);
    if (!alphabet.empty()) {
        for (auto alpha=alphabet.begin(); alpha!=alphabet.end(); ++alpha) {
            in_alphabet[*alpha] = 1;
        }
    }
    size_t num_sequences       = 0;
    size_t sequence_length     = 0;
    if (!readFirstLineOfPhylipAlignmentFile(in, num_sequences, sequence_length))
    {
        in.close();
        return false;
    }

    bool   have_read_names     = 0;
    size_t name_length         = 0; //Number of characters to use for sequence name
    size_t sequence_num        = 0; //Ordinal sequence # to read next

    for (; !in.eof(); ++line_num) {
        std::string line;
        safeGetLine(in, line);
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
        else {
            //The start of the line might be a repeat of the sequence name
            std::string seq_prefix = getFormattedName(sequence_num) + " ";
            if (startsWith(line, seq_prefix.c_str())) 
            {
                line = line.substr(seq_prefix.size(), line.size()-seq_prefix.size());
            }
        }
        std::string& seq_string = at(sequence_num).sequenceData();
        if (!processSequenceLine(in_alphabet, unknown_char,
                                 seq_string, line, line_num) ||
            !validateInterleaving(phylipFilePath, line_num, sequence_num) ) {
            in.close();
            return false;
        }
        ++sequence_num;
        sequence_num %= num_sequences;
    }
    in.close();
    return validateLoadFromPhylip(phylipFilePath, num_sequences, 
                                  sequence_length);
}

/**
 * @brief  Check that interleaving in a phylip file is consistent (that is,
 *         each time the sequences are appended, the accumulated lengths of
 *         all the sequences are consistent).
 * @param  phylipFilePath - the name of the interleaved phylip format file
 *                          just read
 * @param  line_num       - the number of the line just read
 * @param  sequence_num   - the number of the sequence that was appended with
 *                          the state data read from the line
 * @return true  - on success
 * @return false - on failure (an error message quoting file path and line
 *         number will be written to std::cerr, before false is returned).
 */
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

/**
 * @brief  Check, after reading from a phylip format alignment file that
 *         (i)  the expected number of sequences were read
 *         (ii) all the sequences have the same length
 * @param phylipFilePath  - the path of the phylip format file
 * @param num_sequences   - the number of sequences
 * @param sequence_length - the length of the first sequence
 * @return true  - on success
 * @return false - on failure (the file path will be quoted in an error
 *         message written to std::cerr).
 */
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

/**
 * @brief  Examine a line, read from a Phylip format alignment file,
 * @param  line_num     - the line number in the file
 * @param  sequence_num - the index of the sequence, into which the line 
 *                        is being read and/or the one into which the
 *                        states, read from the line, are to be appended.
 * @param  line         - the text of the line (this is both input and
 *                        output; if the line has a sequence name in it,
 *                        the sequence name and trailing white space are
 *                        removed from it before it is returned to the caller.
 * @param  name_length  - output - holds the length of the sequence name
 * @return true  - on success
 * @return false - if the length of the sequence name was not as expected.
 * @note   (My understanding is that phylip sequence names are all supposed
 *         to be white-space padded out to the same length.  
 *         Hopefully that's always the case. -James B)
 */
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

/**
 * @brief  Load an alignment form a fasta or phylip format alignment file
 * @param  fastaFilePath   - the file path of a fasta format file (if set)
 * @param  phylipFilePath  - the file path of a phylip format file (if set)
 *                           (fastaFilePath and phylipFilePath shouldn't both
 *                            be non-blank; if they are, fastaFilePath will
 *                            be honoured, and phylipFilePath will be ignored)
 * @param  alphabet        - a string indicating the characters in the alphabet
 * @param  unknown_char    - the character that indicates "unknown"
 * @param  report_progress - true if progress is to be reported, false if not
 * @param  is_site_variant - a vector (of char!), is_site_variant[i] indicates
 *                           if the site with ordinal i is invariant.
 * @return true  - on success
 * @return false - on failure (error messages will be logged to std::cerr)
 *         both file path parameters, were blank, the file load failed,
 *         there were <2 sequences
 * @note   most of code is here is figuring out which sites are variant sites.
 *         but perhaps that ought to be a separate function?! -James B.
 */

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
        intptr_t seqLen   = front().sequenceLength();
        std::vector<SiteInfo> sites;
        sites.resize(seqLen);
        SiteInfo* siteData = sites.data();
        
        size_t seqCount = size();
        for (size_t s=0; s<seqCount; ++s) {
            const char* sequence = at(s).data();
            for (intptr_t i=0; i<seqLen; ++i) {
                siteData[i].handle(unknown_char, s, sequence[i]);
            }
        }
        
        is_site_variant.resize(seqLen, 0);
        sequence_odd_site_count.resize(seqCount, 0);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t i=0; i<seqLen; ++i) {
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

/**
 * @brief  Sets up the "serialized" data that is used for calculating
 *         distances.
 * @note   buffer - sequence characters, with sequences one after another,
 *         in row major order.
 * @note   sequence_data - an array of pointers into the data
 * @note   unk_buffer - a buffer, where the known/unknown status
 *         of characters are recorded in 32-bit (bitfield) integers.
 *         1 bits indicate unknown characters.
 * @note   unkLen - the number of 32-bit integers, per sequence,
 *         in unk_buffer.
 * @note  if _OPENMP is defined, this is parallelized oer sequences.
 */
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
        for (intptr_t col=0; col<(intptr_t)seqLen; ++col) {
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

/**
 * @brief Determine the number of states (needed for correcting distances)
 * @param alphabet - on entry the alphabet supplied by the user
 *                 - on exit, the same (if one was supplied, minus 
 *                   any duplicates), or the distinct characters found
 *                   in the input, if is_DNA is false, and a blank
 *                   alphabet was supplied.
 */
void SequenceLoader::getNumberOfStates(std::string& alphabet) {
    num_states = 0.0;
    if (is_DNA) {
        num_states = 4;
    }
    else {
        std::vector<size_t> char_counts;
        char_counts.resize(256, 0);
        auto char_count_array = char_counts.data();
        if (!alphabet.empty()) {
            //Determine what the distinct characters in the user-supplied
            //alphabet were.
            for (const unsigned char ch : alphabet) {
                ++char_count_array[ch];
            }
        }
        else
        {
            //Determine what the distinct characters in the input were
            //(if we were told, it wasn't DNA but we weren't supplied an alphabet)
            const unsigned char* start_buffer = reinterpret_cast<unsigned char*>(buffer);
            const unsigned char* end_buffer   = start_buffer + rank * seqLen;
            std::cout << "scanning: " << rank << " sequences of length " << seqLen << std::endl;
            for (const unsigned char* scan=start_buffer; scan<end_buffer; ++scan) {
                ++char_count_array[*scan];
            }
        }        
        alphabet.clear();
        for (int i=0; i<256; ++i) {
            if (i!=unknown_char && 0<char_counts[i]) {
                alphabet.push_back(static_cast<char>(i));
            }
        }        
        num_states = alphabet.length();        
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
    , unk_buffer(nullptr), unknown_data(nullptr), num_states(4) {
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

/**
 * @brief  Calculate pairwise distances, between a pair of sequences
 * @param  row - one sequence number (indicating row of the distance
 *               matrix for which we are determining distances)
 * @param  col - another sequence number (the other sequence)
 * @return double 
 * @note   it is assumed that row and col are valid sequence numbers.
 *         but there is NO check that row != col (as you might expect, 
 *         you'll get 0! But you won't get it cheaply, because that 
 *         special case does not get special treatment!)
 * @note   reads: the serialized sequence data in buffer, via sequence_data
 *         (the per sequence pointers into it).
 *         also reads: the "is it unknown" bitfield data, in unk_buffer 
 *         via unknown_data (the per-sequence pointers into it)
 * @note   may mark one of the sequences as problematic (if the distance
 *         between the sequences cannot be determined because they have 
 *         no known characters in common).
 */
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

/**
 * @brief  Calculate pairwise distances, for every pair of sequences
 * @param  m - reference to the FlatMatrix into which distances are to be
 *             written.
 * @param  alphabet - alphabet to use (if blank, will be determined 
 *                    from the sites of the sequences)
 * @return true  - always 
 * @return false - in theory, could return this if it failed (but it won't)
 *         (though it might throw an out of memory exception, I suppose)
 * @note   If _OPENMP is set parallelizes over rows.
 * @note   Actually only calculates the *upper* triangle, top-to-bottom
 *         and left-to-right, but writes to the *lower* triangle too.
 * @note   The upper triangle is calculated so that parallelization will
 *         (we may hope) "load balance better" (because the small bits of
 *         work that won't take so long, occur for the later rows).
 */
bool SequenceLoader::loadSequenceDistances(FlatMatrix& m, std::string& alphabet) {
    m.setSize(rank);
    for (intptr_t row=0; row<rank; ++row) {
        m.addCluster(sequences[row].getName());
    }
    setUpSerializedData();
    getNumberOfStates(alphabet);

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

/**
 * @brief  This writes a calculated distance matrix directly to a file,
 *         without an in-memory distance matrix (which you might want to
 *         do if you have a *hell* of a lot of sequences, and you're
 *         calculating a distance matrix on a memory-challenged machine).
 * @param  numbered_names - true if names are to be numbered,
 *                          false if the existing names are to be used
 * @param  alphabet       - can be passed in (if blank, alphabet will be
 *                          determined from the sequence character data)
 * @param  filePath       - the path of the (phylip format) output file
 * @return true  - on success
 * @return false - on failure (error messages will be written to std::cerr)
 * @note   the output_format is honoured (if it is "lower" or "upper").
 *         if it's neither of those it is assumed the output format is 
 *         "square".
 * @note   this works by setting up a MOCK FlatMatrix that doesn't have
 *         any memory, and "faking" an override implementation of its
 *         appendRowDistancesToLine() function.  
 * @note   if _OPENMP is set, in each row, distance matrix calculations
 *         are parallelized over the columns to be outputted in that row
 */
bool SequenceLoader::writeDistanceMatrixToFile(bool numbered_names,
                                               std::string& alphabet,
                                               const std::string& filePath) {
    setUpSerializedData();
    getNumberOfStates(alphabet);

    #if USE_PROGRESS_DISPLAY
    bool   isTriangle = contains(output_format,"lower") ||
                        contains(output_format,"upper");
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
        FakeMatrix(const FakeMatrix& rhs) = delete;
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

/**
 * @brief Update site information on the basis of the next character
 * @param unknown_char  - the character indicating an unknown state
 * @param sequenceIndex - the sequence number, in the alignment
 * @param state         - the state, for this site, in that sequence
 */
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
