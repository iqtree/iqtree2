#pragma once
#ifndef distancematrix_sequence_h
#define distancematrix_sequence_h
#include <string>
#include <vector>
#include "flatmatrix.h" //for FlatMatrix class (used in SequenceLoader)

bool isNucleotideOrMissing(const char c);
bool isOpeningBracket     (const char c);
bool processSequenceLine  (const std::vector<int> &in_alphabet,
                           char unknown_char, std::string &sequence,
                           std::string &line, size_t line_num);
double correctedDistance  (double char_dist,  double chars_compared,
                           double num_states, double max_distance);
double uncorrectedDistance(double char_dist,  double chars_compared,
                           double max_distance);
void useNumberedNamesIfAskedTo(bool numbered_names, FlatMatrix& m);

class Sequence {
protected:
    std::string name;
    std::string sequence_data;
    bool        is_problematic;

public:
    explicit Sequence(const std::string& seq_name);
    size_t sequenceLength()           const;
    const char* data()                const;
    const std::string& sequenceData() const;
          std::string& sequenceData()      ;
    const std::string& getName()      const;
    void  setName(const std::string& new_name);
    void  setName(const char* new_name)    ;
    bool  isProblematic()             const;
    void  markAsProblematic()              ;
};                                

class Sequences: public std::vector<Sequence> {
    bool numbered_names;
public:
    typedef std::vector<Sequence> super;
    using super::size;

    Sequences() = delete;
    Sequences(const Sequences& rhs) = default;
    Sequences& operator=(const Sequences& rhs) = default;
    explicit Sequences(bool numbered_names);


    bool checkLastTwoSequenceLengths() const;
    size_t countOfProblematicSequences();
    std::string getFormattedName(size_t i);
    intptr_t getSize() const;
    const std::string& getSequenceName(size_t i) const;
    void setSequenceName(size_t i, const std::string& new_name);
    bool loadSequencesFromFasta   (const std::string& fastaFilePath,
                                   const std::string& alphabet,
                                   bool unknown_char, bool report_progress);
    bool loadSequencesFromPhylip  (const std::string& fastaFilePath,
                                   const std::string& alphabet,
                                   bool unknown_char, bool report_progress);
    bool processPhylipSequenceName(int line_num, int sequence_num, 
                                   std::string& line, size_t& name_length);
    bool validateInterleaving     (const std::string& phylipFilePath,
                                   size_t line_num, size_t sequence_num);
    bool validateLoadFromPhylip   (const std::string& phylipFilePath,
                                   size_t num_sequences, size_t sequence_length);
    bool loadAlignment            (const std::string& fastaFilePath,
                                   const std::string& phylipFilePath,
                                   const std::string& alphabet, char unknown_char,
                                   bool report_progress, 
                                   std::vector<char>& is_site_variant);
};

class SequenceLoader {
protected:
    char        unknown_char;
    bool        is_DNA;
    double      max_distance;
    bool        correcting_distances;
    std::string output_format;
    int         precision;
    int         compression_level;

    Sequences&               sequences;
    const std::vector<char>& is_site_variant;
    bool                     report_progress;
    intptr_t   rank;
    size_t     rawSeqLen;
    size_t     seqLen    = 0; //# of characters that actually vary between two sequences

    //Serialized data drawn from the sequences (N=rank, P=number of variable sites)
    size_t     unkLen        ; //Call this U
    char*      buffer        ; //All of the sequence data, back to back (NP bytes)
    char**     sequence_data ; //Array of N pointers into buffer (one per sequence)
    uint64_t*  unk_buffer    ; //An array of bits, indicating which sites are unknown (NU bytes)
    uint64_t** unknown_data  ; //Array of N pointers into unk_buffer (one per sequence)
    double     num_states    ; //Number of states

    void setUpSerializedData();
    void getNumberOfStates();

public:
    SequenceLoader(char unknown, bool isDNA, 
                   double maximum_distance,
                   Sequences& sequences_to_load,
                   bool use_corected_distances,
                   int  precision_to_use, int compression,
                   const std::string& output_format_to_use,
                   const std::vector<char>& site_variant,
                   bool report_progress_while_loading);
    ~SequenceLoader();
    double getDistanceBetweenSequences(intptr_t row, intptr_t col) const;
    bool   loadSequenceDistances      (FlatMatrix& m);
    bool   writeDistanceMatrixToFile  (bool numbered_names, 
                                       const std::string& filePath);
};

struct SiteInfo {
public:
    char minState;
    char maxState;
    size_t unknownCount;
    SiteInfo();
    void handle(char unknown_char, size_t sequenceIndex, char state);
};

/**
 * @brief Given a string listing characters to be stripped out / replaced
 *        and another listing those to replace them with, calculate
 *        a 256-element "map" such that indexing the map with an input
 *        character will yield the output character (which might be the same)
 * @param chars_to_strip    (input) list of characters to replace
 * @param replacement_chars (input) list of characters to replace them with
 * @param in_char_to_out_char (output) will be a "replacement" map
 * @note  if replacement_chars is shorter than chars_to_strip, characters
 *        that appear at an index, i, in chars_to_strip will be mapped to
 *        the character at index 0, in replacement_chars, or '-' if 
 *        replacement_chars is zero length.
 */
 inline void setUpReplacementArray(const std::string& chars_to_strip, 
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
        char   ch_out = ( (j<rep_count) ? replacement_chars[j] : '-');
        in_char_to_out_char[ix_in] = ch_out;
    }
}

/**
 * @brief  Truncate sequence names, in an alignment of sequences, 
 *         wherever the contain a character, that is found in a string
 *         containing the truncation characters.
 * @tparam S - the type of something that implements getSequenceName()
 *             and setSequenceName() member functions.
 * @param  truncation_chars - a string of characters to be truncated
 * @param  sequences        - the sequences
 * @note   Access to sequences, via getSequenceName() and getSequenceName()
 *         had better be thread-safe!  Because if _OPENMP is defined,
 *         this function parallelizes over sequences.
 * @note   Does NOT warn the caller/user if, after the truncations
 *         there are now duplicate sequence names. Should. Doesn't. -James B.
 */
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
}

/**
 * @brief  Change the names, of sequences, when those names might cause 
 *         trouble for some distance calculation or distance matrix 
 *         phylogenetic inference algorithm.
 * @tparam S - the type of something that implements getSequenceName()
 *             and setSequenceName() member functions.
 * @param  truncation_chars  - a string of characters to be truncated
 * @param  chars_to_strip    - a string of characters to strip
 * @param  replacement_chars - replacment characters, for those to strip
 * @param  sequences         - the sequences
 * @note   Access to sequences, via getSequenceName() and getSequenceName()
 *         had better be thread-safe!  Because if _OPENMP is defined,
 *         this function parallelizes over sequences.
 * @note   Does NOT warn the caller/user if, after the truncations
 *         there are now duplicate sequence names. Should. Doesn't. -James B.
 */
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

#endif

