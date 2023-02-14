//
//  DecentTree.cpp
//
//  Copyright (C) 2020-2022, James Barbetti.
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
#include <utils/argument.h>        //for Argument, ArgumentMap, et al.
#include <utils/progress.h>        //for progress_display::setProgressDisplay()
#include <utils/operatingsystem.h> //for getOSName
#include <utils/hammingdistance.h> //for hammingDistance
#include <utils/stringfunctions.h> //for contains
#include "starttree.h"       //for StartTree::Registry
#include "flatmatrix.h"      //for FlatMatrix
#include "distancematrix.h"  //for loadDistanceMatrixInto
#include "sequence.h"        //for 
#if USE_GZSTREAM
#include <utils/gzstream.h>
#endif

#define PROBLEM(x) { problems << x << ".\n"; }

namespace {
    bool   correcting_distances     = true;  //Use Jukes-Cantor distance correction
                                             //when calculating distances between
                                             //taxons (for an alignment).
    bool   is_DNA                   = true;  //Sequence data, for an alignment 
                                             //is (A,G,C,T) DNA.
    bool   numbered_names           = false; //True if all sequence names are to be
                                             //replaced with names based on the 
                                             //sequence number (e.g. for when 
                                             //sequences in an input alignment 
                                             //have very long names, or names 
                                             //that contain characters that will be
                                             //problematic for the distance matrix 
                                             //program that will be passed the
                                             //distance matrix decentTree calculates
                                             //from the alignment).
    bool   filter_problem_sequences = false; //When processing alignments to generate
                                             //a distance matrix, and there are pairs of
                                             //sequences for which the distance cannot
                                             //be calculated (because there aren't *any* 
                                             //comparable characters, remove the sequence
                                             //with the greatest number of unknown sites).
                                             //(if set).
    char   unknown_char             = 'N';   //The character that indicates an 
                                             //site with an unknown state 
                                             //(e.g. for DNA, unknown nucleotide).
    double max_distance             = 10.0;  //An upper bound on (corrected) distance.
                                             //(The upper bound on uncorrected distance 
                                             // is always: 1).
    int    precision                = 8;     //(in distance matrix or tree output files)
    int    compression_level        = 9;     //(if outputting gzipped files)
    std::string msaOutputPath;     //write .msa formatted version of .fasta (or other) input here
    std::string phylipOutputPath;  //write .phylip formatted version of .fasta (or other) input here
    std::string alphabet;      //defaults to ACGT
    std::string unknown_chars; //defaults to .~_-?N
    std::string format              = "square.interleaved";
    bool        linewrap_format     = false; //if true, when generating a Phylip output file
                                             //add line-feeds within long sequences
                                             //rather than interleaving them,
    bool        interleaved_format  = true;  //if true, when generating an MSA format
                                             //file (an input for FASTME)
                                             //generate the interleaved format 
                                             //(which is the only format FASTME accepts).
                                             //(has no effect unless -msa-out option
                                             //is requested).
                                             //Likewise, when generating a Phylip output file
                                             //(unless linewrap_format is true).
    size_t      interleaving_width = 60;     //how many sites to output per line, if
                                             //interleaved_format is true

    std::string stripName;              //characters to strip from names
    std::string nameReplace("_");       //characters to replace stripepd chars with, in names
    std::string truncateName;           //truncate names when you see one of these characters
                                        //e.g. to make IQTree happy, truncate at space " ".
};

/** Show a banner for decentTree */
void showBanner() {
    std::cout << "\ndecentTree for " << getOSName() << "\n";
    std::cout << "Based on algorithms (UPGMA, NJ, BIONJ) proposed by Sokal & Michener [1958],\n";
    std::cout << "Saitou & Nei [1987], Gascuel [1997] and [2009]\n";
    std::cout << "Incorporating (in NJ-R and BIONJ-R) techniques proposed by Simonson, Mailund, and Pedersen [2011]\n";
    std::cout << "Developed by Olivier Gascuel [2009], Hoa Sien Cuong [2009], James Barbetti [2020-22]\n";
    std::cout << "(To suppress this banner pass -no-banner)\n";
}

/** Show a summary of decentTree's command line syntax and options */
void showUsage() {
    std::cout << "Usage: decenttree (-fasta [fastapath] | -phylip [phypath])\n";
    std::cout << "       (-msa-out [msapath]) (-strip-name [stripped]) \n";
    std::cout << "       (-name-replace [reps]) (-truncate-name-at [chars])\n";
    std::cout << "       (-uncorrected) (-no-matrix) (-dist-out [distout]) (-out-format [shape])\n";
    std::cout << "       (-alphabet [states]) (-unknown [chars]) (-not-dna))\n";
    std::cout << "       -in [mldist] (-c [level]) (-f [prec]) -out [newick] -t [algorithm]\n";
    std::cout << "       (-nt [threads]) (-gz) (-no-banner) (-q)\n";
    std::cout << "Arguments in parentheses () are optional.\n";
    std::cout << "[msapath]    is the path of a .msa format file, to which the input\n";
    std::cout << "             alignment is to be written (in .msa format).\n";
    std::cout << "[fastapath]  is the path of a .fasta format file specifying an alignment\n";
    std::cout << "             of genetic sequences (file may be in .gz format)\n";
    std::cout << "[phypath]    is the path of a .phy format file specifying an alignment\n";
    std::cout << "             genetic sequences (file may be in .gz format)\n";
    std::cout << "             (by default, the character indicating an unknown state is 'N')\n";
    std::cout << "[stripped]   is a list of characters to replace in taxon names, e.g. \" /\"\n";
    std::cout << "[rep]        is a list of characters to replace them with e.g. \"_\"\n";
    std::cout << "             (may be shorter than [strippped]; if so first character is the default.\n";
    std::cout << "[distout]    is the path, of a file, into which the distance matrix is to be written\n";
    std::cout << "             (possibly in a .gz format)\n";
    std::cout << "[shape]      is the shape of a distance matrix output\n";
    std::cout << "             (square for square, upper or lower for triangular)\n";
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

    std::cout << StartTree::Registry::getInstance().getListOfTreeBuilders();
}

/**
 * @brief Load sequences into a distance matrix
 * 
 * @param sequences        (input) the sequences to load
 * @param is_site_variant  a vector, with a size equal to the number of sites,
 *                         which indicates (with a non-zero) value which sites
 *                         vary.
 * @param report_progress  whether to display a progress bar
 * @param m                (output) the distance matrix
 * @param alphabet         the alphabet to use (if blank will be determined by
 *                         reading the sequences)
 * 
 * @return true  (if the sequences could be loaded, their distances calculated)
 * @return false (if an error occurred)
 */

bool loadSequenceDistancesIntoMatrix(Sequences& sequences,
                                     const std::vector<char>& is_site_variant,
                                     bool report_progress, FlatMatrix& m,
                                     std::string& alphabet) {
    SequenceLoader loader(unknown_char, is_DNA, max_distance,
                          sequences, correcting_distances,
                          precision, compression_level, format, 
                          is_site_variant, report_progress);
    bool success = loader.loadSequenceDistances(m, alphabet);
    return success;
}

/**
 * @brief Write sequences as an MSA output file
 *
 * if interleaved_format is true, an interleaved MSA output format 
 * will be used.  The first interleaving_width sites for each sequence will
 * be listed, followed by the next interleaving_width, and so on.
 * 
 * if interleaved_format is false, each sequence (no matter how long)
 * will be written on a single line.
 * 
 * @param sequences the sequences to write to the output file
 * @param msaPath   the path to the file
 * @return true     if the sequences could be written to the file
 * @return false    if an error occurred
 */
bool writeMSAOutputFile(const Sequences& sequences, const std::string& msaPath)
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
            for (size_t pos=0; pos<len; pos+=interleaving_width) {
                for (size_t i=0; i<sequences.size(); ++i) {
                    const std::string& seq_data = sequences[i].sequenceData();
                    if (pos==0) {
                        out << sequences.getFormattedName(i);
                    } else {
                        out << "      ";
                    }
                    size_t posEnd = pos + interleaving_width;
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
    catch (const std::ios::failure& f) {
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

/**
 * @brief Write interleaved sequences to a Phylip output file stream.
 *        e.g.
 * 
 * Taxon1 CATGAGCAT
 * Taxon2 CATCATCAT
 * Taxon1 TACTACTAC
 * Taxon2 ACTACTACT
 * 
 * @param sequences the sequences to write to the output file stream
 * @param len       how long each sequence is (it is assumed all are 
 *                  of the same length, but this is not checked).
 * @param out       the output stream to append them to
 */
void writeInterleavedSequences(const Sequences& sequences, size_t len, 
                            std::ofstream& out) {
    size_t longest_name = sequences.getLengthOfLongestFormattedName();
    for (size_t pos=0; pos<len; pos+=interleaving_width) {
        for (size_t i=0; i<sequences.size(); ++i) {
            const std::string& seq_data = sequences[i].sequenceData();
            //Interleaved phylip always repeats sequence names  
            out << sequences.getFormattedName(i, longest_name);
            size_t posEnd = pos + interleaving_width;
            if (len<posEnd) {
                posEnd = len;
            }
            out << " " << seq_data.substr(pos, posEnd-pos) << "\n";
        }
        out << "\n";
    }
}

/**
 * @brief Write line-wrapped sequences to a Phylip output file stream
 * 
 * @param sequences the sequences to write to the output file stream
 * @param len       how long each sequence is (it is assumed all are 
 *                  of the same length, but this is not checked).
 * @param out       the output stream to append them to
 */
void writeLineWrappedSequences(const Sequences& sequences, size_t len, 
                               std::ofstream& out) {
    size_t longest_name = sequences.getLengthOfLongestFormattedName();
    for (size_t i=0; i<sequences.size(); ++i) {
        const std::string& seq_data = sequences[i].sequenceData();
        out << sequences.getFormattedName(i, longest_name);
        for (size_t pos=0; pos<len; pos+=interleaving_width) {
            size_t posEnd = pos + interleaving_width;
            if (len<posEnd) {
                posEnd = len;
            }
            out << " " << seq_data.substr(pos, posEnd-pos) << "\n";
        }
    }
}

/**
 * @brief Write sequences as an Phylip output file
 *
 * if linewrap_format is true, the Phylip sequences will be line-wrapped.
 * 
 * otherwise, if interleaved_format is true, an interleaved Phylip output format 
 * will be used.  The first interleaving_width sites for each sequence will
 * be listed, followed by the next interleaving_width, and so on.
 * 
 * otherwise, if both linewrap_format and interleaved_format are false, 
 * each sequence (no matter how long) will be written on a single line
 * (one line per sequence).
 * 
 * @param sequences the sequences to write to the output file
 * @param msaPath   the path to the file
 * @return true     if the sequences could be written to the file
 * @return false    if an error occurred
 */
bool writePhylipOutputFile(const Sequences& sequences, const std::string& outPath)
{
    //This is writing an interleaved Phylip file
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    bool success = false;
    try {
        out.open(outPath.c_str());
        size_t len = sequences.size()==0 ? 0 : sequences[0].sequenceLength();
        out << sequences.size() << " " << len << "\n";
        
        //Todo: The following is for the interleaved sequence format
        //      that fastme likes.
        if (linewrap_format) {
            writeLineWrappedSequences(sequences, len, out);
        }
        if (interleaved_format) {
            writeInterleavedSequences(sequences, len, out);
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
            << " Phylip format file: " << outPath << std::endl;
    }
    catch (...) {
        std::cerr << "Unexpected error trying to write"
            << " Phylip format file: " << outPath << std::endl;
    }
    out.close();
    return success;
}

/**
 * @brief Remove, from Sequences (and a FlatMatrix) holding the distances
 *        between those seuqences, all of the sequences that have been 
 *        taggged as problematic (the sequences for which isProblematic()
 *        returns true).
 *        Problematic sequences (which are sequences for which distances 
 *        to at least some other sequences cannot be calculated, because 
 *        they have no known characters in common), are detected by
 *        SequenceLoader::getDistanceBetweenSequences().
 * @param sequences the sequences in the alignment.
 * @param m         the FlatMatrix, containing the distances 
 *                  between the sequences in the alignment,
 *                  which will already have been calculated, by
 *                  SequenceLoader::getDistanceBetweenSequences().
 */
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
    size_t rNew = 0; //index in sequences *in the output*.
                     //r will be the index in old_sequences (the *input*).
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

/**
 * @brief Process input files, and (possibly) load a distance matrix 
 *        (either by reading a distance matrix file, or by calculating -
 *        uncorrected or corrected - hamming distances between sequences)
 * 
 * @param fastaFilePath  the file path for a fasta format sequence file 
 *                       (if one is set, blank if not)
 * @param phylipFilePath the file path for a phylip format sequence file
 *                       (if one is set, blank if not)
 * @param matrixInputFilePath    the file path for a distance matrix file
 * @param reportProgress         true if a progress bar is to be displayed
 * @param distanceOutputFilePath the file path to write a copy of the 
 *                               distance matrix to (if set, blank if not)
 * @param sequences      (output) the sequences
 * @param loadMatrix     true if the distance matrix is to be loaded
 * @param m              (output) the distance matrix to load
 * @return true   (if everything succeeds)
 * @return false  (if no input is specified, or there is any sort of error)
 *
 * @note If fastaFilePath and phylipFilePath are both non-blank, 
 *       fastaFilePath will be honoured, and phylipFilePath ignored.
 * @note if either fastaFilePath and phylipFilePath are set, 
 *       distanceOutputFilePath will be ignored.
 * @note if loadMatrix is false, and a distanceOutputFilePath is supplied,
 *       SequenceLoader's writeDistanceMatrixToFile is used to write
 *       distances, row by row, without loading a distance matrix.
 */
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
                            reportProgress, m, alphabet)) {
                return false;
            }
        }
        if (filter_problem_sequences) {
            removeProblematicSequences(sequences, m);
        }
        if (!msaOutputPath.empty()) {
            writeMSAOutputFile(sequences, msaOutputPath);
        }
        if (!phylipOutputPath.empty()) {
            writePhylipOutputFile(sequences, phylipOutputPath);
        }

        if (!distanceOutputFilePath.empty()) {
            if (loadMatrix) {
                useNumberedNamesIfAskedTo(numbered_names, m);
                return m.writeToDistanceFile(format, precision,
                                 compression_level, reportProgress,
                                 distanceOutputFilePath );
            }
            else {
                SequenceLoader loader(unknown_char, is_DNA, 
                                      max_distance, sequences, 
                                      correcting_distances, precision, 
                                      compression_level, format,
                                      is_site_variant, reportProgress);
                bool success = loader.writeDistanceMatrixToFile(numbered_names, alphabet,
                                                                distanceOutputFilePath);
                return success;
            }
        }
    } else if (!matrixInputFilePath.empty()) {
        if (loadDistanceMatrixInto(matrixInputFilePath,
                                    reportProgress, m)) {
            fixUpSequenceNames(truncateName, stripName, nameReplace, m);
            if (!distanceOutputFilePath.empty()) {
                return m.writeToDistanceFile(format, precision,
                                             compression_level, reportProgress,
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

/**
 * @brief    restrict the range of an (in/out) thing to a range
 * @tparam T the type (for which operator< and operator= must exist)
 * @param lo the minimum value for the range (if restrict_me is less than *lo*,
 *           restrict_me will be set to *lo*).
 * @param hi the maixmum value for the range (if *hi* is less than restrict_me,
 *           restrict_me will be set to *hi*).
 * @param restrict_me 
 *
 * @note it is assumed that lo <= hi. In effect:
 *       if hi<lo, restrict_me will be set to min(restrict_me,lo)...
 */
template <class T> void range_restrict(const T lo, const T hi, T& restrict_me) {
    if (restrict_me < lo) {
        restrict_me = lo;
    } else if (hi < restrict_me) {
        restrict_me = hi;
    }
}

/**
 * @brief Tracks the settings of command-line options passed to 
 *        decentTree
 */
class DecentTreeOptions {
public:
    std::stringstream problems;         //accumulates any error messages
    std::string algorithmName;          //the name of the distance matrix algorithm
    std::string fastaFilePath;          //file path of fasta file, to load equence alignment from
    std::string phylipFilePath;         //file path of phylip file, to load sequence alignment from
    std::string inputFilePath;          //file path of phylip distance matrix file
                                        //(square and triangular) formats are supported
    std::string outputFilePath;         //file path of newick tree format output file
    std::string distanceOutputFilePath; //file path of phylip distance matrix format output file

    bool isOutputZipped            = false; //true if output is to be zipped
    bool isOutputSuppressed        = false; //true if no output is to be written
                                            //e.g. if benchmarking phylogenetic inference
                                            //     algorithms against one another.
    bool isOutputToStandardOutput  = false; //caller asked for newick tree to go to std::cout
    bool isBannerSuppressed        = false; //true if the program's banner is to be suppressed
    int  threadCount               = 0;     //number of threads to use (0 is "unset", or...
                                            //as many as are available)
    bool beSilent                  = false; //true if minimal output is to be written
    bool beVerbose                 = false; //true if additional explanatory output is 
                                            //to be writen to standard output
    bool isMatrixToBeLoaded        = true;  //set to false if caller passes -no-matrix
    bool isTreeConstructionSkipped = false; //becomes true if phylogenetic inference (or 
                                            //"tree construction") is to be skipped.
                                            //(e.g. if decentTree is being used to generate
                                            //an MSA file for FastME, or a distance matrix 
                                            //input for some other distance matrix 
                                            //phylogenetic inference program)
    bool showBar                   = false; //true if a progress bar is to be displayed

    ArgumentMap arg_map;

    /**
     * @brief Construct a new Decent Tree Options object
     */

    DecentTreeOptions()
        : algorithmName(StartTree::Registry::getNameOfDefaultTreeBuilder()) {
        isOutputZipped            = false;
        isOutputSuppressed        = false;
        isOutputToStandardOutput  = false; //caller asked for newick tree to go to std::cout
        isBannerSuppressed        = false;
        threadCount               = 0;
        beSilent                  = false;
        showBar                   = false;
        isMatrixToBeLoaded        = true;  //set to false if caller passes -no-matrix
        isTreeConstructionSkipped = false;
        initializeArgumentMap();
    }

    /**
     * @brief Set up a map, between recognized arguments and string,
     *        integer, double, and boolean variables.
     *        When the command line is parsed, actual arguments will
     *        be compared with recognized arguments, in the map, and
     *        the corresponding std::string, int, double, or bool 
     *        variables set.
     */
    void initializeArgumentMap() {
        arg_map << new StringArgument("-fasta", "fasta file path",             fastaFilePath);
        arg_map << new StringArgument("-phylip", "phylip alignment file path", phylipFilePath);
        arg_map << new StringArgument("-in",    "distance matrix file path",   inputFilePath);
        arg_map << new StringArgument("-dist",  "distance matrix file path",   distanceOutputFilePath);
        arg_map << new IntArgument   ("-c",     "compression level between 1 and 9", 
                                    compression_level);
        arg_map << new IntArgument   ("-f",     "precision level between 4 and 15",
                                    precision);
        arg_map << new StringArgument("-out-format", "output format" 
                                    "(e.g. square, upper, or lower)", format);
        arg_map << new StringArgument("-msa-out",  "msa format file path", msaOutputPath);
        arg_map << new StringArgument("-aln-out",  "phylip alignment file path", phylipOutputPath);
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
        arg_map << new DoubleArgument("-max-dist",    "maximum distance", max_distance);
        arg_map << new IntArgument   ("-nt", "thread count", threadCount);
        arg_map << new SwitchArgument("-q",           beSilent,                 true);
        arg_map << new SwitchArgument("-filter",      filter_problem_sequences, true);
        arg_map << new StringArgument("-alphabet", "list of characters", alphabet);
        arg_map << new StringArgument("-unknown",  "list of characters", unknown_chars);
        arg_map << new SwitchArgument("-num",         numbered_names,           true);
        arg_map << new SwitchArgument("-not-dna",     is_DNA,                   false);
        arg_map << new SwitchArgument("-v",           beVerbose,                true);
        arg_map << new SwitchArgument("-bar",         showBar,                  true);
    }

    /**
     * @brief Process command line options. Restrict maximum distance,
     *        compression level, and precision, to sensible ranges.
     * 
     * @param argc the number of command line options
     * @param argv an array of command line options (string constants)
     *
     * @note  The only command-line argument that needs special 
     *        treatment is -t (which specifies a distance-matrix
     *        algorithm to be used for phylogenetic inference).
     */
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
                    PROBLEM(StartTree::Registry::getInstance().getListOfTreeBuilders());
                }
                ++argNum;
            }
            else {
                PROBLEM("Unrecognized command-line argument, " + arg);
                break;
            }
        }
        if (max_distance<=0) {
            PROBLEM("Maximum distance (" << max_distance << ") too low");
        }
        range_restrict(0, 9,  compression_level );
        range_restrict(1, 15, precision );
        format = string_to_lower(format);
        if (isOutputZipped && !contains(format, ".gz")) {
            //Ensure that distance file will be compressed
            format += ".gz"; 
        }
        if (isOutputToStandardOutput) {
            outputFilePath = "STDOUT";
        }
        isBannerSuppressed |= beSilent;
        if (!alphabet.empty() && is_DNA) {
            is_DNA = false;
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
    }

    /**
     * @brief Check command-line options for consistency.
     * 
     * @return true  (if they're okay)
     * @return false (if they're self-contradictory, e.g. more than one
     *                input file is specified).
     */
    bool checkCommandLineOptions() {
        if (inputFilePath.empty() && fastaFilePath.empty() && phylipFilePath.empty()) {
            PROBLEM("Input (mldist) file should be specified via -in [filepath.mldist]");
            PROBLEM("or alignment (fasta) file may be specified via -fasta [filepath.fasta]");
            PROBLEM("or alignment (phylip) file, via -phylip [filepath.phy]");
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

    /**
     * @brief Set the number of threads to use (if Open MP is supported)
     * @note  If the number of threads hasn't been set, or 0 has been requested,
     *        the maximum thread count will be used.     
     */
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

/**
 * @brief entry point for decenttree.cpp
 * 
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @return int (0 if the command-line arguments are valid and no errors are
 *             reported; 1 if the command-line arguments are invalid or an
 *             error occurs))
 */
int main(int argc, char* argv[]) {
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

/**
 * @brief  Obey the pre-parsed command line options.
 *         (i)  load input files (if instructed to) 
 *              and/or calculate distance matrix from
 *              alignment (if instructed to)
 *         (ii) run the indicated distance matrix 
 *              phylogenetic algorithm (unless told not to)
 *              either in-memory
 * 
 * @param  options the options read from the command line
 * @return int     1 if an error occurred, 0 otherwise
 */
int obeyCommandLineOptions(DecentTreeOptions& options) {

    #if USE_PROGRESS_DISPLAY
    progress_display::setProgressDisplay(options.showBar);
    #endif

    StartTree::BuilderInterface* algorithm = 
        StartTree::Registry::getTreeBuilderByName(options.algorithmName);
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
