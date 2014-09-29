//
// C++ Interface: alignment
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <vector>
#include "pattern.h"
#include "ncl/ncl.h"
#include "tools.h"


const char STATE_UNKNOWN = 126;
const char STATE_INVALID = 127;
const int NUM_CHAR = 256;

enum SeqType {
    SEQ_DNA, SEQ_PROTEIN, SEQ_BINARY, SEQ_MORPH, SEQ_MULTISTATE, SEQ_CODON, SEQ_UNKNOWN
};


#ifdef USE_HASH_MAP
typedef unordered_map<string, int> StringIntMap;
typedef unordered_map<string, double> StringDoubleHashMap;
typedef unordered_map<string, int> PatternIntMap;
//typedef map<string, int> PatternIntMap;
#else
typedef map<string, double> StringDoubleHashMap;
typedef map<string, int> PatternIntMap;
#endif

/**
Multiple Sequence Alignment. Stored by a vector of site-patterns

        @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class Alignment : public vector<Pattern> {
    friend class SuperAlignment;

public:

    /**
            constructor
     */
    Alignment();

    /**
            constructor
            @param filename file name
            @param sequence_type type of the sequence, either "BIN", "DNA", "AA", or NULL
            @param intype (OUT) input format of the file
     */
    Alignment(char *filename, char *sequence_type, InputType &intype);

    /**
            destructor
     */
    virtual ~Alignment();


    /****************************************************************************
            input alignment reader
     ****************************************************************************/

    /**
            add a pattern into the alignment
            @param pat the pattern
            @param site the site index of the pattern from the alignment
            @param freq frequency of pattern
            @return TRUE if pattern contains only gaps or unknown char. 
                            In that case, the pattern won't be added.
     */
    bool addPattern(Pattern &pat, int site, int freq = 1);


    /**
            read the alignment in NEXUS format
            @param filename file name
            @return 1 on success, 0 on failure
     */
    int readNexus(char *filename);

    int buildPattern(StrVector &sequences, char *sequence_type, int nseq, int nsite);

    /**
            read the alignment in PHYLIP format
            @param filename file name
            @param sequence_type type of the sequence, either "BIN", "DNA", "AA", or NULL
            @return 1 on success, 0 on failure
     */
    int readPhylip(char *filename, char *sequence_type);

    /**
            read the alignment in FASTA format
            @param filename file name
            @param sequence_type type of the sequence, either "BIN", "DNA", "AA", or NULL
            @return 1 on success, 0 on failure
     */
    int readFasta(char *filename, char *sequence_type);

    /**
            extract the alignment from a nexus data block, called by readNexus()
            @param data_block data block of nexus file
     */
    void extractDataBlock(NxsCharactersBlock *data_block);

    /**
     * un-group site-patterns, i.e., making #sites = #patterns and pattern frequency = 1 for all patterns
     */
    void ungroupSitePattern();


    /**
     * re-group site-patterns
     * @param groups number of groups
     * @param site_group group ID (0, 1, ...ngroups-1; must be continuous) of all sites
     */
    void regroupSitePattern(int groups, IntVector &site_group);


    /****************************************************************************
            output alignment 
     ****************************************************************************/
    SeqType detectSequenceType(StrVector &sequences);

    virtual char convertState(char state, SeqType seq_type);

    /** 
     * convert state if the number of states (num_states is known)
     * @param state input char to convert
     * @return output char from 0 to 0-num_states or STATE_INVALID or STATE_UNKNOWN
     */
    char convertState(char state);

    virtual void convertStateStr(string &str, SeqType seq_type);

	/**
	 * convert from internal state to user-readable state (e.g., to ACGT for DNA)
	 * Note: does not work for codon data
	 * @param state internal state code
	 * @return user-readable state
	 */
    char convertStateBack(char state);

    /**
	 * convert from internal state to user-readable state (e.g., to ACGT for DNA)
	 * Note: work for all data
	 * @param state internal state code
	 * @return user-readable state string
	 */
	string convertStateBackStr(char state);

	/**
            get alignment site range from the residue range relative to a sequence
            @param seq_id reference sequence
            @param residue_left (IN/OUT) left of range
            @param residue_right (IN/OUT) right of range [left,right)
            @return TRUE if success, FALSE if out of range
     */
    bool getSiteFromResidue(int seq_id, int &residue_left, int &residue_right);

    int buildRetainingSites(const char *aln_site_list, IntVector &kept_sites,
            bool exclude_gaps, bool exclude_const_sites, const char *ref_seq_name);

    void printPhylip(const char *filename, bool append = false, const char *aln_site_list = NULL,
    		bool exclude_gaps = false, bool exclude_const_sites = false, const char *ref_seq_name = NULL);

    void printFasta(const char *filename, bool append = false, const char *aln_site_list = NULL,
    		bool exclude_gaps = false, bool exclude_const_sites = false, const char *ref_seq_name = NULL);

    /**
            Print the number of gaps per site
            @param filename output file name
     */
    void printSiteGaps(const char *filename);

    /****************************************************************************
            get general information from alignment
     ****************************************************************************/

    /**
            @return number of sequences
     */
    inline int getNSeq() {
        return seq_names.size();
    }

    /**
            @return number of sites (alignment columns)
     */
    inline int getNSite() {
        return site_pattern.size();
    }

    /**
             @return number of patterns
     */
    inline int getNPattern() {
        return size();
    }

    inline int getPatternID(int site) {
        return site_pattern[site];
    }

    inline Pattern getPattern(int site) {
        return at(site_pattern[site]);
    }

    /**
     * @param pattern_index (OUT) vector of size = alignment length storing pattern index of all sites
     */
    virtual void getSitePatternIndex(IntVector &pattern_index) {
        pattern_index = site_pattern;
    }

    /**
     * @param freq (OUT) vector of site-pattern frequencies
     */
    virtual void getPatternFreq(IntVector &freq);

    /**
            @param i sequence index
            @return sequence name
     */
    string &getSeqName(int i);

    /**
            @param seq_name sequence name
            @return corresponding ID, -1 if not found
     */
    int getSeqID(string &seq_name);

    /**
            @return length of the longest sequence name
     */
    int getMaxSeqNameLength();

    /**
            check proper and undupplicated sequence names
     */
    void checkSeqName();

    /**
     * check identical sequences
     * @return the number of sequences that are identical to one of the sequences
     */
    int checkIdenticalSeq();

    /**
            Quit if some sequences contain only gaps or missing data
     */
	virtual void checkGappySeq(bool force_error = true);

	/**
	 * return a new alignment if some sequence is totally gappy, or this if all sequence are okey
	 */
	Alignment *removeGappySeq();

    /**
            @return TRUE if seq_id contains only gaps or missing characters
            @param seq_id sequence ID
     */
    bool isGapOnlySeq(int seq_id);

    virtual bool isSuperAlignment() {
        return false;
    }

    /****************************************************************************
            alignment general processing
     ****************************************************************************/

    /**
            extract sub-alignment of a sub-set of sequences
            @param aln original input alignment
            @param seq_id ID of sequences to extract from
            @param min_true_cher the minimum number of non-gap characters, true_char<min_true_char -> delete the sequence
     */
    void extractSubAlignment(Alignment *aln, IntVector &seq_id, int min_true_char);

    /**
            extract a sub-set of patterns
            @param aln original input alignment
            @param ptn_id ID of patterns to extract from
     */
    void extractPatterns(Alignment *aln, IntVector &ptn_id);

    /**
            extract a sub-set of patterns
            @param aln original input alignment
            @param ptn_freq pattern frequency to extract from
     */
    void extractPatternFreqs(Alignment *aln, IntVector &ptn_freq);

    /**
            create a non-parametric bootstrap alignment from an input alignment
            @param aln input alignment
            @param pattern_freq (OUT) resampled pattern frequencies if not NULL
            @param spec bootstrap specification of the form "l1:b1,l2:b2,...,lk:bk"
            	to randomly draw b1 sites from the first l1 sites, etc. Note that l1+l2+...+lk
            	must equal m, where m is the alignment length. Otherwise, an error will occur.
            	If spec == NULL, a standard procedure is applied, i.e., randomly draw m sites.
     */
    virtual void createBootstrapAlignment(Alignment *aln, IntVector* pattern_freq = NULL, const char *spec = NULL);

    /**
            resampling pattern frequency by a non-parametric bootstrap 
            @param pattern_freq (OUT) resampled pattern frequencies
            @param spec bootstrap specification, see above
     */
    virtual void createBootstrapAlignment(IntVector &pattern_freq, const char *spec = NULL);

    /**
            resampling pattern frequency by a non-parametric bootstrap
            @param pattern_freq (OUT) resampled pattern frequencies
            @param spec bootstrap specification, see above
     */
    virtual void createBootstrapAlignment(int *pattern_freq, const char *spec = NULL);

    /**
            create a gap masked alignment from an input alignment. Gap patterns of masked_aln 
                    will be superimposed into aln to create the current alignment object.
            @param aln input alignment
            @param masked_aln gappy alignment of the same size with aln
     */
    void createGapMaskedAlignment(Alignment *masked_aln, Alignment *aln);

    /**
	 * shuffle alignment by randomizing the order of sites
	 */
	virtual void shuffleAlignment();

	/**
            concatenate an alignment into the current alignment object
            @param aln an alignment of the same number of sequences and sequence names    
     */
    void concatenateAlignment(Alignment *aln);

    /**
            copy the input alignment into the current alignment object
            @param aln input alignment
     */
    void copyAlignment(Alignment *aln);

    /**
            extract a sub-set of sites
            @param aln original input alignment
            @param ptn_id ID of sites to extract from (starting from 0)
     */
    void extractSites(Alignment *aln, IntVector &site_id);

    /**
            extract a sub-set of sites
            @param aln original input alignment
            @param spec specification of positions, e.g. "1-100,101-200\2"
     */
    void extractSites(Alignment *aln, const char* spec);

    /****************************************************************************
            Distance functions
     ****************************************************************************/


    /**
            compute the observed distance (number of different pairs of positions per site) 
                    between two sequences
            @param seq1 index of sequence 1
            @param seq2 index of sequence 2
            @return the observed distance between seq1 and seq2 (between 0.0 and 1.0)
     */
    virtual double computeObsDist(int seq1, int seq2);

    /**
            @param seq1 index of sequence 1
            @param seq2 index of sequence 2
            @return Juke-Cantor correction distance between seq1 and seq2
     */
    double computeJCDist(int seq1, int seq2);

    /**
            abstract function to compute the distance between 2 sequences. The default return
            Juke-Cantor corrected distance.
            @param seq1 index of sequence 1
            @param seq2 index of sequence 2		
            @return any distance between seq1 and seq2
     */
    virtual double computeDist(int seq1, int seq2) {
        return computeJCDist(seq1, seq2);
    }


    /**
            write distance matrix into a file in PHYLIP distance format
            @param file_name distance file name
            @param dist_mat distance matrix
     */
    void printDist(const char *file_name, double *dist_mat);

    /**
            write distance matrix into a stream in PHYLIP distance format
            @param out output stream
            @param dist_mat distance matrix
     */
    void printDist(ostream &out, double *dist_mat);

    /**
            read distance matrix from a file in PHYLIP distance format
            @param file_name distance file name
            @param dist_mat distance matrix
            @return the longest distance
     */
    double readDist(const char *file_name, double *dist_mat);

    /**
            read distance matrix from a stream in PHYLIP distance format
            @param in input stream
            @param dist_mat distance matrix
     */
    double readDist(istream &in, double *dist_mat);


    /****************************************************************************
            some statistics
     ****************************************************************************/

    /**
            compute empirical state frequencies from the alignment
            @param state_freq (OUT) is filled with state frequencies, assuming state_freq was allocated with 
                    at least num_states entries.
     */
    virtual void computeStateFreq(double *state_freq);

    /**
     * Make all frequencies a little different and non-zero
     * @param stateFrqArr (IN/OUT) state frequencies
     */
    void convfreq(double *stateFrqArr);

    /**
	 * compute special empirical frequencies for codon alignment: 1x4, 3x4, 3x4C
	 * @param state_freq (OUT) is filled with state frequencies, assuming state_freq was allocated with
	 * at least num_states entries.
	 * @param freq either FREQ_CODON_1x4, FREQ_CODON_3x4, or FREQ_CODON_3x4C
	 * @param ntfreq (OUT) nucleotide frequencies, assuming of size 4 for F1x4 and of size 12 for F3x4.
	 */
	void computeCodonFreq(StateFreqType freq, double *state_freq, double *ntfreq);

	/**
            compute empirical rates between state pairs
            @param rates (OUT) vector of size num_states*(num_states-1)/2 for the rates
     */
    virtual void computeEmpiricalRate(double *rates);

    /**
            compute non-reversible empirical rates between state pairs
            @param rates (OUT) vector of size num_states*(num_states-1) for the rates
     */
    virtual void computeEmpiricalRateNonRev(double *rates);

    /**
            count the fraction of constant sites in the alignment, update the variable frac_const_sites
     */
    virtual void countConstSite();

    /**
     * @return unobserved constant patterns, each entry encoding for one constant character
     */
    string getUnobservedConstPatterns();

    /**
            @return the number of ungappy and unambiguous characters from a sequence
            @param seq_id sequence ID
     */
    int countProperChar(int seq_id);

    /**
            @return unconstrained log-likelihood (without a tree)
     */
    virtual double computeUnconstrainedLogL();

    /** either SEQ_BINARY, SEQ_DNA, SEQ_PROTEIN, SEQ_MORPH, or SEQ_CODON */
    SeqType seq_type;

    /**
            number of states
     */
    int num_states;

    /**
            fraction of constant sites
     */
    double frac_const_sites;

	/**
	 *  map from 64 codon to non-stop codon index
	 */
    char *non_stop_codon;

	/**
	 * For codon sequences: index of 61 non-stop codons to 64 codons
	 * For other sequences: NULL
	 */
	char *codon_table;

	/**
	 * For codon_sequences: 64 amino-acid letters for genetic code of AAA,AAC,AAG,AAT,...,TTT
	 * For other sequences: NULL
	 */
	char *genetic_code;



    /** Added by MA
            Compute the probability of this alignment according to the multinomial distribution with parameters determined by the reference alignment
            @param refAlign the reference alignment
            @param prob (OUT) the returned probabilty
		
            The probability is computed as follows:
            - From the reference alignment, we count the relative pattern frequencies p_1 ... p_k (sum = 1)
            - From THIS alignment, we have frequencies d_1 ... d_k (sum = len = nsite)
            - Prob(THIS | refAlign) = nsite!/(d_1! * ... * d_k!) product(p_i^d_i)
     */
    void multinomialProb(Alignment refAlign, double &prob);

    /** Added by MA
            Compute the probability of the `expected alignment' according to the multinomial distribution with parameters determined by the pattern's observed frequencies in THIS alignment.
            The `expected alignment' consists of patterns with log-likelihoods (under some model+tree) given in the input file (logLL).
            Note that order of the log-likelihoods in inputLL must corresponds to patterns in THIS alignment.

            @param inputLL the input patterns log-likelihood vector
            @param prob (OUT) the returned probability
     */
    void multinomialProb(DoubleVector logLL, double &prob);
    void multinomialProb(double *logLL, double &prob);

    /** Adapted from MA
            compute the probability of the alignment defined by pattern_freq given this alignment	
     */
    double multinomialProb(IntVector &pattern_freq);


protected:


    /**
            sequence names
     */
    vector<string> seq_names;

    /**
            Site to pattern index
     */
    IntVector site_pattern;

    /**
            hash map from pattern to index in the vector of patterns (the alignment)
     */
    PatternIntMap pattern_index;


    /**
            get the appearance for a state, helpful for ambigious states
            @param state the state index
            @param state_app (OUT) state appearance
     */
    void getAppearance(char state, double *state_app);

	/**
	 * special initialization for codon sequences, e.g., setting #states, genetic_code
	 * @param sequence_type user-defined sequence type
	 */
	void initCodon(char *sequence_type);

};


void extractSiteID(Alignment *aln, const char* spec, IntVector &site_id);


#endif
