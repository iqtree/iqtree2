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

enum SeqType {SEQ_DNA, SEQ_PROTEIN, SEQ_BINARY, SEQ_UNKNOWN};


#ifdef USE_HASH_MAP
/*
	Define the hash function of Split
*/
#if defined(WIN32) 
namespace stdext {
#else
namespace __gnu_cxx {
#endif

	template<>
	struct hash<string> {
		size_t operator()(string str) const {
			hash<const char*> hash_str;
			return hash_str(str.c_str());
		}
	};
} // namespace __gnu_cxx
#endif


#ifdef USE_HASH_MAP
typedef hash_map<string, int> PatternIntMap;
//typedef map<string, int> PatternIntMap;
#else
typedef map<string, int> PatternIntMap;
#endif
/**
Multiple Sequence Alignment. Stored by a vector of site-patterns

	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class Alignment : public vector<Pattern>
{
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
    Alignment(char *filename,  char *sequence_type, InputType &intype);

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


	/**
		read the alignment in PHYLIP format
		@param filename file name
		@param sequence_type type of the sequence, either "BIN", "DNA", "AA", or NULL
		@return 1 on success, 0 on failure
	*/
	int readPhylip(char *filename, char *sequence_type);

	/**
		extract the alignment from a nexus data block, called by readNexus()
		@param data_block data block of nexus file
	*/
    void extractDataBlock(NxsCharactersBlock *data_block);


/****************************************************************************
	output alignment 
****************************************************************************/
	SeqType detectSequenceType(StrVector &sequences);

	char convertState(char state, SeqType seq_type);
	void convertState(string &str, SeqType seq_type);

	char convertStateBack(char state);

	void printPhylip(const char *filename, bool append = false);

/****************************************************************************
	get general information from alignment
****************************************************************************/

	/**
		@return number of sequences
	*/
	inline int getNSeq() { return seq_names.size(); }


	/**
		@return number of sites (alignment columns)
	*/
	inline int getNSite() { return site_pattern.size(); }


	/**
		 @return number of patterns
	*/
	inline int getNPattern() { return size(); }

	inline int getPatternID(int site) { return site_pattern[site]; }

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

/****************************************************************************
	alignment general processing
****************************************************************************/

	/**
		extract sub-alignment of a sub-set of sequences
		@param aln original input alignment
		@param seq_id ID of sequences to extract from
	*/
	void extractSubAlignment(Alignment *aln, IntVector &seq_id, int min_true_char);

	/**
		extract a sub-set of patterns
		@param aln original input alignment
		@param ptn_id ID of patterns to extract from
	*/
	void extractPatterns(Alignment *aln, IntVector &ptn_id);

	/**
		TODO
		create a non-parametric bootstrap alignment from the current alignment
		@param boot_aln (OUT) created bootstrap alignment
	*/
	void createBootstrapAlignment(Alignment *boot_aln);

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
	double computeObsDist(int seq1, int seq2);

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
	virtual double computeDist(int seq1, int seq2) { return computeJCDist(seq1, seq2); }


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
	*/
	void readDist(const char *file_name, double *dist_mat);

	/**
		read distance matrix from a stream in PHYLIP distance format
		@param in input stream
		@param dist_mat distance matrix
	*/
	void readDist(istream &in, double *dist_mat);


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
		compute empirical rates between state pairs
		@param rates (OUT) vector of size num_states*(num_states-1)/2 for the rates
	*/
	virtual void computeEmpiricalRate (double *rates);

	/**
		compute non-reversible empirical rates between state pairs
		@param rates (OUT) vector of size num_states*(num_states-1) for the rates
	*/
	virtual void computeEmpiricalRateNonRev (double *rates);

	/**
		count the fraction of constant sites in the alignment, update the variable frac_const_sites
	*/
	virtual void countConstSite();

	/**
		@return unconstrained log-likelihood (without a tree)
	*/
	double computeUnconstrainedLogL();

	/**
		number of states
	*/
	int num_states;

	/**
		fraction of constant sites
	*/
	double frac_const_sites;


private:

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
};

#endif
