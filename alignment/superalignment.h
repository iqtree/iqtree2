/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef SUPERALIGNMENT_H
#define SUPERALIGNMENT_H

#include "alignment.h"

/**
Super alignment representing presence/absence of sequences in
k partitions for a total of n sequences. It has the form:
		Site_1 Site_2 ... Site_k
Seq_1     1      0    ...   1
Seq_2     0      1    ...   0
...      ...
Seq_n     1      1    ...   0

Where (i,j)=1 means Seq_i is present in partition j, 0 otherwise

So data is binary.

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/

class SuperAlignment : public Alignment
{
public:
	/** constructor initialize from a supertree */
    SuperAlignment(Params &params);

	/** constructor initialize empty alignment */
    SuperAlignment();

    /** destructor */
    ~SuperAlignment();

    /**
        load partitions from program Params
        @param params program Params
     */
    void readFromParams(Params &params);
    
    /**
     initialize seq_names, taxon_index, buildPattern
     */
    virtual void init(StrVector *sequence_names = NULL);
    
    /** return that this is a super-alignment structure */
	virtual bool isSuperAlignment() { return true; }

    /** read partition model file */
    void readPartition(Params &params);
    
    /** read RAxML-style partition file */
    void readPartitionRaxml(Params &params);
    
    /** read partition model file in NEXUS format into variable info */
    void readPartitionNexus(Params &params);

    /** read partition as files in a directory */
    void readPartitionDir(string partition_dir, char *sequence_type, InputType &intype, string model, bool remove_empty_seq);

    /** read partition as a comma-separated list of files */
    void readPartitionList(string file_list, char *sequence_type, InputType &intype, string model, bool remove_empty_seq);

    void printPartition(const char *filename, const char *aln_file);
    void printPartition(ostream &out, const char *aln_file = NULL, bool append = false);

    void printPartitionRaxml(const char *filename);
    
    void printBestPartition(const char *filename);
    void printBestPartitionRaxml(const char *filename);

	/**
	 * create taxa_index from super-alignment to sub-alignment
	 * @param part index of sub-alignment
	 */
	void linkSubAlignment(int part);

	/**
	 * @param pattern_index (OUT) vector of size = alignment length storing pattern index of all sites
	 * the index of sites in 2nd, 3rd,... genes have to be increased by the number of patterns in previous genes
	 * so that all indices are distinguishable
	*/
	virtual void getSitePatternIndex(IntVector &pattern_index);

	/**
	 * @param freq (OUT) vector of site-pattern frequencies for all sub-alignments
	*/
	virtual void getPatternFreq(IntVector &pattern_freq);

    /**
     * @param[out] freq vector of site-pattern frequencies
     */
    virtual void getPatternFreq(int *freq);

    /**
        Print all site information to a file
        @param filename output file name
    */
    virtual void printSiteInfo(const char* filename);

    /**
     compute empirical substitution counts between state pairs
     @param normalize true to normalize row sum to 1, false otherwise
     @param[out] pair_freq matrix of size num_states*num_states
     @param[out] state_freq vector of size num_states
     */
    virtual void computeDivergenceMatrix(double *pair_freq, double *state_freq, bool normalize = true);

    /**
     perform matched-pair tests of symmetry of Lars Jermiin et al.
     @param[out] sym results of test of symmetry
     @param[out] marsym results of test of marginal symmetry
     @param[out] intsym results of test of internal symmetry
     @param out output stream to print results
     @param rstream random stream to shuffle alignment columns
     @param out_stat output stream to print pairwise statistics
     */
    virtual void doSymTest(size_t vecid, vector<SymTestResult> &sym, vector<SymTestResult> &marsym,
                           vector<SymTestResult> &intsym, int *rstream = NULL, vector<SymTestStat> *stats = NULL);

    /**
            extract sub-alignment of a sub-set of sequences
            @param aln original input alignment
            @param seq_id ID of sequences to extract from
            @param min_true_cher the minimum number of non-gap characters, true_char<min_true_char -> delete the sequence
            @param min_taxa only keep alignment that has >= min_taxa sequences
            @param[out] kept_partitions (for SuperAlignment) indices of kept partitions
     */
    virtual void extractSubAlignment(Alignment *aln, IntVector &seq_id, int min_true_char, int min_taxa = 0, IntVector *kept_partitions = NULL);

    /**
        extract a subset of partitions to form a new SuperAlignment object
        @param part_id vector of partition IDs
        @return new alignment containing only part_id partitions
     */
    SuperAlignment *extractPartitions(IntVector &part_id);

    /**
     remove a subset of partitions
     @param part_id vector of partition IDs
     */
    void removePartitions(set<int> &part_id);

    /**
     * remove identical sequences from alignment
     * @param not_remove name of sequence where removal is avoided
     * @param keep_two TRUE to keep 2 out of k identical sequences, false to keep only 1
     * @param removed_seqs (OUT) name of removed sequences
     * @param target_seqs (OUT) corresponding name of kept sequence that is identical to the removed sequences
     * @return this if no sequences were removed, or new alignment if at least 1 sequence was removed
     */
    virtual Alignment *removeIdenticalSeq(string not_remove, bool keep_two, StrVector &removed_seqs, StrVector &target_seqs);


    /*
        check if some state is absent, which may cause numerical issues
        @param msg additional message into the warning
        @return number of absent states in the alignment
    */
    virtual int checkAbsentStates(string msg);

	/**
		Quit if some sequences contain only gaps or missing data
	*/
	//virtual void checkGappySeq(bool force_error = true);

	/**
		create a non-parametric bootstrap alignment by resampling sites within partitions
		@param aln input alignment
		@param pattern_freq (OUT) if not NULL, will store the resampled pattern frequencies
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
        @param rstream random generator stream, NULL to use the global randstream
	*/
	virtual void createBootstrapAlignment(int *pattern_freq, const char *spec = NULL, int *rstream = NULL);

	/**
	 * shuffle alignment by randomizing the order of sites over all sub-alignments
	 */
	virtual void shuffleAlignment();

	/**
		compute the observed (Hamming) distance (number of different pairs of positions per site)
			between two sequences
		@param seq1 index of sequence 1
		@param seq2 index of sequence 2
		@return the observed distance between seq1 and seq2 (between 0.0 and 1.0)
	*/
	virtual double computeObsDist(int seq1, int seq2);

	/**
		compute the Juke-Cantor corrected distance between 2 sequences over all partitions
		@param seq1 index of sequence 1
		@param seq2 index of sequence 2		
		@return any distance between seq1 and seq2
	*/
	virtual double computeDist(int seq1, int seq2);

	/**
	 * print the super-alignment to a file
	 * @param filename
	 * @param append TRUE to append to this file, false to write new file
	 */
    virtual void printAlignment(InputType format, ostream &out, const char* file_name
                                , bool append = false, const char *aln_site_list = NULL
                                , int exclude_sites = 0, const char *ref_seq_name = NULL);

	/**
	 * print all sub alignments into files with prefix, suffix is the charset name
	 * @param prefix prefix of output files
	 */
	void printSubAlignments(Params &params);

    /**
     @param quartet ID of four taxa
     @param[out] support number of sites supporting 12|34, 13|24 and 14|23
     */
    virtual void computeQuartetSupports(IntVector &quartet, vector<int64_t> &support);
    
	/**
		@return unconstrained log-likelihood (without a tree)
	*/
	virtual double computeUnconstrainedLogL();

	/**
	 * @return proportion of missing data in super alignment
	 */
	double computeMissingData();

	/**
	 * build all patterns of super alignent from partitions and taxa_index
	 * it is in form of a binary alignment, where 0 means absence and 1 means presence
	 * of a gene in a sequence
	 */
	virtual void buildPattern();

    /**
            count the fraction of constant sites in the alignment, update the variable frac_const_sites
     */
    virtual void countConstSite();

    /**
     * 	@return number of states, if it is a partition model, return max num_states across all partitions
     */
    virtual int getMaxNumStates() {
    	return max_num_states;
    }

    /** order pattern by number of character states and return in ptn_order
    */
    virtual void orderPatternByNumChars(int pat_type);

	/**
		actual partition alignments
	*/
	vector<Alignment*> partitions;

	/**
		matrix represents the index of taxon i in partition j, -1 if the taxon is not present
	*/
	vector<IntVector> taxa_index;

	/** maximum number of states across all partitions */
	int max_num_states;

	/**
	 * concatenate subset of alignments
	 * @param ids IDs of sub-alignments
	 * @return concatenated alignment
	 */
    Alignment *concatenateAlignments(set<int> &ids);

	/**
	 * concatenate all alignments
	 * @return concatenated alignment
	 */
    Alignment *concatenateAlignments();


};

#endif
