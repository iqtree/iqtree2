#ifndef TERRACES_OLD_H
#define TERRACES_OLD_H

#include <gmp.h>

#ifdef __cplusplus
#include <cstddef>
using std::size_t;
#define TERRACES_NOEXCEPT noexcept
extern "C" {
#else
#include <stddef.h>
#define TERRACES_NOEXCEPT
#endif

/*
 Error Codes

 the return value of our function is an error code
 we will define these together as the project proceeds, e.g.
  0: successful completion
 -1: problem parsing Newick tree
 -2: #species in Newick tree does not correspond to number of species in data matrix
 -3: entries in data matrix not either 0 or 1
 -4: less than 4 species in input tree
 -5: only one partition in data matrix
 -6: overflow in number of splits to process
 -7: no output file specified
 -8: input tree is not a binary tree
 -9: there is no root species in the data file (a species present in all partitions)
 -10: there is a species with no partition in the data file
 -11: conflict between the set flags; can't perform all actions simultaneously
 */

#define TERRACE_SUCCESS 0
#define TERRACE_NEWICK_ERROR -1
#define TERRACE_SPECIES_ERROR -2
#define TERRACE_MATRIX_ERROR -3
#define TERRACE_NUM_SPECIES_ERROR -4
#define TERRACE_NUM_PARTITIONS_ERROR -5
#define TERRACE_SPLIT_COUNT_OVERFLOW_ERROR -6
#define TERRACE_OUTPUT_FILE_ERROR -7
#define TERRACE_TREE_NOT_BINARY_ERROR -8
#define TERRACE_NO_ROOT_SPECIES_ERROR -9
#define TERRACE_SPECIES_WITHOUT_PARTITION_ERROR -10
#define TERRACE_FLAG_CONFLICT_ERROR -11
#define TERRACE_INTERNAL_ERROR -99
/* to be extended */

/* check for unused return values */
#if defined(_MSC_VER) && (_MSC_VER >= 1700)
#define CHECK_RESULT _Check_return_
#else
#define CHECK_RESULT __attribute__((warn_unused_result))
#endif

/* Argument to control output of terraceAnalysis function (ta_outspec) */

/**
 count unrooted trees on terrace
 */
#define TA_COUNT 1

/**
 print unrooted trees on terrace to file
 @TODO: should TA_ENUMERATE automatically imply TA_COUNT?
 @TODO: Yes it should!
 */
#define TA_ENUMERATE 2

/**
 just detect if the tree is on a terrace. this should run much quicker than TA_COUNT or
 TA_ENUMERATE,
 because we can brake off, as soon as we have found that thera are at least two trees
 on the terrace.
 @TODO: how the output should look like in this case?
 @TODO: return any integer in terraceSize large than 1 if the tree is on a terrace
 */
#define TA_DETECT 4

/**
 print trees on a terrace in compressed form using some external binary tree compression tool.
 optional, does not need to be implemented, only if you want to.
 */
#define TA_ENUMERATE_COMPRESS 8

/**
 take a maximum comprehensive subset of the partitions if no comprehensive taxon exists.
 This leads to an over-estimation of the terrace size, but works with every data set.
 */
#define TA_UPPER_BOUND 16

// data type containing data to be passed to the algorithm we want to implement

typedef struct {
	size_t numberOfSpecies;
	size_t numberOfPartitions;
	unsigned char* missingDataMatrix;
	const char** speciesNames;
	bool allocatedNameArray;
} missingData;

/**
 * Initialize missing data data type
 *
 * @param numberOfSpecies number of species in dataset
 *
 * @param numberOfPartitions number of partitions in dataset
 *
 * @param speciesNames list of species names in dataset, first entry correpsonds to first row in
 * missingDataMatrix etc.
 *
 * @return poitner to missing data data structure
 */

missingData* initializeMissingData(size_t numberOfSpecies, size_t numberOfPartitions,
                                   const char** speciesNames);

/**
 * Free missing data data structure
 *
 * @param m pointer to missing data data structure
 */

void freeMissingData(missingData* m);

/**
 * set entry in missing data matrix
 *
 * @param m pointer to missing data data structure
 *
 * @param speciesNumber species index
 *
 * @param partitionNumber partition index
 *
 * @param value value to be set
 */

void setDataMatrix(missingData* m, size_t speciesNumber, size_t partitionNumber,
                   unsigned char value);

/**
 * get entry from missing data matrix
 *
 * @param m pointer to missing data data structure
 *
 * @param speciesNumber species index
 *
 * @param partitionNumber partition index
 *
 * @return the value at the specified matrix position
 */

unsigned char getDataMatrix(const missingData* m, size_t speciesNumber, size_t partitionNumber);

/**
 * copy one dimensional array containing the missing data matrix to the matrix in the missing data
 * data type
 *
 * @param matrix one-dimensional
 *
 * @param m pointer to missing data data structure
 *
 */

void copyDataMatrix(const unsigned char* matrix, missingData* m);

/**
 * Function that tells us, given a tree, and a missing data matrix as well as its dimensions,
 * if the tree is on a terrace, how many trees there are on the terrace, it also prints trees on the
 * terrace
 * to file, if specified and might compress them as well.
 *
 * We might need to change the data type of variable terraceSize that is being written in this
 * function.
 * Why?
 *
 * @param[in] m struct containing missing data matrix, number of partitions, number of species, and
 * list of
 * species names
 *
 * @param[in] newickTreeString Unrooted strictly binary phylogenetic tree,
 * in Newick format for which we want to find out if it is on a terrace.
 * Denoted by T in the task specification
 *
 * @param[in] ta_outspec bit-masked integer as combination of TA_* constants to control the outputs
 *
 * @param[out] allTreesOnTerrace output file name for unrooted tree enumeration.
 * Trees should be displayed in standard unrooted Newick format, and you should print one tree per
 * line.
 *
 * qparam[out] terraceSize number of unrooted trees on the terrace
 *
 * @return TERRACE_SUCCESS on success, or an error code (see TERRACE_*) on failure
 */
CHECK_RESULT int terraceAnalysis(missingData* m, const char* newickTreeString, int ta_outspec,
                                 const char* allTreesOnTerraceFile, mpz_t terraceSize);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* TERRACES_OLD_H */
