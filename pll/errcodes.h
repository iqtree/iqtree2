/** 
 * PLL (version 1.0.0) a software library for phylogenetic inference
 * Copyright (C) 2013 Tomas Flouri and Alexandros Stamatakis
 *
 * Derived from 
 * RAxML-HPC, a program for sequential and parallel estimation of phylogenetic
 * trees by Alexandros Stamatakis
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For any other enquiries send an Email to Tomas Flouri
 * Tomas.Flouri@h-its.org
 *
 * When publishing work that uses PLL please cite PLL
 * 
 * @file errcodes.h
 */
#ifndef ERRCODES_H
#define ERRCODES_H

#define PLL_ERROR_FILE_OPEN             1               /**< Error while opening file */
#define PLL_ERROR_INVALID_FILETYPE      2               /**< Invalid fileType given at pllParseAlignmeFile */

#define  PLL_NNI_P_TIP                  1 << 0          /**< Node p is a tip */
#define  PLL_NNI_Q_TIP                  1 << 1          /**< Node p->back is a tip */

#define  PLL_PARTITION_OUT_OF_BOUNDS    1 << 0      /**< Trying to access a partition index that is out of bounds */
#define  PLL_BASE_FREQUENCIES_DO_NOT_SUM_TO_1 1 << 1      /**< base frequencies don't sum to 1.0 */

#define PLL_LINKAGE_LIST_OUT_OF_BOUNDS 1 << 0      /**< trying to link a partition index that is out of bounds */

#define PLL_SUBSTITUTION_RATE_OUT_OF_BOUNDS 1 << 0 /**< trying  to set a substitution rate to a value that is out of bounds */
#define PLL_INVALID_Q_MATRIX_SYMMETRY       1 << 1 /**< specifyng an invalid parameter symmetry in the Q matrix */
#define PLL_Q_MATRIX_SYMMETRY_OUT_OF_BOUNDS 1 << 2 /**<specifying a Q matrix symmetry that is out of bounds */

#define PLL_UNKNOWN_MOLECULAR_DATA_TYPE 1 << 0 /**<PLL is trying to do something for an unknown data type */

#define PLL_INCONSISTENT_SUBST_RATE_OPTIMIZATION_SETTING 1 << 0 /**<PLL detected an inconsistent setting for the Q matrix rate optimization */
#define PLL_INCONSISTENT_Q_MATRIX_SYMMETRIES_ACROSS_LINKED_PARTITIONS 1 << 1 /**<Q matrix symmetry vector is not identical for linked partitions */
#define PLL_INCONSISTENT_Q_MATRIX_ENTRIES_ACROSS_LINKED_PARTITIONS 1 << 2 /**<Q matrix entries are not identical for linked partitions */
#define PLL_INCONSISTENT_ALPHA_STATES_ACROSS_LINKED_PARTITIONS 1 << 3 /**<alpha states are not identical across linked partitions */
#define PLL_INCONSISTENT_ALPHA_VALUES_ACROSS_LINKED_PARTITIONS 1 << 4 /**<alpha values are not identical across linked partitions */
#define PLL_INCONSISTENT_FREQUENCY_STATES_ACROSS_LINKED_PARTITIONS 1 << 5 /**<frequency states are not identical across linked partitions */
#define PLL_INCONSISTENT_FREQUENCY_VALUES_ACROSS_LINKED_PARTITIONS 1 << 6 /**<frequency values are not identical across linked partitions */

#define PLL_NEWICK_ROOTED_TREE          1 << 0          /**< @brief Binary root detected */
#define PLL_NEWICK_BAD_STRUCTURE        1 << 1          /**< @brief Errornous tree detected */



#define PLL_ERROR_PHYLIP_HEADER_SYNTAX         5
#define PLL_ERROR_PHYLIP_BODY_SYNTAX           6
#define PLL_ERROR_FASTA_SYNTAX                 7




#endif
