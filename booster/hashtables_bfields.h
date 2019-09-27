/*

BOOSTER: BOOtstrap Support by TransfER: 
BOOSTER is an alternative method to compute bootstrap branch supports 
in large trees. It uses transfer distance between bipartitions, instead
of perfect match.

Copyright (C) 2017 Frederic Lemoine, Jean-Baka Domelevo Entfellner, Olivier Gascuel

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef _HASHTABLES_BFIELDS_H_
#define _HASHTABLES_BFIELDS_H_

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "stats.h"
#include "externs.h" /* gives the extern declaration of ntax, actual number of taxa in the tree(s) dealt with */

/* here we implement bit arrays to store taxon IDs. A taxon ID is an integer, and thus an index in a large bit array.
   A bipartition (== a subset of all the taxa) is a bit array in which the taxa that are present are all the bits set to 1.
   To be efficient in terms of storing the bipartitions, it is essential to have a variable length for our large bitfields.
   The bitfields are allocated at runtime, when we know the value of ntax, the number of taxa in the tree.
*/

/* TYPE DEFINITIONS */

#define MAX_TAXON_ID	USHRT_MAX
typedef unsigned short Taxon_id;	/* this gives us room for at least 65,536 taxa in the tree, maybe more
					   (depending on implementation). Taxon id 0 IS VALID. We can tweak it further here. */


typedef unsigned long* bfield_t;	/* the bitfield type: a series of consecutive unsigned longs. */
#define chunksize (8 * sizeof(unsigned long))	/* number of bits in a bitfield chunk, e.g. sizeof(unsigned long) = 4 means that chunksize = 32 */
#define nbchunks_bitarray (ntax/chunksize + (ntax%chunksize != 0 ? 1 : 0)) /* euclidean division */
/* and then this value never changes, it is the size of a bitarray in longs for this number of taxa. */



typedef struct _id_hash_table_t_ {
    int num_items;		/* the true number of items (ids) stored in this bit field */
    bfield_t bitarray;	      	/* the bit field */
} id_hash_table_t;


/* FUNCTIONS */


/* on id hash tables */
id_hash_table_t* create_id_hash_table(int size);
id_hash_table_t* complement_id_hashtbl(id_hash_table_t* h, int nbtaxa);

int lookup_id(id_hash_table_t *hashtable, Taxon_id my_id);
int add_id(id_hash_table_t *hashtable, Taxon_id my_id);
int delete_id(id_hash_table_t *hashtable, Taxon_id my_id);
void clear_id_hashtable(id_hash_table_t *hashtable);
void fill_id_hashtable(id_hash_table_t *hashtable, int nb_taxa);
void complement_id_hashtable(id_hash_table_t *destination, const id_hash_table_t *source, int nb_taxa);
unsigned int bitCount (unsigned long value);
void update_id_hashtable(id_hash_table_t *source, id_hash_table_t *destination);
int equal_id_hashtables(id_hash_table_t *tbl1, id_hash_table_t *tbl2);
int complement_id_hashtables(id_hash_table_t *tbl1, id_hash_table_t *tbl2,int nb_taxa);
int equal_or_complement_id_hashtables(id_hash_table_t *tbl1, id_hash_table_t *tbl2, int total);
void free_id_hashtable(id_hash_table_t *hashtable);

id_hash_table_t* suffle_hash_table(id_hash_table_t *hashtable, int total);

void print_id_hashtable(FILE* stream, id_hash_table_t *hashtable, int nbtaxa);


#endif /* _HASHTABLES_BFIELDS_H_ */
