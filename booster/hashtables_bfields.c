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

/* This file implements bit arrays to store Taxon_ids, for use in the Edges of the Tree objects. */
#include "hashtables_bfields.h"
/* ntax is defined as an extern int in this header file.
   chunksize is also defined there. */


id_hash_table_t* create_id_hash_table(int size)
{
	/* here we leave the size parameter for compatibility with the old hashtable implementation,
	   but this parameter IS NOT USED in this one. We use the static variable nbchunks_bitarray insead. */
    id_hash_table_t *new_table = (id_hash_table_t*) malloc(sizeof(id_hash_table_t));
    new_table->num_items = 0;

    /* Attempt to allocate and initialize to 0 the memory for the bitfield  */
    if ((new_table->bitarray = (bfield_t) calloc(nbchunks_bitarray, sizeof(unsigned long))) == NULL)
        return NULL;
    else
    	return new_table;
}

id_hash_table_t* complement_id_hashtbl(id_hash_table_t* h, int nbtaxa) {
	/* this creates a new hashtable and populates it with the complement of h */
	id_hash_table_t* c = create_id_hash_table(0);
	int retval;
	Taxon_id my_id;
	for (my_id = 0; my_id < nbtaxa; my_id++) {
		if (!lookup_id(h,my_id)) { retval = add_id(c, my_id); assert(retval == 0); }
	}
	return c;
}


int lookup_id(id_hash_table_t *hashtable, Taxon_id my_id)
{
    /* Returns whether the taxon is in the hashtable */ 
	if(my_id >= ntax) {
	  fprintf(stderr,"Error in %s: taxon ID %d is out of range. Aborting.\n", __FUNCTION__, my_id);
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}	       
	int chunk = my_id / chunksize;
	unsigned long *pointer = hashtable->bitarray + chunk; /* pointer to the long we want to access */
	int bit_index = my_id % chunksize;
	unsigned long mask = 1UL << bit_index; /* within a long, the lsb corresponds to the taxon with lowest TaxonID */
	return ((*pointer & mask) != 0);
}


int add_id(id_hash_table_t *hashtable, Taxon_id my_id)
{
    /* retcodes:
       0 -> no error, insertion has been performed successfully
       1 -> memory allocation failed: no space in memory (impossible in this implementation, though)
       2 -> the id we want to add already exists in the id_hashtable
    */
	int chunk = my_id / chunksize;
	unsigned long *pointer = hashtable->bitarray + chunk; /* pointer to the long we want to access */
	int bit_index = my_id % chunksize;
	unsigned long mask = (1UL << bit_index); /* within a long, the lsb corresponds to the taxon with lowest TaxonID */
	if (*pointer & mask) return 2;
	else {
		*pointer |= mask; /* sets to 1 the bit corresponding to the taxon. */
    		/* and update the total number of items in the hashtable */
		hashtable->num_items++;
		return 0;
	}
}

int delete_id(id_hash_table_t *hashtable, Taxon_id my_id)
{
    /* retcodes:
       0 -> no error, deletion has been performed successfully
       2 -> the id we are asked to delete was already set at 0 in the id_hashtable
    */
	int chunk = my_id / chunksize;
	unsigned long *pointer = hashtable->bitarray + chunk; /* pointer to the long we want to access */
	int bit_index = my_id % chunksize;
	unsigned long mask = (1UL << bit_index); /* within a long, the lsb corresponds to the taxon with lowest TaxonID */
	if (!(*pointer & mask)) return 2;
	else {
		*pointer &= ~mask; /* sets to 0 the bit corresponding to the taxon. */
    		/* and update the total number of items in the hashtable */
		hashtable->num_items--;
		return 0;
	}
}


void clear_id_hashtable(id_hash_table_t *hashtable) { /* clears completely the hashtable (no taxa) */
	int chunk;
	for (chunk = 0; chunk < nbchunks_bitarray; chunk++) hashtable->bitarray[chunk] = 0UL;
	hashtable->num_items = 0;
}


void fill_id_hashtable(id_hash_table_t *hashtable, int nb_taxa) { /* sets all bits to 1 in the whole hashtable (all taxa) */
	int chunk;
	unsigned long full_one = ~(0UL);
	for (chunk = 0; chunk < nbchunks_bitarray; chunk++) hashtable->bitarray[chunk] = full_one;
	/* the last bits of the last chunk are MEANINGLESS when chunksize is not a divisor of nb_taxa. */
	hashtable->num_items = nb_taxa;
}

void complement_id_hashtable(id_hash_table_t *destination, const id_hash_table_t *source, int nb_taxa) {
	/* transforms destination into the complement of source */
	int chunk;
	for (chunk = 0; chunk < nbchunks_bitarray; chunk++) destination->bitarray[chunk] = ~(source->bitarray[chunk]);
	destination->num_items = nb_taxa - source->num_items;
}

unsigned int bitCount (unsigned long value) {
    unsigned int count = 0;
    while (value) {           // until all bits are zero
        if (value & 0x1)     // check LSB
            count++;
        value >>= 1;              // shift bits, deleting LSB
    }
    return count;
}

void update_id_hashtable(id_hash_table_t *source, id_hash_table_t *destination) {
	/* copies all the items from source into destination. Doesn't erase anything anywhere.
	   Doesn't produce duplicate entries in the destination. */
	int chunk;
	unsigned int added;

	for (chunk = 0; chunk < nbchunks_bitarray; chunk++) {
		/* we first need to know how many new taxa we are going to add in destination */
		added = bitCount(source->bitarray[chunk] & ~destination->bitarray[chunk]); /* 1 in source AND O in dest */
		if (added) {
			/* copy all items from source->bitarray[chunk] into destination */
			destination->bitarray[chunk] = (destination->bitarray[chunk] | source->bitarray[chunk]);
			destination->num_items += added;
		} /* end if added */
	} /* end of the for loop */
} /* end update_id_hashtable */


int equal_id_hashtables(id_hash_table_t *tbl1, id_hash_table_t *tbl2) {
	/* this function compares the contents of the id_hashtables and returns a non-zero when tables are identical,
	   0 otherwise */
	if(tbl1 == NULL) return (tbl2 == NULL);
	if(tbl2 == NULL) return 0; /* because tbl1 not null */
	if(tbl1->num_items != tbl2->num_items) return 0; /* tables cannot be identical if they don't have the
							    same number of stored elements */
	int chunk;
	/* we simply test the equality of the successive longs */
	for (chunk = 0; chunk < nbchunks_bitarray; chunk++) {
		if (tbl1->bitarray[chunk] != tbl2->bitarray[chunk]) return 0;
	}
	/* here all the ids in tbl1 have been found also in tbl2, and the two tables have same size: */
	return 1;

} /* end equal_id_hashtables */


int complement_id_hashtables(id_hash_table_t *tbl1, id_hash_table_t *tbl2,int nb_taxa){
	/* this function compares the contents of the id_hashtables and returns a non-zero when tables are complement,
	   0 otherwise */
  if(tbl1 == NULL) return (tbl2 == NULL);
  if(tbl2 == NULL) return 0; /* because tbl1 not null */
  
  int chunk;
  /* we simply test the equality of the successive longs ==> Does not work for the last chunk */
  /* If the last long is < nbtaxa : the direct complement does not work!
     Example: 
        n taxa = 5
        chunk1 = 00000000 00000000 00000000 00011010
        chunk2 = 00000000 00000000 00000000 00000101
   ==> ~chunk2 = 11111111 11111111 11111111 11111010
        It does not work directly, we must put a mask depending on (nb_taxa%chunksize) 
	for the last chunk
         chunk1 & mask = 00000000 00000000 00000000 00011010
        ~chunk2 & mask = 00000000 00000000 00000000 00011010
	==> OK
	The mask is (((unsigned long)1 << (nb_taxa%chunksize)) - 1);
   */
  for (chunk = 0; chunk < nbchunks_bitarray; chunk++) {
    /* Initialize Mask with 1111....11*/
    unsigned long mask = -1;
    if(nb_taxa<(chunk+1)*chunksize){
      mask = (((unsigned long)1 << (nb_taxa%chunksize)) - 1);
    }
    if ((tbl1->bitarray[chunk]&mask) != ((~(tbl2->bitarray[chunk]))&mask)) return 0;
  }
  /* here all the ids in tbl1 have been found also in tbl2, and the two tables have same size: */
  return 1;
} /* end equal_id_hashtables */


int equal_or_complement_id_hashtables(id_hash_table_t *tbl1, id_hash_table_t *tbl2, int total) {
  return(complement_id_hashtables(tbl1,tbl2,total) ||
	 equal_id_hashtables(tbl1,tbl2));
} /* end equal_or_complement_id_hashtables */


id_hash_table_t* suffle_hash_table(id_hash_table_t *hashtable, int total){
  id_hash_table_t * output = create_id_hash_table(total);
  Taxon_id* taxid_array = malloc(total*sizeof(Taxon_id));
  Taxon_id i = 0;
  for(i=0;i<total;i++){
    taxid_array[i] = i;
  }
  shuffle(taxid_array, total, sizeof(Taxon_id));

  for(i=0;i<total;i++){
    if(lookup_id(hashtable, i)){
      add_id(output, taxid_array[i]);
    }
  }
  free(taxid_array);
  return(output);
}


void free_id_hashtable(id_hash_table_t *hashtable)
{
    if (hashtable==NULL) return;
    /* Free all the longs
     */
    free(hashtable->bitarray);
    free(hashtable);
}



void print_id_hashtable(FILE* stream, id_hash_table_t *hashtable, int nbtaxa) {
	int i, chunk;
	unsigned long mylong, base = 0, mask = 1, true_index;
	char c;
   	for (chunk = 0; chunk < nbchunks_bitarray; chunk++) {
		mylong = hashtable->bitarray[chunk];
		for (i = 0; i < chunksize; i++) { /* for all the bits in the unsigned long, starting with the LSB */
			true_index = base + i;
			if (true_index == nbtaxa) break; /* end of the last loop */
			if (true_index % 8 == 0 && !(chunk==0 && i == 0)) fputc(' ', stream); /* write blocks of 8 chars for legibility */
			if ((mylong & mask) == 1) c= '1' ; else c = '0';
			fputc(c, stream);
			mylong >>= 1;
		} /* end for on all the bits of the long */
		base += chunksize; /* so that in every loop, base is equal to chunk * chunksize */
	} /* end for on all the chunks (unsigned longs) */
    fputc('\n', stream);
} /* end print_id_hashtable */

