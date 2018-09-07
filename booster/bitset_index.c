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

#include "bitset_index.h"

bitset_hashmap* new_bitset_hashmap(int size, float loadfactor) {
  int i;
  bitset_hashmap* bh = malloc(sizeof(bitset_hashmap));
  bh->capacity = size;
  bh->loadfactor = loadfactor;
  bh->total = 0;
  bh->map_array = malloc(size*sizeof(bitset_bucket));
  for(i=0;i<size;i++){
    bh->map_array[i]=NULL;
  }
  return bh;
}

void free_bitset_hashmap(bitset_hashmap *hm){
  bitset_hash_map_free_map_array(hm->map_array, hm->capacity);
  free(hm);
}

void bitset_hash_map_free_map_array(bitset_bucket **map_array, int total){
  int i;
  for(i=0;i<total;i++){
    if(map_array[i]!=NULL){
      bitset_hash_map_free_buckets(map_array[i]->values, map_array[i]->size);
      free(map_array[i]);
    }
  }
  free(map_array);
}

void bitset_hash_map_free_buckets(bitset_keyvalue ** values, int total){
  int i;
  for(i=0;i<total;i++){
    free(values[i]);
  }
  free(values);
}

// returns the index in the hash map, given a hashcode
int bitset_hashmap_indexfor(int hashcode, int capacity) {
  return hashcode & (capacity - 1);
}

// Returns the count for the given Edge
// If the edge is not present, returns -1
// If the edge is present, returns the value
int bitset_hashmap_value(bitset_hashmap *hm, id_hash_table_t *bitset, int nb_taxa) {
  int index = bitset_hashmap_indexfor(bitset_hashcode(bitset,nb_taxa), hm->capacity);
  int k;
  if(hm->map_array[index] != NULL){
      for (k=0;k<hm->map_array[index]->size;k++){
	if(bitset_hashEquals(hm->map_array[index]->values[k]->key,bitset,nb_taxa)) {
	  return hm->map_array[index]->values[k]->value;
	}
      }
    }
  return -1;
}

void bitset_hashmap_putvalue(bitset_hashmap *hm, id_hash_table_t *bitset, int nb_taxa, int value) {
  int index = bitset_hashmap_indexfor(bitset_hashcode(bitset,nb_taxa), hm->capacity);
  int k;
  if(hm->map_array[index] == NULL) {
    hm->map_array[index] = malloc(sizeof(bitset_bucket));
    hm->map_array[index]->size=1;
    hm->map_array[index]->capacity=3;
    hm->map_array[index]->values=malloc(3*sizeof(bitset_keyvalue*));
    hm->map_array[index]->values[0] = malloc(sizeof(bitset_keyvalue));
    hm->map_array[index]->values[0]->key = bitset;
    hm->map_array[index]->values[0]->value = value;
    hm->total++;
  } else {
    for (k=0;k<hm->map_array[index]->size;k++){
      if(bitset_hashEquals(hm->map_array[index]->values[k]->key,bitset,nb_taxa)) {
	hm->map_array[index]->values[k]->value = value;
	return;
      }
    }
    if(hm->map_array[index]->size>=hm->map_array[index]->capacity){
      hm->map_array[index]->values = realloc(hm->map_array[index]->values,hm->map_array[index]->capacity*2*sizeof(bitset_keyvalue*));
      hm->map_array[index]->capacity *= 2;
    }
    hm->map_array[index]->values[hm->map_array[index]->size] = malloc(sizeof(bitset_keyvalue));
    hm->map_array[index]->values[hm->map_array[index]->size]->key = bitset;
    hm->map_array[index]->values[hm->map_array[index]->size]->value = value;
    hm->map_array[index]->size++;
    hm->total++;
  }
}

// Computes a hash code for the bitset associated to an edge
int bitset_hashcode(id_hash_table_t *hashtable, int nb_taxa){
  int hashCodeSet  = 1;
  int hashCodeUnset  = 1;
  int hashCodeAll = 1;
  int nbset = 0;
  int nbunset = 0;
  int bit;
  for (bit = 0; bit < nb_taxa; bit++) {
    if (lookup_id(hashtable, bit)){
      hashCodeSet = 31*hashCodeSet + bit;
      nbset++;
    } else {
      hashCodeUnset = 31*hashCodeUnset + bit;
      nbunset++;
    }
    hashCodeAll = 31*hashCodeAll + bit;
  }
  // If the number of species on the left is the same
  // than the number of species on the right
  // We return the hashcode of the all species
  // Otherwise, we return the hashcode for the minimum
  // between left and right
  // Allows an edge to be kind of "unique"
  if(nbset == nbunset){
    return hashCodeAll;
  } else if(nbset < nbunset){
    return hashCodeSet;
  }
  return hashCodeUnset;
}

// HashCode for an edge bitset.
// Used for insertion in an EdgeMap
int bitset_hashEquals(id_hash_table_t *tbl1, id_hash_table_t *tbl2, int nb_taxa) {
  return equal_or_complement_id_hashtables(tbl1, tbl2, nb_taxa);
}


// Reconstructs the HashMap if the capacity is almost attained (loadfactor)
void bitset_hashmap_rehash(bitset_hashmap *hm, int nb_taxa) {
  // We rehash everything with a new capacity
  if (((float)hm->total) >= ((float)hm->capacity) * hm->loadfactor) {
    int newcapacity = hm->capacity * 2;
    int i,l,k;
    bitset_bucket **new_map_array = malloc(newcapacity*sizeof(bitset_bucket*));
    for(i=0;i<newcapacity;i++){
      new_map_array[i]=NULL;
    }

    for(k=0;k<hm->capacity;k++){
      if (hm->map_array[k] != NULL) {
	for(l=0;l<hm->map_array[k]->size;l++){
	  int index = bitset_hashmap_indexfor(bitset_hashcode(hm->map_array[k]->values[l]->key,nb_taxa), newcapacity);
	  if (new_map_array[index] == NULL) {
	      new_map_array[index] = malloc(sizeof(bitset_bucket));
	      new_map_array[index]->size=1;
	      new_map_array[index]->capacity=3;
	      new_map_array[index]->values=malloc(3*sizeof(bitset_keyvalue*));
	      new_map_array[index]->values[0] = malloc(sizeof(bitset_keyvalue));
	      new_map_array[index]->values[0]->key = hm->map_array[k]->values[l]->key;
	      new_map_array[index]->values[0]->value = hm->map_array[k]->values[l]->value;
	    } else {
	    if(new_map_array[index]->size>=new_map_array[index]->capacity){
	      new_map_array[index]->values = realloc(new_map_array[index]->values,new_map_array[index]->capacity*2*sizeof(bitset_keyvalue*));
	      new_map_array[index]->capacity *= 2;
	    }
	    new_map_array[index]->values[new_map_array[index]->size] = malloc(sizeof(bitset_keyvalue));
	    new_map_array[index]->values[new_map_array[index]->size]->key = hm->map_array[k]->values[l]->key;
	    new_map_array[index]->values[new_map_array[index]->size]->value = hm->map_array[k]->values[l]->value;
	    new_map_array[index]->size++;
    	  }
	}
      }
    }
    hm->capacity = newcapacity;
    bitset_hash_map_free_map_array(hm->map_array,hm->total);
    hm->map_array = new_map_array;
  }
}
