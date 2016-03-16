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
 * @file hash.c
 */
#include <stdio.h>
#include <string.h>
#include "pll.h"
#include "mem_alloc.h"

static const unsigned int initTable[] = 
  {
    53,         97,         193,       389,       769,    
    1543,       3079,       6151,      12289,     24593, 
    49157,      98317,      196613,    393241,    786433, 
    1572869,    3145739,    6291469,   12582917,  25165843, 
    50331653,   100663319,  201326611, 402653189, 805306457, 
    1610612741, 3221225473, 4294967291
  };
       
/** @brief Generate the hash value for a string 

    Generates the hash value of a string \a s.

    @param s     The string to compute the hash for
    @param size  Size of the hash table
    @return      String hash \a s, i.e. index in hash table
*/
unsigned int pllHashString (const char * s, unsigned int size)
{
  unsigned int hash = 0;

  for (; *s; ++s) hash = (hash << 5) - hash + (unsigned int )*s;

  return (hash % size);
}

/** @brief Add a string and its data to a hashtable
    
    Add an \a item and possibly a string \a s to hashtable \a hTable at position
    \a hash, where \a hash must be a value between 0 and \a hTable->size - 1. If
    string \a s is given and another record with the same computed hash and the
    same associated string exists in the hash table, then the new record will \b not be added and the
    value \b PLL_FALSE is returned. Otherwise, the new item is added at the
    beginning of the corresponding linked list and the value \b PLL_TRUE is
    returned.

    @param hTable Hashtable
    @param hash   Position where to store in hash table
    @param s      String
    @param item   Data associated with \a s
    @return       Returns \b PLL_TRUE if added with success, otherwise \b PLL_FALSE
*/
int pllHashAdd  (pllHashTable * hTable, unsigned int hash, const char * s, void * item)
{
  pllHashItem * hItem;

  hItem = hTable->Items[hash];

  /* If a string was given, check whether the record already exists */
  if (s)
   {
     for (; hItem; hItem = hItem->next)
      {
        if (hItem->str && !strcmp (s, hItem->str)) return (PLL_FALSE);
      }
   }

  hItem = (pllHashItem *) rax_malloc (sizeof (pllHashItem));

  /* store the string together with the element if given */
  if (s)
   {
     hItem->str = (char *) rax_malloc ((strlen(s) + 1) * sizeof (char));
     strcpy (hItem->str, s);
   }
  else
   hItem->str = NULL;

  hItem->data = item;

  hItem->next = hTable->Items[hash];
  hTable->Items[hash] = hItem;
  hTable->entries += 1;

  return (PLL_TRUE);
}

       
/** @brief Initialize hash table
    
    Create a hash table of size at least \a n. The size of the hash table will
    be the first prime number higher or equal to \a n.

    @param n  Minimum size of hash table
    @return   In case of success, returns a pointer to the created hash table, otherwise returns \b NULL
*/
pllHashTable * pllHashInit (unsigned int n)
{ 
  pllHashTable * hTable;
  unsigned int i;
  unsigned int primeTableLength;
       
  hTable = (pllHashTable *) rax_malloc (sizeof (pllHashTable));
  if (!hTable) return (NULL);
  
  primeTableLength = sizeof (initTable) / sizeof(initTable[0]);

  i = 0;
 
  while (initTable[i] < n && i < primeTableLength) ++ i;
 
  n = initTable[i];  
 
  hTable->Items = (pllHashItem **) rax_calloc (n, sizeof (pllHashItem *));
  if (!hTable->Items)
   {
     rax_free (hTable);
     return (NULL);
   }
  hTable->size    = n;
  hTable->entries = 0;
 
  return (hTable);
}

/** @brief Retrieve the data stored in hash table for a given string

    Retrieve the data stored in hash table \a hTable under a given string \a s.
    In case the string is found in the hash table, the associated data are
    stored in \a item and the function returns \b PLL_TRUE. In the opposite
    case, or if \a s is given as \b NULL then \b PLL_FALSE is returned.

    @param hTable   Hash table to be searched
    @param s        String to look for
    @param item     Where to store the retrieved data
    @return         Returns \b PLL_TRUE if the string was found, otherwise \b PLL_FALSE
*/
int pllHashSearch (pllHashTable * hTable, char * s, void ** item)
{
  unsigned int pos;
  pllHashItem * hItem;

  if (!s) return (PLL_FALSE);

  pos   = pllHashString (s, hTable->size);
  hItem = hTable->Items[pos];

  for (; hItem; hItem = hItem->next)
   {
     if (hItem->str && !strcmp (s, hItem->str))
      {
        *item = hItem->data;
        return (PLL_TRUE);
      }
   }

  return (PLL_FALSE);
}

/** @brief Deallocate a hash table

    Deallocates the hash table. A callback function may be specified as \a
    cbDealloc which will be executed upon all \a data elements of the hash
    table, for deallocating custom data. If no deallocation is required for the
    custom data, then \a cbDealloc must be set to \b NULL. The strings
    associated with each hash element are deallocated.

    @param hTable    Hash table to be deallocated
    @pram  cbDealloc Callback function to perform deallocation of each data element of the hash table
    @notes
      Deallocates the structure for the hash table. Note that the 
      data associated with the indexed strings are not deallocated.
*/
void pllHashDestroy (pllHashTable ** hTable, void (*cbDealloc)(void *))
{
  unsigned int i;
  pllHashItem * hItem;
  pllHashItem * tmp;

  for (i = 0; i < (*hTable)->size; ++ i)
  {
    hItem = (*hTable)->Items[i];
    while (hItem)
     {
       tmp   = hItem;
       hItem = hItem->next;
       if (tmp->str)  rax_free (tmp->str);
       if (cbDealloc) cbDealloc (tmp->data);
       rax_free (tmp);
     }
  }
  rax_free ((*hTable)->Items);
  rax_free (*hTable);
  *hTable = NULL;
}
