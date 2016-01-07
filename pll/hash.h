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
 * @file hash.h
 */
#ifndef __pll_HASH__
#define __pll_HASH__

struct pllHashItem
{
  void * data;
  char * str;
  struct pllHashItem * next;
};

struct pllHashTable
{
  unsigned int size;
  struct pllHashItem ** Items;
};

unsigned int pllHashString (const char * s, unsigned int size);
int pllHashAdd  (struct pllHashTable * hTable, const char * s, void * item);
struct pllHashTable * pllHashInit (unsigned int n);
int pllHashSearch (struct pllHashTable * hTable, char * s, void ** item);
void pllHashDestroy (struct pllHashTable ** hTable, int);
#endif
