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

#include "sort.h"

void sort_double(double*tab, int size){
  qsort(tab, size, sizeof(double), comp_double);
}

/* void sort_indexes_double(int * indexes, int size, double * values){ */
/*   #ifdef __APPLE__ */
/*   qsort_r(indexes, size, sizeof(int), values, comp_indexes_apple); */
/*   #else */
/*   qsort_r(indexes, size, sizeof(int), comp_indexes, values); */
/*   #endif */
/* } */

int comp_double(const void * elem1, const void * elem2){
  double f = *((double*)elem1);
  double s = *((double*)elem2);
  if (f > s) return  1;
  if (f < s) return -1;
  return 0;
}

int comp_indexes(const void * elem1, const void * elem2, void * other_array){
  int i1 = *((int*)elem1);
  int i2 = *((int*)elem2);
  
  double * other = (double*)other_array;

  double val1 = other[i1];
  double val2 = other[i2];

  if (val1 > val2) return  1;
  if (val1 < val2) return -1;
  return 0;
}

int comp_indexes_apple(void * other_array, const void * elem1, const void * elem2){
  return(comp_indexes(elem1,elem2,other_array));
}
