/**
    Computing with GPUs template
    Copyright (C) 2012 Nikolaos Alachiotis, Fernando Izquierdo, and 
    Solon P. Pissis.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

#define MAX_NUM 1000


int main ( int argc, char * argv [] )
 {
   unsigned int i;
   unsigned int n = MAX_NUM;
   int * a;
   int * b;
   int * c1;
   int * c2;

   a = calloc ( n, sizeof( int ) );
   b = calloc ( n, sizeof( int ) );
   c1 = calloc ( n, sizeof( int ) );
   c2 = calloc ( n, sizeof( int ) );

   for ( i = 0; i < n; i ++ )
	a[i] = b[i] = i;

add_two_vectors_cpu ( n, a, b, c1 );

///////////////////////////////////////////////////////// GPU CALLS /////////////////////////////////////////////////////////////////////////
#ifdef _USE_GPU

add_two_vectors_gpu ( n, a, b, c2 );

#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   for ( i = 0; i < n; i ++ )
	if ( c1[i] != c2[i] )
          { 
		printf("\nError: c1[%d]=%d != c2[%d]=%d", i, c1[i], i, c2[i] );
                break;
          }
   printf(" The programme finished succesfully.\n" );
     
   return ( 0 );
 }
