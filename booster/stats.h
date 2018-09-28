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

#ifndef _STAT_H
#define _STAT_H

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "prng.h"
#include "io.h"

#define S_PI 3.14159265358979323846264338327950288

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/************************************************/
/*                BASIC FUNCTIONS               */
/************************************************/

int min_int(int a, int b);
int max_int(int a, int b);

int max_int_vec(int* myvec, int length);
short unsigned max_short_unsigned_vec(short unsigned* myvec, int length);

double min_double(double a, double b);
double max_double(double a, double b);

void print_int_vec(FILE* out, int* myvec, int length);
void print_double_vec(FILE* out, double* myvec, int length);

double mean_int_vec(int* myvec, int length);
double mean_double_vec(double* myvec, int length);

int median_int_vec(int* myvec, int length);
double median_double_vec(double* myvec, int length);

void summary_double_vec(double* myvec, int length, double* result);
void summary_double_vec_nocopy(double* myvec, int length, double* result);

int sum_vec_of_ints(int* table, int size);
int sum_vec_of_ints_but_one(int* table, int size, int index_to_ignore);

int swap_ints(int* a, int* b);
int swap_doubles(double* a, double* b);

void merge_sorted_int_vecs(int* myvec, int length1, int length2);
void divide_and_conquer_int_vec(int* vec, int length);

void merge_sorted_double_vecs(double* myvec, int length1, int length2);
void divide_and_conquer_double_vec(double* vec, int length);

/************************************************/
/*               STAT FUNCTIONS                 */
/************************************************/
double unif();
double exponentiel(double lambda);
double gauss();
double normal(double mu, double sig);
int    proba(double p);
int    binomial(double p, int nb);

/* Sample num ints from the data (of length size) 
   if !replace then without replacement
*/
int* sample(int* data, int size, int num, int replace);
/* Shuffles the array */
#define BYTE(X) ((unsigned char *)(X)) 
void shuffle(void *obj, size_t nmemb, size_t size);
/* Samples num values from the ungrouped version of the data array:
   Example: data array: 
   data[0]=3; data[1]=0; data[2]=4
   It will return a sample (of size num ) from :
   0,0,0,2,2,2,2
   num must be <= sum(data) : otherwize returns 0 filled array
   The output is grouped by indice , i.e:
   output[0]=2; output[1]=0; output[2]=3
   AND NOT:
   0,0,2,2,2
   So the output has the same size than data , i.e : length
*/
int* sample_from_counts(int* data, int length, int num, int replace);

/* rand in [0,max[ */
int rand_to(int max);

/* ecart type */
double sigma(double * values, int nb_values);
double sum(double * array, int size);
double qnorm(double x, double mean, double sd);
double pnorm(double x);

/* Computes the factorial of n */
double log_fact(int n);
/* Computes the log of factorial of n using rmnj approximation */
double factorial_log_rmnj(int n);

#endif
