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

#include "stats.h"

/************************************************/
/*                BASIC FUNCTIONS               */
/************************************************/

/* this file contains basic operations on numerical arrays (min, max, mean, sorting, median, debug printing, etc) */
int min_int(int a, int b) {
	return (a<b ? a : b);
}

int max_int(int a, int b) {
	return (a>b ? a : b);
}

int max_int_vec(int* myvec, int length) {
	if (length==0) return -1;
	int i, maximum = myvec[0];
	for(i=1;i<length;i++) if (maximum < myvec[i]) maximum = myvec[i];
	return maximum;
}

short unsigned max_short_unsigned_vec(short unsigned* myvec, int length) {
	if (length==0) return -1;
	int i;
	short unsigned maximum = myvec[0];
	for(i=1;i<length;i++) if (maximum < myvec[i]) maximum = myvec[i];
	return maximum;
}

double min_double(double a, double b) {
	return (a<b ? a : b);
}

double max_double(double a, double b) {
	return (a>b ? a : b);
}


void print_int_vec(FILE* out, int* myvec, int length) {
	int i;
	for(i=0;i<length-1;i++) fprintf(out,"%d ", myvec[i]);
	fprintf(out,"%d\n", myvec[length-1]);
}


void print_double_vec(FILE* out, double* myvec, int length) {
	int i;
	for(i=0;i<length-1;i++) fprintf(out,"%.4g ", myvec[i]);
	fprintf(out,"%.4g\n", myvec[length-1]);
}

double mean_int_vec(int* myvec, int length) {
	int i, accu=0;
	for (i=0;i<length;i++) accu += myvec[i];
	return ((double) accu) / length;
}


double mean_double_vec(double* myvec, int length) {
	int i;
	double accu=0.0;
	for (i=0;i<length;i++) accu += myvec[i];
	return accu / length;
}

int median_int_vec(int* myvec, int length) {
	/* we don't want to modify the original vector, so work on a copy that is going
	   to be sorted: */
	int i, mycopy[length];
	for(i=0;i<length;i++) mycopy[i] = myvec[i];
	divide_and_conquer_int_vec(mycopy, length);
	return mycopy[(int)(floor(length/2))];
}


double median_double_vec(double* myvec, int length) {
	/* we don't want to modify the original vector, so work on a copy that is going
	   to be sorted: */
	int i;
	double mycopy[length];
	for(i=0;i<length;i++) mycopy[i] = myvec[i];
	divide_and_conquer_double_vec(mycopy, length);
	return mycopy[(int)(floor(length/2))];
}

void summary_double_vec(double* myvec, int length, double* result) {
	/* the result vector HAS TO BE ALLOCATED BEFOREHAND, size at least 6.
	   Same as the result function in R:
	   0) minimum
	   1) 1st quartile
	   2) median
	   3) mean
	   4) 3rd quartile
	   5) maximum */

	int i;
	double mycopy[length];
	for(i=0;i<length;i++) mycopy[i] = myvec[i];
	divide_and_conquer_double_vec(mycopy, length);
	result[0] = mycopy[0];				/* min */
	result[1] = mycopy[(int)(floor(length/4))];	/* 1st quart. */
	result[2] = mycopy[(int)(floor(length/2))];	/* median */
	result[3] = mean_double_vec(mycopy, length);	/* mean */ 
	result[4] = mycopy[(int)(floor(3*length/4))];	/* 3rd quart. */
	result[5] = mycopy[length-1];			/* max */
} /* end summary_double_vec */


void summary_double_vec_nocopy(double* myvec, int length, double* result) {
	/* the result vector HAS TO BE ALLOCATED BEFOREHAND, size at least 6.
	   Same as the result function in R:
	   0) minimum
	   1) 1st quartile
	   2) median
	   3) mean
	   4) 3rd quartile
	   5) maximum */
	/* nocopy: same as function above but modifies the array in place */

	divide_and_conquer_double_vec(myvec, length);
	result[0] = myvec[0];				/* min */
	result[1] = myvec[(int)(floor(length/4))];	/* 1st quart. */
	result[2] = myvec[(int)(floor(length/2))];	/* median */
	result[3] = mean_double_vec(myvec, length);	/* mean */ 
	result[4] = myvec[(int)(floor(3*length/4))];	/* 3rd quart. */
	result[5] = myvec[length-1];			/* max */
} /* end summary_double_vec_nocopy */


int sum_vec_of_ints(int* table, int size) {
	/* simply gives the sum of the vector */
	int i, accu = 0;
	for (i=0; i< size; i++) accu += table[i];
	return accu;
} /* end sum_vec_of_ints */



int sum_vec_of_ints_but_one(int* table, int size, int index_to_ignore) {
	/* simply gives the sum of the vector */
	int i, accu = 0;
	for (i=0; i< size; i++) if(i != index_to_ignore) accu += table[i];
	return accu;
} /* end sum_vec_of_ints_but_one */

int swap_ints(int* a, int* b) {
	if (a == NULL || b == NULL) return 1;
	int temp = *b;
	*b = *a;
	*a = temp;
	return 0;
}

int swap_doubles(double* a, double* b) {
	if (a == NULL || b == NULL) return 1;
	double temp = *b;
	*b = *a;
	*a = temp;
	return 0;
}


void merge_sorted_int_vecs(int* myvec, int length1, int length2) {
	/* this function assumes that we have myvec[0..(length1-1)]
	   and myvec[length1..(length1+length2-1)] that are two sorted vectors.
	   It merges the two in place, reusing the initial space. */
	int i, index1=0, index2=0, index_res=0, total_length = length1 + length2;
	int temp[total_length];
	int* vec1 = myvec, *vec2 = myvec+length1; /* pointer arithmetic */
	/* index1 and index2 indicate the next elements of the two subvectors to be processed */
	while(index1 < length1 && index2 < length2) {
		/* there are still elements to treat in both vectors */
		if(vec1[index1] <= vec2[index2]) temp[index_res++] = vec1[index1++];
		else temp[index_res++] = vec2[index2++];
	}
	/* now at least one of the input subvecs is fully processed, remains the other: */
	if (index1 < length1) for (i = index1; i < length1; i++) temp[index_res++] = vec1[i];
	else for (i = index2; i < length2; i++) temp[index_res++] = vec2[i];
	/* sanity check */
	if (index_res != total_length) {
	  fprintf(stderr,"fatal error : input lengths do not sum up to output length. Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}
	/* now we copy the result back into the original vector, to do the thing in place */
	for(i=0;i<total_length;i++) myvec[i] = temp[i];
} /* end of merge_sorted_int_vecs */



void divide_and_conquer_int_vec(int* vec, int length) {
	/* this function works "in place" and does not allocate extra memory.
	   The allocation is done during the merge step, but through a local variable there. */
	if (length < 2) return; /* nothing to do here */
	if (length == 2) {
		if (vec[0] > vec[1]) swap_ints(vec,vec+1); /* swapping with pointer arithmetic */
		return; /* we're done */
	} /* end if length == 2 */

	/* implicit else: here length > 2 */
	int breakpoint = (int) floor(length / 2);
	/* breakpoint is the number of values in the first half */
	int length1 = breakpoint, length2 = length - breakpoint;

	divide_and_conquer_int_vec(vec, length1);
	divide_and_conquer_int_vec(vec+breakpoint, length2);
	merge_sorted_int_vecs(vec, length1, length2);
	return ;

} /* end divide_and_conquer_int_vec */


void merge_sorted_double_vecs(double* myvec, int length1, int length2) {
	/* this function assumes that we have myvec[0..(length1-1)]
	   and myvec[length1..(length1+length2-1)] that are two sorted vectors.
	   It merges the two in place, reusing the initial space. */
	int i, index1=0, index2=0, index_res=0, total_length = length1 + length2;
	double temp[total_length];
	double* vec1 = myvec, *vec2 = myvec+length1; /* pointer arithmetic */
	/* index1 and index2 indicate the next elements of the two subvectors to be processed */
	while(index1 < length1 && index2 < length2) {
		/* there are still elements to treat in both vectors */
		if(vec1[index1] <= vec2[index2]) temp[index_res++] = vec1[index1++];
		else temp[index_res++] = vec2[index2++];
	}
	/* now at least one of the input subvecs is fully processed, remains the other: */
	if (index1 < length1) for (i = index1; i < length1; i++) temp[index_res++] = vec1[i];
	else for (i = index2; i < length2; i++) temp[index_res++] = vec2[i];
	/* sanity check */
	if (index_res != total_length) {
	  fprintf(stderr,"fatal error : input lengths do not sum up to output length. Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}
	/* now we copy the result back into the original vector, to do the thing in place */
	for(i=0;i<total_length;i++) myvec[i] = temp[i];
} /* end of merge_sorted_double_vecs */



void divide_and_conquer_double_vec(double* vec, int length) {
	/* this function works "in place" and does not allocate extra memory.
	   The allocation is done during the merge step, but through a local variable there. */
	if (length < 2) return; /* nothing to do here */
	if (length == 2) {
		if (vec[0] > vec[1]) swap_doubles(vec,vec+1); /* swapping with pointer arithmetic */
		return; /* we're done */
	} /* end if length == 2 */

	/* implicit else: here length > 2 */
	int breakpoint = (int) floor(length / 2);
	/* breakpoint is the number of values in the first half */
	int length1 = breakpoint, length2 = length - breakpoint;

	divide_and_conquer_double_vec(vec, length1);
	divide_and_conquer_double_vec(vec+breakpoint, length2);
	merge_sorted_double_vecs(vec, length1, length2);
	return ;

} /* end divide_and_conquer_double_vec */


/************************************************/
/*               STAT FUNCTIONS                 */
/************************************************/


double unif(){
  double unif = 0.5;
  unif = (unif + prng_get_int())/ INT_MAX;
  return(unif);
}

double exponentiel(double lambda){
  double exponentiel = unif();
  exponentiel = -log(1 - exponentiel) / lambda;
  return(exponentiel);
}

double gauss(){
  double unif1 = unif();
  double unif2 = unif();
  double gauss = sqrt(-2*log(unif1))*sin(2 * S_PI * (unif2));
  return(gauss);
}

double normal(double mu, double sig){
  return(mu + (sig*gauss()));
}

int proba(double p){
  return(unif()<p);
}

int binomial(double p, int nb){
  int binom = 0;
  int i = 0;
  for(i = 0; i < nb; i++){
    binom+=unif() < p;
  }
  return(binom);
}

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
int* sample_from_counts(int* data, int length, int num, int replace){
  int total = 0;
  int * values;
  int * counts;
  int i,j;
  int current;
  int* sampled;

  counts = malloc( length * sizeof(int));
  for(i=0; i < length; i++){
    total += data[i];
    counts[i] = 0;
  }
  
  if( total < num ){
    return(counts);
  }

  values = malloc( total * sizeof(int));
  current=0;
  for(i=0;i<length;i++){
    for(j=0;j<data[i];j++){
      values[current] = i;
      current++;
    }
    counts[i]=0;
  }
  sampled = sample(values, total, num, replace);
  for(j=0;j<num;j++){
    counts[sampled[j]]++;
  }
  free(sampled);
  free(values);
  return(counts);
}

/* Sample num ints from the input of length length, with or without replacement */
int* sample(int* data, int length, int num, int replace){
  int * output = malloc(num * sizeof(int));
  int i=0;

  /* Without replacement */
  if(!replace){
    int * temp  =  malloc(length * sizeof(int));
    for(i=0; i < length; i++){
      temp[i] = data[i];
    }
    shuffle(temp,length,sizeof(int));
    for(i=0;i<num;i++){
      output[i] = temp[i];
    }
    free(temp);
  }else{
    /* With replacement */
    for(i=0;i<num;i++){
      output[i] = data[rand_to(length)];
    }
  } 
  return output;
}

/* Shuffles the data in the array of length size */
void shuffle(void *obj, size_t nmemb, size_t size){
  void *temp = malloc(size);
  size_t n = nmemb;
  while ( n > 1 ) {
    size_t k = rand_to(n--);
    memcpy(temp, BYTE(obj) + n*size, size);
    memcpy(BYTE(obj) + n*size, BYTE(obj) + k*size, size);
    memcpy(BYTE(obj) + k*size, temp, size);
  }
  free(temp);
} 

/* take a random int from [0,max[ */
int rand_to(int max){
  return(prng_get_int()%max);
}

double sigma(double * values, int nb_values){
  double mean = 0.0;
  double var = 0.0;
  int i;
  for(i = 0; i < nb_values; i++){
    mean += values[i];
  }

  for(i = 0; i < nb_values; i++){
    var += pow((values[i] - mean),2);
  }
  return(sqrt(var));
}

double sum(double * array, int size){
  int i;
  double sum = 0;
  for(i = 0; i < size; i++){
    sum += array[i];
  }
  return(sum);
}

/* Original C++ implementation found at http://www.wilmott.com/messageview.cfm?catid=10&threadid=38771 */
/* C# implementation found at http://weblogs.asp.net/esanchez/archive/2010/07/29/a-quick-and-dirty-implementation-of-excel-norminv-function-in-c.aspx*/
/*
 *     Compute the quantile function for the normal distribution.
 *
 *     For small to moderate probabilities, algorithm referenced
 *     below is used to obtain an initial approximation which is
 *     polished with a final Newton step.
 *
 *     For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *     Beasley, J. D. and S. G. Springer (1977).
 *     Algorithm AS 111: The percentage points of the normal distribution,
 *     Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */
/* Taken from https://gist.github.com/kmpm/1211922/ */
double qnorm(double p, double mu, double sigma){
    double q, r, val;

    if (p < 0 || p > 1){
      fprintf(stderr,"Warning: p is < 0 or > 1 : returning DBL_MIN\n");
      return NAN;
    }
    if (sigma < 0){
      fprintf(stderr,"Warning: sigma is < 0 : returning NaN\n");
      return NAN;
    }
    if (p == 0){
        return -INFINITY;
    }
    if (p == 1){
        return INFINITY;
    }

    if (sigma == 0){
        return mu;
    }
    q = p - 0.5;
    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
            Produces the normal deviate Z corresponding to a given lower
            tail area of P; Z is accurate to about 1 part in 10**16.
    */
    if (fabs(q) <= .425){/* 0.075 <= p <= 0.925 */
      r = .180625 - q * q;
      val =
	q * (((((((r * 2509.0809287301226727 +
		   33430.575583588128105) * r + 67265.770927008700853) * r +
		 45921.953931549871457) * r + 13731.693765509461125) * r +
	       1971.5909503065514427) * r + 133.14166789178437745) * r +
	     3.387132872796366608)
	/ (((((((r * 5226.495278852854561 +
		 28729.085735721942674) * r + 39307.89580009271061) * r +
	       21213.794301586595867) * r + 5394.1960214247511077) * r +
	     687.1870074920579083) * r + 42.313330701600911252) * r + 1);
    } else { /* closer than 0.075 from {0,1} boundary */
      /* r = min(p, 1-p) < 0.075 */
      if (q > 0)
	r = 1 - p;
      else
	r = p;
      r = sqrt(-log(r));
      /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
      if (r <= 5){ /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
	r += -1.6;
	val = (((((((r * 7.7454501427834140764e-4 +
		     .0227238449892691845833) * r + .24178072517745061177) *
		   r + 1.27045825245236838258) * r +
		  3.64784832476320460504) * r + 5.7694972214606914055) *
		r + 4.6303378461565452959) * r +
	       1.42343711074968357734)
	  / (((((((r *
		   1.05075007164441684324e-9 + 5.475938084995344946e-4) *
		  r + .0151986665636164571966) * r +
		 .14810397642748007459) * r + .68976733498510000455) *
	       r + 1.6763848301838038494) * r +
	      2.05319162663775882187) * r + 1);
      } else { /* very close to  0 or 1 */
	r += -5;
	val = (((((((r * 2.01033439929228813265e-7 +
		     2.71155556874348757815e-5) * r +
		    .0012426609473880784386) * r + .026532189526576123093) *
		  r + .29656057182850489123) * r +
		 1.7848265399172913358) * r + 5.4637849111641143699) *
	       r + 6.6579046435011037772)
	  / (((((((r *
		   2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
		  r + 1.8463183175100546818e-5) * r +
		 7.868691311456132591e-4) * r + .0148753612908506148525)
	       * r + .13692988092273580531) * r +
	      .59983220655588793769) * r + 1);
      }
      if (q < 0.0){
	val = -val;
      }
    }
    return mu + sigma * val;
}


/* From https://en.wikipedia.org/wiki/Normal_distribution */
double pnorm(double x){
  double value,sum,result;
  int i;
  sum = x;
  value=x;
  for(i=1;i<=100;i++){
    value=(value*x*x/(2*i+1));
    sum=sum+value;
  }
  result=0.5+(sum/sqrt(2*S_PI))*exp(-(x*x)/2);
  return(result);
}

double log_fact(int n){
  int i;
  double lf = (double) 0.0;
  for (i = 2; i <= n; i++){
    lf = lf + (double) log((double)i);
  }
  return lf;
}

double factorial_log_rmnj(int n){
  if (n==0) {
    return(0.0);
  } else if (n<=100) {
    return(log_fact(n));
  } else {
    double accu = 0.0;
    accu += (double) log((double)n*(1.0+4.0*n*(1.0+2.0*n)) + 1.0/30.0 - 11.0/(240.0*n))/6.0;
    accu += (double) log(S_PI)/ 2.0;
    accu -= (double) n;
    accu += (double) n * log(n);
    return( accu );
  }
}
