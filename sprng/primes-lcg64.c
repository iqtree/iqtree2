#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "primes-lcg64.h"
#include "primelist-lcg64.h"

#define YES 1
#define NO  0
#define NPRIMES 10000

int primes[NPRIMES];

#ifdef __STDC__
int init_prime(void)
#else
int init_prime()
#endif
{
  int i, j, obtained = 0, isprime;
  
  for(i=3; i < MINPRIME; i += 2)
  {
    isprime = YES;
    
    for(j=0; j < obtained; j++)
      if(i%primes[j] == 0)
      {
	isprime = NO;
	break;
      }
    else if(primes[j]*primes[j] > i)
      break;

    if(isprime == YES)
    {
      primes[obtained] = i;
      obtained++;
    }
  }
  
  return obtained;
}




#ifdef __STDC__
int getprime(int need, unsigned int *prime_array, int offset)
#else
int getprime(need, prime_array, offset)
int need, offset;
unsigned int *prime_array;
#endif
{
  static int initiallized = NO, num_prime;
  unsigned int largest;
  int i, isprime, index, obtained = 0;
  
  if(need <= 0)
  {
    fprintf(stderr,"WARNING: Number of primes needed = %d < 1; None returned\n"
	    , need);
    return 0;
  }
  
  if(offset < 0)
  {
    fprintf(stderr,"WARNING: Offset of prime = %d < 1; None returned\n"
	    , offset);
    return 0;
  }
  

  if(offset+need-1<PRIMELISTSIZE1) 
  {
    memcpy(prime_array,prime_list+offset,need*sizeof(unsigned int));
    return need;
  }

  if(!initiallized)
  {
    num_prime = init_prime();

    
    largest = MAXPRIME;
    initiallized = YES;
  }
  
  /* int offset  <=>  unsigned MAXPRIMEOFFSET ??? (HAS) */
  if(offset > MAXPRIMEOFFSET)
  {
    fprintf(stderr,"WARNING: generator has branched maximum number of times;\nindependence of generators no longer guaranteed");
    offset = offset % MAXPRIMEOFFSET;
  }
  
  if(offset < PRIMELISTSIZE1)	/* search table for previous prime */
  {
    largest = prime_list[offset] + 2;
    offset = 0;
  }
  else
  {
    index = (unsigned int) ((offset-PRIMELISTSIZE1+1)/STEP) + PRIMELISTSIZE1 -  1;
    largest = prime_list[index] + 2;
    offset -= (index-PRIMELISTSIZE1+1)*STEP + PRIMELISTSIZE1 - 1;
  }
  
  
  while(need > obtained && largest > MINPRIME)
  {
    isprime = YES;
    largest -= 2;
    for(i=0; i<num_prime; i++)
      if(largest%primes[i] == 0)
      {
	isprime = NO;
	break;
      }
    
    if(isprime == YES && offset > 0)
      offset--;
    else if(isprime == YES)
      prime_array[obtained++] = largest;
  }
  
  if(need > obtained)
    fprintf(stderr,"ERROR: Insufficient number of primes: needed %d, obtained %d\n", need, obtained);
  
  return obtained;
}


#if 0
main()
{
  unsigned int newprimes[1500], np, i;
  
  np = getprime(2,newprimes,0);
  np = getprime(2,newprimes+2,9);
  np = getprime(2,newprimes+4,12);
  
   for(i=0; i<6; i++)
    printf("%u. %u \n", i, newprimes[i]);
  
  /*while(np--)
    printf("New primes: %u\n", newprimes[np]);

  np = getprime(5,newprimes);
  
  printf("%d new primes obtained ...\n", np);
  
  while(np--)
    printf("New primes: %u\n", newprimes[np]);*/
}
#endif
