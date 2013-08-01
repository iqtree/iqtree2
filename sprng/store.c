#include <stdio.h>

#include "store.h"


int store_long(unsigned long l, int nbytes, unsigned char *c)
{
  int i;
  
  for(i=0; i<nbytes; i++)
    c[i] = (l>>(8*(nbytes-i-1)))&0xff;
  
  return nbytes;		/* return number of chars filled */
}


int store_longarray(unsigned long *l, int n, int nbytes, unsigned char *c)
{
  int i;
  
  for(i=0; i<n; i++)
    c += store_long(l[i],nbytes,c);

  return nbytes*n;
}


int load_long(unsigned char *c, int nbytes, unsigned long *l)
{
  int i;
  
  *l = 0;
  
  for(i=0; i<nbytes;i++)
    *l = (*l<<8) + (c[i]&0xff);
 
  return nbytes;
}


int load_longarray(unsigned char *c, int n, int nbytes, unsigned long *l)
{
  int i;
  
  for(i=0; i<n; i++)
    load_long(c+nbytes*i,nbytes,l+i);

  return nbytes*n;
}

#ifdef _LONG_LONG

int store_longlong(unsigned long long l, int nbytes, unsigned char *c)
{
  int i;
  
  for(i=0; i<nbytes; i++)
    c[i] = (l>>(8*(nbytes-i-1)))&0xff;
  
  return nbytes;		/* return number of chars filled */
}


int store_longlongarray(unsigned long long *l, int n, int nbytes, unsigned char *c)
{
  int i;
 
  for(i=0; i<n; i++)
    c += store_longlong(l[i],nbytes,c);

  return nbytes*n;
}


int load_longlong(unsigned char *c, int nbytes, unsigned long long *l)
{
  int i;
 
  *l = 0;
  
  for(i=0; i<nbytes;i++)
    *l = (*l<<8) + (c[i]&0xff);
 
  return nbytes;
}


int load_longlongarray(unsigned char *c, int n, int nbytes, unsigned long long *l)
{
  int i;

  for(i=0; i<n; i++)
    load_longlong(c+nbytes*i,nbytes,l+i);

  return nbytes*n;
}

#endif /* _LONG_LONG */
int store_int(unsigned int l, int nbytes, unsigned char *c)
{
  int i;
  
  for(i=0; i<nbytes; i++)
    c[i] = (l>>(8*(nbytes-i-1)))&0xff;
  
  return nbytes;		/* return number of chars filled */
}


int store_intarray(unsigned int *l, int n, int nbytes, unsigned char *c)
{
  int i;
  
  for(i=0; i<n; i++)
    c += store_int(l[i],nbytes,c);

  return nbytes*n;
}


int load_int(unsigned char *c, int nbytes, unsigned int *l)
{
  int i;
  
  *l = 0;
  
  for(i=0; i<nbytes;i++)
    *l = (*l<<8) + (c[i]&0xff);
 
  return nbytes;
}


int load_intarray(unsigned char *c, int n, int nbytes, unsigned int *l)
{
  int i;
  
  for(i=0; i<n; i++)
    load_int(c+nbytes*i,nbytes,l+i);

  return nbytes*n;
}


  



#if 0

main()
{
  unsigned long long l1[4], l2[4];
  unsigned char c[80], *temp;
  int i;
  
  temp = c;
  
  l1[0] = 0xffaabb01UL;
  l1[1] = 0x01ababUL;
  l1[2] = 0x0b01UL;
  l1[3] = 0xabbffa01UL;
  
  temp += store_longlongarray(l1,2,4,temp);
  temp += store_longlongarray(l1+2,2,4,temp);
  
  temp = c;
  temp += load_longlongarray(temp,3,4,l2);
  temp += load_longlongarray(temp,1,4,l2+3);
  
  for(i=0; i<4; i++)
    printf("%d. \t %llx \t %llx\n", i, l1[i], l2[i]);

  for(i=0; i<4; i++)
    printf("%d. %x %x %x %x\n", i, (unsigned int) c[i*4], (unsigned int) c[i*4+1], (unsigned int) c[i*4+2], (unsigned int) c[i*4+3]);

  
}

#endif /* if 0/1 */
