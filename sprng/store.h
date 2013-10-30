#ifndef _STORE_H
#define _STORE_H

/* Numbers are stored with most significant bit first (left most) */

int store_long(unsigned long l, int nbytes, unsigned char *c);
int store_longarray(unsigned long *l, int n, int nbytes, unsigned char *c);
int load_long(unsigned char *c, int nbytes, unsigned long *l);
int load_longarray(unsigned char *c, int n, int nbytes, unsigned long *l);
int store_int(unsigned int l, int nbytes, unsigned char *c);
int store_intarray(unsigned int *l, int n, int nbytes, unsigned char *c);
int load_int(unsigned char *c, int nbytes, unsigned int *l);
int load_intarray(unsigned char *c, int n, int nbytes, unsigned int *l);


#ifdef _LONG_LONG

int store_longlong(unsigned long long l, int nbytes, unsigned char *c);
int store_longlongarray(unsigned long long *l, int n, int nbytes, unsigned char *c);
int load_longlong(unsigned char *c, int nbytes, unsigned long long *l);
int load_longlongarray(unsigned char *c, int n, int nbytes, unsigned long long *l);

#endif /* _LONG_LONG */


#endif /* _STORE_H */

