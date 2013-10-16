#ifndef SSORT_H
#define SSORT_H

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define SWAP(x,y) do{ __typeof__ (x) _t = x; x = y; y = _t; } while(0)

int * ssort1main (char ** x, int n);

#endif
