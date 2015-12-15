#ifndef MIC_NATIVE_H_
#define MIC_NATIVE_H_

void newviewGTRGAMMA_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling);

double evaluateGTRGAMMA_MIC(int *ex1, int *ex2, int *wptr,
                 double *x1_start, double *x2_start,
                 double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable, const pllBoolean fastScaling);

void sumGTRGAMMA_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

void coreGTRGAMMA_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr);

// protein data
void newviewGTRGAMMAPROT_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling);

double evaluateGTRGAMMAPROT_MIC(int *ex1, int *ex2, int *wptr,
                 double *x1_start, double *x2_start,
                 double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable, const pllBoolean fastScaling);

void sumGTRGAMMAPROT_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

void coreGTRGAMMAPROT_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr);

// protein data - LG4

void newviewGTRGAMMAPROT_LG4_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement);

double evaluateGTRGAMMAPROT_LG4_MIC(int *wptr,
                 double *x1_start, double *x2_start,
                 double *tipVector[4],
                 unsigned char *tipX1, const int n, double *diagptable);

void sumGTRGAMMAPROT_LG4_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector[4],
    unsigned char *tipX1, unsigned char *tipX2, int n);

void coreGTRGAMMAPROT_LG4_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN[4], double *gammaRates, double lz, int *wrptr);


#endif /* MIC_NATIVE_H_ */
