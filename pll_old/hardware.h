#ifndef PLL_HARDWARE
#define PLL_HARDWARE

/* leaf 1 */
/* edx */
#define PLL_HAS_MMX             1 << 23
#define PLL_HAS_SSE             1 << 25
#define PLL_HAS_SSE2            1 << 26

/* ecx */
#define PLL_HAS_SSE3            1
#define PLL_HAS_SSSE3           1 <<  9
#define PLL_HAS_FMA             1 << 12
#define PLL_HAS_SSE41           1 << 19
#define PLL_HAS_SSE42           1 << 20
#define PLL_HAS_AVX             1 << 28


/* leaf 7 */
/* ebx */
#define PLL_HAS_AVX2            1 <<  5

/* leaf 0x80000001 */
/* ecx*/
#define PLL_HAS_SSE4A           1 <<  6
#define PLL_HAS_FMA4            1 << 16

typedef struct
{
  int has_mmx;
  int has_sse;
  int has_sse2;
  int has_sse3;
  int has_ssse3;
  int has_sse41;
  int has_sse42;
  int has_sse4a;
  int has_avx;
  int has_avx2;
  int has_fma;
  int has_fma4;
  int cpu_sockets;
  int cores;
  char vendor[13];

} pllHardwareInfo;

#endif
