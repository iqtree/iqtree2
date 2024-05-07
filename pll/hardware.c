#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "systypes.h"

#include <string.h>
#include "hardware.h"

#define PLL_FEAT_AVAIL(x,y) (((x) & (y)) == (y))
#define PLL_SYS_CPU_DIR_PATH "/sys/devices/system/cpu/"

#ifdef CLANG_UNDER_VS
    //James B. Workaround for Windows builds where these macros might not be defined
    #ifndef S_ISDIR
    #define S_ISDIR(mode) (((mode) & S_IFMT) == S_IFDIR)
    #endif
    #ifndef S_ISREG
    #define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
    #endif
#endif

static __inline void cpuid(unsigned int op, int count,
                         unsigned int *eax, unsigned int *ebx,
                         unsigned int *ecx, unsigned int *edx)
{
#if defined(WIN32) || defined(WIN64)
	__int32 regs[4];
	__cpuid((int*)regs, (int)op);
	*eax = regs[0];
	*ebx = regs[1];
	*ecx = regs[2];
	*edx = regs[3];
#elif !defined (__ARM_NEON)
	*eax = op;
  *ecx = count;
  asm volatile("cpuid"
        : "=a" (*eax),
          "=b" (*ebx),
          "=c" (*ecx),
          "=d" (*edx)

        : "0" (*eax), "2" (*ecx)
        : "memory");
#endif
}


void show_hardware_info(pllHardwareInfo * hw)
{
  printf ("MMX.........: %d\n"
          "SSE.........: %d\n"
          "SSE2........: %d\n"
          "SSE3........: %d\n"
          "SSSE3.......: %d\n"
          "FMA.........: %d\n"
          "SSE4.1......: %d\n"
          "SSE4.2......: %d\n"
          "AVX.........: %d\n"
          "AVX2........: %d\n"
          "SSE4A.......: %d\n"
          "FMA4........: %d\n\n"
          "Core(s).....: %d\n"
          "CPU Sockets.: %d\n",

          hw->has_mmx, hw->has_sse, hw->has_sse2, hw->has_sse3, hw->has_ssse3,
          hw->has_fma, hw->has_sse41, hw->has_sse42, hw->has_avx, hw->has_avx2,
          hw->has_sse4a, hw->has_fma4, hw->cores, hw->cpu_sockets);
}

static int pll_probe_cpu (pllHardwareInfo * hw)
{
  struct stat cpustat;
  char cpu[30];
  char cpupath[100];
  int i, id, max_physical_id = -1;
  const char * physical_id_path = "/topology/physical_package_id";
  FILE * fd;

  /* check whether the sys cpu dir exists */
  if (stat(PLL_SYS_CPU_DIR_PATH, &cpustat)) return (0);
  
  /* and also check whether it is a dir */
  if (!S_ISDIR(cpustat.st_mode)) return (0);

  /* detect number of processors */
  for (i = 0; ; ++i)
   {
     sprintf(cpu, "cpu%d", i);
     strcpy (cpupath, PLL_SYS_CPU_DIR_PATH);
     strcat (cpupath, cpu);
     if (stat(cpupath, &cpustat)) break;

     strcat (cpupath, physical_id_path);
     if (!stat(cpupath, &cpustat))
      {
        fd = fopen (cpupath,"r");
        fscanf (fd, "%d", &id);
        /* printf ("Detected processor %d belonging to package %d\n", i, id); */
        if (id > max_physical_id) max_physical_id = id;
        fclose (fd);
      }
   }
  
  hw->cores       = i;
  hw->cpu_sockets = max_physical_id + 1;

  return (1);
}

static void pll_probe_hardware (pllHardwareInfo * hw)
{
#if defined ( __ARM_NEON )               // SSE4.2 for "sse2neon.h"
  hw->vendor[12] = 0;
  hw->has_mmx    = 1;
  hw->has_sse    = 1;
  hw->has_sse2   = 1;
  hw->has_sse3   = 1;
  hw->has_ssse3  = 1;
  hw->has_fma    = 1;
  hw->has_sse41  = 1;
  hw->has_sse42  = 1;
  hw->has_avx    = 0;
  hw->has_avx2   = 0;
  hw->has_sse4a  = 0;
  hw->has_fma    = 0;
#else
  unsigned int a, b, c, d;
  c = 0;

  cpuid(0,0,&a,&b,&c,&d);
  *((unsigned int *)(hw->vendor)    ) = b;
  *((unsigned int *)(hw->vendor + 4)) = d;
  *((unsigned int *)(hw->vendor + 8)) = c;
  hw->vendor[12] = 0;

  printf ("%s\n", hw->vendor);

  cpuid(1,0,&a,&b,&c,&d);

  hw->has_mmx   = PLL_FEAT_AVAIL(d,PLL_HAS_MMX); 
  hw->has_sse   = PLL_FEAT_AVAIL(d,PLL_HAS_SSE);
  hw->has_sse2  = PLL_FEAT_AVAIL(d,PLL_HAS_SSE2);

  hw->has_sse3  = PLL_FEAT_AVAIL(c,PLL_HAS_SSE3);
  hw->has_ssse3 = PLL_FEAT_AVAIL(c,PLL_HAS_SSSE3);
  hw->has_fma   = PLL_FEAT_AVAIL(c,PLL_HAS_FMA);
  hw->has_sse41 = PLL_FEAT_AVAIL(c,PLL_HAS_SSE41);
  hw->has_sse42 = PLL_FEAT_AVAIL(c,PLL_HAS_SSE42);
  hw->has_avx   = PLL_FEAT_AVAIL(c,PLL_HAS_AVX);

  cpuid(7,0,&a,&b,&c,&d);

  hw->has_avx2  = PLL_FEAT_AVAIL(b,PLL_HAS_AVX2);

  /* TODO: note, here we have to check whether leaf 0x80000001 exists */
  cpuid(0x80000001,0,&a,&b,&c,&d);

  hw->has_sse4a = PLL_FEAT_AVAIL(c,PLL_HAS_SSE4A);
  hw->has_fma4  = PLL_FEAT_AVAIL(c,PLL_HAS_FMA4);
#endif
}

int pllGetHardwareInfo (pllHardwareInfo * hw)
{
  pll_probe_hardware (hw);
  pll_probe_cpu (hw);

  /* TODO: finish failure checks in probe_hardware and probe_cpu */
  return (1);

}

/* TODO: Remove after testing */
/* 
int main (int argc, char * argv[])
{ 
  pllHardwareInfo hw;

  pll_probe_hardware(&hw);
  pll_probe_cpu(&hw);

  show_hardware_info(&hw);
  return (EXIT_SUCCESS);
}
*/
