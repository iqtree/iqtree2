lib_LTLIBRARIES = libpll-generic.la
lib_LIBRARIES = libpll-generic.a
libpll_generic_la_SOURCES = hash.c stack.c ssort.c queue.c utils.c randomTree.c optimizeModel.c trash.c searchAlgo.c topologies.c fastDNAparsimony.c treeIO.c models.c evaluatePartialGenericSpecial.c evaluateGenericSpecial.c newviewGenericSpecial.c makenewzGenericSpecial.c bipartitionList.c restartHashTable.c recom.c lexer.c alignment.c newick.c parsePartition.c
libpll_generic_la_CFLAGS = -c -O2 -fomit-frame-pointer -funroll-loops -Wall -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef 
libpll_generic_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
libpll_generic_a_SOURCES = $(libpll_generic_la_SOURCES)
libpll_generic_a_CFLAGS = $(libpll_generic_la_CFLAGS)

if BUILD_SSE3
lib_LTLIBRARIES += libpll-sse3.la
libpll_sse3_la_SOURCES = $(libpll_generic_la_SOURCES)
libpll_sse3_la_CFLAGS = -D__SSE3 -msse3 $(libpll_generic_la_CFLAGS)
libpll_sse3_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
lib_LIBRARIES += libpll-sse3.a
libpll_sse3_a_SOURCES = $(libpll_sse3_la_SOURCES)
libpll_sse3_a_CFLAGS = $(libpll_sse3_la_CFLAGS)
endif


if BUILD_AVX
lib_LTLIBRARIES += libpll-avx.la
libpll_avx_la_SOURCES = avxLikelihood.c $(libpll_generic_la_SOURCES)
libpll_avx_la_CFLAGS = -D__SSE3 -D__AVX -mavx -msse3 $(libpll_generic_la_CFLAGS)
libpll_avx_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
lib_LIBRARIES += libpll-avx.a
libpll_avx_a_SOURCES = $(libpll_avx_la_SOURCES)
libpll_avx_a_CFLAGS = $(libpll_avx_la_CFLAGS)
endif

if BUILD_SSE3_PTHREADS
lib_LTLIBRARIES += libpll-sse3-pthreads.la
libpll_sse3_pthreads_la_SOURCES = genericParallelization.c $(libpll_generic_la_SOURCES)
libpll_sse3_pthreads_la_CFLAGS = -D_USE_PTHREADS -D__SSE3 -msse3 $(libpll_generic_la_CFLAGS)
libpll_sse3_pthreads_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
lib_LIBRARIES += libpll-sse3-pthreads.a
libpll_sse3_pthreads_a_SOURCES = $(libpll_sse3_pthreads_la_SOURCES)
libpll_sse3_pthreads_a_CFLAGS = $(libpll_sse3_pthreads_la_CFLAGS)
endif

if BUILD_AVX_PTHREADS
lib_LTLIBRARIES += libpll-avx-pthreads.la
libpll_avx_pthreads_la_SOURCES = avxLikelihood.c genericParallelization.c $(libpll_generic_la_SOURCES)
libpll_avx_pthreads_la_CFLAGS = -D_USE_PTHREADS -D__AVX -mavx -D__SSE3 -msse3 $(libpll_generic_la_CFLAGS)
libpll_avx_pthreads_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
lib_LIBRARIES += libpll-avx-pthreads.a
libpll_avx_pthreads_a_SOURCES = $(libpll_avx_pthreads_la_SOURCES)
libpll_avx_pthreads_a_CFLAGS = $(libpll_avx_pthreads_la_CFLAGS)
endif

pkgincludedir=$(includedir)/pll
pkginclude_HEADERS = pll.h newick.h stack.h hash.h errcodes.h globalVariables.h lexer.h genericParallelization.h treeIO.h queue.h parsePartition.h mem_alloc.h cycle.h
