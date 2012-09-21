#define mymalloc(a) (_mymalloc((a), __LINE__, __FILE__))

#ifndef ANSI_ARGS
#ifdef __STDC__
#define ANSI_ARGS(args) args
#else
#define ANSI_ARGS(args) ()
#endif
#endif

void *_mymalloc ANSI_ARGS((long size, int line, char *message));

