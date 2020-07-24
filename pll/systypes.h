#ifdef WIN32
#include <direct.h>
#endif

#if !defined(WIN32) && !defined(WIN64)
#include <unistd.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#endif
