/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef TIMEUTIL_H
#define TIMEUTIL_H

#include <iqtree_config.h>
#include <stdlib.h>

#include <errno.h>
#include <string.h>
#include <stdint.h>
#if !defined(_MSC_VER)
#include <sys/time.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

//#if defined(_MSC_VER)
//#define inline __inline
//#endif

#if (defined _WIN32 || defined __WIN32__ || defined WIN32 || defined WIN64) 
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x500
#endif
#endif

#ifdef HAVE_GETRUSAGE
	#include <sys/resource.h>
#else 
	#if (defined _WIN32 || defined __WIN32__ || defined WIN64) && ! defined __CYGWIN__
	# include <windows.h>
	#else
	# include <sys/times.h>
	# include <unistd.h>
	#endif
#endif /* HAVE_GETRUSAGE */

/*********************************************
 * gettimeofday()
 ********************************************/
#ifndef HAVE_GETTIMEOFDAY
	#if defined WIN32 || defined _WIN32 || defined __WIN32__
	#include <sys/timeb.h>
	#include <sys/types.h>
	#include <winsock.h>

	struct timezone {
		char dummy;
	};

	__inline void gettimeofday(struct timeval* t, void* timezone)
	{       
		struct _timeb timebuffer;
		_ftime( &timebuffer );
		t->tv_sec=timebuffer.time;
		t->tv_usec=1000*timebuffer.millitm;
	}
	#else /* UNIX */
	#include <sys/time.h>
	__inline void gettimeofday(struct timeval* t, void* timezone) {
		time_t cur_time;
		time(&cur_time);
		t->tv_sec = cur_time;
		t->tv_usec = 0;
	}
	#endif
#endif /* HAVE_GETTIMEOFDAY */



/**
 * @return CPU time in seconds since program was started (corrrect up to micro-seconds)
 * with correction for OpenMP
 */
__inline double getCPUTime() {
#ifdef HAVE_GETRUSAGE
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	return (usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1.0e6);
#elif (defined _WIN32 || defined __WIN32__) && ! defined __CYGWIN__
	/* Fill in the ru_utime and ru_stime members.  */
	FILETIME creation_time;
	FILETIME exit_time;
	FILETIME kernel_time;
	FILETIME user_time;

	if (GetProcessTimes (GetCurrentProcess (),
						&creation_time, &exit_time,
						&kernel_time, &user_time))
	{
		/* Convert to microseconds, rounding.  */
		uint64_t user_usec = ((((uint64_t) user_time.dwHighDateTime << 32) | (uint64_t) user_time.dwLowDateTime) + 5) / 10;
		return (double)user_usec / 1.0e6;
	}
#else
	/* Fill in the ru_utime and ru_stime members.  */
	struct tms time;

	if (times (&time) != (clock_t) -1) {
		unsigned int clocks_per_second = sysconf (_SC_CLK_TCK);
		if (clocks_per_second > 0) {
			uint64_t user_usec;
			user_usec =	(((uint64_t) time.tms_utime * (uint64_t) 1000000U) + clocks_per_second / 2) / clocks_per_second;
			return (double)user_usec / 1.0e6;
		}
	}
#endif
	abort();
}

/**
 * @return real wall-clock time in seconds since Epoch (correct up to micro-seconds)
 */
__inline double getRealTime() {
#ifdef _OPENMP
	return omp_get_wtime();
#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	//Tung: the if statement below causes compiling error because gettimeofday() return void not boolean
	//if (gettimeofday(&tv, NULL)) return -1.0; /* error */
	return (tv.tv_sec + (double)tv.tv_usec / 1.0e6);
#endif
}
/*
#if defined _WIN32 || defined __WIN32__ || defined WIN32
#include <windows.h>
#include <winbase.h>
inline uint64_t getTotalSystemMemory()
{
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
}

#elif defined __APPLE__ || defined __MACH__

#include <sys/types.h>
#include <sys/sysctl.h>

inline uint64_t getTotalSystemMemory()
{
	int mib[2];
	uint64_t physical_memory;
	mib[0] = CTL_HW;
	mib[1] = HW_MEMSIZE;
	size_t length = sizeof(uint64_t);
	sysctl(mib, 2, &physical_memory, &length, NULL, 0);
	return physical_memory;
}
#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/sysinfo.h>

inline uint64_t getTotalSystemMemory()
{
    struct sysinfo memInfo;
	sysinfo (&memInfo);
	int64_t totalram = memInfo.totalram;
	return (totalram * memInfo.mem_unit);
}

#endif*/ /* for declaring getTotalSystemMemory() */

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <Windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#else
#error "Unable to define getMemorySize( ) for an unknown OS."
#endif



/**
 * Returns the size of physical memory (RAM) in bytes.
 */
__inline uint64_t getMemorySize( )
{
#if defined(_WIN32) && (defined(__CYGWIN__) || defined(__CYGWIN32__) || !defined(_WIN64))
	/* Cygwin under Windows. ------------------------------------ */
	/* New 64-bit MEMORYSTATUSEX isn't available.  Use old 32.bit */
#warning "getMemorySize() will be wrong if RAM is actually > 4GB"
	MEMORYSTATUS status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatus( &status );
	return (uint64_t)status.dwTotalPhys;

#elif defined(_WIN32)
	/* Windows. ------------------------------------------------- */
	/* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx( &status );
	return (uint64_t)status.ullTotalPhys;

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* UNIX variants. ------------------------------------------- */
	/* Prefer sysctl() over sysconf() except sysctl() HW_REALMEM and HW_PHYSMEM */

#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))
	int mib[2];
	mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
	mib[1] = HW_MEMSIZE;            /* OSX. --------------------- */
#elif defined(HW_PHYSMEM64)
	mib[1] = HW_PHYSMEM64;          /* NetBSD, OpenBSD. --------- */
#endif
	uint64_t size = 0;               /* 64-bit */
	size_t len = sizeof( size );
	if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
		return (uint64_t)size;
	return 0L;			/* Failed? */

#elif defined(_SC_AIX_REALMEM)
	/* AIX. ----------------------------------------------------- */
	return (uint64_t)sysconf( _SC_AIX_REALMEM ) * (uint64_t)1024L;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
	/* FreeBSD, Linux, OpenBSD, and Solaris. -------------------- */
	return (uint64_t)sysconf( _SC_PHYS_PAGES ) *
		(uint64_t)sysconf( _SC_PAGESIZE );

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
	/* Legacy. -------------------------------------------------- */
	return (uint64_t)sysconf( _SC_PHYS_PAGES ) *
		(uint64_t)sysconf( _SC_PAGE_SIZE );

#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))
	/* DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX. -------- */
	int mib[2];
	mib[0] = CTL_HW;
#if defined(HW_REALMEM)
	mib[1] = HW_REALMEM;		/* FreeBSD. ----------------- */
#elif defined(HW_PYSMEM)
	mib[1] = HW_PHYSMEM;		/* Others. ------------------ */
#endif
	uint64_t size = 0;		/* 32-bit */
	size_t len = sizeof( size );
	if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
		return (uint64_t)size;
	return 0L;			/* Failed? */
#endif /* sysctl and sysconf variants */

#else
	return 0L;			/* Unknown OS. */
#endif
}


#define HOW_LONG(x) \
{ std::cout.precision(6); double startTime = getRealTime(); \
x; \
std::cout << #x << " took " << (getRealTime()-startTime) << std::endl;  }

#endif
