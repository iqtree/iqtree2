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
#include <sys/time.h>

#ifdef HAVE_GETRUSAGE
	#include <sys/resource.h>
#else 
	#if (defined _WIN32 || defined __WIN32__) && ! defined __CYGWIN__
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
	void gettimeofday(struct timeval* t, void* timezone)
	{       
		struct _timeb timebuffer;
		_ftime( &timebuffer );
		t->tv_sec=timebuffer.time;
		t->tv_usec=1000*timebuffer.millitm;
	}
	#else /* UNIX */
	#include <sys/time.h>
	void gettimeofday(struct timeval* t, void* timezone) {
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
inline double getCPUTime() {
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
inline double getRealTime() {
	struct timeval tv;
	if (gettimeofday(&tv, NULL)) return -1.0; /* error */
	return (tv.tv_sec + (double)tv.tv_usec / 1.0e6);
}

#endif
