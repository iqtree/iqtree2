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

#ifdef HAVE_GETRUSAGE
#include <sys/resource.h>
#include <sys/time.h>
#endif

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
	#else // UNIX
	#include <sys/time.h>
	void gettimeofday(struct timeval* t, void* timezone) {
		time_t cur_time;
		time(&cur_time);
		t->tv_sec = cur_time;
		t->tv_usec = 0;
	}
	#endif
#endif // HAVE_GETTIMEOFDAY



/*********************************************
 * getrusage()
 ********************************************/
#ifndef HAVE_GETRUSAGE
/* Specification.  */
//#include <sys/resource.h>
//#include <sys/time.h>

#include <errno.h>
#include <string.h>

/* Get uint64_t.  */
#include <stdint.h>

#if (defined _WIN32 || defined __WIN32__) && ! defined __CYGWIN__
# include <windows.h>
#else
# include <sys/times.h>
# include <unistd.h>
#endif

int getrusage (int who, struct rusage *usage_p)
{
  if (who == RUSAGE_SELF || who == RUSAGE_CHILDREN)
    {
      /* Clear all unsupported members of 'struct rusage'.  */
      memset (usage_p, '\0', sizeof (struct rusage));

#if (defined _WIN32 || defined __WIN32__) && ! defined __CYGWIN__
      if (who == RUSAGE_SELF)
        {
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
              uint64_t kernel_usec =
                ((((uint64_t) kernel_time.dwHighDateTime << 32)
                  | (uint64_t) kernel_time.dwLowDateTime)
                 + 5) / 10;
              uint64_t user_usec =
                ((((uint64_t) user_time.dwHighDateTime << 32)
                  | (uint64_t) user_time.dwLowDateTime)
                 + 5) / 10;

              usage_p->ru_utime.tv_sec = user_usec / 1000000U;
              usage_p->ru_utime.tv_usec = user_usec % 1000000U;
              usage_p->ru_stime.tv_sec = kernel_usec / 1000000U;
              usage_p->ru_stime.tv_usec = kernel_usec % 1000000U;
            }
        }
#else // UNIX
      /* Fill in the ru_utime and ru_stime members.  */
      {
        struct tms time;

        if (times (&time) != (clock_t) -1)
          {
            /* Number of clock ticks per second.  */
            unsigned int clocks_per_second = sysconf (_SC_CLK_TCK);

            if (clocks_per_second > 0)
              {
                clock_t user_ticks;
                clock_t system_ticks;

                uint64_t user_usec;
                uint64_t system_usec;

                if (who == RUSAGE_CHILDREN)
                  {
                    user_ticks   = time.tms_cutime;
                    system_ticks = time.tms_cstime;
                  }
                else
                  {
                    user_ticks   = time.tms_utime;
                    system_ticks = time.tms_stime;
                  }

                user_usec =
                  (((uint64_t) user_ticks * (uint64_t) 1000000U)
                   + clocks_per_second / 2) / clocks_per_second;
                system_usec =
                  (((uint64_t) system_ticks * (uint64_t) 1000000U)
                   + clocks_per_second / 2) / clocks_per_second;

                usage_p->ru_utime.tv_sec = user_usec / 1000000U;
                usage_p->ru_utime.tv_usec = user_usec % 1000000U;
                usage_p->ru_stime.tv_sec = system_usec / 1000000U;
                usage_p->ru_stime.tv_usec = system_usec % 1000000U;
              }
          }
      }
#endif

      return 0;
    }
  else
    {
      errno = EINVAL;
      return -1;
    }
}
#endif // HAVE_GETRUSAGE

/**
 * @return CPU time in seconds since program was started (corrrect up to micro-seconds)
 * with correction for OpenMP
 */
inline double getCPUTime() {
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	return (usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1.0e6);
}

/**
 * @return real wall-clock time in seconds since Epoch (correct up to micro-seconds)
 */
inline double getRealTime() {
	struct timeval tv;
	if (gettimeofday(&tv, NULL)) return -1.0; // error
	return (tv.tv_sec + (double)tv.tv_usec / 1.0e6);
}

#endif
