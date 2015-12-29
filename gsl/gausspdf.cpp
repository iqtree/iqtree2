/* randist/gauss.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2006, 2007 James Theiler, Brian Gough
 * Copyright (C) 2006 Charles Karney
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

//#include <config.h>
#include <math.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

double
gsl_ran_gaussian_pdf (const double x, const double sigma)
{
  double u = x / fabs (sigma);
  double p = (1 / (sqrt (2 * 3.14159265358979323846264338327950288) * fabs (sigma))) * exp (-u * u / 2);
  return p;
}

double
gsl_ran_ugaussian_pdf (const double x)
{
  return gsl_ran_gaussian_pdf (x, 1.0);
}

