/***************************************************************************
 *   Copyright (C) 2009 by Gunter Weiss, Bui Quang Minh, Arndt von Haeseler   *
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "random.h"
#include "whtools.h"

unsigned int kiss(void);

/****************************************************************************************/ 
double sexp(void) 
{ 
    /* q[k-1] = sum(aLOG(2.0)**k/k!) k=1,..,n, */ 
    /* The highest n (here 8) is determined by q[n-1] = 1.0 */ 
    /* within standard precision */ 
    static double q[] = 
    { 
		0.6931471805599453, 
			0.9333736875190459, 
			0.9888777961838675, 
			0.9984959252914960, 
			0.9998292811061389, 
			0.9999833164100727, 
			0.9999985691438767, 
			0.9999998906925558, 
			0.9999999924734159, 
			0.9999999995283275, 
			0.9999999999728814, 
			0.9999999999985598, 
			0.9999999999999289, 
			0.9999999999999968, 
			0.9999999999999999, 
			1.0000000000000000 
    }; 
    double a, u, ustar, umin; 
    int i; 
	 
    a = 0.0; 
    u = ranDum(); 
    for (;;) { 
		u = u + u; 
		if (u > 1.0) 
			break; 
		a = a + q[0]; 
    } 
    u = u - 1.0; 
	 
    if (u <= q[0]) 
		return a + u; 
	 
    i = 0; 
    ustar = ranDum(); 
    umin = ustar; 
    do { 
		ustar = ranDum(); 
		if (ustar < umin) 
			umin = ustar; 
		i = i + 1; 
    } 
    while (u > q[i]); 
    return a + umin * q[0]; 
} 
/******************************************************************/
/*  variables for the kiss random number generator                */
/******************************************************************/
unsigned int k,m,x,y,z,w,carry,r;
/******************************************************************/
/*  kickstart the kiss random number generator                    */
/******************************************************************/
void start_kiss(int seed)
{
#  ifdef PARALLEL
	int n;
#  endif

     x=seed;y=102;z=12;w=34535;
     x = x * 69069 + 1;
     y ^= y << 13;
     y ^= y >> 17;
     y ^= y << 5;
     k = (z >> 2) + (w >> 3) + (carry >> 2);
     m = w + w + z + carry;
     z = w;
     w = m;
     carry = k >> 30;

#  ifdef PARALLEL
	for (n=0; n<mpi_myrank; n++)
		kiss();
#  endif

}
/******************************************************************/
void restart_kiss(unsigned int *vals)
{
	k=vals[0];
	m=vals[1];
	x=vals[2];
	y=vals[3];
	z=vals[4];
	w=vals[5];
	carry=vals[6];
	r=vals[7];
}
/******************************************************************/
void kiss_state(unsigned int *vals)
{
	vals[0]=k;
	vals[1]=m;
	vals[2]=x;
	vals[3]=y;
	vals[4]=z;
	vals[5]=w;
	vals[6]=carry;
	vals[7]=r;
}
/******************************************************************/
/*   Keep It Simple Stupid random number generator from George 
     Marsaglia's DIEHARD cdrom                                    */
/******************************************************************/
unsigned int single_kiss(void)
{
     x = x * 69069 + 1;
     y ^= y << 13;
     y ^= y >> 17;
     y ^= y << 5;
     k = (z >> 2) + (w >> 3) + (carry >> 2);
     m = w + w + z + carry;
     z = w;
     w = m;
     carry = k >> 30;
     return x+y+z;
}

unsigned int kiss(void) {
#ifdef PARALLEL
	int i;
	for (i = 1; i < mpi_size; i++)
		single_kiss();
#endif

	return single_kiss();
}

/******************************************************************/
double dkiss(void)
{
    return ((double)kiss()+0.5)/4294967296.0;
}
/******************************************************************/
/*************************************************************************/ 
double stat_normal(void) 
{ 
	static int iset = 0; 
	static double gset; 
	double  fac,rsq,v1,v2; 
	 
	if (iset == 0) { 
		do { 
			v1 = 2.0*ranDum( )-1.0; 
			v2 = 2.0*ranDum( )-1.0; 
			rsq = v1*v1+v2*v2; 
		} while (rsq >= 1.0 || rsq == 0); 
		fac = sqrt(-2.0*log(rsq)/rsq); 
		gset = v1 * fac; 
		iset = 1; 
		return v2*fac; 
	} else { 
		iset = 0; 
		return gset; 
	} 
} 
/****************************************************************/ 
static double a1 = 0.3333333; 
static double a2 = -0.250003; 
static double a3 = 0.2000062; 
static double a4 = -0.1662921; 
static double a5 = 0.1423657; 
static double a6 = -0.1367177; 
static double a7 = 0.1233795; 
static double e1 = 1.0; 
static double e2 = 0.4999897; 
static double e3 = 0.166829; 
static double e4 = 0.0407753; 
static double e5 = 0.010293; 
static double q1 = 0.04166669; 
static double q2 = 0.02083148; 
static double q3 = 0.00801191; 
static double q4 = 0.00144121; 
static double q5 = -7.388e-5; 
static double q6 = 2.4511e-4; 
static double q7 = 2.424e-4; 
static double sqrt32 = 5.656854; 
static double aa = 0.; 
static double aaa = 0.; 
double rgamma(double a, double scale) 
#define repeat for(;;) 
/* Taken from R */ 
{ 
	static double b, c, d, e, p, q, r, s, t, u, v, w, x; 
	static double q0, s2, si; 
	double ret_val; 
	 
	if (a < 1.0) { 
		/* alternate method for parameters a below 1 */ 
		/* 0.36787944117144232159 = exp(-1) */ 
		aa = 0.0; 
		b = 1.0 + 0.36787944117144232159 * a; 
		repeat { 
			p = b * dkiss(); 
			if (p >= 1.0) { 
				ret_val = -log((b - p) / a); 
				if (sexp() >= (1.0 - a) * log(ret_val)) 
					break; 
			} else { 
				ret_val = exp(log(p) / a); 
				if (sexp() >= ret_val) 
					break; 
			} 
		} 
		return scale * ret_val; 
	} 
	/* Step 1: Recalculations of s2, s, d if a has changed */ 
	if (a != aa) { 
		aa = a; 
		s2 = a - 0.5; 
		s = sqrt(s2); 
		d = sqrt32 - s * 12.0; 
	} 
	/* Step 2: t = standard normal deviate, */ 
	/* x = (s,1/2)-normal deviate. */ 
	/* immediate acceptance (i) */ 
	 
	t = stat_normal(); 
	x = s + 0.5 * t; 
	ret_val = x * x; 
	if (t >= 0.0) 
		return scale * ret_val; 
	 
	/* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */ 
	u = dkiss(); 
	if (d * u <= t * t * t) { 
		return scale * ret_val; 
	} 
	/* Step 4: recalculations of q0, b, si, c if necessary */ 
	 
	if (a != aaa) { 
		aaa = a; 
		r = 1.0 / a; 
		q0 = ((((((q7 * r + q6) * r + q5) * r + q4) 
			* r + q3) * r + q2) * r + q1) * r; 
		 
		/* Approximation depending on size of parameter a */ 
		/* The constants in the expressions for b, si and */ 
		/* c were established by numerical experiments */ 
		 
		if (a <= 3.686) { 
			b = 0.463 + s + 0.178 * s2; 
			si = 1.235; 
			c = 0.195 / s - 0.079 + 0.16 * s; 
		} else if (a <= 13.022) { 
			b = 1.654 + 0.0076 * s2; 
			si = 1.68 / s + 0.275; 
			c = 0.062 / s + 0.024; 
		} else { 
			b = 1.77; 
			si = 0.75; 
			c = 0.1515 / s; 
		} 
	} 
	/* Step 5: no quotient test if x not positive */ 
	 
	if (x > 0.0) { 
		/* Step 6: calculation of v and quotient q */ 
		v = t / (s + s); 
		if (fabs(v) <= 0.25) 
			q = q0 + 0.5 * t * t * ((((((a7 * v + a6) 
			* v + a5) * v + a4) * v + a3) 
			* v + a2) * v + a1) * v; 
		else 
			q = q0 - s * t + 0.25 * t * t + (s2 + s2) 
			* log(1.0 + v); 
		 
		 
		/* Step 7: quotient acceptance (q) */ 
		 
		if (log(1.0 - u) <= q) 
			return scale * ret_val; 
	} 
	/* Step 8: e = standard exponential deviate */ 
	/* u= 0,1 -uniform deviate */ 
	/* t=(b,si)-double exponential (laplace) sample */ 
	 
	repeat { 
		e = sexp(); 
		u = dkiss(); 
		u = u + u - 1.0; 
		if (u < 0.0) 
			t = b - si * e; 
		else 
			t = b + si * e; 
		/* Step  9:  rejection if t < tau(1) = -0.71874483771719 */ 
		if (t >= -0.71874483771719) { 
			/* Step 10:  calculation of v and quotient q */ 
			v = t / (s + s); 
			if (fabs(v) <= 0.25) 
				q = q0 + 0.5 * t * t * ((((((a7 * v + a6) 
				* v + a5) * v + a4) * v + a3) 
				* v + a2) * v + a1) * v; 
			else 
				q = q0 - s * t + 0.25 * t * t + (s2 + s2) 
				* log(1.0 + v); 
			/* Step 11:  hat acceptance (h) */ 
			/* (if q not positive go to step 8) */ 
			if (q > 0.0) { 
				if (q <= 0.5) 
					w = ((((e5 * q + e4) * q + e3) 
					* q + e2) * q + e1) * q; 
				else 
					w = exp(q) - 1.0; 
				/* if t is rejected */ 
				/* sample again at step 8 */ 
				if (c * fabs(u) <= w * exp(e - 0.5 * t * t)) 
					break; 
			} 
		} 
	} 
	x = s + 0.5 * t; 
	return scale * x * x; 
} 




/*******************************************************************/
