/*

BOOSTER: BOOtstrap Support by TransfER: 
BOOSTER is an alternative method to compute bootstrap branch supports 
in large trees. It uses transfer distance between bipartitions, instead
of perfect match.

Copyright (C) 2017 Frederic Lemoine, Jean-Baka Domelevo Entfellner, Olivier Gascuel

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef PRNG_H_INCLUDED
#define PRNG_H_INCLUDED

#include <stddef.h>

long prng_seed_time (void);
void prng_seed_bytes (const void *, size_t);
unsigned char prng_get_octet (void);
unsigned char prng_get_byte (void);
void prng_get_bytes (void *, size_t);
unsigned long prng_get_ulong (void);
long prng_get_long (void);
unsigned prng_get_uint (void);
int prng_get_int (void);
double prng_get_double (void);
double prng_get_double_normal (void);

#endif /* prng.h */
