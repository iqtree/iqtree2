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

#ifndef _VERSION_H_
#define _VERSION_H_

#include <stdio.h>
#include <stdlib.h>
#ifndef CLANG_UNDER_VS
#include <libgen.h>
#endif

#define NAME "booster"

/** 
    Prints the version of the tools
    In the output file.
    (may be stdout or stderr)
 */
void version(FILE *out, char *executable);

/** 
    Prints the version of the tools
    In the output file.
    Without the executable name
    (may be stdout or stderr)
 */
void short_version(FILE *out);

#endif
