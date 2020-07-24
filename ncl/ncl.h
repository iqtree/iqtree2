//	Copyright (C) 1999-2003 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library) version 2.0.
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef NCL_NCL_H
#define NCL_NCL_H

#if defined(_MSC_VER) && !defined(CLANG_UNDER_VS)
#	pragma warning(disable:4786)
#	pragma warning(disable:4291)
#	define vsnprintf _vsnprintf
#endif

#if !defined(__DECCXX)
#	include <cassert>
#	include <cctype>
#	include <cmath>
#	include <cstdarg>
#	include <cstdio>
#	include <cstdarg>
#	include <cstdlib>
#	include <ctime>
#	include <cfloat>
#else
#	include <assert.h>
#	include <ctype.h>
#	include <stdarg.h>
#	include <math.h>
#	include <stdarg.h>
#	include <stdio.h>
#	include <stdlib.h>
#	include <time.h>
#	include <float.h>
#endif

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#if defined(__GNUC__)
#	if __GNUC__ < 3
#		include <strstream>
#	else
#		include <sstream>
#	endif
#endif
#include <vector>
using namespace std;

#if defined(__MWERKS__)
#	if __ide_target("Simple-Win Release") || __ide_target("Phorest-Mac-Release")
#		define NDEBUG
#	else
#		undef NDEBUG
#	endif
#endif

#if defined( __BORLANDC__ )
#	include <dos.h>
#endif

#if defined(__MWERKS__)
#	define HAVE_PRAGMA_UNUSED
		// mwerks (and may be other compilers) want return values even if the function throws an exception
		//
#	define DEMANDS_UNREACHABLE_RETURN

#endif

#include "nxsdefs.h"
#include "nxsstring.h"
#include "nxsexception.h"
#include "nxstoken.h"
#include "nxsblock.h"
#include "nxsreader.h"
#include "nxssetreader.h"
#include "nxstaxablock.h"
#include "nxstreesblock.h"
#include "nxsdistancedatum.h"
#include "nxsdistancesblock.h"
#include "nxsdiscretedatum.h"
#include "nxsdiscretematrix.h"
#include "nxscharactersblock.h"
#include "nxsassumptionsblock.h"
#include "nxsdatablock.h"
#include "nxsemptyblock.h"

#endif
