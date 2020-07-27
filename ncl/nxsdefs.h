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
#ifndef NCL_NXSDEFS_H
#define NCL_NXSDEFS_H

#define NCL_NAME_AND_VERSION  "NCL version 2.0"
#define NCL_COPYRIGHT         "Copyright (c) 1999-2003 by Paul O. Lewis"
#define NCL_HOMEPAGEURL       "http://lewis.eeb.uconn.edu/ncl/"

// Maximum number of states that can be stored; the only limitation is that this
// number be less than the maximum size of an int (not likely to be a problem).
// A good number for this is 76, which is 96 (the number of distinct symbols
// able to be input from a standard keyboard) less 20 (the number of symbols
// symbols disallowed by the NEXUS standard for use as state symbols)
//
#define NCL_MAX_STATES         76

#if defined(__MWERKS__) || defined(__DECCXX) || defined(_MSC_VER) 
	typedef long		file_pos;
#else
	typedef streampos	file_pos;
#endif

#define	SUPPORT_OLD_NCL_NAMES

#include <vector>	//for std::vector
#include <set>      //for std::set
#include <map>      //for std::map
#ifdef CLANG_UNDER_VS
#include <xstddef>  //for std::less
#endif 

#include "nxsstring.h"

typedef std::vector<bool>										NxsBoolVector;
typedef std::vector<char>										NxsCharVector;
typedef std::vector<unsigned>									NxsUnsignedVector;
typedef std::vector<NxsStringVector>								NxsAllelesVector;

typedef std::set< unsigned, std::less<unsigned> >						NxsUnsignedSet;

typedef std::map< unsigned, NxsStringVector, std::less<unsigned> >	NxsStringVectorMap;
typedef std::map< NxsString, NxsString, std::less<NxsString> >		NxsStringMap;
typedef std::map< NxsString, NxsUnsignedSet, std::less<NxsString> >	NxsUnsignedSetMap;

// The following typedefs are simply for maintaining compatibility with existing code.
// The names on the right are deprecated and should not be used.
//
typedef	NxsBoolVector		BoolVect;
typedef NxsUnsignedSet		IntSet;
typedef NxsUnsignedSetMap	IntSetMap;
typedef NxsAllelesVector	AllelesVect;
typedef NxsStringVector		LabelList;
typedef NxsStringVector		StrVec;
typedef NxsStringVector		vecStr;
typedef NxsStringVectorMap	LabelListBag;
typedef NxsStringMap		AssocList;

//class NxsTreesBlock;
//class NxsTaxaBlock;
//class NxsAllelesBlock;
//class NxsAssumptionsBlock;
//class NxsCharactersBlock;
//class NxsDistancesBlock;
//class NxsAssumptionsBlock;
//class NxsDiscreteDatum;
//class NxsDiscreteMatrix;

#endif
