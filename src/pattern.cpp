//
// C++ Implementation: pattern
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pattern.h"

Pattern::Pattern()
 : string()
{
	frequency = 0;
	is_const = false;
}


Pattern::~Pattern()
{
}


void Pattern::computeConst() {
	char ch = at(0);
	is_const = true;
	for (iterator i = begin(); i != end(); i++)
		if (*i != ch) {
			is_const = false;
			return;
		}
}
