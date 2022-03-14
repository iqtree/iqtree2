/***************************************************************************
 *   Copyright (C) 2022, by                                                *
 *   James Barbetti <james_barbetti@yahoo.com>                             *
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

#ifndef scoped_assignment_h
#define scoped_assignment_h

template <class C> class scoped_assignment {
    C& variable;
    C  original_value;

public:
    scoped_assignment(C& var, C value)
        : variable(var), original_value(var) {
        variable = value;
    }
    ~scoped_assignment() {
        variable = original_value;
    }
};

//Notes: 
// 1. SCOPED_ASSIGN requires C++11 or later (for decltype).
// 2. The two CONCATENATED_NAME macros generate a unique variable 
//    name, for each invocation of SCOPED_ASSIGN.
// 3. The name is of the form scoped_restore_nnn, and doesn't 
//    include var, because var might be an expression 
//    that evaluates to an l-value rather than a variable name.
// 4. SCOPED_ASSIGN(x,y); is roughly equivaleent to a
//    declaration of the form:
//    scoped_assignment<decltype(x)> scoped_restore_123(x,y);
//
#define CONCATENATED_NAME_1(a,b) a ## b
#define CONCATENATED_NAME(a,b) CONCATENATED_NAME_1(a,b)
#define SCOPED_ASSIGN(var,value) \
    scoped_assignment<decltype(var)> \
        CONCATENATED_NAME(scoped_restore_, __COUNTER__ ) \
        (var, value) 

#endif
