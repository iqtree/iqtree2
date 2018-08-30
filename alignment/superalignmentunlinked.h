/***************************************************************************
 *   Copyright (C) 2018 by BUI Quang Minh   *
 *   m.bui@anu.edu.au   *
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
#ifndef SUPERALIGNMENTUNLINKED_H
#define SUPERALIGNMENTUNLINKED_H

#include "superalignment.h"

class SuperAlignmentUnlinked : public SuperAlignment
{
public:
    /** constructor initialize from a supertree */
    SuperAlignmentUnlinked(Params &params);
    
    /** constructor initialize empty alignment */
    SuperAlignmentUnlinked();
    
    virtual void init(StrVector *sequence_names = NULL);
    
};

#endif

