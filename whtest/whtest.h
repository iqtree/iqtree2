/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
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

#ifndef WHTEST_H
#define WHTEST_H


int WHTest_run ( int argc, char **argv );

void WHT_setAlignmentSize(int ntax, int nsite);

void WHT_allocateMemory();

void WHT_setSequenceSite(int seqid, int siteid, char c);

void WHT_setSequenceName(int seqid, const char *name);

void WHT_setParams(int nsim, double gamma_shape, char *filename, double *dist);

void WHT_getResults(double *delta, double *delta_quantile, double *p_value);

#endif
