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

#ifndef WEISSLAMBDA_SUB_H
#define WEISSLAMBDA_SUB_H

extern int simulation, nr_basen, taxa;
extern char datei_name[100];

extern double alpha;
/*
extern char ausgabe_0[200];
extern char ausgabe_1[200];
extern char ausgabe_2[200];
extern char *ausgabe_null;
extern char *ausgabe_data;
*/


/*********************************/

void ReadDataSize();

void AllocateMemory();

void ComputeWilksLambdafromData();

double ComputeWeissLambdaQ16(double **q_16);
int CountValidPairs(double **q_mats);

void sort ( unsigned long n, double arr[] );


/************************************/

#endif
