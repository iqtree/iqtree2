/***************************************************************************
 *   Copyright (C) 2009 by Gunter Weiss, BUI Quang Minh, Arndt von Haeseler   *
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


#ifndef EIGEN_H
#define EIGEN_H

int xtoy (double x[], double y[], int n);
int matinv( double x[], int n, int m, double space[]);

void matAbyBisC (double A[], double B[], int n, double C[]);

int eigen(int job, double A[], int n, double rr[], double ri[],
          double vr[], double vi[], double w[]);

#endif 
