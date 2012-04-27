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


#ifndef EIGEN_SYM_H
#define EIGEN_SYM_H

#define NUM_STATE 4
#define ZERO 0.000001


typedef double DVec20[NUM_STATE];
typedef DVec20 DMat20[NUM_STATE];



/**
	computing eigenvalues and eigenvectors of matrix Pi_vec^(-1) * H_mat, where H_mat is a
	symmetric matrix
	@param H_mat (IN)
	@param Pi_vec (IN)
	@param n (IN) size of matrix
	@param eval (OUT) eigenvalues
	@param evec (OUT) eigenvectors
	@param inv_evec (OUT) inverse matrix of eigenvectors
*/
int eigen_sym(DMat20 H_mat, DVec20 Pi_vec, int num_state, 
	DVec20 eval, DMat20 evec, DMat20 inv_evec);
int eigen_sym_core(double *mat, int n, double *eval);

#endif
