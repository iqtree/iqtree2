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

#ifndef WHTEST_SUB_H
#define WHTEST_SUB_H

typedef struct _knoten { 

	struct _knoten *left;
	struct _knoten *right;
	struct _knoten *up;

	double	edge_length;
	int	label;
	int	ixlabel;
	char	bezeichnung[100];

} knoten;


extern knoten *baum;
extern int **seqData;

extern int simulation, nr_basen, taxa;

extern double alpha, beta;

extern char datei_name[100];
extern char ausgabe_report[200];
extern char ausgabe_nj_tree[200];

extern double **q_matrizen;

extern double **distance;


/*********************************/

void ReadDataSize ( char *datafile );
void ReadData ( char *datafile );

void AllocateMemory();
void FreeMemory();

void Compute_Hij();

void Compute_Qij_tij();

void Write_Qij(int a);

void ComputeNeighborJoiningTree();

void Save_Tree( knoten *P );

void FixDistance();
void Save_Distance(char *distfile, double **dist);

/*void Compute_Qm();*/

void Simulate_Sequences_q_hat();
void Compute_q_hat_pairwise();

/*void Compute_q_hat();*/


/************************************/

#endif
