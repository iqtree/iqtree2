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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "whtest_sub.h"
#include "eigen.h"
#include "eigen_sym.h"
#include "random.h"
#include "whtools.h"


knoten *baum;
int **seqData;
double ****H;
double **q_matrizen;
double **distance;
double *alpha_rate;
double q_hat_eigen[4], U_q_hat[16], V_q_hat[16], statPi[4];

/*************************************************/

void ReadDataSize ( char *datafile )
{
	char c;

	FILE *ifp;

	if ( ( ifp = fopen ( datafile, "r" ) ) == NULL ) {
		printf ( "\nERROR: Missing input file %s!\n", datafile );
		Finalize ( 1 );
	}


	if ( fscanf ( ifp, "%d", &taxa ) != 1 ) {
		printf ( "\nERROR: Missing number of taxa!\n" );
		Finalize ( 1 );
	}

	if ( fscanf ( ifp, "%d", &nr_basen ) != 1 ) {
		printf ( "\nERROR: Missing number of sites!\n" );
		Finalize ( 1 );
	}

	do	
		c = fgetc ( ifp ); 	/* skip rest of line */
	while ( c != '\n' );

	fclose ( ifp );

	/*fprintf(stderr,"\ntaxa: %d\t basen: %d\n", taxa,nr_basen);*/
}

/*************************************************/


void AllocateMemory()
{

	int i, j, k;

	/* sequence data */

	seqData = ( int ** ) malloc ( ( 2*taxa-1 ) * sizeof ( int * ) );
	for ( i = 0; i < 2*taxa-1; i++ )
		seqData[i] = ( int * ) calloc ( nr_basen, sizeof ( int ) );


	/* baum structure */

	baum = ( knoten * ) malloc ( ( 2*taxa-1 ) * sizeof ( knoten ) );


	/* distance matrix  */

	distance = ( double ** ) malloc ( taxa * sizeof ( double * ) );
	for ( i = 0; i < taxa; i++ )
		distance[i] = ( double * ) calloc ( taxa, sizeof ( double ) );

	/* divergence matrices */

	H = ( double **** ) malloc ( taxa * sizeof ( double *** ) );
	for ( i = 0; i < taxa; i++ ) {
		H[i] = ( double *** ) malloc ( taxa * sizeof ( double ** ) );

		for ( j = 0; j < taxa; j++ ) {
			H[i][j] = ( double ** ) malloc ( 5 * sizeof ( double * ) );

			for ( k = 0; k < 5; k++ )
				H[i][j][k] = ( double * ) calloc ( 5, sizeof ( double ) );
		}
	}

	/* pairwise rate matrices */

	q_matrizen = ( double ** ) malloc ( taxa* ( taxa-1 ) /2 *sizeof ( double * ) );
	for ( i = 0; i < ( int ) ( taxa* ( taxa-1. ) /2. ); i++ ) {
		q_matrizen[i] = ( double * ) calloc ( 16, sizeof ( double ) );
	}


	/* raten_heterogenitaet */

	alpha_rate = ( double * ) calloc ( nr_basen, sizeof ( double ) );

}

void FreeMemory() {
	int i, j, k;

	/* raten_heterogenitaet */
	free(alpha_rate);

	/* pairwise rate matrices */
	for ( i = ( int ) ( taxa* ( taxa-1. ) /2. ) - 1; i >= 0; i-- ) {
		free(q_matrizen[i]);
	}
	free(q_matrizen);

	/* divergence matrices */

	for ( i = taxa-1; i>=0; i-- ) {

		for ( j = taxa-1; j >= 0; j-- ) {
			for ( k = 4; k >= 0; k-- )
				free(H[i][j][k]);
			free(H[i][j]);
		}
		free(H[i]);
	}

	free(H);

	/* distance matrix  */

	for ( i = taxa-1; i >= 0; i-- )
		free(distance[i]);
	free(distance);


	/* TODO baum structure */

	free(baum);
	/* sequence data */

	for ( i = 2*taxa-2; i >= 0; i-- )
		free(seqData[i]);

	free(seqData);


}


/*************************************************/


void ReadData ( char *datafile )
{
	int i, j;
	char c;

	FILE *ifp;

	if ( ( ifp = fopen ( datafile, "r" ) ) == NULL ) {
		if (isMasterProc())
			printf ( "\nERROR: Missing input file!\n" );
	}


	do	
		c = fgetc ( ifp ); 	/* skip 1st line */
	while ( c != '\n' );

	for ( i = 0; i < taxa; i++ ) {

		/*	do	c = fgetc(ifp); */	 /* skip sequence name */
		/*	while ( c != '\n' && c != ' ');
		*/
		j = 0;
		while ( j < 10 ) {
			fscanf ( ifp, "%c", &c );

			if ( c != '\n' && c != ' ' )
				baum[i].bezeichnung[j] = c;
			else	{
				baum[i].bezeichnung[j] = '\0';
				j = 10;
			}

			j++;
		}

		if (isMasterProc())
			printf("%3i\t%s\n", i+1, baum[i].bezeichnung);
		j = 0;
		while ( j < nr_basen ) {
			fscanf ( ifp, "%c", &c );
			c = toupper ( c );

			if ( c == 'A' || c == '0' )	{seqData[i][j] = 0;  j++;}
			else if ( c == 'C' || c== '1' )	{seqData[i][j] = 1;  j++;}
			else if ( c == 'G' || c== '2' )	{seqData[i][j] = 2;  j++;}
			else if ( c == 'T' || c== '3' )	{seqData[i][j] = 3;  j++;}
			else if ( c == '-' )	{seqData[i][j] = 4;  j++;}
			else if ( c == 'N' )     {seqData[i][j] = 4;  j++;}
			else if ( c == ' ' ) ;
			else if ( c == '\n' );
			else {
				if (isMasterProc())
					fprintf ( stderr,"\nERROR: wrong BASE in datafile!   %c\n", c );
				seqData[i][j] = 5;
				j++;
			}

		}

		if (j != nr_basen) {
			if (isMasterProc())
				printf("ERROR: %s has only %i characters\n", baum[i].bezeichnung, j);
			Finalize(1);
		}

		do	
			c = fgetc ( ifp ); 	/* skip rest of line */
		while ( c != '\n' );
	}

	fclose ( ifp );

}

/*********************************************************************/
#define JMAX 20

/* newton raphson method */
float rtnewt(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc)
{
	int j;
	float df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=1;j<=JMAX;j++) {
		(*funcd)(rtn,&f,&df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			printf("Jumped out of brackets in rtnewt");
		if (fabs(dx) < xacc) return rtn;
	}
	printf("Maximum number of iterations exceeded in rtnewt");
	return 0.0;
}


void FixDistance() {
	int i, j, k, m;
	double T1[16], pi[4], coeff[4], coeff_eigen[4], identity;
	double f, df, rtn, x1 = 0.000001, x2 = 10.0, dx, expf;
	printf("Computing corrected distance matrix based on averaged Q\n");
	
	for (i = 0; i < taxa-1; i++) 
		for (j = i+1; j < taxa; j++) {
			/* get the state frequency pi */
			pi[0] = statPi[0];
			for (k = 1; k < 4; k++) {
				pi[k] = statPi[k] - statPi[k-1];
			}
			matAbyBisC ( U_q_hat, V_q_hat, 4, T1 );
			for (k = 0; k < 4; k++) 
				for (m = 0, coeff[k] = 0.0; m < 4; m++) 
					coeff[k] += pi[m] * U_q_hat[k*4+m] * V_q_hat[m*4+k];


			for (k = 0, identity = 0.0; k < 4; k++) {
				identity += H[i][j][k][k];
				coeff_eigen[k] = coeff[k] * q_hat_eigen[k];
			}
			/* with this transformation, we need to solve the equation f(t) = 0 with
			  f(t) = sum_k {coeff[k]*exp(eigen[k]*t)} - identity
			  the derivative:
			  f'(t) = sum_k {eigen[k]*coeff[k]*exp(eigen[k]*t)}
			  In the following we use Newton-Raphson to find the root of f(t)
			  */

			/* first guess of the distance */
			if (distance[i][j] < 10.0) 
				rtn = distance[i][j];
			else 
				rtn = -(3.0/4.0) * log(1.0 - 4.0/3.0 *(1.0 - identity)); /* Juke-Cantor corrected distance */
			for (k = 1; k <= JMAX; k++) {
				/* compute f(x) and f'(x) at x = rtn */
				for (m = 0, f = -identity, df = 0.0; m < 4; m++) {
					expf = exp(q_hat_eigen[m] * rtn);
					f += coeff[m] * expf;
					df += coeff_eigen[m] * expf;
				}
				dx=f/df;
				rtn -= dx;
				if ((x1-rtn)*(rtn-x2) < 0.0)
					printf("Jumped out of brackets in rtnewt");
				if (fabs(dx) < 0.0001) break;
			}
			distance[i][j] = rtn;
			distance[j][i] = rtn;
		}
}
#undef JMAX

/*********************************************************************/
void FixDistance_old() {
	int i, j, k;
	double T1[16], T2[16];
	
	for (i = 0; i < taxa-1; i++) 
		for (j = i+1; j < taxa; j++) 
			if (distance[i][j] >= 100) {
				for ( k = 0; k < 16; k++ )
					T2[k] = 0;
	
				if ( alpha > 10 )	{
					/* keine ratenheterogenitaet */
	
					T2[0] = q_hat_eigen[0];
					T2[5] = q_hat_eigen[1];
					T2[10] = q_hat_eigen[2];
					T2[15] = q_hat_eigen[3];
				} else 	{
					/* ratenheterogenitaet */
					/* something wrong here! */
				T2[0] = alpha* ( 1.- exp(-q_hat_eigen[0]/alpha ) );
				T2[5] = alpha* ( 1.- exp(-q_hat_eigen[1]/alpha ) );
				T2[10] = alpha* ( 1.-exp(-q_hat_eigen[2]/alpha ) );
				T2[15] = alpha* ( 1.-exp(-q_hat_eigen[3]/alpha ) );
				}
	
				/* T2 = U * diag(eigenwert) * V */
	
				matAbyBisC ( U_q_hat, T2, 4, T1 );
				matAbyBisC ( T1, V_q_hat, 4, T2 );
	
				/* normalisieren durch t = -sum (pi_i * r_ii) , Q hat rate 1 */
	
				distance[i][j] = 0;
				for ( k = 0; k < 4; k++ )
					distance[i][j] -= H[i][j][k][4]*T2[k*5];
	
				distance[j][i] = distance[i][j];
/*
				if (isMasterProc())
					printf("Fix distance (%s,%s) -> %f\n", baum[i].bezeichnung, baum[j].bezeichnung, distance[i][j]);*/
				if (distance[i][j] > 100) {
					if (isMasterProc())
						printf("ERROR: too large distance, try higher alpha please\n");
					Finalize(1);
				}
			}
}

/*********************************************************************/

void Save_Distance(char *distfile, double **dist) {
	FILE *fps;
	int i, j;

	if ((fps = fopen(distfile, "w")) == NULL) {
		printf ( "\nERROR: Cannot write to file %s!\n", distfile );
	}

	fprintf(fps, "%d\n", taxa);

	for (i = 0; i < taxa; i++) {
		fprintf(fps, "%-10s", baum[i].bezeichnung);
		for (j = 0; j < taxa; j++) 
			fprintf(fps, " %f", dist[i][j]);
		fprintf(fps, "\n");
	}

	fclose(fps);
}


/*********************************************************************/


void Compute_Hij()
{

	int i, j, k, l;

	for ( i = 0; i < taxa; i++ ) {
		for ( j = 0; j < taxa; j++ ) {
			for ( k = 0; k < 5; k++ ) {
				for ( l = 0; l < 5; l++ )
					H[i][j][k][l] = 0;
			}


			for ( k = 0; k < nr_basen; k++ ) {
				H[i][j][seqData[i][k]][seqData[j][k]]+=1.;
				H[i][j][seqData[j][k]][seqData[i][k]]+=1.;
				/* symmetrisiert, da reversibilitaet vorausgesetzt */
			}
		}
	}


	for ( i = 0; i < taxa; i++ ) {
		for ( j = 0; j < taxa; j++ ) {
			/* H[i][j][k][4] und H[i][j][4][k] enthalten basenhaeufigkeiten */

			for ( k = 0; k < 4; k++ ) {
				H[i][j][k][4] =	H[i][j][k][0]+H[i][j][k][1]+H[i][j][k][2]+H[i][j][k][3];
				H[i][j][4][k] = H[i][j][k][4];
			}

			/* H[i][j][4][4] enthaelt nr_basen ohne gaps */
			H[i][j][4][4] = H[i][j][0][4]+H[i][j][1][4]+H[i][j][2][4]+H[i][j][3][4];

			for ( k = 0; k < 4; k++ ) {	/* normieren */	
				for ( l = 0; l < 4; l++ )
					H[i][j][k][l] /= H[i][j][4][4];

				H[i][j][k][4] /= H[i][j][4][4];
				H[i][j][4][k] /= H[i][j][4][4];
			}
		}
	}

}


/******************************************************************************/
/*
void Write_Qij ( int a )
{

	FILE *fps;

	int i, k, l;

	if ( a == 0 )
	{
		fps = fopen ( ausgabe_0, "w" );
	}

	else if ( a == 1 )
	{
		fps = fopen ( ausgabe_1, "a" );
	}

	else return;


	for ( i = 0; i < taxa* ( taxa-1 ) /2; i++ )
	{
		for ( k = 0; k < 4; k++ )
		{
			for ( l = 0; l < 4; l++ )
			{

				if ( k != l )
					fprintf ( fps,"%f\t",q_matrizen[i][k*4+l] );
			}
		}

		fprintf ( fps,"\n" );
	}

	fclose ( fps );

}
*/


/******************************************************************************/


void Compute_Qij_tij()
{

	int e, i, j, k, l, index_paar;

	double /*P[16],*/ EigenWert[4], T1[16], U[16], V[16], T2[16];
	DMat20 HMat, EigenVec, EigenVecInv;
	DVec20 PiVec;

	for ( i = 0; i < taxa; i++ )
		distance[i][i] = 0.0;

	for ( i = 0; i < taxa-1; i++ ) {
		for ( j = i+1; j < taxa; j++ ) {
			distance[i][j] = 100;
			distance[j][i] = 100;
			index_paar = ( int ) ( i* ( taxa- ( i+3. ) /2. ) +j-1 );
			for ( k = 0; k < 16; k++ )
				q_matrizen[index_paar][k] = 0;

			for ( k = 0; k < 4; k++ ) {
				PiVec[k] = H[i][j][k][4];
				for ( l = 0; l < 4; l++ ) {		/* P_ij(t) = Pi^(-1)*H_ij */
					/*P[4*k+l] = H[i][j][k][l] / H[i][j][k][4];*/
					HMat[k][l] = H[i][j][k][l];
				}
			}

			/*if ( ( e=eigen ( 1, P, 4, EigenWert, T1, U, V, T2 ) ) !=0 )*/
			if ( ( e=eigen_sym (HMat, PiVec, 4, EigenWert, EigenVec, EigenVecInv ) ) !=0 ) 	{
				fprintf ( stderr, "\ncomplex roots in Eigen\n" );
				return;
			} 
			if ( EigenWert[0] <= 0.0001 || EigenWert[1] <= 0.0001 || EigenWert[2] <= 0.0001 ||  EigenWert[3] <= 0.0001 ||
					EigenWert[0] > 1.01 || EigenWert[1] > 1.01 || EigenWert[2] > 1.01 ||  EigenWert[3] > 1.01 )
			{/*
				fprintf ( stderr, "\nbad numerics in estimation of Eigenvalues (%f, %f, %f, %f) of P(t) %d,%d\n",EigenWert[0], EigenWert[1], EigenWert[2], EigenWert[3],i+1,j+1 );
				fprintf(stderr, "H = %f %f %f %f\n", H[i][j][0][0], H[i][j][0][1], H[i][j][0][2], H[i][j][0][3]);
				fprintf(stderr, "    %f %f %f %f\n", H[i][j][1][0], H[i][j][1][1], H[i][j][1][2], H[i][j][1][3]);
				fprintf(stderr, "    %f %f %f %f\n", H[i][j][2][0], H[i][j][2][1], H[i][j][2][2], H[i][j][2][3]);
				fprintf(stderr, "    %f %f %f %f\n", H[i][j][3][0], H[i][j][3][1], H[i][j][3][2], H[i][j][3][3]);
				fprintf(stderr, "Pt= %f %f %f %f\n", P[0], P[1], P[2], P[3]);
				fprintf(stderr, "    %f %f %f %f\n", P[4], P[5], P[6], P[7]);
				fprintf(stderr, "    %f %f %f %f\n", P[8], P[9], P[10], P[11]);
				fprintf(stderr, "    %f %f %f %f\n", P[12], P[13], P[14], P[15]);
				fprintf(stderr, "Pi= %f %f %f %f\n", H[i][j][0][4], H[i][j][1][4], H[i][j][2][4], H[i][j][3][4]);
				*/
				continue;
			}

			for ( k = 0; k < 4; k++ )	{
				for ( l = 0; l < 4; l++ ) {
					U[k*4+l] = EigenVec[k][l];
					V[k*4+l] = EigenVecInv[k][l];
				}
			}

			/*xtoy ( U, V, 16 );
			matinv ( V, 4, 4, T1 );*/

			/* berechne ratenmatrix  */

			/* T2 = diag(eigenwert)*/

			for ( k = 0; k < 16; k++ )
				T2[k] = 0;

			if ( alpha > 10 )	{
				/* keine ratenheterogenitaet */

				T2[0] = log ( EigenWert[0] );
				T2[5] = log ( EigenWert[1] );
				T2[10] = log ( EigenWert[2] );
				T2[15] = log ( EigenWert[3] );
			} else 	{
				/* ratenheterogenitaet */

				T2[0] = alpha* ( 1.-pow ( EigenWert[0],-1./alpha ) );
				T2[5] = alpha* ( 1.-pow ( EigenWert[1],-1./alpha ) );
				T2[10] = alpha* ( 1.-pow ( EigenWert[2],-1./alpha ) );
				T2[15] = alpha* ( 1.-pow ( EigenWert[3],-1./alpha ) );
			}

			/* T2 = U * diag(eigenwert) * V */

			matAbyBisC ( U, T2, 4, T1 );
			matAbyBisC ( T1, V, 4, T2 );

			/* normalisieren durch t = -sum (pi_i * r_ii) , Q hat rate 1 */

			distance[i][j] = 0;
			for ( k = 0; k < 4; k++ )
				distance[i][j] -= H[i][j][k][4]*T2[k*5];

			/* fix ZERO distance */
			if (fabs(distance[i][j]) < 0.00001) {
				if (distance[i][j] >= 0.0)
					distance[i][j] = 0.00001;
				else
					distance[i][j] = -0.00001;
			}

			distance[j][i] = distance[i][j];

			/* fix TOO LARGE distance */
			if (distance[i][j] > 100) {
				/*printf("Distance saturated, please try higher alpha\n");*/
				continue;
			}



			for ( k = 0; k < 16; k++ )
				q_matrizen[index_paar][k] = T2[k]/distance[i][j];

			/*

			int neg_rate = 0;
			for ( k = 0; k < 16; k++ )
				if ((k %5 != 0 && q_matrizen[index_paar][k] < 0.0)) neg_rate = 1;
			if (neg_rate) {
				printf("Negative non-diagonal entry of Q %d,%d\n", i+1, j+1);
				for ( k = 0; k < 16; k++ ) {
					printf("%+f ", q_matrizen[index_paar][k]);
					if (k % 4 == 3) printf("\n");
				}
				
			}*/

			
		}
	}



}
/*********************************************************************/

void Write_Tree (FILE *fp1, knoten *P )
{
	knoten *Q;


	if ( P->left != 0 && P->right != 0 ) {
		fprintf ( fp1,"(" );

		Q = P->left;
		Write_Tree (fp1, Q );
		fprintf ( fp1,"," );

		Q = P->right;
		Write_Tree (fp1, Q );

		/*if ( P->edge_length>0 )*/
			fprintf ( fp1,"):%f",P->edge_length );
		/*else
			fprintf ( fp1,")" );*/
	} else if	( P->left == 0 && P->right == 0 ) {
		fprintf ( fp1,"%s:%f",P->bezeichnung,P->edge_length );
	} 

}

/************************************************************/

void Save_Tree ( knoten *P )
{

	/*knoten *Q;*/
	FILE *fp1;

	if ( ( fp1=fopen ( ausgabe_report,"a" ) ) == 0 ) {
		fprintf ( stderr,"\nERROR writing file %s\n", ausgabe_report );
		Finalize ( 1 );
	}

	fprintf(fp1, "\nNEIGHBOR-JOINING TREE\n\n");

	fprintf ( fp1,"(" );

	Write_Tree (fp1, P->left );

	fprintf ( fp1,"," );
	Write_Tree (fp1, P->right);

	fprintf ( fp1,")" );


	fprintf ( fp1,";\n\n" );

	fclose ( fp1 );

	if ( ( fp1=fopen ( ausgabe_nj_tree,"w" ) ) == 0 ) {
		fprintf ( stderr,"\nERROR writing file %s\n", ausgabe_nj_tree );
		Finalize ( 1 );
	}

	fprintf ( fp1,"(" );

	Write_Tree (fp1, P->left );

	fprintf ( fp1,"," );
	Write_Tree (fp1, P->right);

	fprintf ( fp1,")" );

	fprintf ( fp1,";\n" );


	fclose ( fp1 );
}


/*********************************************************************/
/******************************************************************************/

void ComputeNeighborJoiningTree()
{

	int i, j, c, p1 = 0, p2 = 0, nr_nodes;
	int *cluster_index;
	double **nj_matrix, nj_distance, *hilfsvektor;
	double current_minimum, max = 0;
	/*knoten *P;*/


	cluster_index = ( int * ) malloc ( taxa * sizeof ( int ) );
	for ( i = 0; i < taxa; i++ )
		cluster_index[i] = i;

	hilfsvektor = ( double * ) calloc ( taxa, sizeof ( double ) );

	/* initialize nj_matrix */

	nj_matrix = ( double ** ) malloc ( taxa * sizeof ( double * ) );

	for ( i = 0; i < taxa; i++ ) {
		nj_matrix[i] = ( double * ) calloc ( taxa+1, sizeof ( double ) );
		nj_matrix[i][i] = 0.0;
	}

	for ( i = 0; i < taxa-1; i++ )	{
		for ( j = i+1; j < taxa; j++ )
			nj_matrix[i][j] = nj_matrix[j][i] = distance[i][j];
	}


	/* initialize tree */

	/*  	baum = (knoten *) malloc ( (2*taxa-1) * sizeof ( knoten ) ); schon oben passiert  */

	for ( i = 0; i < 2*taxa-1; i++ )	{
		baum[i].label = baum[i].ixlabel = i;
		baum[i].left = NULL;
		baum[i].right = NULL;
		baum[i].up = NULL;
		baum[i].edge_length = 0;
	}


	/*  build tree  */

	c = nr_nodes = taxa;

	while ( /*c < 2 * taxa - 1 &&*/ nr_nodes > 2 )	{
		for ( i = 0; i < nr_nodes; i++ )	{
			nj_matrix[i][taxa] = 0;

			for ( j = 0; j < nr_nodes; j++ )
				nj_matrix[i][taxa] += nj_matrix[i][j];
		}



		/*	find_minimum_distance();	*/

		current_minimum = max;

		for ( i = 0; i < nr_nodes - 1; i++ )	{
			for ( j = i + 1; j < nr_nodes; j++ )	{
				nj_distance = nj_matrix[i][j] -
				              ( nj_matrix[i][taxa] + nj_matrix[j][taxa] ) / ( nr_nodes-2. );

				if ( nj_distance < current_minimum ) {
					p1 = i;
					p2 = j;
					current_minimum = nj_distance;
				}
			}
		}


		/*	build_next_node();	*/

		baum[c].left = baum + cluster_index[p1];
		baum[c].right = baum + cluster_index[p2];

		baum[cluster_index[p1]].up = baum + c;
		baum[cluster_index[p2]].up = baum + c;

		baum[cluster_index[p1]].edge_length = 
			( nj_matrix[p1][p2] + (nj_matrix[p1][taxa] - nj_matrix[p2][taxa] ) / ( nr_nodes-2. ) ) /2.;


		if ( baum[cluster_index[p1]].edge_length < 0 )	{
			baum[cluster_index[p1]].edge_length = 0;
			baum[cluster_index[p2]].edge_length = nj_matrix[p1][p2];
		}	else
			baum[cluster_index[p2]].edge_length = nj_matrix[p1][p2] -
			                                      baum[cluster_index[p1]].edge_length;

		if ( baum[cluster_index[p2]].edge_length < 0 )	{
			baum[cluster_index[p2]].edge_length = 0;
			baum[cluster_index[p1]].edge_length = nj_matrix[p1][p2];
		}

		/*	update_nj_matrix();	*/
		for ( j = 0; j < nr_nodes; j++ ) {
			if ( j == p1 && j == p2 ) /* MINH: this condition never works since p1 != p2 ! */
				hilfsvektor[j] = 0;
			else
				hilfsvektor[j] = ( nj_matrix[p1][j] + nj_matrix[p2][j] - nj_matrix[p1][p2] ) / 2.;
		}

		for ( j = 0; j < nr_nodes; j++ )
			nj_matrix[j][p1] = nj_matrix[p1][j] = hilfsvektor[j];

		for ( j = 0; j < nr_nodes-1; j++ )	{
			nj_matrix[p2][j] = nj_matrix[nr_nodes-1][j];
			nj_matrix[j][p2] = nj_matrix[p2][j];
		}

		nj_matrix[p2][p2] = 0;

		for ( j = 0; j < nr_nodes-1; j++ )	{
			nj_matrix[nr_nodes-1][j] = nj_matrix[j][nr_nodes-1] = 0;
		}

		for ( j = 0; j < taxa; j++ )
			nj_matrix[j][taxa] = 0;

		cluster_index[p1] = c;
		cluster_index[p2] = cluster_index[nr_nodes-1];
		nr_nodes--;
		c++;

	}

	/*	verbinde zwei letzte knoten mit virtuellem Wurzelknoten*/

	p1 = 0; p2 = 1;
	baum[c].left = baum + cluster_index[p1];
	baum[c].right = baum + cluster_index[p2];

	baum[cluster_index[p1]].up = baum + c;
	baum[cluster_index[p2]].up = baum + c;

	baum[cluster_index[p1]].edge_length = nj_matrix[p1][p2]/2;
	baum[cluster_index[p2]].edge_length = nj_matrix[p1][p2]/2;

	/* release memory */
	for ( i = taxa-1; i >= 0; i-- )
		free(nj_matrix[i]);
	free(nj_matrix);
	free(hilfsvektor);
	free(cluster_index);

/*
	for (i = 0; i <= 2*(taxa-1); i++) {
		P = baum+i;
		if (P->left != 0 && P->right != 0)
			printf("%d -> (%d, %d)\n", P->label, P->left->label, P->right->label);
		else
			printf("%d\n", P->label);
	}*/
}

/******************************************************************************/

void Compute_q_hat_pairwise()
{

	int i, j, k, num_q;
	double Q[16], T1[16], T2[16];
	double rate[6];
	FILE *fps;
	double rate_sum;
	double *q0_matrix = ( double * ) calloc ( 16, sizeof ( double ) );

	num_q = 0;
	for ( i = 0; i <  (taxa-1)*taxa/2; i++ )
		if (q_matrizen[i][0] != 0.0) num_q ++;

/*	if (isMasterProc())
		printf("%d/%d valid Q matrices\n", num_q, (taxa-1)*taxa/2);*/
	
	for ( k = 0; k < 16; k++ ) {
		q0_matrix[k] = 0.0;

		for ( i = 0; i < ( int ) ( ( taxa-1. ) *taxa/2. ); i++ )
			q0_matrix[k] += q_matrizen[i][k];

		q0_matrix[k] /= num_q;
	}

	/* spektral zerlegung von q0_matrix  */

	for ( k = 0; k < 16; k++ )
		Q[k] = q0_matrix[k];

	if ( ( k=eigen ( 1, Q, 4, q_hat_eigen, T1, U_q_hat, V_q_hat, T2 ) ) !=0 ) {
		if (isMasterProc())
			fprintf ( stderr,"\nno spectral decomposition for q0_matrix\n" );
		free(q0_matrix);
		return;
	} else {
		if ( q_hat_eigen[0] > 0.01 || q_hat_eigen[1] > 0.01 || q_hat_eigen[2] > 0.01 ||  q_hat_eigen[3] > 0.01 ) {
			if (isMasterProc()) {
				fprintf ( stderr,"\n%f\t%f\t%f\t%f\n",q_hat_eigen[0],q_hat_eigen[1],q_hat_eigen[2],q_hat_eigen[3] );
				fprintf ( stderr, "\nbad numerics in estimation of Eigenvalues of NULL-Qmatrix\n" );
			}
			Finalize ( 1 );
		} else {
			xtoy ( U_q_hat, V_q_hat, 16 );
			matinv ( V_q_hat, 4, 4, T1 );
			for ( k = 0; k < 4; k++ )
				statPi[k] = V_q_hat[k]/ ( V_q_hat[0]+V_q_hat[1]+V_q_hat[2]+V_q_hat[3] );
			for ( k = 1; k < 4; k++ )
				statPi[k] += statPi[k-1];
		}
	}

	fps = fopen ( ausgabe_report, "a" );

	fprintf(fps, "\nSUBSTITUTION PROCESS OF HOMOGENEOUS MODEL\n\n");

	fprintf ( fps,"Q matrix:\n" );

	for ( k = 0; k < 16; k++ ) {
		if ( k%4 == 0 )
			fprintf ( fps,"\n" );

		fprintf ( fps,"%f\t",q0_matrix[k] );
	}

	fprintf ( fps,"\n\nBase composition:\n\n%f\t",statPi[0] );

	for ( k = 1; k < 4; k++ )
		fprintf ( fps,"%f\t",statPi[k]-statPi[k-1] );
	fprintf ( fps,"\n" );

	/* print individual rates */
	fprintf(fps, "\nRate:\n\n");
	for (i = 0, k = 0; i < 3; i++)
		for (j=i+1; j < 4; j++, k++) {
			rate[k] = q0_matrix[i*4+j] / (statPi[j] - statPi[j-1]);
		}
	for (k=0; k < 6; k++)
		fprintf(fps, "%f\n", rate[k]/rate[5]);
	fprintf ( fps,"\n" );

	fclose ( fps );

	/* check that the scaling is 1 total subst per site */
	rate_sum = 0.0;
	for (i = 0; i < 4; i++) {
		rate_sum -= q0_matrix[i*4+i] * ((i==0) ? statPi[0] : (statPi[i]-statPi[i-1]));
	}
	if (fabs(rate_sum - 1.0) > 1e-3) {
		if (isMasterProc())
			fprintf ( stderr,"\nq0_matrix not scaled to 1 total subst. per site (%f)\n", rate_sum );
		Finalize(1);
	}

	free(q0_matrix);
}




/******************************************************************************/

/******************************************************************************/
/*   subroutine to Simulate_Sequences_q_hat() **/

void EvolveSequences ( knoten *K, int **sim_sequences, double U[], double V[], double QEigenWert[], double *alpha_rate )
{
	int b, i, j, k, k1, l;
	double x, T1[16], T2[16], P_matrix[16];

	i = K->ixlabel;

	/* evolve left child */

	if ( K->left->edge_length > 0 )
	{
		if ( alpha > 10 )
		{
			for ( k = 0; k < 16; k++ )
				T2[k] = 0;

			T2[0] = exp ( QEigenWert[0]*K->left->edge_length );
			T2[5] = exp ( QEigenWert[1]*K->left->edge_length );
			T2[10] = exp ( QEigenWert[2]*K->left->edge_length );
			T2[15] = exp ( QEigenWert[3]*K->left->edge_length );

			matAbyBisC ( U, T2, 4, T1 );
			matAbyBisC ( T1, V, 4, P_matrix );

			for ( k = 0; k < 4; k++ ) {
				for ( k1 = 1; k1 < 4; k1++ )
					P_matrix[k*4+k1] += P_matrix[k*4+k1-1];
			}

			for ( j = 0; j < nr_basen; j++ ) {
				b = sim_sequences[i][j];
				l = 0;
				x = dkiss();

				while ( x > P_matrix[b*4+l] && l < 3)
					l++;

				sim_sequences[K->left->ixlabel][j] = l;
			}
		} else 	{

			for ( j = 0; j < nr_basen; j++ ) {
				for ( k = 0; k < 16; k++ )
					T2[k] = 0;

				T2[0] = exp ( QEigenWert[0]*K->left->edge_length*alpha_rate[j] );
				T2[5] = exp ( QEigenWert[1]*K->left->edge_length*alpha_rate[j] );
				T2[10] = exp ( QEigenWert[2]*K->left->edge_length*alpha_rate[j] );
				T2[15] = exp ( QEigenWert[3]*K->left->edge_length*alpha_rate[j] );

				matAbyBisC ( U, T2, 4, T1 );
				matAbyBisC ( T1, V, 4, P_matrix );

				for ( k = 0; k < 4; k++ ) {
					for ( k1 = 1; k1 < 4; k1++ )
						P_matrix[k*4+k1] += P_matrix[k*4+k1-1];
				}

				b = sim_sequences[i][j];
				l = 0;
				x = dkiss();

				while ( x > P_matrix[b*4+l] && l < 3)
					l++;

				sim_sequences[K->left->ixlabel][j] = l;
			}
		}
	} else {
		for ( j = 0; j < nr_basen; j++ )
			sim_sequences[K->left->ixlabel][j] = sim_sequences[i][j];
	}

	if ( K->left->ixlabel > taxa-1 ) 		/*  K ist kein blatt */
		EvolveSequences ( K->left, sim_sequences, U, V, QEigenWert, alpha_rate );



	/* evolve right child */


	if ( K->right->edge_length > 0 ) {
		if ( alpha > 10 ) {

			for ( k = 0; k < 16; k++ )
				T2[k] = 0;

			T2[0] = exp ( QEigenWert[0]*K->right->edge_length );
			T2[5] = exp ( QEigenWert[1]*K->right->edge_length );
			T2[10] = exp ( QEigenWert[2]*K->right->edge_length );
			T2[15] = exp ( QEigenWert[3]*K->right->edge_length );

			matAbyBisC ( U, T2, 4, T1 );
			matAbyBisC ( T1, V, 4, P_matrix );

			for ( k = 0; k < 4; k++ ) {
				for ( k1 = 1; k1 < 4; k1++ )
					P_matrix[k*4+k1] += P_matrix[k*4+k1-1];
			}

			for ( j = 0; j < nr_basen; j++ ) {
				b = sim_sequences[i][j];
				l = 0;
				x = dkiss();

				while ( x > P_matrix[b*4+l] && l < 3)
					l++;

				sim_sequences[K->right->ixlabel][j] = l;
			}

		} else {
			for ( j = 0; j < nr_basen; j++ ) {
				for ( k = 0; k < 16; k++ )
					T2[k] = 0;

				T2[0] = exp ( QEigenWert[0]*K->right->edge_length*alpha_rate[j] );
				T2[5] = exp ( QEigenWert[1]*K->right->edge_length*alpha_rate[j] );
				T2[10] = exp ( QEigenWert[2]*K->right->edge_length*alpha_rate[j] );
				T2[15] = exp ( QEigenWert[3]*K->right->edge_length*alpha_rate[j] );

				matAbyBisC ( U, T2, 4, T1 );
				matAbyBisC ( T1, V, 4, P_matrix );

				for ( k = 0; k < 4; k++ ) {
					for ( k1 = 1; k1 < 4; k1++ )
						P_matrix[k*4+k1] += P_matrix[k*4+k1-1];
				}

				b = sim_sequences[i][j];

				l = 0;
				x = dkiss();

				while ( x > P_matrix[b*4+l] && l < 3)
					l++;

				sim_sequences[K->right->ixlabel][j] = l;
			}
		}
	} else {
		for ( j = 0; j < nr_basen; j++ )
			sim_sequences[K->right->ixlabel][j] = sim_sequences[i][j];
	}

	if ( K->right->ixlabel > taxa-1 ) 		/*  K->right ist kein blatt */
		EvolveSequences ( K->right, sim_sequences, U, V, QEigenWert, alpha_rate );

}



/******************************************************************************/


void Simulate_Sequences_q_hat()
{

	int j, l;
	double x;

	if ( alpha > 10 )	{
		/*homogene raten */
		;
	} else {
		/* heterogene raten zuweisen */
		for ( j = 0; j < nr_basen; j++ ) {
			alpha_rate[j] = rgamma ( alpha,beta );
		}
	}

	/* generate sequences along the tree ************/
	/*	root sequence	baum[2*taxa-2] */

	for ( j = 0; j < nr_basen; j++ ) {
		l = 0;
		x = dkiss();

		while ( x > statPi[l] && l < 3)
			l++;

		seqData[2*taxa-2][j] = l;
	}

	EvolveSequences ( baum+ ( 2*taxa-2 ), seqData, U_q_hat, V_q_hat, q_hat_eigen, alpha_rate );

	/* end of generate sequences along the tree ************/

}
