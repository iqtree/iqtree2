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
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "weisslambda_sub.h"
#include "eigen.h"
#include "eigen_sym.h"

/*
int simulation, nr_basen, taxa, paare;

char datei_name[15];
char ausgabe_null[20];
char ausgabe_data[20];
char lambdawerte[20];
*/

/*
char *ausgabe_null;
char *ausgabe_data;
int paare;

double **dataQ;

double **nullbetweenQ;

double WeissLambdaData;
*/

/******************************************************/


void Compute_SSbetween_Matrix ( double **data, int s, double SSbetween[] );

double ComputeWeissLambda ( double WeissMatrix[] );

/******************************************************/
/******************************************************/
#define NR_END 1
#define FREE_ARG char*

int *ivector ( long nl, long nh )
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v= ( int * ) malloc ( ( size_t ) ( ( nh-nl+1+NR_END ) *sizeof ( int ) ) );

	return v-nl+NR_END;
}
void free_ivector ( int *v, long nl, long nh )
/* free an int vector allocated with ivector() */
{
	free ( ( FREE_ARG ) ( v+nl-NR_END ) );
}

#define NRANSI

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort ( unsigned long n, double arr[] )
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	double a,temp;

	istack=ivector ( 1,NSTACK );
	for ( ;; )
	{
		if ( ir-l < M )
		{
			for ( j=l+1;j<=ir;j++ )
			{
				a=arr[j];
				for ( i=j-1;i>=1;i-- )
				{
					if ( arr[i] <= a ) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if ( jstack == 0 ) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		}
		else
		{
			k= ( l+ir ) >> 1;
			SWAP ( arr[k],arr[l+1] )
			if ( arr[l+1] > arr[ir] )
			{
				SWAP ( arr[l+1],arr[ir] )
			}
			if ( arr[l] > arr[ir] )
			{
				SWAP ( arr[l],arr[ir] )
			}
			if ( arr[l+1] > arr[l] )
			{
				SWAP ( arr[l+1],arr[l] )
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for ( ;; )
			{
				do i++; while ( arr[i] < a );
				do j--; while ( arr[j] > a );
				if ( j < i ) break;
				SWAP ( arr[i],arr[j] );
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;

			if ( ir-i+1 >= j-l )
			{
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			}
			else
			{
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector ( istack,1,NSTACK );
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

/******************************************************/
/******************************************************/



/*****************************************************************/

void printSSbetween (double SS[]) {
	int i, j;
	for (i = 0; i < 12; i++) {
		for ( j = 0; j < 12; j++)
			printf("%+7.5f ", SS[i*12+j]);
		printf("\n");
	}
}

int CountValidPairs(double **q_mats) {
	int i, num_q;

	num_q = 0;
	for ( i = 0; i < ( int ) ( ( taxa-1. ) *taxa/2. ); i++ )
		if (q_mats[i][0] != 0.0) num_q ++;

	return num_q;
}


/**
	compute the delta statistic from 16*16 matrix
*/
double ComputeWeissLambdaQ16(double **q_16) {
	double SSbetween[144];
	double delta;
	double **q_12;
	int i, j, pair_id;

	int paare = taxa * (taxa-1) / 2;

	q_12 = ( double ** ) malloc ( paare * sizeof ( double * ) );

	for ( i = 0; i < paare; i++ )
		q_12[i] = ( double * ) calloc ( 12, sizeof ( double ) );

	for (pair_id = 0; pair_id < paare; pair_id++)
		for (i = 0, j = 0; i < 16; i++)
			if (i % 5 != 0) q_12[pair_id][j++] = q_16[pair_id][i];

	/*
	for (pair_id = 0; pair_id < paare; pair_id++) {
		for (i = 0; i < 12; i++) 
			printf("%f ", q_16[pair_id][i]);
		printf("\n");
	}*/


	Compute_SSbetween_Matrix ( q_12, 0, SSbetween );

	/*printSSbetween(SSbetween);*/

	


	for ( i = paare-1; i >= 0; i-- )
		free(q_12[i]);
	free(q_12);

	delta = ComputeWeissLambda(SSbetween);

	return delta;
}

/***************************************************************/

/***************************************************************/






/***************************************************************/

void Compute_SSbetween_Matrix ( double **data, int s, double SSbetween[] )
{

	int i, k, l;
	int paare = taxa*(taxa-1)/2;

	double mean[12];

	int true_pair = 0;

	for ( k = 0; k < 12; k++ )
		mean[k] = 0;

	for ( k = 0; k < 144; k++ )
		SSbetween[k] = 0;

	for ( i = 0; i < paare; i++ )
	if (data[s*paare+i][0] != 0.0)
	{
		true_pair++;
		for ( k = 0; k < 12; k++ )
		{
			mean[k] += data[s*paare+i][k];

			for ( l = 0; l < 12; l++ )
				SSbetween[k*12+l] += data[s*paare+i][k]*data[s*paare+i][l];
		}
	} else { 
		/*fprintf(stderr, "one pair discarded\n");*/
	}

	for ( k = 0; k < 12; k++ )
		mean[k] /= ( double ) true_pair;

	for ( k = 0; k < 12; k++ )
	{
		for ( l = 0; l < 12; l++ )
			SSbetween[k*12+l] = SSbetween[k*12+l] - true_pair * mean[k] * mean[l];

	}

	for ( k = 0; k < 144; k++ )
		SSbetween[k] /= ( true_pair-1. );



	/*for( k = 0; k < 144; k++)
	SSbetween[k] *= simulation;
	*/
}


/***************************************************************/

double ComputeWeissLambda ( double WeissMatrix[] )
{

	double EigenWert[12], W[144]; /*T1[12], U[144], V[144], T2[144];*/

	int k;
/*
	double lambda_summe = 0, log_lambda_sum = 0;
	double product_lambda = 1, productinvlambda = 1;*/
	double product_log_lambda = 1;


	for ( k = 0; k < 144; k++ )
		W[k] = WeissMatrix[k];



	if ( ( k=eigen_sym_core ( W, 12, EigenWert ) ) !=0 )
	/*if ( ( k=eigen ( 1, W, 12, EigenWert, T1, U, V, T2 ) ) !=0 )*/
	{
		fprintf ( stderr, "\ncomplex roots in WilksMatrix\n" );

		return 0;
	}

	else
	{
		if ( EigenWert[0] > 100000 || EigenWert[11] < -0.1 )
		{
			fprintf ( stderr, "\nnumerical problems in eigenvalues of WeissMatrix\n" );
			return 0;
		}

		else
		{

			for ( k = 0; k < 12; k++ )
			{
				product_log_lambda += log ( 1.+ EigenWert[k] );

				/*
				lambda_summe += EigenWert[k];

				log_lambda_sum += log ( EigenWert[0] + EigenWert[k] );

				product_lambda *= ( 1.+ EigenWert[k] );


				productinvlambda /= ( 1.+ EigenWert[k] );*/

			}


			/*		fps = fopen(lambdawerte,"a");

			fprintf(fps,"%f\n",product_log_lambda);

					fclose( fps );

			*/
			return product_log_lambda;

		}
	}

}


/***************************************************************/


