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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*#include <iqtree_config.h>*/
#include "utils/timeutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef CLANG_UNDER_VS
	//These are "safe" in the sense that they won't overrun the buffer.  
	//But they will bomb if the buffer isn't big enough
	#define safe_strcpy(dest, source) strcpy_s(&dest[0], sizeof(dest), source)
	#define safe_strcat(dest, source) \
		(strlen(dest)<sizeof(dest)) \
			? strcpy_s(&dest[0], sizeof(dest)-strlen(dest), source) \
			: strcpy_s((&dest[0])+strlen(dest), 0UL, source)
#else
#define safe_strcpy(dest, source) strcpy(dest, source)
#define safe_strcat(dest, source) strcat(dest, source)
#endif
#include <time.h>
#include <math.h>
#include "weisslambda_sub.h"

#include "whtest_sub.h"
#include "random.h"
#include "whtools.h"
/*
#ifdef WIN32
#include <sys/timeb.h>
#include <sys/types.h>
#include <winsock.h>
void gettimeofday(struct timeval* t, void* timezone)
{       struct _timeb timebuffer;
        _ftime( &timebuffer );
        t->tv_sec=timebuffer.time;
        t->tv_usec=1000*timebuffer.millitm;
}
#else
  #include <sys/time.h>
  #ifndef HAVE_GETTIMEOFDAY
	void gettimeofday(struct timeval* t, void* timezone) {
		time_t cur_time;
		time(&cur_time);
		t->tv_sec = cur_time;
		t->tv_usec = 0;
	}
  #endif
#endif
*/

#ifdef PARALLEL
int mpi_myrank;
int mpi_size;
int mpi_master_rank = 0;

long p_randn;
long p_rand;
#endif /*PARALLEL*/

int isMasterProc() {
#ifdef PARALLEL
	return mpi_myrank == mpi_master_rank;
#else
	return 1;
#endif
}

int isSlaveProc() {
#ifdef PARALLEL
	return mpi_myrank != mpi_master_rank;
#else
	return 0;
#endif
}

/*
int isFirstSlaveProc() {
#ifdef PARALLEL
	if (mpi_size == 1)
		return 1;
	return mpi_myrank == mpi_master_rank+1;
#else
	return 1;
#endif
}*/


void Finalize(int exit_code) {
#ifdef PARALLEL
	MPI_Finalize();
#endif
	if (isMasterProc())
		if (exit_code == 0)
			printf("\nFinished successfully.\n");
	exit(exit_code);
}



int simulation, current_sim, nr_basen, taxa;

int random_seed = -1;
int check_times = 10;
double p_value_cutoff;
double alpha, beta;
/*int verbose_mode = 0;*/
int write_sim_result = 0;
int write_dist_matrix = 0;
int fix_distance = 0;

double delta_data;
double *delta_sim;
double p_wert;

char datei_name[100];

/*
char ausgabe_0[200];
char ausgabe_1[200];
char ausgabe_2[200];
*/
char ausgabe_report[200];
char ausgabe_sim_result[200];
char ausgabe_dist[200];
char ausgabe_nj_tree[200];

double *ml_distance = NULL;

void WHT_setAlignmentSize(int ntax, int nsite) {
	taxa = ntax;
	nr_basen = nsite;
}

void WHT_allocateMemory() {
	AllocateMemory();
}

void WHT_setSequenceSite(int seqid, int siteid, char c) {
	if (c>4) c = 4;
	seqData[seqid][siteid] = c;
}

void WHT_setSequenceName(int seqid, const char *name) {
	safe_strcpy ( baum[seqid].bezeichnung, name);	
}

void WHT_setParams(int nsim, double gamma_shape, char *filename, double *dist) {
	simulation = nsim;
	alpha = gamma_shape;
	safe_strcpy(datei_name, filename);
	current_sim = 0;
	p_value_cutoff = 1.0;

	safe_strcpy ( ausgabe_report, datei_name );
	safe_strcat ( ausgabe_report, ".whtest" );

	safe_strcpy ( ausgabe_sim_result, datei_name );
	safe_strcat ( ausgabe_sim_result, ".whsim" );
	safe_strcpy ( ausgabe_dist, datei_name );
	safe_strcat ( ausgabe_dist, ".whdist" );

	safe_strcpy ( ausgabe_nj_tree, datei_name );
	safe_strcat ( ausgabe_nj_tree, ".nj" );

	ml_distance = dist;
	write_dist_matrix = 1;
	write_sim_result = 1;

}

void WHT_getResults(double *delta, double *delta_quantile, double *p_value) {
	*delta = delta_data;
	*delta_quantile = delta_sim[(int)floor(0.95*simulation)];
	*p_value = p_wert;
}


void SetMLDistance() {
	int i;
	for (i=0; i < taxa; i++)
		memcpy(distance[i], ml_distance + (i*taxa), sizeof(double)*taxa);
}

void usage(char *prog_name) {
	if (!isMasterProc()) Finalize(1);
	printf("Usage: %s <alignment> [OPTIONS]\n", prog_name);
	printf("  <alignment>         alignment file name, in standard PHYLIP format\n");
	printf("OPTIONS:\n");
	printf("  -h                  print usage\n");
	printf("  -s <SIMULATION>     #simulations to assess significance, default is 1000\n");
	printf("  -a <ALPHA>          gamma shape parameter, default is 100 (equal site-rates)\n");
	printf("  -t <CUTOFF>         stop the simulations when p-value exceeds the cutoff\n");
	printf("  -i <N>              check p-value N times during simulation, default 10\n");
	printf("  -seed <#>           use <#> as random number seed\n");
	printf("  -wsim               write simulation results to file .whtest.sim\n");
	printf("  -wdist              write distance matrix to file .whtest.dist\n");
	printf("\n");
	Finalize(1);
}

void parseArg( int argc,char **argv ) {
	int i;
	int arg_i;
	/*char *alpha_arg = NULL;*/

	if (isMasterProc()) {

		printf("\nWELCOME TO WH-TEST\n\n");
		printf("G. Weiss and A. von Haeseler (2003) Testing substitution models\n");
		printf("within a phylogenetic tree. Mol. Biol. Evol, 20(4):572-578\n\n");
	
#ifdef PARALLEL
		printf("You are running MPI parallel version with %d processes\n\n", mpi_size);
#endif


		printf("Program was called with:\n");
		for ( i = 0; i < argc; i++ )
			printf ( "%s ",argv[i] );
		printf ( "\n\n" );
	}


	simulation = 1000;
	current_sim = 0;
	alpha = 100;
	datei_name[0] = 0;
	p_value_cutoff = 1.0;

	for (arg_i = 1; arg_i < argc; arg_i++) {
		if (strcmp(argv[arg_i], "-h") == 0) {
			usage(argv[0]);
		} else if (strcmp(argv[arg_i], "-s") == 0) {
			arg_i++;
			simulation = atoi ( argv[arg_i] );
		} else if (strcmp(argv[arg_i], "-t") == 0) {
			arg_i++;
			p_value_cutoff = atof ( argv[arg_i] );
		} else if (strcmp(argv[arg_i], "-a") == 0) {
			arg_i++;
			/*alpha_arg = argv[arg_i];*/
			alpha = atof ( argv[arg_i] );
		} else if (strcmp(argv[arg_i], "-seed") == 0) {
			arg_i++;
			random_seed = atoi ( argv[arg_i] );
		} else if (strcmp(argv[arg_i], "-i") == 0) {
			arg_i++;
			check_times = atoi ( argv[arg_i] );
		} else if (strcmp(argv[arg_i], "-v") == 0) {
			/*verbose_mode = 1;*/
		} else if (strcmp(argv[arg_i], "-wsim") == 0) {
			write_sim_result = 1;
		} else if (strcmp(argv[arg_i], "-wdist") == 0) {
			write_dist_matrix = 1;
		} else if (strcmp(argv[arg_i], "-fdist") == 0) {
			fix_distance = 1;
		} else if (argv[arg_i][0] != '-') {
			safe_strcpy ( datei_name, argv[arg_i] );
		
			safe_strcpy ( ausgabe_report, datei_name );
			safe_strcat ( ausgabe_report, ".whtest" );

			safe_strcpy ( ausgabe_sim_result, ausgabe_report );
			safe_strcat ( ausgabe_sim_result, ".sim" );
			safe_strcpy ( ausgabe_dist, ausgabe_report );
			safe_strcat ( ausgabe_dist, ".dist" );
			
		} else {
			if (isMasterProc()) {
				printf("Unrecognized %s option, run with '-h' for help\n", argv[arg_i]);
			}
			Finalize(1);
		} 
	}

	if (datei_name[0] == 0) {
		printf("ERROR: Missing input alignment file.\n\n");
		usage(argv[0]);
	}

	if (simulation <= 0 || simulation > 10000) {
		if (isMasterProc())
			fprintf ( stderr,"wrong #simulations: %d\nbetween 1 and 10000 please\n", simulation);
		Finalize( 1 );
	}

	if (alpha < 0.01 || alpha > 100.0) {
		if (isMasterProc())
			fprintf ( stderr,"wrong alpha: %f\nbetween 0.01 and 100 please\n", alpha);
		Finalize ( 1 );
	}

	if (check_times < 0) {
		if (isMasterProc())
			fprintf ( stderr,"wrong time interval: %d\npositive number please\n", check_times);
		Finalize(1);
	}

	if (isMasterProc()) {
		printf("Input file: %s\n", datei_name);
		printf("Number of simulations: %d\n", simulation);
		printf("Gamma shape alpha: %f\n", alpha);
	}

}

void StartReport() {
	FILE *fps = fopen( ausgabe_report, "w" );
	fprintf(fps, "WH-TEST\n\n");
	fprintf(fps, "G. Weiss and A. von Haeseler (2003) Testing substitution models\n");
	fprintf(fps, "within a phylogenetic tree. Mol. Biol. Evol, 20(4):572-578\n\n");
	fprintf(fps, "Input file name: %s\n", datei_name);
	fprintf(fps, "Number of simulations: %d\n", simulation);
	fprintf(fps, "Gamma shape parameter: %f\n", alpha);
	fprintf(fps, "Random number seed: %d\n\n", random_seed);
	fprintf(fps, "SEQUENCE ALIGNMENT\n\n");
	fprintf(fps, "Input data: %d sequences with %d nucleotide sites\n", taxa, nr_basen);
	fprintf(fps, "\n");
	fclose(fps);
}

void FinishReport(time_t begin_time) {
	FILE *fps = fopen( ausgabe_report, "a" );
	char *finishedDate_;
	int prog_time;
	int nHour_, nMin_, nSec_;

	time_t end_time;
	time(&end_time);
	finishedDate_ = ctime(&end_time);

	prog_time = difftime (end_time, begin_time);

	nHour_ = prog_time / 3600;
	nMin_ = (prog_time  - nHour_ * 3600) / 60;
	nSec_ = prog_time  - nMin_ * 60 - nHour_ * 3600;

	/*printf("\nDate and time: %s", finishedDate_);*/
	printf("Runtime: %dh:%dm:%ds\n\n", nHour_, nMin_, nSec_);

	fprintf(fps, "\nTIME STAMP\n\n");
	fprintf(fps, "Date and time: %s", finishedDate_);
	fprintf(fps, "Runtime: %dh:%dm:%ds\n", nHour_, nMin_, nSec_);

	fclose(fps);
}

void ReportResults(double delta_data, double delta_95quantile, double p_value) {

	FILE *fps = fopen( ausgabe_report, "a" );

	fprintf(fps, "\nTEST OF HOMOGENEITY ASSUMPTION OVER BRANCHES\n\n");

	fprintf(fps, "Delta of data:                       %f\n", delta_data);
	fprintf(fps, ".95 quantile of Delta distribution:  %f\n", delta_95quantile);
	fprintf(fps, "Number of simulations performed:     %d\n", current_sim);
	if (current_sim == simulation)
		fprintf(fps, "p-value:                             %f\n", p_value);
	else
		fprintf(fps, "p-value:                             >%f\n", p_value);

	fprintf(fps, "\n");
	if (p_value < 0.05) {
		fprintf(fps, "WH-test rejected the assumption of a single model among branches of the tree\n");
	} else {
		fprintf(fps, "WH-test DID NOT reject the assumption of a single model among branches of the tree\n");
	}

	fclose(fps);
}

int WHTest_run ( int argc,char **argv ) {

	int i;
	/*double *global_sim = NULL;*/
	int count_sim;
	int cur_point;

	double prev_p_wert = 0, own_p_wert;
	int *valid_pairs;
	/*int *global_pairs = NULL;*/
	FILE *delta_file = NULL;
	time_t begin_time;
	struct timeval tv;
	int work_single;
#ifdef PARALLEL
	int *displs, *rcounts;
	double mpi_prog_time, mpi_sim_time;
#endif

	int start_sim, end_sim;

	p_wert = 0.0;

	/*knoten *baum;*/

#ifdef PARALLEL
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

	/* start to count the running time */
	time(&begin_time);

	if (argc>0) parseArg(argc, argv);
	/* initialize random seed based on current time */
	if (isMasterProc()) {
#ifndef HAVE_GETTIMEOFDAY
		if (random_seed < 0) {
			printf("WARNING: Random seed may not be well initialized since gettimeofday() is not available.\n");
			printf("         You can use option -seed <NUMBER> to specify your own seed number.\n");
		}
#endif
		gettimeofday(&tv, NULL);
		srand((unsigned) (tv.tv_sec+tv.tv_usec));
		if (random_seed < 0)
			random_seed = rand();
		if (argc > 0)
		printf("Random number seed: %d\n\n", random_seed);
	}

#ifdef PARALLEL
	MPI_Bcast(&random_seed, 1, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
#endif

	start_kiss ( random_seed );

	beta = 1./alpha;

	if (argc > 0) {
		ReadDataSize ( datei_name );
		AllocateMemory();
	}
	delta_sim = ( double* ) calloc ( simulation, sizeof ( double) );
	valid_pairs = ( int* ) calloc ( simulation, sizeof ( int) );
	/*global_sim = ( double* ) calloc ( simulation, sizeof ( double) );
	global_pairs = ( int* ) calloc ( simulation, sizeof ( int) );*/
#ifdef PARALLEL
	displs = (int*) malloc(mpi_size * sizeof(int));
	rcounts = (int*) malloc(mpi_size * sizeof(int));
#endif

	if (isMasterProc() && argc > 0)
		printf("Input data set (%s) contains %d sequences of length %d\n", datei_name, taxa, nr_basen);
	if (argc > 0) ReadData ( datei_name );


	if (isMasterProc())
		printf("\n");
	if (isMasterProc())
		StartReport();
	
#ifdef PARALLEL
	mpi_prog_time = MPI_Wtime();
#endif

	Compute_Hij();
	Compute_Qij_tij();

	/*if (isMasterProc())
		printf("Computing average of Q matrices\n");*/
	Compute_q_hat_pairwise();

	delta_data = ComputeWeissLambdaQ16(q_matrizen); 

	if (fix_distance) 
		FixDistance();
	if (isMasterProc() && write_dist_matrix)
		Save_Distance(ausgabe_dist, distance);

	if (ml_distance) SetMLDistance();

	if (isMasterProc())
		printf("Computing neighbor-joining tree\n");


	ComputeNeighborJoiningTree();

	if (isMasterProc()) {
		Save_Tree ( baum + ( 2*taxa-2 ) );
		printf("\nStart %d simulations\n", simulation);
	}

#ifdef PARALLEL
	mpi_sim_time = MPI_Wtime();
	work_single = (simulation+mpi_size-1) / mpi_size;
	start_sim = work_single * mpi_myrank;
	end_sim = work_single * (mpi_myrank+1);
	if (end_sim > simulation) end_sim = simulation;
	work_single = end_sim - start_sim;
	for (i = 0; i < mpi_size; i++) {
		displs[i] = work_single * i;
		rcounts[i] = work_single;
		if (i == mpi_size-1) rcounts[i] = simulation - displs[i];
		/*if (isMasterProc())
			printf(" %d ", rcounts[i]);*/
	}
#else
	work_single = simulation;
	start_sim = 0;
	end_sim = simulation;
#endif

	int* check_point = NULL;
	if (check_times > 0) {
		check_point = (int*)malloc(check_times * sizeof(int));
		for (i = 0; i < check_times; i++) {
			check_point[i] = work_single * (i + 1) / check_times;
			if (i == check_times - 1)
				check_point[i] = end_sim - start_sim;
		}
	}

	for ( i = start_sim, count_sim = 0, own_p_wert = 0.0, cur_point = 0; i < end_sim; i++) {
		Simulate_Sequences_q_hat();
		Compute_Hij();
		Compute_Qij_tij();
		delta_sim[i] = ComputeWeissLambdaQ16(q_matrizen);
		valid_pairs[i] = CountValidPairs(q_matrizen);
		count_sim++;
		current_sim = count_sim;
		if (delta_sim[i] >= delta_data) own_p_wert += 1.0;
		p_wert = own_p_wert / simulation;
		if (check_point && count_sim == check_point[cur_point]) {
			cur_point++;
#ifdef PARALLEL
			MPI_Allreduce(&own_p_wert, &p_wert, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			p_wert /= simulation;
			MPI_Reduce(&count_sim, &current_sim, 1, MPI_INT, MPI_SUM, mpi_master_rank, MPI_COMM_WORLD);
#endif
			if (isMasterProc()) {
				printf("%5d done", current_sim);
				printf(", current p-value: %5.3f\n", p_wert);
				if (p_wert > 0.05 && prev_p_wert <= 0.05) {
					printf("NOTE: Homogeneity assumption is NOT rejected (p-value > 0.05)\n");
				}
				prev_p_wert = p_wert;
			}
		}
		if (p_wert > p_value_cutoff) 
			break;
	}


#ifdef PARALLEL
	/*printf("Proc %d done.\n", mpi_myrank);
	MPI_Barrier(MPI_COMM_WORLD);*/
	if (mpi_size > 1) {
			MPI_Gatherv(delta_sim + start_sim, end_sim - start_sim, MPI_DOUBLE, 
				delta_sim, rcounts, displs, MPI_DOUBLE, mpi_master_rank, MPI_COMM_WORLD);
			MPI_Gatherv(valid_pairs + start_sim, end_sim - start_sim, MPI_INT, 
				valid_pairs, rcounts, displs, MPI_INT, mpi_master_rank, MPI_COMM_WORLD);
			for (i = 0, current_sim = 0, p_wert = 0.0; i < simulation; i++) {
				if (delta_sim[i] >= delta_data) p_wert += 1.0;
				if (delta_sim[i] != 0.0) current_sim++;
			}
			p_wert /= simulation;
/*
		} else {
			MPI_Reduce(&own_p_wert, &p_wert, 1, MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI_COMM_WORLD);
			p_wert /= simulation;
			MPI_Reduce(&count_sim, &current_sim, 1, MPI_INT, MPI_SUM, mpi_master_rank, MPI_COMM_WORLD);
		}*/
	}
	/*printf("Process %d did %d simulations\n", mpi_myrank, count_sim);*/
#endif

	if (isMasterProc()) {
		printf("%d simulations done\n", current_sim);
	}


	if (isMasterProc() && write_sim_result) {
		delta_file = fopen(ausgabe_sim_result, "w");
		if (!delta_file) {
			printf ( "\nERROR: Cannot write to file %s!\n", ausgabe_sim_result );
		} else {
			fprintf(delta_file, "Sim.    Delta   Valid_Qs\n");
			for (i = 0, count_sim = 1; i < simulation; i++)
				if (delta_sim[i] != 0.0) {
					fprintf(delta_file, "%d\t%f\t%d\n", count_sim++, delta_sim[i], valid_pairs[i]);
				}
			fclose(delta_file);
		}
	}


	if (isMasterProc()) {

#ifdef PARALLEL
		/*if (verbose_mode) {
			printf("Simulation time: %f\n", MPI_Wtime() - mpi_sim_time);
		}*/
#endif

		sort ( simulation, delta_sim-1);
		printf("\nDelta of input data: %f\n", delta_data);
		printf("0.95 quantile:       %f\n", delta_sim[(int)floor(0.95*simulation)]);

		if (current_sim == simulation)
			printf("P-value:             %f\n\n",p_wert);
		else
			printf("P-value:            >%f\n\n",p_wert);

	if (p_wert < 0.05) {
		printf("RESULT: Model homogeneity is rejected (p-value cutoff 0.05)\n");
	} else {
		printf("RESULT: Model homogeneity is NOT rejected (p-value cutoff 0.05)\n");
	}

		ReportResults(delta_data, delta_sim[(int)floor(0.95*simulation)], p_wert);
	if (argc > 0) {
		printf("All results written to disk:\n");
		printf("     WH-test report file:     %s\n", ausgabe_report);
		if (write_sim_result)
			printf("     Simulation results:      %s\n", ausgabe_sim_result);
		if (write_dist_matrix)
			printf("     Pairwise distances:      %s\n", ausgabe_dist);
	}

		FinishReport(begin_time);
#ifdef PARALLEL	
		/*if (verbose_mode) {
			printf("Total time: %f\n", MPI_Wtime() - mpi_prog_time);
		}*/
#endif
	}

#ifdef PARALLEL
	free(rcounts);
	free(displs);
#endif
	if (check_point) free(check_point);
	free(valid_pairs);
	free(delta_sim);
	FreeMemory();

#ifdef PARALLEL
	MPI_Finalize();
#endif
	if (isMasterProc() && argc > 0)
		printf("Finished successfully.\n");
	return 0;
}
