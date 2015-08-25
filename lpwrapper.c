
/*#include "lp/lp_lib.h"*/
#include "lpwrapper.h"

/*
void __WINAPI msgfunction(lprec *lp, void *userhandle, int msg)
{
	switch(msg) {
	case MSG_LPFEASIBLE:
		printf("Feasible solution found\n");
		break;
	case MSG_MILPFEASIBLE:
		printf("Integer feasible solution found\n");
		break;
	case MSG_MILPBETTER:
		printf("Better integer feasible solution found\n");
		break;
	case MSG_MILPEQUAL:
		printf("Equal MILP solution found\n");
		break;
	}
}
*/

void lp_solve_version_info(int *majorversion, int *minorversion, int *release, int *build) {
	/*lp_solve_version(majorversion, minorversion, release, build);*/
}

int lp_solve(char *filename, int ntaxa, double *score, double *variables, int verbose_mode) {
	return 5;
/*	lprec *lp = NULL;
	int ret;
	int Ncol;
	int index, j;
	double *row = NULL;
	char *name;
	char name2[200];

	lp = read_LP(filename, IMPORTANT + ((verbose_mode < 2) ? 0 : verbose_mode-1), "pd");
	//lp = read_LP(filename, NORMAL, "pd");
	//strcpy(name2, filename);
	//strcat(name2,".cnv");
	//write_lp(lp, name2);

	if (lp == NULL) {
		printf("Could not create an LP_SOLVE getInstance!\n");
		return 1;
	}

	set_mip_gap(lp, TRUE, 0.0);
	
	ret = solve(lp);
	

    if(ret == OPTIMAL || ret == PRESOLVED) {
		ret = 0;
	} else {
		ret = 5;
		printf("LP_SOLVE ERROR: %s\n", get_statustext(lp, ret));	
		exit (1);
	}

	if(ret == 0) {
	// a solution is calculated, now lets get some results
	
	// objective value
	*score = get_objective(lp);
	// variable values 
	Ncol = get_Ncolumns(lp);
	row = (double*) malloc(Ncol * sizeof(*row)); 
	get_variables(lp, row);

	
	for(j = 0; j < Ncol; j++) {
		name = get_col_name(lp, j+1);
		if (name[0] == 'x') { // this is for taxa set
			index = -1;
			index = atoi(name+1);
			//printf(name);
			if (index < 0 || index >= ntaxa) {
				printf("Index x_%d is not in the range!\n", index);
				ret = 6;
				break;
			}
			if (row[j] > tolerance && (1.0 - row[j]) > tolerance) {
				if (verbose_mode >= 3) printf("\n%s = %10.8f", name, row[j]);
				ret = 7;
				if (verbose_mode < 3) break;
			}
			variables[index] = row[j];
		}
		//printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
	}
	free(row);
	
	// we are done now 
	}
	if(lp != NULL) {
		// clean up such that all used memory by lpsolve is freed
		delete_lp(lp);
	}

	return ret;*/
}
