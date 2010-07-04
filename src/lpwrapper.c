
#include "lp/lp_lib.h"
#include "lpwrapper.h"


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

int lp_solve(char *filename, int ntaxa, double *score, double *variables, int verbose_mode) {
	lprec *lp = NULL;
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
		printf("Could not create an LP_SOLVE instance!\n");
		return 1;
	}

	set_mip_gap(lp, TRUE, 0.0);
	strcpy(name2, filename);
	strcat(name2,".log");
	set_outputfile(lp, name2);
	set_verbose(lp, NORMAL);
	
	//if (verbose_mode) 
	//	set_verbose(lp, DETAILED);
	//put_msgfunc(lp, msgfunction, NULL, MSG_LPFEASIBLE | MSG_MILPFEASIBLE | MSG_MILPBETTER | MSG_MILPEQUAL);

	//set_presolve(lp, PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP, get_presolveloops(lp));
	ret = solve(lp);
	
	//strcpy(name2, filename);
	//strcat(name2,".presolve");
	//write_lp(lp, name2);
	
	//int spx_status = get_status(lp);
	//if (spx_status != OPTIMAL) {
	//}

    if(ret == OPTIMAL || ret == PRESOLVED) {
		ret = 0;
	} else {
		ret = 5;
		printf("LP_SOLVE ERROR: %s\n", get_statustext(lp, ret));	
		exit (1);
	}

	if(ret == 0) {
	/* a solution is calculated, now lets get some results */
	
	/* objective value */
	*score = get_objective(lp);
	/* variable values */
	Ncol = get_Ncolumns(lp);
	row = (double*) malloc(Ncol * sizeof(*row)); 
	get_variables(lp, row);

	//printf("%f\n", *score);
	
	//for(j = 0; j < Ncol; j++) {
	//	printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
	//}
	
	for(j = 0; j < Ncol; j++) {
		name = get_col_name(lp, j+1);
		if (name[0] == 'x') { /* this is for taxa set */
			index = -1;
			index = atoi(name+1);
			//printf(name);
			if (index < 0 || index >= ntaxa) {
				printf("Index x_%d is not in the range!\n", index);
				ret = 6;
				break;
			}
			if (row[j] > tolerance && (1.0 - row[j]) > tolerance) {
				if (verbose_mode) printf("\n%s = %10.8f", name, row[j]);
				ret = 7;
				if (!verbose_mode) break;
			}
			variables[index] = row[j];
		}
		//printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
	}
	free(row);
	
	/* we are done now */
	}
	if(lp != NULL) {
		/* clean up such that all used memory by lpsolve is freed */
		delete_lp(lp);
	}

	return ret;
}


int lp_demo() {  

  lprec *lp;
  int Ncol, *colno = NULL, j, ret = 0;
  REAL *row = NULL;

  /* We will build the model row by row
     So we start with creating a model with 0 rows and 2 columns */
  Ncol = 2; /* there are two variables in the model */
  lp = make_lp(0, Ncol);
  if(lp == NULL)
    ret = 1; /* couldn't construct a new model... */

  if(ret == 0) {
    /* let us name our variables. Not required, but can be useful for debugging */
    set_col_name(lp, 1, "x");
    set_col_name(lp, 2, "y");

    /* create space large enough for one row */
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL))
      ret = 2;
  }

  if(ret == 0) {
    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

    /* construct first row (120 x + 210 y <= 15000) */
    j = 0;

    colno[j] = 1; /* first column */
    row[j++] = 120;

    colno[j] = 2; /* second column */
    row[j++] = 210;

    /* add the row to lpsolve */
    if(!add_constraintex(lp, j, row, colno, LE, 15000))
      ret = 3;
  }

  if(ret == 0) {
    /* construct second row (110 x + 30 y <= 4000) */
    j = 0;

    colno[j] = 1; /* first column */
    row[j++] = 110;

    colno[j] = 2; /* second column */
    row[j++] = 30;

    /* add the row to lpsolve */
    if(!add_constraintex(lp, j, row, colno, LE, 4000))
      ret = 3;
  }

  if(ret == 0) {
    /* construct third row (x + y <= 75) */
    j = 0;

    colno[j] = 1; /* first column */
    row[j++] = 1;

    colno[j] = 2; /* second column */
    row[j++] = 1;

    /* add the row to lpsolve */
    if(!add_constraintex(lp, j, row, colno, LE, 75))
      ret = 3;
  }

  if(ret == 0) {
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    /* set the objective function (143 x + 60 y) */
    j = 0;

    colno[j] = 1; /* first column */
    row[j++] = 143;

    colno[j] = 2; /* second column */
    row[j++] = 60;

    /* set the objective in lpsolve */
    if(!set_obj_fnex(lp, j, row, colno))
      ret = 4;
  }

  if(ret == 0) {
    /* set the object direction to maximize */
    set_maxim(lp);

    /* just out of curioucity, now show the model in lp format on screen */
    /* this only works if this is a console application. If not, use write_lp and a filename */
    write_LP(lp, stdout);
    /* write_lp(lp, "model.lp"); */

    /* I only want to see important messages on screen while solving */
    set_verbose(lp, IMPORTANT);

    /* Now let lpsolve calculate a solution */
    ret = solve(lp);
    if(ret == OPTIMAL)
      ret = 0;
    else
      ret = 5;
  }

  if(ret == 0) {
    /* a solution is calculated, now lets get some results */

    /* objective value */
    printf("Objective value: %f\n", get_objective(lp));

    /* variable values */
    get_variables(lp, row);
    for(j = 0; j < Ncol; j++)
      printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);

    /* we are done now */
  }

  /* free allocated memory */
  if(row != NULL)
    free(row);
  if(colno != NULL)
    free(colno);

  if(lp != NULL) {
    /* clean up such that all used memory by lpsolve is freed */
    delete_lp(lp);
  }

  return(ret);
}
