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

#include "modeltest/modeltest3_7.h"
#include "modeltest_wrapper.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

extern char* modelIC;
extern char *modelhLRT;

void startLogFile();
void endLogFile();

int modeltest(char *args, const char *score_file, const char *out_file, char *LRT_model, char *IC_model) {
	char *argv[10];
	int i, argc = 1;

	char * pch;
	pch = strtok (args, " ");
	while (pch != NULL)
	{
		/*printf ("%s\n",pch);*/
		argv[argc] = malloc(strlen(pch)+1);
		strcpy(argv[argc], pch);
		argc++;
		pch = strtok (NULL, " ");
	}
	argv[argc] = malloc(strlen(score_file)+4);
	sprintf(argv[argc],"-i%s", score_file);
	argc++;
/*
	for (i = 1; i < argc; i++)
		printf("%s ", argv[i]);
	printf("argc= %d", argc);
*/
    int    fd, fd_in;
    fpos_t pos, pos_in;


	/* CAUTION: Current version will destroy stdin */

	/*endLogFile();*/
	
#ifdef WIN32
	freopen(out_file, "a", stdout);
/*	freopen(score_file, "r", stdin);*/
#else
    fflush(stdout);
    fgetpos(stdout, &pos);
    fd = dup(fileno(stdout));
    freopen(out_file, "a", stdout);

/*    fflush(stdin);
    fgetpos(stdin, &pos_in);
    fd_in = dup(fileno(stdin));
    freopen(score_file, "r", stdin);*/
#endif

	if (!stdin || !stdout) {
		printf("Cannot open %s and %s\n", score_file, out_file);
		return EXIT_FAILURE;
	}

	int status = run_modeltest(argc, argv);

#ifdef WIN32
	freopen("CON", "w", stdout);
	/*freopen("CON", "r", stdin);*/
#else
    fflush(stdout);
    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos);        

/*    fflush(stdin);
    dup2(fd_in, fileno(stdin));
    close(fd_in);
    clearerr(stdin);
    fsetpos(stdin, &pos_in);*/

#endif

	/*startLogFile();*/
	strcpy(IC_model, modelIC);
	strcpy(LRT_model, modelhLRT);
	Free();

/*	for (; argc > 0; argc--)
		free(argv[argc-1]);*/

	return status;
}

double computePValueChiSquare(double chi_square, double df) {
	return (double)ChiSquare((float)chi_square, (float)df);	
}
