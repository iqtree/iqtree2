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

#include "weisslambda_sub.h"

/* following is outdated */

/*
int performDeltaTest ()
{

	FILE *fps;

	printf("Performing the test...\n");

	ausgabe_null = (char*)&ausgabe_1;
	ausgabe_data = (char*)&ausgabe_0;


	printf ( "alpha  %f ",alpha );

	printf ( "%s\n",ausgabe_data );


	fps = fopen ( "p_values","a" );

	fprintf ( fps,"file name\t%s\n\n",datei_name );

	ReadDataSets();

	fclose ( fps );


	ComputeWeissLambdafromData();

	ComputeWeissLambdafromSimulation();


	fprintf ( stderr,"\n" );
}

*/


