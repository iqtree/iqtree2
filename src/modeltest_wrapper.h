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

#ifndef _MODELTEST_WRAPPER
#define _MODELTEST_WRAPPER


#ifdef __cplusplus
extern "C" {
#endif

/**
	run David Posada's modeltest
	@param score_file the input score file
	@param out_file the output file of modeltest
	@param model (OUT) name of the model selected
*/
int modeltest(char *args, const char *score_file, const char *out_file, char *LRT_model, char *IC_model);

#ifdef __cplusplus
}
#endif


#endif