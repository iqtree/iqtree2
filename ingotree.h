/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  BUI Quang Minh <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef INGOTREE_H
#define INGOTREE_H

#include "iqptree.h"

/**
 * this class implements the site-specific state frequency model
 */
class IngoTree : public IQPTree
{

public:
	/** constructor
	 * 
	 */
	IngoTree(Alignment *alignment, char *site_freq_file); 
	
	/** read in site-specific state frequency
	 *  @param site_freq_file input file name
	 */
	void readSiteFreq(char *site_freq_file);
	

	/**
	 * @param site site ID
	 * @return site state frequency vector for a site
	 */
	double *getSiteFreq(int site);
	
    ~IngoTree();
	
protected:
	double *site_freq;
	
};

#endif // INGOTREE_H
