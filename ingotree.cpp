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


#include "ingotree.h"

IngoTree::IngoTree(Alignment *alignment, char *site_freq_file): IQPTree(alignment)
{
	assert(site_freq_file);
	//readSiteFreq(site_freq_file);
	//cout << "Ungroup site-patterns..." << endl;
	//aln->ungroupSitePattern();
}

IngoTree::~IngoTree()
{
	//if (site_freq) delete [] site_freq;
}

/*void IngoTree::readSiteFreq(char* site_freq_file)
{
	cout << "Reading site-specific state frequency file " << site_freq_file << " ..." << endl;
	site_freq = new double[getAlnNSite() * aln->num_states];
	memset(site_freq, 0, sizeof(double) * getAlnNSite() * aln->num_states);
	int i;
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(site_freq_file);
		double freq;
		string site_spec;
		in.exceptions(ios::badbit);
		for (; !in.eof(); ) {
			// remove the failbit
			in >> site_spec;
			if (in.eof()) break;
			IntVector site_id;
			extractSiteID(aln, site_spec.c_str(), site_id);
			if (site_id.size() == 0) throw "No site ID specified";
			int id = site_id.front();
			double *site_freq_entry = site_freq + (id * aln->num_states);
			double sum = 0;
			for (i = 0; i < aln->num_states; i++) {
				in >> freq;
				if (freq <= 0.0 || freq >= 1.0) throw "Invalid frequency entry";
				site_freq_entry[i] = freq;
				sum += freq;
			}
			if (fabs(sum-1.0) > 1e-4) throw "Frequencies do not sum up to 1";
			for (int i = 1; i < site_id.size(); i++) {
				memcpy(site_freq + (site_id[i] * aln->num_states), site_freq_entry, aln->num_states * sizeof(double));
			}
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (const char* str) {
		outError(str);
	} catch (string str) {
		outError(str);
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	}
}

double* IngoTree::getSiteFreq(int site)
{
	return site_freq + (site * aln->num_states);
}
*/