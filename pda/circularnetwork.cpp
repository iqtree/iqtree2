/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
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
#include "circularnetwork.h"

CircularNetwork::CircularNetwork()
 : PDNetwork()
{
}

CircularNetwork::CircularNetwork(Params &params) : PDNetwork(params) {
}

/********************************************************
	MAIN FUNCTION
********************************************************/

void CircularNetwork::findPD(Params &params, vector<SplitSet> &taxa_set, 
	vector<int> &taxa_order) {
	
	if (!isCircular() || params.run_mode == EXHAUSTIVE || params.run_mode == GREEDY 
		|| params.run_mode == LINEAR_PROGRAMMING || isPDArea()) {
		// call inherited findPD if condition not met
		PDNetwork::findPD(params, taxa_set, taxa_order);
		return;
	}
	// call the entering function
	enterFindPD(params);
	params.detected_mode = DYNAMIC_PROGRAMMING;

	int root_id = (initialset.size() > 0) ? initialset[0] : -1;
	
	if (isBudgetConstraint()) {
		// resize the taxa_set
		taxa_set.resize(params.budget - params.min_budget + 1);

		cout << endl << "Dynamic programming on circular split network..." << endl;
		if (root_id < 0)
			findCircularPDBudget(params, taxa_set, taxa_order);
		else
			findCircularRootedPDBudget(params, taxa_set, taxa_order, root_id);
	} else	{
		// resize the taxa_set
		taxa_set.resize(params.sub_size - params.min_size + 1);
		
		cout << endl << "Dynamic programming on circular split network..." << endl;
		if (root_id < 0)
			findCircularPD(params, taxa_set, taxa_order);
		else {
			findCircularRootedPD(params, taxa_set, taxa_order, root_id);
		}
	}
	// call the leaving function
	leaveFindPD(taxa_set);

}




/**
	display the matrix into out (another version)
*/
template <class T>
void reportMyMat(ostream &out, mmatrix(T) &mat) {
	unsigned int i, j;
	for (i = 0; i < mat.size(); i++) {
		for (j = 0; j < mat[i].size(); j++) {
			if (mat[i][j] == 0) 
				out << " - &  "; 
			else if (j < mat[i].size()-1) 
				out << mat[i][j] << " & ";
			else
				out << mat[i][j];
		}
		if (i < mat.size()-1)
			out << " \\\\";
		out << endl;
	}
} 



void CircularNetwork::computePDInfo(Params &params, DoubleMatrix &table, 
		 DoubleMatrix &dist, int root) {
	int ntaxa = getNTaxa();
	int v, k, w;
	// allocate memory to solution table, set everything to ZERO
	table.resize(params.sub_size-1);
	for (k = 0; k < params.sub_size-1; k++) {
		table[k].resize(ntaxa);
		for (v = root+1; v < ntaxa; v++)
			table[k][v] = INT_MIN;
	}
	//table.setZero();

	
	for (v = root+1; v < ntaxa; v++) {
		// initialize cube[0] to distance matrix
		table[0][v] = dist[root][v];
		// now iteratively calculate cube[k]
		for (k = 1; k < params.sub_size-1 && k < v-root; k++) {
			for (w = k+root; w < v; w++) {
				double sum = table[k-1][w] + dist[v][w];
				if (table[k][v] < sum) {
					table[k][v] = sum;
				}
			}
		}
	}
	//cout << table;
}

double CircularNetwork::computePDScore(int sub_size, DoubleMatrix &table, int root) {
	int ntaxa = getNTaxa();

	int v;
	double max_pd = INT_MIN;
	for (v = root+1; v < ntaxa; v++) {
		if (max_pd < table[0][v] + table[sub_size-2][v]) {
			max_pd = table[0][v] + table[sub_size-2][v];
		}
	}
	return max_pd / 2.0;
}

void rotateTaxaOrder(vector<int> &origin_order, vector<int> &new_order, int root) {
	int ntaxa = origin_order.size();
	int i, id = ntaxa;
	for (i = 0; i < ntaxa; i++) 
		if (origin_order[i] == root) { id = i; break; }
	ASSERT(id < ntaxa);
	new_order.resize(ntaxa);
	for (i = 0; i < ntaxa; i++)
		new_order[i] = origin_order[(i+id) % ntaxa];
}


void CircularNetwork::findCircularPD(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order) {

	int ntaxa = getNTaxa();
	DoubleMatrix dist;

	int k;

	// calculate the distance matrix
	DoubleMatrix table;
	calcDistance(dist, taxa_order);

	for (int root = 0; root <= ntaxa-params.min_size; root++) {
		// dynamic programming main procedure
		// compute table information
		computePDInfo(params, table, dist, root);

		// now construction the optimal PD sets
		for (k = params.min_size; k <= params.sub_size; k++) {
			int index = k - params.min_size;
			double pd_score = computePDScore(k, table, root);

			if (taxa_set[index].getWeight() < pd_score)
				taxa_set[index].removeAll();
			else if (taxa_set[index].getWeight() > pd_score || !params.find_all) 
				// if old pd score is better or equal but not find all, continue
				continue;

			constructPD(k, params.find_all, params.pd_limit, table, dist, taxa_set[index], taxa_order, root);
		}
	}
}


void CircularNetwork::findCircularRootedPD(Params &params, vector<SplitSet> &taxa_set, 
	vector<int> &origin_order, int root) {

	DoubleMatrix dist;

	int k;
	DoubleMatrix table;

	vector<int> taxa_order;
	// rotate the position of root to 0 in the taxa_order
	rotateTaxaOrder(origin_order, taxa_order, root);

	// calculate the distance matrix
	calcDistance(dist, taxa_order);

	// dynamic programming main procedure
	computePDInfo(params, table, dist, 0);

	for (k = params.min_size; k <= params.sub_size; k++) {
		constructPD(k, params.find_all, params.pd_limit, table, dist, taxa_set[k-params.min_size], taxa_order, 0);
	}
}


void CircularNetwork::constructPD(int sub_size, bool find_all, int pd_limit, DoubleMatrix &table, 
	DoubleMatrix  &dist, SplitSet &taxa_set, vector<int> &taxa_order, int root) {
	int ntaxa = getNTaxa();
	double max_pd = INT_MIN;
	vector<int> vec_v;
	int max_v, k, v, w, sp;

	// first calculate the PD score
	for (v = root+1; v < ntaxa; v++) {
		if (max_pd < table[0][v] + table[sub_size-2][v]) {
			max_pd = table[0][v] + table[sub_size-2][v];
			max_v = v;
		}
	}

	if (find_all) {
		vec_v.push_back(max_v);
		// find all v with the same max_pd
		for (v = max_v+1; v < ntaxa; v++) {
			if (max_pd == table[0][v] + table[sub_size-2][v]) {
				vec_v.push_back(v);
			}
		}
	} else {
		// otherwise, only use max_v
		vec_v.push_back(max_v);
	}

	// now try all v
	for (sp = 0; sp < vec_v.size(); sp ++) {

		max_v = vec_v[sp];

		Split *pd_set = new Split(ntaxa, max_pd / 2.0);
	
		pd_set->addTaxon(taxa_order[root]);
		pd_set->addTaxon(taxa_order[max_v]);

		if (!find_all) {
			for (k = sub_size-2; k >= 1; k--) {
				double max = INT_MIN;
				int max_w = 0;
				for (w = root+1; w < max_v; w++) 
					if (max < table[k-1][w] + dist[max_v][w]) {
						max = table[k-1][w] + dist[max_v][w];
						max_w = w;
					}
				pd_set->addTaxon(taxa_order[max_w]);
				max_v = max_w;
			}
			taxa_set.push_back(pd_set);
		} else 
			constructPD(sub_size-2, max_v, pd_limit, pd_set, table, dist, taxa_set, taxa_order, root);
	}
}

void CircularNetwork::constructPD(int sub_size, int max_v, int pd_limit, Split *pd_set, DoubleMatrix &table,
	DoubleMatrix &dist, SplitSet &taxa_set, vector<int> &taxa_order, int root) {

	if (sub_size == 0) {
		taxa_set.push_back(pd_set);
		return;
	}

	int k;

	for (k = sub_size; k >= 1; k--) {
		int w, max_w = 0;
		double max = INT_MIN;
	
		for (w = root+1; w < max_v; w++) 
			if (max < table[k-1][w] + dist[max_v][w]) {
				max = table[k-1][w] + dist[max_v][w];
				max_w = w;
			}

		// recursive if find another PD set	
		for (w = max_w+1; w < max_v && taxa_set.size() < pd_limit; w++) 
			if (max == table[k-1][w] + dist[max_v][w]) {
				Split *new_pd = new Split(*pd_set);
				new_pd->addTaxon(taxa_order[w]);
				constructPD(k-1, w, pd_limit, new_pd, table, dist, taxa_set, taxa_order, root);
			}
	
		pd_set->addTaxon(taxa_order[max_w]);
		max_v = max_w;
		//constructPD(k-1, max_u, max_w, pd_set, table, taxa_set, taxa_order);
	}	
	taxa_set.push_back(pd_set);

}

/********************************************************
	CIRCULAR NETWORKS WITH BUDGET CONSTRAINT
********************************************************/

void CircularNetwork::calcMaxBudget(int budget, mmatrix(int) &max_b, vector<int> &taxa_order) {
	int ntaxa = getNTaxa();
	int u, v;
	max_b.resize(ntaxa-1);
	for (u = 0; u < ntaxa-1; u++) {
		max_b[u].resize(ntaxa);
		max_b[u][u] = pda->costs[taxa_order[u]];
		if (max_b[u][u] > budget) 
			max_b[u][u] = budget;
		for (v = u+1; v < ntaxa; v++) {
			max_b[u][v] = max_b[u][v-1] + pda->costs[taxa_order[v]];
			if (max_b[u][v] > budget) 
				max_b[u][v] = budget;
		}
	}
	for (u = 0; u < ntaxa-1; u++)
		for (v = u+1; v < ntaxa; v++)
			max_b[u][v] -= (pda->costs[taxa_order[u]] + pda->costs[taxa_order[v]]);
}



void CircularNetwork::constructPDBudget(int budget, bool find_all, mmatrix(double) &table,
	mmatrix(double) &dist, SplitSet &taxa_set, vector<int> &taxa_order, 
	mmatrix(int) &max_b, int root) {

	int ntaxa = getNTaxa();
	// now trace back to get the maximum pd_k
	double max_pd = INT_MIN;
	int v, s, b, sp;
	int max_v = -1, total_b;
	vector<int> vec_v;
	// reduce the budget
	budget -= pda->costs[taxa_order[root]];

	for (v = root+1; v < ntaxa; v++) {
		total_b = budget - pda->costs[taxa_order[v]];
		if (total_b > max_b[root][v]) total_b = max_b[root][v];
		if (total_b < 0) continue;
		if (max_pd < dist[root][v] + table[v][total_b]) {
			max_pd = dist[root][v] + table[v][total_b];
			max_v = v;
		}
	}

	// check if find something...
	if (max_v < 0)
		return;

	if (find_all) {
		// find all u,v with the same max_pd
		vec_v.push_back(max_v);
		for (v = max_v+1; v < ntaxa; v++) {
			total_b = budget - pda->costs[taxa_order[v]];
			if (total_b > max_b[root][v]) total_b = max_b[root][v];
			if (total_b < 0) continue;
			if (max_pd == dist[root][v] + table[v][total_b]) {
				vec_v.push_back(v);
			}
		}
	} else {
		// otherwise, only use max_v
		vec_v.push_back(max_v);
	}

	// now try all u,v
	for (sp = 0; sp < vec_v.size(); sp ++) {

		max_v = vec_v[sp];

		Split* pd_set = new Split(ntaxa, max_pd / 2.0);
		pd_set->addTaxon(taxa_order[root]);
		pd_set->addTaxon(taxa_order[max_v]);
		if (!find_all) {
			b = budget - pda->costs[taxa_order[max_v]];
			if (b > max_b[root][max_v]) b = max_b[root][max_v];
			// now trace to the minimum budget required
			while (b > 0 && table[max_v][b] == table[max_v][b-1]) b--;

			// iteratively find taxa inbetween
			while (b >= 0) {
				double max = INT_MIN;
				int max_s = -1;
				for (s = root+1; s < max_v; s++) 
					if (b >= pda->costs[taxa_order[s]]) {
						int sub_b = b - pda->costs[taxa_order[s]];
						if (sub_b > max_b[root][s]) sub_b = max_b[root][s];
						if (sub_b < 0) continue;
						if (max < dist[s][max_v] + table[s][sub_b]) {
							max = dist[s][max_v] + table[s][sub_b];
							max_s = s;
						}
					}
				if (max_s == -1) break;
				pd_set->addTaxon(taxa_order[max_s]);
				b -= pda->costs[taxa_order[max_s]];
				if (b > max_b[root][max_s])
					b = max_b[root][max_s];
				max_v = max_s;
			}
			taxa_set.push_back(pd_set);
		} else {
			b = budget - pda->costs[taxa_order[max_v]];
			if (b > max_b[root][max_v]) b = max_b[root][max_v];
			constructPDBudget(b, max_v, pd_set, table, dist, taxa_set, taxa_order, max_b, root);
		}
	}
}


void CircularNetwork::constructPDBudget(int budget, int max_v, Split *pd_set, 
	mmatrix(double) &table, mmatrix(double) &dist, SplitSet &taxa_set, 
	vector<int> &taxa_order, mmatrix(int) &max_b, int root) {

	int b = budget;

	while (b >= 0) {
		int s, max_s = -1;
		double max = INT_MIN;
	
		for (s = root+1; s < max_v; s++) 
			if (b >= pda->costs[taxa_order[s]]) {
				int sub_b = b - pda->costs[taxa_order[s]];
				if (sub_b > max_b[root][s]) sub_b = max_b[root][s];
				if (sub_b < 0) continue;
				if (max < dist[s][max_v] + table[s][sub_b]) {
					max = dist[s][max_v] + table[s][sub_b];
					max_s = s;
				}
			}

		if (max_s < 0) break;

		// recursive if find another PD set	
		for (s = max_s+1; s < max_v; s++) 
			if (b >= pda->costs[taxa_order[s]]) {
				int sub_b = b - pda->costs[taxa_order[s]];
				if (sub_b > max_b[root][s]) sub_b = max_b[root][s];
				if (sub_b < 0) continue;
				if (max == dist[s][max_v] + table[s][sub_b]) {
					Split *new_pd = new Split(*pd_set);
					new_pd->addTaxon(taxa_order[s]);
					constructPDBudget(sub_b, s, new_pd, table, dist, 
						taxa_set, taxa_order, max_b, root);
				}
			}

		pd_set->addTaxon(taxa_order[max_s]);
		b -= pda->costs[taxa_order[max_s]];
		if (b > max_b[root][max_s]) b = max_b[root][max_s];
		max_v = max_s;
	}

	taxa_set.push_back(pd_set);
}

void CircularNetwork::computePDBudgetInfo(Params &params, mmatrix(double) &table, mmatrix(int) &id, 
	mmatrix(double) &dist, vector<int> &taxa_order, mmatrix(int) &max_b, int root)
{
	int ntaxa = getNTaxa();

	int v, s, b, total_b;

	// allocate memory and initialize table
	table.resize(ntaxa);
	if (verbose_mode >= VB_DEBUG) {
		id.resize(ntaxa);
		for (v = 0; v <= root; v++) {
			for (b = 0; b < table[v].size(); b++)
				table[v][b] = 0;			
			for (b = 0; b < id[v].size(); b++)
				id[v][b] = 0;
		}
	}
	for (v = root+1; v < ntaxa; v++) {
		total_b = max_b[root][v];
		if (total_b < 0) continue;
		table[v].resize(total_b + 1, 0);
		// init table[v][b]
		for (b = 0; b <= total_b; b++)
			table[v][b] = 0;
		if (verbose_mode >= VB_DEBUG) {
			id[v].resize(total_b + 1, 0);
			for (b = 0; b <= total_b; b++) 
				id[v][b] = 0;
		}
		for (b = 0; b <= total_b; b++)
			table[v][b] = dist[root][v];
	}

	// dynamic programming
	for (v = root+2; v < ntaxa; v++) {
		total_b = max_b[root][v];
		if (total_b < 0) continue;

		for (s = root+1; s < v; s++)
			for (b = pda->costs[taxa_order[s]]; b <= total_b; b++) {
				int sub_b = b - pda->costs[taxa_order[s]];
				if (sub_b > max_b[root][s]) sub_b = max_b[root][s];
				double sum = dist[v][s] + table[s][sub_b];
				if (table[v][b] < sum) {
					table[v][b] = sum;
					if (verbose_mode >= VB_DEBUG)
						id[v][b] = s+1;
				}
			}
	}

	if (verbose_mode >= VB_DEBUG)	{
		reportMyMat(cout, table);
		reportMyMat(cout, id);
	}

}


double CircularNetwork::computePDBudgetScore(int budget, mmatrix(double) &table,
	mmatrix(double) &dist, vector<int> &taxa_order, mmatrix(int) &max_b, int root) {

	int ntaxa = getNTaxa();
	double max_pd = INT_MIN;
	int v;
	int total_b;

	budget -= pda->costs[taxa_order[root]];
	for (v = root+1; v < ntaxa; v++) {
		total_b = budget - pda->costs[taxa_order[v]];
		if (total_b > max_b[root][v]) total_b = max_b[root][v];
		if (total_b < 0) continue;
		if (max_pd < dist[root][v] + table[v][total_b]) {
			max_pd = dist[root][v] + table[v][total_b];
		}
	}
	return max_pd / 2.0;
}


void CircularNetwork::findCircularRootedPDBudget(Params &params, vector<SplitSet> &taxa_set, 
	vector<int> &origin_order, int root)
{
	//int ntaxa = getNTaxa();
	int b;

	vector<int> taxa_order;
	// rotate the position of root to 0 in the taxa_order
	rotateTaxaOrder(origin_order, taxa_order, root);

	mmatrix(double) dist;
	// calculate the distance matrix
	calcDistance(dist, taxa_order);

	mmatrix(int) max_b;
	// calculate maximum required budget from u to v
	calcMaxBudget(params.budget, max_b, taxa_order);

	mmatrix(double) table;
	mmatrix(int) id;

	// compute table and id information
	computePDBudgetInfo(params, table, id, dist, taxa_order, max_b, 0);

	// now construction the optimal PD sets
	int set_count = 0;
	for (b = params.min_budget; b <= params.budget; b++) {
		constructPDBudget(b, params.find_all, table, dist, taxa_set[b-params.min_budget], taxa_order, max_b, 0);
		if (verbose_mode >= VB_MAX) {
			cout << "budget " << b << ": " << taxa_set.size()-set_count << " set(s)" << endl;
			set_count = taxa_set.size();
		}
	}
}


void CircularNetwork::findCircularPDBudget(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order) {

	int ntaxa = getNTaxa();
	int b;

	mmatrix(double) dist;
	// calculate the distance matrix
	calcDistance(dist, taxa_order);

	if (verbose_mode >= VB_DEBUG)	{
		reportMyMat(cout, dist);
	}


	mmatrix(int) max_b;
	// calculate maximum required budget from u to v
	calcMaxBudget(params.budget, max_b, taxa_order);

	mmatrix(double) table;
	mmatrix(int) id;

	int root;

	for (root = 0; root < ntaxa-1; root++) {
		// compute table and id information
		computePDBudgetInfo(params, table, id, dist, taxa_order, max_b, root);

		// now construction the optimal PD sets
		for (b = params.min_budget; b <= params.budget; b++) {
			int index = b - params.min_budget;
			double pd_score = computePDBudgetScore(b, table, dist, taxa_order, max_b, root);
			// if the current set is already better, continue
			if (taxa_set[index].getWeight() < pd_score)
				taxa_set[index].removeAll();
			else if (taxa_set[index].getWeight() > pd_score || !params.find_all) 
				// if old pd score is better or equal but not find all, continue
				continue;
			constructPDBudget(b, params.find_all, table, dist, taxa_set[index], taxa_order, max_b, root);
		}
	}
}

