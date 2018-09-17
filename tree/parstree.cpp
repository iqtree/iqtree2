/*
 * parstree.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: diep
 */

#include <cstring>
#include "parstree.h"
#include "utils/tools.h"

ParsTree::ParsTree(): IQTree(){
}

ParsTree::ParsTree(Alignment *alignment): IQTree(alignment){
}

ParsTree::~ParsTree() {
}

void ParsTree::loadCostMatrixFile(char * file_name){
    if(cost_matrix){
        aligned_free(cost_matrix);
        cost_matrix = NULL;
    }
//    if(strcmp(file_name, "fitch") == 0)
////    if(file_name == NULL)
//    	cost_matrix = new SankoffCostMatrix(aln->num_states);
//    else
//    	cost_matrix = new SankoffCostMatrix(file_name);

    if(strcmp(file_name, "fitch") == 0 || strcmp(file_name, "e") == 0) { // uniform cost
    	cost_nstates = aln->num_states;
    	cost_matrix = aligned_alloc<unsigned int>(cost_nstates * cost_nstates);
		for(int i = 0; i < cost_nstates; i++)
			for(int j = 0; j < cost_nstates; j++){
				if(j == i) cost_matrix[i * cost_nstates + j] = 0;
				else cost_matrix[i * cost_nstates + j] = 1;
			}
    } else{ // Sankoff cost
        cout << "Loading cost matrix from " << file_name << "..." << endl;
		ifstream fin(file_name);
		if(!fin.is_open()){
			outError("Reading cost matrix file cannot perform. Please check your input file!");
		}
		fin >> cost_nstates;

		// allocate memory for cost_matrix
    	cost_matrix = aligned_alloc<unsigned int>(cost_nstates * cost_nstates);

		// read numbers from file
		for(int i = 0; i < cost_nstates; i++){
			for(int j = 0; j < cost_nstates; j++)
				fin >> cost_matrix[i * cost_nstates + j];
		}

		fin.close();

    }

    int i, j, k;
    bool changed = false;

    for (k = 0; k < cost_nstates; k++)
        for (i = 0; i < cost_nstates; i++)
            for (j = 0; j < cost_nstates; j++)
                if (cost_matrix[i*cost_nstates+j] > cost_matrix[i*cost_nstates+k] + cost_matrix[k*cost_nstates+j]) {
                    changed = true;
                    cost_matrix[i*cost_nstates+j] = cost_matrix[i*cost_nstates+k] + cost_matrix[k*cost_nstates+j];
                }

    if (changed) {
        cout << "WARING: Cost matrix does not satisfy triangular inenquality and is automatically fixed to:" << endl;
        cout << cost_nstates << endl;
        for (i = 0; i < cost_nstates; i++) {
            for (j = 0; j < cost_nstates; j++)
                cout << "  " << cost_matrix[i*cost_nstates+j];
            cout << endl;
        }
    } else {
        cout << "Cost matrix satisfies triangular inenquality" << endl;
    }


}


//inline UINT computeCostFromChild(UINT child_cost, UINT transition_cost){
//    return (child_cost + transition_cost);
//}

//
//void ParsTree::initParsData(Params* pars_params) {
//    if(!pars_params) return;
//    initCostMatrixLinear();
//    if(cost_matrix == NULL) loadCostMatrixFile(pars_params->sankoff_cost_file);
//}

void ParsTree::printPatternScore() {
//    for(int i = 0; i < aln->getNPattern(); i++)
//        cout << _pattern_pars[i] << ", ";
}

// find minimum spanning tree score for a given pattern
UINT ParsTree::findMstScore(int ptn) {

	//--- Initialize site_states to mark which state character is present in the pattern # 'ptn'
	UINT * site_states = new UINT[aln->num_states];
	// site_states[i] = 0 => state i is present, nonzero means it's absent
	for(int i = 0; i < aln->num_states; i++) site_states[i] = UINT_MAX;
	Pattern pat = aln->at(ptn);
	for(int j = 0; j < pat.size(); j++){
		if(pat[j] < aln->num_states) site_states[pat[j]] = 0;
//		else initLeafSiteParsForAmbiguousState(pat[j], site_states)
	}

	int state_count = 0;
	for(int i = 0; i < aln->num_states; i++)
		if(site_states[i] == 0) state_count++;
	if(state_count <= 1) return 0;

//	cout << "state_count = " << state_count << endl;

	//--- Prim algorithm
	UINT * labelled_value = new UINT[aln->num_states];
	bool * added = new bool[aln->num_states];
	for(int i = 0; i < aln->num_states; i++) labelled_value[i] = UINT_MAX;
	for(int i = 0; i < aln->num_states; i++) added[i] = false;

	int add_node;
//	labeled_value[0] = 0;
	int count = 0;

	do{
		if(count == 0){
			for(int c = 0; c < aln->num_states; c++){
				if((added[c] == false) && (site_states[c] == 0)){
					labelled_value[c] = 0;
//					cout << "c = " << c << endl;
					break;
				}
			}
		}
		// find among nodes unadded the one with smallest value
		int min_label = UINT_MAX;
		add_node = -1;
		for(int c = 0; c < aln->num_states; c++){
			if((added[c] == false) && (site_states[c] == 0))
				if(labelled_value[c] < min_label){
					min_label = labelled_value[c];
					add_node = c;
				}
		}

		if(add_node >= 0){
			added[add_node] = true;
			count++;
		}else break;

		// update adjacent list
		for(int c = 0; c < aln->num_states; c++)
			if((site_states[c] == 0) && (added[c] == false)){
				if(labelled_value[c] > cost_matrix[add_node * cost_nstates + c])
					labelled_value[c] = cost_matrix[add_node * cost_nstates + c];
			}
	}while(count < aln->num_states);

	UINT score = 0;
	for(int i = 0; i < aln->num_states; i++)
		if(site_states[i] == 0)
			score += labelled_value[i];

	delete [] site_states;
	delete [] labelled_value;
	delete [] added;
	return score;

}

