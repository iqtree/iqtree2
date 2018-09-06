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
    cost_matrix = NULL;
}

ParsTree::ParsTree(Alignment *alignment): IQTree(alignment){
    cost_matrix = NULL;
}

ParsTree::~ParsTree() {
    if(cost_matrix){
        aligned_free(cost_matrix);
        cost_matrix = NULL;
    }

    if (central_partial_pars)
        delete[] central_partial_pars;
    central_partial_pars = NULL;
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

void ParsTree::initCostMatrix(CostMatrixType cost_type) {
    if(cost_matrix){
        aligned_free(cost_matrix);
        cost_matrix = NULL;
    }
    ASSERT(aln);
    cost_nstates = aln->num_states;
    // allocate memory for cost_matrix
    cost_matrix = aligned_alloc<unsigned int>(cost_nstates * cost_nstates);

    switch (cost_type) {
        case CM_LINEAR:
            for(int i = 0; i < cost_nstates; i++){
                for(int j = 0; j < cost_nstates; j++)
                    cost_matrix[i * cost_nstates + j] = abs(i-j);
            }
            break;
        case CM_UNIFORM:
            for(int i = 0; i < cost_nstates; i++){
                for(int j = 0; j < cost_nstates; j++)
                    cost_matrix[i * cost_nstates + j] = ((i==j) ? 0 : 1);
            }
            break;
    }
    clearAllPartialLH();
}

/**
compute the tree parsimony score
@return parsimony score of the tree
*/
int ParsTree::computeParsimony(){
//	assert(root && root->isLeaf());
//    cout << "nstates = " << aln->num_states << endl;
//	clearAllPartialLH();
    assert(root->isLeaf());
    PhyloNeighbor *nei = ((PhyloNeighbor*) root->neighbors[0]);
    current_it = nei;
    assert(current_it);
    current_it_back = (PhyloNeighbor*) nei->node->findNeighbor(root);
    assert(current_it_back);

    int nptn = aln->size();
//    if(_pattern_pars == NULL) _pattern_pars = aligned_alloc<BootValTypePars>(nptn + VCSIZE_USHORT);

    return computeParsimonyBranch((PhyloNeighbor*) root->neighbors[0], (PhyloNode*) root);
}

//inline UINT computeCostFromChild(UINT child_cost, UINT transition_cost){
//    return (child_cost + transition_cost);
//}

/**
compute partial parsimony score of the subtree rooted at dad
@param dad_branch the branch leading to the subtree
@param dad its dad, used to direct the traversal
*/
void ParsTree::computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad){
	// don't recompute the parsimony
    if (dad_branch->partial_lh_computed & 2)
        return;

    Node *node = dad_branch->node;
    //assert(node->degree() <= 3);
    int ptn;
    if(aln->num_states != cost_nstates){
        cout << "Your cost matrix is not compatible with the alignment"
            << " in terms of number of states. Please check!" << endl;
        exit(1);
    }
    int nstates = aln->num_states;
    assert(dad_branch->partial_pars);

    int pars_block_size = getParsBlockSize();

    if (node->isLeaf() && dad) {
//        cout << "############# leaf!" << endl;
        // external node
        for(int i = 0; i < pars_block_size - 1; i++)
            dad_branch->partial_pars[i] = 1000;
        dad_branch->partial_pars[pars_block_size - 1] = 0; // reserved for corresponding subtree pars
        for (ptn = 0; ptn < aln->size(); ptn++){
            // ignore const ptn because it does not affect pars score
            if (!aln->at(ptn).isConst()) {
                int ptn_start_index = ptn * nstates;
                char state;
                if (node->name == ROOT_NAME) {
                    state = aln->STATE_UNKNOWN;
                } else {
                    assert(node->id < aln->getNSeq());
                    state = (aln->at(ptn))[node->id];
                }

                if (state < nstates) {
                    dad_branch->partial_pars[ptn_start_index + state] = 0;
                } else {
                    // unknown, ambiguous character
//                    cout << "####### ambigous state = " << int(state) << endl;
                    initLeafSiteParsForAmbiguousState(state, dad_branch->partial_pars + ptn_start_index);
                }
            }
        }
    } else {
//        cout << "############# internal!" << endl;
        // internal node
        UINT i, j, ptn, min_child_ptn_pars;

        UINT * partial_pars = dad_branch->partial_pars;
        for(int i = 0; i < pars_block_size; i++)
            partial_pars[i] = 0;
        UINT *left = NULL, *right = NULL;

        FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
            computePartialParsimony((PhyloNeighbor*) (*it), (PhyloNode*) node);
            if (!left)
                left = ((PhyloNeighbor*)*it)->partial_pars;
            else
                right = ((PhyloNeighbor*)*it)->partial_pars;
        }


        if (node->degree() > 3) {
            FOR_NEIGHBOR_IT(node, dad, it) if ((*it)->node->name != ROOT_NAME) {
                UINT *partial_pars_child = ((PhyloNeighbor*) (*it))->partial_pars;
                for (ptn = 0; ptn < aln->size(); ptn++){
                    // ignore const ptn because it does not affect pars score
                    if (aln->at(ptn).isConst()) continue;
                    int ptn_start_index = ptn*nstates;

                    UINT *partial_pars_child_ptr = &partial_pars_child[ptn_start_index];
                    UINT *partial_pars_ptr = &partial_pars[ptn_start_index];
                    UINT *cost_matrix_ptr = cost_matrix;

                    for(i = 0; i < nstates; i++){
                        // min(j->i) from child_branch
                        min_child_ptn_pars = partial_pars_child_ptr[0] + cost_matrix_ptr[0];
                        for(j = 1; j < nstates; j++) {
                            UINT value = partial_pars_child_ptr[j] + cost_matrix_ptr[j];
                            min_child_ptn_pars = min(value, min_child_ptn_pars);
                        }
                        partial_pars_ptr[i] += min_child_ptn_pars;
                        cost_matrix_ptr += nstates;
                    }
                }
            }
        } else {
            // bifurcating node
            assert(node->degree() == 3);

            switch (nstates) {
            case 4:
                for (ptn = 0; ptn < aln->size(); ptn++){
                    // ignore const ptn because it does not affect pars score
                    if (aln->at(ptn).isConst()) continue;
                    int ptn_start_index = ptn*4;

                    UINT *left_ptr = &left[ptn_start_index];
                    UINT *right_ptr = &right[ptn_start_index];
                    UINT *partial_pars_ptr = &partial_pars[ptn_start_index];
                    UINT *cost_matrix_ptr = cost_matrix;
                    UINT left_contrib, right_contrib;

                    for(i = 0; i < 4; i++){
                        // min(j->i) from child_branch
                        left_contrib = left_ptr[0] + cost_matrix_ptr[0];
                        right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                        for(j = 1; j < 4; j++) {
                            UINT value = left_ptr[j] + cost_matrix_ptr[j];
                            left_contrib = min(value, left_contrib);
                            value = right_ptr[j] + cost_matrix_ptr[j];
                            right_contrib = min(value, right_contrib);
                        }
                        partial_pars_ptr[i] = left_contrib+right_contrib;
                        cost_matrix_ptr += 4;
                    }
                }
                break;
            case 20:
                for (ptn = 0; ptn < aln->size(); ptn++){
                    // ignore const ptn because it does not affect pars score
                    if (aln->at(ptn).isConst()) continue;
                    int ptn_start_index = ptn*20;

                    UINT *left_ptr = &left[ptn_start_index];
                    UINT *right_ptr = &right[ptn_start_index];
                    UINT *partial_pars_ptr = &partial_pars[ptn_start_index];
                    UINT *cost_matrix_ptr = cost_matrix;
                    UINT left_contrib, right_contrib;

                    for(i = 0; i < 20; i++){
                        // min(j->i) from child_branch
                        left_contrib = left_ptr[0] + cost_matrix_ptr[0];
                        right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                        for(j = 1; j < 20; j++) {
                            UINT value = left_ptr[j] + cost_matrix_ptr[j];
                            left_contrib = min(value, left_contrib);
                            value = right_ptr[j] + cost_matrix_ptr[j];
                            right_contrib = min(value, right_contrib);
                        }
                        partial_pars_ptr[i] = left_contrib+right_contrib;
                        cost_matrix_ptr += 20;
                    }
                }
                break;
            default:
                for (ptn = 0; ptn < aln->size(); ptn++){
                    // ignore const ptn because it does not affect pars score
                    if (aln->at(ptn).isConst()) continue;
                    int ptn_start_index = ptn*nstates;

                    UINT *left_ptr = &left[ptn_start_index];
                    UINT *right_ptr = &right[ptn_start_index];
                    UINT *partial_pars_ptr = &partial_pars[ptn_start_index];
                    UINT *cost_matrix_ptr = cost_matrix;
                    UINT left_contrib, right_contrib;

                    for(i = 0; i < nstates; i++){
                        // min(j->i) from child_branch
                        left_contrib = left_ptr[0] + cost_matrix_ptr[0];
                        right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                        for(j = 1; j < nstates; j++) {
                            UINT value = left_ptr[j] + cost_matrix_ptr[j];
                            left_contrib = min(value, left_contrib);
                            value = right_ptr[j] + cost_matrix_ptr[j];
                            right_contrib = min(value, right_contrib);
                        }
                        partial_pars_ptr[i] = left_contrib+right_contrib;
                        cost_matrix_ptr += nstates;
                    }
                }
                break;
            }

        }
        /*

        // calc subtree pars
        for (ptn = 0; ptn < aln->size(); ptn++){
            // ignore const ptn because it does not affect pars score
            if (aln->at(ptn).is_const) continue;
            int ptn_start_index = ptn * nstates;
            UINT min_ptn_pars = partial_pars[ptn_start_index];
            for(i = 1; i < nstates; i++){
                if(partial_pars[ptn_start_index + i] < min_ptn_pars)
                    min_ptn_pars = partial_pars[ptn_start_index + i];
            }
            partial_pars[pars_block_size - 1] += min_ptn_pars * aln->at(ptn).frequency;
        }
        */
    }

    dad_branch->partial_lh_computed |= 2;
}

void ParsTree::initLeafSiteParsForAmbiguousState(char state, UINT * site_partial_pars){
    int i, nstates = aln->num_states;
    if(state < nstates) return; // no need for manipulate normal state

    if (state == aln->STATE_UNKNOWN){
        for(i = 0; i < nstates; i++) site_partial_pars[i] = 0;
        return;
    }

    if (state == STATE_INVALID){
    	cout << "nstates = " << nstates << "; state = " << (int) state << endl;
        outError("Alignment contains invalid state. Please check your data!");
    }

//    for(i = 0; i < nstates; i++) site_partial_pars[i] = UINT_MAX;

    switch (nstates) {
        case 2:
        	cout << "nstates = " << nstates << "; state = " << (int) state << endl;
        	outError("Alignment contains invalid state. Please check your data!");
        	break;
        case 4: // DNA
			switch (state) {
				case 1+4+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					return; // A or G, Purine
				case 2+8+3:
					site_partial_pars[aln->convertState('C')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // C or T, Pyrimidine
				case 1+8+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // A or T, Weak
				case 2+4+3:
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					return; // G or C, Strong
				case 1+2+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					return; // A or C, Amino
				case 4+8+3:
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // G or T, Keto
				case 2+4+8+3:
					site_partial_pars[aln->convertState('C')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return;// C or G or T
				case 1+2+8+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // A or C or T
				case 1+4+8+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // A or G or T
				case 1+2+4+3:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					return; // A or G or C
				case 18:
					site_partial_pars[aln->convertState('A')] = 0;
					site_partial_pars[aln->convertState('C')] = 0;
					site_partial_pars[aln->convertState('G')] = 0;
					site_partial_pars[aln->convertState('T')] = 0;
					return; // UNKNOWN for DNA
				default:
					cout << "nstates = " << nstates << "; state = " << (int) state << endl;
					outError("Alignment contains invalid state. Please check your data!");
					return;
			}
			break;
        case 20: // Protein
        	if (state == 4+8+19){
        		site_partial_pars[aln->convertState('D')] = 0;
        		site_partial_pars[aln->convertState('N')] = 0;
        		return; // Aspartic acid (D) or Asparagine (N)
        	}
        	else if (state == 32+64+19){
        		site_partial_pars[aln->convertState('Q')] = 0;
        		site_partial_pars[aln->convertState('E')] = 0;
        		return; // Glutamine (Q) or Glutamic acid (E)
        	}
        	else if (state == 22){
				for(i = 0; i < nstates; i++) site_partial_pars[i] = 0;
        		return; // UNKNOWN for Protein
        	}
        	else{
        		cout << "nstates = " << nstates << "; state = " << (int) state << endl;
        		outError("Alignment contains invalid state. Please check your data!");
        		return;
        	}
        default:
        	// unknown
        	cout << "nstates = " << nstates << "; state = " << (int) state << endl;
        	outError("Alignment contains invalid state. Please check your data!");
        	return;
    }

}

/**
compute tree parsimony score based on a particular branch
@param dad_branch the branch leading to the subtree
@param dad its dad, used to direct the traversal
@param branch_subst (OUT) if not NULL, the number of substitutions on this branch
@return parsimony score of the tree
*/
int ParsTree::computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);

    if (!central_partial_pars)
        initializeAllPartialPars();

    // DTH: I don't really understand what this is for. ###########
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
//        cout << "swapped\n";
    }

    int nptn = aln->size();
//    if(!_pattern_pars) _pattern_pars = aligned_alloc<BootValTypePars>(nptn+VCSIZE_USHORT);
//    memset(_pattern_pars, 0, sizeof(BootValTypePars) * (nptn+VCSIZE_USHORT));

    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(node_branch, node);

    // now combine likelihood at the branch
    tree_pars = 0;
    int nstates = aln->num_states;
    UINT i, j, ptn;

    UINT *ptn_partial_pars = new UINT[nstates];
    if(!ptn_partial_pars){
        outError("Could not allocate for ptn_partial_pars\n");
        exit(1);
    }

    switch (nstates) {
    case 4:
        for (ptn = 0; ptn < aln->size(); ptn++){
            //_pattern_pars[ptn] = 0;
            if (aln->at(ptn).isConst()) continue;
            
            int ptn_start_index = ptn * 4;
            UINT *node_branch_ptr = &node_branch->partial_pars[ptn_start_index];
            UINT *dad_branch_ptr = &dad_branch->partial_pars[ptn_start_index];
            UINT *cost_matrix_ptr = cost_matrix;
            UINT min_ptn_pars = UINT_MAX;
            for(i = 0; i < 4; i++){
                // min(j->i) from node_branch
                UINT min_score = node_branch_ptr[0] + cost_matrix_ptr[0];
                for(j = 1; j < 4; j++) {
                    UINT value = node_branch_ptr[j] + cost_matrix_ptr[j];
                    min_score = min(value, min_score);

                }
                ptn_partial_pars[i] = min_score + dad_branch_ptr[i];
                min_ptn_pars = min(min_ptn_pars, ptn_partial_pars[i]);
                cost_matrix_ptr += 4;
            }
            //_pattern_pars[ptn] = min_ptn_pars;
            tree_pars += min_ptn_pars * aln->at(ptn).frequency;
        }
        break;

    default:
        for (ptn = 0; ptn < aln->size(); ptn++){
            //_pattern_pars[ptn] = 0;
            if (aln->at(ptn).isConst()) continue;
            
            int ptn_start_index = ptn * nstates;
            UINT *node_branch_ptr = &node_branch->partial_pars[ptn_start_index];
            UINT *dad_branch_ptr = &dad_branch->partial_pars[ptn_start_index];
            UINT *cost_matrix_ptr = cost_matrix;
            UINT min_ptn_pars = UINT_MAX;
            for(i = 0; i < nstates; i++){
                // min(j->i) from node_branch
                UINT min_score = node_branch_ptr[0] + cost_matrix_ptr[0];
                for(j = 1; j < nstates; j++) {
                    UINT value = node_branch_ptr[j] + cost_matrix_ptr[j];
                    min_score = min(value, min_score);

                }
                ptn_partial_pars[i] = min_score + dad_branch_ptr[i];
                min_ptn_pars = min(min_ptn_pars, ptn_partial_pars[i]);
                cost_matrix_ptr += nstates;
            }
            //_pattern_pars[ptn] = min_ptn_pars;
            tree_pars += min_ptn_pars * aln->at(ptn).frequency;
        }
        break;
    }
    if (branch_subst)
        *branch_subst = tree_pars;
    if(ptn_partial_pars){
        delete [] ptn_partial_pars;
        ptn_partial_pars = NULL;
    }
    return tree_pars;
}

void ParsTree::initializeAllPartialPars() {
	PhyloTree::initializeAllPartialPars();
//    if(params->maximum_parsimony && (!_pattern_pars))
//    	_pattern_pars = new UINT[aln->size()];
//    int index = 0;
//    initializeAllPartialPars(index);
}

size_t ParsTree::getParsBlockSize(){
    // the extra one is reserved for subtree pars score
    // TODO minus const ptn
    return aln->size() * aln->num_states + 1;
}

UINT* ParsTree::newBitsBlock(){
	return new UINT[getParsBlockSize()];
}


void ParsTree::initializeAllPartialPars(int &index, PhyloNode *node, PhyloNode *dad) {
    size_t pars_block_size = getParsBlockSize();
    if (!node) {
        node = (PhyloNode*) root;
        // allocate the big central partial pars memory
        if (!central_partial_pars) {
            int memsize = (aln->getNSeq() - 1) * 4 * pars_block_size; // Important for handling informal NEXUS format (e.g. output of TNT)
            if (verbose_mode >= VB_MED)
                cout << "Allocating " << memsize * sizeof(UINT) << " bytes for partial parsimony vectors" << endl;
            central_partial_pars = new UINT[memsize];
            if (!central_partial_pars)
                outError("Not enough memory for partial parsimony vectors");
        }
        index = 0;
    }
    if (dad) {
        // make memory alignment 16
        // assign a region in central_partial_lh to both Neihgbors (dad->node, and node->dad)
        PhyloNeighbor *nei = (PhyloNeighbor*) node->findNeighbor(dad);
        nei->partial_pars = central_partial_pars + (index * pars_block_size);
        nei = (PhyloNeighbor*) dad->findNeighbor(node);
        nei->partial_pars = central_partial_pars + ((index + 1) * pars_block_size);
        index += 2;
        //assert(index < nodeNum * 2 - 1);
    }
    FOR_NEIGHBOR_IT(node, dad, it) initializeAllPartialPars(index, (PhyloNode*) (*it)->node, node);
}

UINT* ParsTree::newPartialPars(){
	UINT *ret = new UINT[getParsBlockSize()];
    return ret;
}
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

