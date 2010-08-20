//
// C++ Implementation: phylotree
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "phylotree.h"
#include "bionj.h"
#include "rateheterogeneity.h"
#include "rateinvar.h"
#include "rategamma.h"
#include "rategammainvar.h"
#include "gtrmodel.h"
#include "modeldna.h"
#include "modelprotein.h"

const double MIN_BRANCH_LEN = 0.000001;
const double MAX_BRANCH_LEN = 9.0;
const double TOL_BRANCH_LEN = 0.00001;
const double TOL_LIKELIHOOD = 0.0001;
const double SCALING_THRESHOLD = 1e-150;
const double LOG_SCALING_THRESHOLD = log(SCALING_THRESHOLD);

/****************************************************************************
	SPRMoves class
****************************************************************************/

void SPRMoves::add(PhyloNode *prune_node, PhyloNode *prune_dad,
		PhyloNode *regraft_node, PhyloNode *regraft_dad, double score) {
	if (size() >= MAX_SPR_MOVES && score <= rbegin()->score)
		return;
	if (size() >= MAX_SPR_MOVES) {
		iterator it = end();
		it--;
		erase(it);
	}
	SPRMove spr;
	spr.prune_node = prune_node;
	spr.prune_dad = prune_dad;
	spr.regraft_node = regraft_node;
	spr.regraft_dad = regraft_dad;
	spr.score = score;
	insert(spr);
}

/****************************************************************************
	PhyloTree class
****************************************************************************/

PhyloTree::PhyloTree()
 : MTree()
{
	aln = NULL;
	model = NULL;
	site_rate = NULL;
	optimize_by_newton = false;
	central_partial_lh = NULL;
	central_partial_pars = NULL;
	model_factory = NULL;
}

PhyloTree::~PhyloTree() {
	if (central_partial_lh)
		delete [] central_partial_lh;
	central_partial_lh = NULL;
	if (central_partial_pars)
		delete [] central_partial_pars;
	central_partial_pars = NULL;
	if (model_factory) delete model_factory;
	if (model) delete model;
	if (site_rate) delete site_rate;
	if (root != NULL)
		freeNode();
	root = NULL;
}

void PhyloTree::assignLeafNames(Node *node, Node *dad) {
	if (!node) node = root;
	if (node->isLeaf()) {
		node->id = atoi(node->name.c_str());
		assert(node->id >= 0 && node->id < leafNum);
		node->name = aln->getSeqName(node->id).c_str();
	}
	FOR_NEIGHBOR_IT(node, dad, it)
		assignLeafNames((*it)->node, node);
}


void PhyloTree::copyTree(MTree *tree) {
	MTree::copyTree(tree);
	if (!aln) return;
	// reset the ID with alignment
	setAlignment(aln);
}

void PhyloTree::copyPhyloTree(PhyloTree *tree) {
	MTree::copyTree(tree);
	if (!tree->aln) return;
	setAlignment(tree->aln);
}

void PhyloTree::setAlignment(Alignment *alignment)
{
	aln = alignment;
	int nseq = aln->getNSeq();
	for (int seq = 0; seq < nseq; seq++) {
		string seq_name = aln->getSeqName(seq);
		Node *node = findNodeName(seq_name);
		if (!node) {
			string str = "Alignment has a sequence name ";
			str += seq_name;
			str += " which is not in the tree";
			outError(str);
		}
		assert(node->isLeaf());
		node->id = seq;
	}
}

void PhyloTree::setModel(SubstModel *amodel) {
	model = amodel;
}

void PhyloTree::setModelFactory(ModelFactory *model_fac) {
	model_factory = model_fac;
}

void PhyloTree::setRate(RateHeterogeneity *rate) {
	site_rate = rate;
}

RateHeterogeneity *PhyloTree::getRate() {
	return site_rate;
}


Node* PhyloTree::newNode(int node_id, const char* node_name) {
	return (Node*) (new PhyloNode(node_id, node_name));
}

Node* PhyloTree::newNode(int node_id, int node_name) {
	return (Node*) (new PhyloNode(node_id, node_name));
}

void PhyloTree::clearAllPartialLh() {
	((PhyloNode*)root->neighbors[0]->node)->clearAllPartialLh((PhyloNode*)root);
}


void PhyloTree::createModel(Params &params) {
	assert(aln);
	optimize_by_newton = params.optimize_by_newton;
	string model_str = params.model_name;

	string::size_type pos = model_str.find('+');
;

	/* create site-rate heterogeneity */
	if (pos != string::npos) {
		string rate_str = model_str.substr(pos);
		if (rate_str == "+I")
			site_rate = new RateInvar(this);
		else if (rate_str.substr(0,2) == "+G") {
			if (rate_str.length() > 2) {
				params.num_rate_cats = convert_int(rate_str.substr(2).c_str());
				if (params.num_rate_cats < 1) outError("Wrong number of rate categories");
			}
			site_rate = new RateGamma(params.num_rate_cats, this);
		} else if (rate_str.substr(0,4) == "+I+G" || rate_str == "+G+I") {
			if (rate_str.length() > 4) {
				params.num_rate_cats = convert_int(rate_str.substr(4).c_str());
				if (params.num_rate_cats < 1) outError("Wrong number of rate categories");
			}
			site_rate = new RateGammaInvar(params.num_rate_cats, this);
		} else
			outError("Invalid rate heterogeneity type");
		model_str = model_str.substr(0, pos);
	} else site_rate = new RateHeterogeneity();

	/* create substitution model */

	if (model_str == "JC" /*&& (params.freq_type == FREQ_UNKNOWN || params.freq_type == FREQ_EQUAL)*/) {
		 model = new SubstModel(aln->num_states);
	} else if (model_str == "GTR") {
		model = new GTRModel(this);
		((GTRModel*)model)->init(params.freq_type);
	} else if (aln->num_states == 4) {
		model = new ModelDNA(model_str.c_str(), params.freq_type, this);
	} else if (aln->num_states == 20) {
		model = new ModelProtein(model_str.c_str(), params.freq_type, this);
	} else {
		outError("Unsupported model type");
	}

	model_factory = new ModelFactory(model, params.store_trans_matrix);

}


string PhyloTree::getModelName() {
	return model->name + site_rate->name;
}

/****************************************************************************
	Parsimony function
****************************************************************************/



int PhyloTree::getBitsBlockSize() {
	// reserve the last entry for parsimony score
	return (aln->num_states*aln->size()+ UINT_BITS - 1) / UINT_BITS + 1;
}

UINT *PhyloTree::newBitsBlock() {
	return new UINT[getBitsBlockSize()]; 
}

UINT PhyloTree::getBitsBlock(UINT *bit_vec, int index) {
	int nstates = aln->num_states;
	int myindex = (index * nstates);
	int bit_pos_begin = myindex >> BITS_DIV;
	int bit_off_begin = myindex & BITS_MODULO;
	int bit_pos_end;
	int bit_off_end = bit_off_begin + nstates;
	if (bit_off_end > UINT_BITS) {
		bit_off_end -= UINT_BITS;
		bit_pos_end = bit_pos_begin+1;
	} else
		bit_pos_end = bit_pos_begin;


	if (bit_pos_begin == bit_pos_end)
		return (bit_vec[bit_pos_begin] >> bit_off_begin) & ((1 << nstates) - 1);
	else {
		int part1 = (bit_vec[bit_pos_begin] >> bit_off_begin);
		int part2 = bit_vec[bit_pos_end] & ((1 << bit_off_end) - 1);
		return part1 | (part2 << (UINT_BITS - bit_off_begin));
	}
}

void PhyloTree::setBitsBlock(UINT *bit_vec, int index, UINT value) {
	int nstates = aln->num_states;
	int myindex = (index * nstates);
	int bit_pos_begin = myindex >> BITS_DIV;
	int bit_off_begin = myindex & BITS_MODULO;
	int bit_pos_end;
	int bit_off_end = bit_off_begin + nstates;
	if (bit_off_end > UINT_BITS) {
		bit_off_end -= UINT_BITS;
		bit_pos_end = bit_pos_begin+1;
	} else
		bit_pos_end = bit_pos_begin;

	int allstates = (1 << nstates) - 1;
	assert(value <= allstates);

	if (bit_pos_begin == bit_pos_end) {
		// first clear the bit between bit_off_begin and bit_off_end
		bit_vec[bit_pos_begin] &= ~(allstates << bit_off_begin);
		// now set the bit
		bit_vec[bit_pos_begin] |= value << bit_off_begin;
	} else {
		int len1 = UINT_BITS - bit_off_begin;
		int allbit1 = (1 << len1) - 1;
		// clear bit from bit_off_begin to UINT_BITS
		bit_vec[bit_pos_begin] &= ~(allbit1 << bit_off_begin);
		// set bit  from bit_off_begin to UINT_BITS
		bit_vec[bit_pos_begin] |= ((value & allbit1) << bit_off_begin);

		// clear bit from 0 to bit_off_end 
		bit_vec[bit_pos_end] &= ~((1 << bit_off_end)-1);
		// now set the bit the value
		bit_vec[bit_pos_end] |= (value >> len1);
	
	}
}



void PhyloTree::computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	// don't recompute the likelihood
	if (dad_branch->partial_lh_computed & 2)
		return;
	Node *node = dad_branch->node;
	assert(node->degree() <= 3);
	int ptn;
	int nstates = aln->num_states;
	int pars_size = getBitsBlockSize();
	assert(dad_branch->partial_pars);

	if (node->isLeaf() && dad) {
		// external node 
		memset(dad_branch->partial_pars, 0, pars_size * sizeof(int));
		for (ptn = 0; ptn < aln->size(); ptn++) {
			char state;
			if (node->name == ROOT_NAME) {
				state = STATE_UNKNOWN;
			} else {
				assert(node->id < aln->getNSeq());
				state = (aln->at(ptn))[node->id];
			}
			if (state == STATE_UNKNOWN) {
				// fill all entries with bit 1
				setBitsBlock(dad_branch->partial_pars, ptn, (1 << nstates)-1);
			} else if (state < nstates) {
				setBitsBlock(dad_branch->partial_pars, ptn, 1 << state);
			} else {
				// ambiguous character, for DNA, RNA
				state = state - (nstates-1);
				setBitsBlock(dad_branch->partial_pars, ptn, state);
			}
		}
	} else {
		// internal node 
		//memset(dad_branch->partial_pars, 127, pars_size * sizeof(int));
		UINT *partial_pars_dad = dad_branch->partial_pars;
		UINT *partial_pars_child1 = NULL, *partial_pars_child2 = NULL;
		// take the intersection of two child states (with &= bit operation)
		FOR_NEIGHBOR_IT(node, dad, it)
		if ((*it)->node->name != ROOT_NAME) {
			computePartialParsimony((PhyloNeighbor*)(*it), (PhyloNode*)node);
			if (!partial_pars_child1)
				partial_pars_child1 = ((PhyloNeighbor*)(*it))->partial_pars;
			else
				partial_pars_child2 = ((PhyloNeighbor*)(*it))->partial_pars;
		}
		assert(partial_pars_child1 && partial_pars_child2);
		for (int i = 0; i < pars_size-1; i++)
			partial_pars_dad[i] = partial_pars_child1[i] & partial_pars_child2[i];
		int partial_pars = partial_pars_child1[pars_size-1] + partial_pars_child2[pars_size-1];
		// now check if some intersection is empty, change to union (Fitch algorithm) and increase the parsimony score
		for (ptn = 0; ptn < aln->size(); ptn++) 
			if (getBitsBlock(partial_pars_dad, ptn) == 0) {
				setBitsBlock(partial_pars_dad, ptn, getBitsBlock(partial_pars_child1, ptn) | getBitsBlock(partial_pars_child2, ptn));
				partial_pars += aln->at(ptn).frequency;
			}
		partial_pars_dad[pars_size-1] = partial_pars;
	}
	dad_branch->partial_lh_computed |= 2;
}


int PhyloTree::computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	PhyloNode *node = (PhyloNode*)dad_branch->node;
	PhyloNeighbor *node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
	assert(node_branch);
	if (!central_partial_pars)
		initializeAllPartialLh();
	// swap node and dad if dad is a leaf
	if (node->isLeaf()) {
		PhyloNode *tmp_node = dad;
		dad = node;
		node = tmp_node;
		PhyloNeighbor *tmp_nei = dad_branch;
		dad_branch = node_branch;
		node_branch = tmp_nei;
		//cout << "swapped\n";
	}
	if ((dad_branch->partial_lh_computed & 2) == 0)
		computePartialParsimony(dad_branch, dad);
	if ((node_branch->partial_lh_computed & 2) == 0)
		computePartialParsimony(node_branch, node);
	// now combine likelihood at the branch

	int pars_size = getBitsBlockSize();
	int nstates = aln->num_states;
	int i, ptn;
	int tree_pars = node_branch->partial_pars[pars_size-1] + dad_branch->partial_pars[pars_size-1];
	UINT *partial_pars = newBitsBlock();
	for (i = 0; i < pars_size-1; i++)
		partial_pars[i] = (node_branch->partial_pars[i] & dad_branch->partial_pars[i]);

	for (ptn = 0; ptn < aln->size(); ptn++)  {
		/*
		for (i = 0; i < aln->getNSeq(); i++)
			cout << aln->convertStateBack(aln->at(ptn)[i]);
		cout << endl;*/

		if (getBitsBlock(partial_pars, ptn) == 0)
			tree_pars += aln->at(ptn).frequency;
	}
	delete [] partial_pars;
	return tree_pars;
}

int PhyloTree::computeParsimony() {
	return computeParsimonyBranch((PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root);
}

int PhyloTree::computeParsimonyScore(int ptn, int &states, PhyloNode *node, PhyloNode *dad) {
	int score = 0;
	states = 0;
	if (!node) node = (PhyloNode*)root;
	if (node->degree() > 3)
		outError("Does not work with multifurcating tree");
	if (verbose_mode == VB_DEBUG)
		cout << ptn << " " << node->id  << "  " << node->name << endl;

	if (node->isLeaf()) {
		char state;
		if (node->name == ROOT_NAME) {
			state = STATE_UNKNOWN;
		} else {
			assert(node->id < aln->getNSeq());
			state = (*aln)[ptn][node->id];
		}
		if (state == STATE_UNKNOWN) {
			states = (1 << aln->num_states) - 1;
		} else if (state < aln->num_states)
			states = (1 << state);
		else {
			// ambiguous character, for DNA, RNA
			states = state - 3;
		}
	}
	if (!node->isLeaf() || node == root) {
		int union_states = 0;
		int intersect_states = (1 << aln->num_states) - 1;
		if (states != 0) {
			union_states = states;
			intersect_states = states;
		}
		FOR_NEIGHBOR_IT(node, dad, it) {
			int states_child;
			int score_child = computeParsimonyScore(ptn, states_child, (PhyloNode*)((*it)->node), node);
			union_states |= states_child;
			intersect_states &= states_child;
			score += score_child;
		}
		if (intersect_states)
			states = intersect_states;
		else {
			states = union_states;
			score++;
		}
	}
	return score;
}


int PhyloTree::computeParsimonyScore() {
	assert(root && root->isLeaf());

	int score = 0;
	for (int ptn = 0; ptn < aln->size(); ptn++)  
	if (!aln->at(ptn).is_const) {
		int states;
		score += computeParsimonyScore(ptn, states) * (*aln)[ptn].frequency;
	}
	return score;
}

/****************************************************************************
	Nearest Neighbor Interchange with parsimony
****************************************************************************/

double PhyloTree::swapNNI(double cur_score, PhyloNode *node1, PhyloNode *node2) {
	assert(node1->degree() == 3 && node2->degree() == 3);
	FOR_NEIGHBOR_DECLARE(node1, node2, it1)
		break;
	Node *node1_nei = (*it1)->node;
	FOR_NEIGHBOR_IT(node2, node1, it2) {
		// do the NNI swap
		Node *node2_nei = (*it2)->node;
		node1->updateNeighbor(node1_nei, node2_nei);
		node1_nei->updateNeighbor(node1, node2);
		node2->updateNeighbor(node2_nei, node1_nei);
		node2_nei->updateNeighbor(node2, node1);

		// compute the score of the swapped topology
		double score = computeParsimonyScore();
		// if better: return
		if (score < cur_score) return score;
		// else, swap back
		node1->updateNeighbor(node2_nei, node1_nei);
		node1_nei->updateNeighbor(node2, node1);
		node2->updateNeighbor(node1_nei, node2_nei);
		node2_nei->updateNeighbor(node1, node2);
	}
	return cur_score;
}

double PhyloTree::searchNNI(double cur_score, PhyloNode *node, PhyloNode *dad)  {
	if (!node) node = (PhyloNode*)root;
	if (!node->isLeaf() && dad && !dad->isLeaf()) {
		double score = swapNNI(cur_score, node, dad);
		if (score < cur_score) return score;
	}
	FOR_NEIGHBOR_IT(node, dad, it) {
		double score = searchNNI(cur_score, (PhyloNode*)(*it)->node, node);
		if (score < cur_score) return score;
	}
	return cur_score;
}

void PhyloTree::searchNNI() {
	cout << "Search with Nearest Neighbor Interchange..." << endl;
	double cur_score = computeParsimonyScore();
	do {
		double score = searchNNI(cur_score);
		if (score >= cur_score) break;
		cout << "Better score found: " << score << endl;
		cur_score = score;
	} while (true);
}


/****************************************************************************
	Stepwise addition (greedy) by maximum parsimony
****************************************************************************/

int PhyloTree::addTaxonMP(Node *added_node, Node* &target_node, Node* &target_dad, Node *node, Node *dad) {
	Neighbor *dad_nei = dad->findNeighbor(node);

	// now insert the new node in the middle of the branch node-dad
	double len = dad_nei->length;
	node->updateNeighbor(dad, added_node, len/2.0);
	dad->updateNeighbor(node, added_node, len/2.0);
	added_node->updateNeighbor((Node*)1, node, len/2.0);
	added_node->updateNeighbor((Node*)2, dad, len/2.0);
	// compute the likelihood
	//clearAllPartialLh();
	int best_score = computeParsimonyScore();
	target_node = node;
	target_dad = dad;
	// remove the added node
	node->updateNeighbor(added_node, dad, len);
	dad->updateNeighbor(added_node, node, len);
	added_node->updateNeighbor(node, (Node*)1, len);
	added_node->updateNeighbor(dad, (Node*)2, len);

	// now tranverse the tree downwards
	FOR_NEIGHBOR_IT(node, dad, it) {
		Node *target_node2;
		Node *target_dad2;
		double score = addTaxonMP(added_node, target_node2, target_dad2, (*it)->node, node);
		if (score < best_score) {
			best_score = score;
			target_node = target_node2;
			target_dad = target_dad2;
		}
	}
	return best_score;
}


void PhyloTree::growTreeMP(Alignment *alignment) {

	cout << "Stepwise addition using maximum parsimony..." << endl;
	aln = alignment;
	int size = aln->getNSeq();
	if (size < 3)
		outError(ERR_FEW_TAXA);

	root = newNode();
	Node *new_taxon;

	// create initial tree with 3 taxa
	for (leafNum = 0; leafNum < 3; leafNum++)
	{
		if (verbose_mode >= VB_MAX) 
			cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
		new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
		root->addNeighbor(new_taxon, 1.0);
		new_taxon->addNeighbor(root, 1.0);
	}
	root = findNodeID(0);
	//optimizeAllBranches();

	// stepwise adding the next taxon
	for (leafNum = 3; leafNum < size; leafNum++)
	{
		if (verbose_mode >= VB_MAX) 
			cout << "Add " << aln->getSeqName(leafNum) << " to the tree";
		// allocate a new taxon and a new ajedcent internal node
		new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
		Node *added_node = newNode();
		added_node->addNeighbor(new_taxon, 1.0);
		new_taxon->addNeighbor(added_node, 1.0);

		// preserve two neighbors
		added_node->addNeighbor((Node*)1, 1.0);
		added_node->addNeighbor((Node*)2, 1.0);


		Node *target_node = NULL;
		Node *target_dad = NULL;
		int score = addTaxonMP(added_node, target_node, target_dad, root->neighbors[0]->node, root);
		if (verbose_mode >= VB_MAX) 
			cout << ", score = " << score << endl;
		// now insert the new node in the middle of the branch node-dad
		double len = target_dad->findNeighbor(target_node)->length;
		target_node->updateNeighbor(target_dad, added_node, len/2.0);
		target_dad->updateNeighbor(target_node, added_node, len/2.0);
		added_node->updateNeighbor((Node*)1, target_node, len/2.0);
		added_node->updateNeighbor((Node*)2, target_dad, len/2.0);
		// compute the likelihood
		//clearAllPartialLh();
		//optimizeAllBranches();
		//optimizeNNI();
	}

	nodeNum = 2 * leafNum - 2;
}


/****************************************************************************
	likelihood function
****************************************************************************/

void PhyloTree::initializeAllPartialLh() {
	int index;
	initializeAllPartialLh(index);
	assert(index == (nodeNum-1)*2);
}

void PhyloTree::initializeAllPartialLh(int &index, PhyloNode *node, PhyloNode *dad) {
	int block_size = aln->size() * aln->num_states * site_rate->getNRate();
	int pars_block_size = getBitsBlockSize();
	if (!node) {
		node = (PhyloNode*)root;
		// allocate the big central partial likelihoods memory
		if (verbose_mode >= VB_MED)
			cout << "Allocating " << (leafNum-1)*4*block_size*sizeof(double) << " bytes for partial likelihood vectors" << endl;
		if (!central_partial_lh) 
			central_partial_lh = new double[(leafNum-1)*4*block_size];
		if (!central_partial_lh) {
			outError("Not enough memory for partial likelihood vectors");
		}
		if (verbose_mode >= VB_MED)
			cout << "Allocating " << (leafNum-1)*4*pars_block_size*sizeof(UINT) << " bytes for partial parsimony vectors" << endl;
		if (!central_partial_pars)
			central_partial_pars = new UINT[(leafNum-1)*4*pars_block_size];
		if (!central_partial_pars) {
			outError("Not enough memory for partial parsimony vectors");
		}
		index = 0;
	}
	if (dad) {
		// assign a region in central_partial_lh to both Neihgbors (dad->node, and node->dad)
		PhyloNeighbor *nei = (PhyloNeighbor*)node->findNeighbor(dad);
		assert(!nei->partial_lh);
		nei->partial_lh = central_partial_lh + (index * block_size);
		nei->partial_pars = central_partial_pars + (index * pars_block_size);
		nei = (PhyloNeighbor*)dad->findNeighbor(node);
		assert(!nei->partial_lh);
		nei->partial_lh = central_partial_lh + ((index+1) * block_size);
		nei->partial_pars = central_partial_pars + ((index+1) * pars_block_size);
		index += 2;
		assert(index < nodeNum*2-1);
	}
	FOR_NEIGHBOR_IT(node, dad, it)
		initializeAllPartialLh(index, (PhyloNode*)(*it)->node, node);
}



double *PhyloTree::newPartialLh() {
	return new double[aln->size() * aln->num_states * site_rate->getNRate()];
}


double PhyloTree::computeLikelihood() {
	assert(model);
	assert(site_rate);
	return computeLikelihoodBranch((PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root);
}


double PhyloTree::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	PhyloNode *node = (PhyloNode*)dad_branch->node;
	PhyloNeighbor *node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
	assert(node_branch);
	if (!central_partial_lh)
		initializeAllPartialLh();
	// swap node and dad if dad is a leaf
	if (node->isLeaf()) {
		PhyloNode *tmp_node = dad;
		dad = node;
		node = tmp_node;
		PhyloNeighbor *tmp_nei = dad_branch;
		dad_branch = node_branch;
		node_branch = tmp_nei;
		//cout << "swapped\n";
	}
	if ((dad_branch->partial_lh_computed & 1) == 0)
		computePartialLikelihood(dad_branch, dad);
	if ((node_branch->partial_lh_computed & 1) == 0)
		computePartialLikelihood(node_branch, node);
	// now combine likelihood at the branch

	double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
	int ncat = site_rate->getNRate();
	double p_invar = site_rate->getPInvar();
	double p_var_cat = (1.0 - p_invar) / (double)ncat;
	int nstates = aln->num_states;
	int block = ncat * nstates;
	int trans_size = nstates * nstates;
	int ptn, cat, state1, state2;

	double trans_mat[ncat*trans_size];
	double state_freq[nstates];
	model->getStateFrequency(state_freq);

	for (cat = 0; cat < ncat; cat++) {
		//trans_mat[cat] = model->newTransMatrix();
		double *trans_cat = trans_mat + (cat*trans_size);
		model_factory->computeTransMatrix(dad_branch->length * site_rate->getRate(cat), trans_cat);
		for (state1 = 0; state1 < nstates; state1++) {
			double *trans_mat_state = trans_cat + (state1*nstates);
			for (state2 = 0; state2 < nstates; state2++)
				trans_mat_state[state2] *= state_freq[state1];
		}
	}

	for (ptn = 0; ptn < aln->size(); ptn++) {
		double lh_ptn = 0.0; // likelihood of the pattern
		int dad_state = 1000; // just something big enough
		if (dad->isLeaf())
			dad_state = (*aln)[ptn][dad->id];
		for (cat = 0; cat < ncat; cat++) {
			double *partial_lh_site = node_branch->partial_lh + (ptn*block + cat*nstates);
			double *partial_lh_child = dad_branch->partial_lh + (ptn*block + cat*nstates);
			if (dad_state  < nstates) { // single state
				// external node
				double *trans_state = trans_mat + (cat*trans_size + dad_state*nstates);
				for (state2 = 0; state2 < nstates; state2++)
					lh_ptn += partial_lh_child[state2] * trans_state[state2];
			} else {
				// internal node, or external node but ambiguous character
				for (state1 = 0; state1 < nstates; state1++) {
					double lh_state = 0.0;  // likelihood of state1
					double *trans_state = trans_mat + (cat*trans_size + state1*nstates);

                                        // *************  FOR OPTIMIZATION ******************************
					for (state2 = 0; state2 < nstates; state2++)
						lh_state += partial_lh_child[state2] * trans_state[state2];
                                        // *************  FOR OPTIMIZATION ******************************
					lh_ptn += lh_state * partial_lh_site[state1];
				}
			}
		}
		lh_ptn *= p_var_cat;
		if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
			lh_ptn += p_invar * state_freq[(*aln)[ptn][0]];
		}
		assert(lh_ptn > 0);
		tree_lh += log(lh_ptn) * (*aln)[ptn].frequency;
	}
	//for (cat = ncat-1; cat >= 0; cat--)
		//delete trans_mat[cat];
	//delete state_freq;
	return tree_lh;
}

void PhyloTree::computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	// don't recompute the likelihood
	if (dad_branch->partial_lh_computed & 1)
		return;
	Node *node = dad_branch->node;
	int ptn, cat;
	int ncat = site_rate->getNRate();
	int nstates = aln->num_states;
	int block = nstates * site_rate->getNRate();
	int lh_size = aln->size() * block;
	double *partial_lh_site;

	dad_branch->lh_scale_factor = 0.0;

	assert(dad_branch->partial_lh);
	//if (!dad_branch->partial_lh)
	//	dad_branch->partial_lh = newPartialLh();
	if (node->isLeaf() && dad) {
		/* external node */
		memset(dad_branch->partial_lh, 0, lh_size * sizeof(double));
		for (ptn = 0; ptn < aln->size(); ptn++) {
			char state;
			partial_lh_site = dad_branch->partial_lh + (ptn*block);
			if (node->name == ROOT_NAME) {
				state = STATE_UNKNOWN;
			} else {
				assert(node->id < aln->getNSeq());
				state = (aln->at(ptn))[node->id];
			}
			if (state == STATE_UNKNOWN) {
				// fill all entries (also over rate category) with 1.0
				for (int state2 = 0; state2 < block; state2++) {
					partial_lh_site[state2] = 1.0;
				}
			} else if (state < nstates) {
				for (cat = 0; cat < ncat; cat++)
					partial_lh_site[cat*nstates + state] = 1.0;
			} else {
				// ambiguous character, for DNA, RNA
				state = state - (nstates-1);
				for (int state2 = 0; state2 < nstates; state2++)
					if (state & (1 << state2)) {
						for (cat = 0; cat < ncat; cat++)
							partial_lh_site[cat*nstates + state2] = 1.0;
					}
			}
		}
	} else {
		/* internal node */
		double *trans_mat[ncat];
		for (cat = 0; cat < ncat; cat++) trans_mat[cat] = model->newTransMatrix();
		for (ptn = 0; ptn < lh_size; ptn++)
			dad_branch->partial_lh[ptn] = 1.0;
		FOR_NEIGHBOR_IT(node, dad, it)
		if ((*it)->node->name != ROOT_NAME) {
			computePartialLikelihood((PhyloNeighbor*)(*it), (PhyloNode*)node);

			dad_branch->lh_scale_factor += ((PhyloNeighbor*)(*it))->lh_scale_factor;

			for (cat = 0; cat < ncat; cat++)
				model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), trans_mat[cat]);

			for (ptn = 0; ptn < aln->size(); ptn++) {
				for (cat = 0; cat < ncat; cat++)
				{
					partial_lh_site = dad_branch->partial_lh + (ptn*block + cat*nstates);
					double *partial_lh_child = ((PhyloNeighbor*)(*it))->partial_lh + (ptn*block + cat*nstates);
					for (int state = 0; state < nstates; state++) {
						double lh_child = 0.0;
						double *trans_state = trans_mat[cat] + (state * nstates);


                                                // *************  FOR OPTIMIZATION ******************************
						for (int state2 = 0; state2 < nstates; state2++)
							lh_child += trans_state[state2] * partial_lh_child[state2];
                                                // *************  FOR OPTIMIZATION ******************************



						partial_lh_site[state] *= lh_child;
					}
				}
				// check if one should scale partial likelihoods
				bool do_scale = true;
				partial_lh_site = dad_branch->partial_lh + (ptn*block);
				for (cat = 0; cat < block; cat++)
					if (partial_lh_site[cat] > SCALING_THRESHOLD) {
						do_scale = false; break;
					}
				if (!do_scale) continue;
				// now do the likelihood scaling
				/*
				double lh_max = partial_lh_site[0];
				for (cat = 1; cat < block; cat++)
					if (lh_max < partial_lh_site[cat]) lh_max = partial_lh_site[cat];
				for (cat = 0; cat < block; cat++)
					partial_lh_site[cat] /= lh_max;
				dad_branch->lh_scale_factor += log(lh_max) * (*aln)[ptn].frequency;
				
				*/
				for (cat = 0; cat < block; cat++)
					partial_lh_site[cat] /= SCALING_THRESHOLD;
				dad_branch->lh_scale_factor += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
				
			}
		}
		for (cat = ncat-1; cat >= 0; cat--)
			delete [] trans_mat[cat];
	}
	dad_branch->partial_lh_computed |= 1;
}


double PhyloTree::optimizeModel(bool fixed_len) {
	double cur_lh;
	model_factory->stopStoringTransMatrix();
	if (fixed_len) 
		cur_lh = computeLikelihood();
	else {
		cur_lh = optimizeAllBranches(1);
	}
	assert(model);
	assert(site_rate);
	do {
		double model_lh = model->optimizeParameters();
/*
		if (model_lh != 0.0 && !fixed_len)
			model_lh = optimizeAllBranches(3); */
		double rate_lh = site_rate->optimizeParameters();
/*		if (rate_lh != 0.0 && !fixed_len)
			rate_lh = optimizeAllBranches(2);*/
		if (model_lh == 0.0 && rate_lh == 0.0) {
			if (!fixed_len) cur_lh = optimizeAllBranches();
			break;
		}
		double new_lh = (rate_lh != 0.0) ? rate_lh : model_lh;
		if (verbose_mode > VB_MIN) {
			model->writeInfo(cout);
			site_rate->writeInfo(cout);
		}
		if (new_lh > cur_lh + TOL_LIKELIHOOD) {
			if (!fixed_len)
				cur_lh = optimizeAllBranches(5);  // loop only 5 times in total
			else
				cur_lh = new_lh;
			if (verbose_mode > VB_MIN)
				cout << "Current Log-likelihood: " << cur_lh << endl;
		} else {
			if (!fixed_len) cur_lh = optimizeAllBranches();
			break;
		}
	} while (true);
	model_factory->startStoringTransMatrix();
	return cur_lh;
}


/****************************************************************************
	computing derivatives of likelihood function
****************************************************************************/

double PhyloTree::computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf)
{
	PhyloNode *node = (PhyloNode*)dad_branch->node;
	PhyloNeighbor *node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
	assert(node_branch);
	// swap node and dad if dad is a leaf
	if (node->isLeaf()) {
		PhyloNode *tmp_node = dad;
		dad = node;
		node = tmp_node;
		PhyloNeighbor *tmp_nei = dad_branch;
		dad_branch = node_branch;
		node_branch = tmp_nei;
		//cout << "swapped\n";
	}
	if ((dad_branch->partial_lh_computed & 1) == 0)
		computePartialLikelihood(dad_branch, dad);
	if ((node_branch->partial_lh_computed & 1) == 0)
		computePartialLikelihood(node_branch, node);
	// now combine likelihood at the branch

	double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
	df = ddf = 0.0;
	int ncat = site_rate->getNRate();
	double p_invar = site_rate->getPInvar();
	double p_var_cat = (1.0 - p_invar) / (double)ncat;
	int nstates = aln->num_states;
	int block = ncat * nstates;
	int trans_size = nstates * nstates;
	int ptn, cat, state1, state2;

	double trans_mat[ncat*trans_size];
	double trans_derv1[ncat*trans_size];
	double trans_derv2[ncat*trans_size];
	double state_freq[nstates];
	model->getStateFrequency(state_freq);

	for (cat = 0; cat < ncat; cat++) {
		//trans_mat[cat] = model->newTransMatrix();
		double *trans_cat = trans_mat + (cat*trans_size);
		double *derv1_cat = trans_derv1 + (cat*trans_size);
		double *derv2_cat = trans_derv2 + (cat*trans_size);
		double rate_val = site_rate->getRate(cat);
		double rate_sqr = rate_val * rate_val;
		model_factory->computeTransDerv(dad_branch->length * rate_val, trans_cat, derv1_cat, derv2_cat);
		for (state1 = 0; state1 < nstates; state1++) {
			double *trans_mat_state = trans_cat + (state1*nstates);
			double *trans_derv1_state = derv1_cat + (state1*nstates);
			double *trans_derv2_state = derv2_cat + (state1*nstates);

			for (state2 = 0; state2 < nstates; state2++) {
				trans_mat_state[state2] *= state_freq[state1];
				trans_derv1_state[state2] *= (state_freq[state1] * rate_val);
				trans_derv2_state[state2] *= (state_freq[state1] * rate_sqr);
			}
		}
	}

	for (ptn = 0; ptn < aln->size(); ptn++) {
		double lh_ptn = 0.0; // likelihood of the pattern
		double lh_ptn_derv1 = 0.0;
		double lh_ptn_derv2 = 0.0;
		int dad_state = STATE_UNKNOWN;
		if (dad->isLeaf())
			dad_state = (*aln)[ptn][dad->id];
		for (cat = 0; cat < ncat; cat++) {
			double *partial_lh_site = node_branch->partial_lh + (ptn*block + cat*nstates);
			double *partial_lh_child = dad_branch->partial_lh + (ptn*block + cat*nstates);
			if (dad_state  < nstates) {
				// external node
				double *trans_state = trans_mat + (cat*trans_size + dad_state*nstates);
				double *derv1_state = trans_derv1 + (cat*trans_size + dad_state*nstates);
				double *derv2_state = trans_derv2 + (cat*trans_size + dad_state*nstates);
				for (state2 = 0; state2 < nstates; state2++) {
					lh_ptn += partial_lh_child[state2] * trans_state[state2];
					lh_ptn_derv1 += partial_lh_child[state2] * derv1_state[state2];
					lh_ptn_derv2 += partial_lh_child[state2] * derv2_state[state2];
				}
			} else {
				// internal node, or external node but ambiguous character
				for (state1 = 0; state1 < nstates; state1++) {
					double lh_state = 0.0;  // likelihood of state1
					double lh_state_derv1 = 0.0;
					double lh_state_derv2 = 0.0;
					double *trans_state = trans_mat + (cat*trans_size + state1*nstates);
					double *derv1_state = trans_derv1 + (cat*trans_size + state1*nstates);
					double *derv2_state = trans_derv2 + (cat*trans_size + state1*nstates);
					for (state2 = 0; state2 < nstates; state2++) {
						lh_state += partial_lh_child[state2] * trans_state[state2];
						lh_state_derv1 += partial_lh_child[state2] * derv1_state[state2];
						lh_state_derv2 += partial_lh_child[state2] * derv2_state[state2];
					}
					lh_ptn += lh_state * partial_lh_site[state1];
					lh_ptn_derv1 += lh_state_derv1 * partial_lh_site[state1];
					lh_ptn_derv2 += lh_state_derv2 * partial_lh_site[state1];
				}
			}
		}
		lh_ptn *= p_var_cat;
		lh_ptn_derv1 *= p_var_cat;
		lh_ptn_derv2 *= p_var_cat;
		if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
			lh_ptn += p_invar * state_freq[(*aln)[ptn][0]];
		}
		assert(lh_ptn > 0);
		double derv1_frac = lh_ptn_derv1 / lh_ptn;
		double derv2_frac = lh_ptn_derv2 / lh_ptn;
		tree_lh += log(lh_ptn) * (*aln)[ptn].frequency;
		df += derv1_frac * (*aln)[ptn].frequency;
		ddf += (derv2_frac - derv1_frac*derv1_frac) * (*aln)[ptn].frequency;
	}
	//for (cat = ncat-1; cat >= 0; cat--)
		//delete trans_mat[cat];
	//delete state_freq;
	return tree_lh;
}

/****************************************************************************
	Branch length optimization by maximum likelihood
****************************************************************************/

double PhyloTree::computeFunction(double value) {
	current_it->length = value;
	current_it_back->length = value;
	return -computeLikelihoodBranch(current_it, (PhyloNode*)current_it_back->node);
}

double PhyloTree::computeFuncDerv(double value, double &df, double &ddf) {
	current_it->length = value;
	current_it_back->length = value;
	double lh = -computeLikelihoodDerv(current_it, (PhyloNode*)current_it_back->node, df, ddf);
	df = -df;
	ddf = -ddf;
	return lh;
}

double PhyloTree::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH) {
	double negative_lh;
	current_it = (PhyloNeighbor*)node1->findNeighbor(node2);
	assert(current_it);
	current_it_back = (PhyloNeighbor*)node2->findNeighbor(node1);
	assert(current_it_back);
	double current_len = current_it->length;
	double ferror, optx;
	/*if (verbose_mode == VB_DEBUG) {
		cout << "For branch " << node1->name << "," << node2->name << endl;
	}*/
	if (optimize_by_newton) // Newton-Raphson method
		optx = minimizeNewton(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, negative_lh);
	else // Brent method
		optx = minimizeOneDimen(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, &negative_lh, &ferror);
	if (fabs(current_len-optx) <= TOL_BRANCH_LEN) // if nothing changes, return
		return -negative_lh;
	current_it->length = optx;
	current_it_back->length = optx;
        if (clearLH) {
            node1->clearReversePartialLh(node2);
            node2->clearReversePartialLh(node1);
        }
	return -negative_lh;
}


double PhyloTree::optimizeChildBranches(PhyloNode *node, PhyloNode *dad) {
	double tree_lh = 0.0;
	FOR_NEIGHBOR_IT(node, dad, it) {
		tree_lh = optimizeOneBranch((PhyloNode*)node, (PhyloNode*)(*it)->node);
	}
	return tree_lh;
}


double PhyloTree::optimizeAllBranches(PhyloNode *node, PhyloNode *dad) {
	//double tree_lh = optimizeChildBranches(node, dad);
	double tree_lh;
	FOR_NEIGHBOR_IT(node, dad, it) {
		//if (!(*it)->node->isLeaf())
			tree_lh = optimizeAllBranches((PhyloNode*)(*it)->node, node);
	}
	if (dad) tree_lh = optimizeOneBranch(node, dad);
	return tree_lh;
}

double PhyloTree::optimizeAllBranches(int iterations) {
	if (verbose_mode > VB_MED)
		cout << "Optimizing branch lenths (max " << iterations << " loops)..." << endl;
	double tree_lh = computeLikelihood();
	//cout << tree_lh << endl;
	for (int i = 0; i < iterations; i++) {
		double new_tree_lh = optimizeAllBranches((PhyloNode*)root);
		//clearAllPartialLh();
		//new_tree_lh = computeLikelihood();
		if (verbose_mode > VB_MAX) {
			cout << "BRANCH LEN " << i+1 << " : ";
			cout.precision(10);
			cout << new_tree_lh << endl;
		}
		if (new_tree_lh <= tree_lh + TOL_LIKELIHOOD)
			return (new_tree_lh>tree_lh) ? new_tree_lh : tree_lh;
		tree_lh = new_tree_lh;
	}
	return tree_lh;
}


/****************************************************************************
	Stepwise addition (greedy) by maximum likelihood
****************************************************************************/

double PhyloTree::addTaxonML(Node *added_node, Node* &target_node, Node* &target_dad, Node *node, Node *dad) {
	Neighbor *dad_nei = dad->findNeighbor(node);

	// now insert the new node in the middle of the branch node-dad
	double len = dad_nei->length;
	node->updateNeighbor(dad, added_node, len/2.0);
	dad->updateNeighbor(node, added_node, len/2.0);
	added_node->updateNeighbor((Node*)1, node, len/2.0);
	added_node->updateNeighbor((Node*)2, dad, len/2.0);
	// compute the likelihood
	clearAllPartialLh();
	double best_score = optimizeChildBranches((PhyloNode*)added_node);
	target_node = node;
	target_dad = dad;
	// remove the added node
	node->updateNeighbor(added_node, dad, len);
	dad->updateNeighbor(added_node, node, len);
	added_node->updateNeighbor(node, (Node*)1, len);
	added_node->updateNeighbor(dad, (Node*)2, len);

	// now tranverse the tree downwards
	FOR_NEIGHBOR_IT(node, dad, it) {
		Node *target_node2;
		Node *target_dad2;
		double score = addTaxonML(added_node, target_node2, target_dad2, (*it)->node, node);
		if (score > best_score) {
			best_score = score;
			target_node = target_node2;
			target_dad = target_dad2;
		}
	}
	return best_score;
}


void PhyloTree::growTreeML(Alignment *alignment) {

	cout << "Stepwise addition using ML..." << endl;
	aln = alignment;
	int size = aln->getNSeq();
	if (size < 3)
		outError(ERR_FEW_TAXA);

	root = newNode();
	Node *new_taxon;

	// create initial tree with 3 taxa
	for (leafNum = 0; leafNum < 3; leafNum++)
	{
		cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
		new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
		root->addNeighbor(new_taxon, 1.0);
		new_taxon->addNeighbor(root, 1.0);
	}
	root = findNodeID(0);
	optimizeAllBranches();

	// stepwise adding the next taxon
	for (leafNum = 3; leafNum < size; leafNum++)
	{
		cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
		// allocate a new taxon and a new ajedcent internal node
		new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
		Node *added_node = newNode();
		added_node->addNeighbor(new_taxon, 1.0);
		new_taxon->addNeighbor(added_node, 1.0);

		// preserve two neighbors
		added_node->addNeighbor((Node*)1, 1.0);
		added_node->addNeighbor((Node*)2, 1.0);


		Node *target_node = NULL;
		Node *target_dad = NULL;
		addTaxonML(added_node, target_node, target_dad, root->neighbors[0]->node, root);
		// now insert the new node in the middle of the branch node-dad
		double len = target_dad->findNeighbor(target_node)->length;
		target_node->updateNeighbor(target_dad, added_node, len/2.0);
		target_dad->updateNeighbor(target_node, added_node, len/2.0);
		added_node->updateNeighbor((Node*)1, target_node, len/2.0);
		added_node->updateNeighbor((Node*)2, target_dad, len/2.0);
		// compute the likelihood
		clearAllPartialLh();
		optimizeAllBranches();
		optimizeNNI();
	}

	nodeNum = 2 * leafNum - 2;
}


/****************************************************************************
	compute BioNJ tree, a more accurate extension of Neighbor-Joining
****************************************************************************/

void PhyloTree::computeBioNJ(Params &params, Alignment *alignment, double* &dist_mat) {
	cout << "Computing BioNJ tree..." << endl;
	string dist_file = params.out_prefix;
	string bionj_file = params.out_prefix;
	dist_file += ".dist";
	bionj_file += ".bionj";

	if (!dist_mat) {
		dist_mat = new double[alignment->getNSeq() * alignment->getNSeq()];
	}
	if (!params.dist_file)
		alignment->computeDist(dist_mat);
	else
		alignment->readDist(params.dist_file, dist_mat);

	alignment->printDist(dist_file.c_str(), dist_mat);
	//delete dist_mat;

	BioNj bionj;
	bionj.create(dist_file.c_str(), bionj_file.c_str());
	bool my_rooted = false;
	readTree(bionj_file.c_str(), my_rooted);
	setAlignment(alignment);
}

int PhyloTree::fixNegativeBranch(double fixed_length, Node *node, Node *dad) {
	if (!node) node = root;
	int fixed = 0;
	FOR_NEIGHBOR_IT(node, dad, it) {
		if ((*it)->length < 0.0) { // negative branch length detected
			if (verbose_mode == VB_DEBUG)
				cout << "Negative branch length " << (*it)->length << " was set to ";
			//(*it)->length = ((double)(rand())/RAND_MAX)*0.1+TOL_BRANCH_LEN;
			(*it)->length = fixed_length;				
			if (verbose_mode == VB_DEBUG)
				cout << (*it)->length << endl;
			// set the backward branch length
			(*it)->node->findNeighbor(node)->length = (*it)->length;
			fixed++;
		}
		fixed += fixNegativeBranch(fixed_length, (*it)->node, node);
	}
	return fixed;
}


/****************************************************************************
	Nearest Neighbor Interchange by maximum likelihood
****************************************************************************/

double PhyloTree::swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2) {
	assert(node1->degree() == 3 && node2->degree() == 3);

	PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
	PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);
	double node12_len = node12_it->length;

	// save the likelihood vector at the two ends of node1-node2
	double *node1_lh_save = node12_it->partial_lh;
	double *node2_lh_save = node21_it->partial_lh;
	node12_it->partial_lh = newPartialLh();
	node21_it->partial_lh = newPartialLh();

	// TUNG save the first found neighbor (2 Neighbor total) of node 1 (excluding node2) in node1_it
	FOR_NEIGHBOR_DECLARE(node1, node2, node1_it)
		break;
	Neighbor *node1_nei = *node1_it;
	double node1_len = node1_nei->length;
	NeighborVec::iterator node1_nei_it = node1_nei->node->findNeighborIt(node1);


	// Neighbors of node2 which are not node1
	FOR_NEIGHBOR_IT(node2, node1, node2_it) {

		// do the NNI swap
		Neighbor *node2_nei = *node2_it;
		// TUNG unused variable ?
		NeighborVec::iterator node2_nei_it = node2_nei->node->findNeighborIt(node2);

		double node2_len = node2_nei->length;

		node1->updateNeighbor(node1_it, node2_nei);
		node2_nei->node->updateNeighbor(node2, node1);

		node2->updateNeighbor(node2_it, node1_nei);
		node1_nei->node->updateNeighbor(node1, node2);

		// partial_lhclear partial likelihood vector
		node12_it->clearPartialLh();
		node21_it->clearPartialLh();

		// compute the score of the swapped topology
		double score;
		score = optimizeOneBranch(node1, node2);
		// if better: return
		if (score > cur_score) {
			node2->clearReversePartialLh(node1);
			node1->clearReversePartialLh(node2);
			cur_score = score;
			cout << "Swapped neighbors :" << node1_nei->node->id << " and " << node2_nei->node->id << endl;
			break;

		}
		//cout << "Could not find better score for NNI with Node " << node1->id << "->" << "Node " << node2->id << endl;

		// else, swap back, also recover the branch lengths
		node1->updateNeighbor(node1_it, node1_nei, node1_len);
		node1_nei->node->updateNeighbor(node2, node1, node1_len);
		node2->updateNeighbor(node2_it, node2_nei, node2_len);
		node2_nei->node->updateNeighbor(node1, node2, node2_len);
		node12_it->length = node12_len;
		node21_it->length = node12_len;
	}

	// restore the partial likelihood vector
	delete [] node21_it->partial_lh;
	delete [] node12_it->partial_lh;
	node12_it->partial_lh = node1_lh_save;
	node21_it->partial_lh = node2_lh_save;
	return cur_score;
}

double PhyloTree::optimizeNNI(double cur_score, PhyloNode *node, PhyloNode *dad)  {
	if (!node) node = (PhyloNode*)root;
	if (!node->isLeaf() && dad && !dad->isLeaf()) {
		double score = swapNNIBranch(cur_score, node, dad);
		if (score > cur_score) return score;
	}

	FOR_NEIGHBOR_IT(node, dad, it) {
		double score = optimizeNNI(cur_score, (PhyloNode*)(*it)->node, node);
		if (score > cur_score) return score;
	}
	return cur_score;
}

double PhyloTree::optimizeNNI() {
	double cur_score = computeLikelihood();
	for (int i = 0; i < 100; i++) {
		double score = optimizeNNI(cur_score);
		if (score <= cur_score) break;
		if (verbose_mode > VB_MED)
			cout << "NNI " << i+1 << " : " << score << endl;
		//cur_score = score;
		cur_score=optimizeAllBranches((PhyloNode*)root);
	}
	return optimizeAllBranches();
}

double PhyloTree::optimizeNNIBranches() {
	if (verbose_mode > VB_MIN)
		cout << "Search with Nearest Neighbor Interchange (NNI) using ML..." << endl;
	double cur_score = computeLikelihood();
	for (int i = 0; i < 100; i++) {
		double score = optimizeNNI();
		if (score <= cur_score + TOL_LIKELIHOOD) break;
		cur_score = score;
	}
	return cur_score;
}

/****************************************************************************
	Subtree Pruning and Regrafting by maximum likelihood
****************************************************************************/

double PhyloTree::optimizeSPR(double cur_score, PhyloNode *node, PhyloNode *dad) {
	if (!node) node = (PhyloNode*)root;
	PhyloNeighbor *dad1_nei = NULL;
	PhyloNeighbor *dad2_nei = NULL;
	PhyloNode *sibling1 = NULL;
	PhyloNode *sibling2 = NULL;
	double sibling1_len = 0.0, sibling2_len = 0.0;

	if (dad && !dad->isLeaf()) {
		assert(dad->degree() == 3);
		// assign the sibling of node, with respect to dad
		FOR_NEIGHBOR_DECLARE(dad, node, it) {
			if (!sibling1) {
				dad1_nei = (PhyloNeighbor*)(*it);
				sibling1 = (PhyloNode*)(*it)->node;
				sibling1_len = (*it)->length;
			} else {
				dad2_nei = (PhyloNeighbor*)(*it);
				sibling2 = (PhyloNode*)(*it)->node;
				sibling2_len = (*it)->length;
			}
		}
		// remove the subtree leading to node
		double sum_len = sibling1_len + sibling2_len;
		sibling1->updateNeighbor(dad, sibling2, sum_len);
		sibling2->updateNeighbor(dad, sibling1, sum_len);
		PhyloNeighbor* sibling1_nei = (PhyloNeighbor*)sibling1->findNeighbor(sibling2);
		PhyloNeighbor* sibling2_nei = (PhyloNeighbor*)sibling2->findNeighbor(sibling1);
		sibling1_nei->clearPartialLh();
		sibling2_nei->clearPartialLh();

		// now try to move the subtree to somewhere else
		vector<PhyloNeighbor*> spr_path;
		FOR_NEIGHBOR(sibling1, sibling2, it) {
			spr_path.push_back(sibling1_nei);
			double score = swapSPR(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*)(*it)->node, sibling1, spr_path);
			// if likelihood score improves, return
			if (score > cur_score) return score;
			spr_path.pop_back();
		}
		FOR_NEIGHBOR(sibling2, sibling1, it) {
			spr_path.push_back(sibling2_nei);
			double score = swapSPR(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*)(*it)->node, sibling2, spr_path);
			// if likelihood score improves, return
			if (score > cur_score) return score;
			spr_path.pop_back();
		}
		// if likelihood does not imporve, swap back
		sibling1->updateNeighbor(sibling2, dad, sibling1_len);
		sibling2->updateNeighbor(sibling1, dad, sibling2_len);
		dad1_nei->node = sibling1;
		dad1_nei->length = sibling1_len;
		dad2_nei->node = sibling2;
		dad2_nei->length = sibling2_len;
		clearAllPartialLh();
	}
	FOR_NEIGHBOR_IT(node, dad, it) {
		double score = optimizeSPR(cur_score, (PhyloNode*)(*it)->node, node);
		if (score > cur_score) return score;
	}
	return cur_score;
}

double PhyloTree::swapSPR(double cur_score, int cur_depth, PhyloNode *node1, PhyloNode *dad1,
	PhyloNode *orig_node1, PhyloNode *orig_node2,
	PhyloNode *node2, PhyloNode *dad2, vector<PhyloNeighbor*> &spr_path) {

	PhyloNeighbor *node1_nei = (PhyloNeighbor*)node1->findNeighbor(dad1);
	PhyloNeighbor *dad1_nei = (PhyloNeighbor*)dad1->findNeighbor(node1);
	double node1_dad1_len = node1_nei->length;
	PhyloNeighbor *node2_nei = (PhyloNeighbor*)node2->findNeighbor(dad2);

	if (dad2) {
		// now, connect (node1-dad1) to the branch (node2-dad2)
		bool first = true;
		PhyloNeighbor *node2_nei = (PhyloNeighbor*)node2->findNeighbor(dad2);
		PhyloNeighbor *dad2_nei = (PhyloNeighbor*)dad2->findNeighbor(node2);
		double len2 = node2_nei->length;
		FOR_NEIGHBOR_IT(dad1, node1, it) {
			if (first) {
				(*it)->node = dad2;
				(*it)->length = len2 / 2;
				dad2->updateNeighbor(node2, dad1, len2 / 2);
				first = false;
			} else {
				(*it)->node = node2;
				(*it)->length = len2 / 2;
				node2->updateNeighbor(dad2, dad1, len2 / 2);
			}
			((PhyloNeighbor*)(*it))->clearPartialLh();
		}
		node2_nei->clearPartialLh();
		dad2_nei->clearPartialLh();
		node1_nei->clearPartialLh();
		vector<PhyloNeighbor*>::iterator it2;
		for (it2 = spr_path.begin(); it2 != spr_path.end(); it2++)
			(*it2)->clearPartialLh();
		clearAllPartialLh();
		// optimize relevant branches
		double score;

		/* testing different branch optimization */
		score = optimizeOneBranch(node1, dad1);
		//score = optimizeOneBranch(dad2, dad1);
		//score = optimizeOneBranch(node2, dad1);
		/*
		PhyloNode *cur_node = dad2;
		for (int i = spr_path.size()-1; i >= 0; i--) {
			score = optimizeOneBranch(cur_node, (PhyloNode*)spr_path[i]->node);
			cur_node = (PhyloNode*)spr_path[i]->node;
		}
		*/
		//score = optimizeAllBranches(dad1);

		// if score improves, return
		if (score > cur_score) return score;
		// else, swap back
		node2->updateNeighbor(dad1, dad2, len2);
		dad2->updateNeighbor(dad1, node2, len2);
		node2_nei->clearPartialLh();
		dad2_nei->clearPartialLh();
		node1_nei->length = node1_dad1_len;
		dad1_nei->length = node1_dad1_len;

		// add to candiate SPR moves
		spr_moves.add(node1, dad1, node2, dad2, score);
	}
	if (cur_depth >= spr_radius) return cur_score;
	spr_path.push_back(node2_nei);
	FOR_NEIGHBOR_IT(node2, dad2, it) {
		double score = swapSPR(cur_score, cur_depth+1, node1, dad1, orig_node1, orig_node2, (PhyloNode*)(*it)->node, node2, spr_path);
		if (score > cur_score) return score;
	}
	spr_path.pop_back();
	return cur_score;
}

double PhyloTree::assessSPRMove(double cur_score, const SPRMove &spr) {
	PhyloNode *dad = spr.prune_dad;
	PhyloNode *node = spr.prune_node;
	PhyloNode *dad2 = spr.regraft_dad;
	PhyloNode *node2 = spr.regraft_node;

	PhyloNeighbor *dad_nei1 = NULL;
	PhyloNeighbor *dad_nei2 = NULL;
	PhyloNode *sibling1 = NULL;
	PhyloNode *sibling2 = NULL;
	double sibling1_len = 0.0, sibling2_len = 0.0;

	PhyloNeighbor *node1_nei = (PhyloNeighbor*)node->findNeighbor(dad);
	PhyloNeighbor *dad1_nei = (PhyloNeighbor*)dad->findNeighbor(node);
	double node1_dad1_len = node1_nei->length;

	// assign the sibling of node, with respect to dad
	FOR_NEIGHBOR_DECLARE(dad, node, it) {
		if (!sibling1) {
			dad_nei1 = (PhyloNeighbor*)(*it);
			sibling1 = (PhyloNode*)(*it)->node;
			sibling1_len = (*it)->length;
		} else {
			dad_nei2 = (PhyloNeighbor*)(*it);
			sibling2 = (PhyloNode*)(*it)->node;
			sibling2_len = (*it)->length;
		}
	}
	// remove the subtree leading to node
	double sum_len = sibling1_len + sibling2_len;
	sibling1->updateNeighbor(dad, sibling2, sum_len);
	sibling2->updateNeighbor(dad, sibling1, sum_len);
	// now try to move the subtree to somewhere else

	bool first = true;
	PhyloNeighbor *node2_nei = (PhyloNeighbor*)node2->findNeighbor(dad2);
	PhyloNeighbor *dad2_nei = (PhyloNeighbor*)dad2->findNeighbor(node2);
	double len2 = node2_nei->length;
	FOR_NEIGHBOR(dad, node, it) {
		if (first) {
			(*it)->node = dad2;
			(*it)->length = len2 / 2;
			dad2->updateNeighbor(node2, dad, len2 / 2);
			first = false;
		} else {
			(*it)->node = node2;
			(*it)->length = len2 / 2;
			node2->updateNeighbor(dad2, dad, len2 / 2);
		}
		((PhyloNeighbor*)(*it))->clearPartialLh();
	}

	clearAllPartialLh();
	// optimize branches
	double score;
	score = optimizeAllBranches(dad);

	// if score improves, return
	if (score > cur_score) return score;
	// else, swap back
	node2->updateNeighbor(dad, dad2, len2);
	dad2->updateNeighbor(dad, node2, len2);

	node1_nei->length = node1_dad1_len;
	dad1_nei->length = node1_dad1_len;

	sibling1->updateNeighbor(sibling2, dad, sibling1_len);
	sibling2->updateNeighbor(sibling1, dad, sibling2_len);
	dad_nei1->node = sibling1;
	dad_nei1->length = sibling1_len;
	dad_nei2->node = sibling2;
	dad_nei2->length = sibling2_len;
	clearAllPartialLh();

	return cur_score;

}

double PhyloTree::optimizeSPR() {
	double cur_score = computeLikelihood();
	//spr_radius = leafNum / 5;
	spr_radius = 2;
	for (int i = 0; i < 100; i++) {
		spr_moves.clear();
		double score = optimizeSPR(cur_score, (PhyloNode*)root->neighbors[0]->node);
		clearAllPartialLh();
		if (score <= cur_score) {
			for (SPRMoves::iterator it = spr_moves.begin(); it != spr_moves.end(); it++) {
				//cout << (*it).score << endl;
				score = assessSPRMove(cur_score, *it);
				// if likelihood score improves, apply to SPR
				if (score > cur_score) break;
			}
			if (score <= cur_score) break;
		}
		cout.precision(10);
		cout << "SPR " << i+1 << " : " << score << endl;
		//cur_score = score;
		cur_score = optimizeAllBranches();
	}
	return optimizeAllBranches();
}

double PhyloTree::optimizeSPRBranches() {
	cout << "Search with Subtree Pruning and Regrafting (SPR) using ML..." << endl;
	double cur_score = computeLikelihood();
	for (int i = 0; i < 100; i++) {
		double score = optimizeSPR();
		if (score <= cur_score + TOL_LIKELIHOOD) break;
		cur_score = score;
	}
	return cur_score;
}


void PhyloTree::pruneSubtree(PhyloNode *node, PhyloNode *dad, PruningInfo &info) {
	bool first = true;
	info.node = node;
	info.dad = dad;
	FOR_NEIGHBOR_IT(dad, node, it) {
		if (first) {
			info.dad_it_left = it;
			info.dad_nei_left = (*it);
			info.dad_lh_left = ((PhyloNeighbor*)(*it))->partial_lh;
			info.left_node = (*it)->node;
			info.left_len = (*it)->length;
			first = false;
		} else {
			info.dad_it_right = it;
			info.dad_nei_right = (*it);
			info.dad_lh_right = ((PhyloNeighbor*)(*it))->partial_lh;
			info.right_node = (*it)->node;
			info.right_len = (*it)->length;
		}
	}
	info.left_it = info.left_node->findNeighborIt(dad);
	info.right_it = info.right_node->findNeighborIt(dad);
	info.left_nei = (*info.left_it);
	info.right_nei = (*info.right_it);

	info.left_node->updateNeighbor(info.left_it, info.dad_nei_right);
	info.right_node->updateNeighbor(info.right_it, info.dad_nei_left);
	((PhyloNeighbor*)info.dad_nei_right)->partial_lh = newPartialLh();
	((PhyloNeighbor*)info.dad_nei_left)->partial_lh = newPartialLh();
}

void PhyloTree::regraftSubtree(PruningInfo &info,
	PhyloNode *in_node, PhyloNode *in_dad) {
	NeighborVec::iterator in_node_it = in_node->findNeighborIt(in_dad);
	NeighborVec::iterator in_dad_it = in_dad->findNeighborIt(in_node);
	Neighbor *in_dad_nei = (*in_dad_it);
	Neighbor *in_node_nei = (*in_node_it);
	double in_len = in_dad_nei->length;
	info.dad->updateNeighbor(info.dad_it_right, in_dad_nei);
	info.dad->updateNeighbor(info.dad_it_left, in_node_nei);
	// SOMETHING NEED TO BE DONE
	//in_dad->updateNeighbor(in_dad_it,

}
