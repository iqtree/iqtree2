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
#include "iqtree.h"
#include "phylosupertree.h"
#include "phylosupertreeplen.h"
#include "mexttree.h"
#include "timeutil.h"
#include "gtrmodel.h"
#include "rategamma.h"
#include <numeric>
#include "pll/pllInternal.h"
#include "nnisearch.h"
#include "vectorclass/vectorclass.h"
#include "vectorclass/vectormath_common.h"


Params *globalParam;
Alignment *globalAlignment;
extern StringIntMap pllTreeCounter;


IQTree::IQTree() : PhyloTree() {
    init();
}

void IQTree::init() {
    k_represent = 0;
    k_delete = k_delete_min = k_delete_max = k_delete_stay = 0;
    dist_matrix = NULL;
    var_matrix = NULL;
    nni_count_est = 0.0;
    nni_delta_est = 0;
    curScore = 0.0; // Current score of the tree
    bestScore = -DBL_MAX; // Best score found so far
    curIt = 1;
    cur_pars_score = -1;
//    enable_parsimony = false;
    estimate_nni_cutoff = false;
    nni_cutoff = -1e6;
    nni_sort = false;
    testNNI = false;
    print_tree_lh = false;
    write_intermediate_trees = 0;
    max_candidate_trees = 0;
    logl_cutoff = 0.0;
    len_scale = 10000;
    save_all_br_lens = false;
    duplication_counter = 0;
    pllInst = NULL;
    pllAlignment = NULL;
    pllPartitions = NULL;
    //boot_splits = new SplitGraph;
    pll2iqtree_pattern_index = NULL;
    fastNNI = true;
}

IQTree::IQTree(Alignment *aln) : PhyloTree(aln) {
    init();
}
void IQTree::setParams(Params &params) {
    searchinfo.speednni = params.speednni;
    searchinfo.nni_type = params.nni_type;
    optimize_by_newton = params.optimize_by_newton;
    candidateTrees.aln = aln;
    candidateTrees.max_candidates = params.popSize;
    candidateTrees.limit = params.limitPopSize;

    sse = params.SSE;
//    if (params.maxtime != 1000000) {
//        params.autostop = false;
//    }
    if (params.min_iterations == -1) {
        if (!params.gbo_replicates) {
            if (params.stop_condition == SC_UNSUCCESS_ITERATION) {
                params.min_iterations = aln->getNSeq() * 100;
            } else if (aln->getNSeq() < 100) {
                params.min_iterations = 200;
            } else {
                params.min_iterations = aln->getNSeq() * 2;
            }
            if (params.iteration_multiple > 1)
                params.min_iterations = aln->getNSeq() * params.iteration_multiple;
        } else {
            params.min_iterations = 100;
        }
    }
    if (params.gbo_replicates)
        params.max_iterations = max(params.max_iterations, max(params.min_iterations, 1000));

    k_represent = params.k_representative;

    if (params.p_delete == -1.0) {
        if (aln->getNSeq() < 4)
            params.p_delete = 0.0; // delete nothing
        else if (aln->getNSeq() == 4)
            params.p_delete = 0.25; // just delete 1 leaf
        else if (aln->getNSeq() == 5)
            params.p_delete = 0.4; // just delete 2 leaves
        else if (aln->getNSeq() < 51)
            params.p_delete = 0.5;
        else if (aln->getNSeq() < 100)
            params.p_delete = 0.3;
        else if (aln->getNSeq() < 200)
            params.p_delete = 0.2;
        else if (aln->getNSeq() < 400)
            params.p_delete = 0.1;
        else
            params.p_delete = 0.05;
    }
    //tree.setProbDelete(params.p_delete);
    if (params.p_delete != -1.0) {
        k_delete = k_delete_min = k_delete_max = ceil(params.p_delete * leafNum);
    } else {
        k_delete = k_delete_min = 10;
        k_delete_max = leafNum / 2;
        if (k_delete_max > 100)
            k_delete_max = 100;
        if (k_delete_max < 20)
            k_delete_max = 20;
        k_delete_stay = ceil(leafNum / k_delete);
    }

    //tree.setIQPIterations(params.stop_condition, params.stop_confidence, params.min_iterations, params.max_iterations);

    stop_rule.initialize(params);

    //tree.setIQPAssessQuartet(params.iqp_assess_quartet);
    iqp_assess_quartet = params.iqp_assess_quartet;
    estimate_nni_cutoff = params.estimate_nni_cutoff;
    nni_cutoff = params.nni_cutoff;
    nni_sort = params.nni_sort;
    testNNI = params.testNNI;

    this->params = &params;
    globalParam = &params;
    globalAlignment = aln;

    write_intermediate_trees = params.write_intermediate_trees;

    if (write_intermediate_trees > 2 || params.gbo_replicates > 0) {
        save_all_trees = 1;
    }
    if (params.gbo_replicates > 0) {
        if (params.iqp_assess_quartet != IQP_BOOTSTRAP) {
            save_all_trees = 2;
        }
    }
    if (params.gbo_replicates > 0 && params.do_compression)
        save_all_br_lens = true;
    print_tree_lh = params.print_tree_lh;
    max_candidate_trees = params.max_candidate_trees;
    if (max_candidate_trees == 0)
        max_candidate_trees = aln->getNSeq() * params.step_iterations;
    setRootNode(params.root);

    string bootaln_name = params.out_prefix;
    bootaln_name += ".bootaln";
    if (params.print_bootaln) {
        ofstream bootalnout;
    	bootalnout.open(bootaln_name.c_str());
    	bootalnout.close();
    }
    size_t i;

    if (params.online_bootstrap && params.gbo_replicates > 0) {
        cout << "Generating " << params.gbo_replicates << " samples for ultrafast bootstrap..." << endl;
        // allocate memory for boot_samples
        boot_samples.resize(params.gbo_replicates);
#ifdef BOOT_VAL_FLOAT
        size_t nptn = get_safe_upper_limit_float(getAlnNPattern());
#else
        size_t nptn = get_safe_upper_limit(getAlnNPattern());
#endif
        BootValType *mem = aligned_alloc<BootValType>(nptn * (size_t)(params.gbo_replicates));
        memset(mem, 0, nptn * (size_t)(params.gbo_replicates) * sizeof(BootValType));
        for (i = 0; i < params.gbo_replicates; i++)
        	boot_samples[i] = mem + i*nptn;

        boot_logl.resize(params.gbo_replicates, -DBL_MAX);
        boot_trees.resize(params.gbo_replicates, -1);
        boot_counts.resize(params.gbo_replicates, 0);
        VerboseMode saved_mode = verbose_mode;
        verbose_mode = VB_QUIET;
        for (i = 0; i < params.gbo_replicates; i++) {
        	if (params.print_bootaln) {
    			Alignment* bootstrap_alignment;
    			if (aln->isSuperAlignment())
    				bootstrap_alignment = new SuperAlignment;
    			else
    				bootstrap_alignment = new Alignment;
    			IntVector this_sample;
    			bootstrap_alignment->createBootstrapAlignment(aln, &this_sample, params.bootstrap_spec);
    			for (size_t j = 0; j < nptn; j++)
    				boot_samples[i][j] = this_sample[j];
				bootstrap_alignment->printPhylip(bootaln_name.c_str(), true);
				delete bootstrap_alignment;
        	} else {
    			IntVector this_sample;
        		aln->createBootstrapAlignment(this_sample, params.bootstrap_spec);
    			for (size_t j = 0; j < nptn; j++)
    				boot_samples[i][j] = this_sample[j];
        	}
        }
        verbose_mode = saved_mode;
        if (params.print_bootaln) {
        	cout << "Bootstrap alignments printed to " << bootaln_name << endl;
        }

        cout << "Max candidate trees (tau): " << max_candidate_trees << endl;
    }

    if (params.root_state) {
        if (strlen(params.root_state) != 1)
            outError("Root state must have exactly 1 character");
        root_state = aln->convertState(params.root_state[0]);
        if (root_state < 0 || root_state >= aln->num_states)
            outError("Invalid root state");
    }
}

void myPartitionsDestroy(partitionList *pl) {
	int i;
	for (i = 0; i < pl->numberOfPartitions; i++) {
		rax_free(pl->partitionData[i]->partitionName);
		rax_free(pl->partitionData[i]);
	}
	rax_free(pl->partitionData);
	rax_free(pl);
}

IQTree::~IQTree() {
    //if (bonus_values)
    //delete bonus_values;
    //bonus_values = NULL;
    if (dist_matrix)
        delete[] dist_matrix;
    dist_matrix = NULL;

    if (var_matrix)
        delete[] var_matrix;
    var_matrix = NULL;

    for (vector<double*>::reverse_iterator it = treels_ptnlh.rbegin(); it != treels_ptnlh.rend(); it++)
        delete[] (*it);
    treels_ptnlh.clear();
    for (vector<SplitGraph*>::reverse_iterator it2 = boot_splits.rbegin(); it2 != boot_splits.rend(); it2++)
        delete (*it2);
    //if (boot_splits) delete boot_splits;
    if (pllPartitions)
    	myPartitionsDestroy(pllPartitions);
    if (pllAlignment)
    	pllAlignmentDataDestroy(pllAlignment);
    if (pllInst)
        pllDestroyInstance(pllInst);

    if (!boot_samples.empty())
    	aligned_free(boot_samples[0]); // free memory
}

void IQTree::createPLLPartition(Params &params, ostream &pllPartitionFileHandle) {
    if (isSuperTree()) {
        PhyloSuperTree *siqtree = (PhyloSuperTree*) this;
        // additional check for stupid PLL hard limit
        if (siqtree->size() > PLL_NUM_BRANCHES)
        	outError("Number of partitions exceeds PLL limit, please increase PLL_NUM_BRANCHES constant in pll.h");
        int i = 0;
        int startPos = 1;
        for (PhyloSuperTree::iterator it = siqtree->begin(); it != siqtree->end(); it++) {
            i++;
            int curLen = ((*it))->getAlnNSite();
            if ((*it)->aln->seq_type == SEQ_DNA) {
                pllPartitionFileHandle << "DNA";
            } else if ((*it)->aln->seq_type == SEQ_PROTEIN) {
            	if (siqtree->part_info[i-1].model_name != "" && siqtree->part_info[i-1].model_name.substr(0, 4) != "TEST")
            		pllPartitionFileHandle << siqtree->part_info[i-1].model_name.substr(0, siqtree->part_info[i-1].model_name.find_first_of("+{"));
            	else
            		pllPartitionFileHandle << "WAG";
            } else
            	outError("PLL only works with DNA/protein alignments");
            pllPartitionFileHandle << ", p" << i << " = " << startPos << "-" << startPos + curLen - 1 << endl;
            startPos = startPos + curLen;
        }
    } else {
        /* create a partition file */
        string model;
        if (aln->seq_type == SEQ_DNA) {
            model = "DNA";
        } else if (aln->seq_type == SEQ_PROTEIN) {
        	if (params.pll && params.model_name != "" && params.model_name.substr(0, 4) != "TEST") {
        		model = params.model_name.substr(0, params.model_name.find_first_of("+{"));
        	} else {
        		model = "WAG";
        	}
        } else {
        	model = "WAG";
        	//outError("PLL currently only supports DNA/protein alignments");
        }
        pllPartitionFileHandle << model << ", p1 = " << "1-" << getAlnNSite() << endl;
    }
}

void IQTree::initializePLL(Params &params) {
    pllAttr.rateHetModel = PLL_GAMMA;
    pllAttr.fastScaling = PLL_FALSE;
    pllAttr.saveMemory = PLL_FALSE;
    pllAttr.useRecom = PLL_FALSE;
    pllAttr.randomNumberSeed = params.ran_seed;
#ifdef _OPENMP
    pllAttr.numberOfThreads = params.num_threads; /* This only affects the pthreads version */
#else
    pllAttr.numberOfThreads = 1;
#endif
    if (pllInst != NULL) {
        pllDestroyInstance(pllInst);
    }
    /* Create a PLL instance */
    pllInst = pllCreateInstance(&pllAttr);

    /* Read in the alignment file */
    stringstream pllAln;
	if (aln->isSuperAlignment()) {
		((SuperAlignment*) aln)->printCombinedAlignment(pllAln);
	} else {
		aln->printPhylip(pllAln);
	}
	string pllAlnStr = pllAln.str();
    pllAlignment = pllParsePHYLIPString(pllAlnStr.c_str(), pllAlnStr.length());

    /* Read in the partition information */
    // BQM: to avoid printing file
    stringstream pllPartitionFileHandle;
    createPLLPartition(params, pllPartitionFileHandle);
    pllQueue *partitionInfo = pllPartitionParseString(pllPartitionFileHandle.str().c_str());

    /* Validate the partitions */
    if (!pllPartitionsValidate(partitionInfo, pllAlignment)) {
        outError("pllPartitionsValidate");
    }

    /* Commit the partitions and build a partitions structure */
    pllPartitions = pllPartitionsCommit(partitionInfo, pllAlignment);

    /* We don't need the the intermediate partition queue structure anymore */
    pllQueuePartitionsDestroy(&partitionInfo);

    /* eliminate duplicate sites from the alignment and update weights vector */
    pllAlignmentRemoveDups(pllAlignment, pllPartitions);

    pllTreeInitTopologyForAlignment(pllInst, pllAlignment);

    /* Connect the alignment and partition structure with the tree structure */
    if (!pllLoadAlignment(pllInst, pllAlignment, pllPartitions, PLL_SHALLOW_COPY)) {
        outError("Incompatible tree/alignment combination");
    }
}


void IQTree::initializeModel(Params &params) {
    try {
        if (!getModelFactory()) {
            if (isSuperTree()) {
                if (params.partition_type) {
                    setModelFactory(new PartitionModelPlen(params, (PhyloSuperTreePlen*) this));
                } else
                    setModelFactory(new PartitionModel(params, (PhyloSuperTree*) this));
            } else {
                setModelFactory(new ModelFactory(params, this));
            }
        }
    } catch (string & str) {
        outError(str);
    }
    setModel(getModelFactory()->model);
    setRate(getModelFactory()->site_rate);

    if (params.pll) {
        if (getRate()->getNDiscreteRate() == 1) {
        	outError("Non-Gamma model is not yet supported by PLL.");
            // TODO: change rateHetModel to PLL_CAT in case of non-Gamma model
        }
        if (getRate()->name.substr(0,2) == "+I")
        	outError("+Invar model is not yet supported by PLL.");
        if (aln->seq_type == SEQ_DNA && getModel()->name != "GTR")
        	outError("non GTR model for DNA is not yet supported by PLL.");
    }

}
double IQTree::getProbDelete() {
    return (double) k_delete / leafNum;
}

void IQTree::resetKDelete() {
    k_delete = k_delete_min;
    k_delete_stay = ceil(leafNum / k_delete);
}

void IQTree::increaseKDelete() {
    if (k_delete >= k_delete_max)
        return;
    k_delete_stay--;
    if (k_delete_stay > 0)
        return;
    k_delete++;
    k_delete_stay = ceil(leafNum / k_delete);
    if (verbose_mode >= VB_MED)
        cout << "Increase k_delete to " << k_delete << endl;
}

//void IQTree::setIQPIterations(STOP_CONDITION stop_condition, double stop_confidence, int min_iterations,
//        int max_iterations) {
//    stop_rule.setStopCondition(stop_condition);
//    stop_rule.setConfidenceValue(stop_confidence);
//    stop_rule.setIterationNum(min_iterations, max_iterations);
//}

RepresentLeafSet* IQTree::findRepresentLeaves(vector<RepresentLeafSet*> &leaves_vec, int nei_id, PhyloNode *dad) {
    PhyloNode *node = (PhyloNode*) (dad->neighbors[nei_id]->node);
    int set_id = dad->id * 3 + nei_id;
    if (leaves_vec[set_id])
        return leaves_vec[set_id];
    RepresentLeafSet *leaves = new RepresentLeafSet;
    RepresentLeafSet * leaves_it[2] = { NULL, NULL };
    leaves_vec[set_id] = leaves;
    RepresentLeafSet::iterator last;
    RepresentLeafSet::iterator cur_it;
    int i, j;
    //double admit_height = 1000000;

    leaves->clear();
    if (node->isLeaf()) {
        // set the depth to zero
        //node->height = 0.0;
        leaves->insert(new RepLeaf(node, 0));
    } else {
        for (i = 0, j = 0; i < node->neighbors.size(); i++)
            if (node->neighbors[i]->node != dad) {
                leaves_it[j++] = findRepresentLeaves(leaves_vec, i, node);
            }
        assert(j == 2 && leaves_it[0] && leaves_it[1]);
        if (leaves_it[0]->empty() && leaves_it[1]->empty()) {
            cout << "wrong";
        }
        RepresentLeafSet::iterator lit[2] = { leaves_it[0]->begin(), leaves_it[1]->begin() };
        while (leaves->size() < k_represent) {
            int id = -1;
            if (lit[0] != leaves_it[0]->end() && lit[1] != leaves_it[1]->end()) {
                if ((*lit[0])->height < (*lit[1])->height)
                    id = 0;
                else if ((*lit[0])->height > (*lit[1])->height)
                    id = 1;
                else { // tie, choose at random
                    id = random_int(2);
                }
            } else if (lit[0] != leaves_it[0]->end())
                id = 0;
            else if (lit[1] != leaves_it[1]->end())
                id = 1;
            else
                break;
            assert(id < 2 && id >= 0);
            leaves->insert(new RepLeaf((*lit[id])->leaf, (*lit[id])->height + 1));
            lit[id]++;
        }
    }
    assert(!leaves->empty());
    /*
     if (verbose_mode >= VB_MAX) {
     for (cur_it = leaves->begin(); cur_it != leaves->end(); cur_it++)
     cout << (*cur_it)->leaf->name << " ";
     cout << endl;
     }*/
    return leaves;
}

/*
 void IQPTree::clearRepresentLeaves(vector<RepresentLeafSet*> &leaves_vec, Node *node, Node *dad) {
 int nei_id;
 for (nei_id = 0; nei_id < node->neighbors.size(); nei_id++)
 if (node->neighbors[nei_id]->node == dad) break;
 assert(nei_id < node->neighbors.size());
 int set_id = node->id * 3 + nei_id;
 if (leaves_vec[set_id]) {
 for (RepresentLeafSet::iterator rlit = leaves_vec[set_id]->begin(); rlit != leaves_vec[set_id]->end(); rlit++)
 delete (*rlit);
 delete leaves_vec[set_id];
 leaves_vec[set_id] = NULL;
 }
 FOR_NEIGHBOR_IT(node, dad, it) {
 clearRepresentLeaves(leaves_vec, (*it)->node, node);
 }
 }*/

void IQTree::deleteNonCherryLeaves(PhyloNodeVector &del_leaves) {
    NodeVector cherry_taxa;
    NodeVector noncherry_taxa;
    // get the vector of non cherry taxa
    getNonCherryLeaves(noncherry_taxa, cherry_taxa);
    root = NULL;
    int num_taxa = aln->getNSeq();
    int num_delete = k_delete;
    if (num_delete > num_taxa - 4)
        num_delete = num_taxa - 4;
    if (verbose_mode >= VB_DEBUG) {
        cout << "Deleting " << num_delete << " leaves" << endl;
    }
    vector<unsigned int> indices_noncherry(noncherry_taxa.size());
    //iota(indices_noncherry.begin(), indices_noncherry.end(), 0);
    unsigned int startValue = 0;
    for (vector<unsigned int>::iterator it = indices_noncherry.begin(); it != indices_noncherry.end(); ++it) {
        (*it) = startValue;
        ++startValue;
    }
    random_shuffle(indices_noncherry.begin(), indices_noncherry.end());
    int i;
    for (i = 0; i < num_delete && i < noncherry_taxa.size(); i++) {
        PhyloNode *taxon = (PhyloNode*) noncherry_taxa[indices_noncherry[i]];
        del_leaves.push_back(taxon);
        deleteLeaf(taxon);
        //cout << taxon->id << ", ";
    }
    int j = 0;
    if (i < num_delete) {
        vector<unsigned int> indices_cherry(cherry_taxa.size());
        //iota(indices_cherry.begin(), indices_cherry.end(), 0);
        startValue = 0;
        for (vector<unsigned int>::iterator it = indices_cherry.begin(); it != indices_cherry.end(); ++it) {
            (*it) = startValue;
            ++startValue;
        }
        random_shuffle(indices_cherry.begin(), indices_cherry.end());
        while (i < num_delete) {
            PhyloNode *taxon = (PhyloNode*) cherry_taxa[indices_cherry[j]];
            del_leaves.push_back(taxon);
            deleteLeaf(taxon);
            i++;
            j++;
        }
    }
    root = cherry_taxa[j];
}

void IQTree::deleteLeaves(PhyloNodeVector &del_leaves) {
    NodeVector taxa;
    // get the vector of taxa
    getTaxa(taxa);
    root = NULL;
    //int num_delete = floor(p_delete * taxa.size());
    int num_delete = k_delete;
    int i;
    if (num_delete > taxa.size() - 4)
        num_delete = taxa.size() - 4;
    if (verbose_mode >= VB_DEBUG) {
        cout << "Deleting " << num_delete << " leaves" << endl;
    }
    // now try to randomly delete some taxa of the probability of p_delete
    for (i = 0; i < num_delete;) {
        int id = random_int(taxa.size());
        if (!taxa[id])
            continue;
        else
            i++;
        PhyloNode *taxon = (PhyloNode*) taxa[id];
        del_leaves.push_back(taxon);
        deleteLeaf(taxon);
        taxa[id] = NULL;
    }
    // set root to the first taxon which was not deleted
    for (i = 0; i < taxa.size(); i++)
        if (taxa[i]) {
            root = taxa[i];
            break;
        }
}

int IQTree::assessQuartet(Node *leaf0, Node *leaf1, Node *leaf2, Node *del_leaf) {
    assert(dist_matrix);
    int nseq = aln->getNSeq();
    //int id0 = leaf0->id, id1 = leaf1->id, id2 = leaf2->id;
    double dist0 = dist_matrix[leaf0->id * nseq + del_leaf->id] + dist_matrix[leaf1->id * nseq + leaf2->id];
    double dist1 = dist_matrix[leaf1->id * nseq + del_leaf->id] + dist_matrix[leaf0->id * nseq + leaf2->id];
    double dist2 = dist_matrix[leaf2->id * nseq + del_leaf->id] + dist_matrix[leaf0->id * nseq + leaf1->id];
    if (dist0 < dist1 && dist0 < dist2)
        return 0;
    if (dist1 < dist2)
        return 1;
    return 2;
}

int IQTree::assessQuartetParsimony(Node *leaf0, Node *leaf1, Node *leaf2, Node *del_leaf) {
    int score[3] = { 0, 0, 0 };
    for (Alignment::iterator it = aln->begin(); it != aln->end(); it++) {
        char ch0 = (*it)[leaf0->id];
        char ch1 = (*it)[leaf1->id];
        char ch2 = (*it)[leaf2->id];
        char chd = (*it)[del_leaf->id];
        if (ch0 >= aln->num_states || ch1 >= aln->num_states || ch2 >= aln->num_states || chd >= aln->num_states)
            continue;
        if (chd == ch0 && ch1 == ch2)
            score[0] += (*it).frequency;
        if (chd == ch1 && ch0 == ch2)
            score[1] += (*it).frequency;
        if (chd == ch2 && ch0 == ch1)
            score[2] += (*it).frequency;
    }
    if (score[0] == score[1] && score[0] == score[2]) {
        int id = random_int(3);
        return id;
    }
    if (score[0] > score[1] && score[0] > score[2])
        return 0;
    if (score[1] < score[2])
        return 2;
    return 1;
}

void IQTree::initializeBonus(PhyloNode *node, PhyloNode *dad) {
    if (!node)
        node = (PhyloNode*) root;
    if (dad) {
        PhyloNeighbor *node_nei = (PhyloNeighbor*) node->findNeighbor(dad);
        PhyloNeighbor *dad_nei = (PhyloNeighbor*) dad->findNeighbor(node);
        node_nei->lh_scale_factor = 0.0;
        node_nei->partial_lh_computed = 0;
        dad_nei->lh_scale_factor = 0.0;
        dad_nei->partial_lh_computed = 0;
    }

    FOR_NEIGHBOR_IT(node, dad, it){
    initializeBonus((PhyloNode*) ((*it)->node), node);
}
}

void IQTree::raiseBonus(Neighbor *nei, Node *dad, double bonus) {
    ((PhyloNeighbor*) nei)->lh_scale_factor += bonus;
    if (verbose_mode >= VB_DEBUG)
        cout << dad->id << " - " << nei->node->id << " : " << bonus << endl;

    //  FOR_NEIGHBOR_IT(nei->node, dad, it)
    //	raiseBonus((*it), nei->node, bonus);
}

double IQTree::computePartialBonus(Node *node, Node* dad) {
    PhyloNeighbor *node_nei = (PhyloNeighbor*) node->findNeighbor(dad);
    if (node_nei->partial_lh_computed)
        return node_nei->lh_scale_factor;

    FOR_NEIGHBOR_IT(node, dad, it){
    node_nei->lh_scale_factor += computePartialBonus((*it)->node, node);
}
    node_nei->partial_lh_computed = 1;
    return node_nei->lh_scale_factor;
}


/*
 bool IQTree::containPosNNI(vector<NNIMove> posNNIs) {
 for (vector<NNIMove>::iterator iter = posNNIs.begin(); iter != posNNIs.end(); iter++) {
 if (iter->newloglh > iter->oldloglh)
 return true;
 }
 return false;
 }
 */

void IQTree::findBestBonus(double &best_score, NodeVector &best_nodes, NodeVector &best_dads, Node *node, Node *dad) {
    double score;
    if (!node)
        node = root;
    if (!dad) {
        best_score = 0;
    } else {
        score = computePartialBonus(node, dad) + computePartialBonus(dad, node);
        if (score >= best_score) {
            if (score > best_score) {
                best_score = score;
                best_nodes.clear();
                best_dads.clear();
            }
            best_nodes.push_back(node);
            best_dads.push_back(dad);
        }
        //cout << node->id << " - " << dad->id << " : " << best_score << endl;
    }

    FOR_NEIGHBOR_IT(node, dad, it){
    findBestBonus(best_score, best_nodes, best_dads, (*it)->node, node);
}
}

void IQTree::assessQuartets(vector<RepresentLeafSet*> &leaves_vec, PhyloNode *cur_root, PhyloNode *del_leaf) {
    const int MAX_DEGREE = 3;
    RepresentLeafSet * leaves[MAX_DEGREE];
    double bonus[MAX_DEGREE];
    memset(bonus, 0, MAX_DEGREE * sizeof(double));
    int cnt = 0;

    // only work for birfucating tree
    assert(cur_root->degree() == MAX_DEGREE);

    // find the representative leaf set for three subtrees

    FOR_NEIGHBOR_IT(cur_root, NULL, it){
    leaves[cnt] = findRepresentLeaves(leaves_vec, cnt, cur_root);
    cnt++;
}
    for (RepresentLeafSet::iterator i0 = leaves[0]->begin(); i0 != leaves[0]->end(); i0++)
        for (RepresentLeafSet::iterator i1 = leaves[1]->begin(); i1 != leaves[1]->end(); i1++)
            for (RepresentLeafSet::iterator i2 = leaves[2]->begin(); i2 != leaves[2]->end(); i2++) {
                int best_id;
                if (iqp_assess_quartet == IQP_DISTANCE)
                    best_id = assessQuartet((*i0)->leaf, (*i1)->leaf, (*i2)->leaf, del_leaf);
                else
                    best_id = assessQuartetParsimony((*i0)->leaf, (*i1)->leaf, (*i2)->leaf, del_leaf);
                bonus[best_id] += 1.0;
            }
    for (cnt = 0; cnt < MAX_DEGREE; cnt++)
        if (bonus[cnt] > 0.0)
            raiseBonus(cur_root->neighbors[cnt], cur_root, bonus[cnt]);

}

void IQTree::reinsertLeavesByParsimony(PhyloNodeVector &del_leaves) {
    PhyloNodeVector::iterator it_leaf;
    assert(root->isLeaf());
    for (it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); it_leaf++) {
        //cout << "Add leaf " << (*it_leaf)->id << " to the tree" << endl;
        initializeAllPartialPars();
        clearAllPartialLH();
        Node *target_node = NULL;
        Node *target_dad = NULL;
        Node *added_node = (*it_leaf)->neighbors[0]->node;
        Node *node1 = NULL;
        Node *node2 = NULL;
        //Node *leaf;
        for (int i = 0; i < 3; i++) {
            if (added_node->neighbors[i]->node->id == (*it_leaf)->id) {
                //leaf = added_node->neighbors[i]->node;
            } else if (!node1) {
                node1 = added_node->neighbors[i]->node;
            } else {
                node2 = added_node->neighbors[i]->node;
            }
        }

        //cout << "(" << node1->id << ", " << node2->id << ")" << "----" << "(" << added_node->id << "," << leaf->id << ")" << endl;
        added_node->updateNeighbor(node1, (Node*) 1);
        added_node->updateNeighbor(node2, (Node*) 2);

        addTaxonMPFast(added_node, target_node, target_dad, root->neighbors[0]->node, root);
        target_node->updateNeighbor(target_dad, added_node, -1.0);
        target_dad->updateNeighbor(target_node, added_node, -1.0);
        added_node->updateNeighbor((Node*) 1, target_node, -1.0);
        added_node->updateNeighbor((Node*) 2, target_dad, -1.0);

    }

}

void IQTree::reinsertLeaves(PhyloNodeVector &del_leaves) {
    PhyloNodeVector::iterator it_leaf;

    //int num_del_leaves = del_leaves.size();
    assert(root->isLeaf());

    for (it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); it_leaf++) {
        if (verbose_mode >= VB_DEBUG)
            cout << "Reinserting " << (*it_leaf)->name << " (" << (*it_leaf)->id << ")" << endl;
        vector<RepresentLeafSet*> leaves_vec;
        leaves_vec.resize(nodeNum * 3, NULL);
        initializeBonus();
        NodeVector nodes;
        getInternalNodes(nodes);
        if (verbose_mode >= VB_DEBUG)
            drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
        //printTree(cout, WT_BR_LEN | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
        for (NodeVector::iterator it = nodes.begin(); it != nodes.end(); it++) {
            assessQuartets(leaves_vec, (PhyloNode*) (*it), (*it_leaf));
        }
        NodeVector best_nodes, best_dads;
        double best_bonus;
        findBestBonus(best_bonus, best_nodes, best_dads);
        if (verbose_mode >= VB_DEBUG)
            cout << "Best bonus " << best_bonus << " " << best_nodes[0]->id << " " << best_dads[0]->id << endl;
        assert(best_nodes.size() == best_dads.size());
        int node_id = random_int(best_nodes.size());
        if (best_nodes.size() > 1 && verbose_mode >= VB_DEBUG)
            cout << best_nodes.size() << " branches show the same best bonus, branch nr. " << node_id << " is chosen"
                    << endl;

        reinsertLeaf(*it_leaf, best_nodes[node_id], best_dads[node_id]);
        //clearRepresentLeaves(leaves_vec, *it_node, *it_leaf);
        /*if (verbose_mode >= VB_DEBUG) {
         printTree(cout);
         cout << endl;
         }*/
        for (vector<RepresentLeafSet*>::iterator rit = leaves_vec.begin(); rit != leaves_vec.end(); rit++)
            if ((*rit)) {
                RepresentLeafSet *tit = (*rit);
                for (RepresentLeafSet::iterator rlit = tit->begin(); rlit != tit->end(); rlit++)
                    delete (*rlit);
                delete (*rit);
            }
    }
    initializeTree(); // BQM: re-index nodes and branches s.t. ping-pong neighbors have the same ID

    if (verbose_mode >= VB_DEBUG)
        drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
}

void IQTree::doParsimonyReinsertion() {
    PhyloNodeVector del_leaves;

    deleteLeaves(del_leaves);

    reinsertLeavesByParsimony(del_leaves);
    fixNegativeBranch(false);
}

void IQTree::setBestTree(string treeString, double treeLogl) {
    bestTreeString = treeString;
    bestScore = treeLogl;
}

void IQTree::doRandomNNIs(int numNNI) {
    map<int, Node*> usedNodes;
    NodeVector nodeList1, nodeList2;
    getInternalBranches(nodeList1, nodeList2);
    int numInBran = nodeList1.size();
    assert(numInBran == aln->getNSeq() - 3);
    for (int i = 0; i < numNNI; i++) {
        int index = random_int(numInBran);
        if (usedNodes.find(nodeList1[index]->id) == usedNodes.end()
                && usedNodes.find(nodeList2[index]->id) == usedNodes.end()) {
            doOneRandomNNI(nodeList1[index], nodeList2[index]);
            usedNodes.insert(map<int, Node*>::value_type(nodeList1[index]->id, nodeList1[index]));
            usedNodes.insert(map<int, Node*>::value_type(nodeList2[index]->id, nodeList2[index]));
        } else {
            usedNodes.clear();
            nodeList1.clear();
            nodeList2.clear();
            getInternalBranches(nodeList1, nodeList2);
            doOneRandomNNI(nodeList1[index], nodeList2[index]);
            usedNodes.insert(map<int, Node*>::value_type(nodeList1[index]->id, nodeList1[index]));
            usedNodes.insert(map<int, Node*>::value_type(nodeList2[index]->id, nodeList2[index]));
        }
    }
}

void IQTree::doIQP() {
    if (verbose_mode >= VB_DEBUG)
        drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
    //double time_begin = getCPUTime();
    PhyloNodeVector del_leaves;
    deleteLeaves(del_leaves);
    reinsertLeaves(del_leaves);

    // just to make sure IQP does it right
    setAlignment(aln);

//    if (enable_parsimony) {
//        cur_pars_score = computeParsimony();
//        if (verbose_mode >= VB_MAX) {
//            cout << "IQP Likelihood = " << curScore << "  Parsimony = " << cur_pars_score << endl;
//        }
//    }
}

double IQTree::inputTree2PLL(string treestring, bool computeLH) {
    double res = 0.0;
    // read in the tree string from IQTree kernel
    pllNewickTree *newick = pllNewickParseString(treestring.c_str());
    pllTreeInitTopologyNewick(pllInst, newick, PLL_FALSE);
    pllNewickParseDestroy(&newick);
    if (computeLH) {
        pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_TRUE, PLL_FALSE);
        res = pllInst->likelihood;
    }
    return res;
}

double* IQTree::getModelRatesFromPLL() {
    assert(aln->num_states == 4);
    int numberOfRates = (pllPartitions->partitionData[0]->states * pllPartitions->partitionData[0]->states
            - pllPartitions->partitionData[0]->states) / 2;
    double* rate_params = new double[numberOfRates];
    for (int i = 0; i < numberOfRates; i++) {
        rate_params[i] = pllPartitions->partitionData[0]->substRates[i];
    }
    return rate_params;
}

void IQTree::pllPrintModelParams() {
    cout.precision(6);
    cout << fixed;
    for (int part = 0; part < pllPartitions->numberOfPartitions; part++) {
        cout << "Alpha[" << part << "]" << ": " << pllPartitions->partitionData[part]->alpha << endl;
        if (aln->num_states == 4) {
            int states, rates;
            states = pllPartitions->partitionData[part]->states;
            rates = ((states * states - states) / 2);
            cout << "Rates[" << part << "]: " << " ac ag at cg ct gt: ";
            for (int i = 0; i < rates; i++) {
                cout << pllPartitions->partitionData[part]->substRates[i] << " ";
            }
            cout << endl;
            cout <<  "Frequencies: ";
            for (int i = 0; i < 4; i++) {
                cout << pllPartitions->partitionData[part]->empiricalFrequencies[i] << " ";
            }
            cout << endl;
        }
    }
    cout.precision(3);
    cout << fixed;
}

double IQTree::getAlphaFromPLL() {
    return pllPartitions->partitionData[0]->alpha;
}

void IQTree::inputModelPLL2IQTree() {
    // TODO add support for partitioned model
    dynamic_cast<RateGamma*>(getRate())->setGammaShape(pllPartitions->partitionData[0]->alpha);
    if (aln->num_states == 4) {
        ((GTRModel*) getModel())->setRateMatrix(pllPartitions->partitionData[0]->substRates);
        getModel()->decomposeRateMatrix();
    }
    ((GTRModel*) getModel())->setStateFrequency(pllPartitions->partitionData[0]->empiricalFrequencies);
}

void IQTree::inputModelIQTree2PLL() {
    // get the alpha parameter
    double alpha = getRate()->getGammaShape();
    if (alpha == 0.0)
        alpha = PLL_ALPHA_MAX;
    if (aln->num_states == 4) {
        // get the rate parameters
        double *rate_param = new double[6];
        getModel()->getRateMatrix(rate_param);
        // get the state frequencies
        double *state_freqs = new double[aln->num_states];
        getModel()->getStateFrequency(state_freqs);

        /* put them into PLL */
        stringstream linkagePattern;
        int partNr;
        for (partNr = 0; partNr < pllPartitions->numberOfPartitions - 1; partNr++) {
            linkagePattern << partNr << ",";
        }
        linkagePattern << partNr;
        char *pattern = new char[linkagePattern.str().length() + 1];
        strcpy(pattern, linkagePattern.str().c_str());
        pllLinkAlphaParameters(pattern, pllPartitions);
        pllLinkFrequencies(pattern, pllPartitions);
        pllLinkRates(pattern, pllPartitions);
        delete[] pattern;

        for (partNr = 0; partNr < pllPartitions->numberOfPartitions; partNr++) {
            pllSetFixedAlpha(alpha, partNr, pllPartitions, pllInst);
            pllSetFixedBaseFrequencies(state_freqs, 4, partNr, pllPartitions, pllInst);
            pllSetFixedSubstitutionMatrix(rate_param, 6, partNr, pllPartitions, pllInst);
        }
        delete[] rate_param;
        delete[] state_freqs;
    } else if (aln->num_states == 20) {
        double *state_freqs = new double[aln->num_states];
        getModel()->getStateFrequency(state_freqs);
        int partNr;
        for (partNr = 0; partNr < pllPartitions->numberOfPartitions; partNr++) {
            pllSetFixedAlpha(alpha, partNr, pllPartitions, pllInst);
            pllSetFixedBaseFrequencies(state_freqs, 20, partNr, pllPartitions, pllInst);
        }
        delete[] state_freqs;
    } else {
        if (params->pll) {
            outError("Phylogenetic likelihood library current does not support data type other than DNA or Protein");
        }
    }
}

void IQTree::pllBuildIQTreePatternIndex(){
    pll2iqtree_pattern_index = new int[pllAlignment->sequenceLength];
    char ** pll_aln = new char *[pllAlignment->sequenceCount];
    for(int i = 0; i < pllAlignment->sequenceCount; i++)
        pll_aln[i] = new char[pllAlignment->sequenceLength];

    int pos;
    for(int i = 0; i < pllAlignment->sequenceCount; i++){
        pos = 0;
        for(int model = 0; model < pllPartitions->numberOfPartitions; model++){
            memcpy(&pll_aln[i][pos],
                    &pllAlignment->sequenceData[i + 1][pllPartitions->partitionData[model]->lower],
                    pllPartitions->partitionData[model]->width);
            pos += pllPartitions->partitionData[model]->width;
        }
    }

	char * pll_site = new char[pllAlignment->sequenceCount + 1];
	char * site = new char[pllAlignment->sequenceCount + 1];
    for(int i = 0; i < pllAlignment->sequenceLength; i++){
        for(int j = 0; j < pllAlignment->sequenceCount; j++)
            pll_site[j]= pll_aln[j][i];
        pll_site[pllAlignment->sequenceCount] = '\0';

        site[pllAlignment->sequenceCount] = '\0';
        for(int k = 0; k < aln->size(); k++){
            for(int p = 0; p < pllAlignment->sequenceCount; p++)
                site[p] = aln->convertStateBack(aln->at(k)[p]);
            pllBaseSubstitute(site, pllPartitions->partitionData[0]->dataType);
            if(memcmp(pll_site,site, pllAlignment->sequenceCount) == 0){
                pll2iqtree_pattern_index[i] = k;
            }
        }
    }

    delete [] pll_site;
    delete [] site;
    for(int i = 0; i < pllAlignment->sequenceCount; i++)
        delete [] pll_aln[i];
    delete [] pll_aln;
}


/**
 * DTH:
 * Substitute bases in seq according to PLL's rules
 * This function should be updated if PLL's rules change.
 * @param seq: data of some sequence to be substituted
 * @param dataType: PLL_DNA_DATA or PLL_AA_DATA
 */
void IQTree::pllBaseSubstitute (char *seq, int dataType)
{
    char meaningDNA[256];
    char  meaningAA[256];
    char * d;

    for (int i = 0; i < 256; ++ i)
    {
        meaningDNA[i] = -1;
        meaningAA[i]  = -1;
    }

    /* DNA data */

    meaningDNA[(int)'A'] =  1;
    meaningDNA[(int)'B'] = 14;
    meaningDNA[(int)'C'] =  2;
    meaningDNA[(int)'D'] = 13;
    meaningDNA[(int)'G'] =  4;
    meaningDNA[(int)'H'] = 11;
    meaningDNA[(int)'K'] = 12;
    meaningDNA[(int)'M'] =  3;
    meaningDNA[(int)'R'] =  5;
    meaningDNA[(int)'S'] =  6;
    meaningDNA[(int)'T'] =  8;
    meaningDNA[(int)'U'] =  8;
    meaningDNA[(int)'V'] =  7;
    meaningDNA[(int)'W'] =  9;
    meaningDNA[(int)'Y'] = 10;
    meaningDNA[(int)'a'] =  1;
    meaningDNA[(int)'b'] = 14;
    meaningDNA[(int)'c'] =  2;
    meaningDNA[(int)'d'] = 13;
    meaningDNA[(int)'g'] =  4;
    meaningDNA[(int)'h'] = 11;
    meaningDNA[(int)'k'] = 12;
    meaningDNA[(int)'m'] =  3;
    meaningDNA[(int)'r'] =  5;
    meaningDNA[(int)'s'] =  6;
    meaningDNA[(int)'t'] =  8;
    meaningDNA[(int)'u'] =  8;
    meaningDNA[(int)'v'] =  7;
    meaningDNA[(int)'w'] =  9;
    meaningDNA[(int)'y'] = 10;

    meaningDNA[(int)'N'] =
    meaningDNA[(int)'n'] =
    meaningDNA[(int)'O'] =
    meaningDNA[(int)'o'] =
    meaningDNA[(int)'X'] =
    meaningDNA[(int)'x'] =
    meaningDNA[(int)'-'] =
    meaningDNA[(int)'?'] = 15;

    /* AA data */

    meaningAA[(int)'A'] =  0;  /* alanine */
    meaningAA[(int)'R'] =  1;  /* arginine */
    meaningAA[(int)'N'] =  2;  /*  asparagine*/
    meaningAA[(int)'D'] =  3;  /* aspartic */
    meaningAA[(int)'C'] =  4;  /* cysteine */
    meaningAA[(int)'Q'] =  5;  /* glutamine */
    meaningAA[(int)'E'] =  6;  /* glutamic */
    meaningAA[(int)'G'] =  7;  /* glycine */
    meaningAA[(int)'H'] =  8;  /* histidine */
    meaningAA[(int)'I'] =  9;  /* isoleucine */
    meaningAA[(int)'L'] =  10; /* leucine */
    meaningAA[(int)'K'] =  11; /* lysine */
    meaningAA[(int)'M'] =  12; /* methionine */
    meaningAA[(int)'F'] =  13; /* phenylalanine */
    meaningAA[(int)'P'] =  14; /* proline */
    meaningAA[(int)'S'] =  15; /* serine */
    meaningAA[(int)'T'] =  16; /* threonine */
    meaningAA[(int)'W'] =  17; /* tryptophan */
    meaningAA[(int)'Y'] =  18; /* tyrosine */
    meaningAA[(int)'V'] =  19; /* valine */
    meaningAA[(int)'B'] =  20; /* asparagine, aspartic 2 and 3*/
    meaningAA[(int)'Z'] =  21; /*21 glutamine glutamic 5 and 6*/
    meaningAA[(int)'a'] =  0;  /* alanine */
    meaningAA[(int)'r'] =  1;  /* arginine */
    meaningAA[(int)'n'] =  2;  /*  asparagine*/
    meaningAA[(int)'d'] =  3;  /* aspartic */
    meaningAA[(int)'c'] =  4;  /* cysteine */
    meaningAA[(int)'q'] =  5;  /* glutamine */
    meaningAA[(int)'e'] =  6;  /* glutamic */
    meaningAA[(int)'g'] =  7;  /* glycine */
    meaningAA[(int)'h'] =  8;  /* histidine */
    meaningAA[(int)'i'] =  9;  /* isoleucine */
    meaningAA[(int)'l'] =  10; /* leucine */
    meaningAA[(int)'k'] =  11; /* lysine */
    meaningAA[(int)'m'] =  12; /* methionine */
    meaningAA[(int)'f'] =  13; /* phenylalanine */
    meaningAA[(int)'p'] =  14; /* proline */
    meaningAA[(int)'s'] =  15; /* serine */
    meaningAA[(int)'t'] =  16; /* threonine */
    meaningAA[(int)'w'] =  17; /* tryptophan */
    meaningAA[(int)'y'] =  18; /* tyrosine */
    meaningAA[(int)'v'] =  19; /* valine */
    meaningAA[(int)'b'] =  20; /* asparagine, aspartic 2 and 3*/
    meaningAA[(int)'z'] =  21; /*21 glutamine glutamic 5 and 6*/

    meaningAA[(int)'X'] =
    meaningAA[(int)'x'] =
    meaningAA[(int)'?'] =
    meaningAA[(int)'*'] =
    meaningAA[(int)'-'] = 22;

    d = (dataType == PLL_DNA_DATA) ? meaningDNA : meaningAA;
    int seq_len = strlen(seq);
    for (int i = 0; i < seq_len; ++ i)
    {
        seq[i] = d[(int)seq[i]];
    }
}

double IQTree::swapTaxa(PhyloNode *node1, PhyloNode *node2) {
    assert(node1->isLeaf());
    assert(node2->isLeaf());

    PhyloNeighbor *node1nei = (PhyloNeighbor*) *(node1->neighbors.begin());
    PhyloNeighbor *node2nei = (PhyloNeighbor*) *(node2->neighbors.begin());

    node2nei->node->updateNeighbor(node2, node1);
    node1nei->node->updateNeighbor(node1, node2);

    // Update the new neightbors of the 2 nodes
    node1->updateNeighbor(node1->neighbors.begin(), node2nei);
    node2->updateNeighbor(node2->neighbors.begin(), node1nei);

    PhyloNeighbor *node1NewNei = (PhyloNeighbor*) *(node1->neighbors.begin());
    PhyloNeighbor *node2NewNei = (PhyloNeighbor*) *(node2->neighbors.begin());

    // Reoptimize the branch lengths
    optimizeOneBranch(node1, (PhyloNode*) node1NewNei->node);
    this->curScore = optimizeOneBranch(node2, (PhyloNode*) node2NewNei->node);
    //drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    return this->curScore;
}

double IQTree::perturb(int times) {
    while (times > 0) {
        NodeVector taxa;
        // get the vector of taxa
        getTaxa(taxa);
        int taxonid1 = random_int(taxa.size());
        PhyloNode *taxon1 = (PhyloNode*) taxa[taxonid1];
        PhyloNode *taxon2;
        int *dists = new int[taxa.size()];
        int minDist = 1000000;
        for (int i = 0; i < taxa.size(); i++) {
            if (i == taxonid1)
                continue;
            taxon2 = (PhyloNode*) taxa[i];
            int dist = taxon1->calDist(taxon2);
            dists[i] = dist;
            if (dist >= 7 && dist < minDist)
                minDist = dist;
        }

        int taxonid2 = -1;
        for (int i = 0; i < taxa.size(); i++) {
            if (dists[i] == minDist)
                taxonid2 = i;
        }

        taxon2 = (PhyloNode*) taxa[taxonid2];

        cout << "Swapping node " << taxon1->id << " and node " << taxon2->id << endl;
        cout << "Distance " << minDist << endl;
        curScore = swapTaxa(taxon1, taxon2);
        //taxa.erase( taxa.begin() + taxaID1 );
        //taxa.erase( taxa.begin() + taxaID2 -1 );

        times--;
        delete[] dists;
    }
    curScore = optimizeAllBranches(1);
    return curScore;
}

extern "C" pllUFBootData * pllUFBootDataPtr;

string IQTree::optimizeModelParameters(bool printInfo) {
	string newTree;
	if (params->pll) {
		cout << "Estimate model parameters (epsilon = " << params->modeps << ")" << endl;
		double stime = getCPUTime();
		pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_FALSE,
				PLL_FALSE);
		pllOptimizeModelParameters(pllInst, pllPartitions, params->modeps);
		curScore = pllInst->likelihood;
		double etime = getCPUTime();
		cout << etime - stime << " seconds (logl: " << curScore << ")" << endl;
		pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions,
				pllInst->start->back, PLL_TRUE,
				PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
				PLL_FALSE, PLL_FALSE);
		if (printInfo) {
			pllPrintModelParams();
		}
		newTree = string(pllInst->tree_string);
	} else {
		string curTree = getTreeString();
//		double *rate_param_bk = NULL;
//		if (aln->num_states == 4) {
//			rate_param_bk = new double[6];
//			getModel()->getRateMatrix(rate_param_bk);
//		}
//		double alpha_bk = getRate()->getGammaShape();
		cout << "Estimate model parameters (epsilon = " << params->modeps << ")" << endl;
		double modOptScore = getModelFactory()->optimizeParameters(params->fixed_branch_length, printInfo, params->modeps);
		if (isSuperTree()) {
			((PhyloSuperTree*) this)->computeBranchLengths();
		}

		if (modOptScore < curScore - 1.0) {
			cout << "  BUG: Tree logl gets worse after model optimization!" << endl;
			cout << "  Old logl: " << curScore << " / " << "new logl: " << modOptScore << endl;
			abort();
//			readTreeString(curTree);
//			initializeAllPartialLh();
//			newTree = curTree;
//			clearAllPartialLH();
//			if (aln->num_states == 4) {
//				assert(rate_param_bk != NULL);
//				((GTRModel*) getModel())->setRateMatrix(rate_param_bk);
//			}
//			dynamic_cast<RateGamma*>(getRate())->setGammaShape(alpha_bk);
//			getModel()->decomposeRateMatrix();
//			cout << "Reset rate parameters!" << endl;
		} else {
			curScore = modOptScore;
			newTree = getTreeString();
		}
	}
	return newTree;
}

double IQTree::doTreeSearch() {
//    double begin_real_time, cur_real_time;
//    begin_real_time = getRealTime();
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    //printResultTree(params);
    // PLEASE PRINT TREE HERE!
    printResultTree();
    string treels_name = params->out_prefix;
    treels_name += ".treels";
    string out_lh_file = params->out_prefix;
    out_lh_file += ".treelh";
    string site_lh_file = params->out_prefix;
    site_lh_file += ".sitelh";

    if (params->print_tree_lh) {
        out_treelh.open(out_lh_file.c_str());
        out_sitelh.open(site_lh_file.c_str());
    }

    if (params->write_intermediate_trees)
        out_treels.open(treels_name.c_str());

    if (params->write_intermediate_trees && save_all_trees != 2) {
        printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
    }

    setRootNode(params->root);
    // keep the best tree into a string
    stringstream bestTreeStream;
    stringstream bestTopoStream;
    string perturb_tree_string;
    string imd_tree;
    printTree(bestTreeStream, WT_TAXON_ID + WT_BR_LEN);
    printTree(bestTopoStream, WT_TAXON_ID + WT_SORT_TAXA);
    string best_tree_topo = bestTopoStream.str();

    stop_rule.addImprovedIteration(1);
    searchinfo.curPerStrength = params->initPerStrength;

	double cur_correlation = 0.0;

	/*====================================================
	 * MAIN LOOP OF THE IQ-TREE ALGORITHM
	 *====================================================*/
    for (curIt = 2; !stop_rule.meetStopCondition(curIt, cur_correlation); curIt++) {
        searchinfo.curIter = curIt;
        // estimate logl_cutoff for bootstrap
        if (params->avoid_duplicated_trees && max_candidate_trees > 0 && treels_logl.size() > 1000) {
        	int predicted_iteration = ((curIt+params->step_iterations-1)/params->step_iterations)*params->step_iterations;
            int num_entries = floor(max_candidate_trees * ((double) curIt / predicted_iteration));
            if (num_entries < treels_logl.size() * 0.9) {
                DoubleVector logl = treels_logl;
                nth_element(logl.begin(), logl.begin() + (treels_logl.size() - num_entries), logl.end());
                logl_cutoff = logl[treels_logl.size() - num_entries] - 1.0;
            } else
                logl_cutoff = 0.0;
            if (verbose_mode >= VB_MED) {
                if (curIt % 10 == 0) {
                    cout << treels.size() << " trees, " << treels_logl.size() << " logls, logl_cutoff= " << logl_cutoff;
                    if (params->store_candidate_trees)
                        cout << " duplicates= " << duplication_counter << " ("
                                << (int) round(100 * ((double) duplication_counter / treels_logl.size())) << "%)" << endl;
                    else
                        cout << endl;
                }
            }
        }

        if (estimate_nni_cutoff && nni_info.size() >= 500) {
            estimate_nni_cutoff = false;
            estimateNNICutoff(params);
        }

        Alignment *saved_aln = aln;

    	/*----------------------------------------
    	 * Perturb the tree
    	 *---------------------------------------*/
        double perturbScore;
        if (iqp_assess_quartet == IQP_BOOTSTRAP) {
            // create bootstrap sample
            Alignment* bootstrap_alignment;
            if (aln->isSuperAlignment())
                bootstrap_alignment = new SuperAlignment;
            else
                bootstrap_alignment = new Alignment;
            bootstrap_alignment->createBootstrapAlignment(aln, NULL, params->bootstrap_spec);
            setAlignment(bootstrap_alignment);
            initializeAllPartialLh();
            clearAllPartialLH();
            curScore = optimizeAllBranches();
        } else {
            if (params->snni) {
                int numNNI = floor(searchinfo.curPerStrength * (aln->getNSeq() - 3));
                //cout << "candidateTrees.size() = " << candidateTrees.size() << endl;
                string candidateTree = candidateTrees.getRandCandTree();
                readTreeString(candidateTree);
                if (params->iqp) {
                    doIQP();
                } else {
                    doRandomNNIs(numNNI);
                }
            } else {
                doIQP();
            }
            setAlignment(aln);
            setRootNode(params->root);
            perturb_tree_string = getTreeString();
            if (params->count_trees) {
                string perturb_tree_topo = getTopology();
                if (pllTreeCounter.find(perturb_tree_topo) == pllTreeCounter.end()) {
                    // not found in hash_map
                    pllTreeCounter[perturb_tree_topo] = 1;
                } else {
                    // found in hash_map
                    pllTreeCounter[perturb_tree_topo]++;
                }
            }

            if (params->pll) {
                pllNewickTree *perturbTree = pllNewickParseString(perturb_tree_string.c_str());
                assert(perturbTree != NULL);
                pllTreeInitTopologyNewick(pllInst, perturbTree, PLL_FALSE);
                pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_TRUE, PLL_FALSE);
                if (params->numSmoothTree >= 1) {
                    pllOptimizeBranchLengths(pllInst, pllPartitions, params->numSmoothTree);
                }
                pllNewickParseDestroy(&perturbTree);
                curScore = pllInst->likelihood;
                perturbScore = curScore;
            } else {
                initializeAllPartialLh();
                clearAllPartialLH();
                if (isSuperTree()) {
                    ((PhyloSuperTree*) this)->mapTrees();
                }
                curScore = optimizeAllBranches(1, TOL_LIKELIHOOD, PLL_NEWZPERCYCLE);
                perturbScore = curScore;
            }
        }

    	/*----------------------------------------
    	 * Optimize tree with NNI
    	 *---------------------------------------*/
        int nni_count = 0;
        int nni_steps = 0;
        if (params->pll) {
        	if (params->partition_file)
        		outError("Unsupported -pll -sp combination!");
            curScore = pllOptimizeNNI(nni_count, nni_steps, searchinfo);
            pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_TRUE,
                    PLL_TRUE, 0, 0, 0, PLL_SUMMARIZE_LH, 0, 0);
            imd_tree = string(pllInst->tree_string);
            readTreeString(imd_tree);
        } else {
            curScore = optimizeNNI(nni_count, nni_steps);
            //imd_tree = getTreeString();
            if (isSuperTree()) {
                ((PhyloSuperTree*) this)->computeBranchLengths();
            }
            imd_tree = getTreeString();
        }

        if (iqp_assess_quartet == IQP_BOOTSTRAP) {
            // restore alignment
            delete aln;
            setAlignment(saved_aln);
            initializeAllPartialLh();
            clearAllPartialLH();
        }

        if (isSuperTree()) {
            ((PhyloSuperTree*) this)->computeBranchLengths();
        }

    	/*----------------------------------------
    	 * Print information
    	 *---------------------------------------*/
        double realtime_remaining = stop_rule.getRemainingTime(curIt);
        cout.setf(ios::fixed, ios::floatfield);

        cout << ((iqp_assess_quartet == IQP_BOOTSTRAP) ? "Bootstrap " : "Iteration ") << curIt << " / LogL: ";
        if (verbose_mode >= VB_MED)
        	cout << perturbScore << " -> ";
        cout << curScore;
        if (verbose_mode >= VB_MED)
        	cout << " / NNIs: " << nni_count << "," << nni_steps;
        cout << " / Time: " << convert_time(getRealTime() - params->start_real_time);

        if (curIt > 10) {
			cout << " (" << convert_time(realtime_remaining) << " left)";
        }
        cout << endl;

        if (params->write_intermediate_trees && save_all_trees != 2) {
            printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
        }

    	/*----------------------------------------
    	 * Update if better tree is found
    	 *---------------------------------------*/
        if (curScore > bestScore) {
            stringstream cur_tree_topo_ss;
            setRootNode(params->root);
            printTree(cur_tree_topo_ss, WT_TAXON_ID | WT_SORT_TAXA);
            if (cur_tree_topo_ss.str() != best_tree_topo) {
                best_tree_topo = cur_tree_topo_ss.str();
                imd_tree = optimizeModelParameters();
                stop_rule.addImprovedIteration(curIt);
                cout << "BETTER TREE FOUND at iteration " << curIt << ": " << curScore;
                cout << " / CPU time: " << (int) round(getCPUTime() - params->startCPUTime) << "s" << endl << endl;
                if (curScore > bestScore) {
                    searchinfo.curPerStrength = params->initPerStrength;
                }
            } else {
                cout << "UPDATE BEST LOG-LIKELIHOOD: " << curScore << endl;
            }
            setBestTree(imd_tree, curScore);
            if (params->write_best_trees) {
                ostringstream iter_string;
                iter_string << curIt;
                printResultTree(iter_string.str());
            }
            printResultTree();
        }

        // check whether the tree can be put into the reference set
        if (params->snni) {
        	candidateTrees.update(imd_tree, curScore);
        	if (verbose_mode >= VB_MED) {
        		candidateTrees.printBestScores();
        		cout << endl;
        	}
        } else {
            // The IQPNNI algorithm
            readTreeString(bestTreeString);
        }

        // DTH: make pllUFBootData usable in summarizeBootstrap
        if(params->pll && params->online_bootstrap && (params->gbo_replicates > 0))
            pllConvertUFBootData2IQTree();
        // DTH: Carefully watch the -pll case here


    	/*----------------------------------------
    	 * convergence criterion for ultrafast bootstrap
    	 *---------------------------------------*/
        if ((curIt) % (params->step_iterations / 2) == 0 && params->stop_condition == SC_BOOTSTRAP_CORRELATION) {
        	// compute split support every half step
            SplitGraph *sg = new SplitGraph;
            summarizeBootstrap(*sg);
            boot_splits.push_back(sg);
            if (params->max_candidate_trees == 0)
                max_candidate_trees = treels_logl.size() * (curIt + (params->step_iterations / 2)) / curIt;
			cout << "INFO: " << treels_logl.size() << " bootstrap candidate trees evaluated (logl-cutoff: " << logl_cutoff << ")" << endl;

			// check convergence every full step
			if (curIt % params->step_iterations == 0) {
	        	cur_correlation = computeBootstrapCorrelation();
	            cout << "INFO: Bootstrap correlation coefficient of split occurrence frequencies: " << cur_correlation << endl;
	            if (!stop_rule.meetStopCondition(curIt, cur_correlation)) {
	                if (params->max_candidate_trees == 0) {
	                    max_candidate_trees = treels_logl.size() * (curIt + params->step_iterations) / curIt;
	                }
	                cout << "INFO: UFBoot does not converge, continue " << params->step_iterations << " more iterations" << endl;
	            }
	        }
        } // end of bootstrap convergence test
    }

    readTreeString(bestTreeString);

    // DTH: Carefully watch the -pll case here

    // OBSOLETE, consider to remove
//    if (params->stop_condition != SC_UNSUCCESS_ITERATION) {
//    	if (params->snni) {
//			int predicted_iteration = stop_rule.getLastImprovedIteration() + params->unsuccess_iteration;
//			if (predicted_iteration > curIt) {
//				cout << endl << "WARNING: " << predicted_iteration <<
//					" iterations are needed by the default stopping rule" << endl;
//			}
//    	} else {
//			int predicted_iteration = stop_rule.getPredictedIteration();
//			if (predicted_iteration > curIt) {
//				cout << endl << "WARNING: " << predicted_iteration << " iterations are needed to ensure "
//						<< int(floor(params->stop_confidence * 100)) << "% confidence that "
//						<< "         the IQPNNI search will not find a better tree" << endl;
//			}
//    	}
//    }

    if (testNNI)
        outNNI.close();
    if (params->write_intermediate_trees)
        out_treels.close();
    if (params->print_tree_lh) {
        out_treelh.close();
        out_sitelh.close();
    }

    // DTH: pllUFBoot deallocation
    if(params->pll) {
        pllDestroyUFBootData();
    }

    return bestScore;
}

/****************************************************************************
 Fast Nearest Neighbor Interchange by maximum likelihood
 ****************************************************************************/
double IQTree::optimizeNNI(int &nni_count, int &nni_steps) {
    bool rollBack = false;
    nni_count = 0;
    int numNNIs = 0; // number of NNI to be applied in each step
    int MAXSTEPS = 50; // maximum number of NNI steps
    for (nni_steps = 1; nni_steps <= MAXSTEPS; nni_steps++) {
        double oldScore = curScore;
        if (!rollBack) { // tree get improved and was not rollbacked
            if (save_all_trees == 2) {
                saveCurrentTree(curScore); // BQM: for new bootstrap
            }
            if (verbose_mode >= VB_DEBUG) {
                cout << "Doing NNI round " << nni_steps << endl;
                if (isSuperTree()) {
                    ((PhyloSuperTree*) this)->printMapInfo();
                }
            }

            nonConfNNIs.clear(); // Vector containing non-conflicting positive NNIs
            optBrans.clear(); // Vector containing branch length of the positive NNIs
            orgBrans.clear(); // Vector containing all current branch of the tree
            plusNNIs.clear(); // Vector containing all positive NNIs
            saveBranches(); // save all current branch lengths
            initPartitionInfo(); // for super tree
            if (searchinfo.speednni && !brans2Eval.empty())
              evalNNIs(brans2Eval);
            else
              evalNNIs();

//            if (!nni_sort) {
//                evalNNIs(); // generate all positive NNI moves
//            } else {
//                evalNNIsSort(params->approximate_nni);
//            }

            /* sort all positive NNI moves (descending) */
            sort(plusNNIs.begin(), plusNNIs.end());
            if (verbose_mode >= VB_DEBUG) {
                cout << "curScore: " << curScore << endl;
                for (int i = 0; i < plusNNIs.size(); i++) {
                    cout << "Logl of positive NNI " << i << " : " << plusNNIs[i].newloglh << endl;
                }
            }

            if (plusNNIs.size() == 0) {
                break;
            }

            /* remove conflicting NNIs */
            genNonconfNNIs();
            numNNIs = nonConfNNIs.size();
            if (verbose_mode >= VB_DEBUG) {
                for (int i = 0; i < nonConfNNIs.size(); i++) {
                    cout << "Log-likelihood of non-conflicting NNI " << i << " : " << nonConfNNIs[i].newloglh << endl;
                }
            }
        }
        // Apply all non-conflicting positive NNIs
        doNNIs(numNNIs);

        if (searchinfo.speednni) {
            brans2Eval.clear();
            updateBrans2Eval(appliedNNIs);
            appliedNNIs.clear();
        }

        // Re-estimate branch lengths of the new tree
        curScore = optimizeAllBranches(1, TOL_LIKELIHOOD, PLL_NEWZPERCYCLE);

		// curScore should be larger than score of the best NNI
        if (curScore >= nonConfNNIs.at(0).newloglh ) {
            nni_count += numNNIs;
            rollBack = false;
        } else {
            /* tree cannot be worse if only 1 NNI is applied */
//            if (numNNIs == 1) {
//                cout << "ERROR / POSSIBLE BUG: current logl=" << curScore << " < " << nonConfNNIs.at(0).newloglh
//                        << "(best NNI)" << endl;
//                abort();
//            }
            if (verbose_mode >= VB_MED) {
                cout << "New score = " << curScore << " after applying " << numNNIs <<
                        " is worse than score = " << nonConfNNIs.at(0).newloglh
                        << " of the best NNI. Roll back tree ..." << endl;
            }

            // restore the tree by reverting all NNIs
            for (int i = 0; i < numNNIs; i++)
                doNNI(nonConfNNIs.at(i));
            // restore the branch lengths
            restoreAllBrans();
            // This is important because after restoring the branch lengths, all partial
            // likelihood need to be cleared.
            clearAllPartialLH();
            rollBack = true;
            // only apply the best NNI
            numNNIs = 1;
            curScore = oldScore;
        }
    }

    if (nni_count == 0) {
        cout << "INFO: Tree is readily NNI-optimized" << endl;
    }
    brans2Eval.clear();
    return curScore;
}

void IQTree::updateBrans2Eval(vector<NNIMove> nnis) {
    for (vector<NNIMove>::iterator it = nnis.begin(); it != nnis.end(); it++) {
        Branch bran((*it).node1, (*it).node2);
        brans2Eval.insert(pair<string, Branch>(bran.getKey(), bran));
        getInBranches(brans2Eval, 2, (*it).node1, (*it).node2);
        getInBranches(brans2Eval, 2, (*it).node2, (*it).node1);
    }
}


double IQTree::pllOptimizeNNI(int &totalNNICount, int &nniSteps, SearchInfo &searchinfo) {
    if((globalParam->online_bootstrap == PLL_TRUE) && (globalParam->gbo_replicates > 0)) {
        pllInitUFBootData();
    }
    searchinfo.numAppliedNNIs = 0;
    searchinfo.curLogl = curScore;
    //cout << "curLogl: " << searchinfo.curLogl << endl;
    const int MAX_NNI_STEPS = 50;
    totalNNICount = 0;
    for (nniSteps = 1; nniSteps <= MAX_NNI_STEPS; nniSteps++) {
        searchinfo.curNumNNISteps = nniSteps;
        searchinfo.posNNIList.clear();
        double newLH = pllDoNNISearch(pllInst, pllPartitions, searchinfo);
        if (searchinfo.curNumAppliedNNIs == 0) { // no positive NNI was found
            searchinfo.curLogl = newLH;
            break;
        } else {
            searchinfo.curLogl = newLH;
            searchinfo.numAppliedNNIs += searchinfo.curNumAppliedNNIs;
        }
    }

    if (nniSteps == (MAX_NNI_STEPS + 1)) {
        cout << "WARNING: NNI search seems to run unusually too long and thus it was stopped!" << endl;
    }

    totalNNICount = searchinfo.numAppliedNNIs;
    pllInst->likelihood = searchinfo.curLogl;
    return searchinfo.curLogl;
}

void IQTree::pllLogBootSamples(int** pll_boot_samples, int nsamples, int npatterns){
    ofstream bfile("boot_samples.log");
    bfile << "Original freq:" << endl;
    int sum = 0;
    for(int i = 0; i < pllAlignment->sequenceLength; i++){
        bfile << setw(4) << pllInst->aliaswgt[i];
        sum += pllInst->aliaswgt[i];
    }
    bfile << endl << "sum = " << sum << endl;

    bfile << "Bootstrap freq:" << endl;

    for(int i = 0; i < nsamples; i++){
        sum = 0;
        for(int j = 0; j < npatterns; j++){
            bfile << setw(4) << pll_boot_samples[i][j];
            sum += pll_boot_samples[i][j];
        }
        bfile << endl << "sum = "  << sum << endl;
    }
    bfile.close();

}

void IQTree::pllInitUFBootData(){
    if(pllUFBootDataPtr == NULL){
        pllUFBootDataPtr = (pllUFBootData *) malloc(sizeof(pllUFBootData));
        pllUFBootDataPtr->boot_samples = NULL;
        pllUFBootDataPtr->candidate_trees_count = 0;

        if(params->online_bootstrap && params->gbo_replicates > 0){
        	if(!pll2iqtree_pattern_index) pllBuildIQTreePatternIndex();

            pllUFBootDataPtr->treels = pllHashInit(max_candidate_trees);
            pllUFBootDataPtr->treels_size = max_candidate_trees; // track size of treels_logl, treels_newick, treels_ptnlh

            pllUFBootDataPtr->treels_logl =
                (double *) malloc(max_candidate_trees * (sizeof(double)));
            if(!pllUFBootDataPtr->treels_logl) outError("Not enough dynamic memory!");
            //memset(pllUFBootDataPtr->treels_logl, 0, max_candidate_trees * (sizeof(double)));

            pllUFBootDataPtr->treels_newick =
                (char **) malloc(max_candidate_trees * (sizeof(char *)));
            if(!pllUFBootDataPtr->treels_newick) outError("Not enough dynamic memory!");
            memset(pllUFBootDataPtr->treels_newick, 0, max_candidate_trees * (sizeof(char *)));


            pllUFBootDataPtr->treels_ptnlh =
                (double **) malloc(max_candidate_trees * (sizeof(double *)));
            if(!pllUFBootDataPtr->treels_ptnlh) outError("Not enough dynamic memory!");
            memset(pllUFBootDataPtr->treels_ptnlh, 0, max_candidate_trees * (sizeof(double *)));

            // aln->createBootstrapAlignment() must be called before this fragment
            pllUFBootDataPtr->boot_samples =
                (int **) malloc(params->gbo_replicates * sizeof(int *));
            if(!pllUFBootDataPtr->boot_samples) outError("Not enough dynamic memory!");
            for(int i = 0; i < params->gbo_replicates; i++){
                pllUFBootDataPtr->boot_samples[i] =
                    (int *) malloc(pllAlignment->sequenceLength * sizeof(int));
                if(!pllUFBootDataPtr->boot_samples[i]) outError("Not enough dynamic memory!");
                for(int j = 0; j < pllAlignment->sequenceLength; j++){
                    pllUFBootDataPtr->boot_samples[i][j] =
                        boot_samples[i][pll2iqtree_pattern_index[j]];
                }
            }

//            pllLogBootSamples(pllUFBootDataPtr->boot_samples,
//                    params->gbo_replicates, pllAlignment->sequenceLength);

            pllUFBootDataPtr->boot_logl =
                (double *) malloc(params->gbo_replicates * (sizeof(double)));
            if(!pllUFBootDataPtr->boot_logl) outError("Not enough dynamic memory!");
            for(int i = 0; i < params->gbo_replicates; i++)
                pllUFBootDataPtr->boot_logl[i] = -DBL_MAX;

            pllUFBootDataPtr->boot_counts =
                (int *) malloc(params->gbo_replicates * (sizeof(int)));
            if(!pllUFBootDataPtr->boot_counts) outError("Not enough dynamic memory!");
            memset(pllUFBootDataPtr->boot_counts, 0, params->gbo_replicates * (sizeof(int)));

            pllUFBootDataPtr->boot_trees =
                (int *) malloc(params->gbo_replicates * (sizeof(int)));
            if(!pllUFBootDataPtr->boot_trees) outError("Not enough dynamic memory!");

            pllUFBootDataPtr->duplication_counter = 0;
        }
    }
    pllUFBootDataPtr->max_candidate_trees = max_candidate_trees;
    pllUFBootDataPtr->save_all_trees = save_all_trees;
    pllUFBootDataPtr->save_all_br_lens = save_all_br_lens;
    pllUFBootDataPtr->logl_cutoff = logl_cutoff;
    pllUFBootDataPtr->n_patterns = pllAlignment->sequenceLength;
}

void IQTree::pllDestroyUFBootData(){
    if(pll2iqtree_pattern_index){
        delete [] pll2iqtree_pattern_index;
        pll2iqtree_pattern_index = NULL;
    }

    if(params->online_bootstrap && params->gbo_replicates > 0){
        pllHashDestroy(&(pllUFBootDataPtr->treels), PLL_TRUE);

        free(pllUFBootDataPtr->treels_logl);

        for(int i = 0; i < pllUFBootDataPtr->candidate_trees_count; i++)
            if(pllUFBootDataPtr->treels_newick[i])
                free(pllUFBootDataPtr->treels_newick[i]);
        free(pllUFBootDataPtr->treels_newick);

        for(int i = 0; i < pllUFBootDataPtr->treels_size; i++)
            if(pllUFBootDataPtr->treels_ptnlh[i])
                free(pllUFBootDataPtr->treels_ptnlh[i]);
        free(pllUFBootDataPtr->treels_ptnlh);

        for(int i = 0; i < params->gbo_replicates; i++)
            free(pllUFBootDataPtr->boot_samples[i]);
        free(pllUFBootDataPtr->boot_samples);

        free(pllUFBootDataPtr->boot_logl);

        free(pllUFBootDataPtr->boot_counts);

        free(pllUFBootDataPtr->boot_trees);
    }
    free(pllUFBootDataPtr);
    pllUFBootDataPtr = NULL;
}


void IQTree::doNNIs(int nni2apply, bool changeBran) {
    for (int i = 0; i < nni2apply; i++) {
        doNNI(nonConfNNIs.at(i));
        appliedNNIs.push_back(nonConfNNIs.at(i));
        if (!params->leastSquareNNI && changeBran) {
            // apply new branch lengths
            changeNNIBrans(nonConfNNIs.at(i));
        }
    }
}


void IQTree::genNonconfNNIs() {
    for (vector<NNIMove>::iterator iterMove = plusNNIs.begin(); iterMove != plusNNIs.end(); iterMove++) {
        bool choosen = true;
        for (vector<NNIMove>::iterator iterNextMove = nonConfNNIs.begin(); iterNextMove != nonConfNNIs.end();
                iterNextMove++) {
            if ((*iterMove).node1 == (*(iterNextMove)).node1 || (*iterMove).node2 == (*(iterNextMove)).node1
                    || (*iterMove).node1 == (*(iterNextMove)).node2 || (*iterMove).node2 == (*(iterNextMove)).node2) {
                choosen = false;
                break;
            }
        }
        if (choosen) {
            nonConfNNIs.push_back(*iterMove);
        }
    }
}

//double IQTree::estN95() {
//    if (vecNumNNI.size() == 0) {
//        return 0;
//    } else {
//        sort(vecNumNNI.begin(), vecNumNNI.end());
//        int index = floor(vecNumNNI.size() * speed_conf);
//        return vecNumNNI[index];
//    }
//}

double IQTree::getAvgNumNNI() {
    if (vecNumNNI.size() == 0) {
        return 0;
    } else {
        double median;
        size_t size = vecNumNNI.size();
        sort(vecNumNNI.begin(), vecNumNNI.end());
        if (size % 2 == 0) {
            median = (vecNumNNI[size / 2 + 1] + vecNumNNI[size / 2]) / 2;
        } else {
            median = vecNumNNI[size / 2];
        }
        return median;
    }
}

double IQTree::estDeltaMedian() {
    if (vecImpProNNI.size() == 0) {
        return 0;
    } else {
        double median;
        size_t size = vecImpProNNI.size();
        sort(vecImpProNNI.begin(), vecImpProNNI.end());
        if (size % 2 == 0) {
            median = (vecImpProNNI[size / 2 + 1] + vecImpProNNI[size / 2]) / 2;
        } else {
            median = vecImpProNNI[size / 2];
        }
        return median;
    }
}

//inline double IQTree::estDelta95() {
//    if (vecImpProNNI.size() == 0) {
//        return 0;
//    } else {
//        sort(vecImpProNNI.begin(), vecImpProNNI.end());
//        int index = floor(vecImpProNNI.size() * speed_conf);
//        return vecImpProNNI[index];
//    }
//}

int IQTree::getDelete() const {
    return k_delete;
}

void IQTree::setDelete(int _delete) {
    k_delete = _delete;
}

void IQTree::changeBranLen(PhyloNode *node1, PhyloNode *node2, double newlen) {
    node1->findNeighbor(node2)->length = newlen;
    node2->findNeighbor(node1)->length = newlen;
    node1->clearReversePartialLh(node2);
    node2->clearReversePartialLh(node1);
}

double IQTree::getBranLen(PhyloNode *node1, PhyloNode *node2) {
    return  node1->findNeighbor(node2)->length;
}

void IQTree::saveBranches(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    if (dad) {
        double len = getBranLen(node, dad);
        string key = nodePair2String(node, dad);
        orgBrans.insert(mapString2Double::value_type(key, len));
    }

    FOR_NEIGHBOR_IT(node, dad, it){
    saveBranches((PhyloNode*) (*it)->node, node);
}
}

void IQTree::restoreAllBrans(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    if (dad) {
        string key = nodePair2String(node, dad);
        Neighbor* bran_it = node->findNeighbor(dad);
        assert(bran_it);
        Neighbor* bran_it_back = dad->findNeighbor(node);
        assert(bran_it_back);
        assert(orgBrans.count(key));
        bran_it->length = orgBrans[key];
        bran_it_back->length = orgBrans[key];
    }

    FOR_NEIGHBOR_IT(node, dad, it){
    restoreAllBrans((PhyloNode*) (*it)->node, node);
}
}


void IQTree::evalNNIs(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    // internal branch
    if (!node->isLeaf() && dad && !dad->isLeaf()) {
        NNIMove myMove = getBestNNIForBran(node, dad, NULL);
        if (myMove.newloglh > curScore + params->loglh_epsilon) {
            addPositiveNNIMove(myMove);
        }
    }

    FOR_NEIGHBOR_IT(node, dad, it){
        evalNNIs((PhyloNode*) (*it)->node, node);
    }
}

void IQTree::evalNNIs(map<string, Branch> brans) {
    for (map<string, Branch>::iterator it = brans.begin(); it != brans.end(); it++) {
        NNIMove myMove = getBestNNIForBran((PhyloNode*)it->second.node1, (PhyloNode*) it->second.node2, NULL);
        if (myMove.newloglh > curScore + params->loglh_epsilon) {
            addPositiveNNIMove(myMove);
        }
    }
}

/**
 *  Currently not used, commented out to simplify the interface of getBestNNIForBran
void IQTree::evalNNIsSort(bool approx_nni) {
        if (myMove.newloglh > curScore + params->loglh_epsilon) {
        if (myMove.newloglh > curScore + params->loglh_epsilon) {
            addPositiveNNIMove(myMove);
        }
            addPositiveNNIMove(myMove);
        }
    NodeVector nodes1, nodes2;
    int i;
    double cur_lh = curScore;
    vector<IntBranchInfo> int_branches;

    getInternalBranches(nodes1, nodes2);
    assert(nodes1.size() == leafNum - 3 && nodes2.size() == leafNum - 3);

    for (i = 0; i < leafNum - 3; i++) {
        IntBranchInfo int_branch;
        PhyloNeighbor *node12_it = (PhyloNeighbor*) nodes1[i]->findNeighbor(nodes2[i]);
        //PhyloNeighbor *node21_it = (PhyloNeighbor*) nodes2[i]->findNeighbor(nodes1[i]);
        int_branch.lh_contribution = cur_lh - computeLikelihoodZeroBranch(node12_it, (PhyloNode*) nodes1[i]);
        if (int_branch.lh_contribution < 0.0)
            int_branch.lh_contribution = 0.0;
        if (int_branch.lh_contribution < fabs(nni_cutoff)) {
            int_branch.node1 = (PhyloNode*) nodes1[i];
            int_branch.node2 = (PhyloNode*) nodes2[i];
            int_branches.push_back(int_branch);
        }
    }
    std::sort(int_branches.begin(), int_branches.end(), int_branch_cmp);
    for (vector<IntBranchInfo>::iterator it = int_branches.begin(); it != int_branches.end(); it++)
        if (it->lh_contribution >= 0.0) // evaluate NNI if branch contribution is big enough
                {
            NNIMove myMove = getBestNNIForBran(it->node1, it->node2, NULL, approx_nni, it->lh_contribution);
            if (myMove.newloglh > curScore) {
                addPositiveNNIMove(myMove);
                if (!estimate_nni_cutoff)
                    for (vector<IntBranchInfo>::iterator it2 = it + 1; it2 != int_branches.end(); it2++) {
                        if (it2->node1 == it->node1 || it2->node2 == it->node1 || it2->node1 == it->node2
                                || it2->node2 == it->node2)
                            it2->lh_contribution = -1.0; // do not evaluate this branch later on
                    }
            }
        } else { // otherwise, only optimize the branch length
            PhyloNode *node1 = it->node1;
            PhyloNode *node2 = it->node2;
            PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
            PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);
            double stored_len = node12_it->length;
            curScore = optimizeOneBranch(node1, node2, false);
            string key("");
            if (node1->id < node2->id) {
                key += convertIntToString(node1->id) + "->" + convertIntToString(node2->id);
            } else {
                key += convertIntToString(node2->id) + "->" + convertIntToString(node1->id);
            }

            optBrans.insert(mapString2Double::value_type(key, node12_it->length));
            node12_it->length = stored_len;
            node21_it->length = stored_len;
        }
}
*/

void IQTree::estimateNNICutoff(Params* params) {
    double *delta = new double[nni_info.size()];
    int i;
    for (i = 0; i < nni_info.size(); i++) {
        double lh_score[4];
        memmove(lh_score, nni_info[i].lh_score, 4 * sizeof(double));
        std::sort(lh_score + 1, lh_score + 4); // sort in ascending order
        delta[i] = lh_score[0] - lh_score[2];
        if (verbose_mode >= VB_MED)
            cout << i << ": " << lh_score[0] << " " << lh_score[1] << " " << lh_score[2] << " " << lh_score[3] << endl;
    }
    std::sort(delta, delta + nni_info.size());
    nni_cutoff = delta[nni_info.size() / 20];
    cout << endl << "Estimated NNI cutoff: " << nni_cutoff << endl;
    string file_name = params->out_prefix;
    file_name += ".nnidelta";
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(file_name.c_str());
        for (i = 0; i < nni_info.size(); i++) {
            out << delta[i] << endl;
        }
        out.close();
        cout << "NNI delta printed to " << file_name << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, file_name);
    }
    delete[] delta;
}

void IQTree::saveCurrentTree(double cur_logl) {
    ostringstream ostr;
    string tree_str;
    StringIntMap::iterator it = treels.end();
    if (params->store_candidate_trees) {
        printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
        tree_str = ostr.str();
        it = treels.find(tree_str);
    }
    int tree_index = -1;
    if (it != treels.end()) { // already in treels
        duplication_counter++;
        tree_index = it->second;
        if (cur_logl <= treels_logl[it->second] + 1e-4) {
            if (cur_logl < treels_logl[it->second] - 5.0)
                if (verbose_mode >= VB_MED)
                    cout << "Current lh " << cur_logl << " is much worse than expected " << treels_logl[it->second]
                            << endl;
            return;
        }
        if (verbose_mode >= VB_MAX)
            cout << "Updated logl " << treels_logl[it->second] << " to " << cur_logl << endl;
        treels_logl[it->second] = cur_logl;
        if (save_all_br_lens) {
            ostr.seekp(ios::beg);
            printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
            treels_newick[it->second] = ostr.str();
        }
        if (boot_samples.empty()) {
            computePatternLikelihood(treels_ptnlh[it->second], &cur_logl);
            return;
        }
        if (verbose_mode >= VB_MAX)
            cout << "Update treels_logl[" << tree_index << "] := " << cur_logl << endl;
    } else {
        if (logl_cutoff != 0.0 && cur_logl <= logl_cutoff + 1e-4)
            return;
        tree_index = treels_logl.size();
        if (params->store_candidate_trees)
            treels[tree_str] = tree_index;
        treels_logl.push_back(cur_logl);
        if (verbose_mode >= VB_MAX)
            cout << "Add    treels_logl[" << tree_index << "] := " << cur_logl << endl;
    }

    if (write_intermediate_trees)
        printTree(out_treels, WT_NEWLINE | WT_BR_LEN);

    int nptn = getAlnNPattern();
    BootValType *pattern_lh = aligned_alloc<BootValType>(nptn);

#ifdef BOOT_VAL_FLOAT
    double *pattern_lh_orig = aligned_alloc_double(nptn);
    computePatternLikelihood(pattern_lh_orig, &cur_logl);
    for (int i = 0; i < nptn; i++)
    	pattern_lh[i] = pattern_lh_orig[i];
#else
    computePatternLikelihood(pattern_lh, &cur_logl);
#endif


    if (boot_samples.empty()) {
        // for runGuidedBootstrap
#ifdef BOOT_VAL_FLOAT
        treels_ptnlh.push_back(pattern_lh_orig);
#else
        treels_ptnlh.push_back(pattern_lh);
#endif
    } else {
        // online bootstrap
        int ptn, nptn = getAlnNPattern();
        int updated = 0;
        int nsamples = boot_samples.size();

        for (int sample = 0; sample < nsamples; sample++) {
            double rell = 0.0;

            // TODO: The following parallel is not very efficient, should wrap the above loop
//#ifdef _OPENMP
//#pragma omp parallel for reduction(+: rell)
//#endif
//            if (sse == LK_NORMAL || sse == LK_EIGEN) {
            if (false) {
            	BootValType *boot_sample = boot_samples[sample];
				for (ptn = 0; ptn < nptn; ptn++)
					rell += pattern_lh[ptn] * boot_sample[ptn];
            } else {
            	// SSE optimized version of the above loop
				BootValType *boot_sample = boot_samples[sample];
#ifdef BOOT_VAL_FLOAT
				VectorClassFloat vc_rell = 0.0;
				int maxptn = nptn - VCSIZE_FLOAT;
				for (ptn = 0; ptn < maxptn; ptn+=VCSIZE_FLOAT)
					vc_rell = mul_add(VectorClassFloat().load_a(&pattern_lh[ptn]), VectorClassFloat().load_a(&boot_sample[ptn]), vc_rell);
#else
				VectorClassMaster vc_rell = 0.0;
				int maxptn = nptn - VCSIZE_MASTER;
				for (ptn = 0; ptn < maxptn; ptn+=VCSIZE_MASTER)
					vc_rell = mul_add(VectorClassMaster().load_a(&pattern_lh[ptn]), VectorClassMaster().load_a(&boot_sample[ptn]), vc_rell);
#endif

				BootValType res = horizontal_add(vc_rell);
				// add the remaining ptn
				for (; ptn < nptn; ptn++)
					res += pattern_lh[ptn] * boot_sample[ptn];
				rell = res;
            }

            if (rell > boot_logl[sample] + params->ufboot_epsilon
                    || (rell > boot_logl[sample] - params->ufboot_epsilon
                            && random_double() <= 1.0 / (boot_counts[sample] + 1))) {
                if (tree_str == "") {
                    printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
                    tree_str = ostr.str();
                    it = treels.find(tree_str);
                    if (it != treels.end()) {
                        tree_index = it->second;
                    } else {
                        tree_index = treels.size();
                        treels[tree_str] = tree_index;
                    }
                }
                if (rell <= boot_logl[sample] + params->ufboot_epsilon) {
                    boot_counts[sample]++;
                } else {
                    boot_counts[sample] = 1;
                }
                boot_logl[sample] = max(boot_logl[sample], rell);
                boot_trees[sample] = tree_index;
                updated++;
            } /*else if (verbose_mode >= VB_MED && rell > boot_logl[sample] - 0.01) {
             cout << "Info: multiple RELL score trees detected" << endl;
             }*/
        }
        if (updated && verbose_mode >= VB_MAX)
            cout << updated << " boot trees updated" << endl;
        /*
         if (tree_index >= max_candidate_trees/2 && boot_splits->empty()) {
         // summarize split support half way for stopping criterion
         cout << "Summarizing current bootstrap supports..." << endl;
         summarizeBootstrap(*boot_splits);
         }*/
    }
    if (save_all_br_lens) {
        ostr.seekp(ios::beg);
        printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
        treels_newick.push_back(ostr.str());
    }
    if (print_tree_lh) {
        out_treelh << cur_logl;
        double prob;
#ifdef BOOT_VAL_FLOAT
        aln->multinomialProb(pattern_lh_orig, prob);
#else
        aln->multinomialProb(pattern_lh, prob);
#endif
        out_treelh << "\t" << prob << endl;

        IntVector pattern_index;
        aln->getSitePatternIndex(pattern_index);
        out_sitelh << "Site_Lh   ";
        for (int i = 0; i < getAlnNSite(); i++)
            out_sitelh << " " << pattern_lh[pattern_index[i]];
        out_sitelh << endl;
    }

    if (!boot_samples.empty()) {
#ifdef BOOT_VAL_FLOAT
    	aligned_free(pattern_lh_orig);
#endif
    	aligned_free(pattern_lh);
    } else {
#ifdef BOOT_VAL_FLOAT
    	aligned_free(pattern_lh);
#endif
    }

}

void IQTree::saveNNITrees(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    if (dad && !node->isLeaf() && !dad->isLeaf()) {
        double *pat_lh1 = new double[aln->getNPattern()];
        double *pat_lh2 = new double[aln->getNPattern()];
        double lh1, lh2;
        computeNNIPatternLh(curScore, lh1, pat_lh1, lh2, pat_lh2, node, dad);
        delete[] pat_lh2;
        delete[] pat_lh1;
    }
    FOR_NEIGHBOR_IT(node, dad, it)saveNNITrees((PhyloNode*) (*it)->node, node);
}

void IQTree::summarizeBootstrap(Params &params, MTreeSet &trees) {
    int sum_weights = trees.sumTreeWeights();
    int i, j;
    if (verbose_mode >= VB_MAX) {
        for (i = 0; i < trees.size(); i++)
            if (trees.tree_weights[i] > 0)
                cout << "Tree " << i + 1 << " weight= " << (double) trees.tree_weights[i] * 100 / sum_weights << endl;
    }
    int max_tree_id = max_element(trees.tree_weights.begin(), trees.tree_weights.end()) - trees.tree_weights.begin();
    if (verbose_mode >= VB_MED) {
		cout << "max_tree_id = " << max_tree_id + 1 << "   max_weight = " << trees.tree_weights[max_tree_id];
		cout << " (" << (double) trees.tree_weights[max_tree_id] * 100 / sum_weights << "%)" << endl;
    }
    // assign bootstrap support
    SplitGraph sg;
    SplitIntMap hash_ss;
    // make the taxa name
    vector<string> taxname;
    taxname.resize(leafNum);
    if (boot_splits.empty()) {
        getTaxaName(taxname);
    } else {
        boot_splits.back()->getTaxaName(taxname);
    }
    /*if (!tree.save_all_trees)
     trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
     else
     trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false);
     */
    trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false); // do not sort taxa

    if (verbose_mode >= VB_MED)
    	cout << sg.size() << " splits found" << endl;

    if (!boot_splits.empty()) {
        // check the stopping criterion for ultra-fast bootstrap
        if (computeBootstrapCorrelation() < params.min_correlation)
            cout << "WARNING: bootstrap analysis did not converge. You should rerun with higher number of iterations (-nm option)" << endl;

    }

    sg.scaleWeight(1.0 / trees.sumTreeWeights(), false, 4);
    string out_file;
    out_file = params.out_prefix;
    out_file += ".splits";
    if (verbose_mode >= VB_MED) {
		sg.saveFile(out_file.c_str(), IN_OTHER, true);
		cout << "Split supports printed to star-dot file " << out_file << endl;
    }
    // compute the percentage of appearance
    sg.scaleWeight(100.0, true);
    //	printSplitSet(sg, hash_ss);
    //sg.report(cout);
    cout << "Creating bootstrap support values..." << endl;
    stringstream tree_stream;
    printTree(tree_stream, WT_TAXON_ID | WT_BR_LEN);
    MExtTree mytree;
    mytree.readTree(tree_stream, rooted);
    mytree.assignLeafID();
    mytree.createBootstrapSupport(taxname, trees, sg, hash_ss);

    // now write resulting tree with supports
    tree_stream.seekp(0, ios::beg);
    mytree.printTree(tree_stream);

    // now read resulting tree
    tree_stream.seekg(0, ios::beg);
    freeNode();
    readTree(tree_stream, rooted);
    assignLeafNames();
    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    } else {
		initializeAllPartialLh();
		clearAllPartialLH();
    }

    if (!save_all_trees) {
        out_file = params.out_prefix;
        out_file += ".suptree";

        printTree(out_file.c_str());
        cout << "Tree with assigned bootstrap support written to " << out_file << endl;
    }

    out_file = params.out_prefix;
    out_file += ".splits.nex";
    sg.saveFile(out_file.c_str(), IN_NEXUS, false);
    cout << "Split supports printed to NEXUS file " << out_file << endl;

    /*
     out_file = params.out_prefix;
     out_file += ".supval";
     writeInternalNodeNames(out_file);

     cout << "Support values written to " << out_file << endl;
     */

    if (params.print_ufboot_trees) {
        string filename = params.out_prefix;
        filename += ".ufboot";
        ofstream out(filename.c_str());
        for (i = 0; i < trees.size(); i++) {
            NodeVector taxa;
            // change the taxa name from ID to real name
            trees[i]->getOrderedTaxa(taxa);
            for (j = 0; j < taxa.size(); j++)
                taxa[j]->name = aln->getSeqName(taxa[j]->id);
            // now print to file
            for (j = 0; j < trees.tree_weights[i]; j++)
                trees[i]->printTree(out, WT_NEWLINE);
        }
        out.close();
        cout << "UFBoot trees printed to " << filename << endl;
    }

}

void IQTree::summarizeBootstrap(Params &params) {
	if (verbose_mode >= VB_MED)
		cout << "Summarizing from " << treels.size() << " candidate trees..." << endl;
    MTreeSet trees;
    IntVector tree_weights;
    int sample;
    tree_weights.resize(treels_logl.size(), 0);
    for (sample = 0; sample < boot_trees.size(); sample++)
        tree_weights[boot_trees[sample]]++;
    trees.init(treels, rooted, tree_weights);
    summarizeBootstrap(params, trees);
}

void IQTree::summarizeBootstrap(SplitGraph &sg) {
    MTreeSet trees;
    IntVector tree_weights;
    tree_weights.resize(treels_logl.size(), 0);
    for (int sample = 0; sample < boot_trees.size(); sample++)
        tree_weights[boot_trees[sample]]++;
    trees.init(treels, rooted, tree_weights);
    //SplitGraph sg;
    SplitIntMap hash_ss;
    // make the taxa name
    vector<string> taxname;
    taxname.resize(leafNum);
    getTaxaName(taxname);

    /*if (!tree.save_all_trees)
     trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
     else
     trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false);
     */
    trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false); // do not sort taxa
}

void IQTree::pllConvertUFBootData2IQTree(){
    // duplication_counter
    duplication_counter = pllUFBootDataPtr->duplication_counter;
    //treels_logl
    treels_logl.clear();
    for(int i = 0; i < pllUFBootDataPtr->candidate_trees_count; i++)
        treels_logl.push_back(pllUFBootDataPtr->treels_logl[i]);

    //boot_trees
    boot_trees.clear();
    for(int i = 0; i < params->gbo_replicates; i++)
        boot_trees.push_back(pllUFBootDataPtr->boot_trees[i]);

    //treels
    treels.clear();
    if(pllUFBootDataPtr->candidate_trees_count > 0){
        struct pllHashItem * hItem;
        struct pllHashTable * hTable = pllUFBootDataPtr->treels;
        for (int i = 0; i < hTable->size; ++ i){
            hItem = hTable->Items[i];
            while (hItem){
                string k(hItem->str);
                treels[k] = *((int *)hItem->data);
                hItem = hItem->next;
            }
        }
    }
}

double computeCorrelation(IntVector &ix, IntVector &iy) {

    assert(ix.size() == iy.size());
    DoubleVector x;
    DoubleVector y;

    double mx = 0.0, my = 0.0; // mean value
    int i;
    x.resize(ix.size());
    y.resize(iy.size());
    for (i = 0; i < x.size(); i++) {
        x[i] = ix[i];
        y[i] = iy[i];
        mx += x[i];
        my += y[i];
    }
    mx /= x.size();
    my /= y.size();
    for (i = 0; i < x.size(); i++) {
        x[i] = x[i] / mx - 1.0;
        y[i] = y[i] / my - 1.0;
    }

    double f1 = 0.0, f2 = 0.0, f3 = 0.0;
    for (i = 0; i < x.size(); i++) {
        f1 += (x[i]) * (y[i]);
        f2 += (x[i]) * (x[i]);
        f3 += (y[i]) * (y[i]);
    }
    if (f2 == 0.0 || f3 == 0.0)
        return 1.0;
    return f1 / (sqrt(f2) * sqrt(f3));
}

double IQTree::computeBootstrapCorrelation() {
    if (boot_splits.size() < 2)
        return 0.0;
    IntVector split_supports;
    SplitIntMap split_map;
    int i;
    // collect split supports
    SplitGraph *sg = boot_splits.back();
    SplitGraph *half = boot_splits[(boot_splits.size() - 1) / 2];
    for (i = 0; i < half->size(); i++)
        if (half->at(i)->trivial() == -1) {
            split_map.insertSplit(half->at(i), split_supports.size());
            split_supports.push_back((int) (half->at(i)->getWeight()));
        }

    // collect split supports for new tree collection
    IntVector split_supports_new;
    split_supports_new.resize(split_supports.size(), 0);
    for (i = 0; i < sg->size(); i++)
        if ((*sg)[i]->trivial() == -1) {
            int index;
            Split *sp = split_map.findSplit((*sg)[i], index);
            if (sp) {
                // split found
                split_supports_new[index] = (int) ((*sg)[i]->getWeight());
            } else {
                // new split
                split_supports_new.push_back((int) ((*sg)[i]->getWeight()));
            }
        }
    if (verbose_mode >= VB_MED)
    	cout << split_supports_new.size() - split_supports.size() << " new splits compared to old boot_splits" << endl;
    if (split_supports_new.size() > split_supports.size())
        split_supports.resize(split_supports_new.size(), 0);

    // now compute correlation coefficient
    double corr = computeCorrelation(split_supports, split_supports_new);
    // printing supports into file
    /*
     string outfile = params->out_prefix;
     outfile += ".splitsup";
     try {
     ofstream out;
     out.exceptions(ios::failbit | ios::badbit);
     out.open(outfile.c_str());
     out << "tau=" << max_candidate_trees / 2 << "\ttau="
     << treels_logl.size() << endl;
     for (int i = 0; i < split_supports.size(); i++)
     out << split_supports[i] << "\t" << split_supports_new[i] << endl;
     out.close();
     cout << "Split support values printed to " << outfile << endl;
     } catch (ios::failure) {
     outError(ERR_WRITE_OUTPUT, outfile);
     }
     */
    return corr;
}

void IQTree::addPositiveNNIMove(NNIMove myMove) {
    plusNNIs.push_back(myMove);
}

void IQTree::printResultTree(string suffix) {
    setRootNode(params->root);
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    if (suffix.compare("") != 0) {
        string iter_tree_name = tree_file_name + "." + suffix;
        printTree(iter_tree_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    } else {
        printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    }
    //printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH);
}

void IQTree::printResultTree(ostream &out) {
    setRootNode(params->root);
    printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
}

/*
void IQTree::printPhylolibModelParams(const char* suffix) {
    char phyloliModelFile[1024];
    strcpy(phyloliModelFile, params->out_prefix);
    strcat(phyloliModelFile, suffix);
    ofstream modelfile;
    modelfile.open(phyloliModelFile);
    for (int model = 0; model < pllInst->NumberOfModels; model++) {
        cout << "Rate parameters: ";
        for (int i = 0; i < 6; i++) {
            cout << pllInst->partitionData[model].substRates[i] << " ";
            modelfile << pllInst->partitionData[model].substRates[i] << " ";
        }
        cout << endl;
        modelfile << endl;
        cout << "Base frequencies: ";
        for (int i = 0; i < aln->num_states; i++) {
            cout << pll_tree->partitionData[model].frequencies[i] << " ";
            modelfile << pll_tree->partitionData[model].frequencies[i] << " ";
        }
        cout << endl;
        modelfile << endl;
        cout << "Gamma shape :" << pll_tree->partitionData[model].alpha << endl;
        modelfile << pll_tree->partitionData[model].alpha << endl;
    }
}
*/

void IQTree::printPhylolibTree(const char* suffix) {
    pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_TRUE, 1, 0, 0, 0,
            PLL_SUMMARIZE_LH, 0, 0);
    char phylolibTree[1024];
    strcpy(phylolibTree, params->out_prefix);
    strcat(phylolibTree, suffix);
    FILE *phylolib_tree = fopen(phylolibTree, "w");
    fprintf(phylolib_tree, "%s", pllInst->tree_string);
    cout << "Tree optimized by Phylolib was written to " << phylolibTree << endl;
}

void IQTree::printIntermediateTree(int brtype) {
    setRootNode(params->root);
    bool duplicated_tree = false;
    double *pattern_lh = NULL;
    double logl = curScore;
    if (params->avoid_duplicated_trees) {
        // estimate logl_cutoff
        stringstream ostr;
        printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
        string tree_str = ostr.str();
        StringIntMap::iterator it = treels.find(tree_str);
        if (it != treels.end()) { // already in treels
            duplicated_tree = true;
            if (curScore > treels_logl[it->second] + 1e-4) {
                if (verbose_mode >= VB_MAX)
                    cout << "Updated logl " << treels_logl[it->second] << " to " << curScore << endl;
                treels_logl[it->second] = curScore;
                computeLikelihood(treels_ptnlh[it->second]);
                if (save_all_br_lens) {
                    ostr.seekp(ios::beg);
                    printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
                    treels_newick[it->second] = ostr.str();
                }
            }
            //pattern_lh = treels_ptnlh[treels[tree_str]];
        } else {
            //cout << __func__ << ": new tree" << endl;
            if (logl_cutoff != 0.0 && curScore <= logl_cutoff + 1e-4)
                duplicated_tree = true;
            else {
                treels[tree_str] = treels_ptnlh.size();
                pattern_lh = new double[aln->getNPattern()];
                computePatternLikelihood(pattern_lh, &logl);
                treels_ptnlh.push_back(pattern_lh);
                treels_logl.push_back(logl);
                if (save_all_br_lens) {
                    ostr.seekp(ios::beg);
                    printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
                    treels_newick.push_back(ostr.str());
                }
            }
        }
        //cout << tree_str << endl;
    } else {
        if (params->print_tree_lh) {
            pattern_lh = new double[aln->getNPattern()];
            computePatternLikelihood(pattern_lh, &logl);
        }
    }

    if (!duplicated_tree) {
        if (write_intermediate_trees)
            printTree(out_treels, brtype);
        if (params->print_tree_lh) {
            out_treelh.precision(10);
            out_treelh << logl;
            double prob;
            aln->multinomialProb(pattern_lh, prob);
            out_treelh << "\t" << prob << endl;
            if (!(brtype & WT_APPEND))
                out_sitelh << aln->getNSite() << endl;
            out_sitelh << "Site_Lh   ";
            for (int i = 0; i < aln->getNSite(); i++)
                out_sitelh << "\t" << pattern_lh[aln->getPatternID(i)];
            out_sitelh << endl;
            if (!params->avoid_duplicated_trees)
                delete[] pattern_lh;
        }
    }
    if (params->write_intermediate_trees == 1 && save_all_trees != 1) {
        return;
    }
    int x = save_all_trees;
    save_all_trees = 2;
    evalNNIs();
    save_all_trees = x;
}

void IQTree::removeIdenticalSeqs(Params &params, StrVector &removed_seqs, StrVector &twin_seqs) {
	Alignment *new_aln;
	if (params.root)
		new_aln = aln->removeIdenticalSeq((string)params.root, removed_seqs, twin_seqs);
	else
		new_aln = aln->removeIdenticalSeq("", removed_seqs, twin_seqs);
	if (removed_seqs.size() > 0)
		cout << "INFO: " << removed_seqs.size() << " identical sequences are removed for tree reconstruction" << endl;
	if (new_aln != aln) {
		aln = new_aln;
	}
}

void IQTree::reinsertIdenticalSeqs(Alignment *orig_aln, StrVector &removed_seqs, StrVector &twin_seqs) {
	if (removed_seqs.empty()) return;
	IntVector id;
	int i;
	id.resize(removed_seqs.size());
	for (i = 0; i < id.size(); i++)
		id[i] = i;
	// randomize order before reinsert back into tree
//	random_shuffle(id.begin(), id.end(), random_int);
	random_shuffle(id.begin(), id.end());

	for (int i = 0; i < removed_seqs.size(); i++) {
		Node *old_taxon = findLeafName(twin_seqs[id[i]]);
		assert(old_taxon);
		double len = old_taxon->neighbors[0]->length;
		Node *old_node = old_taxon->neighbors[0]->node;
		Node *new_taxon = newNode(leafNum+i, removed_seqs[id[i]].c_str());
		Node *new_node = newNode();
		// link new_taxon - new_node
		new_taxon->addNeighbor(new_node, 0.0);
		new_node->addNeighbor(new_taxon, 0.0);
		// link old_taxon - new_node
		new_node->addNeighbor(old_taxon, 0.0);
		old_taxon->updateNeighbor(old_node, new_node, 0.0);
		// link old_node - new_node
		new_node->addNeighbor(old_node, len);
		old_node->updateNeighbor(old_taxon, new_node, len);
	}

    leafNum = leafNum + removed_seqs.size();
    initializeTree();
//    delete iqtree.aln;
    setAlignment(orig_aln);
    // delete all partial_lh, which will be automatically recreated later
    deleteAllPartialLh();
}
