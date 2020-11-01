/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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
#include "phylotreemixlen.h"
#include "phylosupertree.h"
#include "phylosupertreeplen.h"
#include "model/partitionmodelplen.h"
#include "model/modelfactorymixlen.h"
#include "mexttree.h"
#include "utils/timeutil.h"
#include "model/modelmarkov.h"
#include "model/rategamma.h"
//#include "phylotreemixlen.h"
//#include "model/modelfactorymixlen.h"
#include <numeric>
#include "utils/tools.h"
#include "utils/MPIHelper.h"
#include "utils/pllnni.h"

Params *globalParams;
Alignment *globalAlignment;
extern StringIntMap pllTreeCounter;

IQTree::IQTree() : PhyloTree() {
    IQTree::init();
}

void IQTree::init() {
//    PhyloTree::init();
    k_represent = 0;
    k_delete = k_delete_min = k_delete_max = k_delete_stay = 0;
    dist_matrix = NULL;
    var_matrix = NULL;
//    curScore = 0.0; // Current score of the tree
    cur_pars_score = -1;
//    enable_parsimony = false;
    estimate_nni_cutoff = false;
    nni_cutoff = -1e6;
    nni_sort = false;
    testNNI = false;
//    print_tree_lh = false;
//    write_intermediate_trees = 0;
//    max_candidate_trees = 0;
    logl_cutoff = 0.0;
    len_scale = 10000;
//    save_all_br_lens = false;
    duplication_counter = 0;
    //boot_splits = new SplitGraph;
    pll2iqtree_pattern_index = NULL;

    treels_name = Params::getInstance().out_prefix;
    treels_name += ".treels";
    out_lh_file = Params::getInstance().out_prefix;
    out_lh_file += ".treelh";
    site_lh_file = Params::getInstance().out_prefix;
    site_lh_file += ".sitelh";

    if (Params::getInstance().print_tree_lh) {
        out_treelh.open(out_lh_file.c_str());
        out_sitelh.open(site_lh_file.c_str());
    }

    if (Params::getInstance().write_intermediate_trees)
        out_treels.open(treels_name.c_str());
    on_refine_btree = false;
    contree_rfdist = -1;
    boot_consense_logl = 0.0;

}

IQTree::IQTree(Alignment *aln) : PhyloTree(aln) {
    IQTree::init();
}

void IQTree::setCheckpoint(Checkpoint *checkpoint) {
    PhyloTree::setCheckpoint(checkpoint);
    stop_rule.setCheckpoint(checkpoint);
    candidateTrees.setCheckpoint(checkpoint);
    for (auto it = boot_splits.begin(); it != boot_splits.end(); it++)
        (*it)->setCheckpoint(checkpoint);
}

void IQTree::saveUFBoot(Checkpoint *checkpoint) {
    checkpoint->startStruct("UFBoot");
    if (MPIHelper::getInstance().isWorker()) {
        CKP_SAVE(sample_start);
        CKP_SAVE(sample_end);
        checkpoint->startList(boot_samples.size());
        checkpoint->setListElement(sample_start-1);
        for (int id = sample_start; id != sample_end; id++) {
            checkpoint->addListElement();
            stringstream ss;
            ss.precision(10);
            ss << boot_counts[id] << " " << boot_logl[id] << " " << boot_orig_logl[id] << " " << boot_trees[id];
            checkpoint->put("", ss.str());
        }
        checkpoint->endList();
    } else {
        CKP_SAVE(logl_cutoff);
        int boot_splits_size = boot_splits.size();
        CKP_SAVE(boot_splits_size);
        checkpoint->startList(boot_samples.size());
        for (int id = 0; id != boot_samples.size(); id++) {
            checkpoint->addListElement();
            stringstream ss;
            ss.precision(10);
            ss << boot_counts[id] << " " << boot_logl[id] << " " << boot_orig_logl[id] << " " << boot_trees[id];
            checkpoint->put("", ss.str());
        }
        checkpoint->endList();
    }
    checkpoint->endStruct();
}

void IQTree::saveCheckpoint() {
    stop_rule.saveCheckpoint();
    candidateTrees.saveCheckpoint();
    
    if (boot_samples.size() > 0 && !boot_trees.front().empty()) {
        saveUFBoot(checkpoint);
        // boot_splits
        int id = 0;
        for (vector<SplitGraph*>::iterator sit = boot_splits.begin(); sit != boot_splits.end(); sit++, id++) {
            checkpoint->startStruct("UFBootSplit" + convertIntToString(id));
            (*sit)->saveCheckpoint();
            checkpoint->endStruct();
        }
    }
    
    PhyloTree::saveCheckpoint();
    CKP_SAVE(boot_consense_logl);
    CKP_SAVE(contree_rfdist);
}

void IQTree::restoreUFBoot(Checkpoint *checkpoint) {
    checkpoint->startStruct("UFBoot");
    // save boot_samples and boot_trees
    int id;
    checkpoint->startList(params->gbo_replicates);
    int sample_start, sample_end;
    CKP_RESTORE(sample_start);
    CKP_RESTORE(sample_end);
    checkpoint->setListElement(sample_start-1);
    for (id = sample_start; id != sample_end; id++) {
        checkpoint->addListElement();
        string str;
        checkpoint->getString("", str);
        ASSERT(!str.empty());
        stringstream ss(str);
        ss >> boot_counts[id] >> boot_logl[id] >> boot_orig_logl[id] >> boot_trees[id];
    }
    checkpoint->endList();
    checkpoint->endStruct();
}

void IQTree::restoreCheckpoint() {
    PhyloTree::restoreCheckpoint();
    stop_rule.restoreCheckpoint();
    candidateTrees.restoreCheckpoint();

    if (params->gbo_replicates > 0 && checkpoint->hasKeyPrefix("UFBoot")) {
        checkpoint->startStruct("UFBoot");
//        CKP_RESTORE(max_candidate_trees);
        CKP_RESTORE(logl_cutoff);
        // save boot_samples and boot_trees
        int id = 0;
        checkpoint->startList(params->gbo_replicates);
        boot_trees.resize(params->gbo_replicates);
        boot_logl.resize(params->gbo_replicates);
        boot_orig_logl.resize(params->gbo_replicates);
        boot_counts.resize(params->gbo_replicates);
        for (id = 0; id < params->gbo_replicates; id++) {
            checkpoint->addListElement();
            string str;
            checkpoint->getString("", str);
            stringstream ss(str);
            ss >> boot_counts[id] >> boot_logl[id] >> boot_orig_logl[id] >> boot_trees[id];
        }
        checkpoint->endList();
        int boot_splits_size = 0;
        CKP_RESTORE(boot_splits_size);
        checkpoint->endStruct();

        // boot_splits
        for (id = 0; id < boot_splits_size; id++) {
            checkpoint->startStruct("UFBootSplit" + convertIntToString(id));
            SplitGraph *sg = new SplitGraph;
            sg->createBlocks();
            StrVector taxname;
            getTaxaName(taxname);
            for (auto its = taxname.begin(); its != taxname.end(); its++)
                sg->getTaxa()->AddTaxonLabel(NxsString(its->c_str()));
            sg->setCheckpoint(checkpoint);
            sg->restoreCheckpoint();
            boot_splits.push_back(sg);
            checkpoint->endStruct();
        }
    }

    CKP_RESTORE(boot_consense_logl);
    CKP_RESTORE(contree_rfdist);
}

void IQTree::initSettings(Params &params) {

    searchinfo.nni_type = params.nni_type;
    optimize_by_newton = params.optimize_by_newton;
    setLikelihoodKernel(params.SSE);
    if (num_threads > 0) {
        setNumThreads(num_threads);
    } else {
        setNumThreads(params.num_threads);
    }
    candidateTrees.init(this->aln, 200);
    intermediateTrees.init(this->aln, 200000);

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
    if (!params.treeset_file.empty() && params.min_iterations == -1) {
        params.min_iterations = 1;
        params.stop_condition = SC_FIXED_ITERATION;
        params.numInitTrees = 1;
    }
    if (params.gbo_replicates)
        params.max_iterations = max(params.max_iterations, max(params.min_iterations, params.step_iterations));

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

    globalParams = &params;
    globalAlignment = aln;

    //write_intermediate_trees = params.write_intermediate_trees;

    if (Params::getInstance().write_intermediate_trees > 2 || params.gbo_replicates > 0) {
        save_all_trees = 1;
    }
    if (params.gbo_replicates > 0) {
        if (params.iqp_assess_quartet != IQP_BOOTSTRAP) {
            save_all_trees = 2;
        }
    }
//    if (params.gbo_replicates > 0 && params.do_compression)
//        save_all_br_lens = true;
//    print_tree_lh = params.print_tree_lh;
//    max_candidate_trees = params.max_candidate_trees;
//    if (max_candidate_trees == 0)
//        max_candidate_trees = aln->getNSeq() * params.step_iterations;
    if (!params.compute_ml_tree_only) {
        setRootNode(params.root);
    }

    size_t i;

    if (params.online_bootstrap && params.gbo_replicates > 0 && !isSuperTreeUnlinked()) {
        if (aln->getNSeq() < 4)
            outError("It makes no sense to perform bootstrap with less than 4 sequences.");

        string bootaln_name = params.out_prefix;
        bootaln_name += ".bootaln";
        if (params.print_bootaln) {
            ofstream bootalnout;
            bootalnout.open(bootaln_name.c_str());
            bootalnout.close();
        }

        // 2015-12-17: initialize random stream for creating bootstrap samples
        // mainly so that checkpointing does not need to save bootstrap samples
        int *saved_randstream = randstream;
        init_random(params.ran_seed);
        
//        cout << "Generating " << params.gbo_replicates << " samples for ultrafast "
//             << RESAMPLE_NAME << " (seed: " << params.ran_seed << ")..." << endl;
        // allocate memory for boot_samples
        boot_samples.resize(params.gbo_replicates);
        sample_start = 0;
        sample_end = boot_samples.size();

        // compute the sample_start and sample_end
        if (MPIHelper::getInstance().getNumProcesses() > 1) {
            int num_samples = boot_samples.size() / MPIHelper::getInstance().getNumProcesses();
            if (boot_samples.size() % MPIHelper::getInstance().getNumProcesses() != 0)
                num_samples++;
            sample_start = MPIHelper::getInstance().getProcessID() * num_samples;
            sample_end = sample_start + num_samples;
            if (sample_end > boot_samples.size())
                sample_end = boot_samples.size();
        }

        size_t orig_nptn = getAlnNPattern();
#ifdef BOOT_VAL_FLOAT
        size_t nptn = get_safe_upper_limit_float(orig_nptn);
#else
        size_t nptn = get_safe_upper_limit(orig_nptn);
#endif
        BootValType *mem = aligned_alloc<BootValType>(nptn * (size_t)(params.gbo_replicates));
        memset(mem, 0, nptn * (size_t)(params.gbo_replicates) * sizeof(BootValType));
        for (i = 0; i < params.gbo_replicates; i++) {
            boot_samples[i] = mem + i*nptn;
        }
        if (boot_trees.empty()) {
            boot_logl.resize(params.gbo_replicates, -DBL_MAX);
            boot_orig_logl.resize(params.gbo_replicates, -DBL_MAX);
            boot_trees.resize(params.gbo_replicates, "");
            boot_counts.resize(params.gbo_replicates, 0);
        } else {
            cout << "CHECKPOINT: " << boot_trees.size() << " UFBoot trees and " << boot_splits.size() << " UFBootSplits restored" << endl;
        }
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
                for (size_t j = 0; j < orig_nptn; j++)
                    boot_samples[i][j] = this_sample[j];
                bootstrap_alignment->printAlignment(params.aln_output_format, bootaln_name.c_str(), true);
                delete bootstrap_alignment;
            } else {
                IntVector this_sample;
                aln->createBootstrapAlignment(this_sample, params.bootstrap_spec);
                for (size_t j = 0; j < orig_nptn; j++)
                    boot_samples[i][j] = this_sample[j];
            }
        }
        verbose_mode = saved_mode;
        if (params.print_bootaln) {
            cout << "Bootstrap alignments printed to " << bootaln_name << endl;
        }

//        cout << "Max candidate trees (tau): " << max_candidate_trees << endl;
        
        // restore randstream
        finish_random();
        randstream = saved_randstream;

        // Diep: initialize data members to be used in the Refinement Step
        on_refine_btree = false;
        saved_aln_on_refine_btree = NULL;
        if(params.ufboot2corr){
            boot_samples_int.resize(params.gbo_replicates);
            for (size_t i = 0; i < params.gbo_replicates; i++) {
                boot_samples_int[i].resize(nptn, 0);
                for (size_t j = 0; j < orig_nptn; j++)
                    boot_samples_int[i][j] = boot_samples[i][j];
               }
        }
    }
    if (params.root_state) {
        if (strlen(params.root_state) != 1)
            outError("Root state must have exactly 1 character");
        root_state = aln->convertState(params.root_state[0]);
        if (root_state < 0 || root_state >= aln->num_states)
            outError("Invalid root state");
    }
}

IQTree::~IQTree() {
    //if (bonus_values)
    //delete bonus_values;
    //bonus_values = NULL;

    for (vector<SplitGraph*>::reverse_iterator it2 = boot_splits.rbegin(); it2 != boot_splits.rend(); it2++)
        delete (*it2);
    boot_splits.clear();
    //if (boot_splits) delete boot_splits;

    if (!boot_samples.empty()) {
        aligned_free(boot_samples[0]); // free memory
        boot_samples.clear();
    }
}

extern const char *aa_model_names_rax[];

void IQTree::createPLLPartition(Params &params, ostream &pllPartitionFileHandle) {
    if (isSuperTree()) {
        PhyloSuperTree *siqtree = (PhyloSuperTree*) this;
        // additional check for PLL hard limit
        if (params.pll) {
            if (siqtree->size() > PLL_NUM_BRANCHES)
                outError("Number of partitions exceeds PLL limit, please increase PLL_NUM_BRANCHES constant in pll.h");
            int i = 0;
            int startPos = 1;
            
            // prepare proper partition file 
            for (PhyloSuperTree::iterator it = siqtree->begin(); it != siqtree->end(); it++) {
                i++;
                int curLen = ((*it)->aln->seq_type == SEQ_CODON) ? (*it)->getAlnNSite()*3 : (*it)->getAlnNSite();
                if ((*it)->aln->seq_type == SEQ_DNA || (*it)->aln->seq_type == SEQ_CODON) {
                    pllPartitionFileHandle << "DNA";
                } else if ((*it)->aln->seq_type == SEQ_PROTEIN) {
                    if ((*it)->aln->model_name != "" && (*it)->aln->model_name.substr(0, 4) != "TEST" && (*it)->aln->model_name.substr(0, 2) != "MF") {
                        string modelStr = (*it)->aln->model_name.
                                substr(0, (*it)->aln->model_name.find_first_of("+{"));
                        if (modelStr == "LG4")
                            modelStr = "LG4M";
                        bool name_ok = false;
                        for (int j = 0; j < 18; j++)
                            if (modelStr == aa_model_names_rax[j]) {
                                name_ok = true;
                                break;
                            }
                        if (name_ok)
                            pllPartitionFileHandle << modelStr;
                        else
                            pllPartitionFileHandle << "WAG";                    
                    } else {
                        pllPartitionFileHandle << "WAG";
                    }
                } else
                    outError("PLL only works with DNA/protein alignments");
                pllPartitionFileHandle << ", p" << i << " = " << startPos << "-" << startPos + curLen - 1 << endl;
                startPos = startPos + curLen;
            }
        } else {
            // only prepare partition file for computing parsimony trees
            SeqType datatype[] = {SEQ_DNA, SEQ_CODON, SEQ_PROTEIN};
            PhyloSuperTree::iterator it;
            
            for (int i = 0; i < sizeof(datatype)/sizeof(SeqType); i++) {
                bool first = true;
                int startPos = 1;
                for (it = siqtree->begin(); it != siqtree->end(); it++) 
                    if ((*it)->aln->seq_type == datatype[i]) {
                        if (first) {
                        if (datatype[i] == SEQ_DNA || datatype[i] == SEQ_CODON)
                            pllPartitionFileHandle << "DNA";
                        else
                            pllPartitionFileHandle << "WAG";
                        }
                        int curLen = (datatype[i] == SEQ_CODON) ? (*it)->getAlnNSite()*3 : (*it)->getAlnNSite();
                        if (first)
                            pllPartitionFileHandle << ", p" << i << " = ";
                        else
                            pllPartitionFileHandle << ", ";
                            
                        pllPartitionFileHandle << startPos << "-" << startPos + curLen - 1;
                        startPos = startPos + curLen;
                        first = false;
                    } else {
                        startPos = startPos + (*it)->getAlnNSite();
                    }
                if (!first) pllPartitionFileHandle << endl;
            }
        }
    } else {
        /* create a partition file */
        string model;
        if (aln->seq_type == SEQ_DNA || aln->seq_type == SEQ_CODON) {
            model = "DNA";
        } else if (aln->seq_type == SEQ_PROTEIN) {
            if (params.pll && params.model_name != "" && params.model_name.substr(0, 4) != "TEST" && params.model_name.substr(0, 2) != "MF") {
                model = params.model_name.substr(0, params.model_name.find_first_of("+{"));
            } else {
                model = "WAG";
            }
        } else {
            model = "WAG";
            //outError("PLL currently only supports DNA/protein alignments");
        }
        size_t nsite = (aln->seq_type == SEQ_CODON) ? getAlnNSite()*3 : getAlnNSite();
        pllPartitionFileHandle << model << ", p1 = " << "1-" << nsite << endl;
    }
}

void IQTree::computeInitialTree(LikelihoodKernel kernel) {
    double start = getRealTime();
    string initTree;
    string out_file = params->out_prefix;
    int score;
    if ( params->stop_condition == SC_FIXED_ITERATION
         && params->numNNITrees > params->min_iterations) {
        params->numNNITrees = max(params->min_iterations, 1);
    }
    int fixed_number = 0;
    if (params->sankoff_cost_file && !isUsingSankoffParsimony()) {
        loadCostMatrixFile(params->sankoff_cost_file);
    }
    if (aln->ordered_pattern.empty()) {
        aln->orderPatternByNumChars(PAT_VARIANT);
    }
    setParsimonyKernel(kernel);
    bool noisy = !progress_display::getProgressDisplay();
    if (params->user_file) {
        // start the search with user-defined tree
        if (noisy) {
            cout << "Reading input tree file "
                << params->user_file << " ...";
        }
        bool myrooted = params->is_rooted;
        readTree(params->user_file, myrooted);
        if (myrooted && !isSuperTreeUnlinked() && noisy) {
            cout << " rooted tree";
        }
        if (noisy) {
            cout << endl;
        }
        //
        //Todo: In incremental mode, if the user-defined tree
        //      has leaves with names *not* found in the alignment, they
        //      need to be removed.  If it is *missing* leaves for sequences
        //      in the alignment, they will need to be added somehow.
        //      Perhaps by treating the user-defined tree as a constraint
        //      tree? The code for all that would go... here.
        //
        if (params->incremental) {
            updateToMatchAlignment(aln);
        } else {
            setAlignment(aln);
        }
        if (isSuperTree()) {
            wrapperFixNegativeBranch(params->fixed_branch_length == BRLEN_OPTIMIZE &&
                                     params->partition_type == BRLEN_OPTIMIZE);
        }
        else {
            fixed_number = wrapperFixNegativeBranch(false);
        }
        params->numInitTrees = 1;
        params->numNNITrees = 1;
        if (params->pll) {
            pllReadNewick(getTreeString());
        }
        initTree = getTreeString();
        CKP_SAVE(initTree);
        saveCheckpoint();
    } else if (CKP_RESTORE(initTree)) {
        bool was_restored = true;
        readTreeString(initTree);
        if (params->incremental) {
            was_restored = !updateToMatchAlignment(aln);
        }
        if (noisy && was_restored) {
            cout << endl << "CHECKPOINT: Initial tree restored" << endl;
        }
    } else {
        START_TREE_TYPE start_tree = params->start_tree;
        // only own parsimony kernel supports constraint tree
        if (!constraintTree.empty()) {
            start_tree = STT_PARSIMONY;
        }
        switch (start_tree) {
        case STT_PARSIMONY:
            //initCostMatrix(CM_UNIFORM);
            // Create parsimony tree using IQ-Tree kernel
            cout << "Creating fast initial parsimony tree by random order stepwise addition..." << endl;
            start = getRealTime();
            score = computeParsimonyTree(params->out_prefix, aln, randstream);
            cout << "Step-wise parsimony addition took " << getRealTime() - start << " seconds"
                <<", parsimony score: " << score
                << " (based on " << aln->num_parsimony_sites << " sites)"<< endl;
            break;
        case STT_PARSIMONY_JOINING:
            cout << "Creating parsimony tree by parsimony joining..." << endl;
            start = getRealTime();
            score = joinParsimonyTree(params->out_prefix, aln);
            cout << "Parsimony joining took " << getRealTime() - start << " seconds"
                << ", parsimony score: " << score
                << " (based on " << aln->num_parsimony_sites << " sites)"<< endl;
            break;

        case STT_RANDOM_TREE:
        case STT_PLL_PARSIMONY:
            cout << endl;
            cout << "Create initial parsimony tree by phylogenetic likelihood library (PLL)... ";
            pllInst->randomNumberSeed = params->ran_seed + MPIHelper::getInstance().getProcessID();
            pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInst, pllPartitions, params->sprDist);
            resetBranches(pllInst);
            pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back,
                    PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
            PhyloTree::readTreeStringSeqName(string(pllInst->tree_string));
            cout << getRealTime() - start << " seconds" << endl;
            wrapperFixNegativeBranch(true);
            break;
        case STT_BIONJ:
            // This is the old default option: using BIONJ as starting tree
            if (this->dist_matrix == nullptr) {
                //Really, we want to re-use the distance matrix from the tree,
                //two levels higher in the call stack (in runModelFinder)
                //(and this is just an easy workaround to avert a null de-reference).
                //Easier than adding parameters to computeFastMLTree, and to
                //IQTree::computeInitialTree, but not ideal. -James B. 17-Aug-2020
                computeDist(*params, aln, dist_matrix, var_matrix);
            }
            computeBioNJ(*params);
            if (verbose_mode >= VB_MED) {
                cout << "Computing initial tree took " << getRealTime() - start
                << " wall-clock seconds" << endl;
            }
            params->numInitTrees = 1;
            if (isSuperTree()) {
                wrapperFixNegativeBranch(true);
            } else {
                fixed_number = wrapperFixNegativeBranch(false);
            }
            break;
        case STT_USER_TREE:
            ASSERT(0 && "User tree should be handled already");
            break;
        }
        initTree = getTreeString();
        CKP_SAVE(initTree);
        saveCheckpoint();
        checkpoint->dump();
    }
    if (!constraintTree.isCompatible(this)) {
        outError("Initial tree is not compatible with constraint tree");
    }
    if (fixed_number) {
        cout << "WARNING: " << fixed_number << " undefined/negative"
            << " branch lengths are initialized with parsimony" << endl;
    }
    if (params->root) {
        StrVector outgroup_names;
        convert_string_vec(params->root, outgroup_names);
        for (auto it = outgroup_names.begin(); it != outgroup_names.end(); it++) {
            if (aln->getSeqID(*it) < 0) {
                outError("Specified outgroup taxon " + *it + " not found");
            }
        }
    }
    if (params->write_init_tree) {
        out_file += ".init_tree";
        printTree(out_file.c_str(), WT_NEWLINE);
    }
}

int IQTree::addTreeToCandidateSet(string treeString, double score, bool updateStopRule, int sourceProcID) {
    double curBestScore = candidateTrees.getBestScore();
    int pos = candidateTrees.update(treeString, score);
    if (updateStopRule) {
        stop_rule.setCurIt(stop_rule.getCurIt() + 1);
        if (score > curBestScore) {
            hideProgress();
            if (pos != -1) {
                stop_rule.addImprovedIteration(stop_rule.getCurIt());
                cout << "BETTER TREE FOUND at iteration " << stop_rule.getCurIt() << ": " << score << endl;
            } else {
                cout << "UPDATE BEST LOG-LIKELIHOOD: " << score << endl;
            }
            showProgress();
            bestcandidate_changed = true;
            // COMMENT OUT: not safe with MPI version
//            printResultTree();
        }

        curScore = score;
        printIterationInfo(sourceProcID);
    }
    return pos;
}

void IQTree::initCandidateTreeSet(int nParTrees, int nNNITrees) {
    if (nParTrees > 0 && !progress_display::getProgressDisplay()) {
        if (params->start_tree == STT_RANDOM_TREE) {
            cout << "Generating " << nParTrees  << " random trees... ";
            cout.flush();
        }
        else if (nParTrees > 1) {
            cout << "Generating " << nParTrees  << " parsimony trees... ";
            cout.flush();
        }
    }
    double startTime = getRealTime();
//    int numDupPars = 0;
//    bool orig_rooted = rooted;
//    rooted = false;
    int processID = MPIHelper::getInstance().getProcessID();

    int init_size = candidateTrees.size();
    std::string whatAmIDoingHere;
    const char* verb = "loaded";
    switch (params->start_tree) {
        case STT_PLL_PARSIMONY:
            whatAmIDoingHere = "Constructing parsimony trees via PLL";
            verb = "constructed";
            break;
        case STT_RANDOM_TREE:
            whatAmIDoingHere = "Generating Yule-Harding random trees";
            verb = "generated";
            break;
        case STT_PARSIMONY:
            whatAmIDoingHere = "Constructing parsimony trees";
            verb = "constructed";
            break;
        case STT_PARSIMONY_JOINING:
            whatAmIDoingHere = "Constructing parsimony tree via parsimony joining";
            verb = "constructed";
            break;
        case STT_BIONJ:
            whatAmIDoingHere = std::string("Loading ") + params->start_tree_subtype_name + " tree ";
            break;
        case STT_USER_TREE:
            /* fall-through */
        default:
            whatAmIDoingHere = "Loading user tree";
            break;
    }
    if (0<nParTrees) {
        initProgress(nParTrees, whatAmIDoingHere, verb, "tree");
        for (int treeNr = 1; treeNr <= nParTrees; ++treeNr) {
            int parRandSeed = Params::getInstance().ran_seed + processID * nParTrees + treeNr;
            string curParsTree;
            
            /********* Create parsimony tree using PLL *********/
            if (params->start_tree == STT_PLL_PARSIMONY) {
                pllInst->randomNumberSeed = parRandSeed;
                pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInst, pllPartitions, params->sprDist);
                resetBranches(pllInst);
                pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions,
                                pllInst->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE,
                                PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
                curParsTree = string(pllInst->tree_string);
                PhyloTree::readTreeStringSeqName(curParsTree);
                wrapperFixNegativeBranch(true);
            } else if (params->start_tree == STT_RANDOM_TREE) {
                generateRandomTree(YULE_HARDING);
                wrapperFixNegativeBranch(true);
                if (rooted) {
                    rooted = false;
                    convertToRooted();
                }
            } else if (params->start_tree == STT_PARSIMONY) {
                int *rstream;
                init_random(parRandSeed, false, &rstream);
                PhyloTree tree;
                if (!constraintTree.empty()) {
                    tree.constraintTree.readConstraint(constraintTree);
                }
                tree.showNoProgress();
                tree.setParams(params);
                tree.setParsimonyKernel(params->SSE);
                tree.rooted = rooted;
                tree.computeParsimonyTree(NULL, aln, rstream);
                finish_random(rstream);
                PhyloTree::readTreeString(tree.getTreeString());
            } else if (params->start_tree == STT_PARSIMONY_JOINING) {
                joinParsimonyTree(nullptr, aln);
            } else {
                //Use the tree we've already got!
            }
            curParsTree = getTreeString();
            int pos = addTreeToCandidateSet(curParsTree, -DBL_MAX, false, MPIHelper::getInstance().getProcessID());
            // if a duplicated tree is generated, then randomize the tree
            if (pos == -1) {
                readTreeString(curParsTree);
                doRandomNNIs();
                wrapperFixNegativeBranch(true);
                string randTree = getTreeString();
                addTreeToCandidateSet(randTree, -DBL_MAX, false, MPIHelper::getInstance().getProcessID());
            }
            trackProgress(1);
        }
        doneProgress();
    }
    if (!progress_display::getProgressDisplay()) {
        cout << getRealTime() - startTime << " second" << endl;
    }
    /****************************************************************************************
                          Compute logl of all initial trees
    *****************************************************************************************/

    vector<string> initTreeStrings = candidateTrees.getBestTreeStrings();
    candidateTrees.clear();
    
    size_t candidateCount = initTreeStrings.size();
    bool   sayHowLong = init_size < initTreeStrings.size() && !progress_display::getProgressDisplay();
    if (sayHowLong) {
        cout << "Computing log-likelihood of " << initTreeStrings.size() - init_size
            << " initial trees ... ";
    }
    startTime = getRealTime();

    //(b) Initialize candidate trees
    initProgress(candidateCount, "Assessing candidate trees", "assessed", "tree");
    for (vector<string>::iterator it = initTreeStrings.begin(); it != initTreeStrings.end(); ++it) {
        string treeString;
        double score;
        readTreeString(*it);
        if (it-initTreeStrings.begin() >= init_size) {
            treeString = optimizeBranches(params->brlen_num_traversal);
        }
        else {
            computeLogL();
            treeString = getTreeString();
        }
        score = getCurScore();
        candidateTrees.update(treeString,score);
        trackProgress(1.0);
    }
    doneProgress();

    if (Params::getInstance().writeDistImdTrees) {
        intermediateTrees.initTrees(candidateTrees);
    }
    if (sayHowLong) {
        cout << getRealTime() - startTime << " seconds" << endl;
    }
    if (nParTrees > 0) {
        cout << "Current best score: " << candidateTrees.getBestScore() << endl;
    }

    //---- BLOCKING COMMUNICATION
    syncCandidateTrees(nNNITrees, false);

    vector<string> bestInitTrees; // Set of best initial trees for doing NNIs
    bestInitTrees = candidateTrees.getBestTreeStringsForProcess(nNNITrees);

    cout << endl;

    //(c) Run NNI search on candidate trees
    size_t initialTreeCount = bestInitTrees.size();
    if (!progress_display::getProgressDisplay()) {
        cout << "Do NNI search on " << initialTreeCount << " best initial trees" << endl;
    }
    stop_rule.setCurIt(0);
    if (candidateTrees.size() > Params::getInstance().numSupportTrees) {
        candidateTrees.clear();// TODO why clear here???
    }
    candidateTrees.setMaxSize(Params::getInstance().numSupportTrees);
    vector<string>::iterator it;

    candidateCount = bestInitTrees.size();
    size_t candidate = 1;
    for (it = bestInitTrees.begin(); it != bestInitTrees.end(); ++it, ++candidate) {
        readTreeString(*it);
        stringstream treeDescription;
        treeDescription << "candidate tree " << candidate << " (of " << candidateCount << ")";
        string context = treeDescription.str();
        doNNISearch(true, context.c_str());
        string treeString = getTreeString();
        addTreeToCandidateSet(treeString, curScore, true, MPIHelper::getInstance().getProcessID());
        if (Params::getInstance().writeDistImdTrees) {
            intermediateTrees.update(treeString, curScore);
        }
    }

    // TODO turning this
    if (isMixlen()) {
        //(d) Optimizing model parameters for the best of the
        //    candidate trees.
        LOG_LINE(VB_MIN, "Optimizing model parameters for top "
            << min((int)candidateTrees.size(), params->popSize)
            << " candidate trees... ");
        
        startTime = getRealTime();
        bestInitTrees = candidateTrees.getBestTreeStrings(params->popSize);
        size_t countBestTrees = bestInitTrees.size();
        initProgress(countBestTrees, "Optimizing parameters", "optimized parameters for", "tree");
        for (it = bestInitTrees.begin(); it != bestInitTrees.end(); it++) {
            string tree;
            readTreeString(*it);
            //tree = optimizeBranches();
            tree = optimizeModelParameters();
//            cout << "Tree after brlen opt: " << tree << endl;
            hideProgress();
            cout << "Tree " << distance(bestInitTrees.begin(), it)+1 << " / LogL: " << getCurScore() << endl;
            candidateTrees.update(tree, getCurScore());
            showProgress();
            trackProgress(1.0);
        }
        doneProgress();
        LOG_LINE(VB_MIN, "... " << getRealTime() - startTime << " seconds");
    }
    //---- BLOCKING COMMUNICATION
    syncCandidateTrees(Params::getInstance().numSupportTrees, true);
}

string IQTree::generateParsimonyTree(int randomSeed) {
    string parsimonyTreeString;
    if (params->start_tree == STT_PLL_PARSIMONY) {
        pllInst->randomNumberSeed = randomSeed;
        pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInst,
                                                          pllPartitions, params->sprDist);
        resetBranches(pllInst);
        pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions,
                        pllInst->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE,
                        PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
        parsimonyTreeString = string(pllInst->tree_string);
        PhyloTree::readTreeString(parsimonyTreeString);
        wrapperFixNegativeBranch(true);
        parsimonyTreeString = getTreeString();
    } else if (params->start_tree == STT_RANDOM_TREE) {
        generateRandomTree(YULE_HARDING);
        wrapperFixNegativeBranch(true);
        parsimonyTreeString = getTreeString();
    } else {
        computeParsimonyTree(NULL, aln, randstream);
        parsimonyTreeString = getTreeString();
    }
    return parsimonyTreeString;
}

void IQTree::initializePLL(Params &params) {
    pllAttr.rateHetModel = PLL_GAMMA;
    pllAttr.fastScaling = PLL_FALSE;
    pllAttr.saveMemory = PLL_FALSE;
    pllAttr.useRecom = PLL_FALSE;
    pllAttr.randomNumberSeed = params.ran_seed;
    pllAttr.numberOfThreads = max(params.num_threads, 1); /* This only affects the pthreads version */
    if (pllInst != NULL) {
        pllDestroyInstance(pllInst);
    }
    /* Create a PLL getInstance */
    pllInst = pllCreateInstance(&pllAttr);

    /* Read in the alignment file */
    stringstream pllAln;
    aln->printAlignment(IN_PHYLIP, pllAln, "");
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
    if (!pllLoadAlignment(pllInst, pllAlignment, pllPartitions)) {
        outError("Incompatible tree/alignment combination");
    }
}

bool IQTree::isInitializedPLL() {
    return pllInst != nullptr;
}

void IQTree::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    try {
        if (!getModelFactory()) {
            if (isSuperTree()) {
                if (params.partition_type == BRLEN_OPTIMIZE || params.partition_type == TOPO_UNLINKED) {
                    setModelFactory(new PartitionModel(params, (PhyloSuperTree*) this, models_block));
                } else
                    setModelFactory(new PartitionModelPlen(params, (PhyloSuperTreePlen*) this, models_block));

                // mapTrees again in case of different rooting
                if (root)
                    ((PhyloSuperTree*)this)->mapTrees();

            } else {
                setModelFactory(new ModelFactory(params, model_name, this, models_block));
            }
        }
    } catch (string & str) {
        outError(str);
    }
    setModel(getModelFactory()->model);
    setRate(getModelFactory()->site_rate);
    getModelFactory()->setCheckpoint(checkpoint);
    /*
     * MDW: I don't understand why/how checkpointing is being used,
     * however I'm having problems because a restoreCheckpoint
     * happens before anything has been saved (in 
     * phyloanalysis.cpp:runTreeReconstruction). This fixes that
     * problem, but as I don't understand what is going on, I am
     * not at all confident this is the correct solution.
     */

     //!!! BQM: Because iqtree restore the checkpoint from the previous run! (not this current run)
     // That's why saveCheckpoint should NOT be called now! I comment this line out.

//    model->saveCheckpoint();

    if (params.pll) {
        if (getRate()->getNDiscreteRate() == 1) {
            outError("Non-Gamma model is not yet supported by PLL.");
            // TODO: change rateHetModel to PLL_CAT in case of non-Gamma model
        }
        if (getRate()->name.substr(0,2) == "+I")
            outError("+Invar model is not yet supported by PLL.");
        if (aln->seq_type == SEQ_DNA && getModel()->name != "GTR")
            outError("non GTR model for DNA is not yet supported by PLL.");
        pllInitModel(pllInst, pllPartitions);
    }

    if (aln->ordered_pattern.empty())
        aln->orderPatternByNumChars(PAT_VARIANT);

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
    PhyloNode *node = dad->getNeighborByIndex(nei_id)->getNode();
    int set_id = dad->id * 3 + nei_id;
    if (leaves_vec[set_id])
        return leaves_vec[set_id];
    RepresentLeafSet* leaves = new RepresentLeafSet;
    RepresentLeafSet* leaves_it[2] = { NULL, NULL };
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
        ASSERT(j == 2 && leaves_it[0] && leaves_it[1]);
        if ( 2 <= j && leaves_it[0]->empty() && leaves_it[1]->empty()) {
            cout << "wrong";
        }
        if (j < 2) {
            return leaves;
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
            ASSERT(id < 2 && id >= 0);
            leaves->insert(new RepLeaf((*lit[id])->leaf, (*lit[id])->height + 1));
            lit[id]++;
        }
    }
    ASSERT(!leaves->empty());
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
    PhyloNodeVector cherry_taxa;
    PhyloNodeVector noncherry_taxa;
    // get the vector of non cherry taxa
    getNonCherryLeaves(noncherry_taxa, cherry_taxa);
    root = NULL;
    size_t num_taxa = aln->getNSeq();
    int num_delete = k_delete;
    if (num_delete > num_taxa - 4)
        num_delete = num_taxa - 4;
    LOG_LINE(VB_DEBUG, "Deleting " << num_delete << " leaves");
    vector<unsigned int> indices_noncherry(noncherry_taxa.size());
    //iota(indices_noncherry.begin(), indices_noncherry.end(), 0);
    unsigned int startValue = 0;
    for (vector<unsigned int>::iterator it = indices_noncherry.begin(); it != indices_noncherry.end(); ++it) {
        (*it) = startValue;
        ++startValue;
    }
    my_random_shuffle(indices_noncherry.begin(), indices_noncherry.end());
    int i;
    for (i = 0; i < num_delete && i < noncherry_taxa.size(); i++) {
        PhyloNode *taxon = noncherry_taxa[indices_noncherry[i]];
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
        my_random_shuffle(indices_cherry.begin(), indices_cherry.end());
        while (i < num_delete) {
            PhyloNode *taxon = cherry_taxa[indices_cherry[j]];
            del_leaves.push_back(taxon);
            deleteLeaf(taxon);
            i++;
            j++;
        }
    }
    root = cherry_taxa[j];
}

void IQTree::deleteLeaves(PhyloNodeVector &del_leaves) {
    PhyloNodeVector taxa;
    // get the vector of taxa
    getTaxa(taxa);
    root = nullptr;
    //int num_delete = floor(p_delete * taxa.size());
    int num_delete = k_delete;
    int i;
    if (num_delete > taxa.size() - 4) {
        num_delete = taxa.size() - 4;
    }
    LOG_LINE(VB_DEBUG, "Deleting " << num_delete << " leaves");
    // now try to randomly delete some taxa of the probability of p_delete
    for (i = 0; i < num_delete;) {
        int id = random_int(taxa.size());
        if (!taxa[id])
            continue;
        else
            i++;
        PhyloNode *taxon = taxa[id];
        del_leaves.push_back(taxon);
        deleteLeaf(taxon);
        taxa[id] = nullptr;
    }
    // set root to the first taxon which was not deleted
    for (i = 0; i < taxa.size(); i++)
        if (taxa[i]) {
            root = taxa[i];
            break;
        }
}

int IQTree::assessQuartet(Node *leaf0, Node *leaf1, Node *leaf2, Node *del_leaf) {
    ASSERT(dist_matrix);
    size_t nseq = aln->getNSeq();
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
    if (!node) {
        node = getRoot();
    }
    if (dad!=nullptr) {
        PhyloNeighbor *node_nei = node->findNeighbor(dad);
        PhyloNeighbor *dad_nei  = dad->findNeighbor(node);
        node_nei->lh_scale_factor = 0.0;
        node_nei->clearComputedFlags();
        dad_nei->lh_scale_factor = 0.0;
        node_nei->clearComputedFlags();
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        initializeBonus(child, node);
    }
}

void IQTree::raiseBonus(PhyloNeighbor *nei, PhyloNode *dad, double bonus) {
    nei->lh_scale_factor += bonus;
    LOG_LINE(VB_DEBUG, dad->id << " - " << nei->node->id << " : " << bonus);
}

double IQTree::computePartialBonus(PhyloNode *node, PhyloNode* dad) {
    PhyloNeighbor *node_nei = node->findNeighbor(dad);
    if (node_nei->isLikelihoodComputed()) {
        return node_nei->lh_scale_factor;
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        node_nei->lh_scale_factor += computePartialBonus(child, node);
    }
    //Note: formerly this cleared the parsimony computed flag.
    //Now it only sets the likelihood computed flag (James B. 11-Sep-2020).
    node_nei->setLikelihoodComputed(true);
    return node_nei->lh_scale_factor;
}

void IQTree::findBestBonus(double &best_score, NodeVector &best_nodes, NodeVector &best_dads, PhyloNode *node, PhyloNode *dad) {
    double score;
    if (!node)
        node = getRoot();
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

    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        findBestBonus(best_score, best_nodes, best_dads, child, node);
    }
}

void IQTree::assessQuartets(vector<RepresentLeafSet*> &leaves_vec, PhyloNode *cur_root, PhyloNode *del_leaf) {
    const int MAX_DEGREE = 3;
    RepresentLeafSet * leaves[MAX_DEGREE];
    double bonus[MAX_DEGREE];
    memset(bonus, 0, MAX_DEGREE * sizeof(double));
    int cnt = 0;

    // only work for birfucating tree
    ASSERT(cur_root->degree() == MAX_DEGREE);

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
            raiseBonus(cur_root->getNeighborByIndex(cnt), cur_root, bonus[cnt]);

}

void IQTree::reinsertLeavesByParsimony(PhyloNodeVector &del_leaves) {
    ASSERT(0 && "this function is obsolete");
    ASSERT(root->isLeaf());
    for (auto it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); ++it_leaf) {
        //cout << "Add leaf " << (*it_leaf)->id << " to the tree" << endl;
        initializeAllPartialPars();
        clearAllPartialLH();
        Node *target_node = NULL;
        Node *target_dad = NULL;
        Node *added_node = (*it_leaf)->firstNeighbor()->node;
        Node *node1 = NULL;
        Node *node2 = NULL;
        //Node *leaf;
        for (int i = 0; i < 3; i++) {
            auto nei = added_node->neighbors[i];
            if (nei->node->id == (*it_leaf)->id) {
                //leaf = nei->node;
            } else if (!node1) {
                node1 = nei->node;
            } else {
                node2 = nei->node;
            }
        }

        //cout << "(" << node1->id << ", " << node2->id << ")" << "----" << "(" << added_node->id << "," << leaf->id << ")" << endl;
        added_node->updateNeighbor(node1, DUMMY_NODE_1);
        added_node->updateNeighbor(node2, DUMMY_NODE_2);

        best_pars_score = INT_MAX;
        // TODO: this needs to be adapted
//        addTaxonMPFast(added_node, target_node, target_dad, NULL, root->neighbors[0]->node, root);
        target_node->updateNeighbor(target_dad, added_node, -1.0);
        target_dad->updateNeighbor(target_node, added_node, -1.0);
        added_node->updateNeighbor(DUMMY_NODE_1, target_node, -1.0);
        added_node->updateNeighbor(DUMMY_NODE_2, target_dad, -1.0);
    }
}

void IQTree::reinsertLeaves(PhyloNodeVector &del_leaves) {

    //int num_del_leaves = del_leaves.size();
    ASSERT(root->isLeaf());

    for (auto it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); it_leaf++) {
        LOG_LINE(VB_DEBUG, "Reinserting " << (*it_leaf)->name << " (" << (*it_leaf)->id << ")");
        vector<RepresentLeafSet*> leaves_vec;
        leaves_vec.resize(nodeNum * 3, NULL);
        initializeBonus();
        PhyloNodeVector nodes;
        getInternalNodes(nodes);
        if (verbose_mode >= VB_DEBUG) {
            hideProgress();
            drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
            showProgress();
        }
        //printTree(cout, WT_BR_LEN | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
        for (auto it = nodes.begin(); it != nodes.end(); it++) {
            assessQuartets(leaves_vec, *it, *it_leaf);
        }
        NodeVector best_nodes, best_dads;
        double best_bonus;
        findBestBonus(best_bonus, best_nodes, best_dads);
        LOG_LINE(VB_DEBUG, "Best bonus " << best_bonus 
            << " " << best_nodes[0]->id << " " << best_dads[0]->id);
        ASSERT(best_nodes.size() == best_dads.size());
        int node_id = random_int(best_nodes.size());
        if (best_nodes.size() > 1) {
            LOG_LINE(VB_DEBUG, best_nodes.size() << " branches show the same best bonus"
                << ", branch nr. " << node_id << " is chosen");
        }
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

    if (verbose_mode >= VB_DEBUG) {
        hideProgress();
        drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
        showProgress();
    }
}

void IQTree::doParsimonyReinsertion() {
    PhyloNodeVector del_leaves;

    deleteLeaves(del_leaves);

    reinsertLeavesByParsimony(del_leaves);
    fixNegativeBranch(false);
}

void IQTree::getNonTabuBranches(Branches& allBranches, SplitGraph& tabuSplits, Branches& nonTabuBranches, Branches* tabuBranches) {
    if (tabuSplits.size() == 0) {
        return;
    }
    for (Branches::iterator it = allBranches.begin(); it != allBranches.end(); it++) {
        if (isInnerBranch(it->second.first, it->second.second)) {
            int nodeID1 = it->second.first->id;
            int nodeID2 = it->second.second->id;
            Branch curBranch = it->second;
            Split* sp = getSplit(it->second.first, it->second.second);
            if (!tabuSplits.containSplit(*sp)) {
                nonTabuBranches.insert(pair<int,Branch>(pairInteger(nodeID1, nodeID2), curBranch));
            } else {
                if (tabuBranches != NULL) {
                    tabuBranches->insert(pair<int,Branch>(pairInteger(nodeID1, nodeID2), curBranch));
                }
            }
            delete sp;
        }

    }
}

void IQTree::getSplitBranches(Branches &branches, SplitIntMap &splits, Node *node, Node *dad) {
    if (!node) {
        node = root;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
            if (isInnerBranch((*it)->node, node)) {
                Branch curBranch;
                curBranch.first = (*it)->node;
                curBranch.second = node;
                Split* curSplit;
                Split *sp = (*it)->split;
                ASSERT(sp != NULL);
                curSplit = new Split(*sp);
                if (curSplit->shouldInvert())
                    curSplit->invert();
                if (splits.findSplit(curSplit) != NULL) {
                    //curSplit->report(cout);
                    branches.insert(pair<int,Branch>(pairInteger(curBranch.first->id, curBranch.second->id), curBranch));
                }
                delete curSplit;
            }
            getSplitBranches(branches, splits, (*it)->node, node);
        }
}

bool IQTree::shouldEvaluate(Split *curSplit, SplitIntMap &tabuSplits, SplitIntMap &candSplits) {
    bool answer = true;
    /******************** CHECK TABU SPLIT **************************/
    if (tabuSplits.findSplit(curSplit) != NULL) {
        answer = false;
    } else if (!candSplits.empty()) {
        Split *_curSplit;
        /******************** CHECK STABLE SPLIT **************************/
        int value;
        _curSplit = candSplits.findSplit(curSplit, value);
        if (_curSplit == NULL || _curSplit->getWeight() <= params->stableSplitThreshold) {
            answer = true;
        } else { // add a stable branch with a certain probability
            double rndDbl = random_double();
            if (rndDbl > params->stableSplitThreshold) {
                answer = true;
            } else {
                answer = false;
            }
        }
    } else {
        answer = true;
    }
    return answer;
}


void IQTree::getNNIBranches(SplitIntMap &tabuSplits, SplitIntMap &candSplits,Branches &nonNNIBranches, Branches &nniBranches, Node *node, Node *dad) {
    if (!node) {
        node = root;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
            if (isInnerBranch((*it)->node, node)) {
                Branch curBranch;
                curBranch.first = (*it)->node;
                curBranch.second = node;
                int branchID = pairInteger(curBranch.first->id, curBranch.second->id);

                if (params->fixStableSplits) {
                    Split *curSplit;
                    Split *sp = (*it)->split;
                    ASSERT(sp != NULL);
                    curSplit = new Split(*sp);
                    if (curSplit->shouldInvert())
                        curSplit->invert();
                    if (shouldEvaluate(curSplit, tabuSplits, candSplits)) {
                        nniBranches.insert(pair<int, Branch>(branchID, curBranch));
                    } else {
                        nonNNIBranches.insert(pair<int, Branch>(branchID, curBranch));
                    }
                    delete curSplit;
                } else {
                    nniBranches.insert(pair<int, Branch>(branchID, curBranch));
                }
            }
            getNNIBranches(tabuSplits, candSplits, nonNNIBranches, nniBranches, (*it)->node, node);
        }
}

void IQTree::getStableBranches(SplitIntMap &candSplits, double supportValue, Branches &stableBranches, Node *node, Node *dad) {
    if (!node) {
        node = root;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
            if (isInnerBranch((*it)->node, node)) {
                Branch curBranch;
                curBranch.first = (*it)->node;
                curBranch.second = node;
                Split *curSplit;
                Split *sp = (*it)->split;
                ASSERT(sp != NULL);
                curSplit = new Split(*sp);
                if (curSplit->shouldInvert())
                    curSplit->invert();
                int occurences;
                sp = candSplits.findSplit(curSplit, occurences);
                if (sp != NULL) {
                    if ( sp->getWeight() >= supportValue) {
                        stableBranches.insert(
                                pair<int, Branch>(pairInteger(curBranch.first->id, curBranch.second->id), curBranch));
                    }
                }
                delete curSplit;
            }
            getStableBranches(candSplits, supportValue, stableBranches, (*it)->node, node);
        }
}

string IQTree::perturbStableSplits(double suppValue) {
    int numRandNNI = 0;
    Branches stableBranches;
//    initTabuSplits.clear();
//    stableBranches = getStableBranches(candidateTrees.getCandSplits(), suppValue);
//    int maxRandNNI = stableBranches.size() / 2;
    do {
        getStableBranches(candidateTrees.getCandSplits(), suppValue, stableBranches);
        vector<NNIMove> randomNNIs;
        vector<NNIMove> compatibleNNIs;
        for (map<int, Branch>::iterator it = stableBranches.begin(); it != stableBranches.end(); it++) {
            NNIMove randNNI = getRandomNNI(it->second);
            if (constraintTree.isCompatible(randNNI))
                randomNNIs.push_back(randNNI);
        }
        getCompatibleNNIs(randomNNIs, compatibleNNIs);
        for (vector<NNIMove>::iterator it = compatibleNNIs.begin(); it != compatibleNNIs.end(); it++) {
            doNNI(*it);
            numRandNNI++;
//            Split *sp = getSplit(it->node1, it->node2);
//            Split *tabuSplit = new Split(*sp);
//            if (tabuSplit->shouldInvert()) {
//                tabuSplit->invert();
//            }
//            initTabuSplits.insertSplit(tabuSplit, 1);
        }
    } while (stableBranches.size() > 0);

    if (verbose_mode >= VB_MAX) {
        cout << "Tree perturbation: number of random NNI performed = " << numRandNNI << endl;
    }
    setAlignment(aln);
    setRootNode(params->root);

    clearAllPartialLH();

    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    }
    if (params->pll) {
        pllReadNewick(getTreeString());
    }

    resetCurScore();
    return getTreeString();
}

string IQTree::doRandomNNIs(bool storeTabu) {
    int cntNNI = 0;
    int numRandomNNI;
    Branches nniBranches;
    Branches nonNNIBranches;
    if (storeTabu) {
        Branches stableBranches;
        getStableBranches(candidateTrees.getCandSplits(), Params::getInstance().stableSplitThreshold, stableBranches);
        int numNonStableBranches = leafNum - 3 - stableBranches.size();
        numRandomNNI = numNonStableBranches;
    } else {
        numRandomNNI = floor((leafNum - 3) * Params::getInstance().initPS);
        if (leafNum >= 4 && numRandomNNI == 0)
            numRandomNNI = 1;
    }

    initTabuSplits.clear();
    while (cntNNI < numRandomNNI) {
        nniBranches.clear();
        nonNNIBranches.clear();
        getNNIBranches(initTabuSplits, candidateTrees.getCandSplits(), nonNNIBranches, nniBranches);
        if (nniBranches.size() == 0) break;
        // Convert the map data structure Branches to vector of Branch
        vector<Branch> vectorNNIBranches;
        for (Branches::iterator it = nniBranches.begin(); it != nniBranches.end(); ++it) {
            vectorNNIBranches.push_back(it->second);
        }
        int randInt = random_int((int) vectorNNIBranches.size());
        NNIMove randNNI = getRandomNNI(vectorNNIBranches[randInt]);
        if (constraintTree.isCompatible(randNNI)) {
            // only if random NNI satisfies constraintTree
            doNNI(randNNI);
            if (storeTabu) {
                Split *sp = getSplit(randNNI.node1, randNNI.node2);
                Split *tabuSplit = new Split(*sp);
                if (tabuSplit->shouldInvert()) {
                    tabuSplit->invert();
                }
                initTabuSplits.insertSplit(tabuSplit, 1);
            }
        }
        cntNNI++;
    }
    if (verbose_mode >= VB_MAX)
        cout << "Tree perturbation: number of random NNI performed = " << cntNNI << endl;
    setAlignment(aln);
    setRootNode(params->root);

    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    }
    if (params->pll) {
        pllReadNewick(getTreeString());
    }

    clearAllPartialLH();
    resetCurScore();
    return getTreeString();
}


void IQTree::doIQP() {
    if (verbose_mode >= VB_DEBUG) {
        hideProgress();
        drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
        showProgress();
    }
    PhyloNodeVector del_leaves;
    deleteLeaves(del_leaves);
    reinsertLeaves(del_leaves);

    // just to make sure IQP does it right
    setAlignment(aln);
    if (params->pll) {
        pllReadNewick(getTreeString());
    }

    resetCurScore();

    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    }
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
    ASSERT(aln->num_states == 4);
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
    getRate()->setGammaShape(pllPartitions->partitionData[0]->alpha);
    if (aln->num_states == 4) {
        ((ModelMarkov*) getModel())->setRateMatrix(pllPartitions->partitionData[0]->substRates);
        getModel()->decomposeRateMatrix();
    }
    ((ModelMarkov*) getModel())->setStateFrequency(pllPartitions->partitionData[0]->empiricalFrequencies);
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
    ASSERT(node1->isLeaf());
    ASSERT(node2->isLeaf());

    PhyloNeighbor *node1nei = node1->firstNeighbor();
    PhyloNeighbor *node2nei = node2->firstNeighbor();

    node2nei->node->updateNeighbor(node2, node1);
    node1nei->node->updateNeighbor(node1, node2);

    // Update the new neightbors of the 2 nodes
    node1->updateNeighbor(node1->neighbors.begin(), node2nei);
    node2->updateNeighbor(node2->neighbors.begin(), node1nei);

    PhyloNeighbor *node1NewNei = node1->firstNeighbor();
    PhyloNeighbor *node2NewNei = node2->firstNeighbor();

    // Reoptimize the branch lengths
    optimizeOneBranch(node1, node1NewNei->getNode());
    //this->curScore = optimizeOneBranch(node2, node2NewNei->getNode());
    optimizeOneBranch(node2, node2NewNei->getNode());
    //drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    this->curScore = computeLikelihoodFromBuffer();
    return this->curScore;
}

double IQTree::perturb(int times) {
    while (times > 0) {
        PhyloNodeVector taxa;
        // get the vector of taxa
        getTaxa(taxa);
        int taxonid1 = random_int(taxa.size());
        PhyloNode *taxon1 = taxa[taxonid1];
        PhyloNode *taxon2;
        int *dists = new int[taxa.size()];
        int minDist = 1000000;
        for (int i = 0; i < taxa.size(); i++) {
            if (i == taxonid1) {
                continue;
            }
            taxon2 = taxa[i];
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

        taxon2 =  taxa[taxonid2];

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

//extern "C" pllUFBootData * pllUFBootDataPtr;
extern pllUFBootData * pllUFBootDataPtr;

string IQTree::optimizeModelParameters(bool printInfo, double logl_epsilon) {
    prepareToComputeDistances();
    if (logl_epsilon == -1) {
        logl_epsilon = params->modelEps;
    }
    cout << "Estimate model parameters (epsilon = " << logl_epsilon << ")" << endl;
    double stime = getRealTime();
    string newTree;
    if (params->pll) {
        if (curScore == -DBL_MAX) {
            pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_TRUE, PLL_FALSE);
        } else {
            pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_FALSE, PLL_FALSE);
        }
        pllOptimizeModelParameters(pllInst, pllPartitions, logl_epsilon);
        curScore = pllInst->likelihood;
        pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions,
                pllInst->start->back, PLL_TRUE,
                PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH,
                PLL_FALSE, PLL_FALSE);
        if (printInfo) {
            pllPrintModelParams();
        }
        newTree = string(pllInst->tree_string);
        double etime = getRealTime();
        if (printInfo)
            cout << etime - stime << " seconds (logl: " << curScore << ")" << endl;
    } else {
        double modOptScore;
        if (params->opt_gammai) { // DO RESTART ON ALPHA AND P_INVAR
            modOptScore = getModelFactory()->optimizeParametersGammaInvar(params->fixed_branch_length, printInfo, logl_epsilon);
            params->opt_gammai = false;
        } else {
            modOptScore = getModelFactory()->optimizeParameters(params->fixed_branch_length, printInfo, logl_epsilon);
        }
        if (isSuperTree()) {
            ((PhyloSuperTree*) this)->computeBranchLengths();
        }
        if (getModelFactory()->isUnstableParameters() && aln->seq_type != SEQ_CODON) {
            cout << endl;
            outWarning("Estimated model parameters are at boundary that can cause numerical instability!");
            cout << endl;
        }
        if (modOptScore < curScore - 1.0 && params->partition_type != TOPO_UNLINKED) {
            cout << "  BUG: Tree logl gets worse after model optimization!" << endl;
            cout << "  Old logl: " << curScore << " / " << "new logl: " << modOptScore << endl;
            printTree("debug.tree");
            abort();
        } else {
            curScore = modOptScore;
            newTree = getTreeString();
        }
        if (params->print_trees_site_posterior)
            computePatternCategories();
    }
    doneComputingDistances(); //DISABLED
    return newTree;
}

string IQTree::ensureModelParametersAreSet(double initEpsilon) {
    string initTree;
    getModelFactory()->restoreCheckpoint();
    if (getCheckpoint()->getBool("finishedModelInit")) {
        // model optimization already done: ignore this step
        if (!candidateTrees.empty()) {
            readTreeString(getBestTrees()[0]);
        }
        setCurScore(computeLikelihood());
        initTree = getTreeString();
        cout << "CHECKPOINT: Model parameters restored, LogL: " << getCurScore() << endl;
    } else {
        initTree = optimizeModelParameters(true, initEpsilon);
        if (isMixlen()) {
            initTree = ((ModelFactoryMixlen*)getModelFactory())->sortClassesByTreeLength();
        }
        saveCheckpoint();
        getModelFactory()->saveCheckpoint();
        getCheckpoint()->putBool("finishedModelInit", true);
        getCheckpoint()->dump();
    }
    return initTree;
}

void IQTree::printBestScores() {
    vector<double> bestScores = candidateTrees.getBestScores(params->popSize);
    for (vector<double>::iterator it = bestScores.begin(); it != bestScores.end(); it++)
        cout << (*it) << " ";
    cout << endl;
}

double IQTree::computeLogL() {
    if (params->pll) {
        if (curScore == -DBL_MAX) {
            pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_TRUE, PLL_FALSE);
        } else {
            pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_FALSE, PLL_FALSE);
        }
        curScore = pllInst->likelihood;
    } else {
//        if (!lhComputed) {
//            initializeAllPartialLh();
//            clearAllPartialLH();
//        }
        curScore = computeLikelihood();
    }
    return curScore;
}

string IQTree::optimizeBranches(int maxTraversal) {
    string tree;
    if (params->pll) {
        if (curScore == -DBL_MAX) {
            pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_TRUE, PLL_FALSE);
//            lhComputed = true;
        }
        pllOptimizeBranchLengths(pllInst, pllPartitions, maxTraversal);
        curScore = pllInst->likelihood;
        pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE,
                PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
        tree = string(pllInst->tree_string);
    } else {
//        if (!lhComputed) {
//            initializeAllPartialLh();
//            clearAllPartialLH();
//            lhComputed = true;
//        }
        curScore = optimizeAllBranches(maxTraversal, params->loglh_epsilon, PLL_NEWZPERCYCLE);
        tree = getTreeString();
    }
    return tree;
}

double IQTree::doTreeSearch() {
    
    if (params->numInitTrees > 1) {
        cout << "--------------------------------------------------------------------" << endl;
        cout << "|             INITIALIZING CANDIDATE TREE SET                      |" << endl;
        cout << "--------------------------------------------------------------------" << endl;
    }

    double initCPUTime = getRealTime();
    int treesPerProc = (params->numInitTrees) / MPIHelper::getInstance().getNumProcesses() - candidateTrees.size();
    if (params->numInitTrees % MPIHelper::getInstance().getNumProcesses() != 0) {
        treesPerProc++;
    }
    if (treesPerProc < 0) {
        treesPerProc = 0;
    }
    // Master node does one tree less because it already created the BIONJ tree
//    if (MPIHelper::getInstance().isMaster()) {
//        treesPerProc--;
//    }

    // Make sure to get at least 1 tree
    if (treesPerProc < 1 && params->numInitTrees > candidateTrees.size()) {
        treesPerProc = 1;
    }

    /* Initialize candidate tree set */
    if (!getCheckpoint()->getBool("finishedCandidateSet")) {
        initCandidateTreeSet(treesPerProc, params->numNNITrees);
        // write best tree to disk
        printBestCandidateTree();
        saveCheckpoint();
        getCheckpoint()->putBool("finishedCandidateSet", true);
        getCheckpoint()->dump();
    } else {
        cout << "CHECKPOINT: Candidate tree set restored, best LogL: " << candidateTrees.getBestScore() << endl;
    }
    ASSERT(candidateTrees.size() != 0);
    cout << "Finish initializing candidate tree set (" << candidateTrees.size() << ")" << endl;


    cout << "Current best tree score: " << candidateTrees.getBestScore() << " / CPU time: " <<
    getRealTime() - initCPUTime << endl;
    cout << "Number of iterations: " << stop_rule.getCurIt() << endl;

//    string treels_name = params->out_prefix;
//    treels_name += ".treels";
//    string out_lh_file = params->out_prefix;
//    out_lh_file += ".treelh";
//    string site_lh_file = params->out_prefix;
//    site_lh_file += ".sitelh";
//
//    if (params->print_tree_lh) {
//        out_treelh.open(out_lh_file.c_str());
//        out_sitelh.open(site_lh_file.c_str());
//    }

//    if (params->write_intermediate_trees)
//        out_treels.open(treels_name.c_str());

//    if (params->write_intermediate_trees && save_all_trees != 2) {
//        printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
//    }

    setRootNode(params->root);

    if (!getCheckpoint()->getBool("finishedCandidateSet"))
        cout << "CHECKPOINT: " << stop_rule.getCurIt() << " search iterations restored" << endl;

    searchinfo.curPerStrength = params->initPS;
    double cur_correlation = 0.0;


    if ((Params::getInstance().fixStableSplits || Params::getInstance().adaptPertubation) && candidateTrees.size() > 1) {
        candidateTrees.computeSplitOccurences(Params::getInstance().stableSplitThreshold);
    }

    // tracking of worker candidate set is changed from master candidate set
    candidateset_changed.resize(MPIHelper::getInstance().getNumProcesses(), false);
    bestcandidate_changed = false;

    /*==============================================================================================================
                                           MAIN LOOP OF THE IQ-TREE ALGORITHM
     *=============================================================================================================*/

    bool early_stop = stop_rule.meetStopCondition(stop_rule.getCurIt(), cur_correlation);
    if (!early_stop) {
        cout << "--------------------------------------------------------------------" << endl;
        cout << "|               OPTIMIZING CANDIDATE TREE SET                      |" << endl;
        cout << "--------------------------------------------------------------------" << endl;
    }

    // count threshold for computing bootstrap correlation
    int ufboot_count, ufboot_count_check;
    stop_rule.getUFBootCountCheck(ufboot_count, ufboot_count_check);
    showNoProgress();
    
    while (!stop_rule.meetStopCondition(stop_rule.getCurIt(), cur_correlation)) {
        searchinfo.curIter = stop_rule.getCurIt();
        // estimate logl_cutoff for bootstrap
        if (!boot_orig_logl.empty())
            logl_cutoff = *min_element(boot_orig_logl.begin(), boot_orig_logl.end());

        if (estimate_nni_cutoff && nni_info.size() >= 500) {
            estimate_nni_cutoff = false;
            estimateNNICutoff(params);
        }

        Alignment *saved_aln = aln;

        string curTree;
        /*----------------------------------------
         * Perturb the tree
         *---------------------------------------*/
        doTreePerturbation();

        /*----------------------------------------
         * Optimize tree with NNI
         *----------------------------------------*/
        pair<int, int> nniInfos; // <num_NNIs, num_steps>
        nniInfos = doNNISearch(true, "");
        curTree = getTreeString();
        int pos = addTreeToCandidateSet(curTree, curScore, true, MPIHelper::getInstance().getProcessID());
        if (pos != -2 && pos != -1 && (Params::getInstance().fixStableSplits || Params::getInstance().adaptPertubation)) {
            candidateTrees.computeSplitOccurences(Params::getInstance().stableSplitThreshold);
        }
        if (MPIHelper::getInstance().isWorker() || MPIHelper::getInstance().gotMessage()) {
            syncCurrentTree();
        }
        // TODO: cannot check yet, need to somehow return treechanged
//        if (nni_count == 0 && params->snni && numPerturb > 0 && treechanged) {
//            assert(0 && "BUG: NNI could not improved perturbed tree");
//        }
        if (iqp_assess_quartet == IQP_BOOTSTRAP) {
            // restore alignment
            delete aln;
            setAlignment(saved_aln);
            initializeAllPartialLh();
            clearAllPartialLH();
        }
        if (isSuperTree()) {
            ((PhyloSuperTree *) this)->computeBranchLengths();
        }
        /*----------------------------------------
         * Print information
         *---------------------------------------*/
        //printInterationInfo();

//        if (params->write_intermediate_trees && save_all_trees != 2) {
//            printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
//        }

        if (params->snni && verbose_mode >= VB_DEBUG) {
            hideProgress();
            printBestScores();
            showProgress();
        }
        // DTH: make pllUFBootData usable in summarizeBootstrap
        if (params->pll && params->online_bootstrap && (params->gbo_replicates > 0))
            pllConvertUFBootData2IQTree();
        // DTH: Carefully watch the -pll case here

        /*----------------------------------------
         * convergence criterion for ultrafast bootstrap
         *---------------------------------------*/
         
        // MASTER receives bootstrap trees and perform stop convergence test
        if ((stop_rule.getCurIt()) >= ufboot_count &&
            params->stop_condition == SC_BOOTSTRAP_CORRELATION && MPIHelper::getInstance().isMaster()) {
            ufboot_count += params->step_iterations/2;
            // compute split support every half step
            SplitGraph *sg = new SplitGraph;
            summarizeBootstrap(*sg);
            sg->removeTrivialSplits();
            sg->setCheckpoint(checkpoint);
            boot_splits.push_back(sg);
            cout << "Log-likelihood cutoff on original alignment: " << logl_cutoff << endl;
//            MPIHelper::getInstance().sendMsg(LOGL_CUTOFF_TAG, convertDoubleToString(logl_cutoff));

            // check convergence every full step
            if (stop_rule.getCurIt() >= ufboot_count_check) {
                ufboot_count_check += params->step_iterations;
                cur_correlation = computeBootstrapCorrelation();
                cout << "NOTE: Bootstrap correlation coefficient of split occurrence frequencies: " <<
                cur_correlation << endl;
                if (!stop_rule.meetCorrelation(cur_correlation)) {
                    cout << "NOTE: UFBoot does not converge, continue at least " << params->step_iterations << " more iterations" << endl;
                }
            }
            if (params->gbo_replicates && params->online_bootstrap && params->print_ufboot_trees) {
                writeUFBootTrees(*params);
            }
        } // end of bootstrap convergence test

        // print UFBoot trees every 10 iterations

        saveCheckpoint();
        checkpoint->dump();

        if (bestcandidate_changed) {
            printBestCandidateTree();
            bestcandidate_changed = false;
        }
        //if (params->partition_type)
        //     ((PhyloSuperTreePlen*)this)->printNNIcasesNUM();
    }

    // 2019-06-03: check convergence here to avoid effect of refineBootTrees
    if (boot_splits.size() >= 2 && MPIHelper::getInstance().isMaster()) {
        // check the stopping criterion for ultra-fast bootstrap
        if (computeBootstrapCorrelation() < params->min_correlation) {
            cout << "WARNING: bootstrap analysis did not converge."
                << " You should rerun with higher number of iterations (-nm option)" << endl;
        }
    }
    if(params->ufboot2corr) {
        refineBootTrees();
    }
    if (!early_stop) {
        sendStopMessage();
    }
    readTreeString(candidateTrees.getBestTreeStrings()[0]);

    if (testNNI) {
        outNNI.close();
    }
    if (params->write_intermediate_trees) {
        out_treels.close();
    }
    if (params->print_tree_lh) {
        out_treelh.close();
        out_sitelh.close();
    }
    // DTH: pllUFBoot deallocation
    if (params->pll) {
        pllDestroyUFBootData();
    }

#ifdef _IQTREE_MPI
    cout << "Total number of trees received: " << MPIHelper::getInstance().getNumTreeReceived() << endl;
    cout << "Total number of trees sent: " << MPIHelper::getInstance().getNumTreeSent() << endl;
    cout << "Total number of NNI searches done by myself: " << MPIHelper::getInstance().getNumNNISearch() << endl;
    MPIHelper::getInstance().resetNumbers();
#endif

    cout << "TREE SEARCH COMPLETED AFTER " << stop_rule.getCurIt() << " ITERATIONS"
    << " / Time: " << convert_time(getRealTime() - params->start_real_time) << endl << endl;

    return candidateTrees.getBestScore();

}


/*
void IQTree::refineBootTrees(){
    on_refine_btree = true;
    int num_boot_rep = params->gbo_replicates;
    params->gbo_replicates = 0;
    save_all_trees = 0;

    NNI_Type saved_nni_type = params->nni_type;
    if(params->u2c_nni5 == false){
        params->nni5 = false;
        params->nni_type = NNI1;
    }else{
        params->nni5 = true;
        params->nni_type = NNI5;
    }

    string saved_tree = getTreeString();
    saved_aln_on_refine_btree = aln;

    int nptn = getAlnNPattern();
    cout << "npn = " << nptn << endl;

    string tree;
    Alignment * bootstrap_aln;

    int *saved_randstream = randstream;
    init_random(params->ran_seed);

    string file_name = params->out_prefix;
    file_name += ".binfo";
    ofstream binfo(file_name.c_str());

    file_name = params->out_prefix;
    file_name += ".btree.pre";
    ofstream btreep(file_name.c_str(), ios::app);

    file_name = params->out_prefix;
    file_name += ".btree.aft";
    ofstream btreea(file_name.c_str(), ios::app);

    file_name = params->out_prefix;
    file_name += ".btree.pre.brlen";
    ofstream btreep_brlen(file_name.c_str(), ios::app);

    file_name = params->out_prefix;
    file_name += ".btree.aft.brlen";
    ofstream btreea_brlen(file_name.c_str(), ios::app);

    file_name = params->out_prefix;
    file_name += ".refine.aln";
    ofstream refine_aln(file_name.c_str(), ios::app);


    if(!isSuperTree()){
        getModelFactory()->model->writeInfo(binfo);
        getModelFactory()->site_rate->writeInfo(binfo);
    }else{
        ((PhyloSuperTree *)this)->at(0)->getModelFactory()->model->writeInfo(binfo);
        ((PhyloSuperTree *)this)->at(0)->getModelFactory()->site_rate->writeInfo(binfo);
    }

    PhyloSuperTree *stree = (isSuperTree()) ? (PhyloSuperTree*)this : NULL;

    for(int sample = 0; sample < num_boot_rep; sample++){
        if ((sample+1) % 100 == 0)
            cout << sample+1 << " replicates done" << endl;

        if (saved_aln_on_refine_btree->isSuperAlignment())
            bootstrap_aln = new SuperAlignment;
        else
            bootstrap_aln = new Alignment;

//        bootstrap_aln->buildFromPatternFreq(*saved_aln_on_refine_btree, boot_samples_int[sample]);
        IntVector this_sample;
        bootstrap_aln->createBootstrapAlignment(saved_aln_on_refine_btree, &this_sample, params->bootstrap_spec);

        if (isSuperTree())
            stree->setSuperAlignment(bootstrap_aln);
        else
            setAlignment(bootstrap_aln);

        for(int k = 0; k < nptn; k++){
            if(boot_samples_int[sample][k] != this_sample[k]){
                assert(0 && "Not identical");
            }
        }

//        bootstrap_aln->createBootstrapAlignment(saved_aln_on_refine_btree, NULL, params->bootstrap_spec);
        ptn_freq_computed = false;
        computePtnFreq();
        if (isSuperTree()) {
            for (auto it = stree->begin(); it != stree->end(); it++) {
                (*it)->ptn_freq_computed = false;
                (*it)->computePtnFreq();
            }
        }


        tree = boot_trees[sample];
        // Read the bootstrap tree
//        cout << tree << endl;

        readTreeString(tree);


        if(!isSuperTree())
            aln->printPhylip(refine_aln, true);
        else
            ((SuperAlignment *) aln)->printCombinedAlignment(refine_aln, true);

        setRootNode(params->root);

        if (isSuperTree()){
            wrapperFixNegativeBranch(true);
            ((PhyloSuperTree*) this)->mapTrees();
        }
        else{
            initializeAllPartialLh();
            clearAllPartialLH();
            wrapperFixNegativeBranch(false);
        }


//        printTree(cout);
//        cout << endl;
        optimizeBranches();
//        printTree(cout);
//        cout << endl;

        curScore = computeLogL();
        binfo << sample << "\t" << curScore;

        printTree(btreep, WT_NEWLINE | WT_SORT_TAXA);
        printTree(btreep_brlen, WT_NEWLINE | WT_SORT_TAXA | WT_BR_LEN);

        pair<int, int> nniInfos; // <num_NNIs, num_steps>
        nniInfos = doNNISearch(true, "");
        curScore = computeLogL();
        binfo << "\t" << curScore << endl;

        stringstream ostr;
        printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
        tree = ostr.str();
        boot_trees[sample] = getTreeString();
        boot_logl[sample] = curScore;

        printTree(btreea, WT_NEWLINE | WT_SORT_TAXA);
        printTree(btreea_brlen, WT_NEWLINE | WT_SORT_TAXA | WT_BR_LEN);
        delete aln;
    }

    binfo.close();
    btreep.close();
    btreea.close();

    // restore randstream
    finish_random();
    randstream = saved_randstream;

    // Recover the last status of IQTREE

    params->gbo_replicates = num_boot_rep;

    if (isSuperTree())
        stree->setSuperAlignment(saved_aln_on_refine_btree);
    else
        setAlignment(saved_aln_on_refine_btree);

    readTreeString(saved_tree);

    ptn_freq_computed = false;
    computePtnFreq();
    if (isSuperTree()) {
        for (auto it = stree->begin(); it != stree->end(); it++) {
            (*it)->ptn_freq_computed = false;
            (*it)->computePtnFreq();
        }
    }

    if(params->u2c_nni5 == false){
        params->nni_type = saved_nni_type;
        if(params->nni_type == NNI5)
            params->nni5 = true;
    }

    initializeAllPartialLh();
    clearAllPartialLH();
    curScore = optimizeAllBranches();

    save_all_trees = 2;
    on_refine_btree = false;
}
*/

/**********************************************************
 * STANDARD NON-PARAMETRIC BOOTSTRAP
 ***********************************************************/
void IQTree::refineBootTrees() {

    int *saved_randstream = randstream;
    init_random(params->ran_seed);

    params->gbo_replicates = 0;

    NNI_Type saved_nni_type = params->nni_type;
    // TODO: A bug in PhyloSuperTreePlen::swapNNIBranch by nni1
    // Thus always turn on -nni5 by PhyloSuperTreePlen
    if(params->u2c_nni5 == false && (!isSuperTree() || params->partition_type == BRLEN_OPTIMIZE)) {
        params->nni5 = false;
        params->nni_type = NNI1;
    }else{
        params->nni5 = true;
        params->nni_type = NNI5;
    }

    cout << "Refining ufboot trees with NNI ";
    if (params->nni5)
        cout << "5 branches..." << endl;
    else
        cout << "1 branch..." << endl;

    int refined_trees = 0;

    int refined_samples = 0;

    checkpoint->startStruct("UFBoot");
    if (CKP_RESTORE(refined_samples))
        cout << "CHECKPOINT: " << refined_samples << " refined samples restored" << endl;
    checkpoint->endStruct();
    
    // 2018-08-17: delete duplicated memory
    deleteAllPartialLh();

    ModelsBlock *models_block = readModelsDefinition(*params);
    
	// do bootstrap analysis
	for (int sample = refined_samples; sample < boot_trees.size(); sample++) {
        // create bootstrap alignment
        Alignment* bootstrap_alignment;
        if (aln->isSuperAlignment())
            bootstrap_alignment = new SuperAlignment;
        else
            bootstrap_alignment = new Alignment;
        bootstrap_alignment->createBootstrapAlignment(aln, NULL, params->bootstrap_spec);

        // create bootstrap tree
        IQTree *boot_tree;
        if (aln->isSuperAlignment()){
            if(params->partition_type != BRLEN_OPTIMIZE){
                boot_tree = new PhyloSuperTreePlen((SuperAlignment*) bootstrap_alignment, (PhyloSuperTree*) this);
            } else {
                boot_tree = new PhyloSuperTree((SuperAlignment*) bootstrap_alignment, (PhyloSuperTree*) this);
            }
        } else {
            // allocate heterotachy tree if neccessary
            int pos = posRateHeterotachy(aln->model_name);
            
            if (params->num_mixlen > 1) {
                boot_tree = new PhyloTreeMixlen(bootstrap_alignment, params->num_mixlen);
            } else if (pos != string::npos) {
                boot_tree = new PhyloTreeMixlen(bootstrap_alignment, 0);
            } else
                boot_tree = new IQTree(bootstrap_alignment);
        }

        boot_tree->on_refine_btree = true;
        boot_tree->save_all_trees = 0;

        // initialize constraint tree
        if (!constraintTree.empty()) {
            boot_tree->constraintTree.readConstraint(constraintTree);
        }

        boot_tree->setParams(params);

        // 2019-06-03: bug fix setting part_info properly
        if (boot_tree->isSuperTree())
            ((PhyloSuperTree*)boot_tree)->setPartInfo((PhyloSuperTree*)this);

        // copy model
        // BQM 2019-05-31: bug fix with -bsam option
        boot_tree->initializeModel(*params, aln->model_name, models_block);
        boot_tree->getModelFactory()->setCheckpoint(getCheckpoint());
        if (isSuperTree())
            ((PartitionModel*)boot_tree->getModelFactory())->PartitionModel::restoreCheckpoint();
        else
            boot_tree->getModelFactory()->restoreCheckpoint();

        // set likelihood kernel
        boot_tree->setParams(params);
        boot_tree->setLikelihoodKernel(sse);
        boot_tree->setNumThreads(num_threads);

        // load the current ufboot tree
        // 2019-02-06: fix crash with -sp and -bnni
        if (isSuperTree())
            boot_tree->PhyloTree::readTreeString(boot_trees[sample]);
        else
            boot_tree->readTreeString(boot_trees[sample]);
        
        if (boot_tree->isSuperTree() && params->partition_type == BRLEN_OPTIMIZE) {
            if (((PhyloSuperTree*)boot_tree)->size() > 1) {
                // re-initialize branch lengths for unlinked model
                boot_tree->wrapperFixNegativeBranch(true);
            }
        }
        
        // TODO: check if this resolves the crash in reorientPartialLh()
        boot_tree->initializeAllPartialLh();

        // just in case some branch lengths are negative
        if (int num_neg = boot_tree->wrapperFixNegativeBranch(false))
            outWarning("Bootstrap tree " + convertIntToString(sample+1) + " has " +
                convertIntToString(num_neg) + "non-positive branch lengths");

        // REMARK: branch lengths were estimated from original alignments
        // for bootstrap_alignment, they still thus need to be reoptimized a bit
        boot_tree->optimizeBranches(2);

        stringstream sampleDescription;
        sampleDescription << "bootstrap tree " << sample
            << "( of " << boot_trees.size() << ")";
        string context = sampleDescription.str();
        auto num_nnis = boot_tree->doNNISearch(true, context.c_str());
        if (num_nnis.second != 0) {
            refined_trees++;
        }
        if (verbose_mode >= VB_MED) {
            cout << "UFBoot tree " << sample+1 << ": " << boot_logl[sample] << " -> " << boot_tree->getCurScore() << endl;
        }

        stringstream ostr;
        if (params->print_ufboot_trees == 2)
            boot_tree->printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_LEN_SHORT);
        else
            boot_tree->printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
        boot_trees[sample] = ostr.str();
        boot_logl[sample] = boot_tree->curScore;


        // delete memory
        //boot_tree->setModelFactory(NULL);
        boot_tree->save_all_trees = 2;

        bootstrap_alignment = boot_tree->aln;
        delete boot_tree;
        // fix bug: bootstrap_alignment might be changed
        delete bootstrap_alignment;


        if ((sample+1) % 100 == 0)
            cout << sample+1 << " samples done" << endl;

        saveCheckpoint();
        checkpoint->startStruct("UFBoot");
        refined_samples = sample;
        CKP_SAVE(refined_samples);
        checkpoint->endStruct();

        checkpoint->dump();

	}
    
    delete models_block;

    cout << "Total " << refined_trees << " ufboot trees refined" << endl;

    // restore randstream
    finish_random();
    randstream = saved_randstream;

    SplitGraph *sg = new SplitGraph;
    summarizeBootstrap(*sg);
    sg->removeTrivialSplits();
    sg->setCheckpoint(checkpoint);
    boot_splits.push_back(sg);

    saveCheckpoint();
    checkpoint->dump();
    
    // restore
    params->gbo_replicates = boot_trees.size();
    params->nni_type = saved_nni_type;
    if(params->nni_type == NNI5) {
        params->nni5 = true;
    } else
        params->nni5 = false;

    // 2018-08-17: recover memory
    initializeAllPartialLh();

}


void IQTree::printIterationInfo(int sourceProcID) {
    double realtime_remaining = stop_rule.getRemainingTime(stop_rule.getCurIt());
    cout.setf(ios_base::fixed, ios_base::floatfield);

    // only print every 10th iteration or high verbose mode
    if (stop_rule.getCurIt() % 10 == 0 || verbose_mode >= VB_MED) {
        hideProgress();
        cout << ((iqp_assess_quartet == IQP_BOOTSTRAP) ? "Bootstrap " : "Iteration ") << stop_rule.getCurIt() <<
        " / LogL: ";
        cout << curScore;
        cout << " / Time: " << convert_time(getRealTime() - params->start_real_time);
        
        if (stop_rule.getCurIt() > 20) {
            cout << " (" << convert_time(realtime_remaining) << " left)";
        }
        if (MPIHelper::getInstance().getNumProcesses() > 1)
            cout << " / Process: " << sourceProcID;
        cout << endl;
        showProgress();
    }
}

//void IQTree::estimateLoglCutoffBS() {
//    if (params->avoid_duplicated_trees && max_candidate_trees > 0 && treels_logl.size() > 1000) {
//        int predicted_iteration;
//        predicted_iteration = ((stop_rule.getCurIt() + params->step_iterations - 1) / params->step_iterations)
//                              * params->step_iterations;
//        int num_entries = (int) floor(max_candidate_trees * ((double) stop_rule.getCurIt() / predicted_iteration));
//        if (num_entries < treels_logl.size() * 0.9) {
//            DoubleVector logl = treels_logl;
//            nth_element(logl.begin(), logl.begin() + (treels_logl.size() - num_entries), logl.end());
//            logl_cutoff = logl[treels_logl.size() - num_entries] - 1.0;
//        } else
//            logl_cutoff = 0.0;
//        if (verbose_mode >= VB_MED) {
//            if (stop_rule.getCurIt() % 10 == 0) {
//                cout << treels.size() << " trees, " << treels_logl.size() << " logls, logl_cutoff= " << logl_cutoff;
//                if (params->store_candidate_trees)
//                    cout << " duplicates= " << duplication_counter << " ("
//                    << (int) round(100 * ((double) duplication_counter / treels_logl.size())) << "%)" << endl;
//                else
//                    cout << endl;
//            }
//        }
//    }
//}


double IQTree::doTreePerturbation() {
    if (iqp_assess_quartet == IQP_BOOTSTRAP) {
        // create bootstrap sample
        Alignment *bootstrap_alignment;
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
            if (Params::getInstance().five_plus_five) {
                readTreeString(candidateTrees.getNextCandTree());
            } else {
                readTreeString(candidateTrees.getRandTopTree(Params::getInstance().popSize));
            }
            if (Params::getInstance().iqp) {
                doIQP();
            } else if (Params::getInstance().adaptPertubation) {
                perturbStableSplits(Params::getInstance().stableSplitThreshold);
            } else {
                doRandomNNIs(Params::getInstance().tabu);
            }
        } else {
            // Using the IQPNNI algorithm (best tree is selected)
            readTreeString(getBestTrees()[0]);
            doIQP();
        }
        if (params->count_trees) {
            string perturb_tree_topo = getTopologyString(false);
            if (pllTreeCounter.find(perturb_tree_topo) == pllTreeCounter.end()) {
                // not found in hash_map
                pllTreeCounter[perturb_tree_topo] = 1;
            } else {
                // found in hash_map
                pllTreeCounter[perturb_tree_topo]++;
            }
        }
        //optimizeBranches(1);
        curScore = computeLogL();
    }
    return curScore;
}

/****************************************************************************
 Fast Nearest Neighbor Interchange by maximum likelihood
 ****************************************************************************/
pair<int, int> IQTree::doNNISearch(bool write_info, const char* context) {

    computeLogL();
    double curBestScore = getBestScore();

    if (Params::getInstance().write_intermediate_trees && save_all_trees != 2) {
        printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
    }

    pair<int, int> nniInfos; // Number of NNIs and number of steps
    if (params->pll) {
        if (params->partition_file)
            outError("Unsupported -pll -sp combination!");
        curScore = pllOptimizeNNI(nniInfos.first, nniInfos.second, searchinfo);
        pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_TRUE,
                PLL_TRUE, 0, 0, 0, PLL_SUMMARIZE_LH, 0, 0);
        readTreeString(string(pllInst->tree_string));
    } else {
        nniInfos = optimizeNNI(Params::getInstance().speednni, context);
        if (isSuperTree()) {
            ((PhyloSuperTree*) this)->computeBranchLengths();
        }
        if (params->print_trees_site_posterior)
            computePatternCategories();
    }

    if(!on_refine_btree){ // Diep add (IF in Refinement Step, do not optimize model parameters)
        // Better tree or score is found
        if (getCurScore() > curBestScore + params->modelEps) {
            // Re-optimize model parameters (the sNNI algorithm)
            optimizeModelParameters(write_info, params->modelEps * 10);
            getModelFactory()->saveCheckpoint();

            // 2018-01-09: additional optimize root position
            // TODO: does not work with SuperTree yet
            if (rooted && !isSuperTree() && params->root_move_dist > 0)
            {
                optimizeRootPosition(params->root_move_dist, true, params->modelEps * 10);
            }
        }
    }
    MPIHelper::getInstance().setNumNNISearch(MPIHelper::getInstance().getNumNNISearch() + 1);

    return nniInfos;
}

pair<int, int> IQTree::optimizeNNI(bool speedNNI, const char* context) {
    unsigned int totalNNIApplied = 0;
    unsigned int numSteps = 0;
    const int MAXSTEPS = leafNum;
//    unsigned int numInnerBranches = leafNum - 3;
    double curBestScore = candidateTrees.getBestScore();

//    if (isMixlen())
//        optimizeBranches();

    Branches nniBranches;
    Branches nonNNIBranches;
    vector<NNIMove> positiveNNIs;
    vector<NNIMove> appliedNNIs;
    SplitIntMap tabuSplits;
    if (!initTabuSplits.empty()) {
        tabuSplits = initTabuSplits;
    }
    std::string task = "Optimizing NNI";
    if (0<strlen(context)) {
        task += " for ";
        task += context;
    }
    initProgress(MAXSTEPS, task, "done", "step", true);
    double originalScore = curScore;
    for (numSteps = 1; numSteps <= MAXSTEPS; numSteps++) {

//        cout << "numSteps = " << numSteps << endl;
        double oldScore = curScore;
        if (save_all_trees == 2) {
            saveCurrentTree(curScore); // BQM: for new bootstrap
        }
        if (verbose_mode >= VB_DEBUG && !progress_display::getProgressDisplay()) {
            LOG_LINE(VB_DEBUG, "Doing NNI round " << numSteps);
            if (isSuperTree()) {
                ((PhyloSuperTree*) this)->printMapInfo();
            }
        }

        // save all current branch lengths
        DoubleVector lenvec;
        saveBranchLengths(lenvec);

        // for super tree
        initPartitionInfo();

        nniBranches.clear();
        nonNNIBranches.clear();

        bool startSpeedNNI;
        // When tabu and speednni are combined, speednni is only start from third steps
        if (!initTabuSplits.empty() && numSteps < 3) {
            startSpeedNNI = false;
        } else if (speedNNI && !appliedNNIs.empty()) {
            startSpeedNNI = true;
        } else {
            startSpeedNNI = false;
        }

        if (startSpeedNNI) {
            // speedNNI option: only evaluate NNIs that are 2 branches away from the previously applied NNI
            Branches filteredNNIBranches;
            filterNNIBranches(appliedNNIs, filteredNNIBranches);
            for (Branches::iterator it = filteredNNIBranches.begin(); it != filteredNNIBranches.end(); it++) {
                Branch         curBranch = it->second;
                PhyloNeighbor* nei       = (PhyloNeighbor*) curBranch.first->findNeighbor(curBranch.second);
                Split*         curSplit  = nei->split;
                bool tabu = false;
                bool stable = false;
                if (!tabuSplits.empty()) {
                    int value;
                    if (tabuSplits.findSplit(curSplit, value) != NULL)
                        tabu = true;
                }
                if (!candidateTrees.getCandSplits().empty()) {
                    int value;
                    if (candidateTrees.getCandSplits().findSplit(curSplit, value) != NULL)
                        stable = true;

                }
                if (!tabu && !stable) {
                    int branchID =  pairInteger(curBranch.first->id, curBranch.second->id);
                    nniBranches.insert(pair<int, Branch>(branchID, curBranch));
                }
            }
        } else {
            getNNIBranches(tabuSplits, candidateTrees.getCandSplits(), nonNNIBranches, nniBranches);
        }

        if (!tabuSplits.empty()) {
            tabuSplits.clear();
        }

        positiveNNIs.clear();
        evaluateNNIs(nniBranches, positiveNNIs);

        if (positiveNNIs.size() == 0) {
            if (!nonNNIBranches.empty() && totalNNIApplied == 0) {
                evaluateNNIs(nonNNIBranches, positiveNNIs);
                if (positiveNNIs.size() == 0) {
                    break;
                }
            } else {
                break;
            }
        }

        /* sort all positive NNI moves (ASCENDING) */
        sort(positiveNNIs.begin(), positiveNNIs.end());

        /* remove conflicting NNIs */
        appliedNNIs.clear();
        getCompatibleNNIs(positiveNNIs, appliedNNIs);

        // do non-conflicting positive NNIs
        doNNIs(appliedNNIs);
        curScore = optimizeAllBranches(1, params->loglh_epsilon, PLL_NEWZPERCYCLE);
        const double loglh_tolerance = 0.1;
        auto expected_score = appliedNNIs.at(0).newloglh;

        if (curScore < expected_score - params->loglh_epsilon) {
            LOG_LINE( VB_MAX, "Tree getting worse. curScore = " << curScore
                     << " / best score = " << expected_score);
            // tree cannot be worse if only 1 NNI is applied
            double dodgy_score = curScore;
            if (appliedNNIs.size() > 1) {
                // revert all applied NNIs
                doNNIs(appliedNNIs);
                restoreBranchLengths(lenvec);
                clearAllPartialLH();
                // only do the best NNI
                appliedNNIs.resize(1);
                doNNIs(appliedNNIs);
                curScore = optimizeAllBranches(1, params->loglh_epsilon, PLL_NEWZPERCYCLE);
            }
            if (curScore < expected_score - loglh_tolerance ) {
                //In large trees this is possible, because optimizeAllBranches
                //might not have behaved as expected.
                LOG_LINE(VB_MIN, "Current likelihood (" << curScore << ")"
                        << " was less than expected (" << expected_score << ")"
                        << " by more (" << (expected_score - curScore) << ")"
                        << " than " << loglh_tolerance << ".");
            } else {
                LOG_LINE(VB_MED, "Post-NNI Likelihood score (" << dodgy_score << ")"
                    << " was less than expected (" << expected_score << ")"
                    << " by more (" << (expected_score - dodgy_score) << ")"
                    << " than epsilon (" << params->loglh_epsilon << ")."
                    << " After applying only the 1st NNI,"
                    << " Likelihood score is now: " << curScore);
            }
            if (!params->ignore_any_errors) {
                ASSERT(curScore >= expected_score - loglh_tolerance);
            }
            totalNNIApplied++;
        } else {
            totalNNIApplied += appliedNNIs.size();
        }
        if (curScore - oldScore <  params->loglh_epsilon)
        {
            break;
        }
        if (params->snni && (curScore > curBestScore + 0.1)) {
            curBestScore = curScore;
        }
        if (Params::getInstance().write_intermediate_trees && save_all_trees != 2) {
            printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
        }
        if (Params::getInstance().writeDistImdTrees) {
            intermediateTrees.update(getTreeString(), curScore);
        }
        trackProgress(1);
    }
    doneProgress();

    if (totalNNIApplied == 0 && verbose_mode >= VB_MED) {
        cout << "NOTE: Input tree is already NNI-optimal" << endl;
    }

    if (numSteps == MAXSTEPS) {
        cout << "WARNING: NNI search needs unusual large number of steps (" << numSteps << ") to converge!" << endl;
    }

    if(curScore < originalScore - params->loglh_epsilon){
        cout << "AAAAAAAAAAAAAAAAAAA: " << curScore << "\t" << originalScore << "\t" << curScore - originalScore << endl;

    }
    return make_pair(numSteps, totalNNIApplied);
}

void IQTree::filterNNIBranches(vector<NNIMove> &appliedNNIs, Branches &nniBranches) {
    for (vector<NNIMove>::iterator it = appliedNNIs.begin(); it != appliedNNIs.end(); it++) {
        Branch curBranch;
        curBranch.first = it->node1;
        curBranch.second = it->node2;
        int branchID = pairInteger(it->node1->id, it->node2->id);
        if (nniBranches.find(branchID) == nniBranches.end())
            nniBranches.insert(pair<int,Branch>(branchID, curBranch));
        getSurroundingInnerBranches(it->node1, it->node2, 2, nniBranches);
        getSurroundingInnerBranches(it->node2, it->node1, 2, nniBranches);
    }
}

double IQTree::pllOptimizeNNI(int &totalNNICount, int &nniSteps, SearchInfo &searchinfo) {
    if((globalParams->online_bootstrap == PLL_TRUE) && (globalParams->gbo_replicates > 0)) {
        pllInitUFBootData();
    }
    searchinfo.numAppliedNNIs = 0;
    searchinfo.curLogl = curScore;
    //cout << "curLogl: " << searchinfo.curLogl << endl;
    const int MAX_NNI_STEPS = aln->getNSeq();
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
        cout << "WARNING: NNI search needs unusual large number of steps (" << MAX_NNI_STEPS << ") to converge!" << endl;
    }

    if (searchinfo.numAppliedNNIs == 0) {
        cout << "NOTE: Tree is already NNI-optimized" << endl;
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

//            pllUFBootDataPtr->treels = pllHashInit(max_candidate_trees);
//            pllUFBootDataPtr->treels_size = max_candidate_trees; // track size of treels_logl, treels_newick, treels_ptnlh

//            pllUFBootDataPtr->treels_logl =
//                (double *) malloc(max_candidate_trees * (sizeof(double)));
//            if(!pllUFBootDataPtr->treels_logl) outError("Not enough dynamic memory!");



//            pllUFBootDataPtr->treels_ptnlh =
//                (double **) malloc(max_candidate_trees * (sizeof(double *)));
//            if(!pllUFBootDataPtr->treels_ptnlh) outError("Not enough dynamic memory!");
//            memset(pllUFBootDataPtr->treels_ptnlh, 0, max_candidate_trees * (sizeof(double *)));

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


            pllUFBootDataPtr->boot_logl =
                (double *) malloc(params->gbo_replicates * (sizeof(double)));
            if(!pllUFBootDataPtr->boot_logl) outError("Not enough dynamic memory!");
            for(int i = 0; i < params->gbo_replicates; i++)
                pllUFBootDataPtr->boot_logl[i] = -DBL_MAX;

            pllUFBootDataPtr->boot_counts =
                (int *) malloc(params->gbo_replicates * (sizeof(int)));
            if(!pllUFBootDataPtr->boot_counts) outError("Not enough dynamic memory!");
            memset(pllUFBootDataPtr->boot_counts, 0, params->gbo_replicates * (sizeof(int)));

            pllUFBootDataPtr->boot_trees.resize(params->gbo_replicates, "");
            pllUFBootDataPtr->duplication_counter = 0;
        }
    }
//    pllUFBootDataPtr->max_candidate_trees = max_candidate_trees;
    pllUFBootDataPtr->save_all_trees = save_all_trees;
    pllUFBootDataPtr->logl_cutoff = logl_cutoff;
    pllUFBootDataPtr->n_patterns = pllAlignment->sequenceLength;
}

void IQTree::pllDestroyUFBootData(){
    if(pll2iqtree_pattern_index){
        delete [] pll2iqtree_pattern_index;
        pll2iqtree_pattern_index = NULL;
    }

    if(params->online_bootstrap && params->gbo_replicates > 0){
        pllHashDestroy(&(pllUFBootDataPtr->treels), rax_free);

        free(pllUFBootDataPtr->treels_logl);


        for(int i = 0; i < pllUFBootDataPtr->treels_size; i++)
            if(pllUFBootDataPtr->treels_ptnlh[i])
                free(pllUFBootDataPtr->treels_ptnlh[i]);
        free(pllUFBootDataPtr->treels_ptnlh);

        for(int i = 0; i < params->gbo_replicates; i++)
            free(pllUFBootDataPtr->boot_samples[i]);
        free(pllUFBootDataPtr->boot_samples);

        free(pllUFBootDataPtr->boot_logl);

        free(pllUFBootDataPtr->boot_counts);

    }
    free(pllUFBootDataPtr);
    pllUFBootDataPtr = NULL;
}


void IQTree::doNNIs(const vector<NNIMove> &compatibleNNIs, bool changeBran) {
    for (auto it = compatibleNNIs.begin(); it != compatibleNNIs.end(); it++) {
        doNNI(*it);
        if (!params->leastSquareNNI && changeBran) {
            // apply new branch lengths
            changeNNIBrans(*it);
        }
    }
    // 2015-10-14: has to reset this pointer when read in
    current_it = current_it_back = NULL;

}


void IQTree::getCompatibleNNIs(vector<NNIMove> &nniMoves, vector<NNIMove> &compatibleNNIs) {
    compatibleNNIs.clear();
    for (vector<NNIMove>::iterator it1 = nniMoves.begin(); it1 != nniMoves.end(); it1++) {
        bool select = true;
        for (vector<NNIMove>::iterator it2 = compatibleNNIs.begin(); it2 != compatibleNNIs.end(); it2++) {
            if ((*it1).node1 == (*(it2)).node1
                    || (*it1).node2 == (*(it2)).node1
                    || (*it1).node1 == (*(it2)).node2
                    || (*it1).node2 == (*(it2)).node2) {
                select = false;
                break;
            }
        }
        if (select) {
            compatibleNNIs.push_back(*it1);
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

void IQTree::evaluateNNIs(Branches &nniBranches, vector<NNIMove>  &positiveNNIs) {
    for (Branches::iterator it = nniBranches.begin(); it != nniBranches.end(); it++) {
        NNIMove nni = getBestNNIForBran((PhyloNode*) it->second.first, (PhyloNode*) it->second.second, nullptr);
        if (nni.newloglh > curScore) {
            positiveNNIs.push_back(nni);
        }

        // synchronize tree during optimization step
        if (MPIHelper::getInstance().isMaster() && candidateset_changed.size() > 0
            && MPIHelper::getInstance().gotMessage()) {
            syncCurrentTree();
        }
    }
}

//Branches IQTree::getReducedListOfNNIBranches(Branches &previousNNIBranches) {
//    Branches resBranches;
//    for (Branches::iterator it = previousNNIBranches.begin(); it != previousNNIBranches.end(); it++) {
//        getSurroundingInnerBranches(it->second.first, it->second.second, 2, resBranches);
//        getSurroundingInnerBranches(it->second.second, it->second.first, 2, resBranches);
//    }
//}

double IQTree::optimizeNNIBranches(Branches &nniBranches) {
    for (Branches::iterator it = nniBranches.begin(); it != nniBranches.end(); it++) {
        optimizeOneBranch((PhyloNode*) it->second.first, (PhyloNode*) it->second.second, true, PLL_NEWZPERCYCLE);
    }
    curScore = computeLikelihoodFromBuffer();
    return curScore;
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
    PhyloNodeVector nodes1, nodes2;
    int i;
    double cur_lh = curScore;
    vector<IntBranchInfo> int_branches;

    getInternalBranches(nodes1, nodes2);
    assert(nodes1.size() == leafNum - 3 && nodes2.size() == leafNum - 3);

    for (i = 0; i < leafNum - 3; i++) {
        IntBranchInfo int_branch;
        PhyloNeighbor *node12_it = nodes1[i]->findNeighbor(nodes2[i]);
        //PhyloNeighbor *node21_it = nodes2[i]->findNeighbor(nodes1[i]);
        int_branch.lh_contribution = cur_lh - computeLikelihoodZeroBranch(node12_it, nodes1[i]);
        if (int_branch.lh_contribution < 0.0)
            int_branch.lh_contribution = 0.0;
        if (int_branch.lh_contribution < fabs(nni_cutoff)) {
            int_branch.node1 =  nodes1[i];
            int_branch.node2 =  nodes2[i];
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
            PhyloNeighbor *node12_it = node1->findNeighbor(node2);
            PhyloNeighbor *node21_it = node2->findNeighbor(node1);
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

    if (logl_cutoff != 0.0 && cur_logl < logl_cutoff - 1.0)
        return;
//    treels_logl.push_back(cur_logl);
//    num_trees_for_rell++;

    if (Params::getInstance().write_intermediate_trees)
        printTree(out_treels, WT_NEWLINE | WT_BR_LEN);

    int nptn = getAlnNPattern();

#ifdef BOOT_VAL_FLOAT
    int maxnptn = get_safe_upper_limit_float(nptn);
    BootValType *pattern_lh = aligned_alloc<BootValType>(maxnptn);
    memset(pattern_lh, 0, maxnptn*sizeof(BootValType));
    double *pattern_lh_orig = aligned_alloc<double>(nptn);
    computePatternLikelihood(pattern_lh_orig, &cur_logl);
    for (int i = 0; i < nptn; i++)
        pattern_lh[i] = (float)pattern_lh_orig[i];
#else
    int maxnptn = get_safe_upper_limit(nptn);
    BootValType *pattern_lh = aligned_alloc<BootValType>(maxnptn);
    memset(pattern_lh, 0, maxnptn*sizeof(BootValType));
    computePatternLikelihood(pattern_lh, &cur_logl);
#endif


    if (boot_samples.empty()) {
        // for runGuidedBootstrap
    } else {
        // online bootstrap
//        int ptn;
//        int updated = 0;
//        int nsamples = boot_samples.size();
        ostringstream ostr;
        string tree_str;
        setRootNode(params->root);
        if (params->print_ufboot_trees == 2)
            printTree(ostr, WT_TAXON_ID + WT_SORT_TAXA + WT_BR_LEN + WT_BR_LEN_SHORT);
        else
            printTree(ostr, WT_TAXON_ID + WT_SORT_TAXA);
        tree_str = ostr.str();

    #ifdef _OPENMP
        int rand_seed = random_int(1000);
        #pragma omp parallel
        {
        int *rstream;
        init_random(rand_seed + omp_get_thread_num(), false, &rstream);
        #pragma omp for
    #else
        int *rstream = randstream;
    #endif
        for (int sample = sample_start; sample < sample_end; sample++) {
            double rell = 0.0;

            {
                // SSE optimized version of the above loop
                BootValType *boot_sample = boot_samples[sample];

                BootValType res = (this->*dotProduct)(pattern_lh, boot_sample, nptn);

                rell = res;
            }

            bool better = rell > boot_logl[sample] + params->ufboot_epsilon;
            if (!better && rell > boot_logl[sample] - params->ufboot_epsilon) {
                better = (random_double(rstream) <= 1.0 / (boot_counts[sample] + 1));
            }
            if (better) {
                if (rell <= boot_logl[sample] + params->ufboot_epsilon) {
                    boot_counts[sample]++;
                } else {
                    boot_counts[sample] = 1;
                }
                boot_logl[sample] = max(boot_logl[sample], rell);
                boot_orig_logl[sample] = cur_logl;
                boot_trees[sample] = tree_str;
            }
        }
    #ifdef _OPENMP
        finish_random(rstream);
        }
    #endif
    }
    if (Params::getInstance().print_tree_lh) {
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
        for (size_t i = 0; i < getAlnNSite(); ++i)
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
        node = getRoot();
    }
    if (dad && !node->isLeaf() && !dad->isLeaf()) {
        double *pat_lh1 = new double[aln->getNPattern()];
        double *pat_lh2 = new double[aln->getNPattern()];
        double lh1, lh2;
        computeNNIPatternLh(curScore, lh1, pat_lh1, lh2, pat_lh2, node, dad);
        delete[] pat_lh2;
        delete[] pat_lh1;
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        saveNNITrees(child, node);
    }
}

void IQTree::summarizeBootstrap(Params &params, MTreeSet &trees) {
    int sum_weights = trees.sumTreeWeights();
    int i;
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
    trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, NULL, false); // do not sort taxa

    if (verbose_mode >= VB_MED)
    	cout << sg.size() << " splits found" << endl;

    sg.scaleWeight(1.0 / trees.sumTreeWeights(), false, 4);
    string out_file;
    out_file = params.out_prefix;
    out_file += ".splits";
    if (params.print_splits_file) {
        sg.saveFile(out_file.c_str(), IN_OTHER, true);
        cout << "Split supports printed to star-dot file " << out_file << endl;
    }
    // compute the percentage of appearance
    sg.scaleWeight(100.0, true);
    //    printSplitSet(sg, hash_ss);
    //sg.report(cout);
    //cout << "Creating " << RESAMPLE_NAME << " support values..." << endl;
//    stringstream tree_stream;
//    printTree(tree_stream, WT_TAXON_ID | WT_BR_LEN);
//    MExtTree mytree;
//    mytree.readTree(tree_stream, rooted);
//    mytree.assignLeafID();
    assignLeafNameByID();
    createBootstrapSupport(taxname, trees, hash_ss, NULL);

    // now write resulting tree with supports
//    tree_stream.seekp(0, ios::beg);
//    mytree.printTree(tree_stream);
//
//    // now read resulting tree
//    tree_stream.seekg(0, ios::beg);
//    freeNode();
//    // RARE BUG FIX: to avoid cases that identical seqs were removed and leaf name happens to be IDs
//    MTree::readTree(tree_stream, rooted);

    assignLeafNames();

    // 2017-12-05: commented out, not sure why doing this, which disort branch lengths
//    if (isSuperTree()) {
//        ((PhyloSuperTree*) this)->mapTrees();
//    } else {
//        initializeAllPartialLh();
//        clearAllPartialLH();
//    }

    if (!save_all_trees) {
        out_file = params.out_prefix;
        out_file += ".suptree";

        printTree(out_file.c_str());
        cout << "Tree with assigned support written to " << out_file << endl;
    }
    
    if (params.print_splits_nex_file) {
        out_file = params.out_prefix;
        out_file += ".splits.nex";
        sg.saveFile(out_file.c_str(), IN_NEXUS, false);
        cout << "Split supports printed to NEXUS file " << out_file << endl;
    }

    /*
     out_file = params.out_prefix;
     out_file += ".supval";
     writeInternalNodeNames(out_file);

     cout << "Support values written to " << out_file << endl;
     */


}

void IQTree::writeUFBootTrees(Params &params) {
    MTreeSet trees;
//    IntVector tree_weights;
    int i, j;
    string filename = params.out_prefix;
    filename += ".ufboot";
    ofstream out(filename.c_str());

    trees.init(boot_trees, rooted);
    for (i = 0; i < trees.size(); i++) {
        NodeVector taxa;
        // change the taxa name from ID to real name
        trees[i]->getOrderedTaxa(taxa);
        for (j = 0; j < taxa.size() - (int)rooted; j++)
            taxa[j]->name = aln->getSeqName(taxa[j]->id);
        if (removed_seqs.size() > 0) {
            // reinsert removed seqs into each tree
            trees[i]->insertTaxa(removed_seqs, twin_seqs);
        }
        // now print to file
        for (j = 0; j < trees.tree_weights[i]; j++)
            if (params.print_ufboot_trees == 1)
                trees[i]->printTree(out, WT_NEWLINE);
            else
                trees[i]->printTree(out, WT_NEWLINE + WT_BR_LEN);
    }
    cout << "UFBoot trees printed to " << filename << endl;
    out.close();
}

void IQTree::summarizeBootstrap(Params &params) {
    setRootNode(params.root);
    MTreeSet trees;
    trees.init(boot_trees, rooted);
    summarizeBootstrap(params, trees);
}

void IQTree::summarizeBootstrap(SplitGraph &sg) {
    MTreeSet trees;
    //SplitGraph sg;
    trees.init(boot_trees, rooted);
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
    trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, NULL, false); // do not sort taxa
}

void IQTree::pllConvertUFBootData2IQTree(){
    // duplication_counter
    duplication_counter = pllUFBootDataPtr->duplication_counter;
    //treels_logl
//    treels_logl.clear();
//    for(int i = 0; i < pllUFBootDataPtr->candidate_trees_count; i++)
//        treels_logl.push_back(pllUFBootDataPtr->treels_logl[i]);

    //boot_trees
    boot_trees.clear();
    for(int i = 0; i < params->gbo_replicates; i++)
        boot_trees.push_back(pllUFBootDataPtr->boot_trees[i]);

}

double computeCorrelation(IntVector &ix, IntVector &iy) {

    ASSERT(ix.size() == iy.size());
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

    return corr;
}

//void IQTree::addPositiveNNIMove(NNIMove myMove) {
//    plusNNIs.push_back(myMove);
//}

void IQTree::printResultTree(string suffix) {
    if (MPIHelper::getInstance().isWorker()) {
        return;
//        stringstream processTreeFile;
//        processTreeFile << tree_file_name << "." << MPIHelper::getInstance().getProcessID();
//        tree_file_name = processTreeFile.str();
    }
    if (params->suppress_output_flags & OUT_TREEFILE)
        return;
    setRootNode(params->root, true);
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    if (suffix.compare("") != 0) {
        tree_file_name += "." + suffix;
    }
    printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    if (verbose_mode >= VB_MED)
        cout << "Best tree printed to " << tree_file_name << endl;
    setRootNode(params->root, false);
}

void IQTree::printResultTree(ostream &out) {
    setRootNode(params->root);
    printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
}

void IQTree::printBestCandidateTree() {
    if (MPIHelper::getInstance().isWorker())
        return;
    if (params->suppress_output_flags & OUT_TREEFILE)
        return;
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    readTreeString(candidateTrees.getBestTreeStrings(1)[0]);
    setRootNode(params->root);
    printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    if (verbose_mode >= VB_MED)
        cout << "Best tree printed to " << tree_file_name << endl;
}


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
    double *pattern_lh = NULL;
    double logl = curScore;

    if (params->print_tree_lh) {
        pattern_lh = new double[getAlnNPattern()];
        computePatternLikelihood(pattern_lh, &logl);
    }

    if (Params::getInstance().write_intermediate_trees)
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
        for (size_t i = 0; i < aln->getNSite(); ++i)
            out_sitelh << "\t" << pattern_lh[aln->getPatternID(i)];
        out_sitelh << endl;
        delete[] pattern_lh;
    }
    if (params->write_intermediate_trees == 1 && save_all_trees != 1) {
        return;
    }
    int x = save_all_trees;
    save_all_trees = 2;
    // TODO Why is evalNNI() is called in this function?
    //evalNNIs();
    Branches innerBranches;
    vector<NNIMove> positiveNNIs;
    getInnerBranches(innerBranches);
    evaluateNNIs(innerBranches, positiveNNIs);
    save_all_trees = x;
}


void IQTree::convertNNI2Splits(SplitIntMap &nniSplits, int numNNIs, vector<NNIMove> &compatibleNNIs) {
    for (int i = 0; i < numNNIs; i++) {
        Split *sp = new Split(*getSplit(compatibleNNIs[i].node1, compatibleNNIs[i].node2));
        if (sp->shouldInvert()) {
            sp->invert();
        }
        nniSplits.insertSplit(sp, 1);
    }
}

double IQTree::getBestScore() {
    return candidateTrees.getBestScore();
}

vector<string> IQTree::getBestTrees(int numTrees) {
    return candidateTrees.getBestTreeStrings(numTrees);
}


/*******************************************
    MPI stuffs
*******************************************/

void IQTree::syncCandidateTrees(int nTrees, bool updateStopRule) {
    if (MPIHelper::getInstance().getNumProcesses() == 1)
        return;

#ifdef _IQTREE_MPI
    // gather trees to Master

    Checkpoint *ckp = new Checkpoint;

    if (MPIHelper::getInstance().isMaster()) {
        // update candidate set at master
        int trees = 0;
        for (int w = 1; w < MPIHelper::getInstance().getNumProcesses(); w++) {
            int worker = MPIHelper::getInstance().recvCheckpoint(ckp);
            CandidateSet cset;
            cset.setCheckpoint(ckp);
            cset.restoreCheckpoint();
            for (CandidateSet::iterator it = cset.begin(); it != cset.end(); it++)
                addTreeToCandidateSet(it->second.tree, it->second.score, updateStopRule, worker);
            trees += ckp->size();
            ckp->clear();
        }
        cout << "Master: " << trees << " candidate trees gathered from workers" << endl;
        // get the best candidate trees
        int numTrees = max(nTrees, MPIHelper::getInstance().getNumProcesses());
        CandidateSet bestCandidates = candidateTrees.getBestCandidateTrees(numTrees);
        int saved_numNNITrees = params->numNNITrees;
        params->numNNITrees = numTrees;
        bestCandidates.setCheckpoint(ckp);
        bestCandidates.saveCheckpoint();
        params->numNNITrees = saved_numNNITrees;
    } else {
        // send candidate set to master
        CandidateSet cset = candidateTrees.getBestCandidateTrees();
        cset.setCheckpoint(ckp);
        cset.saveCheckpoint();
        MPIHelper::getInstance().sendCheckpoint(ckp, PROC_MASTER);
        cout << "Worker " << MPIHelper::getInstance().getProcessID() << ": " << ckp->size() << " candidate trees sent to master" << endl;
        ckp->clear();
    }

    if (updateStopRule && stop_rule.meetStopCondition(stop_rule.getCurIt(), 0.0)) {
        // 2020-04-30: send stop signal
        ckp->putBool("stop", true);
    }
    
    // broadcast candidate trees from master to worker
    MPIHelper::getInstance().broadcastCheckpoint(ckp);
    cout << ckp->size() << " trees broadcasted to workers" << endl;

    if (MPIHelper::getInstance().isWorker()) {
        // update candidate set at worker
        CandidateSet cset;
        cset.setCheckpoint(ckp);
        cset.restoreCheckpoint();
        for (CandidateSet::iterator it = cset.begin(); it != cset.end(); it++)
            addTreeToCandidateSet(it->second.tree, it->second.score, false, PROC_MASTER);
        
        // 2020-04-40: check stop signal
        if (ckp->getBool("stop")) {
            cout << "Worker " << MPIHelper::getInstance().getProcessID() << " gets STOP message!" << endl;
            stop_rule.shouldStop();
        }
    }

    delete ckp;
#endif
}

void IQTree::syncCurrentTree() {
    if (MPIHelper::getInstance().getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    //------ BLOCKING COMMUNICATION ------//
    Checkpoint *checkpoint = new Checkpoint;
    string tree;
    double score;

    if (MPIHelper::getInstance().isMaster()) {
        // master: receive tree from WORKERS
        int worker = MPIHelper::getInstance().recvCheckpoint(checkpoint);
        MPIHelper::getInstance().increaseTreeReceived();
        CKP_RESTORE(tree);
        CKP_RESTORE(score);
        int pos = addTreeToCandidateSet(tree, score, true, worker);
        if (pos >= 0 && pos < params->popSize) {
            // candidate set is changed, update for other workers
            for (int w = 0; w < candidateset_changed.size(); w++)
                if (w != worker)
                    candidateset_changed[w] = true;
        }

        if (boot_samples.size() > 0) {
            restoreUFBoot(checkpoint);
        }

        // send candidate trees to worker
        checkpoint->clear();
        if (boot_samples.size() > 0)
            CKP_SAVE(logl_cutoff);
        if (candidateset_changed[worker]) {
            CandidateSet cset = candidateTrees.getBestCandidateTrees(Params::getInstance().popSize);
            cset.setCheckpoint(checkpoint);
            cset.saveCheckpoint();
            candidateset_changed[worker] = false;
            MPIHelper::getInstance().increaseTreeSent(Params::getInstance().popSize);
        }
        MPIHelper::getInstance().sendCheckpoint(checkpoint, worker);
    } else {
        // worker: always send tree to MASTER
        tree = getTreeString();
        score = curScore;
        CKP_SAVE(tree);
        CKP_SAVE(score);
        if (boot_samples.size() > 0) {
            saveUFBoot(checkpoint);
        }
        MPIHelper::getInstance().sendCheckpoint(checkpoint, PROC_MASTER);
        MPIHelper::getInstance().increaseTreeSent();

        // now receive the candidate set
        MPIHelper::getInstance().recvCheckpoint(checkpoint, PROC_MASTER);
        if (checkpoint->getBool("stop")) {
            cout << "Worker " << MPIHelper::getInstance().getProcessID() << " gets STOP message!" << endl;
            stop_rule.shouldStop();
        } else {
            CandidateSet cset;
            cset.setCheckpoint(checkpoint);
            cset.restoreCheckpoint();
            for (CandidateSet::iterator it = cset.begin(); it != cset.end(); it++)
                addTreeToCandidateSet(it->second.tree, it->second.score, false, MPIHelper::getInstance().getProcessID());
            MPIHelper::getInstance().increaseTreeReceived(cset.size());
            if (boot_samples.size() > 0)
                CKP_RESTORE(logl_cutoff);
        }
    }

    delete checkpoint;

#endif
}

void IQTree::sendStopMessage() {
    if (MPIHelper::getInstance().getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI

    Checkpoint *checkpoint = new Checkpoint;
    checkpoint->putBool("stop", true);
    stringstream ss;
    checkpoint->dump(ss);
    string str = ss.str();
    string tree;
    double score;

    cout << "Sending STOP message to workers" << endl;

    // send STOP message to all processes
    if (MPIHelper::getInstance().isMaster()) {
        // repeatedly send stop message to all workers
        for (int w = 1; w < MPIHelper::getInstance().getNumProcesses(); w++) {
//            string buf;
//            int worker = MPIHelper::getInstance().recvString(buf);
            checkpoint->clear();
            int worker = MPIHelper::getInstance().recvCheckpoint(checkpoint);
            MPIHelper::getInstance().increaseTreeReceived();
            CKP_RESTORE(tree);
            CKP_RESTORE(score);
            addTreeToCandidateSet(tree, score, true, worker);
            MPIHelper::getInstance().sendString(str, worker, TREE_TAG);
        }
    }

    delete checkpoint;

    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void PhyloTree::warnNumThreads() const {
    if (num_threads <= 1)
        return;
    // return if -nt AUTO
    if (params->num_threads == 0)
        return;
    size_t nptn = getAlnNPattern();
    if (nptn < num_threads*vector_size) {
        outError("Too many threads for short alignments, please reduce number of threads or use -T AUTO to determine it.");
    }
    if (nptn < num_threads*400/aln->num_states) {
        outWarning("Number of threads seems too high for short alignments. Use -T AUTO to determine best number of threads.");
    }
}

int PhyloTree::ensureNumberOfThreadsIsSet(Params *params, bool suppressAnyThreadCountWarnings) {
    #ifdef _OPENMP
        if (num_threads <= 0 ) {
            int bestThreads = testNumThreads();
            omp_set_num_threads(bestThreads);
            if (params!=nullptr) {
                params->num_threads = bestThreads;
            }
        } else if (!suppressAnyThreadCountWarnings) {
            warnNumThreads();
        }
        return num_threads;
    #else
        return 1;
    #endif
}

int PhyloTree::testNumThreads() {
#ifndef _OPENMP
    return 1;
#else
	int max_procs = min(countPhysicalCPUCores()/MPIHelper::getInstance().countSameHost(), params->num_threads_max);
    cout << "Measuring multi-threading efficiency up to " << max_procs << " CPU cores" << endl;
    DoubleVector runTimes;
    int bestProc = 0;
    double saved_curScore = curScore;
    int num_iter = 1;

    // generate different trees
    int tree;
    double min_time = max_procs; // minimum time in seconds
    StrVector trees;
    trees.push_back(getTreeString());
    setLikelihoodKernel(sse);

    initProgress(max_procs, "Determining AUTO threadcount", "tried", "threadcount", true);
    for (int proc = 1; proc <= max_procs; ++proc) {
        omp_set_num_threads(proc);
        setNumThreads(proc);
        initializeAllPartialLh();

        double beginTime = getRealTime();
        double runTime, logl;

        for (tree = 0; tree < trees.size(); tree++) {
            readTreeString(trees[tree]);
            logl = optimizeAllBranches(num_iter);
            runTime = getRealTime() - beginTime;

            // too fast, increase number of iterations
            if (runTime*10 < min_time && proc == 1 && tree == 0) {
                int new_num_iter = 10;
                hideProgress();
                cout << "Increase to " << new_num_iter << " rounds for branch lengths" << endl;
                showProgress();
                logl = optimizeAllBranches(new_num_iter - num_iter);
                num_iter = new_num_iter;
                runTime = getRealTime() - beginTime;
            }

            // considering at least 2 trees
            if ((runTime < min_time && proc == 1) || trees.size() == 1) {
                // time not reached, add more tree
//                readTreeString(trees[0]);
//                doRandomNNIs();
                generateRandomTree(YULE_HARDING);
                wrapperFixNegativeBranch(true);
                trees.push_back(getTreeString());
            }
            curScore = saved_curScore;
        }

        if (proc == 1) {
            hideProgress();
            cout << trees.size() << " trees examined" << endl;
            showProgress();
        }

        deleteAllPartialLh();

        runTimes.push_back(runTime);
        double speedup = runTimes[0] / runTime;

        hideProgress();
        cout << "Threads: " << proc << " / Time: " << runTime << " sec / Speedup: " << speedup
            << " / Efficiency: " << (int)round(speedup*100/proc) << "% / LogL: " << (int)logl << endl;
        showProgress();

        // break if too bad efficiency ( < 50%) or worse than than 10% of the best run time
        if (speedup*2 <= proc || (runTime > runTimes[bestProc]*1.1 && proc>1)) {
            break;
        }

        // update best threads if sufficient
        if (runTime <= runTimes[bestProc]*0.95) {
            bestProc = proc-1;
        }
        trackProgress(1.0);
    }
    doneProgress();
    readTreeString(trees[0]);

    cout << "BEST NUMBER OF THREADS: " << bestProc+1 << endl << endl;
    setNumThreads(bestProc+1);

    return bestProc+1;
#endif
}
