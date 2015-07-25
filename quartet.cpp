//
//  quartet.cpp
//  iqtree
//
//  Created by Minh Bui on 24/07/15.
//
//

#include <stdio.h>

#include "phylotree.h"
#include "phylosupertree.h"

void PhyloTree::computeQuartetLikelihoods(vector<QuartetInfo> &quartet_info) {

    if (leafNum <= 4) 
        outError("Tree must have 5 or more taxa");
        
    quartet_info.resize(params->num_quartets);
    
    int qc[] = {0, 1, 2, 3,  0, 2, 1, 3,  0, 3, 1, 2};
    
#ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
#endif
    for (int qid = 0; qid < params->num_quartets; qid++) {
        // uniformly draw 4 taxa
        quartet_info[qid].seqID[0] = random_int(leafNum);
        do {
            quartet_info[qid].seqID[1] = random_int(leafNum);
        } while (quartet_info[qid].seqID[1] == quartet_info[qid].seqID[0]);
        do {
            quartet_info[qid].seqID[2] = random_int(leafNum);
        } while (quartet_info[qid].seqID[2] == quartet_info[qid].seqID[0] || quartet_info[qid].seqID[2] == quartet_info[qid].seqID[1]);
        do {
            quartet_info[qid].seqID[3] = random_int(leafNum);
        } while (quartet_info[qid].seqID[3] == quartet_info[qid].seqID[0] || quartet_info[qid].seqID[3] == quartet_info[qid].seqID[1]
            || quartet_info[qid].seqID[3] == quartet_info[qid].seqID[2]);
            
        sort(quartet_info[qid].seqID, quartet_info[qid].seqID+4);
            
        // initialize sub-alignment and sub-tree
        Alignment *quartet_aln;
        PhyloTree *quartet_tree;
        if (aln->isSuperAlignment()) {
            quartet_aln = new SuperAlignment;
        } else {
            quartet_aln = new Alignment;
        }
        IntVector seq_id;
        seq_id.insert(seq_id.begin(), quartet_info[qid].seqID, quartet_info[qid].seqID+4);
        quartet_aln->extractSubAlignment(aln, seq_id, 0);
        if (isSuperTree()) {
            quartet_tree = new PhyloSuperTree((SuperAlignment*)quartet_aln, (PhyloSuperTree*)this);
        } else {
            quartet_tree = new PhyloTree(quartet_aln);
        }

        // set up parameters
        quartet_tree->setParams(params);
        quartet_tree->optimize_by_newton = params->optimize_by_newton;
        quartet_tree->setLikelihoodKernel(params->SSE);
        
        // set model and rate
        quartet_tree->setModelFactory(model_factory);
        quartet_tree->setModel(getModel());
        quartet_tree->setRate(getRate());
        // NOTE: we don't need to set phylo_tree in model and rate because parameters are not reoptimized
        
        
        
        // loop over 3 quartets to compute likelihood
        for (int k = 0; k < 3; k++) {
            string quartet_tree_str;
            quartet_tree_str = "(" + quartet_aln->getSeqName(qc[k*4]) + "," + quartet_aln->getSeqName(qc[k*4+1]) + ",(" + 
                quartet_aln->getSeqName(qc[k*4+2]) + "," + quartet_aln->getSeqName(qc[k*4+3]) + "));";
            quartet_tree->readTreeString(quartet_tree_str);
            quartet_tree->initializeAllPartialLh();
            quartet_tree->wrapperFixNegativeBranch(true);
            // optimize branch lengths with logl_epsilon=0.1 accuracy
            quartet_info[qid].logl[k] = quartet_tree->optimizeAllBranches(10, 0.1);
        }
        // reset model & rate so that they are not deleted
        quartet_tree->setModel(NULL);
        quartet_tree->setModelFactory(NULL);
        quartet_tree->setRate(NULL);

        delete quartet_tree;
        delete quartet_aln;
    }

}

void PhyloTree::doLikelihoodMapping() {
    // TODO For Heiko: Please add code here
    vector<QuartetInfo> quartet_info;
    computeQuartetLikelihoods(quartet_info);
    
    // print quartet file
    string filename = (string)params->out_prefix + ".quartetlh";
    ofstream out;
    out.open(filename.c_str());
    for (int qid = 0; qid < params->num_quartets; qid++) {
        out << "(" << quartet_info[qid].seqID[0] << ","
            << quartet_info[qid].seqID[1] << ","
            << quartet_info[qid].seqID[2] << ","
            << quartet_info[qid].seqID[3] << ")"
            << "   " << quartet_info[qid].logl[0] 
            << "   " << quartet_info[qid].logl[1] 
            << "   " << quartet_info[qid].logl[2] << endl;
    }
    out.close();
    
}
