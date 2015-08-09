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
        outError("Tree must have 5 or more taxa with unique sequences!");
        
    quartet_info.resize(params->num_quartets);
    
    int qc[] = {0, 1, 2, 3,  0, 2, 1, 3,  0, 3, 1, 2};
    
    double onethird = 1.0/3.0;
    unsigned char treebits[] = {1, 2, 4};


#ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
#endif
int xxx=1;
    for (int qid = 0; qid < params->num_quartets; qid++) {
// fprintf(stderr, "%d\n", qid); 
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
            
        sort(quartet_info[qid].seqID, quartet_info[qid].seqID+4); // why do you sort them?!? HAS ;^)
            
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

	// determine likelihood order
	int qworder[3]; // local (thread-safe) vector for sorting

	if (quartet_info[qid].logl[0] > quartet_info[qid].logl[1]) {
		if(quartet_info[qid].logl[2] > quartet_info[qid].logl[0]) {
			qworder[0] = 2;
			qworder[1] = 0;
			qworder[2] = 1;		
		} else if (quartet_info[qid].logl[2] < quartet_info[qid].logl[1]) {
			qworder[0] = 0;
			qworder[1] = 1;
			qworder[2] = 2;		
		} else {
			qworder[0] = 0;
			qworder[1] = 2;
			qworder[2] = 1;		
		}
	} else {
		if(quartet_info[qid].logl[2] > quartet_info[qid].logl[1]) {
			qworder[0] = 2;
			qworder[1] = 1;
			qworder[2] = 0;		
		} else if (quartet_info[qid].logl[2] < quartet_info[qid].logl[0]) {
			qworder[0] = 1;
			qworder[1] = 0;
			qworder[2] = 2;		
		} else {
			qworder[0] = 1;
			qworder[1] = 2;
			qworder[2] = 0;		
		}
	}

	// compute Bayesian weights
	double temp;

	quartet_info[qid].qweight[0] = quartet_info[qid].logl[0];
	quartet_info[qid].qweight[1] = quartet_info[qid].logl[1];
	quartet_info[qid].qweight[2] = quartet_info[qid].logl[2];

	temp = quartet_info[qid].qweight[qworder[1]]-quartet_info[qid].qweight[qworder[0]];
	if(temp < -TP_MAX_EXP_DIFF)	/* possible, since 1.0+exp(>36) == 1.0 */
	   quartet_info[qid].qweight[qworder[1]] = 0.0;
	else
	   quartet_info[qid].qweight[qworder[1]] = exp(temp);

        temp = quartet_info[qid].qweight[qworder[2]]-quartet_info[qid].qweight[qworder[0]];
	if(temp < -TP_MAX_EXP_DIFF)	/* possible, since 1.0+exp(>36) == 1.0 */
	   quartet_info[qid].qweight[qworder[2]] = 0.0;
	else
	   quartet_info[qid].qweight[qworder[2]] = exp(temp);

	quartet_info[qid].qweight[qworder[0]] = 1.0;

	temp = quartet_info[qid].qweight[0] + quartet_info[qid].qweight[1] + quartet_info[qid].qweight[2];
	quartet_info[qid].qweight[0] = quartet_info[qid].qweight[0]/temp;
	quartet_info[qid].qweight[1] = quartet_info[qid].qweight[1]/temp;
	quartet_info[qid].qweight[2] = quartet_info[qid].qweight[2]/temp;

	// determine which of the three corners (only meaningful if seqIDs NOT sorted)
	if (treebits[qworder[0]] == 1) {
		quartet_info[qid].corner=0;
	} else {
		if (treebits[qworder[0]] == 2) {
			quartet_info[qid].corner=1;
		} else {
			quartet_info[qid].corner=2;
		}
	}

	// determine which of the 7 regions (only meaningful if seqIDs NOT sorted)
	double temp1, temp2, temp3;
	unsigned char discreteweight[3];
	double sqdiff[3];

	/* 100 distribution */
	temp1 = 1.0 - quartet_info[qid].qweight[qworder[0]];
	sqdiff[0] = temp1*temp1 +
		quartet_info[qid].qweight[qworder[1]]*quartet_info[qid].qweight[qworder[1]] +
		quartet_info[qid].qweight[qworder[2]]*quartet_info[qid].qweight[qworder[2]];
	discreteweight[0] = treebits[qworder[0]];

	/* 110 distribution */
	temp1 = 0.5 - quartet_info[qid].qweight[qworder[0]];
	temp2 = 0.5 - quartet_info[qid].qweight[qworder[1]];
	sqdiff[1] = temp1*temp1 + temp2*temp2 +
		quartet_info[qid].qweight[qworder[2]]*quartet_info[qid].qweight[qworder[2]];
	discreteweight[1] = treebits[qworder[0]] + treebits[qworder[1]];

	/* 111 distribution */
	temp1 = onethird - quartet_info[qid].qweight[qworder[0]];
	temp2 = onethird - quartet_info[qid].qweight[qworder[1]];
	temp3 = onethird - quartet_info[qid].qweight[qworder[2]];
	sqdiff[2] = temp1 * temp1 + temp2 * temp2 + temp3 * temp3;
	discreteweight[2] = (unsigned char) 7;

	/* sort in descending order */
	int sqorder[3]; // local (thread-safe) vector for sorting
	if (sqdiff[0] > sqdiff[1]) {
		if(sqdiff[2] > sqdiff[0]) {
			sqorder[0] = 2;
			sqorder[1] = 0;
			sqorder[2] = 1;		
		} else if (sqdiff[2] < sqdiff[1]) {
			sqorder[0] = 0;
			sqorder[1] = 1;
			sqorder[2] = 2;		
		} else {
			sqorder[0] = 0;
			sqorder[1] = 2;
			sqorder[2] = 1;		
		}
	} else {
		if(sqdiff[2] > sqdiff[1]) {
			sqorder[0] = 2;
			sqorder[1] = 1;
			sqorder[2] = 0;		
		} else if (sqdiff[2] < sqdiff[0]) {
			sqorder[0] = 1;
			sqorder[1] = 0;
			sqorder[2] = 2;		
		} else {
			sqorder[0] = 1;
			sqorder[1] = 2;
			sqorder[2] = 0;		
		}
	}


	// determine which of the 7 regions (only meaningful if seqIDs NOT sorted)
	unsigned char qpbranching = (unsigned char) discreteweight[sqorder[2]];

	if (qpbranching == 1) {
		quartet_info[qid].area=0; // LM_REG1 - top
	}
	if (qpbranching == 2) {
		quartet_info[qid].area=1; // LM_REG2 - right
	}
	if (qpbranching == 4) {
		quartet_info[qid].area=2; // LM_REG3 - left
	}

	if (qpbranching == 3) {
		quartet_info[qid].area=3; // LM_REG4
	}
	if (qpbranching == 6) {
		quartet_info[qid].area=4; // LM_REG5
	}
	if (qpbranching == 5) {
		quartet_info[qid].area=5; // LM_REG6
	}

	if (qpbranching == 7) {
		quartet_info[qid].area=6; // LM_REG7 - center 
	}

	
	
    } // end for num_quartets

} // end PhyloTree::computeQuartetLikelihoods


void PhyloTree::doLikelihoodMapping() {
    // TODO For Heiko: Please add code here
    vector<QuartetInfo> quartet_info;
    int areacount[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int cornercount[4] = {0, 0, 0, 0};
    int resolved, partly, unresolved;
    int qid;

    computeQuartetLikelihoods(quartet_info);
    for (qid = 0; qid < params->num_quartets; qid++) {
        areacount[quartet_info[qid].area]++;
        cornercount[quartet_info[qid].corner]++;
    }
    
    if (params->print_quartet_lh) {
        // print quartet file
        string filename = (string)params->out_prefix + ".quartetlh";
        ofstream out;
        out.open(filename.c_str());
        for (qid = 0; qid < params->num_quartets; qid++) {
            out << "(" << quartet_info[qid].seqID[0] << ","
                << quartet_info[qid].seqID[1] << ","
                << quartet_info[qid].seqID[2] << ","
                << quartet_info[qid].seqID[3] << ")"
                << "   " << quartet_info[qid].logl[0] 
                << "   " << quartet_info[qid].logl[1] 
                << "   " << quartet_info[qid].logl[2] << endl;
        }
        out.close();        
        cout << "Quartet log-likelihoods printed to " << filename << endl;
    }

    resolved   = areacount[0] + areacount[1] + areacount[2];
    partly     = areacount[3] + areacount[4] + areacount[5];
    unresolved = areacount[6];
	
//    fprintf(stdout, "LIKELIHOOD MAPPING ANALYSIS\n\n");
//    fprintf(stdout, "Number of quartets: %d (randomly drawn with replacement)\n\n", (resolved+partly+unresolved));
//    fprintf(stdout, "Overall quartet resolution:\n");
//    fprintf(stdout, "Number of fully resolved  quartets: %6d (= %.2f%%)\n", resolved, 100.0 * resolved/(resolved+partly+unresolved));
//    fprintf(stdout, "Number of partly resolved quartets: %6d (= %.2f%%)\n", partly, 100.0 * partly/(resolved+partly+unresolved));
//    fprintf(stdout, "Number of unresolved      quartets: %6d (= %.2f%%)\n\n", unresolved, 100.0 * unresolved/(resolved+partly+unresolved));

    cout << "\nOverall quartet resolution: (from " << (resolved+partly+unresolved) << " randomly drawn quartets)" << endl;
    cout << "Fully resolved quartets:  " << resolved   << " (= "
        << (double) resolved * 100.0   / (resolved+partly+unresolved) << "%)" << endl;
    cout << "Partly resolved quartets: " << partly     << " (= "
        << (double) partly * 100.0     / (resolved+partly+unresolved) << "%)" << endl;
    cout << "Unresolved quartets:      " << unresolved << " (= "
        << (double) unresolved * 100.0 / (resolved+partly+unresolved) << "%)" << endl << endl;

} // end PhyloTree::doLikelihoodMapping
