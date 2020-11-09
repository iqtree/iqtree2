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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>
#include "tree/phylotree.h"
#include "tree/phylosupertree.h"
#include "tree/phylosupertreeplen.h"
#include "tree/phylosupertreeunlinked.h"
#include "phyloanalysis.h"
#include "alignment/alignment.h"
#include "alignment/superalignment.h"
#include "alignment/superalignmentunlinked.h"
#include "tree/iqtree.h"
#include "tree/phylotreemixlen.h"
#include "model/modelmarkov.h"
#include "model/modeldna.h"
#include "model/modelpomo.h"
#include "nclextra/myreader.h"
#include "model/rateheterogeneity.h"
#include "model/rategamma.h"
#include "model/rateinvar.h"
#include "model/rategammainvar.h"
//#include "modeltest_wrapper.h"
#include "model/modelprotein.h"
#include "model/modelbin.h"
#include "model/modelcodon.h"
#include "utils/stoprule.h"

#include "tree/mtreeset.h"
#include "tree/mexttree.h"
#include "model/ratemeyerhaeseler.h"
#include "whtest/whtest_wrapper.h"
#include "model/partitionmodel.h"
#include "model/modelmixture.h"
#include "model/modelfactorymixlen.h"
//#include "guidedbootstrap.h"
#include "model/modelset.h"
#include "utils/timeutil.h"
#include "tree/upperbounds.h"
#include "utils/MPIHelper.h"
#include "timetree.h"

#ifdef USE_BOOSTER
extern "C" {
#include "booster/booster.h"
}
#endif

#ifdef IQTREE_TERRAPHAST
    #include "terrace/terrace.h"
#endif

void reportReferences(Params &params, ofstream &out) {

    out << "To cite IQ-TREE please use:" << endl << endl
    << "Bui Quang Minh, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf," << endl
    << "Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear (2020)" << endl
    << "IQ-TREE 2: New models and efficient methods for phylogenetic inference" << endl
    << "in the genomic era. Mol. Biol. Evol., in press." << endl
    << "https://doi.org/10.1093/molbev/msaa015" << endl << endl;
    
    bool modelfinder_only = false;
    if (params.model_name.substr(0,4) == "TEST" || params.model_name.substr(0, 2) == "MF" || params.model_name.empty()) {
        out << "To cite ModelFinder please use: " << endl << endl
            << "Subha Kalyaanamoorthy, Bui Quang Minh, Thomas KF Wong, Arndt von Haeseler," << endl
            << "and Lars S Jermiin (2017) ModelFinder: Fast model selection for" << endl
            << "accurate phylogenetic estimates. Nature Methods, 14:587–589." << endl
            << "https://doi.org/10.1038/nmeth.4285" << endl << endl;
        if (params.model_name.find("ONLY") != string::npos || (params.model_name.substr(0,2)=="MF" && params.model_name.substr(0,3)!="MFP"))
            modelfinder_only = true;
    }
    if (posPOMO(params.model_name) != string::npos) {
        out << "For polymorphism-aware models please cite:" << endl << endl
            << "Dominik Schrempf, Bui Quang Minh, Nicola De Maio, Arndt von Haeseler, and Carolin Kosiol" << endl
            << "(2016) Reversible polymorphism-aware phylogenetic models and their application to" << endl
            << "tree inference. J. Theor. Biol., 407:362–370." << endl
            << "https://doi.org/10.1016/j.jtbi.2016.07.042" << endl << endl;
    }

    if (params.site_freq_file || params.tree_freq_file)
    out << "Since you used site-specific frequency model please also cite: " << endl << endl
        << "Huai-Chun Wang, Edward Susko, Bui Quang Minh, and Andrew J. Roger (2018)" << endl
        << "Modeling site heterogeneity with posterior mean site frequency profiles" << endl
        << "accelerates accurate phylogenomic estimation. Syst. Biol., 67:216–235." << endl
        << "https://doi.org/10.1093/sysbio/syx068" << endl << endl;


    if (params.gbo_replicates)
    out << "Since you used ultrafast bootstrap (UFBoot) please also cite: " << endl << endl
        << "Diep Thi Hoang, Olga Chernomor, Arndt von Haeseler, Bui Quang Minh," << endl
        << "and Le Sy Vinh (2018) UFBoot2: Improving the ultrafast bootstrap" << endl
        << "approximation. Mol. Biol. Evol., 35:518–522." << endl
        << "https://doi.org/10.1093/molbev/msx281" << endl << endl;

    if (params.partition_file)
    out << "Since you used partition models please also cite:" << endl << endl
        << "Olga Chernomor, Arndt von Haeseler, and Bui Quang Minh (2016)" << endl
        << "Terrace aware data structure for phylogenomic inference from" << endl
        << "supermatrices. Syst. Biol., 65:997-1008." << endl
        << "https://doi.org/10.1093/sysbio/syw037" << endl << endl;

    if (params.terrace_analysis)
    out << "Since you used terrace analysis please also cite:" << endl << endl
        << "Biczok R, Bozsoky P, Eisenmann P, Ernst J, Ribizel T, Scholz F," << endl
        << "Trefzer A, Weber F, Hamann M, Stamatakis A. (2018)" << endl
        << "Two C++ libraries for counting trees on a phylogenetic" << endl
        << "terrace. Bioinformatics 34:3399–3401." << endl
        << "https://doi.org/10.1093/bioinformatics/bty384" << endl << endl;

    if (params.dating_method == "LSD")
        out << "Since you used least square dating (LSD) please also cite: " << endl << endl
        << "Thu-Hien To, Matthieu Jung, Samantha Lycett, Olivier Gascuel (2016)" << endl
        << "Fast dating using least-squares criteria and algorithms. Syst. Biol. 65:82-97." << endl
        << "https://doi.org/10.1093/sysbio/syv068" << endl << endl;
}

void reportAlignment(ofstream &out, Alignment &alignment, int nremoved_seqs) {
    out << "Input data: " << alignment.getNSeq()+nremoved_seqs << " sequences with "
            << alignment.getNSite() << " ";
    switch (alignment.seq_type) {
    case SEQ_BINARY: out << "binary"; break;
    case SEQ_DNA: out << "nucleotide"; break;
    case SEQ_PROTEIN: out << "amino-acid"; break;
    case SEQ_CODON: out << "codon"; break;
    case SEQ_MORPH: out << "morphological"; break;
    case SEQ_POMO: out << "PoMo"; break;
    default: out << "unknown"; break;
    }
    out << " sites" << endl << "Number of constant sites: "
        << round(alignment.frac_const_sites * alignment.getNSite())
        << " (= " << alignment.frac_const_sites * 100 << "% of all sites)" << endl

        << "Number of invariant (constant or ambiguous constant) sites: "
        << round(alignment.frac_invariant_sites * alignment.getNSite())
        << " (= " << alignment.frac_invariant_sites * 100 << "% of all sites)" << endl

        << "Number of parsimony informative sites: " << alignment.num_informative_sites << endl

        << "Number of distinct site patterns: " << alignment.size() << endl
        << endl;
}

/*
void pruneModelInfo(ModelCheckpoint &model_info, PhyloSuperTree *tree) {
    ModelCheckpoint res_info;
    for (vector<PartitionInfo>::iterator it = tree->part_info.begin(); it != tree->part_info.end(); it++) {
        for (ModelCheckpoint::iterator mit = model_info.begin(); mit != model_info.end(); mit++)
            if (mit->set_name == it->name)
                res_info.push_back(*mit);
    }
    model_info = res_info;

}
*/
void reportModelSelection(ofstream &out, Params &params, ModelCheckpoint *model_info, PhyloTree *tree) {
    out << "Best-fit model according to " << criterionName(params.model_test_criterion) << ": ";
//    ModelCheckpoint::iterator it;
    string best_model;
    PhyloSuperTree *stree = (tree->isSuperTree()) ? ((PhyloSuperTree*)tree) : NULL;
    if (tree->isSuperTree()) {
        SuperAlignment *saln = (SuperAlignment*)stree->aln;
        for (int part = 0; part != stree->size(); part++) {
            if (part != 0)
                out << ",";
            out << saln->partitions[part]->model_name << ":" << saln->partitions[part]->name;
        }
//        string set_name = "";
//        for (it = model_info.begin(); it != model_info.end(); it++) {
//            if (it->set_name != set_name) {
//                if (set_name != "")
//                    out << ",";
//                out << it->name << ":" << it->set_name;
//                set_name = it->set_name;
//            }
//        }
    } else {
//        out << model_info[0].name;
        model_info->getBestModel(best_model);
        out << best_model;
    }

    if (tree->isSuperTree()) {
        out << endl << endl << "List of best-fit models per partition:" << endl << endl;
    } else {
        out << endl << endl << "List of models sorted by "
            << ((params.model_test_criterion == MTC_BIC) ? "BIC" :
                ((params.model_test_criterion == MTC_AIC) ? "AIC" : "AICc"))
            << " scores: " << endl << endl;
    }
    if (tree->isSuperTree())
        out << "  ID  ";
    out << "Model                  LogL         AIC      w-AIC        AICc     w-AICc         BIC      w-BIC" << endl;
    /*
    if (is_partitioned)
        out << "----------";

    out << "----------------------------------------------------------------------------------------" << endl;
    */
    int setid = 1;
    out.precision(3);

    CandidateModelSet models;
    model_info->getOrderedModels(tree, models);
    for (auto it = models.begin(); it != models.end(); it++) {
        if (tree->isSuperTree()) {
            out.width(4);
            out << right << setid << "  ";
            setid++;
        }
        out.width(15);
        out << left << it->getName() << " ";
        out.width(11);
        out << right << it->logl << " ";
        out.width(11);
        out    << it->AIC_score << ((it->AIC_conf) ? " + " : " - ");
        out.unsetf(ios::fixed);
        out.width(8);
        out << it->AIC_weight << " ";
        out.setf(ios::fixed);
        out.width(11);
        out << it->AICc_score << ((it->AICc_conf) ? " + " : " - ");
        out.unsetf(ios::fixed);
        out.width(8);
        out << it->AICc_weight << " ";
        out.setf(ios::fixed);
        out.width(11);
        out << it->BIC_score  << ((it->BIC_conf) ? " + " : " - ");
        out.unsetf(ios::fixed);
        out.width(8);
        out << it->BIC_weight;
        out.setf(ios::fixed);
        out << endl;
    }
    out.precision(4);

    /* TODO
    for (it = model_info.begin(); it != model_info.end(); it++) {
        if (it->AIC_score == DBL_MAX) continue;
        if (it != model_info.begin() && it->set_name != (it-1)->set_name)
            setid++;
        if (is_partitioned && it != model_info.begin() && it->set_name == (it-1)->set_name)
            continue;
        if (is_partitioned) {
            out.width(4);
            out << right << setid << "  ";
        }
        out.width(15);
        out << left << it->name << " ";
        out.width(11);
        out << right << it->logl << " ";
        out.width(11);
        out    << it->AIC_score << ((it->AIC_conf) ? " + " : " - ") << it->AIC_weight << " ";
        out.width(11);
        out << it->AICc_score << ((it->AICc_conf) ? " + " : " - ") << it->AICc_weight << " ";
        out.width(11);
        out << it->BIC_score  << ((it->BIC_conf) ? " + " : " - ") << it->BIC_weight;
        out << endl;
    }
    */
    out << endl;
    out <<  "AIC, w-AIC   : Akaike information criterion scores and weights." << endl
         << "AICc, w-AICc : Corrected AIC scores and weights." << endl
         << "BIC, w-BIC   : Bayesian information criterion scores and weights." << endl << endl

         << "Plus signs denote the 95% confidence sets." << endl
         << "Minus signs denote significant exclusion." <<endl;
    out << endl;
}

void reportModel(ofstream &out, Alignment *aln, ModelSubst *m) {
    int i, j, k;
    ASSERT(aln->num_states == m->num_states);
    double *rate_mat = new double[m->num_states * m->num_states];
    if (!m->isSiteSpecificModel())
        m->getRateMatrix(rate_mat);
    else
        ((ModelSet*)m)->front()->getRateMatrix(rate_mat);

    if (m->num_states <= 4) {
        out << "Rate parameter R:" << endl << endl;

        if (m->num_states > 4)
            out << fixed;
        if (m->isReversible()) {
            for (i = 0, k = 0; i < m->num_states - 1; i++)
                for (j = i + 1; j < m->num_states; j++, k++) {
                    out << "  " << aln->convertStateBackStr(i) << "-" << aln->convertStateBackStr(j) << ": "
                            << rate_mat[k];
                    if (m->num_states <= 4)
                        out << endl;
                    else if (k % 5 == 4)
                        out << endl;
                }

        } else { // non-reversible model
            for (i = 0, k = 0; i < m->num_states; i++)
                for (j = 0; j < m->num_states; j++)
                    if (i != j) {
                        out << "  " << aln->convertStateBackStr(i) << "-" << aln->convertStateBackStr(j)
                                << ": " << rate_mat[k];
                        if (m->num_states <= 4)
                            out << endl;
                        else if (k % 5 == 4)
                            out << endl;
                        k++;
                    }

        }

        //if (tree.aln->num_states > 4)
        out << endl;
        out.unsetf(ios_base::fixed);
    } else if (aln->seq_type == SEQ_PROTEIN && m->getNDim() > 20) {
        ASSERT(m->num_states == 20);
        out << "WARNING: This model has " << m->getNDim() + m->getNDimFreq() << " parameters that may be overfitting. Please use with caution!" << endl << endl;
        double full_mat[400];
        
        out.precision(6);

        if (m->isReversible()) {
            for (i = 0, k = 0; i < m->num_states - 1; i++)
                for (j = i + 1; j < m->num_states; j++, k++) {
                    full_mat[i*m->num_states+j] = rate_mat[k];
                }
            out << "Substitution parameters (lower-diagonal) and state frequencies in PAML format (can be used as input for IQ-TREE): " << endl << endl;
            for (i = 1; i < m->num_states; i++) {
                for (j = 0; j < i; j++)
                    out << " " << full_mat[j*m->num_states+i];
                out << endl;
            }
        } else {
            // non-reversible model
            m->getQMatrix(full_mat);
            out << "Full Q matrix and state frequencies (can be used as input for IQ-TREE): " << endl << endl;
            for (i = 0; i < m->num_states; i++) {
                for (j = 0; j < m->num_states; j++)
                    out << " " << full_mat[i*m->num_states+j];
                out << endl;
            }
        }
        double state_freq[20];
        m->getStateFrequency(state_freq);
        for (i = 0; i < m->num_states; i++)
            out << " " << state_freq[i];
        out << endl << endl;
        out.precision(4);
    }
    
    delete[] rate_mat;

    if (aln->seq_type == SEQ_POMO) {
        m->report(out);
        return;
    }

    out << "State frequencies: ";
    if (m->isSiteSpecificModel())
        out << "(site specific frequencies)" << endl << endl;
    else {
        // 2016-11-03: commented out as this is not correct anymore
//        if (!m->isReversible())
//            out << "(inferred from Q matrix)" << endl;
//        else
            switch (m->getFreqType()) {
            case FREQ_EMPIRICAL:
                out << "(empirical counts from alignment)" << endl;
                break;
            case FREQ_ESTIMATE:
                out << "(estimated with maximum likelihood)" << endl;
                break;
            case FREQ_USER_DEFINED:
                out << ((aln->seq_type == SEQ_PROTEIN) ? "(model)" : "(user-defined)") << endl;
                break;
            case FREQ_EQUAL:
                out << "(equal frequencies)" << endl;
                break;
            default:
                break;
            }
        out << endl;

        if ((m->getFreqType() != FREQ_USER_DEFINED || aln->seq_type == SEQ_DNA) && m->getFreqType() != FREQ_EQUAL) {
            double *state_freqs = new double[m->num_states];
            m->getStateFrequency(state_freqs);
            int ncols=(aln->seq_type == SEQ_CODON) ? 4 : 1;
            for (i = 0; i < m->num_states; i++) {
                out << "  pi(" << aln->convertStateBackStr(i) << ") = " << state_freqs[i];
                if (i % ncols == ncols-1)
                    out << endl;
            }
            delete[] state_freqs;
            out << endl;
        }
        if (m->num_states <= 4 || verbose_mode >= VB_MED) {
            // report Q matrix
            if (verbose_mode >= VB_MED)
                out.precision(6);
            double *q_mat = new double[m->num_states * m->num_states];
            m->getQMatrix(q_mat);

            out << "Rate matrix Q:" << endl << endl;
            for (i = 0, k = 0; i < m->num_states; i++) {
                out << "  " << aln->convertStateBackStr(i);
                for (j = 0; j < m->num_states; j++, k++) {
                    out << "  ";
                    out.width(8);
                    out << q_mat[k];
                }
                out << endl;
            }
            out << endl;
            delete[] q_mat;
        }
    }
}

void reportModel(ofstream &out, PhyloTree &tree) {
//    int i, j, k;
    int i;

    if (tree.getModel()->isMixture() && !tree.getModel()->isPolymorphismAware()) {
        out << "Mixture model of substitution: " << tree.getModelName() << endl;
//        out << "Full name: " << tree.getModelName() << endl;
        ModelSubst *mmodel = tree.getModel();
        out << endl << "  No  Component      Rate    Weight   Parameters" << endl;
        i = 0;
        int nmix = mmodel->getNMixtures();
        for (i = 0; i < nmix; i++) {
            ModelMarkov *m = (ModelMarkov*)mmodel->getMixtureClass(i);
            out.width(4);
            out << right << i+1 << "  ";
            out.width(12);
            out << left << (m)->name << "  ";
            out.width(7);
            out << (m)->total_num_subst << "  ";
            out.width(7);
            out << mmodel->getMixtureWeight(i) << "  " << (m)->getNameParams() << endl;

            if (tree.aln->seq_type == SEQ_POMO) {
                out << endl << "Model for mixture component "  << i+1 << ": " << (m)->name << endl;
                reportModel(out, tree.aln, m);
            }
        }
        if (tree.aln->seq_type != SEQ_POMO && tree.aln->seq_type != SEQ_DNA)
        for (i = 0; i < nmix; i++) {
            ModelMarkov *m = (ModelMarkov*)mmodel->getMixtureClass(i);
            if (m->getFreqType() == FREQ_EQUAL || m->getFreqType() == FREQ_USER_DEFINED)
                continue;
            out << endl << "Model for mixture component "  << i+1 << ": " << (m)->name << endl;
            reportModel(out, tree.aln, m);
        }
        out << endl;
    } else {
        out << "Model of substitution: " << tree.getModelName() << endl << endl;
        reportModel(out, tree.aln, tree.getModel());
    }
}

void reportRate(ostream &out, PhyloTree &tree) {
    RateHeterogeneity *rate_model = tree.getRate();
    out << "Model of rate heterogeneity: " << rate_model->full_name << endl;
    rate_model->writeInfo(out);

    if (rate_model->getNDiscreteRate() > 1 || rate_model->getPInvar() > 0.0) {
        out << endl << " Category  Relative_rate  Proportion" << endl;
        if (rate_model->getPInvar() > 0.0)
            out << "  0         0              " << rate_model->getPInvar()
                    << endl;
        int cats = rate_model->getNDiscreteRate();
        DoubleVector prop;
        if (rate_model->getGammaShape() > 0 || rate_model->getPtnCat(0) < 0) {
            //            prop.resize(cats, (1.0 - rate_model->getPInvar()) / rate_model->getNRate());
            prop.resize(cats);
            for (size_t i = 0; i < cats; i++) {
                prop[i] = rate_model->getProp(i);
            }
        } else {
            prop.resize(cats, 0.0);
            auto   frequencies  = tree.getConvertedSequenceFrequencies();
            size_t num_patterns = tree.aln->getNPattern();
            if (frequencies!=nullptr) {
                for (size_t i = 0; i < num_patterns; i++) {
                    prop[rate_model->getPtnCat(i)] += frequencies[i];
                }
            } else {
                for (size_t i = 0; i < num_patterns; i++) {
                    prop[rate_model->getPtnCat(i)] += tree.aln->at(i).frequency;
                }
            }
            for (size_t i = 0; i < cats; i++) {
                prop[i] /= tree.aln->getNSite();
            }
        }
        for (size_t i = 0; i < cats; i++) {
            out << "  " << i + 1 << "         ";
            out.width(14);
            out << left << rate_model->getRate(i) << " " << prop[i];
            out << endl;
        }
        if (rate_model->isGammaRate()) {
            out << "Relative rates are computed as " << ((rate_model->isGammaRate() == GAMMA_CUT_MEDIAN) ? "MEDIAN" : "MEAN") <<
                " of the portion of the Gamma distribution falling in the category." << endl;
        }
    }
    /*
     if (rate_model->getNDiscreteRate() > 1 || rate_model->isSiteSpecificRate())
     out << endl << "See file " << rate_file << " for site-specific rates and categories" << endl;*/
    out << endl;
}

void reportTree(ofstream &out, Params &params, PhyloTree &tree, double tree_lh, double lh_variance, double main_tree) {
    size_t ssize = tree.getAlnNSite();
    double epsilon = 1.0 / ssize;
    double totalLen = tree.treeLength();
    int df = tree.getModelFactory()->getNParameters(BRLEN_OPTIMIZE);
    double AIC_score, AICc_score, BIC_score;
    computeInformationScores(tree_lh, df, ssize, AIC_score, AICc_score, BIC_score);

    out << "Log-likelihood of the tree: " << fixed << tree_lh;
    if (lh_variance > 0.0)
        out << " (s.e. " << sqrt(lh_variance) << ")";
    out << endl;
    out    << "Unconstrained log-likelihood (without tree): " << tree.aln->computeUnconstrainedLogL() << endl;

    out << "Number of free parameters (#branches + #model parameters): " << df << endl;
//    if (ssize > df) {
//        if (ssize > 40*df)
//            out    << "Akaike information criterion (AIC) score: " << AIC_score << endl;
//        else
//            out << "Corrected Akaike information criterion (AICc) score: " << AICc_score << endl;
//
//        out << "Bayesian information criterion (BIC) score: " << BIC_score << endl;
//    } else
    out    << "Akaike information criterion (AIC) score: " << AIC_score << endl;
    out << "Corrected Akaike information criterion (AICc) score: " << AICc_score << endl;
    out << "Bayesian information criterion (BIC) score: " << BIC_score << endl;

    if (ssize <= df && main_tree) {

        out << endl
            << "**************************** WARNING ****************************" << endl
            << "Number of parameters (K, model parameters and branch lengths): " << df << endl
            << "Sample size (n, alignment length): " << ssize << endl << endl
            << "Given that K>=n, the parameter estimates might be inaccurate." << endl
            << "Thus, phylogenetic estimates should be interpreted with caution." << endl << endl

            << "Ideally, it is desirable that n >> K. When selecting optimal models," << endl
            << "1. use AIC or BIC if n > 40K;" << endl
            << "2. use AICc or BIC if 40K >= n > K;" << endl
            << "3. be extremely cautious if n <= K" << endl << endl

            << "To improve the situation (3), consider the following options:" << endl
            << "  1. Increase the sample size (n)" << endl
            << "  2. Decrease the number of parameters (K) to be estimated. If" << endl
            << "     possible:" << endl
            << "     a. Remove the least important sequences from the alignment" << endl
            << "     b. Specify some of the parameter values for the substitution"<< endl
            << "        model (e.g., the nucleotide or amino acid frequencies)" << endl
            << "     c. Specify some of the parameter values for the rates-across-" << endl
            << "        sites model (e.g., the shape parameter for the discrete" << endl
            << "        Gamma distribution, the proportion of invariable sites, or" << endl
            << "        the rates of change for different rate categories under" << endl
            << "        the FreeRate model)" << endl << endl
            << "Reference:" << endl
            << "Burnham KR, Anderson DR (2002). Model Selection and Multimodel" << endl
            << "Inference: A Practical Information-Theoretic Approach. Springer," << endl
            << "New York." << endl
            << "************************ END OF WARNING ***********************" << endl;
    }
    out << endl;

    if (tree.aln->seq_type == SEQ_POMO) {
        int N = tree.aln->virtual_pop_size;
        out << "NOTE: The branch lengths of PoMo measure mutations and frequency shifts." << endl;
        out << "To compare PoMo branch lengths to DNA substitution models use the tree length" << endl;
        out << "measured in substitutions per site." << endl << endl;
        out << "PoMo branch length = Substitution model branch length * N * N." << endl << endl;
        out << "Total tree length (sum of branch lengths)" << endl;
        out << " - measured in number of mutations and frequency shifts per site: " << totalLen << endl;
        out << " - measured in number of substitutions per site (divided by N^2): " << totalLen / (N * N) << endl;
    }
    else out << "Total tree length (sum of branch lengths): " << totalLen << endl;

    double totalLenInternal = tree.treeLengthInternal(epsilon);
    double totalLenInternalP = totalLenInternal*100.0 / totalLen;
    if (tree.aln->seq_type == SEQ_POMO) {
      int N = tree.aln->virtual_pop_size;
      double totLenIntSub = totalLenInternal/(N * N);
        out << "Sum of internal branch lengths" << endl;
        out << "- measured in mutations and frequency shifts per site: " << totalLenInternal << " (" << totalLenInternalP << "% of tree length)" << endl;
        out << "- measured in substitutions per site: " << totLenIntSub << " (" << totalLenInternalP << "% of tree length)" << endl;
        out << endl;
    }
    else {
        out << "Sum of internal branch lengths: " << totalLenInternal << " (" << totalLenInternalP << "% of tree length)" << endl;
        //    out << "Sum of internal branch lengths divided by total tree length: "
        //            << totalLenInternal / totalLen << endl;
        out << endl;
    }

    if (tree.isMixlen()) {
        DoubleVector lenvec;
        tree.treeLengths(lenvec);
        out << "Class tree lengths: ";
        for (int i = 0; i < lenvec.size(); i++)
            out << " " << lenvec[i];
        out << endl;
    }

    if (params.partition_type == TOPO_UNLINKED) {
        out << "Tree topologies are unlinked across partitions, thus no drawing will be displayed here" << endl;
        out << endl;
        return;
    }
    
    //out << "ZERO BRANCH EPSILON = " << epsilon << endl;
    int zero_internal_branches = tree.countZeroInternalBranches(NULL, NULL, epsilon);
    if (zero_internal_branches > 0) {
        //int zero_internal_branches = tree.countZeroInternalBranches(NULL, NULL, epsilon);
        /*
        out << "WARNING: " << zero_branches
                << " branches of near-zero lengths (<" << epsilon << ") and should be treated with caution!"
                << endl;
        */
        out << "WARNING: " << zero_internal_branches
                << " near-zero internal branches (<" << epsilon << ") should be treated with caution"
                << endl;
        /*
        cout << endl << "WARNING: " << zero_branches
                << " branches of near-zero lengths (<" << epsilon << ") and should be treated with caution!"
                << endl;
        */
        out << "         Such branches are denoted by '**' in the figure below"
                << endl << endl;
    }
    int long_branches = tree.countLongBranches(NULL, NULL, params.max_branch_length-0.2);
    if (long_branches > 0) {
        //stringstream sstr;
        out << "WARNING: " << long_branches << " too long branches (>"
            << params.max_branch_length-0.2 << ") should be treated with caution!" << endl;
        //out << sstr.str();
        //cout << sstr.str();
    }

            //<< "Total tree length: " << tree.treeLength() << endl << endl
    tree.sortTaxa();
    if (tree.rooted)
        out << "NOTE: Tree is ROOTED at virtual root '" << tree.root->name << "'" << endl;
    else
        out << "NOTE: Tree is UNROOTED although outgroup taxon '" << tree.root->name << "' is drawn at root" << endl;

    if (tree.isSuperTree() && params.partition_type == BRLEN_OPTIMIZE)
        out    << "NOTE: Branch lengths are weighted average over all partitions" << endl
            << "      (weighted by the number of sites in the partitions)" << endl;
    if (tree.isMixlen())
        out << "NOTE: Branch lengths are weighted average over heterotachy classes" << endl;

    bool is_codon = tree.aln->seq_type == SEQ_CODON;
    if (tree.isSuperTree()) {
        PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
        is_codon = true;
        for (PhyloSuperTree::iterator sit = stree->begin(); sit != stree->end(); sit++)
            if ((*sit)->aln->seq_type != SEQ_CODON) {
                is_codon = false;
                break;
            }
    }
    if (is_codon)
        out << endl << "NOTE: Branch lengths are interpreted as number of nucleotide substitutions per codon site!"
                << endl << "      Rescale them by 1/3 if you want to have #nt substitutions per nt site" << endl;
    if (main_tree) 
    if (params.aLRT_replicates > 0 || params.gbo_replicates || (params.num_bootstrap_samples && params.compute_ml_tree)) {
        out << "Numbers in parentheses are ";
        if (params.aLRT_replicates > 0) {
            out << "SH-aLRT support (%)";
            if (params.localbp_replicates)
                out << " / local bootstrap support (%)";
        }
        if (params.aLRT_test)
            out << " / parametric aLRT support";
        if (params.aBayes_test)
            out << " / aBayes support";
        if (params.num_bootstrap_samples && params.compute_ml_tree) {
            if (params.aLRT_replicates > 0 || params.aLRT_test || params.aBayes_test)
                out << " /";
            out << " standard " << RESAMPLE_NAME << " support (%)";
        }
        if (params.gbo_replicates) {
            if (params.aLRT_replicates > 0 || params.aLRT_test || params.aBayes_test)
                out << " /";
            out << " ultrafast " << RESAMPLE_NAME << " support (%)";
        }
        out << endl;
    }
    out << endl;

    //tree.setExtendedFigChar();
    tree.setRootNode(params.root, true);
    tree.drawTree(out, WT_BR_SCALE, epsilon);

    out << "Tree in newick format:";
    if (tree.isMixlen())
        out << " (class branch lengths are given in [...] and separated by '/' )";
    if (tree.aln->seq_type == SEQ_POMO)
        out << " (measured in mutations and frequency shifts)";
    out << endl << endl;

    tree.printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
    out << endl << endl;

  if (tree.aln->seq_type == SEQ_POMO) {
    out << "Tree in newick format (measured in substitutions, see above):" << endl;
    out << "WARNING: Only for comparison with substitution models." << endl;
    out << "         These are NOT the branch lengths inferred by PoMo." << endl << endl;
    double len_scale_old = tree.len_scale;
    int N = tree.aln->virtual_pop_size;
    tree.len_scale = 1.0/(N*N);
    tree.printTree(out, WT_BR_SCALE | WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
    tree.len_scale = len_scale_old;
    out << endl << endl;
    }

  tree.setRootNode(params.root, false);

}

void reportCredits(ofstream &out) {
    out << "CREDITS" << endl << "-------" << endl << endl
            << "Some parts of the code were taken from the following packages/libraries:"
            << endl << endl
            << "Schmidt HA, Strimmer K, Vingron M, and von Haeseler A (2002)" << endl
            << "TREE-PUZZLE: maximum likelihood phylogenetic analysis using quartets" << endl
            << "and parallel computing. Bioinformatics, 18(3):502-504." << endl << endl

            //<< "The source code to construct the BIONJ tree were taken from BIONJ software:"
            //<< endl << endl
            << "Gascuel O (1997) BIONJ: an improved version of the NJ algorithm" << endl
            << "based on a simple model of sequence data. Mol. Bio. Evol., 14:685-695." << endl << endl

            //<< "The Nexus file parser was taken from the Nexus Class Library:"
            //<< endl << endl
            << "Paul O. Lewis (2003) NCL: a C++ class library for interpreting data files in" << endl
            << "NEXUS format. Bioinformatics, 19(17):2330-2331." << endl << endl

            << "Mascagni M and Srinivasan A (2000) Algorithm 806: SPRNG: A Scalable Library" << endl
            << "for Pseudorandom Number Generation. ACM Transactions on Mathematical Software," << endl
            << "26: 436-461." << endl << endl

            << "Guennebaud G, Jacob B, et al. (2010) Eigen v3. http://eigen.tuxfamily.org" << endl << endl;
            /*
            << "The Modeltest 3.7 source codes were taken from:" << endl << endl
            << "David Posada and Keith A. Crandall (1998) MODELTEST: testing the model of"
            << endl << "DNA substitution. Bioinformatics, 14(9):817-8." << endl
            */
}

/***********************************************************
 * CREATE REPORT FILE
 ***********************************************************/
extern StringIntMap pllTreeCounter;

void exhaustiveSearchGAMMAInvar(Params &params, IQTree &iqtree);

void searchGAMMAInvarByRestarting(IQTree &iqtree);

void computeLoglFromUserInputGAMMAInvar(Params &params, IQTree &iqtree);

void printOutfilesInfo(Params &params, IQTree &tree) {

    cout << endl << "Analysis results written to: " << endl;
    if (!(params.suppress_output_flags & OUT_IQTREE))
        cout<< "  IQ-TREE report:                " << params.out_prefix << ".iqtree"
            << endl;
    if (params.compute_ml_tree) {
        if (!(params.suppress_output_flags & OUT_TREEFILE)) {
            if (params.model_name.find("ONLY") != string::npos || (params.model_name.substr(0,2)=="MF" && params.model_name.substr(0,3)!="MFP"))
                cout << "  Tree used for ModelFinder:     " << params.out_prefix << ".treefile" << endl;
            else {
                cout << "  Maximum-likelihood tree:       " << params.out_prefix << ".treefile" << endl;
                if (params.partition_type == BRLEN_OPTIMIZE && tree.isSuperTree())
                    cout << "  Partition trees:               " << params.out_prefix << ".parttrees" << endl;
            }
        }
//        if (params.snni && params.write_local_optimal_trees) {
//            cout << "  Locally optimal trees (" << tree.candidateTrees.getNumLocalOptTrees() << "):    " << params.out_prefix << ".suboptimal_trees" << endl;
//        }
    }
    if (params.num_runs > 1)
        cout << "  Trees from independent runs:   " << params.out_prefix << ".runtrees" << endl;

    if (!params.user_file && params.start_tree == STT_BIONJ) {
        cout << "  BIONJ tree:                    " << params.out_prefix << ".bionj"
                << endl;
    }
    if (!params.dist_file) {
        //cout << "  Juke-Cantor distances:    " << params.out_prefix << ".jcdist" << endl;
        if (params.compute_ml_dist)
        cout << "  Likelihood distances:          " << params.out_prefix
                    << ".mldist" << endl;
        if (params.print_conaln)
        cout << "  Concatenated alignment:        " << params.out_prefix
                    << ".conaln" << endl;
    }
    if ((params.model_name.find("TEST") != string::npos || params.model_name.substr(0,2) == "MF") && tree.isSuperTree()) {
        cout << "  Best partitioning scheme:      " << params.out_prefix << ".best_scheme.nex" << endl;
        bool raxml_format_printed = true;

        for (auto it = ((SuperAlignment*)tree.aln)->partitions.begin();
                it != ((SuperAlignment*)tree.aln)->partitions.end(); it++)
            if (!(*it)->aln_file.empty()) {
                raxml_format_printed = false;
                break;
            }
        if (raxml_format_printed)
             cout << "           in RAxML format:      " << params.out_prefix << ".best_scheme" << endl;
    }
    if ((tree.getRate()->getGammaShape() > 0 || params.partition_file) && params.print_site_rate)
        cout << "  Site-specific rates:           " << params.out_prefix << ".rate"
                << endl;

    if ((tree.getRate()->isSiteSpecificRate() || tree.getRate()->getPtnCat(0) >= 0) && params.print_site_rate)
        cout << "  Site-rates by MH model:        " << params.out_prefix << ".rate"
                << endl;

    if (params.print_site_lh)
        cout << "  Site log-likelihoods:          " << params.out_prefix << ".sitelh"
                << endl;

    if (params.print_partition_lh)
        cout << "  Partition log-likelihoods:     " << params.out_prefix << ".partlh"
                << endl;

    if (params.print_site_prob)
        cout << "  Site probability per rate/mix: " << params.out_prefix << ".siteprob"
                << endl;

    if (params.print_ancestral_sequence) {
        cout << "  Ancestral state:               " << params.out_prefix << ".state" << endl;
//        cout << "  Ancestral sequences:           " << params.out_prefix << ".aseq" << endl;
    }

    if (params.write_intermediate_trees)
        cout << "  All intermediate trees:        " << params.out_prefix << ".treels"
                << endl;

    if (params.writeDistImdTrees) {
        tree.intermediateTrees.printTrees(string("ditrees"));
        cout << "  Distinct intermediate trees:   " << params.out_prefix <<  ".ditrees" << endl;
        cout << "  Logl of intermediate trees:    " << params.out_prefix <<  ".ditrees_lh" << endl;
    }

    if (params.gbo_replicates) {
        cout << endl << "Ultrafast " << RESAMPLE_NAME << " approximation results written to:" << endl;
        if (!tree.isSuperTreeUnlinked())
            cout << "  Split support values:          " << params.out_prefix << ".splits.nex" << endl
             << "  Consensus tree:                " << params.out_prefix << ".contree" << endl;
        if (params.print_ufboot_trees)
        cout << "  UFBoot trees:                  " << params.out_prefix << ".ufboot" << endl;

    }

    if (!params.treeset_file.empty()) {
        cout << "  Evaluated user trees:          " << params.out_prefix << ".trees" << endl;

        if (params.print_tree_lh) {
        cout << "  Tree log-likelihoods:          " << params.out_prefix << ".treelh" << endl;
        }
    }
    if (params.lmap_num_quartets >= 0) {
        cout << "  Likelihood mapping plot (SVG): " << params.out_prefix << ".lmap.svg" << endl;
        cout << "  Likelihood mapping plot (EPS): " << params.out_prefix << ".lmap.eps" << endl;
    }
    if (!(params.suppress_output_flags & OUT_LOG))
        cout << "  Screen log file:               " << params.out_prefix << ".log" << endl;
    /*    if (params.model_name == "WHTEST")
     cout <<"  WH-TEST report:           " << params.out_prefix << ".whtest" << endl;*/

    cout << endl;

}

void reportPhyloAnalysis(Params &params, IQTree &tree, ModelCheckpoint &model_info) {
    if (!MPIHelper::getInstance().isMaster()) {
        return;
    }
    if (params.suppress_output_flags & OUT_IQTREE) {
        printOutfilesInfo(params, tree);
        return;
    }
        
    if (params.count_trees) {
        // addon: print #distinct trees
        cout << endl << "NOTE: " << pllTreeCounter.size() << " distinct trees evaluated during whole tree search" << endl;

        IntVector counts;
        for (StringIntMap::iterator i = pllTreeCounter.begin(); i != pllTreeCounter.end(); i++) {
            if (i->second > counts.size())
                counts.resize(i->second+1, 0);
            counts[i->second]++;
        }
        for (IntVector::iterator i2 = counts.begin(); i2 != counts.end(); i2++) {
            if (*i2 != 0) {
                cout << "#Trees occurring " << (i2-counts.begin()) << " times: " << *i2 << endl;
            }
        }
    }
    string outfile = params.out_prefix;

    outfile += ".iqtree";
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(outfile.c_str());
        out << "IQ-TREE " << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR
            << iqtree_VERSION_PATCH << " COVID-edition built " << __DATE__ << endl
                << endl;
        if (params.partition_file)
            out << "Partition file name: " << params.partition_file << endl;
        if (params.aln_file)
            out << "Input file name: " << params.aln_file << endl;

        if (params.user_file)
            out << "User tree file name: " << params.user_file << endl;
        out << "Type of analysis: ";
        bool modelfinder = params.model_name.substr(0,4)=="TEST" || params.model_name.substr(0,2) == "MF" || params.model_name.empty();
        if (modelfinder)
            out << "ModelFinder";
        if (params.compute_ml_tree) {
            if (modelfinder)
                out << " + ";
            out << "tree reconstruction";
        }
        if (params.num_bootstrap_samples > 0) {
            if (params.compute_ml_tree)
                out << " + ";
            out << "non-parametric " << RESAMPLE_NAME << " (" << params.num_bootstrap_samples
                    << " replicates)";
        }
        if (params.gbo_replicates > 0) {
            out << " + ultrafast " << RESAMPLE_NAME << " (" << params.gbo_replicates << " replicates)";
        }
        out << endl;
        out << "Random seed number: " << params.ran_seed << endl << endl;
        out << "REFERENCES" << endl << "----------" << endl << endl;
        reportReferences(params, out);

        out << "SEQUENCE ALIGNMENT" << endl << "------------------" << endl
                << endl;
        if (tree.isSuperTree()) {
      // TODO DS: Changes may be needed here for PoMo.
            out << "Input data: " << tree.aln->getNSeq()+tree.removed_seqs.size() << " taxa with "
                    << tree.aln->getNSite() << " partitions and "
                    << tree.getAlnNSite() << " total sites ("
                    << ((SuperAlignment*)tree.aln)->computeMissingData()*100 << "% missing data)" << endl << endl;

            PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
            int namelen = stree->getMaxPartNameLength();
            int part;
            out.width(max(namelen+6,10));
            out << left << "  ID  Name" << "  Type\tSeq\tSite\tUnique\tInfor\tInvar\tConst" << endl;
            //out << string(namelen+54, '-') << endl;
            part = 0;
            for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
                //out << "FOR PARTITION " << stree->part_info[part].name << ":" << endl << endl;
                //reportAlignment(out, *((*it)->aln));
                out.width(4);
                out << right << part+1 << "  ";
                out.width(max(namelen,4));
                out << left << (*it)->aln->name << "  ";
                switch ((*it)->aln->seq_type) {
                case SEQ_BINARY: out << "BIN"; break;
                case SEQ_CODON: out << "CODON"; break;
                case SEQ_DNA: out << "DNA"; break;
                case SEQ_MORPH: out << "MORPH"; break;
                case SEQ_MULTISTATE: out << "MULTI"; break;
                case SEQ_PROTEIN: out << "AA"; break;
                case SEQ_POMO: out << "POMO"; break;
                case SEQ_UNKNOWN: out << "???"; break;
                }
                out << "\t" << (*it)->aln->getNSeq() << "\t" << (*it)->aln->getNSite()
                    << "\t" << (*it)->aln->getNPattern() << "\t" << (*it)->aln->num_informative_sites
                    << "\t" << (*it)->getAlnNSite() - (*it)->aln->num_variant_sites
                    << "\t" << int((*it)->aln->frac_const_sites*(*it)->getAlnNSite()) << endl;
            }
            out << endl << "Column meanings:" << endl
                << "  Unique: Number of unique site patterns" << endl
                << "  Infor:  Number of parsimony-informative sites" << endl
                << "  Invar:  Number of invariant sites" << endl
                << "  Const:  Number of constant sites (can be subset of invariant sites)" << endl << endl;

        } else
            reportAlignment(out, *(tree.aln), tree.removed_seqs.size());

        out.precision(4);
        out << fixed;

        if (!model_info.empty()) {
            out << "ModelFinder" << endl << "-----------" << endl << endl;
//            if (tree.isSuperTree())
//                pruneModelInfo(model_info, (PhyloSuperTree*)&tree);
            reportModelSelection(out, params, &model_info, &tree);
        }

        out << "SUBSTITUTION PROCESS" << endl << "--------------------" << endl
                << endl;
        if (tree.isSuperTree()) {
            if(params.partition_type == BRLEN_SCALE)
                out << "Edge-linked-proportional partition model with ";
            else if(params.partition_type == BRLEN_FIX)
                out << "Edge-linked-equal partition model with ";
            else if (params.partition_type == BRLEN_OPTIMIZE)
                out << "Edge-unlinked partition model with ";
            else
                out << "Topology-unlinked partition model with ";
            
            if (params.model_joint)
                out << "joint substitution model ";
            else
                out << "separate substitution models ";
            if (params.link_alpha)
                out << "and joint gamma shape";
            else
                out << "and separate rates across sites";
            out << endl << endl;

            PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
            PhyloSuperTree::iterator it;
            int part;
            if(params.partition_type == BRLEN_OPTIMIZE || params.partition_type == TOPO_UNLINKED)
                out << "  ID  Model         TreeLen  Parameters" << endl;
            else
                out << "  ID  Model           Speed  Parameters" << endl;
            //out << "-------------------------------------" << endl;
            for (it = stree->begin(), part = 0; it != stree->end(); it++, part++) {
                out.width(4);
                out << right << (part+1) << "  ";
                out.width(14);
                if(params.partition_type == BRLEN_OPTIMIZE || params.partition_type == TOPO_UNLINKED)
                    out << left << (*it)->getModelName() << " " << (*it)->treeLength() << "  " << (*it)->getModelNameParams() << endl;
                else
                    out << left << (*it)->getModelName() << " " << stree->part_info[part].part_rate  << "  " << (*it)->getModelNameParams() << endl;
            }
            out << endl;
            /*
            for (it = stree->begin(), part = 0; it != stree->end(); it++, part++) {
                reportModel(out, *(*it));
                reportRate(out, *(*it));
            }*/
            PartitionModel *part_model = (PartitionModel*)tree.getModelFactory();
            for (auto itm = part_model->linked_models.begin(); itm != part_model->linked_models.end(); itm++) {
                for (it = stree->begin(); it != stree->end(); it++)
                    if ((*it)->getModel() == itm->second) {
                        out << "Linked model of substitution: " << itm->second->getName() << endl << endl;
                        bool fixed = (*it)->getModel()->fixParameters(false);
                        reportModel(out, (*it)->aln, (*it)->getModel());
                        (*it)->getModel()->fixParameters(fixed);
                        break;
                    }
            }
        } else {
            reportModel(out, tree);
            reportRate(out, tree);
        }
        if (params.lmap_num_quartets >= 0) {
            tree.reportLikelihoodMapping(out);
        }

        /*
        out << "RATE HETEROGENEITY" << endl << "------------------" << endl
                << endl;
        if (tree.isSuperTree()) {
            PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
            int part = 0;
            for (PhyloSuperTree::iterator it = stree->begin();
                    it != stree->end(); it++, part++) {
                out << "FOR PARTITION " << stree->part_info[part].name << ":"
                        << endl << endl;
                reportRate(out, *(*it));
            }
        } else
            reportRate(out, tree);
        */
        // Bootstrap analysis:
        //Display as outgroup: a

        if (params.model_name == "WHTEST") {
            out << "TEST OF MODEL HOMOGENEITY" << endl
                    << "-------------------------" << endl << endl;
            out << "Delta of input data:                 "
                    << params.whtest_delta << endl;
            out << ".95 quantile of Delta distribution:  "
                    << params.whtest_delta_quantile << endl;
            out << "Number of simulations performed:     "
                    << params.whtest_simulations << endl;
            out << "P-value:                             "
                    << params.whtest_p_value << endl;
            if (params.whtest_p_value < 0.05) {
                out
                        << "RESULT: Homogeneity assumption is rejected (p-value cutoff 0.05)"
                        << endl;
            } else {
                out
                        << "RESULT: Homogeneity assumption is NOT rejected (p-value cutoff 0.05)"
                        << endl;
            }
            out << endl << "*** For this result please cite:" << endl << endl;
            out
                    << "G. Weiss and A. von Haeseler (2003) Testing substitution models"
                    << endl
                    << "within a phylogenetic tree. Mol. Biol. Evol, 20(4):572-578"
                    << endl << endl;
        }
        
        if (params.num_runs > 1) {
            out << "MULTIPLE RUNS" << endl
                << "-------------" << endl << endl;
            out << "Run     logL" << endl;
            DoubleVector runLnL;
            tree.getCheckpoint()->getVector("runLnL", runLnL);
            for (int run = 0; run < runLnL.size(); run++)
                out << run+1 << "\t" << fixed << runLnL[run] << endl;
            out << endl;
        }
/*
        out << "TREE SEARCH" << endl << "-----------" << endl << endl
                << "Stopping rule: "
                << ((params.stop_condition == SC_STOP_PREDICT) ? "Yes" : "No")
                << endl << "Number of iterations: "
                << tree.stop_rule.getNumIterations() << endl
                << "Probability of deleting sequences: " << params.p_delete
                << endl << "Number of representative leaves: "
                << params.k_representative << endl
                << "NNI log-likelihood cutoff: " << tree.getNNICutoff() << endl
                << endl;
*/
        if (params.compute_ml_tree) {
            if (params.model_name.find("ONLY") != string::npos || (params.model_name.substr(0,2) == "MF" && params.model_name.substr(0,3) != "MFP")) {
                out << "TREE USED FOR ModelFinder" << endl
                    << "-------------------------" << endl << endl;
            } else if (params.min_iterations == 0) {
                if (params.user_file)
                    out << "USER TREE" << endl
                        << "---------" << endl << endl;
                else
                    out << "STARTING TREE" << endl
                        << "-------------" << endl << endl;
            } else {
                out << "MAXIMUM LIKELIHOOD TREE" << endl
                    << "-----------------------" << endl << endl;
            }

            tree.setRootNode(params.root);

            if (params.gbo_replicates) {
                if (tree.boot_consense_logl > tree.getBestScore() + 0.1 && !tree.isSuperTreeUnlinked()) {
                    out << endl << "**NOTE**: Consensus tree has higher likelihood than ML tree found! Please use consensus tree below." << endl;
                }
            }

            reportTree(out, params, tree, tree.getBestScore(), tree.logl_variance, true);

            if (tree.isSuperTree() && verbose_mode >= VB_MED) {
                PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
//                stree->mapTrees();
//                int empty_branches = stree->countEmptyBranches();
//                if (empty_branches) {
//                    stringstream ss;
//                    ss << empty_branches << " branches in the overall tree with no phylogenetic information due to missing data!";
//                    outWarning(ss.str());
//                }

                int part = 0;
                for (PhyloSuperTree::iterator it = stree->begin();
                        it != stree->end(); it++, part++) {
                    out << "FOR PARTITION " << (*it)->aln->name
                            << ":" << endl << endl;
                    (*it)->setRootNode(params.root);
//                    reportTree(out, params, *(*it), (*it)->computeLikelihood(), (*it)->computeLogLVariance(), false);
                    reportTree(out, params, *(*it), stree->part_info[part].cur_score, 0.0, false);
                }
            }

        }
        /*
         if (params.write_intermediate_trees) {
         out << endl << "CONSENSUS OF INTERMEDIATE TREES" << endl << "-----------------------" << endl << endl
         << "Number of intermediate trees: " << tree.stop_rule.getNumIterations() << endl
         << "Split threshold: " << params.split_threshold << endl
         << "Burn-in: " << params.tree_burnin << endl << endl;
         }*/

        if (params.consensus_type == CT_CONSENSUS_TREE && !tree.isSuperTreeUnlinked()) {
            out << "CONSENSUS TREE" << endl << "--------------" << endl << endl;
            out << "Consensus tree is constructed from "
                    << (params.num_bootstrap_samples ? params.num_bootstrap_samples : params.gbo_replicates)
                    << " " << RESAMPLE_NAME << " trees";
            if (params.gbo_replicates || params.num_bootstrap_samples) {
                out << endl << "Log-likelihood of consensus tree: " << fixed << tree.boot_consense_logl;
            }
            string con_file = params.out_prefix;
            con_file += ".contree";

            // -- Mon Apr 17 21:14:53 BST 2017
            // DONE Minh: merged correctly
            if (params.compute_ml_tree)
                out << endl << "Robinson-Foulds distance between ML tree and consensus tree: "
                    << tree.contree_rfdist << endl;
            // --
            
            out << endl << "Branches with support >"
                    << floor(params.split_threshold * 1000) / 10 << "% are kept";
            if (params.split_threshold == 0.0)
                out << " (extended consensus)";
            if (params.split_threshold == 0.5)
                out << " (majority-rule consensus)";
            if (params.split_threshold >= 0.99)
                out << " (strict consensus)";

            out << endl << "Branch lengths are optimized by maximum likelihood on original alignment" << endl;
            out << "Numbers in parentheses are " << RESAMPLE_NAME << " supports (%)" << endl << endl;

            bool rooted = false;
            MTree contree;
            contree.readTree(con_file.c_str(), rooted);
            contree.drawTree(out, WT_BR_SCALE);
            out << endl << "Consensus tree in newick format: " << endl << endl;
            contree.printTree(out);
            out << endl << endl;
//            tree.freeNode();
//            tree.root = NULL;
//            tree.readTree(con_file.c_str(), rooted);
//            if (removed_seqs.size() > 0) {
//                tree.reinsertIdenticalSeqs(tree.aln, removed_seqs, twin_seqs);
//            }
//            tree.setAlignment(tree.aln);

            // bug fix
//            if ((tree.sse == LK_EIGEN || tree.sse == LK_EIGEN_SSE) && !tree.isBifurcating()) {
//                cout << "NOTE: Changing to old kernel as consensus tree is multifurcating" << endl;
//                tree.changeLikelihoodKernel(LK_SSE);
//            }

//            tree.initializeAllPartialLh();
//            tree.fixNegativeBranch(false);
//            if (tree.isSuperTree())
//                ((PhyloSuperTree*) &tree)->mapTrees();
//            tree.optimizeAllBranches();
//            tree.printTree(con_file.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
//            tree.sortTaxa();
//            tree.drawTree(out, WT_BR_SCALE);
//            out << endl << "Consensus tree in newick format: " << endl << endl;
//            tree.printResultTree(out);
//            out << endl << endl;
        }
#ifdef IQTREE_TERRAPHAST
        if (params.terrace_analysis && params.compute_ml_tree) {
            
            out << "TERRACE ANALYSIS" << endl << "----------------" << endl << endl;
            cout << "Running additional analysis: Phylogenetic Terraces ..."<< endl;
            
            string filename = params.out_prefix;
            filename += ".terrace";

            try
            {
                Terrace terrace(tree, (SuperAlignment*)(tree.aln));

                uint64_t terrace_size = terrace.getSize();

                if (terrace_size == 1) {
                    out << "The tree does not lie on a terrace." << endl;
                } else {
                    out << "The tree lies on a terrace of size ";

                    if (terrace_size == UINT64_MAX) {
                        out << "at least " << terrace_size << " (integer overflow)";
                    } else {
                        out << terrace_size;
                    }

                    out << endl;

                    ofstream terraceout;
                    terraceout.open(filename.c_str());

                    terrace.printTreesCompressed(terraceout);

                    terraceout.close();

                    out << "Terrace trees written (in compressed Newick format) to " << filename << endl;
                }
            }
            catch (std::exception& e)
            {
                out << "ERROR: Terrace analysis using Terraphast failed: " << e.what() << endl << endl;
            }

            out << endl;
            out << "For documentation, see the technical supplement to Biczok et al. (2018)" << endl;
            out << "https://doi.org/10.1093/bioinformatics/bty384";

            out << endl << endl;
            cout<< "Done. Results are written in "<<params.out_prefix<<".iqtree file."<<endl;
        }
#endif
        /* evaluate user trees */
        vector<TreeInfo> info;
        IntVector distinct_trees;
        if (!params.treeset_file.empty()) {
            evaluateTrees(params.treeset_file, params, &tree, info, distinct_trees);
            out.precision(4);
            out.setf(ios_base::fixed);

			out << endl << "USER TREES" << endl << "----------" << endl << endl;
			out << "See " << params.out_prefix << ".trees for trees with branch lengths." << endl << endl;
			if (params.topotest_replicates && info.size() > 1) {
                if (params.do_au_test && params.topotest_replicates < 10000)
                    out << "WARNING: Too few replicates for AU test. At least -zb 10000 for reliable results!" << endl << endl;
                out << "Tree      logL    deltaL  bp-RELL    p-KH     p-SH    ";
				if (params.do_weighted_test)
					out << "p-WKH    p-WSH    ";
                out << "   c-ELW";
                if (params.do_au_test)
                    out << "       p-AU";

                out << endl << "------------------------------------------------------------------";
                if (params.do_weighted_test)
                    out << "------------------";
                if (params.do_au_test)
                    out << "-------";
                out << endl;
			} else {
				out << "Tree      logL    deltaL" << endl;
				out << "-------------------------" << endl;

			}
			double maxL = -DBL_MAX;
			int tid, orig_id;
			for (tid = 0; tid < info.size(); tid++)
				if (info[tid].logl > maxL) maxL = info[tid].logl;
			for (orig_id = 0, tid = 0; orig_id < distinct_trees.size(); orig_id++) {
				out.width(3);
				out << right << orig_id+1 << " ";
				if (distinct_trees[orig_id] >= 0) {
					out << " = tree " << distinct_trees[orig_id]+1 << endl;
					continue;
				}
                out.unsetf(ios::fixed);
				out.precision(10);
				out.width(12);
				out << info[tid].logl << " ";
				out.width(7);
                out.precision(5);
				out << maxL - info[tid].logl;
				if (!params.topotest_replicates || info.size() <= 1) {
					out << endl;
					tid++;
					continue;
				}
				out.precision(3);
				out << "  ";
				out.width(6);
				out << info[tid].rell_bp;
				if (info[tid].rell_confident)
					out << " + ";
				else
					out << " - ";
				out.width(6);
				out << right << info[tid].kh_pvalue;
				if (info[tid].kh_pvalue < 0.05)
					out << " - ";
				else
					out << " + ";
				out.width(6);
				out << right << info[tid].sh_pvalue;
				if (info[tid].sh_pvalue < 0.05)
					out << " - ";
				else
					out << " + ";

				if (params.do_weighted_test) {
					out.width(6);
					out << right << info[tid].wkh_pvalue;
					if (info[tid].wkh_pvalue < 0.05)
						out << " - ";
					else
						out << " + ";
					out.width(6);
					out << right << info[tid].wsh_pvalue;
					if (info[tid].wsh_pvalue < 0.05)
						out << " - ";
					else
						out << " + ";
				}
				out.width(9);
				out << right << info[tid].elw_value;
				if (info[tid].elw_confident)
					out << " + ";
				else
					out << " - ";

                if (params.do_au_test) {
                    out.width(8);
                    out << right << info[tid].au_pvalue;
                    if (info[tid].au_pvalue < 0.05)
                        out << " - ";
                    else
                        out << " + ";
                }
                out.setf(ios::fixed);

                out << endl;
                tid++;
            }
            out << endl;

            if (params.topotest_replicates) {
                out <<  "deltaL  : logL difference from the maximal logl in the set." << endl
                     << "bp-RELL : bootstrap proportion using RELL method (Kishino et al. 1990)." << endl
                     << "p-KH    : p-value of one sided Kishino-Hasegawa test (1989)." << endl
                     << "p-SH    : p-value of Shimodaira-Hasegawa test (2000)." << endl;
                if (params.do_weighted_test) {
                    out << "p-WKH   : p-value of weighted KH test." << endl
                     << "p-WSH   : p-value of weighted SH test." << endl;
                }
                out     << "c-ELW   : Expected Likelihood Weight (Strimmer & Rambaut 2002)." << endl;
                if (params.do_au_test) {
                    out << "p-AU    : p-value of approximately unbiased (AU) test (Shimodaira, 2002)." << endl;
                }
                out  << endl
                     << "Plus signs denote the 95% confidence sets." << endl
                     << "Minus signs denote significant exclusion."  << endl
                     << "All tests performed "
                     << params.topotest_replicates << " resamplings using the RELL method."<<endl;
            }
            out << endl;
        }


        time_t cur_time;
        time(&cur_time);

        char *date_str;
        date_str = ctime(&cur_time);
        out.unsetf(ios_base::fixed);
        out << "TIME STAMP" << endl << "----------" << endl << endl
                << "Date and time: " << date_str << "Total CPU time used: "
                << (double) params.run_time << " seconds (" << convert_time(params.run_time) << ")" << endl
                << "Total wall-clock time used: " << getRealTime() - params.start_real_time
                << " seconds (" << convert_time(getRealTime() - params.start_real_time) << ")" << endl << endl;

        //reportCredits(out); // not needed, now in the manual
        out.close();

    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, outfile);
    }
    
    printOutfilesInfo(params, tree);
}

void checkZeroDist(Alignment *aln, double *dist) {
    size_t ntaxa = aln->getNSeq();
    IntVector checked;
    checked.resize(ntaxa, 0);
    auto minLen = Params::getInstance().min_branch_length;
    for (size_t i = 0; i < ntaxa - 1; ++i) {
        if (checked[i])
            continue;
        string str = "";
        bool first = true;
        auto distRow = dist + i*ntaxa;
        for (size_t j = i + 1; j < ntaxa; ++j)
            if (distRow[j] <= minLen ) {
                if (first)
                    str = "ZERO distance between sequences "
                            + aln->getSeqName(i);
                str += ", " + aln->getSeqName(j);
                checked[j] = 1;
                first = false;
            }
        checked[i] = 1;
        if (str != "") {
            outWarning(str);
        }
    }
}


void printAnalysisInfo(int model_df, IQTree& iqtree, Params& params) {
//    if (!params.raxmllib) {
    cout << "Model of evolution: ";
    if (iqtree.isSuperTree()) {
        cout << iqtree.getModelName() << " (" << model_df << " free parameters)" << endl;
    } else {
        cout << iqtree.getModelName() << " with ";
        switch (iqtree.getModel()->getFreqType()) {
        case FREQ_EQUAL:
            cout << "equal";
            break;
        case FREQ_EMPIRICAL:
            cout << "counted";
            break;
        case FREQ_USER_DEFINED:
            cout << "user-defined";
            break;
        case FREQ_ESTIMATE:
            cout << "optimized";
            break;
        case FREQ_CODON_1x4:
            cout << "counted 1x4";
            break;
        case FREQ_CODON_3x4:
            cout << "counted 3x4";
            break;
        case FREQ_CODON_3x4C:
            cout << "counted 3x4-corrected";
            break;
        case FREQ_DNA_RY:
            cout << "constrained A+G=C+T";
            break;
        case FREQ_DNA_WS:
            cout << "constrained A+T=C+G";
            break;
        case FREQ_DNA_MK:
            cout << "constrained A+C=G+T";
            break;
        case FREQ_DNA_1112:
            cout << "constrained A=C=G";
            break;
        case FREQ_DNA_1121:
            cout << "constrained A=C=T";
            break;
        case FREQ_DNA_1211:
            cout << "constrained A=G=T";
            break;
        case FREQ_DNA_2111:
            cout << "constrained C=G=T";
            break;
        case FREQ_DNA_1122:
            cout << "constrained A=C,G=T";
            break;
        case FREQ_DNA_1212:
            cout << "constrained A=G,C=T";
            break;
        case FREQ_DNA_1221:
            cout << "constrained A=T,C=G";
            break;
        case FREQ_DNA_1123:
            cout << "constrained A=C";
            break;
        case FREQ_DNA_1213:
            cout << "constrained A=G";
            break;
        case FREQ_DNA_1231:
            cout << "constrained A=T";
            break;
        case FREQ_DNA_2113:
            cout << "constrained C=G";
            break;
        case FREQ_DNA_2131:
            cout << "constrained C=T";
            break;
        case FREQ_DNA_2311:
            cout << "constrained G=T";
            break;
        default:
            outError("Wrong specified state frequencies");
        }
        cout << " frequencies (" << model_df << " free parameters)" << endl;
    }
    cout << "Fixed branch lengths: "
            << ((params.fixed_branch_length) ? "Yes" : "No") << endl;

    if (params.min_iterations > 0) {
        cout << "Tree search algorithm: " << (params.snni ? "Stochastic nearest neighbor interchange" : "IQPNNI") << endl;
        cout << "Termination condition: ";
        if (params.stop_condition == SC_REAL_TIME) {
            cout << "after " << params.maxtime << " minutes" << endl;
        } else if (params.stop_condition == SC_UNSUCCESS_ITERATION) {
            cout << "after " << params.unsuccess_iteration << " unsuccessful iterations" << endl;
        } else if (params.stop_condition == SC_FIXED_ITERATION) {
                cout << params.min_iterations << " iterations" << endl;
        } else if(params.stop_condition == SC_WEIBULL) {
                cout << "predicted in [" << params.min_iterations << ","
                        << params.max_iterations << "] (confidence "
                        << params.stop_confidence << ")" << endl;
        } else if (params.stop_condition == SC_BOOTSTRAP_CORRELATION) {
            cout << "min " << params.min_correlation << " correlation coefficient" << endl;
        }

        if (!params.snni) {
            cout << "Number of representative leaves  : " << params.k_representative << endl;
            cout << "Probability of deleting sequences: " << iqtree.getProbDelete() << endl;
            cout << "Number of leaves to be deleted   : " << iqtree.getDelete() << endl;
            cout << "Important quartets assessed on: "
                    << ((params.iqp_assess_quartet == IQP_DISTANCE) ?
                            "Distance" : ((params.iqp_assess_quartet == IQP_PARSIMONY) ? "Parsimony" : "Bootstrap"))
                    << endl;
        }
        cout << "NNI assessed on: " << ((params.nni5) ? "5 branches" : "1 branch") << endl;
    }
    cout << "Phylogenetic likelihood library: " << (params.pll ? "Yes" : "No") << endl;
    if (params.fixed_branch_length != BRLEN_FIX)
        cout << "Branch length optimization method: "
            << ((iqtree.optimize_by_newton) ? "Newton" : "Brent") << endl;
    cout << "Number of Newton-Raphson steps in NNI evaluation and branch length optimization: " << NNI_MAX_NR_STEP
            << " / " << PLL_NEWZPERCYCLE << endl;
    cout << "SSE instructions: "
            << ((iqtree.sse) ? "Yes" : "No") << endl;
    cout << endl;
}

void computeMLDist ( Params& params, IQTree& iqtree
                   , double begin_wallclock_time, double begin_cpu_time) {
    double longest_dist;
    cout << "Computing ML distances based on estimated model parameters..." << endl;
    double *ml_dist = nullptr;
    double *ml_var  = nullptr;
    iqtree.decideDistanceFilePath(params);
    longest_dist = iqtree.computeDist(params, iqtree.aln, ml_dist, ml_var);
    cout << "Computing ML distances took "
        << (getRealTime() - begin_wallclock_time) << " sec (of wall-clock time) "
        << (getCPUTime() - begin_cpu_time) << " sec (of CPU time)" << endl;
    size_t n = iqtree.aln->getNSeq();
    size_t nSquared = n*n;
    if ( iqtree.dist_matrix == nullptr ) {
        iqtree.dist_matrix = ml_dist;
        ml_dist = nullptr;
    } else {
        memmove(iqtree.dist_matrix, ml_dist,
                sizeof(double) * nSquared);
        delete[] ml_dist;
    }
    if ( iqtree.var_matrix == nullptr ) {
        iqtree.var_matrix = ml_var;
        ml_var = nullptr;
    } else {
        memmove(iqtree.var_matrix, ml_var,
                sizeof(double) * nSquared);
        delete[] ml_var;
    }
    if (!params.dist_file)
    {
        iqtree.printDistanceFile();
    }
    double max_genetic_dist = MAX_GENETIC_DIST;
    if (iqtree.aln->seq_type == SEQ_POMO) {
        int N = iqtree.aln->virtual_pop_size;
        max_genetic_dist *= N * N;
    }
    if (longest_dist > max_genetic_dist * 0.99) {
        outWarning("Some pairwise ML distances are too long (saturated)");
    }
}

void computeInitialDist(Params &params, IQTree &iqtree) {
    double longest_dist;
    if (params.dist_file) {
        cout << "Reading distance matrix file " << params.dist_file << " ..." << endl;
    } else if (params.compute_jc_dist) {
        cout << "Computing Jukes-Cantor distances..." << endl;
    } else if (params.compute_obs_dist) {
        cout << "Computing observed distances..." << endl;
    }
    if (params.compute_jc_dist || params.compute_obs_dist || params.partition_file) {
        longest_dist = iqtree.computeDist(params, iqtree.aln, iqtree.dist_matrix, iqtree.var_matrix);
        //if (!params.suppress_zero_distance_warnings) {
        //  checkZeroDist(iqtree.aln, iqtree.dist_matrix);
        //}
        double max_genetic_dist = MAX_GENETIC_DIST;
        if (iqtree.aln->seq_type == SEQ_POMO) {
            int N = iqtree.aln->virtual_pop_size;
            max_genetic_dist *= N * N;
        }
        if (longest_dist > max_genetic_dist * 0.99) {
            outWarning("Some pairwise distances are too long (saturated)");
        }
    }
}

void initializeParams(Params &params, IQTree &iqtree)
{
//    iqtree.setCurScore(-DBL_MAX);
    bool ok_tree = iqtree.root;
    if (iqtree.isSuperTreeUnlinked())
        ok_tree = ((PhyloSuperTree*)&iqtree)->front()->root;
    if (!ok_tree)
    {
        // compute initial tree
        if (!params.compute_ml_tree_only) {
            iqtree.computeInitialTree(params.SSE);
        }
    }
    ASSERT(iqtree.aln);

    if (iqtree.aln->model_name == "WHTEST") {
        if (iqtree.aln->seq_type != SEQ_DNA)
            outError("Weiss & von Haeseler test of model homogeneity only works for DNA");
        iqtree.aln->model_name = "GTR+G";
    }
    if (params.gbo_replicates)
    {
        params.speed_conf = 1.0;
    }

    // TODO: check if necessary
//    if (iqtree.isSuperTree())
//        ((PhyloSuperTree*) &iqtree)->mapTrees();

    // set parameter for the current tree
//    iqtree.setParams(params);
}


void pruneTaxa(Params &params, IQTree &iqtree, double *pattern_lh, NodeVector &pruned_taxa, StrVector &linked_name) {
    int num_low_support;
    double mytime;

    if (params.aLRT_threshold <= 100 && (params.aLRT_replicates > 0 || params.localbp_replicates > 0)) {
        mytime = getCPUTime();
        cout << "Testing tree branches by SH-like aLRT with " << params.aLRT_replicates << " replicates..." << endl;
        iqtree.setRootNode(params.root);
        double curScore =  iqtree.getCurScore();
        iqtree.computePatternLikelihood(pattern_lh, &curScore);
        num_low_support = iqtree.testAllBranches(params.aLRT_threshold, curScore,
                pattern_lh, params.aLRT_replicates, params.localbp_replicates, params.aLRT_test, params.aBayes_test);
        iqtree.printResultTree();
        cout << "  " << getCPUTime() - mytime << " sec." << endl;
        cout << num_low_support << " branches show low support values (<= " << params.aLRT_threshold << "%)" << endl;

        //tree.drawTree(cout);
        cout << "Collapsing stable clades..." << endl;
        iqtree.collapseStableClade(params.aLRT_threshold, pruned_taxa, linked_name, iqtree.dist_matrix);
        cout << pruned_taxa.size() << " taxa were pruned from stable clades" << endl;
    }

    if (!pruned_taxa.empty()) {
        cout << "Pruned alignment contains " << iqtree.aln->getNSeq()
                << " sequences and " << iqtree.aln->getNSite() << " sites and "
                << iqtree.aln->getNPattern() << " patterns" << endl;
        iqtree.initializeAllPartialLh();
        iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.optimizeAllBranches());
        //cout << "Log-likelihood    after reoptimizing model parameters: " << tree.curScore << endl;
        iqtree.optimizeNNI(true, "");
        cout << "Log-likelihood after optimizing partial tree: "
                << iqtree.getCurScore() << endl;
    }

}

void restoreTaxa(IQTree &iqtree, double *saved_dist_mat, NodeVector &pruned_taxa, StrVector &linked_name) {
    if (!pruned_taxa.empty()) {
        cout << "Restoring full tree..." << endl;
        iqtree.restoreStableClade(iqtree.aln, pruned_taxa, linked_name);
        delete[] iqtree.dist_matrix;
        iqtree.dist_matrix = saved_dist_mat;
        iqtree.initializeAllPartialLh();
        iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.optimizeAllBranches());
        //cout << "Log-likelihood    after reoptimizing model parameters: " << tree.curScore << endl;
        pair<int, int> nniInfo;
        nniInfo = iqtree.optimizeNNI(true, "");
        cout << "Log-likelihood    after reoptimizing full tree: " << iqtree.getCurScore() << endl;
        //iqtree.setBestScore(iqtree.getModelFactory()->optimizeParameters(params.fixed_branch_length, true, params.model_eps));

    }
}
void runApproximateBranchLengths(Params &params, IQTree &iqtree) {
    if (!params.fixed_branch_length && params.leastSquareBranch) {
        cout << endl << "Computing Least Square branch lengths..." << endl;
        iqtree.optimizeAllBranchesLS();
        iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.computeLikelihood());
        string filename = params.out_prefix;
        filename += ".lstree";
        iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        cout << "Logl of tree with LS branch lengths: " << iqtree.getCurScore() << endl;
        cout << "Tree with LS branch lengths written to " << filename << endl;
        if (params.print_branch_lengths) {
            if (params.manuel_analytic_approx) {
                cout << "Applying Manuel's analytic approximation.." << endl;
                iqtree.approxAllBranches();
            }
            ofstream out;
            filename = params.out_prefix;
            filename += ".lsbrlen";
            out.open(filename.c_str());
            iqtree.printBranchLengths(out);
            out.close();
            cout << "LS Branch lengths written to " << filename << endl;
        }
        cout << "Total LS tree length: " << iqtree.treeLength() << endl;
    }

    if (params.pars_branch_length) {
        cout << endl << "Computing parsimony branch lengths..." << endl;
        iqtree.fixNegativeBranch(true);
        iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.computeLikelihood());
        string filename = params.out_prefix;
        filename += ".mptree";
        iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        cout << "Logl of tree with MP branch lengths: " << iqtree.getCurScore() << endl;
        cout << "Tree with MP branch lengths written to " << filename << endl;
        if (params.print_branch_lengths) {
            ofstream out;
            filename = params.out_prefix;
            filename += ".mpbrlen";
            out.open(filename.c_str());
            iqtree.printBranchLengths(out);
            out.close();
            cout << "MP Branch lengths written to " << filename << endl;
        }
        cout << "Total MP tree length: " << iqtree.treeLength() << endl;

    }

    if (params.bayes_branch_length) {
        cout << endl << "Computing Bayesian branch lengths..." << endl;
        iqtree.computeAllBayesianBranchLengths();
        iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.computeLikelihood());
        string filename = params.out_prefix;
        filename += ".batree";
        iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        cout << "Logl of tree with Bayesian branch lengths: " << iqtree.getCurScore() << endl;
        cout << "Tree with Bayesian branch lengths written to " << filename << endl;
        if (params.print_branch_lengths) {
            ofstream out;
            filename = params.out_prefix;
            filename += ".babrlen";
            out.open(filename.c_str());
            iqtree.printBranchLengths(out);
            out.close();
            cout << "Bayesian Branch lengths written to " << filename << endl;
        }
        cout << "Total Bayesian tree length: " << iqtree.treeLength() << endl;

    }

}

void printSiteRates(IQTree &iqtree, const char *rate_file, bool bayes) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(rate_file);
        out << "# Site-specific subtitution rates determined by ";
        if (bayes)
            out<< "empirical Bayesian method" << endl;
        else
            out<< "maximum likelihood" << endl;
        out << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" <<  rate_file << "',header=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl;
        if (iqtree.isSuperTree()) {
            out << "#   Part:   Partition ID (1=" << ((PhyloSuperTree*)&iqtree)->front()->aln->name << ", etc)" << endl
            << "#   Site:   Site ID within partition (starting from 1 for each partition)" << endl;
        } else
            out << "#   Site:   Alignment site ID" << endl;
        
        if (bayes)
            out << "#   Rate:   Posterior mean site rate weighted by posterior probability" << endl
                << "#   Cat:    Category with highest posterior (0=invariable, 1=slow, etc)" << endl
                << "#   C_Rate: Corresponding rate of highest category" << endl;
        else
            out << "#   Rate:   Site rate estimated by maximum likelihood" << endl;
        if (iqtree.isSuperTree())
            out << "Part\t";
        out << "Site\tRate";
        if (bayes)
            out << "\tCat\tC_Rate" << endl;
        else
            out << endl;
        iqtree.writeSiteRates(out, bayes);
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, rate_file);
    }
    cout << "Site rates printed to " << rate_file << endl;
}

void printMiscInfo(Params &params, IQTree &iqtree, double *pattern_lh) {
    if (params.print_site_lh && !params.pll) {
        string site_lh_file = params.out_prefix;
        site_lh_file += ".sitelh";
        if (params.print_site_lh == WSL_SITE)
            printSiteLh(site_lh_file.c_str(), &iqtree, pattern_lh);
        else
            printSiteLhCategory(site_lh_file.c_str(), &iqtree, params.print_site_lh);
    }

    if (params.print_partition_lh && !iqtree.isSuperTree()) {
        outWarning("-wpl does not work with non-partition model");
        params.print_partition_lh = false;
    }
    if (params.print_partition_lh && !params.pll) {
        string part_lh_file = (string)params.out_prefix + ".partlh";
        printPartitionLh(part_lh_file.c_str(), &iqtree, pattern_lh);
    }

    if (params.print_site_prob && !params.pll) {
        printSiteProbCategory(((string)params.out_prefix + ".siteprob").c_str(), &iqtree, params.print_site_prob);
    }
    
    if (params.print_ancestral_sequence) {
        printAncestralSequences(params.out_prefix, &iqtree, params.print_ancestral_sequence);
    }
    
    if (params.print_site_state_freq != WSF_NONE && !params.site_freq_file && !params.tree_freq_file) {
        string site_freq_file = params.out_prefix;
        site_freq_file += ".sitesf";
        printSiteStateFreq(site_freq_file.c_str(), &iqtree);
    }

    if (params.print_trees_site_posterior) {
        cout << "Computing mixture posterior probabilities" << endl;
        IntVector pattern_cat;
        int num_mix = iqtree.computePatternCategories(&pattern_cat);
        cout << num_mix << " mixture components are necessary" << endl;
        string site_mix_file = (string)params.out_prefix + ".sitemix";
        ofstream out(site_mix_file.c_str());
        if (!out.is_open())
            outError("File " + site_mix_file + " could not be opened");
        out << "Ptn\tFreq\tNumMix" << endl;
        int ptn;
        for (ptn = 0; ptn < pattern_cat.size(); ptn++)
            out << ptn << "\t" << (int)iqtree.ptn_freq[ptn] << "\t" << pattern_cat[ptn] << endl;
        out.close();
        cout << "Pattern mixtures printed to " << site_mix_file << endl;

        site_mix_file = (string)params.out_prefix + ".sitemixall";
        out.open(site_mix_file.c_str());
        int ncat = iqtree.getRate()->getNRate();
        if (iqtree.getModel()->isMixture() && !iqtree.getModelFactory()->fused_mix_rate)
            ncat = iqtree.getModel()->getNMixtures();
        out << "Ptn\tFreq\tNumMix\tCat" << endl;

        int c;
        for (ptn = 0; ptn < iqtree.ptn_cat_mask.size(); ptn++) {
            int num_cat = popcount_lauradoux((unsigned*)&iqtree.ptn_cat_mask[ptn], 2);
            out << ptn << "\t" << (int)iqtree.ptn_freq[ptn] << "\t" << num_cat << "\t";
            for (c = 0; c < ncat; c++)
                if (iqtree.ptn_cat_mask[ptn] & ((uint64_t)1<<c))
                    out << "1";
                else
                    out << "0";
            out << endl;
        }
        out.close();
    }

    if (params.print_branch_lengths) {
        if (params.manuel_analytic_approx) {
            cout << "Applying Manuel's analytic approximation.." << endl;
            iqtree.approxAllBranches();
        }
        string brlen_file = params.out_prefix;
        brlen_file += ".brlen";
        ofstream out;
        out.open(brlen_file.c_str());
        iqtree.printBranchLengths(out);
        out.close();
        cout << "Branch lengths written to " << brlen_file << endl;
    }

    if (params.write_branches) {
        string filename = string(params.out_prefix) + ".branches.csv";
        ofstream out;
        out.open(filename.c_str());
        iqtree.writeBranches(out);
        out.close();
        cout << "Branch lengths written to " << filename << endl;
    }
    
    if (params.print_conaln && iqtree.isSuperTree()) {
        string str = params.out_prefix;
        str = params.out_prefix;
        str += ".conaln";
        iqtree.aln->printAlignment(params.aln_output_format, str.c_str());
    }
    
    if (params.print_partition_info && iqtree.isSuperTree()) {
        ASSERT(params.print_conaln);
        string aln_file = (string)params.out_prefix + ".conaln";
        string partition_info = params.out_prefix;
        partition_info += ".partinfo.nex";
        ((SuperAlignment*)(iqtree.aln))->printPartition(partition_info.c_str(), aln_file.c_str());
        partition_info = (string)params.out_prefix + ".partitions";
        ((SuperAlignment*)(iqtree.aln))->printPartitionRaxml(partition_info.c_str());
    }

    if (params.mvh_site_rate) {
        RateMeyerHaeseler *rate_mvh = new RateMeyerHaeseler(params.rate_file,
                &iqtree, params.rate_mh_type);
        cout << endl << "Computing site-specific rates by "
                << rate_mvh->full_name << "..." << endl;
        rate_mvh->runIterativeProc(params, iqtree);
        cout << endl << "BEST SCORE FOUND : " << iqtree.getBestScore()<< endl;
        string mhrate_file = params.out_prefix;
        mhrate_file += ".mhrate";
        try {
            ofstream out;
            out.exceptions(ios::failbit | ios::badbit);
            out.open(mhrate_file.c_str());
            iqtree.writeSiteRates(out, true);
            out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, mhrate_file);
        }

        if (params.print_site_lh) {
            string site_lh_file = params.out_prefix;
            site_lh_file += ".mhsitelh";
            printSiteLh(site_lh_file.c_str(), &iqtree);
        }
    }

    if (params.print_site_rate & 1) {
        string rate_file = params.out_prefix;
        rate_file += ".rate";
        printSiteRates(iqtree, rate_file.c_str(), true);
    }

    if (params.print_site_rate & 2) {
        string rate_file = params.out_prefix;
        rate_file += ".mlrate";
        printSiteRates(iqtree, rate_file.c_str(), false);
    }

    if (params.fixed_branch_length == BRLEN_SCALE) {
        string filename = (string)params.out_prefix + ".blscale";
        iqtree.printTreeLengthScaling(filename.c_str());
        cout << "Scaled tree length and model parameters printed to " << filename << endl;
    }

}

void printFinalSearchInfo(Params &params, IQTree &iqtree, double search_cpu_time, double search_real_time) {
    cout << "Total tree length: " << iqtree.treeLength() << endl;

    if (iqtree.isSuperTree() && verbose_mode >= VB_MAX) {
        PhyloSuperTree *stree = (PhyloSuperTree*) &iqtree;
        cout << stree->evalNNIs << " NNIs evaluated from " << stree->totalNNIs << " all possible NNIs ( " <<
                (int)(((stree->evalNNIs+1.0)/(stree->totalNNIs+1.0))*100.0) << " %)" << endl;
        cout<<"Details for subtrees:"<<endl;
        int part = 0;
        for(auto it = stree->begin(); it != stree->end(); it++,part++){
            cout << part+1 << ". " << (*it)->aln->name << ": " << stree->part_info[part].evalNNIs<<" ( "
                << (int)(((stree->part_info[part].evalNNIs+1.0)/((stree->totalNNIs+1.0) / stree->size()))*100.0)
                << " %)" << endl;
        }
    }

    params.run_time = (getCPUTime() - params.startCPUTime);
    cout << endl;
    cout << "Total number of iterations: " << iqtree.stop_rule.getCurIt() << endl;
//    cout << "Total number of partial likelihood vector computations: " << iqtree.num_partial_lh_computations << endl;
    cout << "CPU time used for tree search: " << search_cpu_time
            << " sec (" << convert_time(search_cpu_time) << ")" << endl;
    cout << "Wall-clock time used for tree search: " << search_real_time
            << " sec (" << convert_time(search_real_time) << ")" << endl;
    cout << "Total CPU time used: " << (double) params.run_time << " sec ("
            << convert_time((double) params.run_time) << ")" << endl;
    cout << "Total wall-clock time used: "
            << getRealTime() - params.start_real_time << " sec ("
            << convert_time(getRealTime() - params.start_real_time) << ")" << endl;

}

void printTrees(vector<string> trees, Params &params, string suffix) {
    ofstream treesOut((string(params.out_prefix) + suffix).c_str(),
            ofstream::out);
    for (vector<string>::iterator it = trees.begin(); it != trees.end(); it++) {
        treesOut << (*it);
        treesOut << endl;
    }
    treesOut.close();
}

/************************************************************
 *  MAIN TREE RECONSTRUCTION
 ***********************************************************/

void startTreeReconstruction(Params &params, IQTree* &iqtree, ModelCheckpoint &model_info) {
    if (params.root) {
        StrVector outgroup_names;
        convert_string_vec(params.root, outgroup_names);
        for (auto it = outgroup_names.begin(); it != outgroup_names.end(); it++)
            if (iqtree->aln->getSeqID(*it) < 0)
                outError("Alignment does not have specified outgroup taxon ", *it);
    }

//    if (params.count_trees && pllTreeCounter == NULL)
//        pllTreeCounter = new StringIntMap;

    // Temporary fix since PLL only supports DNA/Protein: switch to IQ-TREE parsimony kernel
    if (params.start_tree == STT_PLL_PARSIMONY) {
        if (iqtree->isSuperTreeUnlinked()) {
            params.start_tree = STT_PARSIMONY;
        } else if (iqtree->isSuperTree()) {
            PhyloSuperTree *stree = (PhyloSuperTree*)iqtree;
            for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++)
                if ((*it)->aln->seq_type != SEQ_DNA && (*it)->aln->seq_type != SEQ_PROTEIN)
                    params.start_tree = STT_PARSIMONY;
        } else if (iqtree->aln->seq_type != SEQ_DNA && iqtree->aln->seq_type != SEQ_PROTEIN)
            params.start_tree = STT_PARSIMONY;
    }

    /***************** Initialization for PLL and sNNI ******************/
    if (params.start_tree == STT_PLL_PARSIMONY || params.start_tree == STT_RANDOM_TREE || params.pll) {
        /* Initialized all data structure for PLL*/
        iqtree->initializePLL(params);
    }
    
    /********************* Compute pairwise distances *******************/
    if (params.start_tree == STT_BIONJ || params.iqp || params.leastSquareBranch) {
        computeInitialDist(params, *iqtree);
    }
    
    /******************** Pass the parameter object params to IQTree *******************/
    iqtree->setParams(&params);

    /*************** SET UP PARAMETERS and model testing ****************/

       // FOR TUNG: swapping the order cause bug for -m TESTLINK
//    iqtree.initSettings(params);

    runModelFinder(params, *iqtree, model_info);
}
        
/**
 optimize branch lengths of consensus tree
 */
void optimizeConTree(Params &params, IQTree *tree) {
    string contree_file = string(params.out_prefix) + ".contree";
    
    DoubleVector rfdist;
    tree->computeRFDist(contree_file.c_str(), rfdist);
    tree->contree_rfdist = (int)rfdist[0];
    
    tree->readTreeFile(contree_file);
    
    tree->initializeAllPartialLh();
    tree->fixNegativeBranch(false);
    
    tree->boot_consense_logl = tree->optimizeAllBranches();
    cout << "Log-likelihood of consensus tree: " << tree->boot_consense_logl << endl;
    tree->setRootNode(params.root);
    tree->insertTaxa(tree->removed_seqs, tree->twin_seqs);
    tree->printTree(contree_file.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    string contree = tree->getTreeString();
    tree->getCheckpoint()->put("contree", contree);
}

void handleGammaInvariantOptions(Params &params, IQTree &iqtree);

void handleQuartetLikelihoodMapping(Params &params, IQTree &iqtree);

void runTreeReconstruction(Params &params, IQTree* &iqtree) {

    //    string dist_file;
    params.startCPUTime = getCPUTime();
    params.start_real_time = getRealTime();
    
    int absent_states = 0;
    if (iqtree->isSuperTree()) {
        PhyloSuperTree *stree = (PhyloSuperTree*)iqtree;
        for (auto i = stree->begin(); i != stree->end(); i++)
            absent_states += (*i)->aln->checkAbsentStates("partition " + (*i)->aln->name);
    } else {
        absent_states = iqtree->aln->checkAbsentStates("alignment");
    }
    if (absent_states > 0) {
        cout << "NOTE: " << absent_states << " states (see above) are not present and thus removed from Markov process to prevent numerical problems" << endl;
    }
    
    // Make sure that no partial likelihood of IQ-TREE is initialized when PLL is used to save memory
    if (params.pll) {
        iqtree->deleteAllPartialLh();
    }
    
    /***************** Initialization for PLL and sNNI ******************/
    if ((params.start_tree == STT_PLL_PARSIMONY || params.start_tree == STT_RANDOM_TREE || params.pll) && !iqtree->isInitializedPLL()) {
        /* Initialized all data structure for PLL*/
        iqtree->initializePLL(params);
    }
    
    
    /********************* Compute pairwise distances *******************/
    if ((params.start_tree == STT_BIONJ || params.iqp || params.leastSquareBranch) && !iqtree->root) {
        computeInitialDist(params, *iqtree);
    }
    
    /******************** Pass the parameter object params to IQTree *******************/
    iqtree->setParams(&params);
    
    ModelsBlock *models_block = readModelsDefinition(params);

    initializeParams(params, *iqtree);

    if (posRateHeterotachy(iqtree->aln->model_name) != string::npos && !iqtree->isMixlen()) {
        // create a new instance
        IQTree* iqtree_new = new PhyloTreeMixlen(iqtree->aln, 0);
        iqtree_new->setCheckpoint(iqtree->getCheckpoint());
        if (!iqtree->constraintTree.empty())
            iqtree_new->constraintTree.readConstraint(iqtree->constraintTree);
        iqtree_new->removed_seqs = iqtree->removed_seqs;
        iqtree_new->twin_seqs = iqtree->twin_seqs;
        if (params.start_tree == STT_PLL_PARSIMONY || params.start_tree == STT_RANDOM_TREE || params.pll) {
            /* Initialized all data structure for PLL*/
            iqtree_new->initializePLL(params);
        }
        iqtree_new->setParams(&params);
        iqtree_new->copyPhyloTree(iqtree, false);

        // replace iqtree object
        delete iqtree;
        iqtree = iqtree_new;
    }

    if (!params.compute_ml_tree_only) {
        iqtree->setRootNode(params.root);
    }

    iqtree->restoreCheckpoint();

    if (params.online_bootstrap && params.gbo_replicates > 0) {
        cout << "Generating " << params.gbo_replicates << " samples for ultrafast "
        << RESAMPLE_NAME << " (seed: " << params.ran_seed << ")..." << endl;
    }

    iqtree->initSettings(params);

    /*********************** INITIAL MODEL OPTIMIZATION *****************/


    if (!iqtree->getModelFactory()) {
        iqtree->initializeModel(params, iqtree->aln->model_name, models_block);
    }
    if (iqtree->getRate()->isHeterotachy() && !iqtree->isMixlen()) {
        ASSERT(0 && "Heterotachy tree not properly created");
    }
//    iqtree.restoreCheckpoint();

    delete models_block;

    // UpperBounds analysis. Here, to analyse the initial tree without any tree search or optimization
    /*
    if (params.upper_bound) {
        iqtree.setCurScore(iqtree.computeLikelihood());
        cout<<iqtree.getCurScore()<<endl;
        UpperBounds(&params, iqtree.aln, iqtree);
        exit(0);
    }
    */

    // degree of freedom
    cout << endl;
    if (verbose_mode >= VB_MED) {
        cout << "ML-TREE SEARCH START WITH THE FOLLOWING PARAMETERS:" << endl;
        int model_df = iqtree->getModelFactory()->getNParameters(BRLEN_OPTIMIZE);
        printAnalysisInfo(model_df, *iqtree, params);
    }

    if (!params.pll) {
        uint64_t total_mem = getMemorySize();
        if (params.lh_mem_save == LM_MEM_SAVE && params.max_mem_size > total_mem)
            params.max_mem_size = total_mem;

        uint64_t mem_required = iqtree->getMemoryRequired();

        if (mem_required >= total_mem*0.95 && !iqtree->isSuperTree()) {
            // switch to memory saving mode
            if (params.lh_mem_save != LM_MEM_SAVE) {
                params.max_mem_size = (total_mem*0.95)/mem_required;
                params.lh_mem_save = LM_MEM_SAVE;
                mem_required = iqtree->getMemoryRequired();
                cout << "NOTE: Switching to memory saving mode using " << (mem_required / 1073741824.0) << " GB ("
                    <<  (mem_required*100/total_mem) << "% of normal mode)" << endl;
                cout << "NOTE: Use -mem option if you want to restrict RAM usage further" << endl;
            }
            if (mem_required >= total_mem) {
                params.lh_mem_save = LM_MEM_SAVE;
                params.max_mem_size = 0.0;
                mem_required = iqtree->getMemoryRequired();
            }
        }
        if (mem_required >= total_mem) {
            cerr << "ERROR: Your RAM is below minimum requirement of " << (mem_required / 1073741824.0) << " GB RAM" << endl;
            outError("Memory saving mode cannot work, switch to another computer!!!");
        }

//#if defined __APPLE__ || defined __MACH__
        cout << "NOTE: " << (mem_required / 1048576) << " MB RAM (" << (mem_required / 1073741824) << " GB) is required!" << endl;
//#else
//        cout << "NOTE: " << ((double) mem_size / 1000.0) / 1000 << " MB RAM is required!" << endl;
//#endif
        if (params.memCheck)
            exit(0);
#ifdef BINARY32
        if (mem_required >= 2000000000) {
            outError("Memory required exceeds 2GB limit of 32-bit executable");
        }
#endif
        int max_procs = countPhysicalCPUCores();
        if (mem_required * max_procs > total_mem * iqtree->num_threads && iqtree->num_threads > 0) {
            outWarning("Memory required per CPU-core (" + convertDoubleToString((double)mem_required/iqtree->num_threads/1024/1024/1024)+
            " GB) is higher than your computer RAM per CPU-core ("+convertIntToString(total_mem/max_procs/1024/1024/1024)+
            " GB), thus multiple runs may exceed RAM!");
        }
    }

    bool   finishedInitTree = false;
    double initEpsilon = params.min_iterations == 0 ? params.modelEps : (params.modelEps*10);
    string initTree;
    iqtree->prepareToComputeDistances();
    //None of his will work until there are actually taxa tree
    //(we cannot do it until we *have* that).
    if (!params.compute_ml_tree_only) {
        iqtree->ensureNumberOfThreadsIsSet(&params, false);
        iqtree->initializeAllPartialLh();
        handleGammaInvariantOptions(params, *iqtree);
        
        // Optimize model parameters and branch lengths using ML for the initial tree
        iqtree->clearAllPartialLH();
        initTree = iqtree->ensureModelParametersAreSet(initEpsilon);
        handleQuartetLikelihoodMapping(params, *iqtree);
        finishedInitTree = iqtree->getCheckpoint()->getBool("finishedInitTree");
        
        // now overwrite with random tree
        if (params.start_tree == STT_RANDOM_TREE && !finishedInitTree) {
            cout << "Generate random initial Yule-Harding tree..." << endl;
            iqtree->generateRandomTree(YULE_HARDING);
            iqtree->wrapperFixNegativeBranch(true);
            iqtree->initializeAllPartialLh();
            initTree = iqtree->optimizeBranches(params.brlen_num_traversal);
            cout << "Log-likelihood of random tree: " << iqtree->getCurScore() << endl;
        }
        
        /****************** NOW PERFORM MAXIMUM LIKELIHOOD TREE RECONSTRUCTION ******************/
        
        // Update best tree
        if (!finishedInitTree) {
            iqtree->addTreeToCandidateSet(initTree, iqtree->getCurScore(), false, MPIHelper::getInstance().getProcessID());
            iqtree->printResultTree();
            iqtree->intermediateTrees.update(iqtree->getTreeString(), iqtree->getCurScore());
            if (iqtree->isSuperTreeUnlinked()) {
                PhyloSuperTree* stree = (PhyloSuperTree*)iqtree;
                for (auto it = stree->begin(); it != stree->end(); it++)
                    ((IQTree*)(*it))->addTreeToCandidateSet((*it)->getTreeString(),
                                                            (*it)->getCurScore(), false, MPIHelper::getInstance().getProcessID());
            }
        }
        
        if (params.min_iterations && !iqtree->isBifurcating()) {
            outError("Tree search does not work with initial multifurcating tree. Please specify `-n 0` to avoid this.");
        }
        
        // Compute maximum likelihood distance
        // ML distance is only needed for NNI/IQP
        
        if ((params.min_iterations <= 1 || params.numInitTrees <= 1) && params.start_tree != STT_BIONJ) {
            params.compute_ml_dist = false;
        }
        if ((params.user_file || params.start_tree == STT_RANDOM_TREE) && params.snni && !params.iqp) {
            params.compute_ml_dist = false;
        }
        if (params.constraint_tree_file) {
            params.compute_ml_dist = false;
        }
        if (iqtree->isSuperTreeUnlinked()) {
            params.compute_ml_dist = false;
        }
        std::string distFilePath =  iqtree->getDistanceFileWritten();
        if (!distFilePath.empty()) {
            cout << "Wrote distance file to... "
                << iqtree->getDistanceFileWritten() << endl;
        }
    }
    bool wantMLDistances = MPIHelper::getInstance().isMaster() && !iqtree->getCheckpoint()->getBool("finishedCandidateSet");
    if (wantMLDistances) {
        wantMLDistances = !finishedInitTree && ((!params.dist_file && params.compute_ml_dist) || params.leastSquareBranch);
    }
        
    //Compute ML distances, and generate BIONJ tree from those
    if (wantMLDistances || params.compute_ml_tree_only) {
        computeMLDist(params, *iqtree, getRealTime(), getCPUTime());
        bool wasMLDistanceWrittenToFile = false;
        if (!params.user_file) {
            if (params.start_tree != STT_RANDOM_TREE) {
                if (!params.compute_ml_tree_only) {
                    iqtree->resetCurScore();
                }
                double start_bionj = getRealTime();
                iqtree->computeBioNJ(params);
                if (verbose_mode >= VB_MED) {
                    cout << "Wall-clock time spent creating initial tree was "
                    << getRealTime() - start_bionj << " seconds" << endl;
                }
                wasMLDistanceWrittenToFile  = !params.dist_file;
                if (params.compute_ml_tree_only) {
                    iqtree->initializeAllPartialPars();
                    iqtree->clearAllPartialLH();
                    iqtree->fixNegativeBranch(iqtree->isSuperTree());
                    //Because wrapperFixNegativeBranch resets the score,
                    //and that complains if the number of threads isn't set.
                    //(but... if you set the number of threads first, it
                    //complains about the negative path lengths).
                    
                    //This fails if there are any lengths <=0 (so it has to
                    //go after the fix-up for negative branch lengths).
                    iqtree->ensureNumberOfThreadsIsSet(&params, false);
                    iqtree->initializeAllPartialLh();
                    handleGammaInvariantOptions(params, *iqtree);
                    initTree = iqtree->ensureModelParametersAreSet(initEpsilon);
                    handleQuartetLikelihoodMapping(params, *iqtree);
                } else {
                    iqtree->wrapperFixNegativeBranch(iqtree->isSuperTree());
                    iqtree->initializeAllPartialLh();
                    if (params.start_tree == STT_BIONJ) {
                        initTree = iqtree->optimizeModelParameters(params.min_iterations==0, initEpsilon);
                    } else {
                        initTree = iqtree->optimizeBranches();
                    }
                }
                cout << "Log-likelihood of " << params.start_tree_subtype_name
                    << " tree: " << iqtree->getCurScore() << endl;
                iqtree->candidateTrees.update(initTree, iqtree->getCurScore());
            }
        }
        if (!wasMLDistanceWrittenToFile && !params.dist_file) {
            double write_begin_time = getRealTime();
            iqtree->printDistanceFile();
            if (verbose_mode >= VB_MED) {
                #ifdef _OPENMP
                #pragma omp critical (io)
                #endif
                cout << "Time taken to write distance file ( " << iqtree->dist_file << ") : "
                << getRealTime() - write_begin_time << " seconds " << endl;
            }
        }
    }
    //iqtree->saveCheckpoint();

    double cputime_search_start = getCPUTime();
    double realtime_search_start = getRealTime();

    if (params.leastSquareNNI) {
        iqtree->computeSubtreeDists();
    }
    if (params.model_name == "WHTEST") {
        cout << endl << "Testing model homogeneity by Weiss & von Haeseler (2003)..." << endl;
        WHTest(params, *iqtree);
    }
    NodeVector pruned_taxa;
    StrVector linked_name;
    double *saved_dist_mat = iqtree->dist_matrix;
    double *pattern_lh = new double[iqtree->getAlnNPattern()];

    // prune stable taxa
    iqtree->doneComputingDistances();
    pruneTaxa(params, *iqtree, pattern_lh, pruned_taxa, linked_name);

    /***************************************** DO STOCHASTIC TREE SEARCH *******************************************/
    if (params.min_iterations > 0 && !params.tree_spr) {
        iqtree->prepareToComputeDistances();
        iqtree->doTreeSearch();
        iqtree->doneComputingDistances();
        iqtree->setAlignment(iqtree->aln);
    } else {
        iqtree->candidateTrees.saveCheckpoint();
        /* do SPR with likelihood function */
        if (params.tree_spr) {
            //tree.optimizeSPRBranches();
            cout << "Doing SPR Search" << endl;
            cout << "Start tree.optimizeSPR()" << endl;
            double spr_score = iqtree->optimizeSPR();
            cout << "Finish tree.optimizeSPR()" << endl;
            //PhyloNode* nextToRoot = tree.getRoot()->firstNeighbor()->getNode();
            //double spr_score = tree.optimizeSPR(tree.curScore, nextToRoot);
            if (spr_score <= iqtree->getCurScore()) {
                cout << "SPR search did not found any better tree" << endl;
            }
        }
    }
    // restore pruned taxa
    restoreTaxa(*iqtree, saved_dist_mat, pruned_taxa, linked_name);

    double search_cpu_time = getCPUTime() - cputime_search_start;
    double search_real_time = getRealTime() - realtime_search_start;

    if (!MPIHelper::getInstance().isMaster()) {
        delete[] pattern_lh;
        return;
    }
    if (params.snni && params.min_iterations && verbose_mode >= VB_MED) {
        cout << "Log-likelihoods of " << params.popSize << " best candidate trees: " << endl;
        iqtree->printBestScores();
        cout << endl;
    }
    if (!params.final_model_opt) {
        iqtree->setCurScore(iqtree->computeLikelihood());
    } else if (params.min_iterations) {
        iqtree->readTreeString(iqtree->getBestTrees()[0]);
        iqtree->initializeAllPartialLh();
        iqtree->clearAllPartialLH();
        cout << "--------------------------------------------------------------------" << endl;
        cout << "|                    FINALIZING TREE SEARCH                        |" << endl;
        cout << "--------------------------------------------------------------------" << endl;

        if (iqtree->getCheckpoint()->getBool("finishedModelFinal")) {
            iqtree->setCurScore(iqtree->computeLikelihood());
            cout << "CHECKPOINT: Final model parameters restored" << endl;
        } else {
            cout << "Performs final model parameters optimization" << endl;
            string tree;
            Params::getInstance().fixStableSplits = false;
            Params::getInstance().tabu = false;
            tree = iqtree->optimizeModelParameters(true);
            iqtree->addTreeToCandidateSet(tree, iqtree->getCurScore(), false, MPIHelper::getInstance().getProcessID());
            iqtree->getCheckpoint()->putBool("finishedModelFinal", true);
            iqtree->saveCheckpoint();
        }
    }

    if (iqtree->isSuperTree()) {
        ((PhyloSuperTree*) iqtree)->computeBranchLengths();
        ((PhyloSuperTree*) iqtree)->printBestPartitionParams((string(params.out_prefix) + ".best_model.nex").c_str());
    }

    cout << "BEST SCORE FOUND : " << iqtree->getCurScore() << endl;

    if (params.write_candidate_trees) {
        printTrees(iqtree->getBestTrees(), params, ".imd_trees");
    }

    if (params.pll)
        iqtree->inputModelPLL2IQTree();

    /* root the tree at the first sequence */
    // BQM: WHY SETTING THIS ROOT NODE????
//    iqtree->root = iqtree->findLeafName(iqtree->aln->getSeqName(0));
//    assert(iqtree->root);
    iqtree->setRootNode(params.root);


    if (!params.pll) {
        iqtree->computeLikelihood(pattern_lh);
        // compute logl variance
        iqtree->logl_variance = iqtree->computeLogLVariance();
    }

    printMiscInfo(params, *iqtree, pattern_lh);

    if (params.root_test) {
        cout << "Testing root positions..." << endl;
        iqtree->testRootPosition(true, params.loglh_epsilon);
    }
    
    /****** perform SH-aLRT test ******************/
    if ((params.aLRT_replicates > 0 || params.localbp_replicates > 0 || params.aLRT_test || params.aBayes_test) && !params.pll) {
        double mytime = getRealTime();
        params.aLRT_replicates = max(params.aLRT_replicates, params.localbp_replicates);
        cout << endl;
        if (params.aLRT_replicates > 0)
            cout << "Testing tree branches by SH-like aLRT with "
                << params.aLRT_replicates << " replicates..." << endl;
        if (params.localbp_replicates)
            cout << "Testing tree branches by local-BP test with " << params.localbp_replicates << " replicates..." << endl;
        if (params.aLRT_test)
            cout << "Testing tree branches by aLRT parametric test..." << endl;
        if (params.aBayes_test)
            cout << "Testing tree branches by aBayes parametric test..." << endl;
        iqtree->setRootNode(params.root);
        if (iqtree->isBifurcating()) {
            iqtree->testAllBranches(params.aLRT_threshold, iqtree->getCurScore(),
                    pattern_lh, params.aLRT_replicates, params.localbp_replicates, params.aLRT_test, params.aBayes_test);
            cout << getRealTime() - mytime << " sec." << endl;
        } else {
            outWarning("Tree is multifurcating and such test is not applicable");
            params.aLRT_replicates = params.localbp_replicates = params.aLRT_test = params.aBayes_test = 0;
        }
    }

    if (params.gbo_replicates > 0) {
        cout << "Creating " << RESAMPLE_NAME << " support values..." << endl;
        if (!params.online_bootstrap)
            outError("Obsolete feature");
//            runGuidedBootstrap(params, iqtree->aln, iqtree);
        else
            iqtree->summarizeBootstrap(params);
    }

    if (params.collapse_zero_branch) {
        cout << "Collapsing near-zero internal branches... ";
        cout << iqtree->collapseInternalBranches(NULL, NULL, params.min_branch_length*4);
        cout << " collapsed" << endl;
    }

    printFinalSearchInfo(params, *iqtree, search_cpu_time, search_real_time);

    if (params.gbo_replicates && params.online_bootstrap && params.print_ufboot_trees)
        iqtree->writeUFBootTrees(params);

    if (params.gbo_replicates && params.online_bootstrap && !iqtree->isSuperTreeUnlinked()) {
        
        cout << endl << "Computing " << RESAMPLE_NAME << " consensus tree..." << endl;
        string splitsfile = params.out_prefix;
        splitsfile += ".splits.nex";
        double weight_threshold = (params.split_threshold<1) ? params.split_threshold : (params.gbo_replicates-1.0)/params.gbo_replicates;
        weight_threshold *= 100.0;
        computeConsensusTree(splitsfile.c_str(), 0, 1e6, -1,
                             weight_threshold, NULL, params.out_prefix, NULL, &params);
        // now optimize branch lengths of the consensus tree
        string current_tree = iqtree->getTreeString();
        optimizeConTree(params, iqtree);
        // revert the best tree
        iqtree->readTreeString(current_tree);
    }
    if (Params::getInstance().writeDistImdTrees) {
        cout << endl;
        cout << "Recomputing the log-likelihood of the intermediate trees ... " << endl;
        iqtree->intermediateTrees.recomputeLoglOfAllTrees(*iqtree);
    }
    
    // BUG FIX: readTreeString(bestTreeString) not needed before this line
    iqtree->printResultTree();
    iqtree->saveCheckpoint();

    if (params.upper_bound_NNI) {
        string out_file_UB = params.out_prefix;
        out_file_UB += ".UB.NNI.main";
        ofstream out_UB;
        out_UB.exceptions(ios::failbit | ios::badbit);
        out_UB.open((char *) out_file_UB.c_str(), std::ofstream::out | std::ofstream::app);
        out_UB << iqtree->leafNum << "\t" << iqtree->aln->getNSite() << "\t" << iqtree->params->upper_bound_frac << "\t"
        << iqtree->skippedNNIub << "\t" << iqtree->totalNNIub << "\t" << iqtree->getBestScore() << endl;
        //iqtree->minUB << "\t" << iqtree->meanUB/iqtree->skippedNNIub << "\t" << iqtree->maxUB << endl;
        out_UB.close();
    }

    if (params.out_file) {
        iqtree->printTree(params.out_file);
    }
    delete[] pattern_lh;

    runApproximateBranchLengths(params, *iqtree);
}

void handleGammaInvariantOptions(Params &params, IQTree &iqtree) {
    if (iqtree.getRate()->name.find("+I+G") != string::npos) {
        if (params.alpha_invar_file != NULL) { // COMPUTE TREE LIKELIHOOD BASED ON THE INPUT ALPHA AND P_INVAR VALUE
            computeLoglFromUserInputGAMMAInvar(params, iqtree);
            exit(0);
        }
        if (params.exh_ai) {
            exhaustiveSearchGAMMAInvar(params, iqtree);
            exit(0);
        }
    }
}

void handleQuartetLikelihoodMapping(Params &params, IQTree &iqtree) {
    if (params.lmap_num_quartets >= 0) {
        cout << endl << "Performing likelihood mapping with ";
        if (params.lmap_num_quartets > 0)
            cout << params.lmap_num_quartets;
        else
            cout << "all";
        cout << " quartets..." << endl;
        double lkmap_time = getRealTime();
        iqtree.doLikelihoodMapping();
        cout << "Likelihood mapping needed " << getRealTime()-lkmap_time << " seconds" << endl << endl;
    }
}

/**********************************************************
 * MULTIPLE TREE RECONSTRUCTION
 ***********************************************************/
void runMultipleTreeReconstruction(Params &params, Alignment *alignment, IQTree *tree) {
    ModelCheckpoint *model_info = new ModelCheckpoint;
    
    if (params.suppress_output_flags & OUT_TREEFILE)
        outError("Suppress .treefile not allowed with -runs option");
    string treefile_name = params.out_prefix;
    treefile_name += ".treefile";
    string runtrees_name = params.out_prefix;
    runtrees_name += ".runtrees";
    DoubleVector runLnL;
    
    if (tree->getCheckpoint()->getVector("runLnL", runLnL)) {
        cout << endl << "CHECKPOINT: " << runLnL.size() << " independent run(s) restored" << endl;
    } else if (MPIHelper::getInstance().isMaster()) {
        // first empty the runtrees file
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(runtrees_name.c_str());
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, runtrees_name);
        }
    }
    
    double start_time = getCPUTime();
    double start_real_time = getRealTime();
    
    int orig_seed = params.ran_seed;
    int run;
    int best_run = 0;
    for (run = 0; run < runLnL.size(); run++)
        if (runLnL[run] > runLnL[best_run])
            best_run = run;

    // do multiple tree reconstruction
    for (run = runLnL.size(); run < params.num_runs; run++) {

        tree->getCheckpoint()->startStruct("run" + convertIntToString(run+1));
        
        params.ran_seed = orig_seed + run*1000 + MPIHelper::getInstance().getProcessID();
        
        cout << endl << "---> START RUN NUMBER " << run + 1 << " (seed: " << params.ran_seed << ")" << endl;
        
        tree->getCheckpoint()->put("seed", params.ran_seed);
        
        // initialize random stream for replicating the run
        
        int *saved_randstream = randstream;
        init_random(params.ran_seed);
        
        IQTree *iqtree;
        if (alignment->isSuperAlignment()){
            if(params.partition_type != BRLEN_OPTIMIZE){
                iqtree = new PhyloSuperTreePlen((SuperAlignment*) alignment, (PhyloSuperTree*) tree);
            } else {
                iqtree = new PhyloSuperTree((SuperAlignment*) alignment, (PhyloSuperTree*) tree);
            }
        } else {
            // allocate heterotachy tree if neccessary
            int pos = posRateHeterotachy(alignment->model_name);
            
            if (params.num_mixlen > 1) {
                iqtree = new PhyloTreeMixlen(alignment, params.num_mixlen);
            } else if (pos != string::npos) {
                iqtree = new PhyloTreeMixlen(alignment, 0);
            } else
                iqtree = new IQTree(alignment);
        }
        
        if (!tree->constraintTree.empty()) {
            iqtree->constraintTree.readConstraint(tree->constraintTree);
        }
        
        // set checkpoint
        iqtree->setCheckpoint(tree->getCheckpoint());
        iqtree->num_precision = tree->num_precision;
        
        runTreeReconstruction(params, iqtree);
        // read in the output tree file
        stringstream ss;
        iqtree->printTree(ss);
        if (MPIHelper::getInstance().isMaster())
            try {
                ofstream tree_out;
                tree_out.exceptions(ios::failbit | ios::badbit);
                tree_out.open(runtrees_name.c_str(), ios_base::out | ios_base::app);
                tree_out.precision(10);
                tree_out << "[ lh=" << iqtree->getBestScore() << " ]";
                tree_out << ss.str() << endl;
                tree_out.close();
            } catch (ios::failure) {
                outError(ERR_WRITE_OUTPUT, runtrees_name);
            }
        // fix bug: set the model for original tree after testing
        if ((params.model_name.substr(0,4) == "TEST" || params.model_name.substr(0,2) == "MF") && tree->isSuperTree()) {
            PhyloSuperTree *stree = ((PhyloSuperTree*)tree);
            stree->part_info =  ((PhyloSuperTree*)iqtree)->part_info;
        }
        runLnL.push_back(iqtree->getBestScore());

        if (MPIHelper::getInstance().isMaster()) {
            if (params.num_bootstrap_samples > 0 && params.consensus_type == CT_CONSENSUS_TREE &&
                (run == 0 || iqtree->getBestScore() > runLnL[best_run])) {
                // 2017-12-08: optimize branch lengths of consensus tree
                // now optimize branch lengths of the consensus tree
                string current_tree = iqtree->getTreeString();
                optimizeConTree(params, iqtree);
                // revert the best tree
                iqtree->readTreeString(current_tree);
                iqtree->saveCheckpoint();
            }
        }
        if (iqtree->getBestScore() > runLnL[best_run])
            best_run = run;
            
        if (params.num_runs == 1)
            reportPhyloAnalysis(params, *iqtree, *model_info);
        delete iqtree;
        
        tree->getCheckpoint()->endStruct();
        // clear all checkpointed information
//        tree->getCheckpoint()->keepKeyPrefix("iqtree");
        tree->getCheckpoint()->putVector("runLnL", runLnL);
//        tree->getCheckpoint()->putBool("finished", false);
        tree->getCheckpoint()->dump(true);
        // restore randstream
        finish_random();
        randstream = saved_randstream;
    }
    
    cout << endl << "---> SUMMARIZE RESULTS FROM " << runLnL.size() << " RUNS" << endl << endl;
    
    cout << "Run " << best_run+1 <<  " gave best log-likelihood: " << runLnL[best_run] << endl;

    // initialize tree and model structure
    ModelsBlock *models_block = readModelsDefinition(params);
    tree->setParams(&params);
    tree->setNumThreads(params.num_threads);
    if (!tree->getModelFactory()) {
        tree->initializeModel(params, tree->aln->model_name, models_block);
    }
    if (tree->getRate()->isHeterotachy() && !tree->isMixlen()) {
        ASSERT(0 && "Heterotachy tree not properly created");
    }
    delete models_block;

    // restore the tree and model from the best run
    tree->getCheckpoint()->startStruct("run" + convertIntToString(best_run+1));
    tree->restoreCheckpoint();
    tree->getModelFactory()->restoreCheckpoint();
    tree->setCurScore(runLnL[best_run]);
    if (params.gbo_replicates && !tree->isSuperTreeUnlinked()) {
        
        string out_file = (string)params.out_prefix + ".splits";
        if (params.print_splits_file) {
            tree->boot_splits.back()->saveFile(out_file.c_str(), IN_OTHER, true);
            cout << "Split supports printed to star-dot file " << out_file << endl;
        }

        if (params.print_splits_nex_file) {
            out_file = (string)params.out_prefix + ".splits.nex";
            tree->boot_splits.back()->saveFile(out_file.c_str(), IN_NEXUS, false);
            cout << "Split supports printed to NEXUS file " << out_file << endl;
        }

        // overwrite .ufboot trees
        if (params.print_ufboot_trees)
            tree->writeUFBootTrees(params);

        // overwrite .contree
        string contree;
        if (!tree->getCheckpoint()->getString("contree", contree))
            ASSERT(0 && "Couldn't restore contree");
        string contree_file = string(params.out_prefix) + ".contree";
        string current_tree = tree->getTreeString();
        tree->readTreeString(contree);
        tree->setRootNode(params.root);
        tree->insertTaxa(tree->removed_seqs, tree->twin_seqs);
        tree->printTree(contree_file.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        tree->readTreeString(current_tree);
        cout << "Consensus tree written to " << contree_file << endl;
    }
    tree->getCheckpoint()->endStruct();
    
    // overwrite .treefile
    tree->printResultTree();

    if (MPIHelper::getInstance().isMaster()) {
        cout << "Total CPU time for " << params.num_runs << " runs: " << (getCPUTime() - start_time) << " seconds." << endl;
        cout << "Total wall-clock time for " << params.num_runs << " runs: " << (getRealTime() - start_real_time) << " seconds." << endl << endl;
    }
    delete model_info;
}

void computeLoglFromUserInputGAMMAInvar(Params &params, IQTree &iqtree) {
    RateHeterogeneity *site_rates = iqtree.getRate();
    site_rates->setFixPInvar(true);
    site_rates->setFixGammaShape(true);
    vector<double> alphas, p_invars, logl;
    ifstream aiFile;
    aiFile.open(params.alpha_invar_file, ios_base::in);
    if (aiFile.good()) {
        double alpha, p_invar;
        while (aiFile >> alpha >> p_invar) {
            alphas.push_back(alpha);
            p_invars.push_back(p_invar);
        }
        aiFile.close();
        cout << "Computing tree logl based on the alpha and p_invar values in " << params.alpha_invar_file << " ..." <<
        endl;
    } else {
        stringstream errMsg;
        errMsg << "Could not find file: " << params.alpha_invar_file;
        outError(errMsg.str().c_str());
    }
    string aiResultsFileName = string(params.out_prefix) + "_" + string(params.alpha_invar_file) + ".results";
    ofstream aiFileResults;
    aiFileResults.open(aiResultsFileName.c_str());
    aiFileResults << fixed;
    aiFileResults.precision(4);
    DoubleVector lenvec;
    aiFileResults << "Alpha P_Invar Logl TreeLength\n";
    for (int i = 0; i < alphas.size(); i++) {
        iqtree.saveBranchLengths(lenvec);
        aiFileResults << alphas.at(i) << " " << p_invars.at(i) << " ";
        site_rates->setGammaShape(alphas.at(i));
        site_rates->setPInvar(p_invars.at(i));
        iqtree.clearAllPartialLH();
        double lh = iqtree.getModelFactory()->optimizeParameters(params.fixed_branch_length, false, 0.001);
        aiFileResults << lh << " " << iqtree.treeLength() << "\n";
        iqtree.restoreBranchLengths(lenvec);
    }
    aiFileResults.close();
    cout << "Results were written to: " << aiResultsFileName << endl;
    cout << "Wall clock time used: " << getRealTime() - params.start_real_time << endl;
}

void searchGAMMAInvarByRestarting(IQTree &iqtree) {
    if (!Params::getInstance().fixed_branch_length)
        iqtree.setCurScore(iqtree.optimizeAllBranches(1));
    else
        iqtree.setCurScore(iqtree.computeLikelihood());
    RateHeterogeneity* site_rates = (iqtree.getRate());
    double values[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    vector<double> initAlphas;
    if (Params::getInstance().randomAlpha) {
        while (initAlphas.size() < 10) {
            double initAlpha = random_double();
            initAlphas.push_back(initAlpha + iqtree.params->min_gamma_shape*2);
        }
    } else {
        initAlphas.assign(values, values+10);
    }
    double bestLogl = iqtree.getCurScore();
    double bestAlpha = 0.0;
    double bestPInvar = 0.0;
    double initPInvar = iqtree.getRate()->getPInvar();

    /* Back up branch lengths and substitutional rates */
    DoubleVector lenvec;
    DoubleVector bestLens;
    iqtree.saveBranchLengths(lenvec);
    int numRateEntries = iqtree.getModel()->getNumRateEntries();
    double *rates = new double[numRateEntries];
    double *bestRates = new double[numRateEntries];
    iqtree.getModel()->getRateMatrix(rates);
    int numStates = iqtree.aln->num_states;
    double *state_freqs = new double[numStates];
    iqtree.getModel()->getStateFrequency(state_freqs);
    double *bestStateFreqs =  new double[numStates];

    for (int i = 0; i < 10; i++) {
        cout << endl;
        cout << "Testing alpha: " << initAlphas[i] << endl;
        // Initialize model parameters
        iqtree.restoreBranchLengths(lenvec);
        ((ModelMarkov*) iqtree.getModel())->setRateMatrix(rates);
        ((ModelMarkov*) iqtree.getModel())->setStateFrequency(state_freqs);
        iqtree.getModel()->decomposeRateMatrix();
        site_rates->setGammaShape(initAlphas[i]);
        site_rates->setPInvar(initPInvar);
        iqtree.clearAllPartialLH();
        iqtree.optimizeModelParameters(verbose_mode >= VB_MED, Params::getInstance().testAlphaEps);
        double estAlpha = iqtree.getRate()->getGammaShape();
        double estPInv = iqtree.getRate()->getPInvar();
        double logl = iqtree.getCurScore();
        cout << "Est. alpha: " << estAlpha << " / Est. pinv: " << estPInv
        << " / Logl: " << logl << endl;

        if (iqtree.getCurScore() > bestLogl) {
            bestLogl = logl;
            bestAlpha = estAlpha;
            bestPInvar = estPInv;
            bestLens.clear();
            iqtree.saveBranchLengths(bestLens);
            iqtree.getModel()->getRateMatrix(bestRates);
            iqtree.getModel()->getStateFrequency(bestStateFreqs);
        }
    }
    site_rates->setGammaShape(bestAlpha);
    site_rates->setFixGammaShape(false);
    site_rates->setPInvar(bestPInvar);
    site_rates->setFixPInvar(false);
    ((ModelMarkov*) iqtree.getModel())->setRateMatrix(bestRates);
    ((ModelMarkov*) iqtree.getModel())->setStateFrequency(bestStateFreqs);
    iqtree.restoreBranchLengths(bestLens);
    iqtree.getModel()->decomposeRateMatrix();
    iqtree.clearAllPartialLH();
    iqtree.setCurScore(iqtree.computeLikelihood());
    cout << endl;
    cout << "Best initial alpha: " << bestAlpha << " / initial pinv: " << bestPInvar << " / ";
    cout << "Logl: " << iqtree.getCurScore() << endl;

    delete [] rates;
    delete [] state_freqs;
    delete [] bestRates;
    delete [] bestStateFreqs;
}

// Test alpha fom 0.1 to 15 and p_invar from 0.1 to 0.99, stepsize = 0.01
void exhaustiveSearchGAMMAInvar(Params &params, IQTree &iqtree) {
    double alphaMin = 0.01;
    double alphaMax = 10.00;
    double p_invarMin = 0.01;
    double p_invarMax = 1.00;
//    double p_invarMax = iqtree.aln->frac_const_sites;
    double stepSize = 0.01;
    int numAlpha = (int) floor((alphaMax - alphaMin)/stepSize);
    int numInvar = (int) floor((p_invarMax - p_invarMin)/stepSize);

    cout << "EVALUATING " << numAlpha*numInvar << " COMBINATIONS OF " << " alpha=" << alphaMin << ".." << alphaMax
         << " AND " << " p-invar=" << p_invarMin << ".." << p_invarMax
         << " (epsilon: " << params.modelEps << ")" << endl;

//    vector<string> results;
//    results.reserve((unsigned long) (numAlpha * numInvar));
    DoubleVector lenvec;
    iqtree.saveBranchLengths(lenvec);

    RateHeterogeneity* site_rates = (iqtree.getRate());
    site_rates->setFixPInvar(true);
    site_rates->setFixGammaShape(true);

    string aiResultsFileName = string(params.out_prefix) + ".ai_results";
    ofstream aiFileResults;
    aiFileResults.open(aiResultsFileName.c_str());
    aiFileResults << fixed;
    aiFileResults.precision(4);
    aiFileResults << "alpha p_invar logl tree_len\n";

    for (double alpha = alphaMin; alpha < alphaMax; alpha = alpha + stepSize) {
        cout << "alpha = " << alpha << endl;
        for (double p_invar = p_invarMin; p_invar < p_invarMax; p_invar = p_invar + stepSize) {
            site_rates->setGammaShape(alpha);
            site_rates->setPInvar(p_invar);
            iqtree.clearAllPartialLH();
            double lh = iqtree.getModelFactory()->optimizeParameters(params.fixed_branch_length, false, params.modelEps);
//            stringstream ss;
//            ss << fixed << setprecision(2) << alpha << " " << p_invar << " " << lh << " " << iqtree.treeLength();
            aiFileResults << alpha << " " << p_invar << " " << lh << " " << iqtree.treeLength() << endl;
            //cout << ss.str() << endl;
//            results.push_back(ss.str());
            iqtree.restoreBranchLengths(lenvec);
        }
    }
//    for (vector<string>::iterator it = results.begin(); it != results.end(); it++) {
//                aiFileResults << (*it) << endl;
//            }
    aiFileResults.close();
    cout << "Results were written to: " << aiResultsFileName << endl;
    cout << "Wall clock time used: " << getRealTime() - params.start_real_time << endl;
}

/**********************************************************
 * STANDARD NON-PARAMETRIC BOOTSTRAP
 ***********************************************************/
void runStandardBootstrap(Params &params, Alignment *alignment, IQTree *tree) {
    ModelCheckpoint *model_info = new ModelCheckpoint;
    StrVector removed_seqs, twin_seqs;

    // turn off all branch tests
    int saved_aLRT_replicates = params.aLRT_replicates;
    int saved_localbp_replicates = params.localbp_replicates;
    bool saved_aLRT_test = params.aLRT_test;
    bool saved_aBayes_test = params.aBayes_test;
    params.aLRT_replicates = 0;
    params.localbp_replicates = 0;
    params.aLRT_test = false;
    params.aBayes_test = false;
    
    if (params.suppress_output_flags & OUT_TREEFILE)
        outError("Suppress .treefile not allowed for standard bootstrap");
    string treefile_name = params.out_prefix;
    treefile_name += ".treefile";
    string boottrees_name = params.out_prefix;
    boottrees_name += ".boottrees";
    string bootaln_name = params.out_prefix;
    bootaln_name += ".bootaln";
    string bootlh_name = params.out_prefix;
    bootlh_name += ".bootlh";
    int bootSample = 0;
    if (tree->getCheckpoint()->get("bootSample", bootSample)) {
        cout << "CHECKPOINT: " << bootSample << " bootstrap analyses restored" << endl;
    } else if (MPIHelper::getInstance().isMaster()) {
        // first empty the boottrees file
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(boottrees_name.c_str());
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, boottrees_name);
        }

        // empty the bootaln file
        if (params.print_bootaln)
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(bootaln_name.c_str());
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, bootaln_name);
        }
    }

    double start_time = getCPUTime();
    double start_real_time = getRealTime();

    startTreeReconstruction(params, tree, *model_info);
    
    // 2018-06-21: bug fix: alignment might be changed by -m ...MERGE
    alignment = tree->aln;
    
    // do bootstrap analysis
    for (int sample = bootSample; sample < params.num_bootstrap_samples; sample++) {
        cout << endl << "===> START " << RESAMPLE_NAME_UPPER << " REPLICATE NUMBER "
                << sample + 1 << endl << endl;

        // 2015-12-17: initialize random stream for creating bootstrap samples
        // mainly so that checkpointing does not need to save bootstrap samples
        int *saved_randstream = randstream;
        init_random(params.ran_seed + sample);

        Alignment* bootstrap_alignment;
        cout << "Creating " << RESAMPLE_NAME << " alignment (seed: " << params.ran_seed+sample << ")..." << endl;

        if (alignment->isSuperAlignment())
            bootstrap_alignment = new SuperAlignment;
        else
            bootstrap_alignment = new Alignment;
        bootstrap_alignment->createBootstrapAlignment(alignment, NULL, params.bootstrap_spec);

        // restore randstream
        finish_random();
        randstream = saved_randstream;

        if (params.print_tree_lh && MPIHelper::getInstance().isMaster()) {
            double prob;
            bootstrap_alignment->multinomialProb(*alignment, prob);
            ofstream boot_lh;
            if (sample == 0)
                boot_lh.open(bootlh_name.c_str());
            else
                boot_lh.open(bootlh_name.c_str(), ios_base::out | ios_base::app);
            boot_lh << "0\t" << prob << endl;
            boot_lh.close();
        }
        IQTree *boot_tree;
        if (alignment->isSuperAlignment()){
            if(params.partition_type != BRLEN_OPTIMIZE){
                boot_tree = new PhyloSuperTreePlen((SuperAlignment*) bootstrap_alignment, (PhyloSuperTree*) tree);
            } else {
                boot_tree = new PhyloSuperTree((SuperAlignment*) bootstrap_alignment, (PhyloSuperTree*) tree);
            }
        } else {
            // allocate heterotachy tree if neccessary
            int pos = posRateHeterotachy(alignment->model_name);
            
            if (params.num_mixlen > 1) {
                boot_tree = new PhyloTreeMixlen(bootstrap_alignment, params.num_mixlen);
            } else if (pos != string::npos) {
                boot_tree = new PhyloTreeMixlen(bootstrap_alignment, 0);
            } else
                boot_tree = new IQTree(bootstrap_alignment);
        }
        if (params.print_bootaln && MPIHelper::getInstance().isMaster()) {
            bootstrap_alignment->printAlignment(params.aln_output_format, bootaln_name.c_str(), true);
        }

        if (params.print_boot_site_freq && MPIHelper::getInstance().isMaster()) {
            printSiteStateFreq((((string)params.out_prefix)+"."+convertIntToString(sample)+".bootsitefreq").c_str(), bootstrap_alignment);
                bootstrap_alignment->printAlignment(params.aln_output_format, (((string)params.out_prefix)+"."+convertIntToString(sample)+".bootaln").c_str());
        }

        if (!tree->constraintTree.empty()) {
            boot_tree->constraintTree.readConstraint(tree->constraintTree);
        }

        // set checkpoint
        boot_tree->setCheckpoint(tree->getCheckpoint());
        boot_tree->num_precision = tree->num_precision;

        runTreeReconstruction(params, boot_tree);
        // read in the output tree file
        stringstream ss;
        boot_tree->printTree(ss);
//        try {
//            ifstream tree_in;
//            tree_in.exceptions(ios::failbit | ios::badbit);
//            tree_in.open(treefile_name.c_str());
//            tree_in >> tree_str;
//            tree_in.close();
//        } catch (ios::failure) {
//            outError(ERR_READ_INPUT, treefile_name);
//        }
        // write the tree into .boottrees file
        if (MPIHelper::getInstance().isMaster())
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(boottrees_name.c_str(), ios_base::out | ios_base::app);
            tree_out << ss.str() << endl;
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, boottrees_name);
        }
        // OBSOLETE fix bug: set the model for original tree after testing
//        if ((params.model_name.substr(0,4) == "TEST" || params.model_name.substr(0,2) == "MF") && tree->isSuperTree()) {
//            PhyloSuperTree *stree = ((PhyloSuperTree*)tree);
//            stree->part_info =  ((PhyloSuperTree*)boot_tree)->part_info;
//        }
        if (params.num_bootstrap_samples == 1)
            reportPhyloAnalysis(params, *boot_tree, *model_info);
        // WHY was the following line missing, which caused memory leak?
        bootstrap_alignment = boot_tree->aln;
        delete boot_tree;
        // fix bug: bootstrap_alignment might be changed
        delete bootstrap_alignment;

        // clear all checkpointed information
        tree->getCheckpoint()->keepKeyPrefix("iqtree");
        tree->getCheckpoint()->put("bootSample", sample+1);
        tree->getCheckpoint()->putBool("finished", false);
        tree->getCheckpoint()->dump(true);
    }


    if (params.consensus_type == CT_CONSENSUS_TREE && MPIHelper::getInstance().isMaster()) {

        cout << endl << "===> COMPUTE CONSENSUS TREE FROM " << params.num_bootstrap_samples
            << RESAMPLE_NAME_UPPER << " TREES" << endl << endl;
        string root_name = (params.root) ? params.root : alignment->getSeqName(0);
        const char* saved_root = params.root;
        params.root = root_name.c_str();
        computeConsensusTree(boottrees_name.c_str(), 0, 1e6, -1,
                params.split_threshold, NULL, params.out_prefix, NULL, &params);
        params.root = saved_root;
    }

    if (params.compute_ml_tree) {
        cout << endl << "===> START ANALYSIS ON THE ORIGINAL ALIGNMENT" << endl << endl;
        // restore branch tests
        params.aLRT_replicates = saved_aLRT_replicates;
        params.localbp_replicates = saved_localbp_replicates;
        params.aLRT_test = saved_aLRT_test;
        params.aBayes_test = saved_aBayes_test;

        if (params.num_runs == 1)
            runTreeReconstruction(params, tree);
        else
            runMultipleTreeReconstruction(params, tree->aln, tree);

        if (MPIHelper::getInstance().isMaster()) {
            if (params.consensus_type == CT_CONSENSUS_TREE && params.num_runs == 1) {
                // 2017-12-08: optimize branch lengths of consensus tree
                optimizeConTree(params, tree);
            }

            cout << endl << "===> ASSIGN " << RESAMPLE_NAME_UPPER
                << " SUPPORTS TO THE TREE FROM ORIGINAL ALIGNMENT" << endl << endl;
            MExtTree ext_tree;
            assignBootstrapSupport(boottrees_name.c_str(), 0, 1e6,
                    treefile_name.c_str(), false, treefile_name.c_str(),
                    params.out_prefix, ext_tree, NULL, &params);
            tree->copyTree(&ext_tree);
            reportPhyloAnalysis(params, *tree, *model_info);
        }
    } else if (params.consensus_type == CT_CONSENSUS_TREE && MPIHelper::getInstance().isMaster()) {
        int mi = params.min_iterations;
        STOP_CONDITION sc = params.stop_condition;
        params.min_iterations = 0;
        params.stop_condition = SC_FIXED_ITERATION;
        runTreeReconstruction(params, tree);
        params.min_iterations = mi;
        params.stop_condition = sc;
        tree->stop_rule.initialize(params);
        optimizeConTree(params, tree);
        reportPhyloAnalysis(params, *tree, *model_info);
    } else
        cout << endl;

#ifdef USE_BOOSTER
    if (params.transfer_bootstrap) {
        // transfer bootstrap expectation (TBE)
        cout << "Performing transfer bootstrap expectation..." << endl;
        string input_tree = (string)params.out_prefix + ".treefile";
        string boot_trees = (string)params.out_prefix + ".boottrees";
        string out_tree = (string)params.out_prefix + ".tbe.tree";
        string out_raw_tree = (string)params.out_prefix + ".tbe.rawtree";
        string stat_out = (string)params.out_prefix + ".tbe.stat";
        main_booster(input_tree.c_str(), boot_trees.c_str(), out_tree.c_str(),
                     (params.transfer_bootstrap==2) ? out_raw_tree.c_str() : NULL,
                     stat_out.c_str(), (verbose_mode >= VB_MED) ? 0 : 1);
        cout << "TBE tree written to " << out_tree << endl;
        if (params.transfer_bootstrap == 2)
            cout << "TBE raw tree written to " << out_raw_tree << endl;
        cout << "TBE statistic written to " << stat_out << endl;
        cout << endl;
    }
#endif
    
    if (MPIHelper::getInstance().isMaster()) {
        cout << "Total CPU time for " << RESAMPLE_NAME << ": " << (getCPUTime() - start_time) << " seconds." << endl;
    cout << "Total wall-clock time for " << RESAMPLE_NAME << ": " << (getRealTime() - start_real_time) << " seconds." << endl << endl;
    cout << "Non-parametric " << RESAMPLE_NAME << " results written to:" << endl;
    if (params.print_bootaln)
        cout << RESAMPLE_NAME_I << " alignments:     " << params.out_prefix << ".bootaln" << endl;
    cout << RESAMPLE_NAME_I << " trees:          " << params.out_prefix << ".boottrees" << endl;
    if (params.consensus_type == CT_CONSENSUS_TREE)
        cout << "  Consensus tree:           " << params.out_prefix << ".contree" << endl;
    cout << endl;
    }
    delete model_info;
}

void convertAlignment(Params &params, IQTree *iqtree) {
    Alignment *alignment = iqtree->aln;
    if (params.num_bootstrap_samples || params.print_bootaln) {
        // create bootstrap alignment
        Alignment* bootstrap_alignment;
        cout << "Creating " << RESAMPLE_NAME << " alignment..." << endl;
        if (alignment->isSuperAlignment())
            bootstrap_alignment = new SuperAlignment;
        else
            bootstrap_alignment = new Alignment;
        bootstrap_alignment->createBootstrapAlignment(alignment, NULL, params.bootstrap_spec);
        delete alignment;
        alignment = bootstrap_alignment;
        iqtree->aln = alignment;
    }

    int exclude_sites = 0;
    if (params.aln_nogaps)
        exclude_sites += EXCLUDE_GAP;
    if (params.aln_no_const_sites)
        exclude_sites += EXCLUDE_INVAR;

    if (alignment->isSuperAlignment()) {
        alignment->printAlignment(params.aln_output_format, params.aln_output, false, params.aln_site_list,
                                  exclude_sites, params.ref_seq_name);
        if (params.print_subaln)
            ((SuperAlignment*)alignment)->printSubAlignments(params);
        if (params.aln_output_format != IN_NEXUS) {
            string partition_info = string(params.aln_output) + ".nex";
            ((SuperAlignment*)alignment)->printPartition(partition_info.c_str(), params.aln_output);
            partition_info = (string)params.aln_output + ".partitions";
            ((SuperAlignment*)alignment)->printPartitionRaxml(partition_info.c_str());
        }
    } else if (params.gap_masked_aln) {
        Alignment out_aln;
        Alignment masked_aln(params.gap_masked_aln, params.sequence_type, params.intype, params.model_name);
        out_aln.createGapMaskedAlignment(&masked_aln, alignment);
        out_aln.printAlignment(params.aln_output_format, params.aln_output, false, params.aln_site_list,
                exclude_sites, params.ref_seq_name);
        string str = params.gap_masked_aln;
        str += ".sitegaps";
        out_aln.printSiteGaps(str.c_str());
    } else  {
        alignment->printAlignment(params.aln_output_format, params.aln_output, false, params.aln_site_list,
                exclude_sites, params.ref_seq_name);
    }
}

/**
    2016-08-04: compute a site frequency model for profile mixture model
*/
void computeSiteFrequencyModel(Params &params, Alignment *alignment) {

    cout << endl << "===> COMPUTING SITE FREQUENCY MODEL BASED ON TREE FILE " << params.tree_freq_file << endl;
    ASSERT(params.tree_freq_file);
    PhyloTree *tree = new PhyloTree(alignment);
    tree->setParams(&params);
    bool myrooted = params.is_rooted;
    tree->readTree(params.tree_freq_file, myrooted);
    tree->setAlignment(alignment);
    tree->setRootNode(params.root);
    
    ModelsBlock *models_block = readModelsDefinition(params);
    tree->setModelFactory(new ModelFactory(params, alignment->model_name, tree, models_block));
    delete models_block;
    tree->setModel(tree->getModelFactory()->model);
    tree->setRate(tree->getModelFactory()->site_rate);
    tree->setLikelihoodKernel(params.SSE);
    tree->setNumThreads(params.num_threads);
    
    if (!tree->getModel()->isMixture())
        outError("No mixture model was specified!");
    uint64_t mem_size = tree->getMemoryRequired();
    uint64_t total_mem = getMemorySize();
    cout << "NOTE: " << (mem_size / 1024) / 1024 << " MB RAM is required!" << endl;
    if (mem_size >= total_mem) {
        outError("Memory required exceeds your computer RAM size!");
    }
#ifdef BINARY32
    if (mem_size >= 2000000000) {
        outError("Memory required exceeds 2GB limit of 32-bit executable");
    }
#endif

    tree->ensureNumberOfThreadsIsSet(nullptr, false);

    tree->initializeAllPartialLh();
    // 2017-12-07: Increase espilon ten times (0.01 -> 0.1) to speedup PMSF computation
    tree->getModelFactory()->optimizeParameters(params.fixed_branch_length, true, params.modelEps*10);

    size_t nptn = alignment->getNPattern(), nstates = alignment->num_states;
    double *ptn_state_freq = new double[nptn*nstates];
    tree->computePatternStateFreq(ptn_state_freq);
    alignment->site_state_freq.resize(nptn);
    for (size_t ptn = 0; ptn < nptn; ptn++) {
        double *f = new double[nstates];
        memcpy(f, ptn_state_freq+ptn*nstates, sizeof(double)*nstates);
        alignment->site_state_freq[ptn] = f;
    }
    alignment->getSitePatternIndex(alignment->site_model);
    printSiteStateFreq(((string)params.out_prefix+".sitefreq").c_str(), tree, ptn_state_freq);
    params.print_site_state_freq = WSF_NONE;
    
    delete [] ptn_state_freq;
    delete tree;
    
    cout << endl << "===> CONTINUE ANALYSIS USING THE INFERRED SITE FREQUENCY MODEL" << endl;
}


/**********************************************************
 * TOP-LEVEL FUNCTION
 ***********************************************************/

IQTree *newIQTree(Params &params, Alignment *alignment) {
    IQTree *tree;
    if (alignment->isSuperAlignment()) {
        if (params.partition_type == TOPO_UNLINKED) {
            tree = new PhyloSuperTreeUnlinked((SuperAlignment*)alignment);
        } else if(params.partition_type != BRLEN_OPTIMIZE){
            // initialize supertree - Proportional Edges case
            tree = new PhyloSuperTreePlen((SuperAlignment*)alignment, params.partition_type);
        } else {
            // initialize supertree stuff if user specifies partition file with -sp option
            tree = new PhyloSuperTree((SuperAlignment*)alignment);
        }
        // this alignment will actually be of type SuperAlignment
        //        alignment = tree->aln;
        if (((PhyloSuperTree*)tree)->rescale_codon_brlen)
            cout << "NOTE: Mixed codon and other data, branch lengths of codon partitions are rescaled by 3!" << endl;
        
    } else {
        // allocate heterotachy tree if neccessary
        int pos = posRateHeterotachy(alignment->model_name);
        
        if (params.num_mixlen > 1) {
            tree = new PhyloTreeMixlen(alignment, params.num_mixlen);
        } else if (pos != string::npos) {
            tree = new PhyloTreeMixlen(alignment, 0);
        } else
            tree = new IQTree(alignment);
    }

    return tree;
}

/** get ID of bad or good symtest results */
void getSymTestID(vector<SymTestResult> &res, set<int> &id, bool bad_res) {
    if (bad_res) {
        // get significant test ID
        switch (Params::getInstance().symtest) {
            case SYMTEST_BINOM:
                for (auto i = res.begin(); i != res.end(); i++)
                    if (i->pvalue_binom < Params::getInstance().symtest_pcutoff)
                        id.insert(i - res.begin());
                break;
            case SYMTEST_MAXDIV:
                for (auto i = res.begin(); i != res.end(); i++)
                    if (i->pvalue_maxdiv < Params::getInstance().symtest_pcutoff)
                        id.insert(i - res.begin());
                break;
            default:
                break;
        }
    } else {
        // get non-significant test ID
        switch (Params::getInstance().symtest) {
            case SYMTEST_BINOM:
                for (auto i = res.begin(); i != res.end(); i++)
                    if (i->pvalue_binom >= Params::getInstance().symtest_pcutoff)
                        id.insert(i - res.begin());
                break;
            case SYMTEST_MAXDIV:
                for (auto i = res.begin(); i != res.end(); i++)
                    if (i->pvalue_maxdiv >= Params::getInstance().symtest_pcutoff)
                        id.insert(i - res.begin());
                break;
            default:
                break;
        }
    }
}

double computePValueSMax(vector<SymTestResult> &sym, int start, int step) {
    double orig_max = sym[start].max_stat;
    int count = 0, num = 0;
    for (size_t i = start; i < sym.size(); i += step, num++)
        if (sym[i].max_stat >= orig_max)
            count++;
    return double(count)/num;
    
}

void doSymTest(Alignment *alignment, Params &params) {
    double start_time = getRealTime();
    cout << "Performing matched-pair tests of symmetry...";
    vector<SymTestResult> sym, marsym, intsym;

    size_t num_parts = 1;
    if (alignment->isSuperAlignment())
        num_parts = ((SuperAlignment*)alignment)->partitions.size();
    
    string filename_stat = string(params.out_prefix) + ".symstat.csv";
    ofstream *out_stat = NULL;
    if (params.symtest_stat) {
        out_stat = new ofstream;
        out_stat->open(filename_stat);
        *out_stat
        << "# Statistic values for matched-pair tests of symmetry" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#    dat=read.csv('" <<  filename_stat << "',comment.char='#')" << endl
        << "# Columns are comma-separated with following meanings:" << endl
        << "#    ID:     Partition ID" << endl
        << "#    Seq1:   ID of sequence 1 within partition" << endl
        << "#    Seq1:   ID of sequence 2 within partition" << endl
        << "#    Sym:    Statistic for test of symmetry" << endl
        << "#    SymChi: Chi-square p-value for test of symmetry" << endl
        << "#    Mar:    Statistic for test of marginal symmetry" << endl
        << "#    MarChi: Chi-square p-value for marginal test of symmetry" << endl
        << "#    Int:    Statistic for test of internal symmetry" << endl
        << "#    MarChi: Chi-square p-value for internal test of symmetry" << endl
        << "ID,Seq1,Seq2,Sym,SymChi,Mar,MarChi,Int,IntChi" << endl;
        
    }

    sym.resize(num_parts*params.symtest_shuffle);
    marsym.resize(num_parts*params.symtest_shuffle);
    intsym.resize(num_parts*params.symtest_shuffle);

    for (int i = 0; i < params.symtest_shuffle; i++) {
        vector<SymTestStat> *stats = NULL;
        if (params.symtest_stat)
            stats = new vector<SymTestStat>;
        if (i == 0) // original alignment
            alignment->doSymTest(i*num_parts, sym, marsym, intsym, NULL, stats);
        else {
            int *rstream;
            init_random(params.ran_seed+i+1, false, &rstream);
            alignment->doSymTest(i*num_parts, sym, marsym, intsym, rstream, stats);
            finish_random(rstream);
        }
        if ((i+1)*10 % params.symtest_shuffle == 0) {
            cout << " " << (i+1)*100 / params.symtest_shuffle << "%";
            cout.flush();
        }
        if (!stats)
            continue;
        for (auto it = stats->begin(); it != stats->end(); it++) {
            *out_stat << it->part << ',' << it->seq1 << ',' << it->seq2 << ','
            << it->chi2_sym << ',' << it->pval_sym << ','
            << it->chi2_marsym << ',' << it->pval_marsym << ','
            << it->chi2_intsym << ',' << it->pval_intsym << endl;
        }
        delete stats;
    }

    if (out_stat) {
        out_stat->close();
        delete out_stat;
    }
    
    if (params.symtest_shuffle > 1) {
        // compute p-value for s-max approach
        for (int part = 0; part < num_parts; part++) {
            sym[part].pvalue_perm = computePValueSMax(sym, part, num_parts);
            marsym[part].pvalue_perm = computePValueSMax(marsym, part, num_parts);
            intsym[part].pvalue_perm = computePValueSMax(intsym, part, num_parts);
        }
    }

    string filename = string(params.out_prefix) + ".symtest.csv";
    ofstream out;
    out.open(filename);
    out << "# Matched-pair tests of symmetry" << endl
    << "# This file can be read in MS Excel or in R with command:" << endl
    << "#    dat=read.csv('" <<  filename << "',comment.char='#')" << endl
    << "# Columns are comma-separated with following meanings:" << endl
    << "#    Name:    Partition name" << endl
    << "#    SymSig:  Number of significant sequence pairs by test of symmetry" << endl
    << "#    SymNon:  Number of non-significant sequence pairs by test of symmetry" << endl
    << ((Params::getInstance().symtest == SYMTEST_BINOM) ? "#    SymBi:   P-value for binomial test of symmetry" : "#    SymPval: P-value for maximum test of symmetry") << endl;
    if (params.symtest_shuffle > 1)
        out << "#    SymMax:  Maximum of pair statistics by test of symmetry" << endl
            << "#    SymPerm: P-value for permutation test of symmetry" << endl;
    
    out << "#    MarSig:  Number of significant sequence pairs by test of marginal symmetry" << endl
    << "#    MarNon:  Number of non-significant sequence pairs by test of marginal symmetry" << endl
    << ((Params::getInstance().symtest == SYMTEST_BINOM) ? "#    MarBi:   P-value for binomial test of marginal symmetry" : "#    MarPval: P-value for maximum test of marginal symmetry") << endl;
    if (params.symtest_shuffle > 1)
        out << "#    MarMax:  Maximum of pair statistics by test of marginal symmetry" << endl
            << "#    MarPerm: P-value for permutation test of marginal symmetry" << endl;
    out << "#    IntSig:  Number of significant sequence pairs by test of internal symmetry" << endl
    << "#    IntNon:  Number of non-significant sequence pairs by test of internal symmetry" << endl
    << ((Params::getInstance().symtest == SYMTEST_BINOM) ? "#    IntBi:   P-value for binomial test of symmetry" : "#    IntPval: P-value for maximum test of internal symmetry") << endl;
    if (params.symtest_shuffle > 1)
        out << "#    IntMax:  Maximum of pair statistics by test of internal symmetry" << endl
        << "#    IntPerm: P-value for permutation test of internal symmetry" << endl;
    
    out << "Name,SymSig,SymNon," << ((Params::getInstance().symtest == SYMTEST_BINOM) ? "SymBi" : "SymPval")
        << ((params.symtest_shuffle > 1) ? ",SymMax,SymPerm" : "")
        << ",MarSig,MarNon," << ((Params::getInstance().symtest == SYMTEST_BINOM) ? "MarBi" : "MarPval")
        << ((params.symtest_shuffle > 1) ? ",MarMax,MarPerm" : "")
        << ",IntSig,IntNon," << ((Params::getInstance().symtest == SYMTEST_BINOM) ? "IntBi" : "IntPval")
        << ((params.symtest_shuffle > 1) ? ",IntMax,IntPerm" : "") << endl;
    
    if (alignment->isSuperAlignment()) {
        SuperAlignment *saln = (SuperAlignment*)alignment;
        for (int part = 0; part < saln->partitions.size(); part++)
            out << saln->partitions[part]->name << ','
                << sym[part] << ',' << marsym[part] << ','  << intsym[part] << endl;
    } else {
        out << alignment->name << ',' << sym[0] << ',' << marsym[0] << ','  << intsym[0] << endl;
    }

    if (params.symtest_shuffle > 1) {
        for (int part = num_parts; part < sym.size(); part++) {
            sym[part].pvalue_perm = marsym[part].pvalue_perm = intsym[part].pvalue_perm = -1.0;
            out << part % num_parts << ','
            << sym[part] << ',' << marsym[part] << ','  << intsym[part] << endl;
        }
        // erase the rest
        sym.erase(sym.begin()+num_parts, sym.end());
        marsym.erase(marsym.begin()+num_parts, marsym.end());
        intsym.erase(intsym.begin()+num_parts, intsym.end());
    }
    
    out.close();
    cout << " " << getRealTime() - start_time << " seconds" << endl;
    if (params.symtest_stat)
        cout << "SymTest statistics written to " << filename_stat << endl;
    cout << "SymTest results written to " << filename << endl;

    // now filter out partitions
    if (alignment->isSuperAlignment()) {
        set<int> part_id;
        if (params.symtest_remove == 1) {
            // remove bad loci
            if (params.symtest_type == 0)
                getSymTestID(sym, part_id, true);
            else if (params.symtest_type == 1)
                getSymTestID(marsym, part_id, true);
            else
                getSymTestID(intsym, part_id, true);
        } else if (params.symtest_remove == 2) {
            // remove good loci
            if (params.symtest_type == 0)
                getSymTestID(sym, part_id, false);
            else if (params.symtest_type == 1)
                getSymTestID(marsym, part_id, false);
            else
                getSymTestID(intsym, part_id, false);
        }
        if (!part_id.empty()) {
            SuperAlignment *saln = (SuperAlignment*)alignment;
            cout << "Removing " << part_id.size()
            << ((params.symtest_remove == 1)? " bad" : " good") << " partitions (pvalue cutoff = "
            << params.symtest_pcutoff << ")..." << endl;
            if (part_id.size() < alignment->getNSite())
                saln->removePartitions(part_id);
            else
                outError("Can't remove all partitions");
            if (params.aln_output_format == IN_NEXUS) {
                string aln_file = (string)params.out_prefix + ((params.symtest_remove == 1)? ".good.nex" : ".bad.nex");
                alignment->printAlignment(params.aln_output_format, aln_file.c_str());
            } else {
                string aln_file = (string)params.out_prefix + ((params.symtest_remove == 1)? ".good.phy" : ".bad.phy");
                alignment->printAlignment(params.aln_output_format, aln_file.c_str());
                string nexus_filename = (string)params.out_prefix + ((params.symtest_remove == 2)? ".good.nex" : ".bad.nex");
                saln->printPartition(nexus_filename.c_str(), aln_file.c_str());
            }
        }
    }
    if (params.symtest_only)
        exit(EXIT_SUCCESS);
}

void runPhyloAnalysis(Params &params, Checkpoint *checkpoint) {
    Alignment *alignment;

    checkpoint->putBool("finished", false);
    checkpoint->setDumpInterval(params.checkpoint_dump_interval);

    /****************** read in alignment **********************/
    if (params.partition_file) {
        // Partition model analysis
        if (params.partition_type == TOPO_UNLINKED)
            alignment = new SuperAlignmentUnlinked(params);
        else
            alignment = new SuperAlignment(params);
    } else {
        alignment = createAlignment(params.aln_file, params.sequence_type, params.intype, params.model_name);

        if (params.freq_const_patterns) {
            int orig_nsite = alignment->getNSite();
            alignment->addConstPatterns(params.freq_const_patterns);
            cout << "INFO: " << alignment->getNSite() - orig_nsite << " const sites added into alignment" << endl;
        }

        // Initialize site-frequency model
        if (params.tree_freq_file) {
            if (checkpoint->getBool("finishedSiteFreqFile")) {
                alignment->readSiteStateFreq(((string)params.out_prefix + ".sitefreq").c_str());
                params.print_site_state_freq = WSF_NONE;
                cout << "CHECKPOINT: Site frequency model restored" << endl;
            } else {
                computeSiteFrequencyModel(params, alignment);
                checkpoint->putBool("finishedSiteFreqFile", true);
                checkpoint->dump();
            }
        }
        if (params.site_freq_file) {
            alignment->readSiteStateFreq(params.site_freq_file);
        }
    }

    if (params.symtest) {
        doSymTest(alignment, params);
    }

    if (params.print_aln_info) {
        string site_info_file = string(params.out_prefix) + ".alninfo";
        alignment->printSiteInfo(site_info_file.c_str());
        cout << "Alignment sites statistics printed to " << site_info_file << endl;
    }

    /*************** initialize tree ********************/
    IQTree *tree = newIQTree(params, alignment);
    
    tree->setCheckpoint(checkpoint);
    if (params.min_branch_length <= 0.0) {
        params.min_branch_length = 1e-6;
        if (!tree->isSuperTree() && tree->getAlnNSite() >= 100000) {
            params.min_branch_length = 0.1 / (tree->getAlnNSite());
            tree->num_precision = max((int)ceil(-log10(Params::getInstance().min_branch_length))+1, 6);
            cout.precision(12);
            cout << "NOTE: minimal branch length is reduced to " << params.min_branch_length << " for long alignment" << endl;
            cout.precision(3);
        }
        // Increase the minimum branch length if PoMo is used.
        if (alignment->seq_type == SEQ_POMO) {
            params.min_branch_length *= alignment->virtual_pop_size * alignment->virtual_pop_size;
            cout.precision(12);
            cout << "NOTE: minimal branch length is increased to " << params.min_branch_length << " because PoMo infers number of mutations and frequency shifts" << endl;
            cout.precision(3);
        }
    }
    // Increase the minimum branch length if PoMo is used.
    if (alignment->seq_type == SEQ_POMO) {
        params.max_branch_length *= alignment->virtual_pop_size * alignment->virtual_pop_size;
        cout.precision(1);
        cout << "NOTE: maximal branch length is increased to " << params.max_branch_length << " because PoMo infers number of mutations and frequency shifts" << endl;
        cout.precision(3);
    }


    if (params.concatenate_aln) {
        Alignment aln(params.concatenate_aln, params.sequence_type, params.intype, params.model_name);
        cout << "Concatenating " << params.aln_file << " with " << params.concatenate_aln << " ..." << endl;
        alignment->concatenateAlignment(&aln);
    }

    if (params.constraint_tree_file) {
        cout << "Reading constraint tree " << params.constraint_tree_file << "..." << endl;
        tree->constraintTree.readConstraint(params.constraint_tree_file, alignment->getSeqNames());
        if (params.start_tree == STT_PLL_PARSIMONY)
            params.start_tree = STT_PARSIMONY;
        else if (params.start_tree == STT_BIONJ)
            outError("Constraint tree does not work with -t BIONJ");
        if (params.num_bootstrap_samples || params.gbo_replicates)
            cout << "INFO: Constraint tree will be applied to ML tree and all bootstrap trees." << endl;
    }

    if (params.compute_seq_identity_along_tree) {
        if (!params.user_file)
            outError("Please supply a user tree file!");
        tree->readTree(params.user_file, params.is_rooted);
        if (!tree->rooted && !params.root) {
            outError("Tree is unrooted, thus you have to specify a root with -o option");
        }
        tree->setAlignment(tree->aln);
        if (!tree->rooted)
            tree->setRootNode(params.root);
        tree->computeSeqIdentityAlongTree();
        if (verbose_mode >= VB_MED)
            tree->drawTree(cout);
        string out_tree = (string)params.out_prefix + ".seqident_tree";
        tree->printTree(out_tree.c_str());
        cout << "Tree with sequence identity printed to " << out_tree << endl;
    } else if (params.aln_output) {
        /************ convert alignment to other format and write to output file *************/
        convertAlignment(params, tree);
    } else if (params.gbo_replicates > 0 && params.user_file && params.second_tree) {
        // run one of the UFBoot analysis
//        runGuidedBootstrap(params, alignment, *tree);
        outError("Obsolete feature");
    } else if (params.avh_test) {
        // run one of the wondering test for Arndt
//        runAvHTest(params, alignment, *tree);
        outError("Obsolete feature");
    } else if (params.bootlh_test) {
        // run Arndt's plot of tree likelihoods against bootstrap alignments
//        runBootLhTest(params, alignment, *tree);
        outError("Obsolete feature");
    } else if (params.num_bootstrap_samples == 0) {
    /********************************************************************************
                    THE MAIN MAXIMUM LIKELIHOOD TREE RECONSTRUCTION
     ********************************************************************************/
        ModelCheckpoint *model_info = new ModelCheckpoint;
        alignment->checkGappySeq(params.remove_empty_seq);

        // remove identical sequences
        if (params.ignore_identical_seqs) {
            tree->removeIdenticalSeqs(params);
            if (tree->removed_seqs.size() > 0 && MPIHelper::getInstance().isMaster() && (params.suppress_output_flags & OUT_UNIQUESEQ) == 0) {
                string filename = (string)params.out_prefix + ".uniqueseq.phy";
                tree->aln->printAlignment(params.aln_output_format, filename.c_str());
                cout << endl << "For your convenience alignment with unique sequences printed to " << filename << endl;
            }
        }
        alignment = NULL; // from now on use tree->aln instead

        startTreeReconstruction(params, tree, *model_info);
        // call main tree reconstruction
        if (params.num_runs == 1)
            runTreeReconstruction(params, tree);
        else
            runMultipleTreeReconstruction(params, tree->aln, tree);
        
        if (MPIHelper::getInstance().isMaster()) {
            reportPhyloAnalysis(params, *tree, *model_info);
        }

        // reinsert identical sequences
        if (tree->removed_seqs.size() > 0) {
            // BUG FIX: dont use reinsertIdenticalSeqs anymore
            tree->insertTaxa(tree->removed_seqs, tree->twin_seqs);
            tree->printResultTree();
        }
        delete model_info;
        
        if (params.dating_method != "") {
            doTimeTree(tree);
        }

    } else {
        // the classical non-parameter bootstrap (SBS)
//        if (params.model_name.find("LINK") != string::npos || params.model_name.find("MERGE") != string::npos)
//            outError("-m TESTMERGE is not allowed when doing standard bootstrap. Please first\nfind partition scheme on the original alignment and use it for bootstrap analysis");
        if (alignment->getNSeq() < 4)
            outError("It makes no sense to perform bootstrap with less than 4 sequences.");
        runStandardBootstrap(params, alignment, tree);
    }

//    if (params.upper_bound) {
//            UpperBounds(&params, alignment, tree);
//    }

    if(verbose_mode >= VB_MED){
        if(tree->isSuperTree() && params.partition_type != BRLEN_OPTIMIZE){
            ((PhyloSuperTreePlen*) tree)->printNNIcasesNUM();
        }
    }
    // 2015-09-22: bug fix, move this line to before deleting tree
    alignment = tree->aln;
    delete tree;
    // BUG FIX: alignment can be changed, should delete tree->aln instead
    // 2015-09-22: THIS IS STUPID: after deleting tree, one cannot access tree->aln anymore
//    alignment = tree->aln;
    delete alignment;

    checkpoint->putBool("finished", true);
    checkpoint->dump(true);
}

/**
    Perform separate tree reconstruction when tree topologies
    are unlinked between partitions
 */
void runUnlinkedPhyloAnalysis(Params &params, Checkpoint *checkpoint) {
    SuperAlignment *super_aln;
    
    ASSERT(params.partition_file);
    
    /****************** read in alignment **********************/
    // Partition model analysis
    super_aln = new SuperAlignmentUnlinked(params);
    PhyloSuperTree *super_tree = new PhyloSuperTree(super_aln);
 
    /**** do separate tree reconstruction for each partition ***/
    
    MTreeSet part_trees;
    
    if (params.user_file) {
        // reading user tree file for all partitions
        bool is_rooted = false;
        part_trees.readTrees(params.user_file, is_rooted, 0, super_aln->partitions.size());
        if (is_rooted)
            outError("Rooted trees not allowed: ", params.user_file);
        if (part_trees.size() != super_aln->partitions.size())
            outError("User tree file does not have the same number of trees as partitions");
        params.user_file = NULL;
    }

    ModelCheckpoint *model_info = new ModelCheckpoint;
    int part = 0;
    for (auto alnit = super_aln->partitions.begin(); alnit != super_aln->partitions.end(); alnit++, part++) {
        
        checkpoint->startStruct((*alnit)->name);

        // allocate heterotachy tree if neccessary
        int pos = posRateHeterotachy((*alnit)->model_name);
        IQTree *tree;
        
        if (params.num_mixlen > 1) {
            tree = new PhyloTreeMixlen((*alnit), params.num_mixlen);
        } else if (pos != string::npos) {
            tree = new PhyloTreeMixlen((*alnit), 0);
        } else
            tree = new IQTree((*alnit));

        tree->setCheckpoint(checkpoint);
        if (checkpoint->getBool("finished")) {
            tree->restoreCheckpoint();
        } else {
            if (!part_trees.empty())
                tree->copyTree(part_trees[part]);

            startTreeReconstruction(params, tree, *model_info);
            // call main tree reconstruction
            if (params.num_runs == 1)
                runTreeReconstruction(params, tree);
            else
                runMultipleTreeReconstruction(params, tree->aln, tree);
            checkpoint->putBool("finished", true);
            checkpoint->dump();
        }

        super_tree->at(part)->copyTree(tree);
        
        delete tree;
        checkpoint->endStruct();
    }
    
    IQTree *iqtree = super_tree;
    super_tree->setCheckpoint(checkpoint);
    startTreeReconstruction(params, iqtree, *model_info);
    runTreeReconstruction(params, iqtree);
    if (MPIHelper::getInstance().isMaster())
        reportPhyloAnalysis(params, *iqtree, *model_info);

    delete super_tree;
    delete super_aln;
    delete model_info;
}

void assignBranchSupportNew(Params &params) {
    if (!params.user_file)
        outError("No target tree file provided");
    if (params.num_threads == 0)
        outError("-nt AUTO is not supported for concordance factor analysis, please specify no. cores");
    PhyloTree *tree;
    Alignment *aln = NULL;
    if (params.site_concordance) {
        if (!params.aln_file && !params.partition_file)
            outError("Please provide an alignment (-s) or partition file");
        if (params.partition_file) {
            params.compute_seq_composition = false;
            aln = new SuperAlignment(params);
            tree = new PhyloSuperTree((SuperAlignment*)aln);
        } else {
            aln = createAlignment(params.aln_file, params.sequence_type, params.intype, params.model_name);
            tree = new PhyloTree;
        }
    } else {
        tree = new PhyloTree;
    }
    tree->setParams(&params);

    cout << "Reading tree " << params.user_file << " ..." << endl;
    bool rooted = params.is_rooted;
    tree->readTree(params.user_file, rooted);
    cout << ((tree->rooted) ? "rooted" : "un-rooted") << " tree with "
        << tree->leafNum - tree->rooted << " taxa and " << tree->branchNum << " branches" << endl;

    // 2018-12-13: move initialisation to fix rooted vs unrooted tree
    if (params.site_concordance) {
        tree->setAlignment(aln);
        if (tree->isSuperTree())
            ((PhyloSuperTree*)tree)->mapTrees();
    }
    
    BranchVector branches;
    tree->getInnerBranches(branches);
    BranchVector::iterator brit;
    for (brit = branches.begin(); brit != branches.end(); brit++) {
        Neighbor *branch = brit->second->findNeighbor(brit->first);
        string label = brit->second->name;
        if (!label.empty())
            PUT_ATTR(branch, label);
    }
    
    map<string,string> meanings;
    
    if (!params.treeset_file.empty()) {
        MTreeSet trees(params.treeset_file.c_str(), rooted, params.tree_burnin, params.tree_max_count);
        double start_time = getRealTime();
        cout << "Computing gene concordance factor..." << endl;
        tree->computeGeneConcordance(trees, meanings);
        if (params.internode_certainty)
            tree->computeQuartetConcordance(trees);
        cout << getRealTime() - start_time << " sec" << endl;
    }
    if (params.site_concordance) {
        cout << "Computing site concordance factor..." << endl;
        double start_time = getRealTime();
        tree->computeSiteConcordance(meanings);
        cout << getRealTime() - start_time << " sec" << endl;
        delete aln;
    }
    string prefix = (params.out_prefix) ? params.out_prefix : params.user_file;
    string str = prefix + ".cf.tree";
    tree->printTree(str.c_str());
    cout << "Tree with concordance factors written to " << str << endl;
    str = prefix + ".cf.tree.nex";
    string filename = prefix + ".cf.stat";
    tree->printNexus(str, WT_BR_LEN, "See " + filename + " for branch annotation meanings." +
                     " This file is best viewed in FigTree.");
    cout << "Annotated tree (best viewed in FigTree) written to " << str << endl;
    if (verbose_mode >= VB_DEBUG)
        tree->drawTree(cout);
    str = prefix + ".cf.branch";
    tree->printTree(str.c_str(), WT_BR_LEN + WT_INT_NODE + WT_NEWLINE);
    cout << "Tree with branch IDs written to " << str << endl;
    ofstream out;
    out.open(filename.c_str());
    out << "# Concordance factor statistics" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" <<  filename << "',header=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl
        << "#   ID: Branch ID" << endl;
    map<string,string>::iterator mit;
    for (mit = meanings.begin(); mit != meanings.end(); mit++)
        if (mit->first[0] != '*')
            out << "#   " << mit->first << ": " << mit->second << endl;
    out << "#   Label: Existing branch label" << endl;
    out << "#   Length: Branch length" << endl;
    for (mit = meanings.begin(); mit != meanings.end(); mit++)
        if (mit->first[0] == '*')
            out << "# " << mit->first << ": " << mit->second << endl;
    out << "ID";
    for (mit = meanings.begin(); mit != meanings.end(); mit++)
        if (mit->first[0] != '*')
            out << "\t" << mit->first;
    out << "\tLabel\tLength" << endl;
    for (brit = branches.begin(); brit != branches.end(); brit++) {
        Neighbor *branch = brit->second->findNeighbor(brit->first);
        int ID = brit->second->id;
        out << ID;
        for (mit = meanings.begin(); mit != meanings.end(); mit++) {
            if (mit->first[0] == '*')
                continue; // ignore NOTES
            out << '\t';
            string val;
            if (branch->getAttr(mit->first, val))
                out << val;
            else
                out << "NA";
        }
        double length = branch->length;
        string label;
        GET_ATTR(branch, label);
        out << '\t' << label << '\t' << length << endl;
    }
    out.close();
    cout << "Concordance factors per branch printed to " << filename << endl;
    
    if (params.print_cf_quartets) {
        filename = prefix + ".cf.quartet";
        out.open(filename);
        out << "# Site concordance factor for all resampled quartets (with replacement)" << endl
            << "# This file can be read in MS Excel or in R with command:" << endl
            << "#   tab=read.table('" <<  filename << "',header=TRUE)" << endl
            << "# Columns are tab-separated with following meaning:" << endl
            << "#   ID: Branch ID" << endl
            << "#   QuartID: Quartet ID" << endl
            << "#   Seq1: ID of sequence 1 on 'left' side of the branch" << endl
            << "#   Seq2: ID of sequence 2 on 'left' side of the branch" << endl
            << "#   Seq3: ID of sequence 3 on 'right' side of the branch" << endl
            << "#   Seq4: ID of sequence 4 on 'right' side of the branch" << endl
            << "#   qCF: Fraction of concordant sites supporting quartet Seq1,Seq2|Seq3,Seq4 (=qCF_N/qN)" << endl
            << "#   qCF_N: Number of concordant sites supporting quartet Seq1,Seq2|Seq3,Seq4" << endl
            << "#   qDF1: Fraction of discordant sites supporting quartet Seq1,Seq3|Seq2,Seq4  (=qDF1_N/qN)" << endl
            << "#   qDF1_N: Number of discordant sites supporting quartet Seq1,Seq3|Seq2,Seq4" << endl
            << "#   qDF2: Fraction of discordant sites supporting quartet Seq1,Seq4|Seq2,Seq3 (=qDF2_N/qN)" << endl
            << "#   qDF2_N: Number of discordant sites supporting quartet Seq1,Seq4|Seq2,Seq3" << endl
            << "#   qN: Number of decisive sites with four taxa Seq1,Seq2,Seq3,Seq4 (=qCF_N+qDF1_N+qDF2_N)" << endl
            << "ID\tQuartID\tSeq1\tSeq2\tSeq3\tSeq4\tqCF\tqCF_N\tqDF1\tqDF1_N\tqDF2\tqDF2_N\tqN" << endl;
        for (brit = branches.begin(); brit != branches.end(); brit++) {
            Neighbor *branch = brit->second->findNeighbor(brit->first);
            int ID = brit->second->id;
            for (int qid = 0; ; qid++) {
                string qstr;
                if (branch->attributes.find("q" + convertIntToString(qid)) == branch->attributes.end())
                    break;
                out << ID << '\t' << qid+1 << '\t' << branch->attributes["q" + convertIntToString(qid)] << endl;
            }
        }
        out.close();
        cout << "Site concordance factors for quartets printed to " << filename << endl;
    }
    
    if (!params.site_concordance_partition)
        return;
    
    // print concordant/discordant gene trees
    filename = prefix + ".cf.stat_tree";
    out.open(filename);
    out << "# Concordance factor statistics for decisive trees" << endl
    << "# This file can be read in MS Excel or in R with command:" << endl
    << "#   tab2=read.table('" <<  filename << "',header=TRUE)" << endl
    << "# Columns are tab-separated with following meaning:" << endl
    << "#   ID: Branch ID" << endl
    << "#   TreeID: Tree ID" << endl
    << "#   gC: 1/0 if tree is concordant/discordant with branch" << endl
    << "#   gD1: 1/0 if NNI-1 tree is concordant/discordant with branch" << endl
    << "#   gD2: 1/0 if NNI-2 tree is concordant/discordant with branch" << endl
    << "# NOTE: NA means that tree is not decisive for branch" << endl
    << "ID\tTreeID\tgC\tgD1\tgD2" << endl;
    for (brit = branches.begin(); brit != branches.end(); brit++) {
        Neighbor *branch = brit->second->findNeighbor(brit->first);
        int ID = brit->second->id;
        for (int part = 1; ; part++) {
            string gC, gD1, gD2;
            if (!branch->getAttr("gC" + convertIntToString(part), gC))
                break;
            branch->getAttr("gD1" + convertIntToString(part), gD1);
            branch->getAttr("gD2" + convertIntToString(part), gD2);
            out << ID << '\t' << part << '\t' << gC << '\t' << gD1 << '\t' << gD2 << endl;
        }
    }
    out.close();
    cout << "Concordance factors per branch and tree printed to " << filename << endl;
    
    if (!params.site_concordance_partition || !tree->isSuperTree())
        return;
    // print partition-wise concordant/discordant sites
    filename = prefix + ".cf.stat_loci";
    out.open(filename);
    out << "# Concordance factor statistics for loci" << endl
    << "# This file can be read in MS Excel or in R with command:" << endl
    << "#   tab2=read.table('" <<  filename << "',header=TRUE)" << endl
    << "# Columns are tab-separated with following meaning:" << endl
    << "#   ID: Branch ID" << endl
    << "#   PartID: Locus ID" << endl
    << "#   sC: Number of concordant sites averaged over " << params.site_concordance << " quartets" << endl
    << "#   sD1: Number of discordant sites for alternative quartet 1" << endl
    << "#   sD2: Number of discordant sites for alternative quartet 2" << endl
    << "# NOTE: NA means that locus is not decisive for branch" << endl
    << "ID\tPartID\tsC\tsD1\tsD2" << endl;
    for (brit = branches.begin(); brit != branches.end(); brit++) {
        Neighbor *branch = brit->second->findNeighbor(brit->first);
        int ID = brit->second->id;
        for (int part = 1; ; part++) {
            string sC, sD1, sD2;
            if (!branch->getAttr("sC" + convertIntToString(part), sC))
                break;
            if (!branch->getAttr("sD1" + convertIntToString(part), sD1))
                break;
            if (!branch->getAttr("sD2" + convertIntToString(part), sD2))
                break;
            out << ID << '\t' << part << '\t' << sC << '\t' << sD1 << '\t' << sD2 << endl;
        }
    }
    out.close();
    cout << "Concordance factors per branch and locus printed to " << filename << endl;
}



/**
 * assign split occurence frequencies from a set of input trees onto a target tree
 * NOTE: input trees must have the same taxon set
 * @param input_trees file containing NEWICK tree strings
 * @param burnin number of beginning trees to discard
 * @param max_count max number of trees to read in
 * @param target_tree the target tree
 * @param rooted TRUE if trees are rooted, false for unrooted trees
 * @param output_file file name to write output tree with assigned support values
 * @param out_prefix prefix of output file
 * @param mytree (OUT) resulting tree with support values assigned from target_tree
 * @param tree_weight_file file containing INTEGER weights of input trees
 * @param params program parameters
 */
void assignBootstrapSupport(const char *input_trees, int burnin, int max_count,
        const char *target_tree, bool rooted, const char *output_tree,
        const char *out_prefix, MExtTree &mytree, const char* tree_weight_file,
        Params *params) {
    bool myrooted = rooted;
    // read the tree file
    cout << "Reading tree " << target_tree << " ..." << endl;
    mytree.init(target_tree, myrooted);
    if (mytree.rooted)
        cout << "rooted tree detected" << endl;
    else
        cout << "unrooted tree detected" << endl;
    // reindex the taxa in the tree to aphabetical names
    NodeVector taxa;
    mytree.getTaxa(taxa);
    sort(taxa.begin(), taxa.end(), nodenamecmp);
    int i = 0;
    for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++) {
        (*it)->id = i++;
    }

    /*
     string filename = params.boot_trees;
     filename += ".nolen";
     boot_trees.printTrees(filename.c_str(), false);
     return;
     */
    SplitGraph sg;
    SplitIntMap hash_ss;
    // make the taxa name
    vector<string> taxname;
    taxname.resize(mytree.leafNum);
    mytree.getTaxaName(taxname);

    // read the bootstrap tree file
    double scale = 100.0;
    if (params->scaling_factor > 0)
        scale = params->scaling_factor;

    MTreeSet boot_trees;
    if (params && detectInputFile(input_trees) == IN_NEXUS) {
        sg.init(*params);
        for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++)
            hash_ss.insertSplit((*it), (*it)->getWeight());
        StrVector sgtaxname;
        sg.getTaxaName(sgtaxname);
        i = 0;
        for (StrVector::iterator sit = sgtaxname.begin();
                sit != sgtaxname.end(); sit++, i++) {
            Node *leaf = mytree.findLeafName(*sit);
            if (!leaf)
                outError("Tree does not contain taxon ", *sit);
            leaf->id = i;
        }
        scale /= sg.maxWeight();
    } else {
        myrooted = rooted;
        boot_trees.init(input_trees, myrooted, burnin, max_count,
                        tree_weight_file);
        if (mytree.rooted != boot_trees.isRooted())
            outError("Target tree and tree set have different rooting");
        if (boot_trees.equal_taxon_set) {
            boot_trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, params->support_tag);
            scale /= boot_trees.sumTreeWeights();
        }
    }
    //sg.report(cout);
    if (!sg.empty()) {
        cout << "Rescaling split weights by " << scale << endl;
        if (params->scaling_factor < 0)
            sg.scaleWeight(scale, true);
        else {
            sg.scaleWeight(scale, false, params->numeric_precision);
        }

        cout << sg.size() << " splits found" << endl;
    }
    // compute the percentage of appearance
    //    printSplitSet(sg, hash_ss);
    //sg.report(cout);
    cout << "Creating " << RESAMPLE_NAME << " support values..." << endl;
    if (!sg.empty())
        mytree.createBootstrapSupport(taxname, boot_trees, hash_ss, params->support_tag);
    else {
        //mytree.createBootstrapSupport(boot_trees);
        cout << "Unequal taxon sets, rereading trees..." << endl;
        DoubleVector rfdist;
        mytree.computeRFDist(input_trees, rfdist, 1);
    }
    
    //mytree.scaleLength(100.0/boot_trees.size(), true);
    string out_file;
    if (output_tree)
        out_file = output_tree;
    else {
        if (out_prefix)
            out_file = out_prefix;
        else
            out_file = target_tree;
        out_file += ".suptree";
    }

    mytree.printTree(out_file.c_str());
    cout << "Tree with assigned support written to " << out_file
            << endl;
    /*
    if (out_prefix)
        out_file = out_prefix;
    else
        out_file = target_tree;
    out_file += ".supval";
    mytree.writeInternalNodeNames(out_file);

    cout << "Support values written to " << out_file << endl;
    */
}

void computeConsensusTree(const char *input_trees, int burnin, int max_count,
        double cutoff, double weight_threshold, const char *output_tree,
        const char *out_prefix, const char *tree_weight_file, Params *params) {
    bool rooted = false;

    // read the bootstrap tree file
    /*
     MTreeSet boot_trees(input_trees, rooted, burnin, tree_weight_file);
     string first_taxname = boot_trees.front()->root->name;
     //if (params.root) first_taxname = params.root;

     SplitGraph sg;

     boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);*/

    //sg.report(cout);
    SplitGraph sg;
    SplitIntMap hash_ss;
    // make the taxa name
    //vector<string> taxname;
    //taxname.resize(mytree.leafNum);
    //mytree.getTaxaName(taxname);

    // read the bootstrap tree file
    double scale = 100.0;
    if (params->scaling_factor > 0)
        scale = params->scaling_factor;

    MTreeSet boot_trees;
    if (params && detectInputFile(input_trees) == IN_NEXUS) {
        char *user_file = params->user_file;
        params->user_file = (char*) input_trees;
        params->split_weight_summary = SW_COUNT; // count number of splits
        sg.init(*params);
        params->user_file = user_file;
        for (SplitGraph::iterator it = sg.begin(); it != sg.end();)
            if ((*it)->getWeight() > weight_threshold) {
                hash_ss.insertSplit((*it), (*it)->getWeight());
                it++;
            } else {
                // delete the split
                if (it != sg.end()-1) {
                    *(*it) = (*sg.back());
                }
                delete sg.back();
                sg.pop_back();
            }
        /*        StrVector sgtaxname;
         sg.getTaxaName(sgtaxname);
         i = 0;
         for (StrVector::iterator sit = sgtaxname.begin(); sit != sgtaxname.end(); sit++, i++) {
         Node *leaf = mytree.findLeafName(*sit);
         if (!leaf) outError("Tree does not contain taxon ", *sit);
         leaf->id = i;
         }*/
        scale /= sg.maxWeight();
    } else {
        boot_trees.init(input_trees, rooted, burnin, max_count,
                tree_weight_file);
        boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);
        scale /= boot_trees.sumTreeWeights();
        cout << sg.size() << " splits found" << endl;
    }
    //sg.report(cout);
    if (verbose_mode >= VB_MED)
        cout << "Rescaling split weights by " << scale << endl;
    if (params->scaling_factor < 0)
        sg.scaleWeight(scale, true);
    else {
        sg.scaleWeight(scale, false, params->numeric_precision);
    }



    //cout << "Creating greedy consensus tree..." << endl;
    MTree mytree;
    SplitGraph maxsg;
    sg.findMaxCompatibleSplits(maxsg);

    if (verbose_mode >= VB_MAX)
        maxsg.saveFileStarDot(cout);
    //cout << "convert compatible split system into tree..." << endl;
    mytree.convertToTree(maxsg);
    //cout << "done" << endl;
    if (!mytree.rooted) {
        string taxname;
        if (params->root)
            taxname = params->root;
        else
            taxname = sg.getTaxa()->GetTaxonLabel(0);
        Node *node = mytree.findLeafName(taxname);
        if (node)
            mytree.root = node;
    }
    // mytree.scaleLength(100.0 / boot_trees.sumTreeWeights(), true);

    // mytree.getTaxaID(maxsg.getSplitsBlock()->getCycle());
    //maxsg.saveFile(cout);

    string out_file;

    if (output_tree)
        out_file = output_tree;
    else {
        if (out_prefix)
            out_file = out_prefix;
        else
            out_file = input_trees;
        out_file += ".contree";
    }

//    if (removed_seqs.size() > 0)
//        mytree.insertTaxa(removed_seqs, twin_seqs);

    mytree.printTree(out_file.c_str(), WT_BR_CLADE);
    cout << "Consensus tree written to " << out_file << endl;

    if (output_tree)
        out_file = output_tree;
    else {
        if (out_prefix)
            out_file = out_prefix;
        else
            out_file = input_trees;
        out_file += ".splits";
    }

    //sg.scaleWeight(0.01, false, 4);
    if (params->print_splits_file) {
        sg.saveFile(out_file.c_str(), IN_OTHER, true);
        cout << "Non-trivial split supports printed to star-dot file " << out_file << endl;
    }

}

void computeConsensusNetwork(const char *input_trees, int burnin, int max_count,
        double cutoff, int weight_summary, double weight_threshold, const char *output_tree,
        const char *out_prefix, const char* tree_weight_file) {
    bool rooted = false;

    // read the bootstrap tree file
    MTreeSet boot_trees(input_trees, rooted, burnin, max_count,
            tree_weight_file);

    SplitGraph sg;
    //SplitIntMap hash_ss;

    boot_trees.convertSplits(sg, cutoff, weight_summary, weight_threshold);

    string out_file;

    if (output_tree)
        out_file = output_tree;
    else {
        if (out_prefix)
            out_file = out_prefix;
        else
            out_file = input_trees;
        out_file += ".nex";
    }

    sg.saveFile(out_file.c_str(), IN_NEXUS);
    cout << "Consensus network printed to " << out_file << endl;

    if (output_tree)
        out_file = output_tree;
    else {
        if (out_prefix)
            out_file = out_prefix;
        else
            out_file = input_trees;
        out_file += ".splits";
    }
    if (verbose_mode >= VB_MED) {
        sg.saveFile(out_file.c_str(), IN_OTHER, true);
        cout << "Non-trivial split supports printed to star-dot file " << out_file << endl;
    }

}
