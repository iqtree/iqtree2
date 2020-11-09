/*
 * phylotesting.h
 *
 *  Created on: Aug 23, 2013
 *      Author: minh
 */

#ifndef PHYLOTESTING_H_
#define PHYLOTESTING_H_

#include "utils/tools.h"
#include "utils/checkpoint.h"
#include "nclextra/modelsblock.h"
#include "alignment/superalignment.h"

class PhyloTree;
class IQTree;
class ModelCheckpoint;

const int MF_SAMPLE_SIZE_TRIPLE = 1;
const int MF_IGNORED            = 2;
const int MF_RUNNING            = 4;
const int MF_WAITING            = 8;
const int MF_DONE               = 16;

/**
    Candidate model under testing
 */
class CandidateModel {
    
public:
    
    /** constructor */
    CandidateModel(int flag = 0) {
        logl = 0.0;
        df = 0;
        tree_len = 0.0;
        aln = NULL;
        AIC_score = DBL_MAX;
        AICc_score = DBL_MAX;
        BIC_score = DBL_MAX;
        this->flag = flag;
    }
    
    CandidateModel(string subst_name, string rate_name, Alignment *aln, int flag = 0) : CandidateModel(flag) {
        this->subst_name = orig_subst_name = subst_name;
        this->rate_name = orig_rate_name = rate_name;
        this->aln = aln;
    }
    
    CandidateModel(Alignment *aln, int flag = 0) : CandidateModel(flag) {
        this->aln = aln;
        getUsualModel(aln);
    }
    
    string getName() {
        return subst_name + rate_name;
    }
    
    /**
     get usual model for a given alignment
     @param aln input alignment
     @return length of the alignment
     */
    size_t getUsualModel(Alignment *aln);
    
    /**
     evaluate this model
     @param params program parameters
     @param in_aln input alignment
     @param[in] in_model_info input checkpointing information
     @param[out] out_model_info output checkpointing information
     @param models_block models block
     @param num_thread number of threads
     @param brlen_type BRLEN_OPTIMIZE | BRLEN_FIX | BRLEN_SCALE | TOPO_UNLINKED
     @return tree string
     */
    string evaluate(Params &params,
                    ModelCheckpoint &in_model_info, ModelCheckpoint &out_model_info,
                    ModelsBlock *models_block, int &num_threads, int brlen_type);
    
    /**
     evaluate concatenated alignment
     */
    string evaluateConcatenation(Params &params, SuperAlignment *super_aln,
                                 ModelCheckpoint &model_info, ModelsBlock *models_block, int num_threads);

    /**
     compute information criterion scores (AIC, AICc, BIC)
     */
    void computeICScores(size_t sample_size);
    void computeICScores();

    /**
     compute information criterion scores (AIC, AICc, BIC)
     */
    double computeICScore(size_t sample_size);
    
    /** @return model score */
    double getScore();

    /** @return model score */
    double getScore(ModelTestCriterion mtc);

    /**
     save model into checkpoint
     */
    void saveCheckpoint(Checkpoint *ckp) {
        stringstream ostr;
        ostr.precision(10);
        ostr << logl << " " << df << " " << tree_len;
        if (!tree.empty())
            ostr << " " << tree;
        ckp->put(getName(), ostr.str());
    }
    
    /**
     restore model from checkpoint
     */
    bool restoreCheckpoint(Checkpoint *ckp) {
        string val;
        if (ckp->getString(getName(), val)) {
            stringstream str(val);
            str >> logl >> df >> tree_len;
            return true;
        }
        return false;
    }
    
    /**
     restore model from checkpoint
     */
    bool restoreCheckpointRminus1(Checkpoint *ckp, CandidateModel *model) {
        size_t posR;
        const char *rates[] = {"+R", "*R", "+H", "*H"};
        for (int i = 0; i < sizeof(rates)/sizeof(char*); i++) {
            if ((posR = model->rate_name.find(rates[i])) != string::npos) {
                int cat = convert_int(model->rate_name.substr(posR+2).c_str());
                subst_name = model->subst_name;
                rate_name = model->rate_name.substr(0, posR+2) + convertIntToString(cat-1);
                return restoreCheckpoint(ckp);
            }
        }
        return false;
    }
    
    /** turn on some flag with OR operator */
    bool setFlag(int flagsToSet) {
        bool flagChanged = (this->flag & flagsToSet) != flagsToSet;
        this->flag |= flagsToSet;
        return flagChanged;
    }

    bool hasFlag(int flag_to_check) {
        return (this->flag & flag_to_check) != 0;
    }
    
    string set_name; // subset name
    string subst_name; // substitution matrix name
    string orig_subst_name; // original substitution name
    string rate_name; // rate heterogeneity name
    string orig_rate_name; // original rate heterogeneity name
    double logl; // tree log likelihood
    int df;      // #parameters
    double tree_len; // tree length, added 2015-06-24 for rcluster algorithm
    string tree; // added 2015-04-28: tree string
    double AIC_score, AICc_score, BIC_score;    // scores
    double AIC_weight, AICc_weight, BIC_weight; // weights
    bool AIC_conf, AICc_conf, BIC_conf;         // in confidence set?

    Alignment *aln; // associated alignment
    
protected:
    
    /** flag */
    int flag;
};

/**
 set of candidate models
 */
class CandidateModelSet : public vector<CandidateModel> {
public:

    CandidateModelSet() : vector<CandidateModel>() {
        current_model = -1;
    }
    
    /** get ID of the best model */
    int getBestModelID(ModelTestCriterion mtc);
    
    /**
     * get the list of model
     * @param params program parameters
     * @param aln alignment
     * param separate_rate true to separate rates from models
     * @param merge_phase true to consider models for merging phase
     * @return maximum number of rate categories
     */
    int generate(Params &params, Alignment *aln, bool separate_rate, bool merge_phase);

    /**
     Filter out all "non-promising" rate models
                @return the number of "non-promising" models that were filtered out.
     */
    int filterRates(int finished_model);

    /**
     Filter out all "non-promising" substitution models
                @return the number of "non-promising" models that were filtered out.
     */
    int filterSubst(int finished_model);

    /**
     testing the best-fit model
     return in params.freq_type and params.rate_type
     @param params global program parameters
     @param in_tree phylogenetic tree
     @param model_info (IN/OUT) information for all models considered
     @param models_block global model definition
     @param num_threads number of threads
     @param brlen_type BRLEN_OPTIMIZE | BRLEN_FIX | BRLEN_SCALE | TOPO_UNLINK
     @param set_name for partition model selection
     @param in_model_name a specific model name if testing one model
     @param adjust model adjustment for modelomatic
     @param merge_phase true to consider models for merging phase
     @return name of best-fit-model
     */
    CandidateModel test(Params &params, PhyloTree* in_tree, ModelCheckpoint &model_info,
                ModelsBlock *models_block, int num_threads, int brlen_type,
                string set_name = "", string in_model_name = "",
                bool merge_phase = false);

    /**
     for a rate model XXX+R[k], return XXX+R[k-j] that finished
     @return the index of fewer category +R model that finished
     */
    int getLowerKModel(int model) {
        size_t posR;
        const char *rates[] = {"+R", "*R", "+H", "*H"};
        for (int i = 0; i < sizeof(rates)/sizeof(char*); i++) {
            if ((posR = at(model).rate_name.find(rates[i])) == string::npos)
                continue;
            int cat = convert_int(at(model).rate_name.substr(posR+2).c_str());
            for (int prev_model = model-1; prev_model >= 0; prev_model--, cat--) {
                string name = at(model).rate_name.substr(0, posR+2) + convertIntToString(cat-1);
                if (at(prev_model).rate_name != name)
                    break;
                if (!at(prev_model).hasFlag(MF_DONE))
                    continue;
                return prev_model;
            }
        }
        return -1;
    }

    int getHigherKModel(int model) {
        size_t posR;
        const char *rates[] = {"+R", "*R", "+H", "*H"};
        for (int i = 0; i < sizeof(rates)/sizeof(char*); i++) {
            if ((posR = at(model).rate_name.find(rates[i])) == string::npos)
                continue;
            size_t this_posR = at(model).rate_name.find(rates[i]);
            ASSERT(this_posR != string::npos);
            int cat = convert_int(at(model).rate_name.substr(this_posR+2).c_str());
            for (int next_model = model+1; next_model < size(); next_model++, cat++) {
//                if (at(next_model).name.substr(0, posR) != orig_name.substr(0, posR))
//                    break;
                string rate_name = at(model).rate_name.substr(posR, 2) + convertIntToString(cat+1);
                if (at(next_model).rate_name.find(rate_name) == string::npos)
                    break;
                return next_model;
            }
        }
        return -1;
    }

    /** get the next model to evaluate in parallel */
    int64_t getNextModel();

    /**
     evaluate all models in parallel
     */
    CandidateModel evaluateAll(Params &params, PhyloTree* in_tree, ModelCheckpoint &model_info,
                     ModelsBlock *models_block, int num_threads, int brlen_type,
                     string in_model_name = "", bool merge_phase = false, bool write_info = true);
    
private:
    
    /** current model */
    int64_t current_model;
};

//typedef vector<ModelInfo> ModelCheckpoint;

class ModelCheckpoint : public Checkpoint {

public:

    /*
        get the best model
        @param[out] best_model name of the best model
        @return TRUE if best model found, FALSE otherwise (unfinished job)
    */
    bool getBestModel(string &best_model);

    /*
     get the best model list
     @param[out] best_model_list list of the best model
     @return TRUE if best model found, FALSE otherwise (unfinished job)
     */
    bool getBestModelList(string &best_model_list);

    /*
     put the best model list
     @param best_model_list list of the best model
     @return TRUE if best model found, FALSE otherwise (unfinished job)
     */
    void putBestModelList(string &best_model_list);

    /*
        get the ordered model list according to AIC, AICc or BIC
        @param tree associated tree
        @param[out] ordered_models list of models ordered by specified criterion
        @return TRUE if ordered_models found, FALSE otherwise (unfinished job)
    */
    bool getOrderedModels(PhyloTree *tree, CandidateModelSet &ordered_models);

    /*
        get the best tree
        @param[out] best_tree NEWICK string of the best tree
        @return TRUE if best tree found, FALSE otherwise (unfinished job)
    */
    bool getBestTree(string &best_tree);

};

/**
 * computing AIC, AICc, and BIC scores
 */
void computeInformationScores(double tree_lh, int df, int ssize, double &AIC, double &AICc, double &BIC);

double computeInformationScore(double tree_lh, int df, int ssize, ModelTestCriterion mtc);

string criterionName(ModelTestCriterion mtc);

/**
 perform ModelFinder to find the best-fit model
 @param params program parameters
 @param iqtree phylogenetic tree
 @param model_info (IN/OUT) information for all models considered
 */
void runModelFinder(Params &params, IQTree &iqtree, ModelCheckpoint &model_info);

/**
 testing the best-fit model
 return in params.freq_type and params.rate_type
 @param set_name for partitioned analysis
 @param in_tree phylogenetic tree
 @param model_info (IN/OUT) information for all models considered
 @param set_name for partition model selection
 @param print_mem_usage true to print RAM memory used (default: false) 
 @return name of best-fit-model
 */
//string testModel(Params &params, PhyloTree* in_tree, ModelCheckpoint &model_info,
//		ModelsBlock *models_block, int num_threads, int brlen_type,
//        string set_name = "", bool print_mem_usage = false, string in_model_name = "");


/**
    get sequence type for a model name
    @param model_name model name string
    @param seq_type (OUT) sequence type, SEQ_UNKNOWN if is not determined
    @return 1 for parametric model, 2 for empirical model
*/
int detectSeqType(const char *model_name, SeqType &seq_type);

string detectSeqTypeName(string model_name);


#endif /* PHYLOTESTING_H_ */
