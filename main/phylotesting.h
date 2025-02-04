/*
 * phylotesting.h
 *
 *  Created on: Aug 23, 2013
 *      Author: minh
 */

#ifndef PHYLOTESTING_H_
#define PHYLOTESTING_H_

#ifdef _IQTREE_MPI
    #include <mpi.h>
#endif

#include "utils/tools.h"
#include "utils/checkpoint.h"
#include "nclextra/modelsblock.h"
#include "alignment/superalignment.h"
#include "utils/MPIHelper.h"
#include <set>

class PhyloTree;
class IQTree;
class ModelCheckpoint;
class SyncChkPoint;
class PhyloSuperTree;
class SubsetPair;

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
        syncChkPoint = nullptr;
        //init_first_mix = false;
        model_selection_action = 0;
    }
    
    CandidateModel(string subst_name, string rate_name, Alignment *aln, int flag = 0) : CandidateModel(flag) {
        this->subst_name = orig_subst_name = subst_name;
        this->rate_name = orig_rate_name = rate_name;
        this->aln = aln;
        syncChkPoint = nullptr;
        //init_first_mix = false;
        model_selection_action = 0;
    }
    
    CandidateModel(Alignment *aln, int flag = 0) : CandidateModel(flag) {
        this->aln = aln;
        getUsualModel(aln);
        syncChkPoint = nullptr;
        //init_first_mix = false;
        model_selection_action = 0;
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
    void setFlag(int flag) {
        this->flag |= flag;
    }

    bool hasFlag(int flag) {
        return (this->flag & flag) != 0;
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
    
    // indicate whether it is the first k-class mixture model
    // if so, then the parameters will be initialized from the previous (k-1)-class mixture model
    //bool init_first_mix;

    /** the nest relationships of all candidate Q matrices */
    map<string, vector<string> > nest_network;

    /** the value of the action in function runModelSelection*/
    int model_selection_action;

    Alignment *aln; // associated alignment

    /**
     Synchronization of check point for MPI
     */
    SyncChkPoint* syncChkPoint;

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
        syncChkPoint = nullptr;
        under_mix_finder = false;
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
     Filter out all "non-promissing" rate models
     */
    void filterRates(int finished_model);

    /**
     Filter out all "non-promissing" substitution models
     */
    void filterSubst(int finished_model);
    
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
     @param generate_candidates true to generate candidates in the beginning
     @param skip_all_when_drop true to skip the testing of all the subsequent models when the current model becomes less favorable
     @return name of best-fit-model
     */
    CandidateModel test(Params &params, PhyloTree* in_tree, ModelCheckpoint &model_info,
                ModelsBlock *models_block, int num_threads, int brlen_type,
                string set_name = "", string in_model_name = "",
                bool merge_phase = false,
                bool generate_candidates = true,
                bool skip_all_when_drop = false);

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

    /**
     Synchronization of check point for MPI
     */
    SyncChkPoint* syncChkPoint;

    /** the nest relationships of all candidate Q matrices */
    map<string, vector<string> > nest_network;

    /** whether it is under the process of mixture finder */
    bool under_mix_finder;
    
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

/** model information by merging two partitions */
struct ModelPair {
    /** score after merging */
    double score;
    /** ID of partition 1 */
    int part1;
    /** ID of partition 2 */
    int part2;
    /** log-likelihood */
    double logl;
    /** degree of freedom */
    int df;
    /** tree length */
    double tree_len;
    /** IDs of merged partitions */
    set<int> merged_set;
    /** set name */
    string set_name;
    /* best model name */
    string model_name;
};

class ModelPairSet : public multimap<double, ModelPair> {

public:

    /** insert a partition pair */
    void insertPair(ModelPair &pair) {
        insert(value_type(pair.score, pair));
    }

    /**
        find the maximum compatible partition pairs
        @param num max number of pairs to return
    */
    void getCompatiblePairs(int num, ModelPairSet &res) {
        set<int> part_ids;

        for (auto it = begin(); it != end() && res.size() < num; it++) {

            // check for compatibility
            vector<int> overlap;
            set_intersection(part_ids.begin(), part_ids.end(),
                it->second.merged_set.begin(), it->second.merged_set.end(),
                std::back_inserter(overlap));

            if (!overlap.empty()) continue;

            // take the union
            part_ids.insert(it->second.merged_set.begin(), it->second.merged_set.end());

            // put the compatible pair to the set
            res.insertPair(it->second);
        }
    }

};

#ifdef _IQTREE_MPI

/*
 * This class is designed for a job to perform merging between two partitions (for MPI)
 */
class MergeJob {
public:
    int id1;
    int id2;
    set<int> geneset1;
    set<int> geneset2;
    double treelen1;
    double treelen2;
    // constructors
    MergeJob();
    MergeJob(int id_1, int id_2, set<int>& geneset_1, set<int>& geneset_2, double treelen_1, double treelen_2);
    bool isEmpty();
    void copyFrom(MergeJob* anotherMergeJob);
    void setEmpty();
    void toString(string& str);
    void loadString(string& str);
};

#endif

/*
 * This class is designed for partition finder
 */
class PartitionFinder {

private:
    int brlen_type;
    bool test_merge;
    SuperAlignment *super_aln;


    // retreive the answers from checkpoint
    // and remove those jobs from the array jobIDs
    void retreiveAnsFrChkpt(vector<pair<int,double> >& jobs, int job_type);

    /**
     * compute and process the best model for partitions (without MPI)
     * nthreads : the number of threads available for these jobs
     */
    void getBestModelforPartitionsNoMPI(int nthreads, vector<pair<int,double> >& jobs);

    /**
     * compute and process the best model for merges (without MPI)
     * nthreads : the number of threads available for these jobs
     */
    void getBestModelforMergesNoMPI(int nthreads, vector<pair<int,double> >& jobs);

    /**
     * compute the best model
     * job_type = 1 : for all partitions
     * job_type = 2 : for all merges
     */
    void getBestModel(int job_type);

#ifdef _IQTREE_MPI
    
    // The following functions are for MPI

    /**
     * Process the computation of the best model for a single partition with MPI
     *
     * nthreads : number of threads available for this job
     * need_next_treeID : whether it is needed to get the next tree ID
     *
     * if need_next_treeID, then
     *    if WORKER and IS_ASYN_COMM = 1 (i.e. asynchronous communication)
     *        return the index of the array storing MPI_Request
     *    else
     *        return the next Job ID from master
     * else
     *    return -1
     */
    int computeBestModelforOnePartitionMPI(int tree_id, int nthreads, bool need_next_treeID, SyncChkPoint& syncChkPt, double& run_time, double& wait_time);

    /**
     * Process the computation of the best model for a merge with MPI
     *
     * nthreads : number of threads available for this job
     * need_next_job : whether it is needed to get the next job
     *
     * if need_next_job
     *    job will be updated to the next job
     */
    void getBestModelForOneMergeMPI(MergeJob* job, int nthreads, bool need_next_job, SyncChkPoint& syncChkPt, double& run_time, double& wait_time);

    /**
	 * compute and process the best model for partitions (for MPI)
	 */
	void getBestModelforPartitionsMPI(int nthreads, vector<int> &jobs, double* run_time, double* wait_time, double* fstep_time, int* partNum, double& cpu_time, double& wall_time);

	/**
	 * compute and process the best model for merges (for MPI)
	 */
	void getBestModelforMergesMPI(int nthreads, vector<MergeJob* >& jobs, double* run_time, double* wait_time, double* fstep_time, int* partNum, double& cpu_time, double& wall_time);
    
    /*
     * Consolidate the partition results (for MPI)
     */
    void consolidPartitionResults();
    
    /*
     * Consolidate the merge results (for MPI)
     */
    void consolidMergeResults();

#endif

public:
    ModelCheckpoint *model_info;
    DoubleVector lhvec; // log-likelihood for each partition
    IntVector dfvec; // number of parameters for each partition
    DoubleVector lenvec; // tree length for each partition
    double lhsum;
    int dfsum;
    double start_time;
    int64_t total_num_model;
    int64_t num_model;
    vector<SubsetPair> closest_pairs;
    vector<set<int> > gene_sets;
    PhyloSuperTree* in_tree;
    size_t  ssize;
    Params *params;
    double inf_score;
    ModelPairSet better_pairs; // list of all better pairs of partitions than current partitioning scheme

    ModelsBlock *models_block;
    int num_threads;
    int num_processes;

    int nextjob;
    int jobdone;
    int tot_job_num;
    vector<int> remain_job_list;
    
#ifdef _IQTREE_MPI
    vector<MergeJob*> remain_mergejobs;
#endif
    
    int base;
    
    // for ONE-SIDE communication
    double last_syn_time;
    vector<int> tree_id_vec;
    vector<double> tree_len_vec;
    vector<string> model_name_vec;
    vector<double> score_vec;
    vector<string> set_name_vec;
    vector<int> tag_vec;
    int tot_jobs_done;
    ModelCheckpoint process_model_info;

#ifdef _IQTREE_MPI
    // shared memory space to be accessed by the other processors
    int   *val_ptr;
    MPI_Win   win;
#endif

    /* Constructor
     */
    PartitionFinder(Params *inparams, PhyloSuperTree* intree, ModelCheckpoint *modelinfo,
                    ModelsBlock *modelsblock, int numthreads);

    /* Destructor
     */
    ~PartitionFinder();
        
    /*
     * Perform the computation
     */
    void test_PartitionModel();

#ifdef _IQTREE_MPI
    /*
     *  initialize the shared memory space to be accessed by the other processors
     */
    void initialMPIShareMemory();

    /*
     *  free the shared memory space
     */
    void freeMPIShareMemory();
    
    /*
     * assign initial partition jobs to processors
     * input: a set of partition jobs ordered by the estimated computational costs
     * output: number of items in currJobs
     */
    int partjobAssignment(vector<pair<int,double> > &job_ids, vector<int> &currJobs);
    
    /*
     * assign initial merge jobs to processors
     * input: a set of merge jobs ordered by the estimated computational costs
     * output: number of items in currJobs
     */
    int mergejobAssignment(vector<pair<int,double> > &job_ids, vector<MergeJob* >&currJobs);

#endif // _IQTREE_MPI

    /*
     * Show the result of best model for the partition
     */
    void showPartitionResult(ModelCheckpoint& part_model_info, int tree_id, double tree_len, const string& model_name, double score, int tag);

    /*
     * Show the result of best model for the partition
     */
    void showPartitionResults(ModelCheckpoint& part_model_info, vector<int>& tree_id, vector<double>& tree_len, vector<string>& model_name, vector<double>& score, vector<int>& tag);
    
    /*
     * Show the the other worker's result of best model for the merge
     */
    void showMergeResult(ModelCheckpoint& part_model_info, double tree_len, const string& model_name, string& set_name, bool done_before, int tag);

    /*
     * Show the the other worker's result of best model for the merge
     */
    void showMergeResults(ModelCheckpoint& part_model_info, vector<double>& tree_len, vector<string>& model_name, vector<string>& set_name, vector<int>& tag, int tot_jobs_done);
};

#ifdef _IQTREE_MPI

/*
 * This class is designed for synchronization of checkpoints for partition finder (for MPI)
 */
class SyncChkPoint {

private:
    // shared among threads
    PartitionFinder* pfinder;

public:

    int mytag;

    /*  constructor
     */
    SyncChkPoint(PartitionFinder* pf, int thres_id);
    
    /*
     * Show the result of best model
     */
    void showResult(ModelCheckpoint& part_model_info, int tag);

    /*
     * FOR MASTER - synchronize the checkpoints from the other processors
     * Receive checkpoint from worker and send the next Job ID to workers
     * increase the value of next_job and job_done by 1
     * update the master's checkpoint: model_info
     */
    void masterSyncOtherChkpts(bool chk_gotMessage = true);
    
    /*
     * FOR WORKER
     * send checkpoint to master
     * clear the checkpoint
     *
     * if need_nextJobID and NOT ASYN_COMM (i.e. non-asynchronous communication)
     * return the next Job ID from master
     * else -1
     */
    int sendChkptToMaster(ModelCheckpoint &model_info, bool need_nextJobID, int job_type, MergeJob* mergeJob = nullptr, bool forceToSyn = false);

    /*
     * receive an integer from the master (for synchronous communication)
     */
    // int recvInt(int tag);

    /*
     * get the next Job ID
     */
    int getNextJobID();

    /*
     * get the next Merge Job
     */
    void getNextMergeJob(MergeJob* mergejob);

    void sendCheckpoint(Checkpoint *ckp, int dest, int tag);
    
    // void recvCheckpoint(Checkpoint *ckp, int src, int tag);
    
    void recvAnyCheckpoint(Checkpoint *ckp, int& src, int& tag);
    
    void recvAnyString(string &str, int& src, int& tag);

    /*
     * Check for incoming messages
     * if there is a message, collect the tag value and the source
     */
    bool gotMessage(int& tag, int& source);
    
    void sendMergeJobToWorker(MergeJob& mergeJob, int dest, int tag);
    
    void recMergeJobFrMaster(MergeJob& mergeJob, int tag);
    
    void broadcastVecSetInt(vector<set<int> >& gene_sets);
    
    void broadcastVecStr(vector<string>& model_names);
    
    int* toIntArr(vector<set<int> >& gene_sets, int& buffsize);
    
    void loadFrIntArr(vector<set<int> >& gene_sets, int* buff, int buffsize);
    
    char* toCharArr(vector<string>& model_names, int& buffsize);
    
    void loadFrCharArr(vector<string>& model_names, char* buff);
};

#endif

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
 @param best_subst_name (OUT) information for all models considered
 @param best_rate_name (OUT) information for all models considered
 @param nest_network (IN) nest relationships of all DNA models considered
 @param under_mix_finder (IN) whether MixtureFinder is being used
 */
void runModelFinder(Params &params, IQTree &iqtree, ModelCheckpoint &model_info, string &best_subst_name, string &best_rate_name, map<string, vector<string> > nest_network, bool under_mix_finder = false);

/**
 optimisation of Q-Mixture model, including estimation of best number of classes in the mixture
 @param params program parameters
 @param iqtree phylogenetic tree
 @param model_info (IN/OUT) information for all models considered
 */
void optimiseQMixModel(Params &params, IQTree* &iqtree, ModelCheckpoint &model_info);

/**
 perform ModelFinderNN to find the best-fit model (uses neural network for model inference)
 @param params program parameters
 @param iqtree phylogenetic tree
 @param model_info (IN/OUT) information for all models considered
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

string convertSeqTypeToSeqTypeName(SeqType seq_type);

string detectSeqTypeName(string model_name);

/****************************************************/
/*    Q MATRICES NESTING CHECK                      */
/****************************************************/

/**
 * get the index of a DNA model in dna_model_names
 */
int findModelIndex(const string& model, const char* model_set[], size_t size);

/**
 * reorder the input dna models or RHAS models as default (from simple to flexible)
 */
void reorderModelNames(StrVector& model_names, const char* model_set[], size_t size);

/**
 * check whether rate_type2 is nested in rate_type1
 */
bool isRateTypeNested(string rate_type1, string rate_type2);

/**
 * build the nest relationships of all candidate Q matrices
 */
map<string, vector<string> > generateNestNetwork(StrVector model_names, StrVector freq_names);

#endif /* PHYLOTESTING_H_ */
