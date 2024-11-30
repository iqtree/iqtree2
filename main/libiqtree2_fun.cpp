#include "libiqtree2_fun.h"

#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
#include <winsock2.h>
#endif

class input_options {
public:
    vector<string> flags;
    vector<string> values;

    void insert(string flag, string value="") {
        flags.push_back(flag);
        values.push_back(value);
    }
    
    // set Params according to the input options from PiQTREE
    // only invoke this function after the default values of parameters are set
    // this function defines which IQ-TREE options are availble for PiQTREE
    void set_params(Params& params);
};

void cleanup(Params& params) {
    if (params.state_freq_set != NULL) {
        delete[] params.state_freq_set;
        params.state_freq_set = NULL;
    }
}

string build_phylogenetic(vector<string>& names, vector<string>& seqs, string model, string intree,
                          int rand_seed, string prog, input_options* in_options);

// Calculates the robinson fould distance between two trees
int robinson_fould(const string& tree1, const string& tree2) {

    int output;

    try {
        MTree first_tree;
        bool is_rooted = false;
        std::vector<double> rfdist;

        // read in the first tree
        first_tree.read_TreeString(tree1, is_rooted);
        
        // second tree
        stringstream second_tree_str;
        second_tree_str << tree2;
        second_tree_str.seekg(0, ios::beg);
        
        // compute the RF distance
        first_tree.computeRFDist(second_tree_str, rfdist);
        
        output = (int)rfdist[0];
    } catch (std::runtime_error& e) {
        // reset the output and error buffers
        funcExit();
        throw e;
    }

    return output;
}

// Generates a set of random phylogenetic trees
// tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
string random_tree(int num_taxa, string tree_gen_mode, int num_trees, int rand_seed) {
    string output;
    
    try {
        PhyloTree ptree;
        int seed = rand_seed;
        if (seed == 0)
            seed = make_new_seed();
        cout << "seed: " << seed << endl;
        init_random(seed);
        
        TreeGenType tree_mode;
        if (tree_gen_mode == "YULE_HARDING") {
            tree_mode = YULE_HARDING;
        } else if (tree_gen_mode == "UNIFORM") {
            tree_mode = UNIFORM;
        } else if (tree_gen_mode == "CATERPILLAR") {
            tree_mode = CATERPILLAR;
        } else if (tree_gen_mode == "BALANCED") {
            tree_mode = BALANCED;
        } else if (tree_gen_mode == "BIRTH_DEATH") {
            tree_mode = BIRTH_DEATH;
        } else if (tree_gen_mode == "STAR_TREE") {
            tree_mode = STAR_TREE;
        } else {
            outError("Unknown mode: " + tree_gen_mode);
        }
        
        Params params = Params::getInstance();
        params.setDefault();
        params.sub_size = num_taxa;
        params.tree_gen = tree_mode;
        params.repeated_time = num_trees;
        params.ignore_checkpoint = true; // overrid the output file if exists
        params.user_file = (char*) "";
        
        ostringstream ostring;
        generateRandomTree(params, ostring);
        output = ostring.str();
    } catch (std::runtime_error& e) {
        // reset the output and error buffers
        funcExit();
        throw e;
    }

    return output;
}

// Perform phylogenetic analysis on the input alignment (in string format)
// With estimation of the best topology
// num_thres -- number of cpu threads to be used, default: 1
string build_tree(vector<string>& names, vector<string>& seqs, string model, int rand_seed, int bootstrap_rep, int num_thres) {
    string intree = "";
    string output;
    try {
        input_options* in_options = NULL;
        if (bootstrap_rep > 0 || num_thres > 1) {
            in_options = new input_options();
            if (bootstrap_rep > 0)
                in_options->insert("-bb", convertIntToString(bootstrap_rep));
            if (num_thres > 1)
                in_options->insert("-nt", convertIntToString(num_thres));
        }
        output = build_phylogenetic(names, seqs, model, intree, rand_seed, "build_tree", in_options);
        if (in_options != NULL)
            delete in_options;
    } catch (std::runtime_error& e) {
        // reset the output and error buffers
        funcExit();
        throw e;
    }
    return output;
}

// Perform phylogenetic analysis on the input alignment (in string format)
// With restriction to the input toplogy
// num_thres -- number of cpu threads to be used, default: 1
string fit_tree(vector<string>& names, vector<string>& seqs, string model, string intree, int rand_seed, int num_thres) {
    string output;
    try {
        input_options* in_options = NULL;
        if (num_thres > 1) {
            in_options = new input_options();
            in_options->insert("-nt", convertIntToString(num_thres));
        }
        output = build_phylogenetic(names, seqs, model, intree, rand_seed, "fit_tree", in_options);
        if (in_options != NULL)
            delete in_options;
    } catch (std::runtime_error& e) {
        // reset the output and error buffers
        funcExit();
        throw e;
    }
    return output;
}

// Perform phylogenetic analysis with ModelFinder
// on the input alignment (in string format)
// model_set -- a set of models to consider
// freq_set -- a set of frequency types
// rate_set -- a set of RHAS models
// num_thres -- number of cpu threads to be used, default: 1
string modelfinder(vector<string>& names, vector<string>& seqs, int rand_seed, string model_set, string freq_set, string rate_set, int num_thres) {
    
    input_options* in_options = NULL;
    string output;
    string intree = "";
    string model = "MF"; // modelfinder
    int i;

    try {
        in_options = new input_options();
        // handle model_set, freq_set, rate_set
        if (!model_set.empty())
            in_options->insert("-mset", model_set);
        if (!freq_set.empty())
            in_options->insert("-mfreq", freq_set);
        if (!rate_set.empty())
            in_options->insert("-mrate", rate_set);
        if (num_thres > 1)
            in_options->insert("-nt", convertIntToString(num_thres));

        output = build_phylogenetic(names, seqs, model, intree, rand_seed, "modelfinder", in_options);
        
        delete in_options;
    } catch (std::runtime_error& e) {
        // reset the output and error buffers
        funcExit();
        throw e;
    }
    return output;
}

// Build pairwise JC distance matrix
// output: set of distances
// (n * i + j)-th element of the list represents the distance between i-th and j-th sequence,
// where n is the number of sequences
// num_thres -- number of cpu threads to be used, default: 1
vector<double> build_distmatrix(vector<string>& names, vector<string>& seqs, int num_thres) {
    vector<double> output;
    string prog = "build_matrix";
    
    try {
        int n = names.size();
        int n_sq = n * n;
        output.clear();
        if (n == 1) {
            output.push_back(0.0);
        } else {
            output.resize(n_sq);
            extern VerboseMode verbose_mode;
            progress_display::setProgressDisplay(false);
            verbose_mode = VB_QUIET; // (quiet mode)
            Params params = Params::getInstance();
            params.setDefault();

            int rand_seed = make_new_seed();
            string out_prefix_str = prog + "_" + convertIntToString(rand_seed);
            _log_file = out_prefix_str + ".log";
            bool append_log = false;
            startLogFile(append_log);

        #ifdef _OPENMP
            int max_procs = countPhysicalCPUCores();
            if (num_thres > max_procs)
                num_thres = max_procs;
            if (num_thres > 0) {
                Params::getInstance().num_threads = num_thres;
                omp_set_num_threads(num_thres);
            }
         #endif
            
            PhyloTree ptree;
            ptree.aln = new Alignment(names, seqs, params.sequence_type, params.model_name);
            
            // compute the matrix
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int i = 0; i < n; i++) {
                double* dmat = &output[i * n];
                // j = i
                dmat[i] = 0.0;
                // j != i
                for (int j = i+1; j < n; j++) {
                    dmat[j] = ptree.aln->computeJCDist(i, j);
                    output[j * n + i] = dmat[j];
                }
            }
            
            delete ptree.aln;
            funcExit();
        }
    } catch (std::runtime_error& e) {
        // reset the output and error buffers
        funcExit();
        throw e;
    }
    return output;
}

// Using Rapid-NJ to build tree from a distance matrix
string build_njtree(vector<string>& names, vector<double>& distances) {
    string output;

    try {
        // check the size of names and distances
        if (names.size() < 3)
            outError("The size of names must be at least 3");
        size_t n = names.size();
        size_t sq_n = n * n;
        if (distances.size() != sq_n)
            outError("The size of distances must equal to the square of the size of names");
        
        string prog = "build_njtree";
        extern VerboseMode verbose_mode;
        progress_display::setProgressDisplay(false);
        verbose_mode = VB_QUIET; // (quiet mode)
        Params params = Params::getInstance();
        params.setDefault();
        
        int rand_seed = make_new_seed();
        string out_prefix_str = prog + "_" + convertIntToString(rand_seed);
        _log_file = out_prefix_str + ".log";
        bool append_log = false;
        startLogFile(append_log);
        
        string algn_name = "NJ-R"; // Rapid NJ
        StartTree::BuilderInterface* algorithm = StartTree::Factory::getTreeBuilderByName(algn_name);
        stringstream stree;
        double* distMatrix = &distances[0];
        if (!algorithm->constructTreeInMemory2(names, distMatrix, stree)) {
            outError("Tree construction failed.");
        }
        output = stree.str();
        funcExit();
    } catch (std::runtime_error& e) {
        // reset the output and error buffers
        funcExit();
        throw e;
    }
    return output;
}

// ----------------------------------------------
// function for performing plylogenetic analysis
// ----------------------------------------------

// Perform phylogenetic analysis on the input alignment (in string format)
// if intree exists, then the topology will be restricted to the intree
string build_phylogenetic(vector<string>& names, vector<string>& seqs, string model, string intree,
                          int rand_seed, string prog, input_options* in_options) {
    // perform phylogenetic analysis on the input sequences
    // all sequences have to be the same length

    int instruction_set;
    
    // checking whether all seqs are in the same length
    if (seqs.size() > 0) {
        int slen = seqs[0].length();
        for (int i=1; i<seqs.size(); i++) {
            if (seqs[i].length() != slen) {
                outError("The input sequences are not in the same length");
            }
        }
    }

    extern VerboseMode verbose_mode;
    progress_display::setProgressDisplay(false);
    // verbose_mode = VB_MIN;
    verbose_mode = VB_QUIET; // (quiet mode)
    Params::getInstance().setDefault();
    Params::getInstance().num_threads = 1; // only allow single thread at this moment
    Params::getInstance().aln_file = (char*) "";
    Params::getInstance().model_name = model;
    
    if (intree != "") {
        // tree exists, then the resulting phylogenetic tree will be restricted to the input topology
        Params::getInstance().min_iterations = 0;
        Params::getInstance().stop_condition = SC_FIXED_ITERATION;
        Params::getInstance().start_tree = STT_USER_TREE;
        Params::getInstance().intree_str = intree;
    }

    if (in_options != NULL) {
        // assign the input options to Params
        in_options->set_params(Params::getInstance());
    }

    if (rand_seed == 0)
        rand_seed = make_new_seed();
    Params::getInstance().ran_seed = rand_seed;
    // cout << "Seed: " << Params::getInstance().ran_seed << endl << flush;
    init_random(Params::getInstance().ran_seed);

    string out_prefix_str = prog + "_" + convertIntToString(rand_seed);
    Params::getInstance().out_prefix = (char *) out_prefix_str.c_str();

    Checkpoint *checkpoint = new Checkpoint;
    string filename = (string)Params::getInstance().out_prefix +".ckp.gz";
    checkpoint->setFileName(filename);

    bool append_log = false;

    if (!Params::getInstance().ignore_checkpoint && fileExists(filename)) {
        checkpoint->load();
        if (checkpoint->hasKey("finished")) {
            if (checkpoint->getBool("finished")) {
                if (Params::getInstance().force_unfinished) {
                    if (MPIHelper::getInstance().isMaster())
                        cout << "NOTE: Continue analysis although a previous run already finished" << endl;
                } else {
                    delete checkpoint;
                    if (MPIHelper::getInstance().isMaster())
                        outError("Checkpoint (" + filename + ") indicates that a previous run successfully finished\n" +
                            "Use `-redo` option if you really want to redo the analysis and overwrite all output files.\n" +
                            "Use `--redo-tree` option if you want to restore ModelFinder and only redo tree search.\n" +
                            "Use `--undo` option if you want to continue previous run when changing/adding options."
                        );
                    else
                        exit(EXIT_SUCCESS);
                    exit(EXIT_FAILURE);
                }
            } else {
                append_log = true;
            }
        } else {
            if (MPIHelper::getInstance().isMaster())
                outWarning("Ignore invalid checkpoint file " + filename);
            checkpoint->clear();
        }
    }

    if (MPIHelper::getInstance().isWorker())
        checkpoint->setFileName("");

    _log_file = Params::getInstance().out_prefix;
    _log_file += ".log";
    startLogFile(append_log);
    time_t start_time;

    if (append_log) {
        cout << endl << "******************************************************"
             << endl << "CHECKPOINT: Resuming analysis from " << filename << endl << endl;
    }

    MPIHelper::getInstance().syncRandomSeed();

    signal(SIGABRT, &funcAbort);
    signal(SIGFPE, &funcAbort);
    signal(SIGILL, &funcAbort);
    signal(SIGSEGV, &funcAbort);
#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__ && !defined WIN64
    signal(SIGBUS, &funcAbort);
#endif
    printCopyright(cout);

    char hostname[100];
#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
    gethostname(hostname, sizeof(hostname));
    WSACleanup();
#else
    gethostname(hostname, sizeof(hostname));
#endif

    instruction_set = instrset_detect();
#if defined(BINARY32) || defined(__NOAVX__)
    instruction_set = min(instruction_set, (int)LK_SSE42);
#endif
    if (instruction_set < LK_SSE2) outError("Your CPU does not support SSE2!");
    bool has_fma3 = (instruction_set >= LK_AVX) && hasFMA3();

#ifdef __FMA__
    bool has_fma =  has_fma3;
    if (!has_fma) {
        outError("Your CPU does not support FMA instruction, quiting now...");
    }
#endif

    cout << "Host:    " << hostname << " (";
    switch (instruction_set) {
    case 0: cout << "x86, "; break;
    case 1: cout << "SSE, "; break;
    case 2: cout << "SSE2, "; break;
    case 3: cout << "SSE3, "; break;
    case 4: cout << "SSSE3, "; break;
    case 5: cout << "SSE4.1, "; break;
    case 6: cout << "SSE4.2, "; break;
    case 7: cout << "AVX, "; break;
    case 8: cout << "AVX2, "; break;
    default: cout << "AVX512, "; break;
    }
    if (has_fma3) cout << "FMA3, ";
    cout << (int)(((getMemorySize()/1024.0)/1024)/1024) << " GB RAM)" << endl;

    time(&start_time);
    cout << "Time:    " << ctime(&start_time);

    // increase instruction set level with FMA
    if (has_fma3 && instruction_set < LK_AVX_FMA)
        instruction_set = LK_AVX_FMA;

    Params::getInstance().SSE = min(Params::getInstance().SSE, (LikelihoodKernel)instruction_set);

    cout << "Kernel:  ";

    if (Params::getInstance().lk_safe_scaling) {
        cout << "Safe ";
    }

    if (Params::getInstance().pll) {
#ifdef __AVX__
        cout << "PLL-AVX";
#else
        cout << "PLL-SSE3";
#endif
    } else {
        if (Params::getInstance().SSE >= LK_AVX512)
            cout << "AVX-512";
        else if (Params::getInstance().SSE >= LK_AVX_FMA) {
            cout << "AVX+FMA";
        } else if (Params::getInstance().SSE >= LK_AVX) {
            cout << "AVX";
        } else if (Params::getInstance().SSE >= LK_SSE2){
            cout << "SSE2";
        } else
            cout << "x86";
    }

#ifdef _OPENMP
    if (Params::getInstance().num_threads >= 1) {
        omp_set_num_threads(Params::getInstance().num_threads);
        Params::getInstance().num_threads = omp_get_max_threads();
    }
//    int max_threads = omp_get_max_threads();
    int max_procs = countPhysicalCPUCores();
    cout << " - ";
    if (Params::getInstance().num_threads > 0)
        cout << Params::getInstance().num_threads  << " threads";
    else
        cout << "auto-detect threads";
    cout << " (" << max_procs << " CPU cores detected)";
    if (Params::getInstance().num_threads  > max_procs) {
        cout << endl;
        outError("You have specified more threads than CPU cores available");
    }
    // omp_set_nested(false); // don't allow nested OpenMP parallelism
    omp_set_max_active_levels(1);
#else
    if (Params::getInstance().num_threads != 1) {
        cout << endl << endl;
        outError("Number of threads must be 1 for sequential version.");
    }
#endif

    int num_procs = countPhysicalCPUCores();

    //cout << "sizeof(int)=" << sizeof(int) << endl;
    cout << endl << endl;
    
    // show msgs which are delayed to show
    cout << Params::getInstance().delay_msgs;

    cout.precision(3);
    cout.setf(ios::fixed);
    
    // checkpoint general run information
    checkpoint->startStruct("iqtree");
    int seed = Params::getInstance().ran_seed;
    CKP_SAVE(seed);
    CKP_SAVE(start_time);

    // check for incompatible version
    string version;
    stringstream sversion;
    sversion << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR << iqtree_VERSION_PATCH;
    version = sversion.str();
    CKP_SAVE(version);
    checkpoint->endStruct();
    
    Params params = Params::getInstance();
    IQTree *tree;
    Alignment *alignment = new Alignment(names, seqs, params.sequence_type, params.model_name);
    bool align_is_given = true;
    ModelCheckpoint* model_info = NULL;
    if (model == "MF") {
        model_info = new ModelCheckpoint;
    }

    runPhyloAnalysis(params, checkpoint, tree, alignment, align_is_given, model_info);
    
    stringstream ss;
    if (model_info != NULL) {
        // output the modelfinder results in YAML format
        model_info->dump(ss);
    } else {
        // output the checkpoint in YAML format
        checkpoint->dump(ss);
    }
    
    alignment = tree->aln;
    delete tree;
    delete alignment;
    cleanup(params);

    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    try{
        delete checkpoint;
        if (model_info != NULL)
            delete model_info;
    }catch(int err_num){}

    finish_random();
    funcExit();
    
    return ss.str();
}

// --------------------------------------------------
// Handle the input options of PiQTREE
// --------------------------------------------------

void input_options::set_params(Params& params) {
    ASSERT(flags.size() == values.size());
    int n = flags.size();
    for (int i = 0; i < n; i++) {
        if (flags[i] == "-keep-indent") {
            params.ignore_identical_seqs = false;
            cout << "params.ignore_identical_seqs = " << params.ignore_identical_seqs << endl;
        }
        else if (flags[i] == "-mset") {
            params.model_set = values[i];
            cout << "params.model_set = " << params.model_set << endl;
        }
        else if (flags[i] == "-mfreq") {
            int clen = values[i].length();
            if (clen > 0) {
                params.state_freq_set = new char[clen + 1];
                strcpy(params.state_freq_set, values[i].c_str());
            }
        }
        else if (flags[i] == "-mrate") {
            params.ratehet_set = values[i];
            cout << "params.ratehet_set = " << params.ratehet_set << endl;
        }
        else if (flags[i] == "-bb") {
            params.gbo_replicates = atoi(values[i].c_str());
            if (params.gbo_replicates < 1000)
                outError("#replicates must be >= 1000");
            params.consensus_type = CT_CONSENSUS_TREE;
            params.stop_condition = SC_BOOTSTRAP_CORRELATION;
        }
        else if (flags[0] == "-nt") {
            params.num_threads = atoi(values[i].c_str());
        }
    }
}
