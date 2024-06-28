#include "libiqtree2_fun.h"

string build_phylogenetic(vector<string> names, vector<string> seqs, string model, string intree, int rand_seed);

// Calculates the robinson fould distance between two trees
int robinson_fould(const string& tree1, const string& tree2) {
    // Placeholder implementation
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

    return (int)rfdist[0];
}

// Generates a set of random phylogenetic trees
// tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
string random_tree(int num_taxa, string tree_gen_mode, int num_trees, int rand_seed) {
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
        cerr << "Unknown mode: " << tree_gen_mode << endl;
        exit(1);
    }
    
    Params params = Params::getInstance();
    params.setDefault();
    params.sub_size = num_taxa;
    params.tree_gen = tree_mode;
    params.repeated_time = num_trees;
    params.ignore_checkpoint = true; // overrid the output file if exists
    params.user_file = "";

    ostringstream ostring;
    generateRandomTree(params, ostring);
    return ostring.str();
    // generateRandomTree(params);
    // return "done";
}

// Perform phylogenetic analysis on the input alignment (in string format)
// With estimation of the best topology
string build_tree(vector<string> names, vector<string> seqs, string model, int rand_seed) {
    string intree = "";
    return build_phylogenetic(names, seqs, model, intree, rand_seed);
}

// Perform phylogenetic analysis on the input alignment (in string format)
// With restriction to the input toplogy
string fit_tree(vector<string> names, vector<string> seqs, string model, string intree, int rand_seed) {
    return build_phylogenetic(names, seqs, model, intree, rand_seed);
}



// ----------------------------------------------
// function for performing plylogenetic analysis
// ----------------------------------------------

// Perform phylogenetic analysis on the input alignment (in string format)
// if intree exists, then the topology will be restricted to the intree
string build_phylogenetic(vector<string> names, vector<string> seqs, string model, string intree, int rand_seed) {
    // perform phylogenetic analysis on the input sequences
    // all sequences have to be the same length
    
    // checking whether all seqs are in the same length
    if (seqs.size() > 0) {
        int slen = seqs[0].length();
        for (int i=1; i<seqs.size(); i++) {
            if (seqs[i].length() != slen) {
                outError("The input sequences are not in the same length");
                exit(EXIT_FAILURE);
            }
        }
    }

    extern VerboseMode verbose_mode;
    progress_display::setProgressDisplay(false);
    verbose_mode = VB_MIN; // or VB_QUIET (quiet mode)
    Params::getInstance().setDefault();
    Params::getInstance().num_threads = 1; // only allow single thread at this moment
    Params::getInstance().aln_file = "";
    Params::getInstance().model_name = model;
    
    if (intree != "") {
        // tree exists, then the resulting phylogenetic tree will be restricted to the input topology
        Params::getInstance().min_iterations = 0;
        Params::getInstance().stop_condition = SC_FIXED_ITERATION;
        Params::getInstance().start_tree = STT_USER_TREE;
        Params::getInstance().intree_str = intree;
    }
    
    if (rand_seed == 0)
        rand_seed = make_new_seed();
    Params::getInstance().ran_seed = rand_seed;
    cout << "Seed:    " << Params::getInstance().ran_seed <<  " ";
    init_random(Params::getInstance().ran_seed);

    string out_prefix_str = "build_tree_" + convertIntToString(rand_seed);
    Params::getInstance().out_prefix = (char *) out_prefix_str.c_str();

    Checkpoint *checkpoint = new Checkpoint;
    string filename = (string)Params::getInstance().out_prefix +".ckp.gz";
    checkpoint->setFileName(filename);

    bool append_log = false;
    int instruction_set;

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
    omp_set_nested(false); // don't allow nested OpenMP parallelism
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
    
    runPhyloAnalysis(params, checkpoint, tree, alignment, align_is_given);
    // output the checkpoint in YAML format
    stringstream ss;
    checkpoint->dump(ss);
    
    alignment = tree->aln;
    delete tree;
    delete alignment;

    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    try{
    delete checkpoint;
    }catch(int err_num){}

    finish_random();
    
    return ss.str();
}
