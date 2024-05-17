#include "libiqtree2_fun.h"

// Calculates the RF distance between two trees
int calculate_RF_distance(const string& tree1, const string& tree2) {
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

// Generates a random phylogenetic tree
void generate_random_tree_file(int numtaxa, int seed, string tree_gen_mode, string outfile) {
    PhyloTree ptree;
    
    
    init_random(seed);
    
    TreeGenType genMode;
    if (tree_gen_mode == "YULE_HARDING") {
        genMode = YULE_HARDING;
    } else if (tree_gen_mode == "UNIFORM") {
        genMode = UNIFORM;
    } else if (tree_gen_mode == "CATERPILLAR") {
        genMode = CATERPILLAR;
    } else if (tree_gen_mode == "BALANCED") {
        genMode = BALANCED;
    } else if (tree_gen_mode == "BIRTH_DEATH") {
        genMode = BIRTH_DEATH;
    } else if (tree_gen_mode == "STAR_TREE") {
        genMode = STAR_TREE;
    } else if (tree_gen_mode == "CIRCULAR_SPLIT_GRAPH") {
        genMode = CIRCULAR_SPLIT_GRAPH;
    } else if (tree_gen_mode == "TAXA_SET") {
        genMode = TAXA_SET;
    } else {
        cerr << "Unknown mode: " << tree_gen_mode << endl;
        return;
    }

    Params params = Params::getInstance();
    params.setDefault();
    params.sub_size = numtaxa;
    params.tree_gen = genMode;
    params.user_file = (char *) outfile.c_str();
    params.ignore_checkpoint = true; // overrid the output file if exists

    // default parameter values
    /*
    params.repeated_time = 1; // one random tree to generate
    params.max_len = 1.0; // maximum branch length
    params.min_len = 1e-6; // minimum branch length
    params.mean_len = 0.1; // mean branch length\
    */

    generateRandomTree(params);
}

void phylogenetic_analysis(string& align_file, int ncpus) {
    // perform phylogenetic analysis on the input alignment file

    extern VerboseMode verbose_mode;
    progress_display::setProgressDisplay(false);
    verbose_mode = VB_MIN;
    Params::getInstance().setDefault();
    Params::getInstance().num_threads = ncpus;
    Params::getInstance().aln_file = (char *) align_file.c_str();
    Params::getInstance().out_prefix = (char *) Params::getInstance().aln_file;
    if (Params::getInstance().out_prefix[strlen(Params::getInstance().out_prefix)-1] == '/' || Params::getInstance().out_prefix[strlen(Params::getInstance().out_prefix)-1] == '\\') {
        Params::getInstance().out_prefix[strlen(Params::getInstance().out_prefix)-1] = 0;
    }

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

    /*
    double x=1e-100;
    double y=1e-101;
    if (x > y) cout << "ok!" << endl;
    else cout << "shit!" << endl;
    */
    //FILE *pfile = popen("hostname","r");
    char hostname[100];
#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
    gethostname(hostname, sizeof(hostname));
    WSACleanup();
#else
    gethostname(hostname, sizeof(hostname));
#endif
    //fgets(hostname, sizeof(hostname), pfile);
    //pclose(pfile);

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
//    if (has_fma4) cout << "FMA4, ";
//#if defined __APPLE__ || defined __MACH__
    cout << (int)(((getMemorySize()/1024.0)/1024)/1024) << " GB RAM)" << endl;
//#else
//    cout << (int)(((getMemorySize()/1000.0)/1000)/1000) << " GB RAM)" << endl;
//#endif

    checkpoint->get("iqtree.seed", Params::getInstance().ran_seed);
    cout << "Seed:    " << Params::getInstance().ran_seed <<  " ";
    init_random(Params::getInstance().ran_seed + MPIHelper::getInstance().getProcessID(), true);

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
    
    runPhyloAnalysis(Params::getInstance(), checkpoint);
    
    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    try{
    delete checkpoint;
    }catch(int err_num){}

    finish_random();
}
