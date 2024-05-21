#include "suppFunc.h"

inline void separator(ostream &out, int type) {
    switch (type) {
    case 0:
        out << endl << "==============================================================================" << endl;
        break;
    case 1:
        out << endl << "-----------------------------------------------------------" << endl;
        break;
    default:
        break;
    }
}


void printCopyright(ostream &out) {
#ifdef IQ_TREE
     out << "IQ-TREE";
    #ifdef _IQTREE_MPI
    out << " MPI";
    #endif
    #ifdef _OPENMP
    out << " multicore";
    #endif
    #ifdef __AVX512KNL
    out << " Xeon Phi KNL";
    #endif
     out << " version ";
#else
     out << "PDA - Phylogenetic Diversity Analyzer version ";
#endif
    out << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR << iqtree_VERSION_PATCH << " COVID-edition";
    out << " for " << getOSName();
    out << " built " << __DATE__;
#if defined DEBUG
    out << " - debug mode";
#endif

#ifdef IQ_TREE
    out << endl << "Developed by Bui Quang Minh, James Barbetti, Nguyen Lam Tung,"
        << endl << "Olga Chernomor, Heiko Schmidt, Dominik Schrempf, Michael Woodhams, Ly Trong Nhan." << endl << endl;
#else
    out << endl << "Copyright (c) 2006-2014 Olga Chernomor, Arndt von Haeseler and Bui Quang Minh." << endl << endl;
#endif
}

void printRunMode(ostream &out, RunMode run_mode) {
    switch (run_mode) {
        case DETECTED: out << "Detected"; break;
        case GREEDY: out << "Greedy"; break;
        case PRUNING: out << "Pruning"; break;
        case BOTH_ALG: out << "Greedy and Pruning"; break;
        case EXHAUSTIVE: out << "Exhaustive"; break;
        case DYNAMIC_PROGRAMMING: out << "Dynamic Programming"; break;
        case LINEAR_PROGRAMMING: out << "Integer Linear Programming"; break;
        default: outError(ERR_INTERNAL);
    }
}

/**
    summarize the running with header
*/
void summarizeHeader(ostream &out, Params &params, bool budget_constraint, InputType analysis_type) {
    printCopyright(out);
    out << "Input tree/split network file name: " << params.user_file << endl;
    if(params.eco_dag_file)
        out << "Input food web file name: "<<params.eco_dag_file<<endl;
     out << "Input file format: " << ((params.intype == IN_NEWICK) ? "Newick" : ( (params.intype == IN_NEXUS) ? "Nexus" : "Unknown" )) << endl;
    if (params.initial_file != NULL)
        out << "Initial taxa file: " << params.initial_file << endl;
    if (params.param_file != NULL)
        out << "Parameter file: " << params.param_file << endl;
    out << endl;
    out << "Type of measure: " << ((params.root != NULL || params.is_rooted) ? "Rooted": "Unrooted") <<
            (analysis_type== IN_NEWICK ? " phylogenetic diversity (PD)" : " split diversity (SD)");
    if (params.root != NULL) out << " at " << params.root;
    out << endl;
    if (params.run_mode != CALC_DIST && params.run_mode != PD_USER_SET) {
        out << "Search objective: " << ((params.find_pd_min) ? "Minimum" : "Maximum") << endl;
        out << "Search algorithm: ";
        printRunMode(out, params.run_mode);
        if (params.run_mode == DETECTED) {
            out << " -> ";
            printRunMode(out, params.detected_mode);
        }
        out << endl;
        out << "Search option: " << ((params.find_all) ? "Multiple optimal sets" : "Single optimal set") << endl;
    }
    out << endl;
    out << "Type of analysis: ";
    switch (params.run_mode) {
        case PD_USER_SET: out << "PD/SD of user sets";
            if (params.pdtaxa_file) out << " (" << params.pdtaxa_file << ")"; break;
        case CALC_DIST: out << "Distance matrix computation"; break;
        default:
            out << ((budget_constraint) ? "Budget constraint " : "Subset size k ");
            if (params.intype == IN_NEWICK)
                out << ((analysis_type == IN_NEWICK) ? "on tree" : "on tree -> split network");
            else
                out << "on split network";
    }
    out << endl;
    //out << "Random number seed: " << params.ran_seed << endl;
}

void summarizeFooter(ostream &out, Params &params) {
    separator(out);
    time_t beginTime;
    time (&beginTime);
    char *date;
    date = ctime(&beginTime);

    out << "Time used: " << params.run_time  << " seconds." << endl;
    out << "Finished time: " << date << endl;
}


int getMaxNameLen(vector<string> &setName) {
    int len = 0;
    for (vector<string>::iterator it = setName.begin(); it != setName.end(); it++)
        if (len < (*it).length())
            len = (*it).length();
    return len;
}

void printPDUser(ostream &out, Params &params, PDRelatedMeasures &pd_more) {
    out << "List of user-defined sets of taxa with PD score computed" << endl << endl;
    int maxlen = getMaxNameLen(pd_more.setName)+2;
    out.width(maxlen);
    out << "Name" << "     PD";
    if (params.exclusive_pd) out << "   excl.-PD";
    if (params.endemic_pd) out << "   PD-Endem.";
    if (params.complement_area) out << "   PD-Compl. given area " << params.complement_area;
    out << endl;
    int cnt;
    for (cnt = 0; cnt < pd_more.setName.size(); cnt++) {
        out.width(maxlen);
        out << pd_more.setName[cnt] << " ";
        out.width(7);
        out << pd_more.PDScore[cnt] << "  ";
        if (params.exclusive_pd) {
            out.width(7);
            out << pd_more.exclusivePD[cnt] << "  ";
        }
        if (params.endemic_pd) {
            out.width(7);
            out << pd_more.PDEndemism[cnt] << "  ";
        }
        if (params.complement_area) {
            out.width(8);
            out << pd_more.PDComplementarity[cnt];
        }
        out << endl;
    }
    separator(out, 1);
}

void summarizeTree(Params &params, PDTree &tree, vector<PDTaxaSet> &taxa_set,
    PDRelatedMeasures &pd_more) {
    string filename;
    if (params.out_file == NULL) {
        filename = params.out_prefix;
        filename += ".pda";
    } else
        filename = params.out_file;

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename.c_str());

        summarizeHeader(out, params, false, IN_NEWICK);
        out << "Tree size: " << tree.leafNum-params.is_rooted << " taxa, " <<
            tree.nodeNum-1-params.is_rooted << " branches" << endl;
        separator(out);

        vector<PDTaxaSet>::iterator tid;

        if (params.run_mode == PD_USER_SET) {
            printPDUser(out, params, pd_more);
        }
        else if (taxa_set.size() > 1)
            out << "Optimal PD-sets with k = " << params.min_size-params.is_rooted <<
            " to " << params.sub_size-params.is_rooted << endl << endl;


        int subsize = params.min_size-params.is_rooted;
        if (params.run_mode == PD_USER_SET) subsize = 1;
        for (tid = taxa_set.begin(); tid != taxa_set.end(); tid++, subsize++) {
            if (tid != taxa_set.begin())
                separator(out, 1);
            if (params.run_mode == PD_USER_SET) {
                out << "Set " << subsize << " has PD score of " << tid->score << endl;
            }
            else {
                out << "For k = " << subsize << " the optimal PD score is " << (*tid).score << endl;
                out << "The optimal PD set has " << subsize << " taxa:" << endl;
            }
            for (NodeVector::iterator it = (*tid).begin(); it != (*tid).end(); it++)
                if ((*it)->name != ROOT_NAME){
                    out << (*it)->name << endl;
                }
            if (!tid->tree_str.empty()) {
                out << endl << "Corresponding sub-tree: " << endl;
                out << tid->tree_str << endl;
            }
            tid->clear();
        }
        taxa_set.clear();

        summarizeFooter(out, params);
        out.close();
        cout << endl << "Results are summarized in " << filename << endl << endl;

    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}


void printTaxaSet(Params &params, vector<PDTaxaSet> &taxa_set, RunMode cur_mode) {
    int subsize = params.min_size-params.is_rooted;
    ofstream out;
    ofstream scoreout;
    string filename;
    filename = params.out_prefix;
    filename += ".score";
    scoreout.open(filename.c_str());
    if (!scoreout.is_open())
        outError(ERR_WRITE_OUTPUT, filename);
    cout << "PD scores printed to " << filename << endl;

    if (params.nr_output == 1) {
        filename = params.out_prefix;
        filename += ".pdtaxa";
        out.open(filename.c_str());
        if (!out.is_open())
            outError(ERR_WRITE_OUTPUT, filename);
    }
    for (vector<PDTaxaSet>::iterator tid = taxa_set.begin(); tid != taxa_set.end(); tid++, subsize++) {
        if (params.nr_output > 10) {
            filename = params.out_prefix;
            filename += ".";
            filename += subsize;
            if (params.run_mode == BOTH_ALG) {
                if (cur_mode == GREEDY)
                    filename += ".greedy";
                else
                    filename += ".pruning";
            } else {
                filename += ".pdtree";
            }
            (*tid).printTree((char*)filename.c_str());

            filename = params.out_prefix;
            filename += ".";
            filename += subsize;
            filename += ".pdtaxa";
            (*tid).printTaxa((char*)filename.c_str());
        } else {
            out << subsize << " " << (*tid).score << endl;
            scoreout << subsize << " " << (*tid).score << endl;
            (*tid).printTaxa(out);
        }
    }

    if (params.nr_output == 1) {
        out.close();
        cout << "All taxa list(s) printed to " << filename << endl;
    }

    scoreout.close();
}


/**
    run PD algorithm on trees
*/
void runPDTree(Params &params)
{

    if (params.run_mode == CALC_DIST) {
        bool is_rooted = false;
        MExtTree tree(params.user_file, is_rooted);
        cout << "Tree contains " << tree.leafNum << " taxa." << endl;
        cout << "Calculating distance matrix..." << endl;
        tree.calcDist(params.dist_file);
        cout << "Distances printed to " << params.dist_file << endl;
        return;
    }

    double t_begin, t_end;
    //char filename[300];
    //int idx;

    vector<PDTaxaSet> taxa_set;

    if (params.run_mode == PD_USER_SET) {
        // compute score of user-defined sets
        t_begin = getCPUTime();
        cout << "Computing PD score for user-defined set of taxa..." << endl;
        PDTree tree(params);
        PDRelatedMeasures pd_more;
        tree.computePD(params, taxa_set, pd_more);

        if (params.endemic_pd)
            tree.calcPDEndemism(taxa_set, pd_more.PDEndemism);
        if (params.complement_area != NULL)
            tree.calcPDComplementarity(taxa_set, params.complement_area, pd_more.PDComplementarity);

        t_end = getCPUTime();
        params.run_time = (t_end-t_begin);
        summarizeTree(params, tree, taxa_set, pd_more);
        return;
    }


    /*********************************************
        run greedy algorithm
    *********************************************/

    if (params.sub_size < 2) {
        outError(ERR_NO_K);
    }

    bool detected_greedy = (params.run_mode != PRUNING);

    Greedy test_greedy;

    test_greedy.init(params);

    if (params.root == NULL && !params.is_rooted)
        cout << endl << "Running PD algorithm on UNROOTED tree..." << endl;
    else
        cout << endl << "Running PD algorithm on ROOTED tree..." << endl;

    if (verbose_mode >= VB_DEBUG)
        test_greedy.drawTree(cout, WT_INT_NODE + WT_BR_SCALE + WT_BR_LEN);

    if (params.run_mode == GREEDY || params.run_mode == BOTH_ALG ||
        (params.run_mode == DETECTED)) {

        if (params.run_mode == DETECTED && params.sub_size >= test_greedy.leafNum * 7 / 10
            && params.min_size < 2)
            detected_greedy = false;

        if (detected_greedy) {
            params.detected_mode = GREEDY;
            t_begin=getCPUTime();
            cout << endl << "Greedy Algorithm..." << endl;

            taxa_set.clear();
            test_greedy.run(params, taxa_set);

            t_end=getCPUTime();
            params.run_time = (t_end-t_begin);
            cout << "Time used: " << params.run_time << " seconds." << endl;
            if (params.min_size == params.sub_size)
                cout << "Resulting tree length = " << taxa_set[0].score << endl;

            if (params.nr_output > 0)
                printTaxaSet(params, taxa_set, GREEDY);

            PDRelatedMeasures pd_more;

            summarizeTree(params, test_greedy, taxa_set, pd_more);
        }
    }

    /*********************************************
        run pruning algorithm
    *********************************************/
    if (params.run_mode == PRUNING || params.run_mode == BOTH_ALG ||
        (params.run_mode == DETECTED)) {

        Pruning test_pruning;

        if (params.run_mode == PRUNING || params.run_mode == BOTH_ALG) {
            //Pruning test_pruning(params);
            test_pruning.init(params);
        } else if (!detected_greedy) {
            test_pruning.init(test_greedy);
        } else {
            return;
        }
        params.detected_mode = PRUNING;
        t_begin=getCPUTime();
        cout << endl << "Pruning Algorithm..." << endl;
        taxa_set.clear();
        test_pruning.run(params, taxa_set);

        t_end=getCPUTime();
        params.run_time = (t_end-t_begin) ;
        cout << "Time used: " << params.run_time << " seconds.\n";
        if (params.min_size == params.sub_size)
            cout << "Resulting tree length = " << taxa_set[0].score << endl;

        if (params.nr_output > 0)
            printTaxaSet(params, taxa_set, PRUNING);

        PDRelatedMeasures pd_more;

        summarizeTree(params, test_pruning, taxa_set, pd_more);

    }

}

void checkSplitDistance(ostream &out, PDNetwork &sg) {
    mmatrix(double) dist;
    sg.calcDistance(dist);
    int ntaxa = sg.getNTaxa();
    int i, j;
    bool found = false;
    for (i = 0; i < ntaxa-1; i++) {
        bool first = true;
        for (j = i+1; j < ntaxa; j++)
            if (abs(dist[i][j]) <= 1e-5) {
                if (!found) {
                    out << "The following sets of taxa (each set in a line) have very small split-distance" << endl;
                    out << "( <= 1e-5) as computed from the split system. To avoid a lot of multiple" << endl;
                    out << "optimal PD sets to be reported, one should only keep one taxon from each set" << endl;
                    out << "and exclude the rest from the analysis." << endl << endl;
                }
                if (first)
                    out << sg.getTaxa()->GetTaxonLabel(i);
                found = true;
                first = false;
                out << ", " << sg.getTaxa()->GetTaxonLabel(j);
            }
        if (!first) out << endl;
    }
    if (found)
        separator(out);
}



/**
    check if the set are nested and there are no multiple optimal sets.
    If yes, return the ranking as could be produced by a greedy algorithm
*/
bool makeRanking(vector<SplitSet> &pd_set, IntVector &indices, IntVector &ranking) {
    vector<SplitSet>::iterator it;
    IntVector::iterator inti;
    ranking.clear();
    bool nested = true;
    Split *cur_sp = NULL;
    int id = 1;
    for (it = pd_set.begin(); it != pd_set.end(); it++) {
        if ((*it).empty()) continue;
        if ((*it).size() > 1) {
            nested = false;
            ranking.push_back(-10);
            indices.push_back(0);
        }
        Split *sp = (*it)[0];

        if (!cur_sp) {
            IntVector sp_tax;
            sp->getTaxaList(sp_tax);
            ranking.insert(ranking.end(), sp_tax.begin(), sp_tax.end());
            for (inti = sp_tax.begin(); inti != sp_tax.end(); inti++)
                indices.push_back(id++);
        } else {
            if ( !cur_sp->subsetOf(*sp)) {
                ranking.push_back(-1);
                indices.push_back(0);
                nested = false;
            }
            Split sp_diff(*sp);
            sp_diff -= *cur_sp;
            Split sp_diff2(*cur_sp);
            sp_diff2 -= *sp;
            IntVector sp_tax;
            sp_diff2.getTaxaList(sp_tax);
            ranking.insert(ranking.end(), sp_tax.begin(), sp_tax.end());
            for (inti = sp_tax.begin(); inti != sp_tax.end(); inti++)
                indices.push_back(-id);
            sp_diff.getTaxaList(sp_tax);
            ranking.insert(ranking.end(), sp_tax.begin(), sp_tax.end());
            for (inti = sp_tax.begin(); inti != sp_tax.end(); inti++)
                indices.push_back(id);
            if ( !cur_sp->subsetOf(*sp)) {
                ranking.push_back(-2);
                indices.push_back(0);
            }
            id++;
        }
        cur_sp = sp;
    }
    return nested;
}


void printNexusSets(const char *filename, PDNetwork &sg, vector<SplitSet> &pd_set) {
    try {
        ofstream out;
        out.open(filename);
        out << "#NEXUS" << endl << "BEGIN Sets;" << endl;
        vector<SplitSet>::iterator it;
        for (it = pd_set.begin(); it != pd_set.end(); it++) {
            int id = 1;
            for (SplitSet::iterator sit = (*it).begin(); sit != (*it).end(); sit++, id++) {
                IntVector taxa;
                (*sit)->getTaxaList(taxa);
                out << "   TAXSET Opt_" << taxa.size() << "_" << id << " =";
                for (IntVector::iterator iit = taxa.begin(); iit != taxa.end(); iit++) {
                    if (sg.isPDArea())
                        out << " '" << sg.getSetsBlock()->getSet(*iit)->name << "'";
                    else
                        out << " '" << sg.getTaxa()->GetTaxonLabel(*iit) << "'";
                }
                out << ";" << endl;
            }
        }
        out << "END; [Sets]" << endl;
        out.close();
        cout << endl << "Optimal sets are written to nexus file " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }

}



void computeTaxaFrequency(SplitSet &taxa_set, DoubleVector &freq) {
    ASSERT(taxa_set.size());
    int ntaxa = taxa_set[0]->getNTaxa();
    int i;

    freq.resize(ntaxa, 0);
    for (SplitSet::iterator it2 = taxa_set.begin(); it2 != taxa_set.end(); it2++) {
        for ( i = 0; i < ntaxa; i++)
            if ((*it2)->containTaxon(i)) freq[i] += 1.0;
    }

    for ( i = 0; i < ntaxa; i++)
        freq[i] /= taxa_set.size();

}

/**
    summarize the running results
*/
void summarizeSplit(Params &params, PDNetwork &sg, vector<SplitSet> &pd_set, PDRelatedMeasures &pd_more, bool full_report) {
    int i;


    if (params.nexus_output) {
        string nex_file = params.out_prefix;
        nex_file += ".pdsets.nex";
        printNexusSets(nex_file.c_str(), sg, pd_set);
    }
    string filename;
    if (params.out_file == NULL) {
        filename = params.out_prefix;
        filename += ".pda";
    } else
        filename = params.out_file;

    try {
        ofstream out;
        out.open(filename.c_str());
        /****************************/
        /********** HEADER **********/
        /****************************/
        summarizeHeader(out, params, sg.isBudgetConstraint(), IN_NEXUS);

        out << "Network size: " << sg.getNTaxa()-params.is_rooted << " taxa, " <<
            sg.getNSplits()-params.is_rooted << " splits (of which " <<
            sg.getNTrivialSplits() << " are trivial splits)" << endl;
        out << "Network type: " << ((sg.isCircular()) ? "Circular" : "General") << endl;

        separator(out);

        checkSplitDistance(out, sg);

        int c_num = 0;
        //int subsize = (sg.isBudgetConstraint()) ? params.budget : (params.sub_size-params.is_rooted);
        //subsize -= pd_set.size()-1;
        int subsize = (sg.isBudgetConstraint()) ? params.min_budget : params.min_size-params.is_rooted;
        int stepsize = (sg.isBudgetConstraint()) ? params.step_budget : params.step_size;
        if (params.detected_mode != LINEAR_PROGRAMMING) stepsize = 1;
        vector<SplitSet>::iterator it;
        SplitSet::iterator it2;


        if (params.run_mode == PD_USER_SET) {
            printPDUser(out, params, pd_more);
        }

        /****************************/
        /********** SUMMARY *********/
        /****************************/

        if (params.run_mode != PD_USER_SET && !params.num_bootstrap_samples) {
            out << "Summary of the PD-score and the number of optimal PD-sets with the same " << endl << "optimal PD-score found." << endl;

            if (sg.isBudgetConstraint())
                out << endl << "Budget   PD-score   %PD-score   #PD-sets" << endl;
            else
                out << endl << "Size-k   PD-score   %PD-score   #PD-sets" << endl;

            int sizex = subsize;
            double total = sg.calcWeight();

            for (it = pd_set.begin(); it != pd_set.end(); it++, sizex+=stepsize) {
                out.width(6);
                out << right << sizex << " ";
                out.width(10);
                out << right << (*it).getWeight() << " ";
                out.width(10);
                out << right << ((*it).getWeight()/total)*100.0 << " ";
                out.width(6);
                out << right << (*it).size();
                out << endl;
            }

            out << endl;
            if (!params.find_all)
                out << "Note: You did not choose the option to find multiple optimal PD sets." << endl <<
                    "That's why we only reported one PD-set per size-k or budget. If you want" << endl <<
                    "to determine all multiple PD-sets, use the '-a' option.";
            else {
                out << "Note: The number of multiple optimal PD sets to be reported is limited to " << params.pd_limit << "." << endl <<
                    "There might be cases where the actual #PD-sets exceeds that upper-limit but" << endl <<
                    "won't be listed here. Please refer to the above list to identify such cases." << endl <<
                    "To increase the upper-limit, use the '-lim <limit_number>' option.";
            }
            out << endl;
            separator(out);
        }

        if (!full_report) {
            out.close();
            return;
        }


        /****************************/
        /********* BOOTSTRAP ********/
        /****************************/
        if (params.run_mode != PD_USER_SET && params.num_bootstrap_samples) {
            out << "Summary of the bootstrap analysis " << endl;
            for (it = pd_set.begin(); it != pd_set.end(); it++) {
                DoubleVector freq;
                computeTaxaFrequency((*it), freq);
                out << "For k/budget = " << subsize << " the " << ((sg.isPDArea()) ? "areas" : "taxa")
                    << " supports are: " << endl;
                for (i = 0; i < freq.size(); i++)
                    out << ((sg.isPDArea()) ? sg.getSetsBlock()->getSet(i)->name : sg.getTaxa()->GetTaxonLabel(i))
                        << "\t" << freq[i] << endl;
                if ((it+1) != pd_set.end()) separator(out, 1);
            }
            out << endl;
            separator(out);
        }

        /****************************/
        /********** RANKING *********/
        /****************************/

        if (params.run_mode != PD_USER_SET && !params.num_bootstrap_samples) {


            IntVector ranking;
            IntVector index;

            out << "Ranking based on the optimal sets" << endl;


            if (!makeRanking(pd_set, index, ranking)) {
                out << "WARNING: Optimal sets are not nested, so ranking should not be considered stable" << endl;
            }
            if (subsize > 1) {
                out << "WARNING: The first " << subsize << " ranks should be treated equal" << endl;
            }
            out << endl << "Rank*   ";
            if (!sg.isPDArea())
                out << "Taxon names" << endl;
            else
                out << "Area names" << endl;


            for (IntVector::iterator intv = ranking.begin(), intid = index.begin(); intv != ranking.end(); intv ++, intid++) {
                if (*intv == -10)
                    out << "<--- multiple optimal set here --->" << endl;
                else if (*intv == -1)
                    out << "<--- BEGIN: greedy does not work --->" << endl;
                else if (*intv == -2)
                    out << "<--- END --->" << endl;
                else {
                    out.width(5);
                    out <<  right << *intid << "   ";
                    if (sg.isPDArea())
                        out << sg.getSetsBlock()->getSet(*intv)->name << endl;
                    else
                        out << sg.getTaxa()->GetTaxonLabel(*intv) << endl;
                }
            }
            out << endl;
            out <<  "(*) Negative ranks indicate the point at which the greedy algorithm" << endl <<
                    "    does not work. In that case, the corresponding taxon/area names" << endl <<
                    "    should be deleted from the optimal set of the same size" << endl;
            separator(out);
        }

        int max_len = sg.getTaxa()->GetMaxTaxonLabelLength();

        /****************************/
        /***** DETAILED SETS ********/
        /****************************/

        if (params.run_mode != PD_USER_SET)
            out << "Detailed information of all taxa found in the optimal PD-sets" << endl;

        if (pd_set.size() > 1) {
            if (sg.isBudgetConstraint())
                out << "with budget = " << params.min_budget <<
                    " to " << params.budget << endl << endl;
            else
                out << "with k = " << params.min_size-params.is_rooted <<
                    " to " << params.sub_size-params.is_rooted << endl << endl;
        }

        if (params.run_mode != PD_USER_SET)
            separator(out,1);

        for (it = pd_set.begin(); it != pd_set.end(); it++, subsize+=stepsize) {

            // check if the pd-sets are the same as previous one
            if (sg.isBudgetConstraint() && it != pd_set.begin()) {
                vector<SplitSet>::iterator prev, next;
                for (next=it, prev=it-1; next != pd_set.end() && next->getWeight() == (*prev).getWeight() &&
                    next->size() == (*prev).size(); next++ ) ;
                if (next != it) {
                    // found something in between!
                    out << endl;
                    //out << endl << "**************************************************************" << endl;
                    out << "For budget = " << subsize << " -> " << subsize+(next-it-1)*stepsize <<
                        " the optimal PD score and PD sets" << endl;
                    out << "are identical to the case when budget = " << subsize-stepsize << endl;
                    //out << "**************************************************************" << endl;
                    subsize += (next-it)*stepsize;
                    it = next;
                    if (it == pd_set.end()) break;
                }
            }

            if (it != pd_set.begin()) separator(out, 1);

            int num_sets = (*it).size();
            double weight = (*it).getWeight();

            if (params.run_mode != PD_USER_SET) {
                out << "For " << ((sg.isBudgetConstraint()) ? "budget" : "k") << " = " << subsize;
                out << " the optimal PD score is " << weight << endl;

                if (num_sets == 1) {
                    if (!sg.isBudgetConstraint())
                        out << "The optimal PD set has " << (*it)[0]->countTaxa()-params.is_rooted <<
                            ((sg.isPDArea()) ? " areas" : " taxa");
                    else
                        out << "The optimal PD set has " << (*it)[0]->countTaxa()-params.is_rooted <<
                        ((sg.isPDArea()) ? " areas" : " taxa") << " and requires " << sg.calcCost(*(*it)[0]) << " budget";
                    if (!sg.isPDArea()) out << " and covers " << sg.countSplits(*(*it)[0]) <<
                        " splits (of which " << sg.countInternalSplits(*(*it)[0]) << " are internal splits)";
                    out << endl;
                }
                else
                    out << "Found " << num_sets << " PD sets with the same optimal score." << endl;
            }
            for (it2 = (*it).begin(), c_num=1; it2 != (*it).end(); it2++, c_num++){
                Split *this_set = *it2;

                if (params.run_mode == PD_USER_SET && it2 != (*it).begin())
                    separator(out, 1);

                if (params.run_mode == PD_USER_SET) {
                    if (!sg.isBudgetConstraint())
                        out << "Set " << c_num << " has PD score of " << this_set->getWeight();
                    else
                        out << "Set " << c_num << " has PD score of " << this_set->getWeight() <<
                        " and requires " << sg.calcCost(*this_set) << " budget";
                } else if (num_sets > 1) {
                    if (!sg.isBudgetConstraint())
                        out << endl << "PD set " << c_num;
                    else
                        out << endl << "PD set " << c_num << " has " << this_set->countTaxa()-params.is_rooted <<
                        " taxa and requires " << sg.calcCost(*this_set) << " budget";
                }

                if (!sg.isPDArea() && (num_sets > 1 || params.run_mode == PD_USER_SET ))
                    out << " and covers " << sg.countSplits(*(*it)[0]) << " splits (of which "
                    << sg.countInternalSplits(*(*it)[0]) << " are internal splits)";
                out << endl;

                if (params.run_mode != PD_USER_SET && sg.isPDArea()) {
                    for (i = 0; i < sg.getSetsBlock()->getNSets(); i++)
                        if (this_set->containTaxon(i)) {
                            if (sg.isBudgetConstraint()) {
                                out.width(max_len);
                                out << left << sg.getSetsBlock()->getSet(i)->name << "\t";
                                out.width(10);
                                out << right << sg.getPdaBlock()->getCost(i);
                                out << endl;

                            } else {
                                out << sg.getSetsBlock()->getSet(i)->name << endl;
                            }
                        }

                    Split sp(sg.getNTaxa());
                    for (i = 0; i < sg.getSetsBlock()->getNSets(); i++)
                        if (this_set->containTaxon(i))
                            sp += *(sg.area_taxa[i]);
                    out << endl << "which contains " << sp.countTaxa() - params.is_rooted << " taxa: " << endl;
                    for (i = 0; i < sg.getNTaxa(); i++)
                        if (sg.getTaxa()->GetTaxonLabel(i) != ROOT_NAME && sp.containTaxon(i))
                            out << sg.getTaxa()->GetTaxonLabel(i) << endl;

                } else
                for ( i = 0; i < sg.getNTaxa(); i++)
                    if (sg.getTaxa()->GetTaxonLabel(i) != ROOT_NAME && this_set->containTaxon(i)) {
                        if (sg.isBudgetConstraint()) {
                            out.width(max_len);
                            out << left << sg.getTaxa()->GetTaxonLabel(i) << "\t";
                            out.width(10);
                            out << right << sg.getPdaBlock()->getCost(i);
                            out << endl;

                        } else {
                            out << sg.getTaxa()->GetTaxonLabel(i) << endl;
                        }
                    }
            }
        }

        /****************************/
        /********** FOOTER **********/
        /****************************/

        summarizeFooter(out, params);

        out.close();
        cout << endl << "Results are summarized in " << filename << endl << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

void printGainMatrix(char *filename, mmatrix(double) &delta_gain, int start_k) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        int k = start_k;
        for (mmatrix(double)::iterator it = delta_gain.begin(); it != delta_gain.end(); it++, k++) {
            out << k;
            for (int i = 0; i < (*it).size(); i++)
                out << "  " << (*it)[i];
            out << endl;
        }
        out.close();
        cout << "PD gain matrix printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

/**
    run PD algorithm on split networks
*/
void runPDSplit(Params &params) {

    cout << "Using NCL - Nexus Class Library" << endl << endl;

    // init a split graph class from the parameters
    CircularNetwork sg(params);
    int i;

    // this vector of SplitSet store all the optimal PD sets
    vector<SplitSet> pd_set;
    // this define an order of taxa (circular order in case of circular networks)
    vector<int> taxa_order;
    // this store a particular taxa set
    Split taxa_set;


    if (sg.isCircular()) {
        // is a circular network, get circular order
        for (i = 0; i < sg.getNTaxa(); i++)
            taxa_order.push_back(sg.getCircleId(i));
    } else
        // otherwise, get the incremental order
        for (i = 0; i < sg.getNTaxa(); i++)
            taxa_order.push_back(i);

    PDRelatedMeasures pd_more;

    // begining time of the algorithm run
    double time_begin = getCPUTime();
    //time(&time_begin);
    // check parameters
    if (sg.isPDArea()) {
        if (sg.isBudgetConstraint()) {
            int budget = (params.budget >= 0) ? params.budget : sg.getPdaBlock()->getBudget();
            if (budget < 0 && params.pd_proportion == 0.0) params.run_mode = PD_USER_SET;
        } else {
            int sub_size = (params.sub_size >= 1) ? params.sub_size : sg.getPdaBlock()->getSubSize();
            if (sub_size < 1 && params.pd_proportion == 0.0) params.run_mode = PD_USER_SET;

        }
    }

    if (params.run_mode == PD_USER_SET) {
        // compute score of user-defined sets
        cout << "Computing PD score for user-defined set of taxa..." << endl;
        pd_set.resize(1);
        sg.computePD(params, pd_set[0], pd_more);
        if (params.endemic_pd)
            sg.calcPDEndemism(pd_set[0], pd_more.PDEndemism);

        if (params.complement_area != NULL)
            sg.calcPDComplementarity(pd_set[0], params.complement_area, pd_more.setName, pd_more.PDComplementarity);

    } else {
        // otherwise, call the main function
        if (params.num_bootstrap_samples) {
            cout << endl << "======= START BOOTSTRAP ANALYSIS =======" << endl;
            MTreeSet *mtrees = sg.getMTrees();
            if (mtrees->size() < 100)
                cout << "WARNING: bootstrap may be unstable with less than 100 trees" << endl;
            vector<string> taxname;
            sg.getTaxaName(taxname);
            i = 1;
            for (MTreeSet::iterator it = mtrees->begin(); it != mtrees->end(); it++, i++) {
                cout << "---------- TREE " << i << " ----------" << endl;
                // convert tree into split sytem
                SplitGraph sg2;
                (*it)->convertSplits(taxname, sg2);
                // change the current split system
                for (SplitGraph::reverse_iterator it = sg.rbegin(); it != sg.rend(); it++) {
                    delete *it;
                }
                sg.clear();
                sg.insert(sg.begin(), sg2.begin(), sg2.end());
                sg2.clear();

                // now findPD on the converted tree-split system
                sg.findPD(params, pd_set, taxa_order);
            }
            cout << "======= DONE BOOTSTRAP ANALYSIS =======" << endl << endl;
        } else {
            sg.findPD(params, pd_set, taxa_order);
        }
    }

    // ending time
    double time_end = getCPUTime();
    //time(&time_end);
    params.run_time = time_end - time_begin;

    cout << "Time used: " << (double) (params.run_time) << " seconds." << endl;

    if (verbose_mode >= VB_DEBUG && !sg.isPDArea()) {
        cout << "PD set(s) with score(s): " << endl;
        for (vector<SplitSet>::iterator it = pd_set.begin(); it != pd_set.end(); it++)
        for (SplitSet::iterator it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
            //(*it)->report(cout);
            cout << "  " << (*it2)->getWeight() << "    ";
            for (i = 0; i < sg.getNTaxa(); i++)
                if ((*it2)->containTaxon(i))
                cout << sg.getTaxa()->GetTaxonLabel(i) << "  ";
            if (sg.isBudgetConstraint())
                cout << " (budget = " << sg.calcCost(*(*it2)) << ")";
            cout << endl;
        }
    }

    sg.printOutputSetScore(params, pd_set);


    summarizeSplit(params, sg, pd_set, pd_more, true);

    if (params.calc_pdgain) {
        mmatrix(double) delta_gain;
        sg.calcPDGain(pd_set, delta_gain);
        string filename = params.out_prefix;
        filename += ".pdgain";
        printGainMatrix((char*)filename.c_str(), delta_gain, pd_set.front().front()->countTaxa());
        //cout << delta_gain;
    }


    //for (i = pd_set.size()-1; i >= 0; i--)
    //    delete pd_set[i];

}

void printSplitSet(SplitGraph &sg, SplitIntMap &hash_ss) {
/*
    for (SplitIntMap::iterator it = hash_ss.begin(); it != hash_ss.end(); it++) {
        if ((*it)->getWeight() > 50 && (*it)->countTaxa() > 1)
        (*it)->report(cout);
    }*/
    sg.getTaxa()->Report(cout);
    for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++) {
        if ((*it)->getWeight() > 50 && (*it)->countTaxa() > 1)
        (*it)->report(cout);
    }
}

void readTaxaOrder(char *taxa_order_file, StrVector &taxa_order) {

}

void calcTreeCluster(Params &params) {
    ASSERT(params.taxa_order_file);
    MExtTree tree(params.user_file, params.is_rooted);
//    StrVector taxa_order;
    //readTaxaOrder(params.taxa_order_file, taxa_order);
    NodeVector taxa;
    mmatrix(int) clusters;
    clusters.reserve(tree.leafNum - 3);
    tree.getTaxa(taxa);
    sort(taxa.begin(), taxa.end(), nodenamecmp);
    tree.createCluster(taxa, clusters);
    int cnt = 1;

    string treename = params.out_prefix;
    treename += ".clu-id";
    tree.printTree(treename.c_str());

    for (mmatrix(int)::iterator it = clusters.begin(); it != clusters.end(); it++, cnt++) {
        ostringstream filename;
        filename << params.out_prefix << "." << cnt << ".clu";
        ofstream out(filename.str().c_str());

        ostringstream filename2;
        filename2 << params.out_prefix << "." << cnt << ".name-clu";
        ofstream out2(filename2.str().c_str());

        out << "w" << endl << "c" << endl << "4" << endl << "b" << endl << "g" << endl << 4-params.is_rooted << endl;
        IntVector::iterator it2;
        NodeVector::iterator it3;
        for (it2 = (*it).begin(), it3 = taxa.begin(); it2 != (*it).end(); it2++, it3++)
            if ((*it3)->name != ROOT_NAME) {
                out << char((*it2)+'a') << endl;
                out2 << (*it3)->name << "  " << char((*it2)+'a') << endl;
            }
        out << "y" << endl;
        out.close();
        out2.close();
        cout << "Cluster " << cnt << " printed to " << filename.rdbuf() << " and " << filename2.rdbuf() << endl;
    }
}


void printTaxa(Params &params) {
    MTree mytree(params.user_file, params.is_rooted);
    vector<string> taxname;
    taxname.resize(mytree.leafNum);
    mytree.getTaxaName(taxname);
    sort(taxname.begin(), taxname.end());

    string filename = params.out_prefix;
    filename += ".taxa";

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename.c_str());
        for (vector<string>::iterator it = taxname.begin(); it != taxname.end(); it++) {
            if ((*it) != ROOT_NAME) out << (*it);
            out << endl;
        }
        out.close();
        cout << "All taxa names printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

void printAreaList(Params &params) {
    MSetsBlock *sets;
    sets = new MSetsBlock();
     cout << "Reading input file " << params.user_file << "..." << endl;

    MyReader nexus(params.user_file);

    nexus.Add(sets);

    MyToken token(nexus.inf);
    nexus.Execute(token);

    //sets->Report(cout);

    TaxaSetNameVector *allsets = sets->getSets();

    string filename = params.out_prefix;
    filename += ".names";

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename.c_str());
        for (TaxaSetNameVector::iterator it = allsets->begin(); it != allsets->end(); it++) {
            out << (*it)->name;
            out << endl;
        }
        out.close();
        cout << "All area names printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }

    delete sets;
}

void scaleBranchLength(Params &params) {
    params.is_rooted = true;
    PDTree tree(params);
    if (params.run_mode == SCALE_BRANCH_LEN) {
        cout << "Scaling branch length with a factor of " << params.scaling_factor << " ..." << endl;
        tree.scaleLength(params.scaling_factor, false);
    } else {
        cout << "Scaling clade support with a factor of " << params.scaling_factor << " ..." << endl;
        tree.scaleCladeSupport(params.scaling_factor, false);
    }
    if (params.out_file != NULL)
        tree.printTree(params.out_file);
    else {
        tree.printTree(cout);
        cout << endl;
    }
}

void calcDistribution(Params &params) {

    PDTree mytree(params);

    string filename = params.out_prefix;
    filename += ".randompd";

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename.c_str());
        for (int size = params.min_size; size <= params.sub_size; size += params.step_size) {
            out << size;
            for (int sample = 0; sample < params.sample_size; sample++) {
                Split taxset(mytree.leafNum);
                taxset.randomize(size);
                mytree.calcPD(taxset);
                out << "  " << taxset.getWeight();
            }
            out << endl;
        }
        out.close();
        cout << "PD distribution is printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

void printRFDist(string filename, double *rfdist, int n, int m, int rf_dist_mode, bool print_msg) {
    int i, j;

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        if (Params::getInstance().output_format == FORMAT_CSV) {
            out << "# Robinson-Foulds distances" << endl
            << "# This file can be read in MS Excel or in R with command:" << endl
            << "#    dat=read.csv('" <<  filename << "',comment.char='#')" << endl
            << "# Columns are comma-separated with following meanings:" << endl
            << "#    ID1:     Tree 1 ID" << endl
            << "#    ID2:     Tree 2 ID" << endl
            << "#    Dist:    Robinson-Foulds distance" << endl
            << "ID1,ID2,Dist" << endl;
            if (rf_dist_mode == RF_ADJACENT_PAIR) {
                for (i = 0; i < n; i++)
                    out << i+1 << ',' << i+2 << ',' << rfdist[i] << endl;
            } else if (Params::getInstance().rf_same_pair) {
                for (i = 0; i < n; i++)
                    out << i+1 << ',' << i+1 << ',' << rfdist[i] << endl;
            } else {
                for (i = 0; i < n; i++)  {
                    for (j = 0; j < m; j++)
                        out << i+1 << ',' << j+1 << ',' << rfdist[i*m+j] << endl;
                }
            }
        } else if (rf_dist_mode == RF_ADJACENT_PAIR || Params::getInstance().rf_same_pair) {
            out << "XXX        ";
            out << 1 << " " << n << endl;
            for (i = 0; i < n; i++)
                out << " " << rfdist[i];
            out << endl;
        } else {
            // all pairs
            out << n << " " << m << endl;
            for (i = 0; i < n; i++)  {
                out << "Tree" << i << "      ";
                for (j = 0; j < m; j++)
                    out << " " << rfdist[i*m+j];
                out << endl;
            }
        }
        out.close();
        if (print_msg)
            cout << "Robinson-Foulds distances printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

void computeRFDistExtended(const char *trees1, const char *trees2, const char *filename) {
    cout << "Reading input trees 1 file " << trees1 << endl;
    int ntrees = 0, ntrees2 = 0;
    double *rfdist_raw = NULL;
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(trees1);
        IntVector rfdist;
        for (ntrees = 1; !in.eof(); ntrees++) {
            MTree tree;
            bool is_rooted = false;

            // read in the tree and convert into split system for indexing
            tree.readTree(in, is_rooted);
            if (verbose_mode >= VB_DEBUG)
                cout << ntrees << " " << endl;
            DoubleVector dist;
            tree.computeRFDist(trees2, dist);
            ntrees2 = dist.size();
            rfdist.insert(rfdist.end(), dist.begin(), dist.end());
            char ch;
            in.exceptions(ios::goodbit);
            (in) >> ch;
            if (in.eof()) break;
            in.unget();
            in.exceptions(ios::failbit | ios::badbit);

        }

        in.close();
        ASSERT(ntrees * ntrees2 == rfdist.size());
        rfdist_raw = new double[rfdist.size()];
        copy(rfdist.begin(), rfdist.end(), rfdist_raw);

    } catch (ios::failure) {
        outError(ERR_READ_INPUT, trees1);
    }

    printRFDist(filename, rfdist_raw, ntrees, ntrees2, RF_TWO_TREE_SETS_EXTENDED);
    delete [] rfdist_raw;
}

void computeRFDistSamePair(const char *trees1, const char *trees2, const char *filename) {
    cout << "Reading input trees 1 file " << trees1 << endl;
    int ntrees = 0, ntrees2 = 0;
    double *rfdist_raw = NULL;
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(trees1);

        ifstream in2;
        in2.exceptions(ios::failbit | ios::badbit);
        in2.open(trees2);

        DoubleVector rfdist;
        for (ntrees = 1; !in.eof() && !in2.eof(); ntrees++) {
            MTree tree;
            bool is_rooted = false;
            // read in the tree and convert into split system for indexing
            tree.readTree(in, is_rooted);

            if (verbose_mode >= VB_DEBUG)
                cout << ntrees << " " << endl;
            DoubleVector dist;
            tree.computeRFDist(in2, dist, 0, true);
            ntrees2 = dist.size();
            rfdist.insert(rfdist.end(), dist.begin(), dist.end());
            char ch;
            in.exceptions(ios::goodbit);
            (in) >> ch;
            if (in.eof()) break;
            in.unget();
            in.exceptions(ios::failbit | ios::badbit);
            
        }
        
        in.close();
        in2.close();
        ASSERT(ntrees * ntrees2 == rfdist.size());
        rfdist_raw = new double[rfdist.size()];
        copy(rfdist.begin(), rfdist.end(), rfdist_raw);
        
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, trees1);
    }
    
    printRFDist(filename, rfdist_raw, ntrees, ntrees2, RF_TWO_TREE_SETS_EXTENDED);

    delete [] rfdist_raw;
}

void computeRFDist(Params &params) {

    if (!params.user_file) outError("User tree file not provided");

    string filename = params.out_prefix;
    filename += ".rfdist";

    if (params.rf_dist_mode == RF_TWO_TREE_SETS_EXTENDED) {
        computeRFDistExtended(params.user_file, params.second_tree, filename.c_str());
        return;
    }

    if (params.rf_same_pair) {
        computeRFDistSamePair(params.user_file, params.second_tree, filename.c_str());
        return;
    }

    MTreeSet trees(params.user_file, params.is_rooted, params.tree_burnin, params.tree_max_count);
    int n = trees.size(), m = trees.size();
    double *rfdist;
    double *incomp_splits = NULL;
    string infoname = params.out_prefix;
    infoname += ".rfinfo";
    string treename = params.out_prefix;
    treename += ".rftree";
    if (params.rf_dist_mode == RF_TWO_TREE_SETS) {
        MTreeSet treeset2(params.second_tree, params.is_rooted, params.tree_burnin, params.tree_max_count);
        cout << "Computing Robinson-Foulds distances between two sets of trees" << endl;
        m = treeset2.size();
        size_t size = n*m;
        if (params.rf_same_pair) {
            if (m != n)
                outError("Tree sets has different number of trees");
            size = n;
        }
        rfdist = new double [size];
        memset(rfdist, 0, size*sizeof(double));
        if (verbose_mode >= VB_MAX) {
            incomp_splits = new double [size];
            memset(incomp_splits, 0, size*sizeof(double));
        }
        if (verbose_mode >= VB_MED)
            trees.computeRFDist(rfdist, &treeset2, params.rf_same_pair,
                                infoname.c_str(),treename.c_str(), incomp_splits);
        else
            trees.computeRFDist(rfdist, &treeset2, params.rf_same_pair);
    } else {
        rfdist = new double [n*n];
        memset(rfdist, 0, n*n* sizeof(double));
        trees.computeRFDist(rfdist, params.rf_dist_mode, params.split_weight_threshold);
    }

    //if (verbose_mode >= VB_MED) printRFDist(cout, rfdist, n, m, params.rf_dist_mode);

    printRFDist(filename, rfdist, n, m, params.rf_dist_mode);

    if (incomp_splits) {
        filename = params.out_prefix;
        filename += ".incomp";
        printRFDist(filename, incomp_splits, n, m, params.rf_dist_mode, false);
        cout << "Number of incompatible splits in printed to " << filename << endl;
    }

    if (incomp_splits) delete [] incomp_splits;
    delete [] rfdist;
}


void testInputFile(Params &params) {
    SplitGraph sg(params);
    if (sg.isWeaklyCompatible())
        cout << "The split system is weakly compatible." << endl;
    else
        cout << "The split system is NOT weakly compatible." << endl;

}

/**MINH ANH: for some statistics about the branches on the input tree*/
void branchStats(Params &params){
    MaTree mytree(params.user_file, params.is_rooted);
    mytree.drawTree(cout,WT_TAXON_ID + WT_INT_NODE);
    //report to output file
    string output;
    if (params.out_file)
        output = params.out_file;
    else {
        if (params.out_prefix)
            output = params.out_prefix;
        else
            output = params.user_file;
        output += ".stats";
    }

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(output.c_str());
        mytree.printBrInfo(out);
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, output);
    }
    cout << "Information about branch lengths of the tree is printed to: " << output << endl;
    
    /***** Following added by BQM to print internal branch lengths */
    NodeVector nodes1, nodes2;
    mytree.generateNNIBraches(nodes1, nodes2);
    output = params.out_prefix;
    output += ".inlen";
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(output.c_str());
        for (int i = 0; i < nodes1.size(); i++)
            out << nodes1[i]->findNeighbor(nodes2[i])->length << " ";
        out << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, output);
    }
    cout << "Internal branch lengths printed to: " << output << endl;
}

/**MINH ANH: for comparison between the input tree and each tree in a given set of trees*/
void compare(Params &params){
    MaTree mytree(params.second_tree, params.is_rooted);
    //sort taxon names and update nodeID, to be consistent with MTreeSet
    NodeVector taxa;
    mytree.getTaxa(taxa);
    sort(taxa.begin(), taxa.end(), nodenamecmp);
    int i;
    NodeVector::iterator it;
    for (it = taxa.begin(), i = 0; it != taxa.end(); it++, i++)
            (*it)->id = i;

    string drawFile = params.second_tree;
    drawFile += ".draw";
    try {
        ofstream out1;
        out1.exceptions(ios::failbit | ios::badbit);
        out1.open(drawFile.c_str());
        mytree.drawTree(out1,WT_TAXON_ID + WT_INT_NODE);
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, drawFile);
    }
    cout << "Tree with branchID (nodeID) was printed to: " << drawFile << endl;


    MTreeSet trees(params.user_file,params.is_rooted, params.tree_burnin, params.tree_max_count);
    DoubleMatrix brMatrix;
    DoubleVector BSDs;
    IntVector RFs;
    mytree.comparedTo(trees, brMatrix, RFs, BSDs);
    int numTree = trees.size();
    int numNode = mytree.nodeNum;

    string output;
    if (params.out_file)
        output = params.out_file;
    else {
        if (params.out_prefix)
            output = params.out_prefix;
        else
            output = params.user_file;
        output += ".compare";
    }

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(output.c_str());
        //print the header
        out << "tree  " ;
        for (int nodeID = 0; nodeID < numNode; nodeID++ )
            if ( brMatrix[0][nodeID] != -2 )
                out << "br_" << nodeID << "  ";
        out << "RF  BSD" << endl;
        for ( int treeID = 0; treeID < numTree; treeID++ )
        {
            out << treeID << "  ";
            for (int nodeID = 0; nodeID < numNode; nodeID++ )
                if ( brMatrix[treeID][nodeID] != -2 )
                    out << brMatrix[treeID][nodeID] << "  ";
            out << RFs[treeID] << "  " << BSDs[treeID] << endl;
        }
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, output);
    }
    cout << "Comparison with the given set of trees is printed to: " << output << endl;
}

/**MINH ANH: to compute 'guided bootstrap' alignment*/
void guidedBootstrap(Params &params)
{
    MaAlignment inputAlign(params.aln_file,params.sequence_type, params.intype, params.model_name);
    inputAlign.readLogLL(params.siteLL_file);

    string outFre_name = params.out_prefix;
    outFre_name += ".patInfo";

    inputAlign.printPatObsExpFre(outFre_name.c_str());

    string gboAln_name = params.out_prefix;
    gboAln_name += ".gbo";

    MaAlignment gboAlign;
    double prob;
    gboAlign.generateExpectedAlignment(&inputAlign, prob);
    gboAlign.printAlignment(IN_PHYLIP, gboAln_name.c_str());


    string outProb_name = params.out_prefix;
    outProb_name += ".gbo.logP";
    try {
        ofstream outProb;
        outProb.exceptions(ios::failbit | ios::badbit);
        outProb.open(outProb_name.c_str());
        outProb.precision(10);
        outProb << prob << endl;
        outProb.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, outProb_name);
    }

    cout << "Information about patterns in the input alignment is printed to: " << outFre_name << endl;
    cout << "A 'guided bootstrap' alignment is printed to: " << gboAln_name << endl;
    cout << "Log of the probability of the new alignment is printed to: " << outProb_name << endl;
}

/**MINH ANH: to compute the probability of an alignment given the multinomial distribution of patterns frequencies derived from a reference alignment*/
void computeMulProb(Params &params)
{
    Alignment refAlign(params.second_align, params.sequence_type, params.intype, params.model_name);
    Alignment inputAlign(params.aln_file, params.sequence_type, params.intype, params.model_name);
    double prob;
    inputAlign.multinomialProb(refAlign,prob);
    //Printing
    string outProb_name = params.out_prefix;
    outProb_name += ".mprob";
    try {
        ofstream outProb;
        outProb.exceptions(ios::failbit | ios::badbit);
        outProb.open(outProb_name.c_str());
        outProb.precision(10);
        outProb << prob << endl;
        outProb.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, outProb_name);
    }
    cout << "Probability of alignment " << params.aln_file << " given alignment " << params.second_align << " is: " << prob << endl;
    cout << "The probability is printed to: " << outProb_name << endl;
}

void processNCBITree(Params &params) {
    NCBITree tree;
    Node *dad = tree.readNCBITree(params.user_file, params.ncbi_taxid, params.ncbi_taxon_level, params.ncbi_ignore_level);
    if (params.ncbi_names_file) tree.readNCBINames(params.ncbi_names_file);

    cout << "Dad ID: " << dad->name << " Root ID: " << tree.root->name << endl;
    string str = params.user_file;
    str += ".tree";
    if (params.out_file) str = params.out_file;
    //tree.printTree(str.c_str(), WT_SORT_TAXA | WT_BR_LEN);
    cout << "NCBI tree printed to " << str << endl;
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(str.c_str());
        tree.printTree(out, WT_SORT_TAXA | WT_BR_LEN | WT_TAXON_ID, tree.root, dad);
        out << ";" << endl;
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, str);
    }
}

outstreambuf* outstreambuf::open( const char* name, ios::openmode mode) {
    if (!(Params::getInstance().suppress_output_flags & OUT_LOG) && MPIHelper::getInstance().isMaster()) {
        fout.open(name, mode);
        if (!fout.is_open()) {
            cerr << "ERROR: Could not open " << name << " for logging" << endl;
            exit(EXIT_FAILURE);
            return NULL;
        }
        fout_buf = fout.rdbuf();
    }
    cout_buf = cout.rdbuf();
    cout.rdbuf(this);
    return this;
}

bool outstreambuf::is_open() {
    return fout.is_open();
}

outstreambuf* outstreambuf::close() {
    cout.rdbuf(cout_buf);
    if ( fout.is_open()) {
        sync();
        fout.close();
        return this;
    }
    return NULL;
}

int outstreambuf::overflow( int c) { // used for output buffer only
    if ((verbose_mode >= VB_MIN && MPIHelper::getInstance().isMaster()) || verbose_mode >= VB_MED)
        if (cout_buf->sputc(c) == EOF) return EOF;
    if (Params::getInstance().suppress_output_flags & OUT_LOG)
        return c;
    if (!MPIHelper::getInstance().isMaster())
        return c;
    if (fout_buf->sputc(c) == EOF) return EOF;
    return c;
}



int outstreambuf::sync() { // used for output buffer only
    if ((verbose_mode >= VB_MIN && MPIHelper::getInstance().isMaster()) || verbose_mode >= VB_MED)
        cout_buf->pubsync();
    if ((Params::getInstance().suppress_output_flags & OUT_LOG) || !MPIHelper::getInstance().isMaster())
        return 0;
    return fout_buf->pubsync();
}

extern "C" void startLogFile(bool append_log) {
    if (append_log)
        _out_buf.open(_log_file.c_str(), ios::app);
    else
        _out_buf.open(_log_file.c_str());
    _err_buf.init(_out_buf.get_fout_buf());
    _must_buf.init(_out_buf.get_cout_buf(), _out_buf.get_fout_buf());
}

extern "C" void endLogFile() {
    if (_out_buf.is_open())
        _out_buf.close();
}

void funcExit(void) {
    if(_exit_wait_optn) {
        printf("\npress [return] to finish: ");
        fflush(stdout);
        while (getchar() != '\n');
    }
    
    endLogFile();
    MPIHelper::getInstance().finalize();
}

extern "C" void funcAbort(int signal_number)
{
    /*Your code goes here. You can output debugging info.
      If you return from this function, and it was called
      because abort() was called, your program will exit or crash anyway
      (with a dialog box on Windows).
     */
#if (defined(__GNUC__) || defined(__clang__)) && !defined(WIN32) && !defined(WIN64) && !defined(__CYGWIN__)
    print_stacktrace(cerr);
#endif

    cerr << endl << "*** IQ-TREE CRASHES WITH SIGNAL ";
    switch (signal_number) {
        case SIGABRT: cerr << "ABORTED"; break;
        case SIGFPE:  cerr << "ERRONEOUS NUMERIC"; break;
        case SIGILL:  cerr << "ILLEGAL INSTRUCTION"; break;
        case SIGSEGV: cerr << "SEGMENTATION FAULT"; break;
#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__ && !defined WIN64
        case SIGBUS: cerr << "BUS ERROR"; break;
#endif
    }
    cerr << endl;
    cerr << "*** For bug report please send to developers:" << endl << "***    Log file: " << _log_file;
    cerr << endl << "***    Alignment files (if possible)" << endl;
    funcExit();
    signal(signal_number, SIG_DFL);
}

extern "C" void getintargv(int *argc, char **argv[])
{
    int    done;
    int    count;
    int    n;
    int    l;
    char   ch;
    char  *argtmp;
    char **argstr;

    argtmp = (char  *)calloc(10100, sizeof(char));
    argstr = (char **)calloc(100, sizeof(char*));
    for(n=0; n<100; n++) {
        argstr[n] = &(argtmp[n * 100]);
    }
    n=1;

    fprintf(stdout, "\nYou seem to have click-started this program,");
    fprintf(stdout, "\ndo you want to enter commandline parameters: [y]es, [n]o: ");
    fflush(stdout);

    /* read one char */
    ch = getc(stdin);
    if (ch != '\n') {
        do ;
        while (getc(stdin) != '\n');
    }
    ch = (char) tolower((int) ch);

    if (ch == 'y') {
        done=FALSE;

        fprintf(stdout, "\nEnter single parameter [! for none]: ");
        fflush(stdout);
        count = fscanf(stdin, "%s", argstr[n]);
        do ;
        while (getc(stdin) != '\n');

        if(argstr[0][0] == '!') {
            count = 0;
        } else {
            if (strlen(argstr[n]) > 100) {
                fprintf(stdout, "\nParameter too long!!!\n");
            } else {
                n++;
            }
        }

        while(!done) {
            fprintf(stdout, "\nCurrent commandline: ");
            for(l=1; l<n; l++) {
                fprintf(stdout, "%s ", argstr[l]);
            }
            fprintf(stdout, "\nQuit [q]; confirm [y]%s%s%s: ",
                (n<99 ? ", extend [e]" : ""),
                (n>1 ? ", delete last [l]" : ""),
                (n>1 ? ", delete all [a]" : ""));
            fflush(stdout);

            /* read one char */
            ch = getc(stdin);
            /* ch = getchar(); */
            if (ch != '\n') {
                do ;
                while (getc(stdin) != '\n');
                /* while (getchar() != '\n'); */
            }
            ch = (char) tolower((int) ch);
        
            switch (ch) {
                case 'y':
                    done=TRUE;
                    break;
                case 'e':
                    fprintf(stdout, "\nEnter single parameter [! for none]: ");
                    fflush(stdout);
                    count = fscanf(stdin, "%s", argstr[n]);
                    do ;
                    while (getc(stdin) != '\n');
        
                    if(argstr[0][0] == '!') {
                        count = 0;
                    } else {
                        if (strlen(argstr[n]) > 100) {
                            fprintf(stdout, "\nParameter too long!!!\n");
                        } else {
                            n++;
                        }
                    }
                    break;
                case 'l':
                    if (n>1) n--;
                    break;
                case 'a':
                    n=1;
                    break;
                case 'q':
                       // tp_exit(0, NULL, FALSE, __FILE__, __LINE__, _exit_wait_optn);
                    if(_exit_wait_optn) {
                        printf("\npress [return] to finish: ");
                        fflush(stdout);
                        while (getchar() != '\n');
                    }
                    exit(0);
                    break;
            }
        }
    }

    *argc = n;
    *argv = argstr;
} /* getintargv */

/*********************************************************************************************************************************
    Olga: ECOpd - phylogenetic diversity with ecological constraint: choosing a viable subset of species which maximizes PD/SD
*********************************************************************************************************************************/

void processECOpd(Params &params) {
    double startTime = getCPUTime();
    params.detected_mode = LINEAR_PROGRAMMING;
    cout<<"----------------------------------------------------------------------------------------"<<endl;
    int i;
    double score;
    double *variables;
    int threads = params.gurobi_threads;
    params.gurobi_format=true;

    string model_file,subFoodWeb,outFile;

    model_file = params.out_prefix;
    model_file += ".lp";

    subFoodWeb = params.out_prefix;
    subFoodWeb += ".subFoodWeb";

    outFile = params.out_prefix;
    outFile += ".pda";

    //Defining the input phylo type: t - rooted/unrooted tree, n - split network
    params.intype=detectInputFile(params.user_file);
    if(params.intype == IN_NEWICK){
        params.eco_type = "t";
    } else if(params.intype == IN_NEXUS){
        params.eco_type = "n";
    }

    // Checking whether to treat the food web as weighted or non weighted
    if(params.diet_max == 0){
        params.eco_weighted = false;
    }else if(params.diet_max > 100 || params.diet_max < 0){
        cout<<"The minimum percentage of the diet to be conserved for each predator"<<endl;
        cout<<"d = "<<params.diet_max<<endl;
        cout<<"ERROR: Wrong value of parameter d. It must be within the range 0 <= d <= 100"<<endl;
        exit(0);
    }else{
        params.eco_weighted = true;
    }

    if(strcmp(params.eco_type,"t")==0){
    /*--------------------------------- EcoPD Trees ---------------------------------*/
        ECOpd tree(params.user_file,params.is_rooted);

        // Setting all the information-----------------
        tree.phyloType = "t";
        tree.TaxaNUM = tree.leafNum;
        if(verbose_mode == VB_MAX){
            cout<<"TaxaNUM = "<<tree.TaxaNUM<<endl;
            cout<<"LeafNUM = "<<tree.leafNum<<endl;
            cout<<"root_id = "<<tree.root->id<<" root_name = "<<tree.root->name<<endl;

            for(i=0; i<tree.leafNum; i++){
                cout<<i<<" "<<tree.findNodeID(i)->name <<endl;
            }
        }

        //Getting Species Names from tree
        for(i = 0; i < tree.TaxaNUM; i++)
            (tree.phyloNames).push_back(tree.findNodeID(i)->name);
        //for(i=0;i<tree.phyloNames.size();i++)
        //    cout<<"["<<i<<"] "<<tree.phyloNames[i]<<endl;

        // Full species list including info from tree and food web. Here adding names from phyloInput.
        for(i=0; i<tree.TaxaNUM; i++)
            tree.names.push_back(&(tree.phyloNames[i]));

        // Read the taxa to be included in the final optimal subset
        if(params.initial_file)
            tree.readInitialTaxa(params.initial_file);

        // Read the DAG file, Synchronize species on the Tree and in the Food Web
        tree.weighted = params.eco_weighted;
        tree.T = params.diet_max*0.01;
        tree.readDAG(params.eco_dag_file);
        tree.defineK(params);

        // IP formulation
        cout<<"Formulating an IP problem..."<<endl;
        if(tree.rooted){
            tree.printECOlpRooted(model_file.c_str(),tree);
        } else {
            tree.printECOlpUnrooted(model_file.c_str(),tree);
        }

        // Solve IP problem
        cout<<"Solving the problem..."<<endl;
        variables = new double[tree.nvar];
        int g_return = gurobi_solve((char*)model_file.c_str(), tree.nvar, &score, variables, verbose_mode, threads);
        if(verbose_mode == VB_MAX){
            cout<<"GUROBI finished with "<<g_return<<" return."<<endl;
            for(i=0; i<tree.nvar; i++)
                cout<<"x"<<i<<" = "<<variables[i]<<endl;
            cout<<"score = "<<score<<endl;
        }
        tree.dietConserved(variables);
        params.run_time = getCPUTime() - startTime;
        tree.printResults((char*)outFile.c_str(),variables,score,params);
        tree.printSubFoodWeb((char*)subFoodWeb.c_str(),variables);
        delete[] variables;

    } else if(strcmp(params.eco_type,"n")==0){
    /*----------------------------- EcoPD SplitNetwork ------------------------------*/
        params.intype=detectInputFile(params.user_file);
        PDNetwork splitSYS(params);
        ECOpd ecoInfDAG;

        // Get the species names from SplitNetwork
        splitSYS.speciesList(&(ecoInfDAG.phyloNames));
        //for(i=0;i<ecoInfDAG.phyloNames.size();i++)
        //    cout<<"["<<i<<"] "<<ecoInfDAG.phyloNames[i]<<endl;

        ecoInfDAG.phyloType = "n";
        ecoInfDAG.TaxaNUM = splitSYS.getNTaxa();

        // Full species list including info from tree and food web
        for(i=0; i<ecoInfDAG.TaxaNUM; i++)
            ecoInfDAG.names.push_back(&(ecoInfDAG.phyloNames[i]));

        ecoInfDAG.weighted = params.eco_weighted;
        // Read the taxa to be included in the final optimal subset
        if(params.initial_file)
            ecoInfDAG.readInitialTaxa(params.initial_file);
        ecoInfDAG.T = params.diet_max*0.01;
        ecoInfDAG.readDAG(params.eco_dag_file);
        ecoInfDAG.defineK(params);

        cout<<"Formulating an IP problem..."<<endl;
        splitSYS.transformEcoLP(params, model_file.c_str(), 0);
        /**
         * (subset_size-4) - influences constraints for conserved splits.
         * should be less than taxaNUM in the split system.
         * With 0 prints all the constraints.
         * Values different of 0 reduce the # of constraints.
         **/

        ecoInfDAG.printInfDAG(model_file.c_str(),splitSYS,params);
        cout<<"Solving the problem..."<<endl;
        variables = new double[ecoInfDAG.nvar];
        int g_return = gurobi_solve((char*)model_file.c_str(), ecoInfDAG.nvar, &score, variables, verbose_mode, threads);
        if(verbose_mode == VB_MAX){
            cout<<"GUROBI finished with "<<g_return<<" return."<<endl;
            for(i=0; i<ecoInfDAG.nvar; i++)
                cout<<"x"<<i<<" = "<<variables[i]<<endl;
            cout<<"score = "<<score<<endl;
        }
        ecoInfDAG.splitsNUM = splitSYS.getNSplits();
        ecoInfDAG.totalSD = splitSYS.calcWeight();
        ecoInfDAG.dietConserved(variables);
        params.run_time = getCPUTime() - startTime;
        ecoInfDAG.printResults((char*)outFile.c_str(),variables, score,params);
        ecoInfDAG.printSubFoodWeb((char*)subFoodWeb.c_str(),variables);
        delete[] variables;
    }
}

void collapseLowBranchSupport(char *user_file, char *split_threshold_str) {
    DoubleVector minsup;
    convert_double_vec(split_threshold_str, minsup, '/');
    if (minsup.empty())
        outError("wrong -minsupnew argument, please use back-slash separated string");
    MExtTree tree;
    bool isrooted = false;
    tree.readTree(user_file, isrooted);
    tree.collapseLowBranchSupport(minsup);
    tree.collapseZeroBranches(NULL, NULL, -1.0);
    if (verbose_mode >= VB_MED)
        tree.drawTree(cout);
    string outfile = (string)user_file + ".collapsed";
    tree.printTree(outfile.c_str());
    cout << "Tree with collapsed branches written to " << outfile << endl;
}
