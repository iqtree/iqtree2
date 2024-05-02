/*
 *  terraceanalysis.cpp
 *  Interface to work with terraces
 *  Created on: Sep 10, 2020
 *      Author: Olga
 */

#define VAR_DECLS
#include "terrace/global_vars.h"
#include "terrace/ParallelContext.hpp"
#include "terraceanalysis.h"
#include "terrace/terracenode.hpp"
#include "terrace/terracetree.hpp"
#include "terrace/terrace.hpp"
#include "terrace/presenceabsencematrix.hpp"
#include "tree/mtreeset.h"
#include "alignment/superalignment.h"
#include "utils/timeutil.h"
#include <omp.h>
#include <thread>
#include <memory>
#include <chrono>

#define DYNAMIC_ORDER 1

void terrace_main(Params &params){

    params.startCPUTime = getCPUTime();
    params.start_real_time = getRealTime();

    pContext = make_shared_global<ParallelContext>(params.num_threads);

    if(params.num_threads > 1){

        int number_of_threads = params.num_threads;

        global_terrace_trees = 0;
        global_intermediate_trees = 0;
        global_dead_ends = 0;

        //std::vector<std::thread> vecOfThreads(number_of_threads);

        /* for(int i = 0; i<num_threads; i++){
            vecOfThreads[i] = std::thread(runterraceanalysis_parallel, std::ref(params),i+1);
        } */


        #pragma omp parallel for num_threads(number_of_threads)
        for(int i = 0; i<number_of_threads; i++){
            
            pContext->working[i] = true;
            int tid = omp_get_thread_num();
            runterraceanalysis_parallel(params, tid);
        }


        #pragma omp barrier

        if(!taskQueue.empty()){
            int taskQueueSize = taskQueue.size();
            
            for(int i = taskQueueSize-1; i>=0; i--){
                
                Task *task = taskQueue.front();
                taskQueue.pop();

                delete task->list_taxa_to_insert;
                for(int i = 0; i<task->path_size; i++){
                    delete[] task->path_up_to_now[i];
                }

                delete[] task->path_up_to_now;
            }
        }
        
    } else {

        global_terrace_trees = 0;
        global_intermediate_trees = 0;
        global_dead_ends = 0;

        runterraceanalysis(params);
    }

}

void runterraceanalysis(Params &params){

    cout<<"---------------------------------------------------------\n";
    cout<<"\nBegin analysis with Gentrius...\n\n";
    cout<<"---------------------------------------------------------\n";

    if(params.matrix_order){
        cout<<"Creating presence-absence matrix..\n\n";
        std::shared_ptr<PresenceAbsenceMatrix> matrix = make_shared_global<PresenceAbsenceMatrix>();
        if(params.partition_file && params.aln_file){
            matrix->get_from_alignment(params);
        }else{
            assert(params.pr_ab_matrix && "ERROR: no input presence-absence matrix!");
            matrix->read_pr_ab_matrix(params.pr_ab_matrix);
        }

        string out_file;
        ofstream out;
        out_file = params.out_prefix;
        out_file += ".pr_ab_matrix";
        out.open(out_file);
        matrix->print_pr_ab_matrix(out);
        out.close();

        if(params.print_m_overlap){
            matrix->getPartOverlapComplete();

            string out_file;
            ofstream out;
            out_file = params.out_prefix;
            out_file += ".part_overlap";
            out.open(out_file);
            matrix->print_overlap_matrix(out);
            out.close();
        }
        cout<<"---------------------------------------------------------\n";
        cout<<"Presence-absence matrix was written into "<<out_file<<"\n";
        cout<<"---------------------------------------------------------\n";
        cout<<"Total wall-clock time used: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")\n";
        cout<<"Total CPU time used: "<< getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")\n\n";

        time_t end_time;
        time(&end_time);
        cout << "Date and Time: " << ctime(&end_time);
        exit(0);
    }


    /*
    *  This is an actual terrace:
    *  - main tree - is terrace representative tree
    *  - induced partition trees define the terrace
    */
    Terrace* terrace;
    
    if(params.partition_file && params.user_file){
        std::shared_ptr<PresenceAbsenceMatrix> matrix = make_shared_global<PresenceAbsenceMatrix>();
        matrix->get_from_alignment(params);
        TerraceTree tree;
        tree.readTree(params.user_file, params.is_rooted);
        if(tree.rooted){
            cout<<"WARNING: Gentrius is only used with unrooted trees!\nConverting rooted tree to unrooted...\n";
            tree.convertToUnrooted();
        }
        if(!tree.isBifurcating()){
            cout<<"ERROR: The tree has multifurcations! Only bifurcating trees are allowed. Exiting...\n";
            exit(0);
        }
        terrace = new Terrace(tree,matrix);
        params.print_pr_ab_matrix=true;
        cout<<"\n---------------------------------------------------------\n";
    } else if(params.user_file && params.pr_ab_matrix){
        terrace = new Terrace(params.user_file,params.is_rooted,params.pr_ab_matrix);
        if(!terrace->isBifurcating()){
            cout<<"ERROR: The tree has multifurcations! Only bifurcating trees are allowed. Exiting...\n";
            exit(0);
        }
    } else if(params.user_file){
        vector<TerraceTree*> subtrees;
        read_tree_set(params.user_file, params.is_rooted, subtrees);
        assert(subtrees.size()>1 && "ERROR: the input set of subtrees must be larger than 1 to be considered for the analysis.");
        if(subtrees[0]->rooted){
            cout<<"WARNING: Gentrius is only used with unrooted trees!\n";
        }
        int ti=0;
        for(const auto &t: subtrees){
            ti++;
            if(t->rooted){
                cout<<"Tree "<<ti<<": converting rooted tree to unrooted...\n";
                t->convertToUnrooted();
                //t->printTree(cout);
            }
            if(!t->isBifurcating()){
                cout<<"ERROR: Tree "<<ti<<" has multifurcations! Only bifurcating trees are allowed. Exiting...\n";
                exit(0);
            }
        }
        terrace = new Terrace(subtrees);
    }else{
        throw "ERROR: to start terrace analysis input a tree and either a presence-absence matrix or alignment with partition info or a set of overlapping unrooted trees!";
    }


    if(params.terrace_query_set){
        run_terrace_check(terrace,params);
    
    }else{
        //cout<<"---------------------------------------------------------\n";
        cout<<"Notes:\n";
        cout<<"- Gentrius generates a stand, i.e. a collection of species-trees, which display (also compatible with) the same loci subtrees\n";
        cout<<"- If one of the stopping rules is triggered, the stand was not generated completely. You can change the thresholds or even turn them off to try generating all trees.\n";
        cout<<"- If with default thresholds Gentrius did not generate any species-tree, use alternative exploratory approach (see manual).\n";
        cout<<"---------------------------------------------------------\n";
        
        /*
        *  Create an auxiliary terrace:
        *  - main tree - is the initial tree to be expanded by adding taxa to obtain a tree from the considered terrace (or to reach dead end)
        *  - induced trees - the common subtrees of the main tree and higher level induced partition tree, respectivelly per partition
        *
        *  There will be also a vector of terraces with just one partition. Per partition each terrace is a pair of high and low (common) induced partition trees
        *  - main tree - high level partition tree
        *  - induced tree - a common subtree between "initial to be expanded main tree" and high level induced tree
        */

        //terrace->printInfo();
        //terrace->matrix->print_pr_ab_matrix();

        if(params.print_induced_trees){
            string out_file;
            ofstream out;
            out_file = params.out_prefix;
            out_file += ".stand_info";
            out.open(out_file);
            terrace->printInfo(out);
            out.close();
        }
        if(params.print_pr_ab_matrix){
            string out_file;
            ofstream out;
            out_file = params.out_prefix;
            out_file += ".pr_ab_matrix";
            out.open(out_file);
            terrace->matrix->print_pr_ab_matrix(out);
            out.close();
        }

        terrace->rm_leaves = params.terrace_remove_m_leaves;
        terrace->start_real_time = params.start_real_time;

        /* if(params.num_threads <4){
            terrace->taskThreshold = params.num_threads + 1;
        } else if (params.num_threads >= 4 && params.num_threads < 8){
            terrace->taskThreshold = params.num_threads - 1;
        } else {
            terrace->taskThreshold = params.num_threads/2;
        }
         */
        //terrace->taskThreshold = params.num_threads + 1;

        //Terrace* test = new Terrace(terrace->induced_trees);
        //test->rm_leaves = 0;
        //run_generate_trees(test, params,test->rm_leaves);
        run_generate_trees(terrace, params,terrace->rm_leaves);
    }

    //pthread_cancel(pthread_self());
    
};

void runterraceanalysis_parallel(Params &params, int thread_num){
    
    if(thread_num == 0){
        cout<<"---------------------------------------------------------\n";
        cout<<"\nBegin analysis with Gentrius...\n\n";
        cout<<"---------------------------------------------------------\n";
    }

    if(params.matrix_order){
        
        if(thread_num == 0){
            cout<<"Creating presence-absence matrix..\n\n";
            std::shared_ptr<PresenceAbsenceMatrix> matrix = make_shared_global<PresenceAbsenceMatrix>();
            if(params.partition_file && params.aln_file){
                matrix->get_from_alignment(params);
            }else{
                assert(params.pr_ab_matrix && "ERROR: no input presence-absence matrix!");
                matrix->read_pr_ab_matrix(params.pr_ab_matrix);
            }

            string out_file;
            ofstream out;
            out_file = params.out_prefix;
            out_file += ".pr_ab_matrix";
            out.open(out_file);
            matrix->print_pr_ab_matrix(out);
            out.close();

            if(params.print_m_overlap){
                matrix->getPartOverlapComplete();

                string out_file;
                ofstream out;
                out_file = params.out_prefix;
                out_file += ".part_overlap";
                out.open(out_file);
                matrix->print_overlap_matrix(out);
                out.close();
            }
            cout<<"---------------------------------------------------------\n";
            cout<<"Presence-absence matrix was written into "<<out_file<<"\n";
            cout<<"---------------------------------------------------------\n";
            cout<<"Total wall-clock time used: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")\n";
            cout<<"Total CPU time used: "<< getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")\n\n";

            time_t end_time;
            time(&end_time);
            cout << "Date and Time: " << ctime(&end_time);
        }

       return;

    }


    /*
    *  This is an actual terrace:
    *  - main tree - is terrace representative tree
    *  - induced partition trees define the terrace
    */
    Terrace* terrace;
    
    //terrace->artificial_thread_num = 0;

    if(params.partition_file && params.user_file){
        std::shared_ptr<PresenceAbsenceMatrix> matrix = make_shared_global<PresenceAbsenceMatrix>();
        matrix->get_from_alignment(params);
        TerraceTree tree;
        tree.readTree(params.user_file, params.is_rooted);
        if(tree.rooted){
            if(thread_num == 0)
                cout<<"WARNING: Gentrius is only used with unrooted trees!\nConverting rooted tree to unrooted...\n";
            
            tree.convertToUnrooted();
        }
        if(!tree.isBifurcating()){
            if(thread_num == 0)
                cout<<"ERROR: The tree has multifurcations! Only bifurcating trees are allowed. Exiting...\n";
            
            return;
        }
        terrace = new Terrace(tree,matrix);
        params.print_pr_ab_matrix=true;
        if(thread_num == 0)
            cout<<"\n---------------------------------------------------------\n";
        
    } else if(params.user_file && params.pr_ab_matrix){
        terrace = new Terrace(params.user_file,params.is_rooted,params.pr_ab_matrix);
        if(!terrace->isBifurcating()){
            if(thread_num == 0)    
                cout<<"ERROR: The tree has multifurcations! Only bifurcating trees are allowed. Exiting...\n";
            
            return;
        }
    } else if(params.user_file){
        vector<TerraceTree*> subtrees;
        read_tree_set(params.user_file, params.is_rooted, subtrees);
        assert(subtrees.size()>1 && "ERROR: the input set of subtrees must be larger than 1 to be considered for the analysis.");
        if(subtrees[0]->rooted){
            if(thread_num == 0)
                cout<<"WARNING: Gentrius is only used with unrooted trees!\n";
        }
        int ti=0;
        for(const auto &t: subtrees){
            ti++;
            if(t->rooted){
                if(thread_num == 0)
                    cout<<"Tree "<<ti<<": converting rooted tree to unrooted...\n";
                
                t->convertToUnrooted();
                //t->printTree(cout);
            }
            if(!t->isBifurcating()){
                if(thread_num == 0)
                    cout<<"ERROR: Tree "<<ti<<" has multifurcations! Only bifurcating trees are allowed. Exiting...\n";
                
                return;
            }
        }
        terrace = new Terrace(subtrees);
    }else{
        throw "ERROR: to start terrace analysis input a tree and either a presence-absence matrix or alignment with partition info or a set of overlapping unrooted trees!";
    }


    if(params.terrace_query_set){
        if(thread_num == 0)
            run_terrace_check(terrace,params);
    
    }else{
        //cout<<"---------------------------------------------------------\n";
        if(thread_num == 0){    
            cout<<"Notes:\n";
            cout<<"- Gentrius generates a stand, i.e. a collection of species-trees, which display (also compatible with) the same loci subtrees\n";
            cout<<"- If one of the stopping rules is triggered, the stand was not generated completely. You can change the thresholds or even turn them off to try generating all trees.\n";
            cout<<"- If with default thresholds Gentrius did not generate any species-tree, use alternative exploratory approach (see manual).\n";
            cout<<"---------------------------------------------------------\n";
        }    
        /*
        *  Create an auxiliary terrace:
        *  - main tree - is the initial tree to be expanded by adding taxa to obtain a tree from the considered terrace (or to reach dead end)
        *  - induced trees - the common subtrees of the main tree and higher level induced partition tree, respectivelly per partition
        *
        *  There will be also a vector of terraces with just one partition. Per partition each terrace is a pair of high and low (common) induced partition trees
        *  - main tree - high level partition tree
        *  - induced tree - a common subtree between "initial to be expanded main tree" and high level induced tree
        */

        //terrace->printInfo();
        //terrace->matrix->print_pr_ab_matrix();

        if(params.print_induced_trees && thread_num == 0){
            string out_file;
            ofstream out;
            out_file = params.out_prefix;
            out_file += ".stand_info";
            out.open(out_file);
            terrace->printInfo(out);
            out.close();
        }
        if(params.print_pr_ab_matrix && thread_num == 0){
            string out_file;
            ofstream out;
            out_file = params.out_prefix;
            out_file += ".pr_ab_matrix";
            out.open(out_file);
            terrace->matrix->print_pr_ab_matrix(out);
            out.close();
        }

        terrace->rm_leaves = params.terrace_remove_m_leaves;
        terrace->start_real_time = params.start_real_time;

        if(params.num_threads < 8){
            terrace->taskThreshold = params.num_threads + 1;
        } else {
            terrace->taskThreshold = params.num_threads/2;
        }

        terrace->remaining_threads_to_assign = params.num_threads;
        
        run_generate_trees_parallel(terrace, params,terrace->rm_leaves, thread_num);
    }
}



void run_generate_trees(Terrace* terrace, Params &params,const int m){

    assert(m>=0 && "ERROR: required number m of leaves to be removed from the input tree is incorrect! m must be >=0. Exiting...");
    if(m>=0){
        terrace->rm_leaves = m;
    }
    // GET INFO FOR INITIAL TREE and TAXA TO INSERT
    vector<string> taxa_names_sub;      // taxa to appear on initial tree
    vector<string> list_taxa_to_insert; // list of taxa to be inserted.

    // partition print
    int part_init = terrace->matrix->getINFO_init_tree_taxon_order(taxa_names_sub,list_taxa_to_insert,terrace->rm_leaves);
    if(part_init == -1){
        global_terrace_trees = 1;
        return;
    }

    // INITIAL TERRACE
    // - to be used for generating trees
    // - agile master tree is the initial tree
    // - induced partition trees are common subtrees with induced partition trees of considered terrace

    // Creating a subterrace: submatrix and an initial tree
    // the submatrix contains only rows that have 1 in the selected partition -> (if part_init = 6, column 6 will have only 1's)
    std::shared_ptr<PresenceAbsenceMatrix> submatrix = make_shared_global<PresenceAbsenceMatrix>();
    terrace->matrix->getSubPrAbMatrix(taxa_names_sub, submatrix.get());

    TerraceTree tree_init;
    if(terrace->root){
        tree_init.copyTree_byTaxonNames(terrace,taxa_names_sub);
    }else{
        assert(part_init!=-1);
        tree_init.copyTree(terrace->induced_trees[part_init]);
    }
        //tree_init.drawTree(cout, WT_BR_SCALE | WT_TAXON_ID | WT_NEWLINE);

    Terrace* init_terrace = new Terrace(tree_init, submatrix);
    init_terrace->pContext = pContext;
    init_terrace->real_thread_num = 0;
    init_terrace->artificial_thread_num = 0;
    init_terrace->parallel_exec = false;
    init_terrace->initial_split_done = true;
    init_terrace->taskThreshold = terrace->taskThreshold;
    
    init_terrace->rm_leaves = m;
    init_terrace->master_terrace = terrace;

    init_terrace->terrace_trees_num = 0;
    init_terrace->intermediated_trees_num = 0;
    init_terrace->dead_ends_num = 0;
    
    init_terrace->cur_terrace_trees = 0;
    init_terrace->cur_intermediate_trees = 0;

    // times
    init_terrace->start_real_time = terrace->start_real_time;
    init_terrace->terrace_trees_time = init_terrace->start_real_time;
    init_terrace->intermediate_trees_time = init_terrace->start_real_time;
    init_terrace->dead_ends_time = init_terrace->start_real_time;
     
    //init_terrace->printInfo();
    //init_terrace->matrix->print_pr_ab_matrix();

    init_terrace->out_file = params.out_prefix;
    init_terrace->out_file += ".stand_trees";

    if(params.terrace_non_stop){
        cout<<"All stopping rules for stand generation are turned off.\n";
        cout<<"For the number of intermediate trees and stand trees the hard stopping rule is the maximum value of unsigned integer: "<<UINT_MAX<<"\n";
        init_terrace->intermediate_max_trees = UINT_MAX;
        init_terrace->terrace_max_trees = UINT_MAX;
        init_terrace->seconds_max = -1;
    } else {
        if(params.terrace_stop_intermediate_num > 0){
            init_terrace->intermediate_max_trees = params.terrace_stop_intermediate_num;
        } else if(params.terrace_stop_intermediate_num == 0){
            cout<<"The stopping rule based on the number of intermediate trees: the threshold is set to the maximum value of unsigned integer: "<<UINT_MAX<<"\n";
            init_terrace->intermediate_max_trees = UINT_MAX;
        }
        if(params.terrace_stop_terrace_trees_num > 0){
            init_terrace->terrace_max_trees = params.terrace_stop_terrace_trees_num;
        } else if(params.terrace_stop_terrace_trees_num == 0){
            cout<<"The stopping rule based on the size of a stand: the threshold is set to the maximum value of unsigned integer: "<<UINT_MAX<<"\n";
            init_terrace->terrace_max_trees = UINT_MAX;
        }

        if(params.terrace_stop_time > 0){
            init_terrace->seconds_max = params.terrace_stop_time*3600;
        } else if(params.terrace_stop_time == 0){
            cout<<"The stopping rule based on the CPU time limit is turned off.\n";
            init_terrace->seconds_max = -1;
        }
    }

    cout<<"---------------------------------------------------------"<<"\n";
    cout<<"Current stopping thresholds:\n";
    cout<<"1. Stop if stand size reached: "<<init_terrace->terrace_max_trees<<"\n";
    cout<<"2. Stop if the number of intermediate visited trees reached: "<<init_terrace->intermediate_max_trees<<"\n";
    cout<<"3. Stop if the CPU time reached: "<<init_terrace->seconds_max<<" seconds"<<"\n";


    // links tree_init to the induced trees 
    init_terrace->trees_out_lim = params.terrace_print_lim;
    init_terrace->linkTrees(true, false); // branch_back_map, taxon_back_map; in this case you only want to map branches

    // links induced trees in init_terrace with induced trees in terrace (master terrace)
    vector<Terrace*> part_tree_pairs;
    init_terrace->create_Top_Low_Part_Tree_Pairs(part_tree_pairs, terrace);

    /*cout<<"Taxa to insert on initial tree (initial order): "<<"\n";
    for(i=0; i<list_taxa_to_insert.size(); i++){
        cout<<i<<":"<<list_taxa_to_insert[i]<<"\n";
    }*/

    init_terrace->fillbrNodes();
    cout<<"---------------------------------------------------------"<<"\n";
    cout<<"\n"<<"READY TO GENERATE TREES FROM A STAND"<<"\n";
    if(terrace->rm_leaves>0){
        cout<<"using BACKWARD approach (exploratory only! does not guarantee generating all trees)\n";
    }else{
        cout<<"using FORWARD approach (generates all trees, can be exponentially many)\n";
    }

    cout<<"\n"<<"---------------------------------------------------------"<<"\n";
    cout<<"INPUT INFO:"<<"\n";
    cout<<"---------------------------------------------------------"<<"\n";
    cout<<"Number of taxa: "<<terrace->taxa_num<<"\n";
    cout<<"Number of loci: "<<terrace->part_num<<"\n";
    cout<<"Number of taxa with minimal coverage (row sum = 1): "<<terrace->matrix->uniq_taxa_num<<"\n";
    terrace->matrix->percent_missing();
    cout<<"% of missing entries in taxon per locus presence-absence matrix: "<<terrace->matrix->missing_percent<<"\n";
    cout<<"Number of taxa on initial tree: "<<init_terrace->taxa_num<<"\n";
    cout<<"Number of taxa to be inserted: "<< list_taxa_to_insert.size()<<"\n";
    cout<<"---------------------------------------------------------"<<"\n";

    //init_terrace->print_ALL_DATA(part_tree_pairs);
    if(params.print_terrace_trees){
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(init_terrace->out_file);
        out.close();

        init_terrace->out.exceptions(ios::failbit | ios::badbit);
        init_terrace->out.open(init_terrace->out_file,std::ios_base::app);

    }else{
        init_terrace->terrace_out = false;
    }

    init_terrace->matrix->uniq_taxa_num = terrace->matrix->uniq_taxa_num;
    init_terrace->matrix->uniq_taxa_to_insert_num = terrace->matrix->uniq_taxa_to_insert_num;
    cout<<"\n"<<"Generating trees from a stand...."<<"\n";
    
    
    bool use_dynamic_taxon_order = (bool)DYNAMIC_ORDER;
    bool thread_call = false;
    int taxon_to_insert = 0;
    std::string taxon_name = "";
    int thread_num = 0;
    
    generateTerrace_wrapper(init_terrace,
                            part_tree_pairs, 
                            list_taxa_to_insert,
                            taxon_to_insert, 
                            thread_call, 
                            thread_num, 
                            taxon_name,
                            use_dynamic_taxon_order);
    
    
    std::for_each(part_tree_pairs.begin(), part_tree_pairs.end(), delete_pointer_element<Terrace*>() );
    //(*part_tree_pairs).clear();
    
    //list_taxa_to_insert->clear();
    delete init_terrace;

    //std::for_each( part_tree_pairs.begin(), part_tree_pairs.end(), delete_pointer_element<Terrace*>() );
    
};

void run_generate_trees_parallel(Terrace* terrace, Params &params, const int m, int thread_num){
    
    assert(m>=0 && "ERROR: required number m of leaves to be removed from the input tree is incorrect! m must be >=0. Exiting...");
    if(m>=0){
        terrace->rm_leaves = m;
    }
    // GET INFO FOR INITIAL TREE and TAXA TO INSERT
    vector<string> taxa_names_sub;      // taxa to appear on initial tree
    vector<string> list_taxa_to_insert; // list of taxa to be inserted.

    // partition print
    int part_init = terrace->matrix->getINFO_init_tree_taxon_order(taxa_names_sub, list_taxa_to_insert,terrace->rm_leaves, thread_num);
    if(part_init == -1){
        if(thread_num == 0) global_terrace_trees = 1;
        return;
    }
    // INITIAL TERRACE
    // - to be used for generating trees
    // - agile master tree is the initial tree
    // - induced partition trees are common subtrees with induced partition trees of considered terrace

    // Creating a subterrace: submatrix and an initial tree
    // the submatrix contains only rows that have 1 in the selected partition -> (if part_init = 6, column 6 will have only 1's)
    std::shared_ptr<PresenceAbsenceMatrix> submatrix = make_shared_global<PresenceAbsenceMatrix>();
    terrace->matrix->getSubPrAbMatrix(taxa_names_sub, submatrix.get());

    TerraceTree tree_init;
    if(terrace->root){
        tree_init.copyTree_byTaxonNames(terrace,taxa_names_sub);
    }else{
        assert(part_init!=-1);
        tree_init.copyTree(terrace->induced_trees[part_init]);
    }
        //tree_init.drawTree(cout, WT_BR_SCALE | WT_TAXON_ID | WT_NEWLINE);

    Terrace* init_terrace = new Terrace(tree_init, submatrix);
    init_terrace->start_real_time = terrace->start_real_time;
    init_terrace->pContext = pContext;
    init_terrace->real_thread_num = thread_num;
    init_terrace->artificial_thread_num = thread_num;
    init_terrace->parallel_exec = true;
    init_terrace->initial_split_done = false;
    init_terrace->initial_split_done_index = -1;
    init_terrace->remaining_threads_to_assign = terrace->remaining_threads_to_assign;
    init_terrace->taskThreshold = terrace->taskThreshold;

    init_terrace->rm_leaves = m;
    init_terrace->master_terrace = terrace;

    init_terrace->terrace_trees_num = 0;
    init_terrace->intermediated_trees_num = 0;
    init_terrace->dead_ends_num = 0;
    
    init_terrace->cur_terrace_trees = 0;
    init_terrace->cur_intermediate_trees = 0;
    
    //init_terrace->printInfo();
    //init_terrace->matrix->print_pr_ab_matrix();

    init_terrace->out_file = params.out_prefix;
    init_terrace->out_file += ".stand_trees";

    if(params.terrace_non_stop){
        if(thread_num == 0){
            cout<<"All stopping rules for stand generation are turned off.\n";
            cout<<"For the number of intermediate trees and stand trees the hard stopping rule is the maximum value of unsigned integer: "<<UINT_MAX<<"\n";
        }
        init_terrace->intermediate_max_trees = UINT_MAX;
        init_terrace->terrace_max_trees = UINT_MAX;
        init_terrace->seconds_max = -1;
    } else {
        if(params.terrace_stop_intermediate_num > 0){
            init_terrace->intermediate_max_trees = params.terrace_stop_intermediate_num;
        } else if(params.terrace_stop_intermediate_num == 0){
            if(thread_num == 0)
                cout<<"The stopping rule based on the number of intermediate trees: the threshold is set to the maximum value of unsigned integer: "<<UINT_MAX<<"\n";
            
            init_terrace->intermediate_max_trees = UINT_MAX;
        }
        if(params.terrace_stop_terrace_trees_num > 0){
            init_terrace->terrace_max_trees = params.terrace_stop_terrace_trees_num;
        } else if(params.terrace_stop_terrace_trees_num == 0){
            if(thread_num == 0)
                cout<<"The stopping rule based on the size of a stand: the threshold is set to the maximum value of unsigned integer: "<<UINT_MAX<<"\n";
            
            init_terrace->terrace_max_trees = UINT_MAX;
        }

        if(params.terrace_stop_time > 0){
            init_terrace->seconds_max = params.terrace_stop_time*3600;
        } else if(params.terrace_stop_time == 0){
            if(thread_num == 0)
                cout<<"The stopping rule based on the CPU time limit is turned off.\n";
            
            init_terrace->seconds_max = -1;
        }
    }

    if(thread_num == 0){
        cout<<"---------------------------------------------------------"<<"\n";
        cout<<"Current stopping thresholds:\n";
        cout<<"1. Stop if stand size reached: "<<init_terrace->terrace_max_trees<<"\n";
        cout<<"2. Stop if the number of intermediate visited trees reached: "<<init_terrace->intermediate_max_trees<<"\n";
        cout<<"3. Stop if the CPU time reached: "<<init_terrace->seconds_max<<" seconds"<<"\n";
    }

    // links tree_init to the induced trees 
    init_terrace->trees_out_lim = params.terrace_print_lim;
    init_terrace->linkTrees(true, false); // branch_back_map, taxon_back_map; in this case you only want to map branches

    // links induced trees in init_terrace with induced trees in terrace (master terrace)
    vector<Terrace*> part_tree_pairs;
    init_terrace->create_Top_Low_Part_Tree_Pairs(part_tree_pairs, terrace);


    init_terrace->fillbrNodes();
    if(thread_num == 0){
        cout<<"---------------------------------------------------------"<<"\n";
        cout<<"\n"<<"READY TO GENERATE TREES FROM A STAND"<<"\n";
        if(terrace->rm_leaves>0){
            cout<<"using BACKWARD approach (exploratory only! does not guarantee generating all trees)\n";
        }else{
            cout<<"using FORWARD approach (generates all trees, can be exponentially many)\n";
        }


        cout<<"\n"<<"---------------------------------------------------------"<<"\n";
        cout<<"INPUT INFO:"<<"\n";
        cout<<"---------------------------------------------------------"<<"\n";
        cout<<"Number of taxa: "<<terrace->taxa_num<<"\n";
        cout<<"Number of loci: "<<terrace->part_num<<"\n";
        cout<<"Number of taxa with minimal coverage (row sum = 1): "<<terrace->matrix->uniq_taxa_num<<"\n";
        terrace->matrix->percent_missing();
        cout<<"% of missing entries in taxon per locus presence-absence matrix: "<<terrace->matrix->missing_percent<<"\n";
        cout<<"Number of taxa on initial tree: "<<init_terrace->taxa_num<<"\n";
        cout<<"Number of taxa to be inserted: "<< list_taxa_to_insert.size()<<"\n";
        cout<<"---------------------------------------------------------"<<"\n";
    }

    //init_terrace->print_ALL_DATA(part_tree_pairs);
    if(params.print_terrace_trees){
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(init_terrace->out_file);
        out.close();

        init_terrace->out.exceptions(ios::failbit | ios::badbit);
        init_terrace->out.open(init_terrace->out_file,std::ios_base::app);

    }else{
        init_terrace->terrace_out = false;
    }

    init_terrace->matrix->uniq_taxa_num = terrace->matrix->uniq_taxa_num;
    init_terrace->matrix->uniq_taxa_to_insert_num = terrace->matrix->uniq_taxa_to_insert_num;
    
    if(thread_num == 0)
        cout<<"\n"<<"Generating trees from a stand...."<<"\n";
    
    bool use_dynamic_taxon_order = (bool)DYNAMIC_ORDER;
    bool thread_call = false;
    int taxon_to_insert = 0;
    std::string taxon_name = "";

    generateTerrace_wrapper(init_terrace,
                            part_tree_pairs, 
                            list_taxa_to_insert,
                            taxon_to_insert, 
                            thread_call, 
                            thread_num, 
                            taxon_name,
                            use_dynamic_taxon_order);
    
    
    /* for(int i=list_taxa_to_insert.size()-1; i>=0; i--){
        init_terrace->remove_one_taxon(list_taxa_to_insert[i], part_tree_pairs);
    } */

    startThread(init_terrace, &part_tree_pairs, thread_num);

    std::for_each(part_tree_pairs.begin(), part_tree_pairs.end(), delete_pointer_element<Terrace*>());
    delete init_terrace;

}

void generateTerrace_wrapper(Terrace* init_terrace,
                            vector<Terrace*> &part_tree_pairs, 
                            vector<string> &list_taxa_to_insert,
                            int &taxon_to_insert, 
                            bool &thread_call, 
                            int &new_thread_num, 
                            string &taxon_name,
                            bool &use_dynamic_taxon_order, 
                            vector<int>* ids1,
                            vector<int>* ids2)

{
    //init_terrace->threshold = std::min((0.25*list_taxa_to_insert.size()), 10.);

    if(thread_call){

        init_terrace->generateTerraceTrees(init_terrace->master_terrace, 
                                                part_tree_pairs, 
                                                list_taxa_to_insert, 
                                                taxon_to_insert, 
                                                use_dynamic_taxon_order, 
                                                true,
                                                taxon_name,
                                                ids1, 
                                                ids2);
        
    } else{

        //pContext->working[new_thread_num] = true;
        #pragma omp barrier
        
        // The taxon order is dynamicly adapted using information about allowed branches
        init_terrace->generateTerraceTrees(init_terrace->master_terrace, 
                                            part_tree_pairs, 
                                            list_taxa_to_insert, 
                                            taxon_to_insert, 
                                            use_dynamic_taxon_order, 
                                            thread_call,
                                            "");

    }
    
    global_terrace_trees += init_terrace->terrace_trees_num;
    global_intermediate_trees += init_terrace->intermediated_trees_num;
    global_dead_ends += init_terrace->dead_ends_num;

    {   
        omp_set_lock(&pContext->thread_mutex);
        pContext->working[new_thread_num] = false;
        //std::cout << "Complete thread " << new_thread_num << std::endl;
        
        if(pContext->CheckIfFinished() && pContext->taskCount == 0){
            if(!pContext->complete){
                
                std::cout<<"\n"<<"---------------------------------------------------------"<<"\n";
                std::cout<<"\n"<<"Done!"<<"\n"<<"\n";
                pContext->complete = true;
                
                pContext->condQueue.notify_all();
                init_terrace->write_summary_generation();
            }

        }
        omp_unset_lock(&pContext->thread_mutex);
    }

}

void split_threads(Terrace* init_terrace,
                    vector<string> &list_taxa_to_insert, 
                    int taxon_to_insert, 
                    bool use_dynamic_order,
                    string &taxon_name,
                    vector<int> &ids1,
                    vector<int> &ids2)
{
    
    Task* new_task = new Task;
    
    vector<string> *_list_taxa_to_insert = new vector<string>;
    *_list_taxa_to_insert = list_taxa_to_insert;

    // path up to now
    new_task->path_up_to_now = new int*[init_terrace->path_size];
    
    for(int i = 0; i<init_terrace->path_size; i++){
        new_task->path_up_to_now[i] = new int[2];
    }

    for(int i = 0; i < taxon_to_insert; i++){
        new_task->path_up_to_now[i][0] = init_terrace->path_up_to_now[i][0];
        new_task->path_up_to_now[i][1] = init_terrace->path_up_to_now[i][1];
    }

    vector<int> * _ids1 = new vector<int>;
    vector<int> * _ids2 = new vector<int>;

    *_ids1 = ids1;
    *_ids2 = ids2;

    new_task->list_taxa_to_insert = _list_taxa_to_insert;
    new_task->_taxon_name = taxon_name;
    new_task->path_size = init_terrace->path_size;
    new_task->_taxon_to_insert = taxon_to_insert;
    new_task->ids1 = _ids1;
    new_task->ids2 = _ids2;
    new_task->use_dynamic_order = use_dynamic_order;
    
    new_task->taskFunction = &generateTerrace_wrapper;
    
    submitTask(new_task);

}

void run_terrace_check(Terrace *terrace,Params &params){

    int i, count = 0;
    vector<MTree*> trees_on, trees_off;

    ifstream in;
    in.open(params.terrace_query_set);
    while (!in.eof()) {
        count++;
        MTree *tree = new MTree();
        tree->readTree(in, params.is_rooted);
        if(terrace->check_two_trees(tree)){
            cout<<"Checking query tree "<<count<<"..."<<"ON the stand."<<"\n";
            trees_on.push_back(tree);
        }else{

            cout<<"Checking query tree "<<count<<"..."<<"DOES NOT lie on the stand!"<<"\n";
            trees_off.push_back(tree);
        }
        delete tree;

        char ch;
        (in) >> ch;
        if (in.eof()) break;
        in.unget();
    }
    in.close();

    if(!trees_on.empty() and !trees_off.empty()){
        string out_file_on = params.out_prefix;
        out_file_on += ".on_stand";

        ofstream out_on;
        out_on.exceptions(ios::failbit | ios::badbit);
        out_on.open(out_file_on);

        for(i=0; i<trees_on.size(); i++){
            trees_on[i]->printTree(out_on,WT_BR_SCALE | WT_NEWLINE);
        }

        out_on.close();
    }

    if(!trees_off.empty()){
        string out_file_off = params.out_prefix;
        out_file_off += ".off_stand";

        ofstream out_off;
        out_off.exceptions(ios::failbit | ios::badbit);
        out_off.open(out_file_off);

        for(i=0; i<trees_off.size(); i++){
            trees_off[i]->printTree(out_off,WT_BR_SCALE | WT_NEWLINE);
        }

        out_off.close();
    }

    cout<<"\n"<<"----------------------------------------------------------------------------"<<"\n";
    cout<<"Done."<<"\n";
    cout<<"Checked "<<count<<" trees:"<<"\n";
    if(trees_on.size() == count){
        cout<<"- all trees belong to the same stand"<<"\n";
    }else if(trees_off.size() == count){
        cout<<"- none of the trees belongs to considered stand"<<"\n";
    }else{
        cout<<" - "<<trees_on.size()<<" trees belong to considered stand"<<"\n";
        cout<<" - "<<trees_off.size()<<" trees do not belong to considered stand"<<"\n";
    }
    cout<<"----------------------------------------------------------------------------"<<"\n"<<"\n";

}

void read_tree_set(const char *infile, bool &is_rooted, vector<TerraceTree*> &subtrees){

    ifstream in;
    in.open(infile);
    while (!in.eof()) {
        TerraceTree *tree = new TerraceTree();
        tree->readTree(in, is_rooted);
        subtrees.push_back(tree);

        char ch;
        (in) >> ch;
        if (in.eof()) break;
        in.unget();
    }
    in.close();
}


// ----------- Parallelization functions --------------------------------------

void submitTask(Task* task){

    {    
        std::lock_guard<mutex> guard(pContext->mutexQueue);
        taskQueue.push(task);
        pContext->taskCount++;
    }
    pContext->condQueue.notify_one();

}

void executeTask(Terrace* init_terrace, vector<Terrace*> *part_tree_pairs, Task* task){
    
    init_terrace->initial_split_done = true;

    init_terrace->terrace_trees_num = 0;
    init_terrace->intermediated_trees_num = 0;
    init_terrace->dead_ends_num = 0;

    init_terrace->cur_terrace_trees = 0;
    init_terrace->cur_intermediate_trees = 0;

    init_terrace->artificial_thread_num = 1;

    for(int i = 0; i < task->_taxon_to_insert; i++){
        init_terrace->path_up_to_now[i][0] = task->path_up_to_now[i][0];
        init_terrace->path_up_to_now[i][1] = task->path_up_to_now[i][1];
        
    }

    for(int i = 0; i<init_terrace->path_size; i++){
        delete[] task->path_up_to_now[i];
    }

    delete[] task->path_up_to_now;

    int tti = task->_taxon_to_insert;

    for(int i = init_terrace->initial_split_done_index; i< tti; i++){
        TerraceNode* node1 = (TerraceNode*)init_terrace->findNodeID(init_terrace->path_up_to_now[i][0]);
        TerraceNode* node2 = (TerraceNode*)init_terrace->findNodeID(init_terrace->path_up_to_now[i][1]);
        init_terrace->extendNewTaxon(task->list_taxa_to_insert->at(i),node1, node2, *part_tree_pairs, true);
    }
    
    bool thread_call = true;
    
    task->taskFunction(init_terrace,
                        *part_tree_pairs, 
                        *task->list_taxa_to_insert,
                        task->_taxon_to_insert, 
                        thread_call,
                        init_terrace->real_thread_num,
                        task->_taxon_name,
                        task->use_dynamic_order,
                        task->ids1,
                        task->ids2);
    
    
    for(int i=tti-1; i >=init_terrace->initial_split_done_index; i--){
        init_terrace->remove_one_taxon(task->list_taxa_to_insert->at(i), *part_tree_pairs);
    }

    delete task->list_taxa_to_insert;

    task->ids1->clear();
    task->ids2->clear();
    delete task->ids1;
    delete task->ids2;
    
    delete task;
}


void* startThread(Terrace* init_terrace, vector<Terrace*> *part_tree_pairs, int &thread_num){ 
    
    pContext->working[thread_num] = false;

    while(1){
        
        Task* task;
        
        std::unique_lock<mutex> lk(pContext->mutexQueue);
        pContext->condQueue.wait(lk, []{return (pContext->taskCount > 0 || pContext->complete);});
        
        if(pContext->complete){
            lk.unlock();
            pContext->condQueue.notify_all();
            break;

        }else {
            
            pContext->working[thread_num] = true;

            task = taskQueue.front();
            taskQueue.pop();
            pContext->taskCount--;
            //std::cout << "Task count = " << pContext->taskCount << std::endl;

            lk.unlock();
            executeTask(init_terrace, part_tree_pairs, task);
            
        }
    }

    return 0;
}


// 3770145