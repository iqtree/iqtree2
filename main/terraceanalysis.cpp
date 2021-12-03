/*
 *  terraceanalysis.cpp
 *  Interface to work with terraces
 *  Created on: Sep 10, 2020
 *      Author: Olga
 */

#include "terraceanalysis.h"
#include "terrace/terracenode.hpp"
#include "terrace/terracetree.hpp"
#include "terrace/terrace.hpp"
#include "terrace/presenceabsencematrix.hpp"
#include "tree/mtreeset.h"
#include "alignment/superalignment.h"
#include "utils/timeutil.h"


void runterraceanalysis(Params &params){
	
    params.startCPUTime = getCPUTime();
    params.start_real_time = getRealTime();

    cout<<"---------------------------------------------------------\n";
    cout<<"\nBegin analysis with Gentrius...\n\n";
    cout<<"---------------------------------------------------------\n";
    
    if(params.matrix_order){
        cout<<"Creating presence-absence matrix..\n\n";
        PresenceAbsenceMatrix *matrix = new PresenceAbsenceMatrix();
        if(params.partition_file && params.aln_file){
            matrix->get_from_alignment(params);
        }else{
            assert(params.pr_ab_matrix && "ERROR: no input presence-absence matrix!");
            matrix->read_pr_ab_matrix(params.pr_ab_matrix);
        }
        //matrix->print_pr_ab_matrix();
        
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
    Terrace *terrace;
    
    if(params.partition_file && params.user_file){
        PresenceAbsenceMatrix *matrix = new PresenceAbsenceMatrix();
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
    
    /*terrace->out_file = params.out_prefix;
    terrace->out_file += ".terrace_info";
    
    ofstream out_terrace_info;
    out_terrace_info.exceptions(ios::failbit | ios::badbit);
    out_terrace_info.open(terrace->out_file);
    
    out_terrace_info.close();*/
    
    
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
        
        //Terrace* test = new Terrace(terrace->induced_trees);
        //test->rm_leaves = 0;
        //run_generate_trees(test, params,test->rm_leaves);
        run_generate_trees(terrace, params,terrace->rm_leaves);
    }

};

void run_generate_trees(Terrace *terrace, Params &params,const int m){
    
    assert(m>=0 && "ERROR: required number m of leaves to be removed from the input tree is incorrect! m must be >=0. Exiting...");
    if(m>=0){
        terrace->rm_leaves = m;
    }
    // GET INFO FOR INITIAL TREE and TAXA TO INSERT
    vector<string> taxa_names_sub;      // taxa to appear on initial tree
    vector<string> list_taxa_to_insert; // list of taxa to be inserted.
    int part_init = terrace->matrix->getINFO_init_tree_taxon_order(taxa_names_sub,list_taxa_to_insert,terrace->rm_leaves);
    

    // INITIAL TERRACE
    // - to be used for generating trees
    // - agile master tree is the initial tree
    // - induced partition trees are common subtrees with induced partition trees of considered terrace
    
    // Creating a subterrace: submatrix and an initial tree
    PresenceAbsenceMatrix *submatrix = new PresenceAbsenceMatrix();
    terrace->matrix->getSubPrAbMatrix(taxa_names_sub, submatrix);

    TerraceTree tree_init;
    if(terrace->root){
        tree_init.copyTree_byTaxonNames(terrace,taxa_names_sub);
    }else{
        assert(part_init!=-1);
        tree_init.copyTree(terrace->induced_trees[part_init]);
    }
        //tree_init.drawTree(cout, WT_BR_SCALE | WT_TAXON_ID | WT_NEWLINE);
    
    Terrace *init_terrace = new Terrace(tree_init, submatrix);
    init_terrace->rm_leaves = m;
    init_terrace->master_terrace = terrace;

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
    
    
    init_terrace->trees_out_lim = params.terrace_print_lim;
    init_terrace->linkTrees(true, false); // branch_back_map, taxon_back_map; in this case you only want to map branches

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
    cout<<"Number of taxa to be inserted: "<<list_taxa_to_insert.size()<<"\n";
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
    
    bool use_dynamic_taxon_order = true;
    if(use_dynamic_taxon_order){
        // TAXON ORDER: Based on the number of allowed branches.
        vector<string> ordered_taxa_to_insert;
        ordered_taxa_to_insert = list_taxa_to_insert;
        // The taxon order is dynamicly adapted using information about allowed branches
        init_terrace->generateTerraceTrees(terrace, part_tree_pairs, list_taxa_to_insert, 0, &ordered_taxa_to_insert);
    }else{
        // Taxon order is fixed. Taxa inserted in order they appear in list_taxa_to_insert
        init_terrace->generateTerraceTrees(terrace, part_tree_pairs, list_taxa_to_insert, 0, nullptr);
    }
    cout<<"\n"<<"---------------------------------------------------------"<<"\n";
    cout<<"\n"<<"Done!"<<"\n"<<"\n";

    init_terrace->write_summary_generation();
};



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

