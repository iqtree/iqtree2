/*
 *  presenceabsencematrix.cpp
 *  Class of presence-absence matrices (gene pr-ab)
 *  Created on: Sep 10, 2020
 *      Author: Olga
 */

#include "presenceabsencematrix.hpp"
#include "alignment/superalignment.h"
#include "tree/node.h"
#include "tree/mtree.h"
#include "utils/timeutil.h"

PresenceAbsenceMatrix::PresenceAbsenceMatrix(){
    taxa_num = 0;
    part_num = 0;
    init();
};
PresenceAbsenceMatrix::~PresenceAbsenceMatrix(){};

void PresenceAbsenceMatrix::get_from_alignment(Params &params){
    int i, j;
    
    // Input is an alignment with partition info
    SuperAlignment* alignment = new SuperAlignment(params);
    IntVector pattern_0_1;
    
    taxa_num = alignment->getNSeq();
    part_num = alignment->getNSite();
    
    for(i=0; i<taxa_num; i++){
        pattern_0_1.clear();
        for(j=0; j<part_num; j++){
            if(alignment->taxa_index[i][j]==-1){
                pattern_0_1.push_back(0);
            }else{
                pattern_0_1.push_back(1);
            }
        }
        pr_ab_matrix.push_back(pattern_0_1);
        taxa_names.push_back(alignment->getSeqName(i));
    }

    init();
    
    delete alignment;
}

void PresenceAbsenceMatrix::read_pr_ab_matrix(const char *infile){
    ifstream in;
    //cout<<"\n"<<"-----------------------------------------------------"<<"\n";
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        in.exceptions(ios::badbit);
        read_pr_ab_matrix(in);
        in.close();
    } catch (const char* str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }
}

void PresenceAbsenceMatrix::read_pr_ab_matrix(istream &in){
    string str_rest,name;
    
    if(!(in>>taxa_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    if(!(in>>part_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    
    int i=0,j=0;
    for(i=0; i<taxa_num; i++){
        if(!(in>>name)) throw "Each line should start with a taxon name!";
        //if(name == "0" or name == "1") throw "Each line should start with a taxon name! 0 and 1 are not allowed as taxon names.";
        //cout<<"Taxon "<<i<<": "<<name<<"\n";
        taxa_names.push_back(name);
        vector<int> vec(part_num, -1);
        for(j=0; j<part_num; j++){
            if(!(in >> vec[j])) throw "Could not read a matrix entry! For each species make sure there are as many entries as the number of partitions specified in the first line of the file. Moreover, presence-absence matrix should only contain 0, 1!";
            if(vec[j] < 0) throw "Error: A negative entry! Presence-absence matrix should only contain 0, 1!";
            if(vec[j] > 1) throw "Error: The entry is greater than 1! Presence-absence matrix should only contain 0, 1!";
        }
        pr_ab_matrix.push_back(vec);
    }
    
    init();
};

void PresenceAbsenceMatrix::print_pr_ab_matrix(ostream &out){
    int sum = 0, i, j;
    IntVector part_sum;
    part_sum.resize(part_num,0);
    if(!Params::getInstance().print_pr_ab_matrix){
        out<<"Presence-absence matrix:"<<"\n";
    }
    out<<taxa_num<<" "<<part_num<<"\n";
    for(i=0; i<taxa_num; i++){
        sum=0;
        out<<taxa_names[i];
        for(j=0; j<part_num; j++){
            out<<" "<<pr_ab_matrix[i][j];
            sum+=pr_ab_matrix[i][j];
            part_sum[j]+=pr_ab_matrix[i][j];
        }
        if(!Params::getInstance().print_pr_ab_matrix){
            out<<" | "<<sum;
        }
        out<<"\n";
    }
    if(!Params::getInstance().print_pr_ab_matrix){
        out<<"--------------------"<<"\n";
        out<<"Partition coverage: ";
        for(j=0; j<part_num; j++){
            out<<part_sum[j]<<" ";
        }
        out<<"\n\n";
    }
};


int PresenceAbsenceMatrix::findTaxonID(string taxon_name){
    int id;
    for(id=0; id<taxa_names.size(); id++){
        if(taxa_names[id]==taxon_name){
            //cout<<"GET_TAXON_ID: taxa_name = "<<taxon_name<<" | taxa_id = "<<id<<"\n";
            return id;
        }
    }
    return -1;
};

void PresenceAbsenceMatrix::init(){
    flag_reorderAccordingToTree = false;
    missing_percent = 0.0;
}

void PresenceAbsenceMatrix::getPartTaxa(int part, MTree *tree, MTree *part_tree, NodeVector &part_taxa){
    
    NodeVector taxa_nodes;
    tree->getTaxa(taxa_nodes);
    
    //part_taxa.resize(taxa_num);
    
    //if(!flag_reorderAccordingToTree){
    //    reorderAccordingToTree(taxa_nodes);
    //}
    
    Node *node;
    int taxon_matrix_id;
    for(NodeVector::iterator it=taxa_nodes.begin(); it<taxa_nodes.end(); it++){
        //cout<<(*it)->name<<" id = "<<(*it)->id<<"\n";
        taxon_matrix_id = findTaxonID((*it)->name);
        //cout<<"TAXON_MATRIX_ID:"<<taxon_matrix_id<<"\n";
        if(pr_ab_matrix[taxon_matrix_id][part] == 1){
            node = part_tree->findLeafName((*it)->name);
            //cout<<"PREPARING PART_TAXA: part = "<<part<<"|leaf_id = "<<(*it)->id<<"|leaf_name = "<<(*it)->name<<"\n";
            assert(node && "ERROR: The leaf is not found on partition tree!");
            //part_taxa[taxon_matrix_id]=node;
            part_taxa.push_back(node);
        }
    }
    
    bool check = false;
    if(check){
        //cout<<"GET_taxa_nodes_info...."<<"\n";
        cout<<"Presence-absence info:"<<"\n";
        for(int i=0; i<taxa_num; i++){
            cout<<pr_ab_matrix[i][part]<<" ";
        }
        cout<<"\n";
        cout<<"Taxon names info:"<<"\n";
        for(int i=0; i<taxa_num; i++){
            cout<<taxa_names[i]<<" ";
        }
        cout<<"\n";
        
        cout<<"Partition taxa info:"<<"\n";
        for(int i=0; i<taxa_num; i++){
            if(part_taxa[i]){
                cout<<part_taxa[i]->name<<"("<<part_taxa[i]->id<<") ";
            }else{
                cout<<"NA ";
            }
        }
        cout<<"\n";
    }
}

void PresenceAbsenceMatrix::reorderAccordingToTree(NodeVector taxa_nodes){

    // WARNING: when adding new taxa, this function is not helpful, because the ids of new taxa (at the current setting, as of 06.10.20) are larger than the number of taxa (id of new taxon is set to the number of nodes, which is then increased by 1, when a taxon is added)
    
    //cout<<"BEFORE reordering according to the tree:"<<"\n";
    //print_pr_ab_matrix();
    
    int i=0;
    int id=-1;
    
    vector<IntVector> aux_matrix;
    vector<string> aux_names;
    
    aux_matrix.resize(taxa_num);
    aux_names.resize(taxa_num);
    
    for(i=0; i<taxa_nodes.size();i++){
        id=findTaxonID(taxa_nodes[i]->name);
        //cout<<taxa_nodes[i]->name<<" id = "<<id<<"\n";
        aux_matrix[taxa_nodes[i]->id]=pr_ab_matrix[id];
        aux_names[taxa_nodes[i]->id]=taxa_names[id];
    }
    
    pr_ab_matrix.clear();
    taxa_names.clear();
    
    for(i=0; i<taxa_num; i++){
        pr_ab_matrix.push_back(aux_matrix[i]);
        taxa_names.push_back(aux_names[i]);
    }
    
    //cout<<"AFTER reordering according to the tree:"<<"\n";
    //print_pr_ab_matrix();
}

vector<IntVector> getSubMatrix(vector<IntVector> pr_ab_complete, vector<string> taxa_names, MTree* tree){
    
    //int id;
    vector<IntVector> sub_matrix;
    NodeVector taxa_nodes;
    NodeVector::iterator it;
    string taxon_name;
    tree->getTaxa(taxa_nodes);
    for(it=taxa_nodes.begin(); it<taxa_nodes.end(); it++){
        taxon_name=(*it)->name;
        //id=getTaxonID_in_pr_ab_m(taxon_name);
    }
    
    return sub_matrix;
};

void PresenceAbsenceMatrix::getSubPrAbMatrix(vector<string> taxa_names_subset, PresenceAbsenceMatrix *submatrix, IntVector *parts){
 
    vector<string> not_found_taxon_names;
    bool found = false;
    
    int i,j,h;
    for(i=0; i<taxa_names_subset.size(); i++){
        //cout<<i<<": subset_taxon = "<<taxa_names_subset[i]<<"\n";
        found = false;
        for(j=0; j<taxa_names.size(); j++){
            //cout<<j<<": taxon_name = "<<taxa_names[j]<<"\n";
            if(taxa_names_subset[i]==taxa_names[j]){
                //cout<<"MATCH"<<"\n";
                if(parts){
                    IntVector taxon_coverage;
                    for(IntVector::iterator k=parts->begin(); k<parts->end(); k++){
                        h = (*k);
                        taxon_coverage.push_back(pr_ab_matrix[j][h]);
                    }
                    submatrix->pr_ab_matrix.push_back(taxon_coverage);
                }else{
                    submatrix->pr_ab_matrix.push_back(pr_ab_matrix[j]);
                }
                submatrix->taxa_names.push_back(taxa_names[j]);
                found = true;
                break;
            }
        }
        if(!found){
            cout<<"Taxon "<<taxa_names_subset[i]<<" is not found in the presence-absence matrix..."<<"\n";
            not_found_taxon_names.push_back(taxa_names_subset[i]);
        }
    }
    
    bool print_info = false;
    if(not_found_taxon_names.size()<taxa_names_subset.size()){
        submatrix->taxa_num = submatrix->taxa_names.size();
        submatrix->part_num = submatrix->pr_ab_matrix[0].size();
        if(print_info){
            cout<<"INFO: original matrix."<<"\n";
            print_pr_ab_matrix();
            cout<<"\n"<<"INFO: a submatrix for "<<submatrix->taxa_num<<" taxa was extracted."<<"\n";
            submatrix->print_pr_ab_matrix();
        }
    }
}

void PresenceAbsenceMatrix::getSubPrAbMatrix(NodeVector taxon_nodes, PresenceAbsenceMatrix *submatrix, IntVector *parts){

    vector<string> taxon_names;
    for(NodeVector::iterator it = taxon_nodes.begin(); it<taxon_nodes.end(); it++){
        taxon_names.push_back((*it)->name);
    }
    
    getSubPrAbMatrix(taxon_names,submatrix,parts);
}

void PresenceAbsenceMatrix::extend_by_new_taxa(string taxon_name, IntVector pr_ab_pattern){
    
    taxa_names.push_back(taxon_name);
    pr_ab_matrix.push_back(pr_ab_pattern);
    
    taxa_num+=1;
    
    flag_reorderAccordingToTree = false;
}

void PresenceAbsenceMatrix::remove_taxon(string taxon_name){
    
    //cout<<"REMOVING taxon "<<taxon_name<<" from matrix."<<"\n";
    int id = findTaxonID(taxon_name);
    
    /*cout<<"BEFORE:"<<"\n";
    print_pr_ab_matrix();
    cout<<"MATRIX   dim: "<<pr_ab_matrix.size()<<"x"<<pr_ab_matrix[0].size()<<"\n";
    cout<<"TAXAname dim:"<<taxa_names.size()<<"\n";*/
    
    pr_ab_matrix.erase(pr_ab_matrix.begin()+id);
    taxa_names.erase(taxa_names.begin()+id);
    
    taxa_num-=1;
    
    /*cout<<"AFTER:"<<"\n";
    print_pr_ab_matrix();
    cout<<"MATRIX   dim: "<<pr_ab_matrix.size()<<"x"<<pr_ab_matrix[0].size()<<"\n";
    cout<<"TAXAname dim:"<<taxa_names.size()<<"\n";*/
    
    flag_reorderAccordingToTree = false;
    
}

int PresenceAbsenceMatrix::getINFO_init_tree_taxon_order(vector<string> &taxa_names_sub, vector<string> &list_taxa_to_insert,const int m){
    
    //cout<<"\n"<<"=================================================="<<"\n"<<"INITIAL tree and Taxon Order:"<<"\n"<<"\n";
    int i,j,k;
    
    IntVector part_cov, taxon_cov, uniq_taxa;
    part_cov.resize(part_num,0);
    taxon_cov.resize(taxa_num,0);
    uniq_taxa.resize(part_num,0);
    
    for(i=0; i<part_num;i++){
        for(j=0; j<taxa_num; j++){
            taxon_cov[j]+=pr_ab_matrix[j][i];
        }
    }
    
    for(j=0; j<taxa_num; j++){
        if(taxon_cov[j]==0){
            cout<<"ERROR: taxon "<<j<<" ("<<taxa_names[j]<<")"<<" is not covered by any partition (row sum = 0). Remove it from the dataset.\n";
            assert(taxon_cov[j]!=0);
        }
    }
    
    
    for(i=0; i<part_num; i++){
        for(j=0; j<taxa_num; j++){
            if(taxon_cov[j]!=1){ // taxa covered only by one partition
                part_cov[i]+=pr_ab_matrix[j][i]; // per partition you only sum up taxa, which are present in more than one partition
            }
            if(taxon_cov[j]==1 && pr_ab_matrix[j][i]==1){
                uniq_taxa[i]+=1;
            }
        }
        if(part_cov[i]+uniq_taxa[i] == taxa_num){
            cout<<"\n"<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
            cout<<"INPUT INFO:"<<"\n";
            cout<<"---------------------------------------------------------"<<"\n";
            cout<<"Number of taxa: "<<taxa_num<<"\n";
            cout<<"Number of partitions: "<<part_num<<"\n";
            cout<<"Number of special taxa (row sum = 1): "<<uniq_taxa_num<<"\n";
            percent_missing();
            cout<<"% of missing entries in supermatrix: "<<missing_percent<<"\n";
            cout<<"---------------------------------------------------------"<<"\n";
            cout<<"\n"<<"INFO: "<<"\n";
            cout<<"At least one partition (partition "<<i+1<<") covers all taxa."<<"\n"<<"There are only trivial terraces (contain just 1 tree) for this dataset. Great!"<<"\n"<<"\n";
            cout<<"---------------------------------------------------------"<<"\n";
            cout<<"SUMMARY:"<<"\n";
            cout<<"Number of trees on terrace: "<<1<<"\n";
            cout<<"Number of intermediated trees visited: "<<0<<"\n";
            cout<<"Number of dead ends encountered: "<<0<<"\n";
            cout<<"---------------------------------------------------------"<<"\n";
            cout<<"Total wall-clock time used: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")"<<"\n";
            cout<<"Total CPU time used: "
            << getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")" << "\n";
            cout<<"\n";
            
            exit(0);
        }
    }
    
    uniq_taxa_num = 0;
    for(const auto& u: uniq_taxa){
        uniq_taxa_num += u;
    }
    
    // order partitions by taxon coverage | addition: if partitions have the same number of taxa, order them by the overlap with partitions with larger number of taxa
    vector<int> ordered_partitions;
    
    // BEGIN: ORDER TYPE 1 ====================================================================
    int id;
    ordered_partitions.push_back(0);
    
    bool inserted;
    for(i=1; i<part_num; i++){
        inserted = FALSE;
        for(j=0; j<ordered_partitions.size(); j++){ // TODO: maybe implement something more efficient
            id=ordered_partitions[j];
            //if(part_cov[i]+uniq_taxa[i]>=part_cov[id]+uniq_taxa[id]){
              if(part_cov[i]>=part_cov[id]){
                ordered_partitions.insert(ordered_partitions.begin()+j, i);
                inserted = TRUE;
                break;
            }
        }
        if(!inserted){
            ordered_partitions.push_back(i);
        }
    }
     
    // END: ORDER TYPE 1 ====================================================================
    
    // BEGIN: ORDER TYPE 2 ====================================================================
    /*vector<IntVector> new_order;
    IntVector aux;
    i=0;
    for(const auto& it: part_cov){
        aux.clear();
        aux.push_back(i);
        aux.push_back(it);
        new_order.push_back(aux);
        i++;
    }
    // ADDITIONAL ORDER: refining order within partitions with the same coverage according to their overlap with partitions of higher order
    reordering(new_order, new_order.begin(), new_order.end(), part_cov);
    for(const auto& v: new_order){
        ordered_partitions.push_back(v[0]);
    }
    */
    // END: ORDER TYPE 2 ====================================================================
    
    //cout<<"Partition with maximum coverage based on non-specific taxa is "<<ordered_partitions[0]<<" with "<<part_cov[ordered_partitions[0]]<<" taxa."<<"\n";
    
    // Vector for part id stores its order
    vector<int> part_order;
    part_order.resize(part_num);
    for(i=0; i<ordered_partitions.size(); i++){
        part_order[ordered_partitions[i]]=i;
    }
    
    /*for(i=0; i<part_num; i++){
        cout<<"Partition "<<i<<": ordered as "<<part_order[i]<<" with coverage of "<<part_cov[i]<<"\n";
    }*/
    
    int part_max = ordered_partitions[0];
    
    vector<IntVector> ordered_taxa_ids;
    ordered_taxa_ids.resize(part_num+1);
    
    // collecting taxon names from partition with the largest subtree excluding unique taxa
    
    for(j=0; j<taxa_num; j++){
        //if(pr_ab_matrix[j][part_max]==1 && taxon_cov[j]>1){ // collecting only non-unique taxa
        if(pr_ab_matrix[j][part_max]==1){ // collecting all taxa from the chosen partition also unique
            taxa_names_sub.push_back(taxa_names[j]);
        } else {
            // ordering other taxa by the number of partitions they are covered by
            for(i=1;i<part_num+1;i++){
                if(taxon_cov[j]==i){
                    ordered_taxa_ids[i].push_back(j);
                    break;
                }
            }
        }
    }
    
    int part_init = -1;
    if(m==0){
        cout<<"Partition "<<part_max+1<<" is chosen for the initial tree."<<"\n";
        part_init = part_max;
        // count the number of unique taxa to insert on part_init
        uniq_taxa_to_insert_num=uniq_taxa_num;
        for(j=0; j<taxa_num; j++){
            if(pr_ab_matrix[j][part_init]==1 and taxon_cov[j]==1){
                uniq_taxa_to_insert_num--;
            }
        }
    } else {
        cout<<"The initial tree will be created by removing from the input tree "<<m<<" leaves."<<"\n";
        cout<<"Note, that this procedure does not guarantee generating all trees from a stand! It is only meant to investigate, if at least some trees from a stand can be generated."<<"\n";
    }
    
    IntVector ordered_ids;
    
    for(i=1; i<part_num; i++){
        if(ordered_taxa_ids[i].size()>1){
            //cout<<"\n"<<"ANALYSING COVERAGE SIZE "<<i<<"\n";
            vector<IntVector> cov_within_cat;
            cov_within_cat.resize(ordered_taxa_ids[i].size());
            ordered_ids.clear();
            
            for(j=0; j<ordered_taxa_ids[i].size(); j++){
                for(k=0;k<part_num;k++){
                    if(pr_ab_matrix[ordered_taxa_ids[i][j]][k]==1){
                        cov_within_cat[j].push_back(part_order[k]);
                    }
                }
                sort(cov_within_cat[j].begin(), cov_within_cat[j].end());
                
                //cout<<"TAXON_COVERAGE_INFO: taxon_"<<ordered_taxa_ids[i][j]<<"|"<<taxa_names[ordered_taxa_ids[i][j]];
                //for(k=0; k<cov_within_cat[j].size(); k++){
                //    cout<<" "<<cov_within_cat[j][k];
                //}
                //cout<<"\n";
            }
            
            // REORDER taxa within each coverage category
            orderTaxaByCoverage(ordered_taxa_ids[i],cov_within_cat,ordered_ids);
            
            ordered_taxa_ids[i].clear();
            ordered_taxa_ids[i] = ordered_ids;
            
            //cout<<"REORDERED taxa:"<<"\n";
            //for(k=0; k<ordered_taxa_ids[i].size(); k++){
            //    cout<<k+1<<"|taxon "<<ordered_taxa_ids[i][k]<<"|"<<taxa_names[ordered_taxa_ids[i][k]]<<"\n";
            //}
        }
    }
    
    const int must_insert = taxa_num - m;
    int inserted_num = taxa_names_sub.size();

    // ===================================================================
    // TAXON ORDER 1: first by coverage, then by partition order
    // ===================================================================
    if(m==0){
        for(i=ordered_taxa_ids.size()-1; i>-1; i--){
            if(!ordered_taxa_ids[i].empty()){
                //cout<<"Taxa with "<<i<<" gene coverage: "<<"\n";
                for(j=0; j<ordered_taxa_ids[i].size(); j++){
                    list_taxa_to_insert.push_back(taxa_names[ordered_taxa_ids[i][j]]);
                }
                //cout<<"\n";
            }
        }
    }else{
        for(i=ordered_taxa_ids.size()-1; i>-1; i--){
            if(!ordered_taxa_ids[i].empty()){
                //cout<<"Taxa with "<<i<<" gene coverage: "<<"\n";
                for(j=0; j<ordered_taxa_ids[i].size(); j++){
                    if(must_insert == inserted_num){
                        //cout<<" "<<j<<": "<<ordered_taxa_by_coverage[i][j]<<"\n";
                        list_taxa_to_insert.push_back(taxa_names[ordered_taxa_ids[i][j]]);
                    }else{
                        taxa_names_sub.push_back(taxa_names[ordered_taxa_ids[i][j]]);
                        inserted_num++;
                    }
                }
                //cout<<"\n";
            }
        }
    }
   
    /*cout<<"\n"<<"FINAL taxon order:"<<"\n";
    for(j=0; j<list_taxa_to_insert.size(); j++){
        cout<<j<<": "<<list_taxa_to_insert[j]<<"\n";
    }*/
    
    return part_init;
}

void PresenceAbsenceMatrix::orderTaxaByCoverage(vector<int> &taxon_ids, vector<IntVector> &coverage_info, IntVector &ordered_taxa){
    
    int i,j,k;
    
    IntVector ordered_taxa_ids;
    ordered_taxa_ids.push_back(0);
    bool inserted_check;
    
    /*for(i=0; i<coverage_info.size(); i++){
        for(j=0; j<coverage_info[i].size(); j++){
            cout<<" "<<coverage_info[i][j];
        }
        cout<<"\n";
    }*/
    
    for(i=1; i<taxon_ids.size(); i++){
        //cout<<"CHECKING taxon "<<i<<" with id "<<taxon_ids[i]<<"\n";
        inserted_check = FALSE;
        
        for(j=0; j<ordered_taxa_ids.size(); j++){
            //cout<<"comparing with taxon "<<j<<" with id "<<taxon_ids[ordered_taxa_ids[j]]<<"\n";
            int ident=0;
            for(k=0; k<coverage_info[0].size(); k++){
                //cout<<"partition_order_info:"<<coverage_info[i][k]<<" vs. "<<coverage_info[ordered_taxa_ids[j]][k]<<"\n";
                if(coverage_info[i][k]<coverage_info[ordered_taxa_ids[j]][k]){
                    // insert
                    ordered_taxa_ids.insert(ordered_taxa_ids.begin()+j,i);
                    inserted_check = TRUE;
                    //cout<<"inserted"<<"\n";
                    break;
                } else if(coverage_info[i][k]>coverage_info[ordered_taxa_ids[j]][k]){
                    // do not insert, check next taxon
                    //cout<<"breaking, check next taxon.."<<"\n";
                    break;
                } else {
                    ident+=1;
                    // TODO: maybe do something more elaborate with taxa having identical pattern of partition coverage, such that you do not have to check next taxon with n taxa with the same partition order..
                }
            }
            
            if(inserted_check){
                break;
            }
        }
        
        // if the taxon is last one, insert to the end of the vector
        if(!inserted_check){
            ordered_taxa_ids.push_back(i);
        }
    }
    
    for(i=0; i<taxon_ids.size(); i++){
        ordered_taxa.push_back(taxon_ids[ordered_taxa_ids[i]]);
    }
    
}

void PresenceAbsenceMatrix::percent_missing(){
    
    missing_percent = 0.0;
    int i,j,a = 0;
    for(i=0; i<taxa_num; i++){
        for(j=0; j<part_num; j++){
            a+= pr_ab_matrix[i][j];
        }
    }

    double total = taxa_num*part_num;
    missing_percent=(total-a)/total*100.0;
}

void PresenceAbsenceMatrix::orderPartByOverlap(IntVector &ordered_partitions,IntVector &part_cov){
    
    int i;
    int cov = part_cov[ordered_partitions[0]];
    
    IntVector group;
    int upper_lim=0;
    
    for(i=0;i<part_num; i++){
        if(part_cov[ordered_partitions[i]]==cov){
            group.push_back(i);
        }else{
            // refine within the group of the same cov
            if(group.size()>1){
                // at the moment, there is no higher partition, check overlap between partitions
                if(upper_lim==0){
                    if(group.size()>2){
                        cout<<"BEFORE: ordered_partitions: ";
                        for(int l=0; l<part_num; l++){
                            cout<<" "<<ordered_partitions[l];
                        }
                        cout<<"\n";
                        
                        orderPartByOverlap_within(upper_lim, ordered_partitions, group, part_cov);
                        
                        cout<<"AFTER: ordered_partitions: ";
                        for(int l=0; l<part_num; l++){
                            cout<<" "<<ordered_partitions[l];
                        }
                        cout<<"\n";
                    }
                }else{
                    // there are partitions before considered group
                    orderPartByOverlap_preceding(upper_lim, ordered_partitions, group, part_cov);
                }
            }
            
            // update value for next group
            if(i!=part_num-1){
                cov=part_cov[ordered_partitions[i+1]];
                upper_lim=i+1;
            }
        }
    }
}

void PresenceAbsenceMatrix::orderPartByOverlap_within(int upper_lim, IntVector &ordered_partitions, IntVector &group, IntVector &part_cov){

    int j, k, val;
    bool inserted;
    
    vector<IntVector> overlap;
    IntVector aux, new_order;
    aux.resize(group.size(),0);
    for(j=0; j<group.size(); j++){
        overlap.push_back(aux);
    }
    
    for(j=0; j<group.size()-1; j++){
        for(k=j+1; k<group.size(); k++){
            
            val=0;
            for(int t=0; t<taxa_num; t++){
                if(pr_ab_matrix[t][group[j]]==pr_ab_matrix[t][group[k]] && pr_ab_matrix[t][group[j]]==1){
                    val=val+1;
                    if(val==part_cov[group[j]] or val==part_cov[group[k]]){
                        break;
                    }
                }
            }
            overlap[j][k]=val;
            overlap[k][j]=val;
        }
    }
    cout<<"------------"<<"\n";
    for(j=0; j<group.size(); j++){
        cout<<group[j];
        for(k=0; k<group.size(); k++){
            cout<<" "<<overlap[j][k];
        }
        cout<<"\n";
    }
    cout<<"------------"<<"\n";
    for(j=0; j<group.size(); j++){
        sort(overlap[j].begin(), overlap[j].end(),greater<int>());
    }
    cout<<"------------"<<"\n";
    for(j=0; j<group.size(); j++){
        cout<<group[j];
        for(k=0; k<group.size(); k++){
            cout<<" "<<overlap[j][k];
        }
        cout<<"\n";
    }
    cout<<"------------"<<"\n";
    
    cout<<"NEW ORDER"<<"\n";
    new_order.push_back(0);
    inserted=false;
    for(j=1; j<group.size(); j++){
        inserted=false;
        for(k=0; k<new_order.size();k++){
            if(overlap[j][0]>=overlap[new_order[k]][0]){
                new_order.insert(new_order.begin()+k,j);
                inserted=true;
                break;
            }
        }
        if(!inserted){
            new_order.push_back(j);
        }
    }
    
    int id=upper_lim;
    for(j=0;j<new_order.size();j++){
        ordered_partitions[id]=group[new_order[j]];
        id++;
    }
    
    if(group.size()>2){
        IntVector group_sub;
        upper_lim=upper_lim+1;
        for(j=1;j<new_order.size();j++){
            group_sub.push_back(new_order[j]);
        }
        orderPartByOverlap_preceding(upper_lim, ordered_partitions, group_sub, part_cov);
    }
    
};


void PresenceAbsenceMatrix::orderPartByOverlap_preceding(int upper_lim, IntVector &ordered_partitions, IntVector &group, IntVector &part_cov){
    
    int i,j,k, t,val;
    vector<IntVector> overlap;
    IntVector aux;
    aux.resize(upper_lim,0);
    for(i=0;i<group.size();i++){
        overlap.push_back(aux);
        for(j=0;j<upper_lim; j++){
            val=0;
            for(t=0; t<taxa_num; t++){
                if(pr_ab_matrix[t][ordered_partitions[j]]==pr_ab_matrix[t][group[i]] && pr_ab_matrix[t][group[i]]==1){
                    val=val+1;
                    if(val==part_cov[group[i]] or val==part_cov[ordered_partitions[j]]){
                        break;
                    }
                }
            }
            overlap[i][j]=val;
        }
        //sort(overlap[i].begin(), overlap[i].end(),greater<int>()); // there is no need to sort, because you need to follow the order of preceding partitions
    }
    
    cout<<"------------"<<"\n";
    cout<<"preceding"<<"\n";
    cout<<"------------"<<"\n";
    for(j=0; j<group.size(); j++){
        cout<<group[j];
        for(k=0; k<upper_lim; k++){
            cout<<" "<<overlap[j][k];
        }
        cout<<"\n";
    }
    cout<<"------------"<<"\n";
    
    IntVector new_order;
    new_order.push_back(0);
    bool inserted=false;
    for(j=1; j<group.size(); j++){
        inserted=false;
        for(k=0; k<new_order.size();k++){
            if(overlap[j][0]>=overlap[new_order[k]][0]){
                new_order.insert(new_order.begin()+k,j);
                inserted=true;
                break;
            }
        }
        if(!inserted){
            new_order.push_back(j);
        }
    }
    
    //sort(arr, arr + n, greater<int>());

    int id=upper_lim;
    for(j=0;j<new_order.size();j++){
        ordered_partitions[id]=group[new_order[j]];
        id++;
    }
}


void PresenceAbsenceMatrix::getPartOverlap(int part, IntVector &group, IntVector &part_cov, IntVector &overlap){
    
    int i,id,t,val;
    
    for(i=0;i<group.size();i++){
        val=0;
        id=group[i];
        for(t=0; t<taxa_num; t++){
            if(pr_ab_matrix[t][id]==1){
                if(pr_ab_matrix[t][part]==pr_ab_matrix[t][id]){
                    val++;
                    if(val==part_cov[part] or val==part_cov[id]){
                        break;
                    }
                }
                
            }
        }
        overlap.push_back(val);
    }
};

int PresenceAbsenceMatrix::get2PartOverlap(int part_1, int part_2, int max_overlap){
    int val=0, t;
    
    for(t=0; t<taxa_num; t++){
        if(pr_ab_matrix[t][part_2]==1){
            if(pr_ab_matrix[t][part_1]==pr_ab_matrix[t][part_2]){
                val++;
                if(val==max_overlap or val==max_overlap){
                    break;
                }
            }
        }
    }
    
    return val;
}

void PresenceAbsenceMatrix::reordering(vector<IntVector> &new_order, vector<IntVector>::iterator it_b, vector<IntVector>::iterator it_e, IntVector &part_cov,int i){
  
    //cout<<"======================================="<<"\n";
    //puts("REORDERING:");
    vector<IntVector>::iterator it;
    /*for(it=it_b;it!=it_e; it++){
        cout<<(*it)[0]<<"|"<<(*it)[1]<<"\n";
    }
    cout<<"======================================="<<"\n";*/
    vector<IntVector>::iterator it_end = it_e;
    int j=1,k=0;
    sort(it_b, it_e, [j](vector<int> a,  vector<int> b){return a[j] > b[j];});
    it_b = adjacent_find(it_b, it_e,[j](vector<int> a,  vector<int> b){return a[j] == b[j];});
    
    
    while(it_b != it_end){
        //cout<<"--------------------------------------"<<"\n";
        int val = (*it_b)[j];
        it_e = it_b;
        while((*it_e)[j]==val && it_e!=it_end-1){
            it_e++;
        }
        if(it_e==it_end-1 && (*it_e)[j]==val){
            it_e++;
        }
        
        int dist=distance(new_order.begin(),it_b);
        //cout<<"dist:="<<dist<<" | i:="<<i<<"\n";
        //cout<<"--------------------------------------"<<"\n";
        if(i<dist or (i>dist && i<=part_num)){
            /*cout<<"Case I/II: decision on overlap with preceding/following partitions"<<"\n";
            cout<<"overlap based on i:="<<i<<"|partition_"<<new_order[i][0]<<"\n";
            cout<<"--------------------------------------"<<"\n";
            puts("Considered group:");
            for(it=it_b;it!=it_e; it++){
                cout<<(*it)[0]<<"|"<<(*it)[1]<<"\n";
            }*/
            for(it=it_b; it!=it_e; it++){
                (*it)[1]=get2PartOverlap(new_order[i][0],(*it)[0],min(part_cov[new_order[i][0]],part_cov[(*it)[0]]));
            }
            /*puts("New info:");
            for(it=it_b;it!=it_e; it++){
                cout<<(*it)[0]<<"|"<<(*it)[1]<<"\n";
            }*/
            reordering(new_order, it_b, it_e, part_cov,i+1);
        }else if(i==dist){
            if(distance(it_b,it_e)>2){
                //cout<<"Case III: decision on within group overlap"<<"\n";
                // get_max_overlap_within_group
                
                
                
                
            }else if(it_e!=new_order.end()){
                //cout<<"Case II: decision on overlap with following partition"<<"\n";
                // switch to following partitions
                k=distance(new_order.begin(),it_e);
                /*cout<<"overlap based on k:="<<k<<"|partition_"<<new_order[k][0]<<"\n";
                cout<<"--------------------------------------"<<"\n";
                puts("Considered group:");
                for(it=it_b;it!=it_e; it++){
                    cout<<(*it)[0]<<"|"<<(*it)[1]<<"\n";
                }*/
                for(it=it_b; it!=it_e; it++){
                    (*it)[1]=get2PartOverlap(new_order[k][0],(*it)[0],min(part_cov[new_order[k][0]],part_cov[(*it)[0]]));
                }
                /*puts("New info:");
                for(it=it_b;it!=it_e; it++){
                    cout<<(*it)[0]<<"|"<<(*it)[1]<<"\n";
                }*/
                reordering(new_order, it_b, it_e, part_cov,j+1);
            }
        }
        if(distance(it_e,it_end)>1){
            it_b = adjacent_find(it_e, it_end,[j](vector<int> a,  vector<int> b){return a[j] == b[j];});;
        }else{
            it_b = it_end;
        }
    }
        
};


void PresenceAbsenceMatrix::getPartOverlapComplete(){
    assert(part_num!=0 && "ERROR: assertion part_num!=0 failed in getPartOverlapComplete()..");
    assert(taxa_num!=0 && "ERROR: assertion taxa_num!=0 failed in getPartOverlapComplete()..");
    
    overlap_matrix.resize(part_num);
    for(auto &v: overlap_matrix){
        v.resize(part_num,0);
    }
    int i{0}, j{0};
    for(; i<part_num-1; i++){
        for(j=i+1; j<part_num; j++){
            overlap_matrix[i][j]=get2PartOverlap(i, j, part_num);
            overlap_matrix[j][i]=overlap_matrix[i][j];
        }
    }
};

void PresenceAbsenceMatrix::print_overlap_matrix(ostream &out){
    assert(overlap_matrix.size()==part_num && "ERROR: assertion overlap_matrix.size()==part_num failed in print_overlap_matrix()..");
    out<<"\n"<<"Printing matrix of partition overlap:"<<"\n";
    int i{0};
    /*cout<<" ";
    for(;i<part_num;i++){
        cout<<" "<<i;
    }
    cout<<"\n";*/
    i=0;
    for(const auto& v: overlap_matrix){
        IntVector sum;
        sum.resize(part_num,0);
        out<<"Part "<<i+1<<":";
        for(const auto& u: v){
            sum[i]+=u;
            out<<" "<<u;
        }
        out<<" | "<<sum[i]<<"\n";
        i++;
    }
    out<<"\n";
};

void PresenceAbsenceMatrix::get_from_subtrees(vector<TerraceTree*> subtrees){
    
    part_num = subtrees.size();
    
    string taxon_name;
    IntVector vec;
    vec.resize(part_num,0);

    int p=0, id;
    
    for(const auto &t: subtrees){
        for(const auto &l: t->leafNodes){
            taxon_name = l.first;
            id=findTaxonID(taxon_name);
            if(id==-1){
                taxa_num++;
                taxa_names.push_back(taxon_name);
                pr_ab_matrix.push_back(vec);
                pr_ab_matrix[taxa_num-1][p]=1;
            }else{
                pr_ab_matrix[id][p]=1;
            }
        }
        p++;
    }
}
