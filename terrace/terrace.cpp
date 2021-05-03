/*
 *  terrace.cpp
 *  Terrace class is used to generate trees that have idetical set of induced partition trees.
 *  Created on: Sep 14, 2020
 *      Author: Olga
 */

#include "terrace.hpp"
#include "main/terraceanalysis.h"
#include "terracenode.hpp"
#include "tree/mtreeset.h"
#include "utils/timeutil.h"

Terrace::Terrace(){};
Terrace::~Terrace(){

    for(vector<TerraceTree*>::reverse_iterator it=induced_trees.rbegin(); it<induced_trees.rend();it++){
        delete (*it);
    }
    induced_trees.clear();
    
    delete matrix;
    master_terrace = nullptr;
};

void Terrace::init(){
    
    taxa_num = 0;
    part_num = 0;
    intermediated_trees_num = 0;
    terrace_trees_num = 0;
    dead_ends_num = 0;
    terrace_out = true;
    trees_out_lim = 0;
    
    terrace_max_trees = 1000000;
    intermediate_max_trees = 10000000;
    seconds_max = 604800; // the default value is 7 days
    
}

Terrace::Terrace(const char *infile_tree, bool is_rooted,const char *infile_matrix){
    
    init();
    
    readTree(infile_tree,is_rooted);
    
    matrix = new PresenceAbsenceMatrix();
    matrix->read_pr_ab_matrix(infile_matrix);
    
    taxa_num = matrix->pr_ab_matrix.size();
    part_num = matrix->pr_ab_matrix[0].size();

    //NodeVector taxa_nodes;
    //getTaxa(taxa_nodes);
    //cout<<"BEFORE reodering"<<"\n";
    //matrix->print_pr_ab_matrix();
    get_part_trees();
    //matrix->reorderAccordingToTree(taxa_nodes);
    //printInfo();
};

Terrace::Terrace(TerraceTree tree, PresenceAbsenceMatrix *m){
    
    init();

    if(tree.leafNum>2){
        copyTree(&tree);
    }else{
        string taxa_set = "";
        vector<uint32_t> check_int;
        check_int.resize(2,1);
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        copyTree(&tree,taxa_set);
    }
    
    matrix = m;
    
    taxa_num = leafNum;
    part_num = matrix->pr_ab_matrix[0].size();
    
    //matrix->print_pr_ab_matrix();
    //NodeVector taxa_nodes;
    //getTaxa(taxa_nodes);
    //cout<<"BEFORE reodering"<<"\n";
    //matrix->print_pr_ab_matrix();
    //matrix->reorderAccordingToTree(taxa_nodes);
    //cout<<"AFTER reodering"<<"\n";
    //matrix->print_pr_ab_matrix();
    
    get_part_trees();
    //printInfo();
    
    fillLeafNodes();
}

Terrace::Terrace(TerraceTree tree, PresenceAbsenceMatrix *m, vector<TerraceTree*> input_induced_trees){
 
    init();
    
    if(tree.leafNum>2){
        copyTree(&tree);
    }else{
        string taxa_set = "";
        vector<uint32_t> check_int;
        check_int.resize(2,1);
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        copyTree(&tree,taxa_set);
    }
    
    matrix = m;
    
    taxa_num = leafNum;
    part_num = input_induced_trees.size();
    
    //printTree(cout);
    //matrix->print_pr_ab_matrix();
    
    set_part_trees(input_induced_trees);
    //printInfo();

    fillLeafNodes();
}

void Terrace::get_part_trees(){
    
    assert(matrix!=nullptr && "ERROR: Presence-absence matrix is absent! I cannot get induced partition trees..");
    //cout<<"Preparing induced partition trees..."<<"\n";
    
    int id,k;
    string taxa_set = "";
    string taxon_name;
    NodeVector taxa_nodes;
    
    NodeVector::iterator it2;
    vector<uint32_t> check_int;
    check_int.resize(taxa_num);
    
    getTaxa(taxa_nodes);
    int leafNum_part_tree = 0;
    for(int part=0; part<part_num; part++){
        leafNum_part_tree = 0;
        TerraceTree* induced_tree = new TerraceTree();
        for(it2=taxa_nodes.begin();it2!=taxa_nodes.end();it2++){
            taxon_name=(*it2)->name;
            id=matrix->findTaxonID(taxon_name);
            assert(id != -1 && "Not all of the taxa appear in the pr_ab_matrix!");
            check_int[(*it2)->id] = matrix->pr_ab_matrix[id][part];
            //cout<<"Taxon["<<(*it2)->name<<"|"<<(*it2)->id<<"|"<<id<<"] in partition "<<part<<" is "<<check_int[(*it2)->id]<<"||MATRIX_INFO:"<<matrix->taxa_names[id]<<"-"<<matrix->pr_ab_matrix[id][part]<<"\n";
            if(check_int[(*it2)->id]==1){
                leafNum_part_tree+=1;
            }
        }
        taxa_set.clear();
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        
        if(leafNum_part_tree>1){
            induced_tree->copyTree(this,taxa_set);
            induced_tree->leafNum = leafNum_part_tree;
        } else if(leafNum_part_tree==1){
            for(k=0; k<taxa_num; k++){
                if(check_int[taxa_nodes[k]->id]==1){
                    taxon_name=taxa_nodes[k]->name;
                    break;
                }
            }
            
            induced_tree->root = induced_tree->newNode(0, taxon_name.c_str());
            induced_tree->leafNum = 1;
            induced_tree->nodeNum = 1;
        }
        
        if(induced_tree->leafNum>0){
            induced_tree->fillLeafNodes();
        }
        
        induced_trees.push_back(induced_tree);
    }
    
}

void Terrace::set_part_trees(vector<TerraceTree*> input_induced_trees){
    
    for(int i=0; i<input_induced_trees.size(); i++){
        induced_trees.push_back(input_induced_trees[i]);
    }
    
    // TODO: you need a check, that your input_induced_trees indeed correspond to presence-absence pattern of pr_ab_matrix
}

void Terrace::unset_part_trees(){
    induced_trees.clear();
    part_num=0;
    cleanAllLinkINFO();
};

void Terrace::printInfo(ostream &out){

    //out<<"\n"<<"=================================================="<<"\n";
    out<<"Printing terrace information:"<<"\n"<<"\n";
    out<<"Terrace representative tree:"<<"\n";
    print_terrace_tree(true,out);
    out<<"\n";
    
    /*if(matrix){
        out<<"Presence-absence matrix:"<<"\n";
        matrix->print_pr_ab_matrix(out);
    }*/
    
    if(induced_trees.size()>0){
        int i=0;
        out<<"Induced partition trees:"<<"\n";
        for(vector<TerraceTree*>::iterator it = induced_trees.begin(); it < induced_trees.end(); it++){
            i++;
            out<<"Part["<<i<<"]: ";
            (*it)->print_terrace_tree(false,out);
        }
    }
    //out<<"\n"<<"=================================================="<<"\n"<<"\n";
}

void Terrace::linkTrees(bool back_branch_map, bool back_taxon_map){
    NodeVector part_taxa;
    
    for(int part=0; part<part_num; part++){
        if(induced_trees[part]->leafNum>2){
            //part_taxa.clear();
            //matrix->getPartTaxa(part, this, induced_trees[part], part_taxa);
            linkTree(part, part_taxa, back_branch_map, back_taxon_map);
        }
    }
    
}

void Terrace::linkTree(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node, TerraceNode *dad){
    
    // SEHR WICHTIG! WARNING: do not mix mapping from the parent tree and upper level induced partition trees, because empty branches and empty taxa will be messed up on partition trees.
    
    // INFO/CHECK: make sure that your branches have proper unique ids, further code relies on that.
    
    if (!node) {
        if (!root->isLeaf())
            node = (TerraceNode*) root;
        else
            node = (TerraceNode*) root->neighbors[0]->node;
        ASSERT(node);
        if (node->isLeaf()) // two-taxa parent tree
            dad = (TerraceNode*)node->neighbors[0]->node;
    }
    TerraceNeighbor *nei = NULL;
    TerraceNeighbor *dad_nei = NULL;
    if (dad) {
        nei = (TerraceNeighbor*)node->findNeighbor(dad);
        dad_nei = (TerraceNeighbor*)dad->findNeighbor(node);
        if (nei->link_neighbors.empty()) nei->link_neighbors.resize(part_num);
        if (dad_nei->link_neighbors.empty()) dad_nei->link_neighbors.resize(part_num);
        nei->link_neighbors[part] = NULL;
        dad_nei->link_neighbors[part] = NULL;
    }
    if (node->isLeaf()) {
        ASSERT(dad);
        
        TerraceNode *node_part = nullptr;
        //==================================================
        /*for(NodeVector::iterator it=part_taxa.begin(); it<part_taxa.end(); it++){
            if((*it)){
                if(node->name == (*it)->name){
                    node_part = (TerraceNode*)(*it);
                    break;
                }
            }
        }
        
        TerraceNode *node_part_check = node_part;
        node_part=nullptr;*/
        
        if(induced_trees[part]->leafNodes.find(node->name)!=induced_trees[part]->leafNodes.end()){
            node_part=(TerraceNode*)induced_trees[part]->leafNodes[node->name];
            assert(node->name==node_part->name);
        }
        
        //assert(node_part_check==node_part && "ERROR: linkTree : getPartTaxa != hash_map");*/
        //==================================================
        
        if (node_part) {
            TerraceNode *dad_part = (TerraceNode*)node_part->neighbors[0]->node;
            TerraceNeighbor *dad_part_nei = (TerraceNeighbor*)dad_part->findNeighbor(node_part);
            ASSERT(node_part->isLeaf());
            nei->link_neighbors[part] = (TerraceNeighbor*) node_part->neighbors[0];
            dad_nei->link_neighbors[part] = dad_part_nei;
            
            if(back_branch_map){
                // Back map from induced partition tree onto parent tree
                // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
                ((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors.push_back(nei);
                dad_part_nei->link_neighbors.push_back(dad_nei);
            }
            
            if(back_taxon_map){
                // Back map from induced partition tree onto parent tree. Here used for low to top induced partition trees backward map
                // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
                ((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors_lowtop_back.push_back(nei);
                dad_part_nei->link_neighbors_lowtop_back.push_back(dad_nei);
            }
            
        } else {
            
            // the empty branches are needed for a map from parent to partition
            dad->empty_br_dad_nei.push_back(dad_nei);
            dad->empty_br_node_nei.push_back(nei);
        }
        return;
    }
    
    //cout<<"MASTER NODE:"<<node->id<<"\n";
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        //cout<<"Node->id = "<<(*it)->node->id<<" | Dad->id = "<<node->id<<"\n";
        linkTree(part, part_taxa, back_branch_map, back_taxon_map, (TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
    
    if (!dad) {
        // Check, if the final dad has empty_branches and map them to branch, which is available.
        // Note, that if there are some empty branches/taxa, there will be exactly one branch available for mapping.
        if(!node->empty_br_node_nei.empty()){
            //cout<<"I got here! Partition: "<<part<<"| Node->id:"<<node->id<<"\n";
            FOR_NEIGHBOR_DECLARE(node, NULL, it) {
                if(((TerraceNeighbor*)(*it))->link_neighbors[part]){
                    TerraceNeighbor* node_nei_part = (TerraceNeighbor*)((TerraceNeighbor*)(*it))->link_neighbors[part];
                    TerraceNeighbor* dad_nei_part = (TerraceNeighbor*)((TerraceNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part];
                    int i;
                    for(i=0; i<node->empty_br_dad_nei.size(); i++){
                        ((TerraceNeighbor*)node->empty_br_node_nei[i])->link_neighbors[part] = node_nei_part;
                        ((TerraceNeighbor*)node->empty_br_dad_nei[i])->link_neighbors[part] = dad_nei_part;
            
                        if(back_branch_map){
                            node_nei_part->link_neighbors.push_back(node->empty_br_node_nei[i]);
                            dad_nei_part->link_neighbors.push_back(node->empty_br_dad_nei[i]);
                        }
                        
                        if(back_taxon_map){
                            node_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_node_nei[i]);
                            dad_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_dad_nei[i]);
                        }
                    }
                        
                    node->empty_br_node_nei.clear();
                    node->empty_br_dad_nei.clear();
                    
                    return;
                }
            }
        }
        return;
    }
    linkBranch(part, nei, dad_nei, back_branch_map, back_taxon_map);
}

void Terrace::linkBranch(int part, TerraceNeighbor *nei, TerraceNeighbor *dad_nei, bool back_branch_map, bool back_taxon_map) {
    
    TerraceNode *node = (TerraceNode*)dad_nei->node;
    TerraceNode *dad = (TerraceNode*)nei->node;
    nei->link_neighbors[part] = NULL;
    dad_nei->link_neighbors[part] = NULL;
    vector<TerraceNeighbor*> part_vec;
    vector<TerraceNeighbor*> child_part_vec;
    
    //cout<<"Linking branch ("<<node->id<<","<<dad->id<<");"<<"\n";
    
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        if (((TerraceNeighbor*)*it)->link_neighbors[part]) {
            part_vec.push_back((TerraceNeighbor*)(((TerraceNeighbor*)*it)->link_neighbors[part]));
            child_part_vec.push_back((TerraceNeighbor*)(((TerraceNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]));
            
            // INFO: "folding" of branches, which can be folded to branch and map to "nodeA-nodeB" or "nodeB-nodeA", children neighbors can be different. Therefore, in some cases the assertion won't be fulfilled, but this is not an error and, thus, the assertion is commented out.
            //ASSERT(child_part_vec.back()->node == child_part_vec.front()->node || child_part_vec.back()->id == child_part_vec.front()->id);
        }
    }
    
    if(child_part_vec.size()==2){
        if(child_part_vec[0]->node != child_part_vec[1]->node){
            //cout<<"My weird case I..."<<"\n";
            //cout<<"child_part_vec[0]="<<child_part_vec[0]->node->id<<"\n";
            //cout<<"child_part_vec[1]="<<child_part_vec[1]->node->id<<"\n";
            //cout<<"part_vec[0]="<<part_vec[0]->node->id<<"\n";
            //cout<<"part_vec[1]="<<part_vec[1]->node->id<<"\n";
            TerraceNeighbor * nei_aux;
            if(part_vec[0]->node == part_vec[1]->node){
                nei_aux = part_vec[0];
                part_vec[0]=child_part_vec[0];
                child_part_vec[0]=nei_aux;
                nei_aux=part_vec[1];
                part_vec[1]=child_part_vec[1];
                child_part_vec[1]=nei_aux;
            }else{
                for(int i=0; i<part_vec.size(); i++){
                    for(int j=0; j<child_part_vec.size(); j++){
                        if(part_vec[i]->node == child_part_vec[j]->node and part_vec[1-i]->node != child_part_vec[1-j]->node){
                            nei_aux = part_vec[i];
                            part_vec[i] = child_part_vec[1-j];
                            child_part_vec[1-j] = nei_aux;
                            break;
                        }
                    }
                }
            }
            //cout<<"child_part_vec[0]="<<child_part_vec[0]->node->id<<"\n";
            //cout<<"child_part_vec[1]="<<child_part_vec[1]->node->id<<"\n";
            //cout<<"part_vec[0]="<<part_vec[0]->node->id<<"\n";
            //cout<<"part_vec[1]="<<part_vec[1]->node->id<<"\n";
        }
    }

    
    
    // QUESTION: node or dad empty branch/taxa?
    // QUESTION: if you start mapping the dad, you might encounter situations, when a dad does not know yet, that there are some empty branches.
    // So at the moment you need only empty branches of node and the dad will be mapped at the final state (i.e. final node - the start of the tree traversal)
    
    if (part_vec.empty()){
        int i=0;
        //cout<<"CASE 0: CHILDREN DO NOT HAVE IMAGE"<<"\n"<<"\n";
        if(!node->empty_br_dad_nei.empty()){
            for(i=0; i<node->empty_br_dad_nei.size(); i++){
                dad->empty_br_dad_nei.push_back(node->empty_br_dad_nei[i]);
                dad->empty_br_node_nei.push_back(node->empty_br_node_nei[i]);
            }
            node->empty_br_node_nei.clear();
            node->empty_br_dad_nei.clear();
        }
        
        // since below node the subtrees are empty, add also current branch to an empty list of the dad
        dad->empty_br_dad_nei.push_back(dad_nei);
        dad->empty_br_node_nei.push_back(nei);

        return;
    }
    
    int i=0;
    if (part_vec.size() == 1) {
        
        //cout<<"CASE 1: ONE CHILD HAS IMAGE"<<"\n"<<"\n";
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = part_vec[0];
        
        if(back_branch_map){
            child_part_vec[0]->link_neighbors.push_back(nei);
            part_vec[0]->link_neighbors.push_back(dad_nei);
        }
        
        if(back_taxon_map){
            child_part_vec[0]->link_neighbors_lowtop_back.push_back(nei);
            part_vec[0]->link_neighbors_lowtop_back.push_back(dad_nei);
        }
        
        // Get maps for empty branches
        if(node->empty_br_node_nei.size()>0){
            //cout<<"INFO CHECK (one child image, get maps for the empty branches): ("<<node->id<<","<<dad->id<<")"<<"\n";
            for(i=0; i<node->empty_br_node_nei.size(); i++){
                
                ((TerraceNeighbor*)node->empty_br_node_nei[i])->link_neighbors[part]=child_part_vec[0];
                ((TerraceNeighbor*)node->empty_br_dad_nei[i])->link_neighbors[part] = part_vec[0];
                
                if(back_branch_map){
                    child_part_vec[0]->link_neighbors.push_back(node->empty_br_node_nei[i]);
                    part_vec[0]->link_neighbors.push_back(node->empty_br_dad_nei[i]);
                    //cout<<i<<":"<<"("<<node->empty_br_node_nei[i]->node->id<<","<<node->empty_br_dad_nei[i]->node->id<<") -> ("<<child_part_vec[0]->node->id<<","<<part_vec[0]->node->id<<")"<<"\n";
                }
                
                if(back_taxon_map){
                    child_part_vec[0]->link_neighbors_lowtop_back.push_back(node->empty_br_node_nei[i]);
                    part_vec[0]->link_neighbors_lowtop_back.push_back(node->empty_br_dad_nei[i]);
                }
                
            }
            node->empty_br_node_nei.clear();
            node->empty_br_dad_nei.clear();
        }
        return;
    }
    
    // what's this below? -> When a subtree on the other side is empty. Do nothing for it.
    // In your case, you have have to map empty branches to (child_part_vec[0] ---- part_vec[1]) or (child_part_vec[1] ---- part_vec[0]). Does it matter, which nodes are mapped where? I suspect, that no. Aber there is always aber
    if (part_vec[0] == child_part_vec[1]) {
        // ping-pong, out of sub-tree
        ASSERT(part_vec[1] == child_part_vec[0]);
        
        // Get maps for empty branches
        // the question is if your dad already has the details about all empty branches/taxa below ????? do not touch yet dad's empty branches, they should be considered at a later stage.
        // I think, in this situation you only need to map current branch, because node cannot have empty branches and taxa from its two children
        // again, i'm not sure if the direction of neighbors is important or not...
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = child_part_vec[1];

        
        if(back_branch_map){
            child_part_vec[0]->link_neighbors.push_back(nei);
            child_part_vec[1]->link_neighbors.push_back(dad_nei);
        }
        if(back_taxon_map){
            child_part_vec[0]->link_neighbors_lowtop_back.push_back(nei);
            child_part_vec[1]->link_neighbors_lowtop_back.push_back(dad_nei);
        }
        return;
    }
    
    // WARNING/INFO: FOR THE UPDATE_MAP: CASE 2 (two identical images for child branches), might also happen, when
    // child_part_vec[0]=child_part_vec[1]
    // part_vec[0]=part_vec[1]
    // this is due to the fact, that all branches have a link_neighbors and the direction in which the branch "folds" is not certain
    
    if(child_part_vec[0]==child_part_vec[1]){

        assert(part_vec[0]==part_vec[1]);
        
        nei->link_neighbors[part] = child_part_vec[0];
        dad_nei->link_neighbors[part] = part_vec[0];
        
        if(back_branch_map){
            child_part_vec[0]->link_neighbors.push_back(nei);
            part_vec[0]->link_neighbors.push_back(dad_nei);
        }
        if(back_taxon_map){
            child_part_vec[0]->link_neighbors_lowtop_back.push_back(nei);
            part_vec[0]->link_neighbors_lowtop_back.push_back(dad_nei);
        }
        return;
    }
    
    
    //cout<<"CASE 2: TWO CHILDREN HAVE IMAGEs AND THEY ARE DIFFERENT"<<"\n"<<"\n";
    TerraceNode *node_part = (TerraceNode*) child_part_vec[0]->node;
    TerraceNode *dad_part = nullptr;
    FOR_NEIGHBOR(node_part, NULL, it) {
        bool appear = false;
        for (vector<TerraceNeighbor*>::iterator it2 = part_vec.begin(); it2 != part_vec.end(); it2++){
            if ((*it2) == (*it)) {
                appear = true; break;
            }
        }
        if (!appear) {
            ASSERT(!dad_part);
            dad_part = (TerraceNode*)(*it)->node;
        }
    }
    
    nei->link_neighbors[part] = (TerraceNeighbor*)node_part->findNeighbor(dad_part);
    dad_nei->link_neighbors[part] = (TerraceNeighbor*)dad_part->findNeighbor(node_part);
    
    if(back_branch_map){
        ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.push_back(nei);
        ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors.push_back(dad_nei);
    }
    
    if(back_taxon_map){
        ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors_lowtop_back.push_back(nei);
        ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors_lowtop_back.push_back(dad_nei);
    }
}

void Terrace::update_map(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node, TerraceNode *dad){
    
    if (node->isLeaf()){
        if(!dad){
            dad = (TerraceNode*)node->neighbors[0]->node;
            if(!dad->isLeaf()){
                node = dad;
                dad = nullptr;
            }
        }
    }
    
    if(!dad){
        FOR_NEIGHBOR_DECLARE(node, dad, it) {
            //cout<<"Node->id = "<<(*it)->node->id<<" | Dad->id = "<<node->id<<"\n";
            update_map(part, part_taxa, back_branch_map, back_taxon_map, (TerraceNode*) (*it)->node, (TerraceNode*) node);
        }
        
        // Check, if the final dad has empty_branches and map them to branch, which is available.
        // Note, that if there are some empty branches/taxa, there will be exactly one branch available for mapping.
        if(!node->empty_br_node_nei.empty()){
            //cout<<"I got here! Partition: "<<part<<"| Node->id:"<<node->id<<"\n";
            FOR_NEIGHBOR_DECLARE(node, NULL, it) {
                if(((TerraceNeighbor*)(*it))->link_neighbors[part]){
                    TerraceNeighbor* node_nei_part = (TerraceNeighbor*)((TerraceNeighbor*)(*it))->link_neighbors[part];
                    TerraceNeighbor* dad_nei_part = (TerraceNeighbor*)((TerraceNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part];
                    int i;
                    for(i=0; i<node->empty_br_dad_nei.size(); i++){
                        ((TerraceNeighbor*)node->empty_br_node_nei[i])->link_neighbors[part] = node_nei_part;
                        ((TerraceNeighbor*)node->empty_br_dad_nei[i])->link_neighbors[part] = dad_nei_part;
                        
                        if(back_branch_map){
                            node_nei_part->link_neighbors.push_back(node->empty_br_node_nei[i]);
                            dad_nei_part->link_neighbors.push_back(node->empty_br_dad_nei[i]);
                        }
                        
                        if(back_taxon_map){
                            node_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_node_nei[i]);
                            dad_nei_part->link_neighbors_lowtop_back.push_back(node->empty_br_dad_nei[i]);
                        }
                    }
                    
                    node->empty_br_node_nei.clear();
                    node->empty_br_dad_nei.clear();
                    
                    return;
                }
            }
        }
        return;
        
    } else {
        TerraceNeighbor *nei = (TerraceNeighbor*)node->findNeighbor(dad);
        TerraceNeighbor *dad_nei = (TerraceNeighbor*)dad->findNeighbor(node);
    
        // if there is a map already, return
        if(nei->link_neighbors[part]){
            if(dad_nei->link_neighbors[part])
                return;
        } else {
    
            if (node->isLeaf()) {
                ASSERT(dad);
                TerraceNode *node_part = nullptr;
            
                /*for(NodeVector::iterator it=part_taxa.begin(); it<part_taxa.end(); it++){
                    if((*it)){
                        if(node->name == (*it)->name){
                            node_part = (TerraceNode*)(*it);
                            break;
                        }
                    }
                }
                TerraceNode *node_part_check = node_part;
                node_part=nullptr;*/
                
                if(induced_trees[part]->leafNodes.find(node->name)!=induced_trees[part]->leafNodes.end()){
                    node_part=(TerraceNode*)induced_trees[part]->leafNodes[node->name];
                    assert(node->name==node_part->name);
                }
                
                //assert(node_part_check == node_part && "ERROR: update_map : getPartTaxa != hash_map");
            
                if (node_part) {
                    TerraceNode *dad_part = (TerraceNode*)node_part->neighbors[0]->node;
                    TerraceNeighbor *dad_part_nei = (TerraceNeighbor*)dad_part->findNeighbor(node_part);
                    ASSERT(node_part->isLeaf());
                    nei->link_neighbors[part] = (TerraceNeighbor*) node_part->neighbors[0];
                    dad_nei->link_neighbors[part] = dad_part_nei;
                    
                    if(back_branch_map){
                        // Back map from induced partition tree onto parent tree
                        // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
                        ((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors.push_back(nei);
                        dad_part_nei->link_neighbors.push_back(dad_nei);
                    }
                    
                    if(back_taxon_map){
                        // Back map from induced partition tree onto parent tree (for LOW-TOP backmap)
                        // for partition trees link_neighbours contain all nei's from parent tree, that map on this branch
                        ((TerraceNeighbor*) node_part->neighbors[0])->link_neighbors_lowtop_back.push_back(nei);
                        dad_part_nei->link_neighbors_lowtop_back.push_back(dad_nei);
                    }
                    
                } else {
                    // the empty branches are needed for a map from parent to partition
                    dad->empty_br_dad_nei.push_back(dad_nei);
                    dad->empty_br_node_nei.push_back(nei);
                }
                return;
            }
        
            //cout<<"MASTER NODE:"<<node->id<<"\n";
            FOR_NEIGHBOR_DECLARE(node, dad, it) {
                //cout<<"Node->id = "<<(*it)->node->id<<" | Dad->id = "<<node->id<<"\n";
                update_map(part, part_taxa, back_branch_map, back_taxon_map, (TerraceNode*) (*it)->node, (TerraceNode*) node);
            }

            linkBranch(part, nei, dad_nei, back_branch_map, back_taxon_map);
        }
    }
    
}

void Terrace::printMapInfo(int partition){
    
    cout<<"Mapping info"<<"\n"<<"\n";
    NodeVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    int part = 0;
    drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    if(partition==-1){
        for (vector<TerraceTree*>::iterator it = induced_trees.begin(); it != induced_trees.end(); it++, part++) {
            cout << "\n" << "Subtree for partition " << part << "("<<(*it)->leafNum<<","<<(*it)->nodeNum<<","<<(*it)->branchNum<<")"<< "\n";
            // INFO: drawing of two-taxon tree fails.
            if((*it)->leafNum>2){
                (*it)->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE );
            }
            for (int i = 0; i < nodes1.size(); i++) {
                if(((TerraceNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors.size()>0){
                    Neighbor *nei1 = ((TerraceNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part];
                    Neighbor *nei2 = ((TerraceNeighbor*)nodes2[i]->findNeighbor(nodes1[i]))->link_neighbors[part];
                    cout << nodes1[i]->findNeighbor(nodes2[i])->id << ":";
                    if (nodes1[i]->isLeaf()) cout << nodes1[i]->name; else cout << nodes1[i]->id;
                    cout << ",";
                    if (nodes2[i]->isLeaf()) cout << nodes2[i]->name; else cout << nodes2[i]->id;
                    //cout <<"("<<nodes1[i]->findNeighbor(nodes2[i])->length<<")"<< " -> ";
                    cout<<" -> ";
                    if (nei2) {
                        cout << nei2->id << ":";
                        if (nei2->node->isLeaf())
                            cout << nei2->node->name;
                        else cout << nei2->node->id;
                    }
                    else cout << -1;
                    cout << ",";
                    if (nei1){
                        if (nei1->node->isLeaf())
                            cout << nei1->node->name;
                        else cout << nei1->node->id;
                        //cout <<"("<<nei1->length<<")";
                    }
                    else cout << -1;
                    cout << "\n";
                }
            }
        }
    } else {
        part = partition;
        cout << "\n" << "Subtree for partition " << part << "("<<induced_trees[part]->leafNum<<","<<induced_trees[part]->nodeNum<<","<<induced_trees[part]->branchNum<<")"<< "\n";
        // INFO: drawing of two-taxon tree fails.
        if(induced_trees[part]->leafNum>2){
            induced_trees[part]->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE );
        }
        for (int i = 0; i < nodes1.size(); i++) {
            if(((TerraceNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors.size()>0){
                Neighbor *nei1 = ((TerraceNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part];
                Neighbor *nei2 = ((TerraceNeighbor*)nodes2[i]->findNeighbor(nodes1[i]))->link_neighbors[part];
                cout << nodes1[i]->findNeighbor(nodes2[i])->id << ":";
                if (nodes1[i]->isLeaf()) cout << nodes1[i]->name; else cout << nodes1[i]->id;
                cout << ",";
                if (nodes2[i]->isLeaf()) cout << nodes2[i]->name; else cout << nodes2[i]->id;
                //cout <<"("<<nodes1[i]->findNeighbor(nodes2[i])->length<<")"<< " -> ";
                cout<<" -> ";
                if (nei2) {
                    cout << nei2->id << ":";
                    if (nei2->node->isLeaf())
                        cout << nei2->node->name;
                    else cout << nei2->node->id;
                }
                else cout << -1;
                cout << ",";
                if (nei1){
                    if (nei1->node->isLeaf())
                        cout << nei1->node->name;
                    else cout << nei1->node->id;
                    //cout <<"("<<nei1->length<<")";
                }
                else cout << -1;
                cout << "\n";
            }
        }
    }
    cout << "\n";
}

void Terrace::printBackMapInfo(){

    // INFO: Also for upper level induced partition trees their induced subtrees are the same induced partition trees as for the master tree
    cout<<"\n"<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<"\n";
    cout<<"\n"<<"BACKWARD mapping information:"<<"\n";
    cout<<"\n"<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<"\n";
    
    int i,j,k;
    NodeVector node_1, node_2;
    TerraceNeighbor *nei12, *nei21;
    
   for(i=0; i<part_num; i++){
        cout<<"\n"<<"---------------------------------------------"<<"\n";
        cout<<"\n"<<"Partition "<<i<<":"<< "("<<induced_trees[i]->leafNum<<","<<induced_trees[i]->nodeNum<<","<<induced_trees[i]->branchNum<<")"<<"\n";
        cout<<"\n"<<"---------------------------------------------"<<"\n";
        if(induced_trees[i]->leafNum>1){
            node_1.clear();
            node_2.clear();
            induced_trees[i]->getBranches(node_1, node_2);
            for(j=0; j<node_1.size();j++){
                nei12 = (TerraceNeighbor*) node_1[j]->findNeighbor(node_2[j]);
                nei21 = (TerraceNeighbor*) node_2[j]->findNeighbor(node_1[j]);
                cout<<"\n"<<"* branch "<<nei12->id<<": ";
                if(node_1[j]->isLeaf()){
                    cout<<node_1[j]->name;
                }
                cout<<"("<<node_1[j]->id<<")"<<",";
                if(node_2[j]->isLeaf()){
                    cout<<node_2[j]->name;
                }
                cout<<"("<<node_2[j]->id<<")"<<"\n";
                cout<<"+ link_neighbors:"<<"\n";
                if(!nei12->link_neighbors.empty()){
                    for(k=0; k<nei12->link_neighbors.size();k++){
                        cout<<" - "<<nei21->link_neighbors[k]->id<<":";
                        if(nei21->link_neighbors[k]->node->isLeaf()){
                            cout<<nei21->link_neighbors[k]->node->name;
                        }
                        cout<<"("<<nei21->link_neighbors[k]->node->id<<"),";
                        if(nei12->link_neighbors[k]->node->isLeaf()){
                            cout<<nei12->link_neighbors[k]->node->name;
                        }
                        cout<<"("<<nei12->link_neighbors[k]->node->id<<")"<<"\n";
                    }
                }

                cout<<"+ link_neighbors_lowtop_back:"<<"\n";
                if(!nei12->link_neighbors_lowtop_back.empty()){
                    for(k=0; k<nei12->link_neighbors_lowtop_back.size();k++){
                        cout<<" - "<<nei21->link_neighbors_lowtop_back[k]->id<<":";
                        if(nei21->link_neighbors_lowtop_back[k]->node->isLeaf()){
                            cout<<nei21->link_neighbors_lowtop_back[k]->node->name;
                        }
                        cout<<"("<<nei21->link_neighbors_lowtop_back[k]->node->id<<"),";
                        if(nei12->link_neighbors_lowtop_back[k]->node->isLeaf()){
                            cout<<nei12->link_neighbors_lowtop_back[k]->node->name;
                        }
                        cout<<"("<<nei12->link_neighbors_lowtop_back[k]->node->id<<")"<<"\n";
                    }
                }
            }
        }
    }
}

void Terrace::create_Top_Low_Part_Tree_Pairs(vector<Terrace*> &part_tree_pairs, Terrace *terrace){
    
    int i=0, j=0;
    NodeVector aux_taxon_nodes;
    vector<TerraceTree*> aux_induced_part_trees;
    IntVector parts;
    bool back_branch_map = false, back_taxon_map = true;
    
    for(i=0; i<terrace->part_num; i++){
        /*parts.clear();
        parts.push_back(i);
        aux_taxon_nodes.clear();
        
        // 1. get submatrix for taxa, which are common for the initial tree and for the top level induced partition tree. common taxa, are the leaves of the low level induced tree (assert: they should be with 1's!!!)
        PresenceAbsenceMatrix *aux_submatrix = new PresenceAbsenceMatrix();
        if(induced_trees[i]->leafNum>0){
            induced_trees[i]->getTaxa(aux_taxon_nodes);
            matrix->getSubPrAbMatrix(aux_taxon_nodes, aux_submatrix,&parts);
        }
        // 2. now update matrix with taxa, which occur only on the top level induced partition tree and not on the initial tree. Simply add 0, because it does not occur on the common subtree, i.e. low level induced partition tree
        aux_taxon_nodes.clear();
        terrace->induced_trees[i]->getTaxa(aux_taxon_nodes);
        
        parts.clear(); // auxiliary use of the variable, just need an integer vec with 0 value.
        parts.push_back(0);
        
        for(j=0; j<aux_taxon_nodes.size(); j++){
            if(induced_trees[i]->leafNum>0){
                if(aux_submatrix->findTaxonID(aux_taxon_nodes[j]->name) == -1){
                    aux_submatrix->pr_ab_matrix.push_back(parts);
                    aux_submatrix->taxa_names.push_back(aux_taxon_nodes[j]->name);
                    aux_submatrix->taxa_num +=1;
                }
            }else{
                aux_submatrix->pr_ab_matrix.push_back(parts);
                aux_submatrix->taxa_names.push_back(aux_taxon_nodes[j]->name);
                aux_submatrix->taxa_num +=1;
            }
        }
        aux_submatrix->part_num = 1;
        //aux_submatrix->print_pr_ab_matrix();
        //terrace->induced_trees[i]->printTree(cout,WT_BR_LEN_ROUNDING + WT_NEWLINE);
        //cout<<"\n";
         */
        
        aux_induced_part_trees.clear();
        aux_induced_part_trees.push_back(induced_trees[i]);
        // IF MATRIX is used
        //Terrace *aux_terrace = new Terrace(*(terrace->induced_trees[i]), aux_submatrix, aux_induced_part_trees);
        // IF HASHMAP is used
        Terrace *aux_terrace = new Terrace(*(terrace->induced_trees[i]), nullptr, aux_induced_part_trees);
        aux_terrace->linkTrees(back_branch_map, back_taxon_map);
        //aux_terrace->printMapInfo();
        //aux_terrace->printBackMapInfo();
        part_tree_pairs.push_back(aux_terrace);
        
        //aux_terrace->matrix->print_pr_ab_matrix();
    }
}

//void Terrace::getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, NodeVector *node1_vec, NodeVector *node2_vec, NodeVector *branch_end_1, NodeVector *branch_end_2){
void Terrace::getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, NodeVector *node1_vec, NodeVector *node2_vec){
    // INFO/CHECK: make sure that your branches have proper unique ids, below code relies on that.
    //cout<<"\n"<<"**********************************************"<<"\n"<<"\n";
    //cout<<"IN getALLOWEDbranches: "<<taxon_name<<"\n";
    //cout<<"\n"<<"**********************************************"<<"\n"<<"\n";

    //printMapInfo();
    //printBackMapInfo();
    
    int i, j, h, k;
    TerraceNode *node;
    TerraceNeighbor *nei, *dad_nei, *link_nei, *link_dad_nei;
    
    IntVector branch_ids, intersection;
    bool found = false;
    
    //NodeVector branch_end_1, branch_end_2;
    //Neighbor* aux_nei;
    /*if(branch_end_1->size()==0){
        this->getBranches(*branch_end_1, *branch_end_2);
    }*/
    
    for(i=0;i<aux_terrace.size();i++){
        //node = (TerraceNode*) aux_terrace[i]->findLeafName(taxon_name);
        //if(node){
            //cout<<"Taxon: "<<taxon_name<<"| FOUND in partition "<<i<<"\n";
        if(aux_terrace[i]->leafNodes.find(taxon_name)!=aux_terrace[i]->leafNodes.end()){
            node=(TerraceNode*) aux_terrace[i]->leafNodes[taxon_name];
            assert(node->name==taxon_name);
            //assert(node->isLeaf());
            
            if(induced_trees[i]->leafNum>2){
                nei = (TerraceNeighbor*) node->neighbors[0];
                dad_nei = (TerraceNeighbor*) nei->node->findNeighbor(node);
                
                link_nei = (TerraceNeighbor*) nei->link_neighbors[0];
                link_dad_nei = (TerraceNeighbor*) dad_nei->link_neighbors[0];
                
                assert(link_nei->link_neighbors.size()==link_dad_nei->link_neighbors.size());
            }
            if(branch_ids.empty()){
                if(induced_trees[i]->leafNum>2){
                    for(j=0; j<link_nei->link_neighbors.size(); j++){
                        branch_ids.push_back(link_nei->link_neighbors[j]->id);
                        
                        //intersection.push_back(link_nei->link_neighbors[j]->id);
                        /*node1_vec->push_back((TerraceNode*)link_nei->link_neighbors[j]->node);
                        node2_vec->push_back((TerraceNode*)link_dad_nei->link_neighbors[j]->node);*/
                    }
                    /*if(branch_ids.size()>0 && branch_ids[0]==352){
                        cout<<taxon_name<<"\n";
                        cout<<"INITIAL: branch_ids:\n";
                        for(j=0; j<branch_ids.size(); j++){
                            cout<<" "<<branch_ids[j];
                        }
                        cout<<"\n";
                    }*/
                } else {
                    //cout<<"All branches are allowed for this tree. Push them back"<<"\n";
                    /*for(k=0; k<branch_end_1->size(); k++){
                        aux_nei = branch_end_1->at(k)->findNeighbor(branch_end_2->at(k));
                        
                        branch_ids.push_back(aux_nei->id);
                        //intersection.push_back(aux_nei->id);
                        
                        node1_vec->push_back((TerraceNode*)branch_end_1->at(k));
                        node2_vec->push_back((TerraceNode*)branch_end_2->at(k));
                    }*/
                    for(k=0; k<branchNum;k++){
                        branch_ids.push_back(k);
                    }
                }
            } else if(induced_trees[i]->leafNum>2){ // otherwise, there is no constraint
                
                // TODO: The way to keep track of allowed branches is not optimal, too many unnessary savings and clearings. Optimize in the next step. Also the intersection is O(n^2).
                intersection.clear();
                /*node1_vec->clear();
                node2_vec->clear();*/
                
                // check, if the element occurs in branch ids
                for(j=0; j<link_nei->link_neighbors.size(); j++){
                    //cout<<"PART"<<i<<" | Allowed "<<j<<": "<<link_nei->link_neighbors[j]->id;
                    found = false;
                    for(h=0; h<branch_ids.size(); h++){
                        if(link_nei->link_neighbors[j]->id==branch_ids[h]){
                            found = true;
                            //cout<<" FOUND: "<<branch_ids[h];
                            break;
                        }
                    }
                    if(found){
                        intersection.push_back(link_nei->link_neighbors[j]->id);
                        /*node1_vec->push_back((TerraceNode*)link_nei->link_neighbors[j]->node);
                        node2_vec->push_back((TerraceNode*)link_dad_nei->link_neighbors[j]->node);*/
                        if(intersection.size()==branch_ids.size()){
                            break;
                        }
                    }
                    //cout<<"\n";
                }
                
                branch_ids.clear();
                
                if(intersection.size()>0){
                    //for(j=0; j<intersection.size(); j++){
                    //    branch_ids.push_back(intersection[j]);
                    //}
                    //intersection.clear();
                    /*cout<<"BEFORE:\nintersection:\n";
                    for(j=0; j<intersection.size(); j++){
                        cout<<" "<<intersection[j];
                    }
                    cout<<"\n";
                    cout<<"branch_ids:\n";
                    for(j=0; j<branch_ids.size(); j++){
                        cout<<" "<<branch_ids[j];
                    }
                    cout<<"\n";
                    */
                    branch_ids=std::move(intersection);
                    
                    /*cout<<"AFTER:\nintersection:\n";
                    for(j=0; j<intersection.size(); j++){
                        cout<<" "<<intersection[j];
                    }
                    cout<<"\n";*/
                    /*cout<<"branch_ids:\n";
                    for(j=0; j<branch_ids.size(); j++){
                        cout<<" "<<branch_ids[j];
                    }
                    cout<<"\n";*/
                } else {
                    /*node1_vec->clear();
                    node2_vec->clear();*/
                    break;
                }
            }
        }
    }
    
    // Collect branches
    /*if(branch_ids.size()>0){
        for(i=0; i<branch_ids.size();i++){
            for(k=0; k<branch_end_1.size(); k++){
                aux_nei = branch_end_1[k]->findNeighbor(branch_end_2[k]);
                if(aux_nei->id==branch_ids[i]){
                    node1_vec->push_back((TerraceNode*)branch_end_1[k]);
                    node2_vec->push_back((TerraceNode*)branch_end_2[k]);
                    branch_end_1[k]=nullptr;
                    branch_end_2[k]=nullptr;
                    branch_end_1.erase(branch_end_1.begin()+k);
                    branch_end_2.erase(branch_end_2.begin()+k);
                    break;
                }
            }
        }
    }*/
    if(branch_ids.size()>0){
        for(const auto& b: branch_ids){
            node1_vec->push_back((TerraceNode*)brNodes[b][0]);
            node2_vec->push_back((TerraceNode*)brNodes[b][1]);
        }
    }
    //cout<<"IN getALLOWEDbranches: "<<taxon_name<<":"<<node1_vec->size()<<"\n";
}

void Terrace::extendNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> part_tree_pairs,vector<int> pr_ab_info){
    
    //cout<<"=========================================================\n";
    //cout<<"GOING TO INSERT: "<<node_name<<"\n";
    //cout<<"=========================================================\n";
    
    
    TerraceNode* leaf_node;
    TerraceNeighbor *nei_1, *nei_2;
    nei_1 = (TerraceNeighbor*) node_1_branch->findNeighbor(node_2_branch);
    nei_2 = (TerraceNeighbor*) node_2_branch->findNeighbor(node_1_branch);
    
    NodeVector induced_part_tree_branch_1, induced_part_tree_branch_2;
    induced_part_tree_branch_1.resize(part_num,nullptr);
    induced_part_tree_branch_2.resize(part_num,nullptr);
    TerraceNeighbor *nei_part_1, *nei_part_2;
    NodeVector part_taxa;
    
    int i,j;

    for(i=0; i<part_num; i++){
        
        if(induced_trees[i]->leafNum>2){
            
            induced_part_tree_branch_1[i]=((TerraceNeighbor*)nei_2->link_neighbors[i])->node;
            induced_part_tree_branch_2[i]=((TerraceNeighbor*)nei_1->link_neighbors[i])->node;
            
            nei_part_1 = (TerraceNeighbor*)induced_part_tree_branch_1[i]->findNeighbor(induced_part_tree_branch_2[i]);
            nei_part_2 = (TerraceNeighbor*)induced_part_tree_branch_2[i]->findNeighbor(induced_part_tree_branch_1[i]);
        
            //if(part_tree_pairs[i]->matrix->findTaxonID(node_name)!=-1){
            if(part_tree_pairs[i]->leafNodes.find(node_name)!=part_tree_pairs[i]->leafNodes.end()){
                //cout<<"Partition:"<<i<<"| IN_extendNewTaxon:case>2, taxon will be inserted\n";
                
                // SANITY CHECK:
                //assert(part_tree_pairs[i]->matrix->findTaxonID(node_name)!=-1);
                // --------------- CLEARING MAPS FROM AGILE TREE TO COMMON PART SUBTREES and BACK --------------------------------------------------------------------------
                // STEP 1: clearing pointers for all branches that mapped to the involved branch on induced partition tree (clearing forward maps)
                for(j=0; j<nei_part_1->link_neighbors.size();j++){
                    ((TerraceNeighbor*)nei_part_1->link_neighbors[j])->link_neighbors[i] = nullptr;
                    ((TerraceNeighbor*)nei_part_2->link_neighbors[j])->link_neighbors[i] = nullptr;
                }
                
                // STEP 2: clearing backward maps from induced partition trees to parent tree
                nei_part_1->link_neighbors.clear();
                nei_part_2->link_neighbors.clear();
                
                
                // --------------- CLEARING MAPS FROM INDUCED PART TREE TO COMMON PART SUBTREES and BACK --------------------------------------------------------------------------
                // STEP 1: clearing pointers for all branches that mapped to the involved branch on induced partition tree (clearing forward maps)
                for(j=0; j<nei_part_1->link_neighbors_lowtop_back.size();j++){
                    ((TerraceNeighbor*)nei_part_1->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                    ((TerraceNeighbor*)nei_part_2->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                }
                
                // STEP 2: clearing backward maps from induced partition trees to parent tree
                nei_part_1->link_neighbors_lowtop_back.clear();
                nei_part_2->link_neighbors_lowtop_back.clear();
                
                //part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa();

                // --------------- INSERT A TAXON ON COMMON PARTITION SUBTREE ----------------------------------------------------------------------------------------------------
                leaf_node = induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) induced_part_tree_branch_1[i], (TerraceNode*) induced_part_tree_branch_2[i],true);
                //part_tree_pairs[i]->matrix->pr_ab_matrix[part_tree_pairs[i]->matrix->findTaxonID(node_name)][0]=1;

                // --------------- UPDATE MAPs LOCALLY for TOP LEVEL INDUCED PART TREES ------------------------------------------------------------------------------------------
                //part_taxa.clear();
                //part_tree_pairs[i]->matrix->getPartTaxa(0, part_tree_pairs[i], part_tree_pairs[i]->induced_trees[0], part_taxa);
                //TerraceNode *central_node_part = (TerraceNode*) part_tree_pairs[i]->findLeafName(node_name)->neighbors[0]->node;
                TerraceNode *central_node_part = (TerraceNode*) part_tree_pairs[i]->leafNodes[node_name]->neighbors[0]->node;
                part_tree_pairs[i]->update_map(0,part_taxa, false, true,central_node_part);
                
            } else {
                //cout<<"Partition:"<<i<<"| IN_extendNewTaxon:case>2, taxon will NOT be inserted\n";
                // for all partitions, which do not have the taxon of interest, remove the branch, which will be devided by the insertion of a new taxon from the backward map.
                // in the next step it will be substituted by the three branches: one with a new taxon and two ends of the devided branch
                for(j=0; j<nei_part_1->link_neighbors.size();j++){
                    if((nei_part_1->link_neighbors[j]->node == node_1_branch && nei_part_2->link_neighbors[j]->node == node_2_branch) or (nei_part_1->link_neighbors[j]->node == node_2_branch && nei_part_2->link_neighbors[j]->node == node_1_branch)){
                        
                        nei_part_1->link_neighbors.erase(nei_part_1->link_neighbors.begin()+j);
                        nei_part_2->link_neighbors.erase(nei_part_2->link_neighbors.begin()+j);
                        
                        break;
                    }
                }
            }
        } else if(part_tree_pairs[i]->leafNodes.find(node_name)!=part_tree_pairs[i]->leafNodes.end()){
            // SANITY CHECK:
            //assert(part_tree_pairs[i]->matrix->findTaxonID(node_name)!=-1);
            
            //if(part_tree_pairs[i]->matrix->findTaxonID(node_name)!=-1){
            //if(part_tree_pairs[i]->findLeafName(node_name)){
            
            //part_tree_pairs[i]->matrix->pr_ab_matrix[part_tree_pairs[i]->matrix->findTaxonID(node_name)][0]=1;
            
            if(induced_trees[i]->leafNum == 2){
                //cout<<"Partition:"<<i<<"| IN_extendNewTaxon:case==2, taxon will be inserted\n";
                assert(induced_trees[i]->root->isLeaf() && "ERROR: in extendNewTaxon: root is not a leaf... something is wrong");
                induced_part_tree_branch_1[i]=induced_trees[i]->root;
                induced_part_tree_branch_2[i]=induced_trees[i]->root->neighbors[0]->node;
                leaf_node = induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) induced_part_tree_branch_1[i], (TerraceNode*) induced_part_tree_branch_2[i],true);
                
                // Partitions with less than 3 taxa were not linked before. When you insert 3rd taxon, you should map them
                //NodeVector part_taxa;
                //part_taxa.clear();
                //part_tree_pairs[i]->matrix->getPartTaxa(0, part_tree_pairs[i], part_tree_pairs[i]->induced_trees[0], part_taxa);
                part_tree_pairs[i]->linkTree(0, part_taxa, false, true);
                
            } else {
                //cout<<"Partition:"<<i<<"| IN_extendNewTaxon:case<2, taxon will be inserted\n";
                leaf_node = induced_trees[i]->insertNewTaxon(node_name, nullptr, nullptr,true);
            }
        }
    }
    
    ((TerraceNeighbor*)node_1_branch->findNeighbor(node_2_branch))->link_neighbors.clear();
    ((TerraceNeighbor*)node_2_branch->findNeighbor(node_1_branch))->link_neighbors.clear();
    
    //cleanAllLinkNeighboursAndTaxa();
    //cout<<"=========================================================\n";
    //cout<<"MAIN insertion on agile tree\n";
    leaf_node = insertNewTaxon(node_name,node_1_branch,node_2_branch,true,true);
    //cout<<"=========================================================\n";
    taxa_num += 1;
    //matrix->extend_by_new_taxa(node_name, pr_ab_info);

    TerraceNode* center_node = (TerraceNode*) leaf_node->neighbors[0]->node;
    
    TerraceNeighbor *center_node_nei;
    //TerraceNeighbor *nei_aux;
    // INFO: since you are introducing new branches, make sure the link_neighbor vector is initialised for them
    FOR_NEIGHBOR_IT(center_node, NULL, it){
        center_node_nei=(TerraceNeighbor*)(*it)->node->findNeighbor(center_node);
        center_node_nei->link_neighbors.resize(part_num,nullptr);
        ((TerraceNeighbor*)(*it))->link_neighbors.resize(part_num,nullptr);
    }

    //cout<<"============================================\n";
    //cout<<"YES->INSERTED: "<<node_name<<"\n";
    //cout<<"============================================\n";
    //cout<<"Reached update part"<<"\n";
    for(i=0; i<part_num; i++){
        if(induced_trees[i]->leafNum>2){
            //if(induced_trees[i]->findLeafName(node_name)){
            //if(part_tree_pairs[i]->matrix->findTaxonID(node_name)!=-1){
            if(part_tree_pairs[i]->leafNodes.find(node_name)!=part_tree_pairs[i]->leafNodes.end()){
                // SANITY CHECK:
                //assert(part_tree_pairs[i]->matrix->findTaxonID(node_name)!=-1);
                //cout<<"Partition:"<<i<<"- larger than 2 - leaf occurs"<<"\n";
                // if a taxon was inserted to the induced partition tree, update
                //part_taxa.clear();
                //matrix->getPartTaxa(i, this, induced_trees[i], part_taxa);
                //cout<<"Partition:"<<i<<"| IN_extendNewTaxon_AFTER:case>2, taxon WAS inserted => update_map\n";
                update_map(i,part_taxa, true, false, center_node);

            }else{
                //cout<<"Partition:"<<i<<"| IN_extendNewTaxon_AFTER:case>2, taxon DOES NOT OCCUR => simple update\n";
                //cout<<"Partition:"<<i<<"- larger than 2 - leaf does not occur"<<"\n";
                
                //cout<<"if a taxon does not occur on the induced partition tree"<<"\n";
                nei_part_1 = (TerraceNeighbor*)induced_part_tree_branch_1[i]->findNeighbor(induced_part_tree_branch_2[i]);
                nei_part_2 = (TerraceNeighbor*)induced_part_tree_branch_2[i]->findNeighbor(induced_part_tree_branch_1[i]);
                
                /*//cout<<"First branch: one part of the devided branch"<<"\n";
                nei_aux = (TerraceNeighbor*)node_1_branch->findNeighbor(center_node);
                nei_aux->link_neighbors[i] = nei_part_1;
                nei_part_1->link_neighbors.push_back(nei_aux);
                
                nei_aux = (TerraceNeighbor*)center_node->findNeighbor(node_1_branch);
                nei_aux->link_neighbors[i] = nei_part_2;
                nei_part_2->link_neighbors.push_back(nei_aux);
                
                //cout<<"Second branch: another part of the devided branch"<<"\n";
                nei_aux = (TerraceNeighbor*)center_node->findNeighbor(node_2_branch);
                nei_aux->link_neighbors[i] = nei_part_1;
                nei_part_1->link_neighbors.push_back(nei_aux);
                
                nei_aux = (TerraceNeighbor*)node_2_branch->findNeighbor(center_node);
                nei_aux->link_neighbors[i] = nei_part_2;
                nei_part_2->link_neighbors.push_back(nei_aux);
                
                //cout<<"Third branch: incident to a newly inserted taxon: center->leaf"<<"\n";
                nei_aux = (TerraceNeighbor*)center_node->findNeighbor(leaf_node);
                nei_aux->link_neighbors[i] = nei_part_1;
                nei_part_1->link_neighbors.push_back(nei_aux);
		
                //cout<<"Third branch: incident to a newly inserted taxon: leaf->center"<<"\n";
                nei_aux = (TerraceNeighbor*)leaf_node->findNeighbor(center_node);
                nei_aux->link_neighbors[i] = nei_part_2;
                nei_part_2->link_neighbors.push_back(nei_aux);*/
                
            
                //assert(center_node->neighbors.size()==3 && "ERROR: The central node does not have 3 neighbours! Case leafNum>2 & findLeafName(node_name) = false.");
                FOR_NEIGHBOR_IT(center_node, NULL, it){
                    
                    center_node_nei=(TerraceNeighbor*)(*it)->node->findNeighbor(center_node);

                    // backward map
                    nei_part_1->link_neighbors.push_back(center_node_nei);
                    nei_part_2->link_neighbors.push_back((TerraceNeighbor*)(*it));
                    
                    // forward map
                    center_node_nei->link_neighbors[i] = nei_part_1;
                    ((TerraceNeighbor*)(*it))->link_neighbors[i]=nei_part_2;
                }
            }
        }
        //else{
            //cout<<"Partition:"<<i<<"| IN_extendNewTaxon_AFTER:case<=2 => do nothing\n";
            //cout<<"Partition:"<<i<<"- less than 2"<<"\n";
        //}
        
        /*if(i==1){
            cout<<"============================================\n";
            cout<<"PARTITION INFO: "<<i<<"\n";
            cout<<"============================================\n";
            part_tree_pairs[i]->printMapInfo();
            part_tree_pairs[i]->printBackMapInfo();
            cout<<"============================================\n";
        }*/
    }
    
    
    //cout<<"INTERMEDIATE_INFO_TAXA_"<<leafNum<<"_INSERTED_"<<node_name<<"_TREE_"<<"\n";
    //printTree(cout, WT_BR_SCALE | WT_NEWLINE);
    
    intermediated_trees_num +=1;
    if(intermediated_trees_num % 100000 == 0){
        cout<<"... trees generated - "<<intermediated_trees_num<<"; intermediated - "<<intermediated_trees_num-terrace_trees_num<<"; terrace - "<<terrace_trees_num<<"; dead paths - "<<dead_ends_num<<"\n";
    }
    
    if(intermediated_trees_num - terrace_trees_num == intermediate_max_trees){
        for(const auto &p: part_tree_pairs){
            p->unset_part_trees();
        }
        write_warning_stop(1);
    }
    
    double time= getCPUTime()-Params::getInstance().startCPUTime;
    if(seconds_max!=-1 and time > seconds_max){
        write_warning_stop(3);
    }
    
    // For DEBUGGING:
    /*bool clean_induced_part_maps = true;
    for(i=0; i<part_num; i++){
        part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
    }
    cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
    relinkALL(part_tree_pairs);*/
}

void Terrace::extendNewTaxon_naive(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> part_tree_pairs,vector<int> pr_ab_info){
    
    TerraceNode *leaf_node;
    TerraceNeighbor *nei_1, *nei_2;
    nei_1 = (TerraceNeighbor*) node_1_branch->findNeighbor(node_2_branch);
    nei_2 = (TerraceNeighbor*) node_2_branch->findNeighbor(node_1_branch);
    
    Neighbor *nei_part_1, *nei_part_2;
    
    NodeVector part_taxa;
    
    int i;
    
    for(i=0; i<part_num; i++){
        //if(part_tree_pairs[i]->findLeafName(node_name)){
        //if(part_tree_pairs[i]->matrix->findTaxonID(node_name)!=-1){
        if(part_tree_pairs[i]->leafNodes.find(node_name)!=part_tree_pairs[i]->leafNodes.end()){
            if(induced_trees[i]->leafNum>2){
                nei_part_1 = nei_1->link_neighbors[i];
                nei_part_2 = nei_2->link_neighbors[i];
                
                leaf_node = induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) nei_part_1->node, (TerraceNode*)nei_part_2->node,true);
                
            }else if(induced_trees[i]->leafNum == 2){
                assert(induced_trees[i]->root->isLeaf() && "ERROR: in extendNewTaxon: root is not a leaf... something is wrong");
                leaf_node = induced_trees[i]->insertNewTaxon(node_name, (TerraceNode*) induced_trees[i]->root, (TerraceNode*) induced_trees[i]->root->neighbors[0]->node,true);
            }else if(induced_trees[i]->leafNum == 1){
                TerraceNode* node_new = (TerraceNode*)induced_trees[i]->newNode(1, node_name.c_str());
                induced_trees[i]->root->addNeighbor(node_new, 0.0,0);
                node_new->addNeighbor(induced_trees[i]->root, 0.0,0);
                induced_trees[i]->leafNum += 1;
                induced_trees[i]->nodeNum += 1;
                induced_trees[i]->initializeTree();
            }else{
                assert(induced_trees[i]->leafNum == 0 && "ERROR in extendNewTaxon_naive: a tree should have 0 zero leaves...");
                induced_trees[i]->root = (TerraceNode*)induced_trees[i]->newNode(0, node_name.c_str());
                induced_trees[i]->leafNum += 1;
                induced_trees[i]->nodeNum += 1;
            }
            
            int id = part_tree_pairs[i]->matrix->findTaxonID(node_name);
            part_tree_pairs[i]->matrix->pr_ab_matrix[id][0] = 1;
            
        }
        part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(true);
    }
    
    cleanAllLinkNeighboursAndTaxa(true);
    leaf_node = insertNewTaxon(node_name, node_1_branch, node_2_branch,true);
    matrix->extend_by_new_taxa(node_name, pr_ab_info);
    taxa_num += 1;
    
    relinkALL(part_tree_pairs);
    
    intermediated_trees_num +=1;
    
}

void Terrace::generateTerraceTrees(Terrace *terrace, vector<Terrace*> &part_tree_pairs, vector<string> &list_taxa_to_insert, int taxon_to_insert, vector<string> *ordered_taxa_to_insert){
    
    string taxon_name;
    taxon_name = list_taxa_to_insert[taxon_to_insert];
    //NodeVector branch_node_1, branch_node_2;
    NodeVector node1_vec_branch, node2_vec_branch;
    if(ordered_taxa_to_insert){
        
        if(ordered_taxa_to_insert->size()>1){
            taxon_name = getNextTaxon(part_tree_pairs,ordered_taxa_to_insert,node1_vec_branch,node2_vec_branch);
            if(list_taxa_to_insert[taxon_to_insert]!=taxon_name){
                vector<string>::iterator it_t = std::find(list_taxa_to_insert.begin()+taxon_to_insert+1,list_taxa_to_insert.end(),taxon_name);
                assert(it_t!=list_taxa_to_insert.end() && "ERROR: in generateTerraceTrees: taxon not found in the remainder of the list!");
                list_taxa_to_insert.erase(it_t);
                list_taxa_to_insert.insert(list_taxa_to_insert.begin()+taxon_to_insert, taxon_name);
            }
        }else{
            if(ordered_taxa_to_insert->size()==1){
                ordered_taxa_to_insert->clear();
            }
            
            //getAllowedBranches(taxon_name, part_tree_pairs, &node1_vec_branch, &node2_vec_branch, &branch_node_1, &branch_node_2);
            getAllowedBranches(taxon_name, part_tree_pairs, &node1_vec_branch, &node2_vec_branch);
        }
    }else{
        //getAllowedBranches(taxon_name, part_tree_pairs, &node1_vec_branch, &node2_vec_branch,&branch_node_1, &branch_node_2);
        getAllowedBranches(taxon_name, part_tree_pairs, &node1_vec_branch, &node2_vec_branch);
    }
    
    //cout<<"\n"<<"*******************************************************"<<"\n";
    //cout<<"GENERATE_TERRACE_TREES | TAXON "<<taxon_name<<"\n";
    //cout<<"*******************************************************"<<"\n";
    
    
    int j, id;
    
    if(!node1_vec_branch.empty()){
        //cout<<"NUM_OF_ALLOWED_BRANCHES_"<<taxon_name<<"_"<<node1_vec_branch.size()<<"\n";
        //cout<<"ALL ALLOWED BRANCHES:"<<"\n";
        //for(j=0; j<node1_vec_branch.size(); j++){
        //    cout<<j<<":"<<node1_vec_branch[j]->id<<"-"<<node2_vec_branch[j]->id<<"\n";
        //}
        
        for(j=0; j<node1_vec_branch.size(); j++){
            //cout<<"-----------------------------------"<<"\n"<<"INSERTing taxon "<<taxon_name<<" on branch "<<j+1<<" out of "<<node1_vec_branch.size()<<": "<<node1_vec_branch[j]->id<<"-"<<node2_vec_branch[j]->id<<"\n"<<"-----------------------------------"<<"\n";
            
            id = terrace->matrix->findTaxonID(taxon_name);
            assert(id!=-1);
            extendNewTaxon(taxon_name,(TerraceNode*)node1_vec_branch[j],(TerraceNode*)node2_vec_branch[j],part_tree_pairs,terrace->matrix->pr_ab_matrix[id]);
            //this->printTree(cout,WT_NEWLINE);
            //extendNewTaxon_naive(taxon_name,(TerraceNode*)node1_vec_branch[j],(TerraceNode*)node2_vec_branch[j],part_tree_pairs,terrace->matrix->pr_ab_matrix[id]);
            
            if(taxon_to_insert != list_taxa_to_insert.size()-1){
                
                generateTerraceTrees(terrace, part_tree_pairs, list_taxa_to_insert, taxon_to_insert+1,ordered_taxa_to_insert);
                
                // INFO: IF NEXT TAXON DOES NOT HAVE ALLOWED BRANCHES CURRENT TAXON IS DELETED AND NEXT BRANCH IS EXPLORED.
                //remove_one_taxon_naive(taxon_name,part_tree_pairs);
                remove_one_taxon(taxon_name,part_tree_pairs);
            } else {
                if(terrace_out){
                    if(trees_out_lim==0 or terrace_trees_num<trees_out_lim){
                    //terrace_trees.push_back(getTreeTopologyString(this));
                        
                    //ofstream out;
                    //out.exceptions(ios::failbit | ios::badbit);
                    //out.open(out_file,std::ios_base::app);
                    printTree(out, WT_BR_SCALE | WT_NEWLINE);
                    //out.close();
                    }
                }
                terrace_trees_num+=1;
                //if(terrace_trees_num % 1000 == 0){
                //    cout<<"... generated tree "<<terrace_trees_num<<"\n";
                //}
                //printTree(cout, WT_BR_SCALE | WT_NEWLINE);
                if(terrace_trees_num == terrace_max_trees){
                    write_warning_stop(2);
                }
                //remove_one_taxon_naive(taxon_name,part_tree_pairs);
                remove_one_taxon(taxon_name,part_tree_pairs);
            }
        }
    } else {
        //cout<<"NUM_OF_ALLOWED_BRANCHES_"<<taxon_name<<"_0_dead_end"<<"\n";
        //cout<<"For a given taxon "<<taxon_name<<" there are no allowed branches.. Dead end.."<<"\n";
        dead_ends_num +=1;
    }
    
    if(ordered_taxa_to_insert){
        taxon_name = list_taxa_to_insert[taxon_to_insert];
        ordered_taxa_to_insert->insert(ordered_taxa_to_insert->begin(),taxon_name);
        //cout<<"Added taxon back: "<<taxon_name<<"|ordered_taxa_to_insert->size()="<<ordered_taxa_to_insert->size()<<"\n";
    }

}

void Terrace::remove_one_taxon(string taxon_name, vector<Terrace*> part_tree_pairs){
    
    //cout<<"-----------------------------------"<<"\n"<<"REMOVING TAXON: "<<taxon_name<<"\n"<<"-----------------------------------"<<"\n";
    
    int i, j, h;
    NodeVector induced_part_tree_branch_1, induced_part_tree_branch_2;
    induced_part_tree_branch_1.resize(part_num,nullptr);
    induced_part_tree_branch_2.resize(part_num,nullptr);
    
    Node *node_1, *node_2; // branch nodes of the agile (parent) tree, that will be the ends of the joint branch after taxon removal
    TerraceNeighbor *neiC1, *nei1C;
    
    Node *node_leaf_main, *node_central_main;
    Node *node_leaf, *node_central;
    NodeVector branch_nodes;
    TerraceNeighbor *nei1, *nei2, *nei_aux;
    
    //node_leaf_main = findLeafName(taxon_name);
    node_leaf_main = leafNodes[taxon_name];
    node_central_main = node_leaf_main->neighbors[0]->node;
    FOR_NEIGHBOR_DECLARE(node_central_main, node_leaf_main, it){
        branch_nodes.push_back((*it)->node);
    }
    node_1 = branch_nodes[0];
    node_2 = branch_nodes[1];
    branch_nodes.clear();
    
    // one of the two branches that will be joint on the agile tree
    neiC1 = (TerraceNeighbor*) node_central_main->findNeighbor(node_1);
    nei1C = (TerraceNeighbor*) node_1->findNeighbor(node_central_main);
    
    // Below NeiVectors will store neibours from agile tree and top induced trees that need to be updated
    // Corresponding neibours will have their forward map equal to nei_low_12 on corresponding low induced subtree, i.e. I'm going to "fold" all nei's in the same direction
    // nei_12 ->link = nei_low_12, while nei_21 ->link = nei_low_21
    vector<NeighborVec> update_nei; // neis on agile tree, for each partition its own set of neis to be updated
    update_nei.resize(part_num);
    
    vector<NeighborVec> update_nei_top; // neis on top induced partition tree, for each pair its own set of neis to be updated
    update_nei_top.resize(part_num);
    
    NeighborVec update_nei_aux, update_nei_top_aux;
    bool appear;
    
    for(i=0; i<part_num; i++){
        //cout<<"-------------------------------"<<"\n";
        //cout<<"\n"<<"Partition "<<i<<":"<<"\n";
        
        update_nei_aux.clear();
        update_nei_top_aux.clear();
        
        //part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(true);
        
        /*
         *  If the taxon to be removed is present in partition, clear corresponding maps (forward and backward maps),
         *  which involve three branches (6 neighbors): incident to the taxon and two parts of the branch to be joint after taxon removal.
         *  For all partitions collect branch nodes of common partition subtrees, whose maps (forward or backward) have to
         *  be modified (clear and update).
         */
        
        //if(part_tree_pairs[i]->findLeafName(taxon_name)){
        //if(part_tree_pairs[i]->matrix->findTaxonID(taxon_name)!=-1){
        if(part_tree_pairs[i]->leafNodes.find(taxon_name)!=part_tree_pairs[i]->leafNodes.end()){
            //cout<<"- with Leaf"<<"\n";
            
            // SANITY CHECK:
            //assert(part_tree_pairs[i]->matrix->findTaxonID(taxon_name)!=-1);
            
            //int id = part_tree_pairs[i]->matrix->findTaxonID(taxon_name);
            //part_tree_pairs[i]->matrix->pr_ab_matrix[id][0] = 0;
            
            if(induced_trees[i]->leafNum>2){
                
                //cout<<"- with more than 2 leaves:"<<induced_trees[i]->leafNum<<"\n";
                //node_leaf = induced_trees[i]->findLeafName(taxon_name);
                node_leaf = induced_trees[i]->leafNodes[taxon_name];
                node_central = node_leaf->neighbors[0]->node;
                
                
                branch_nodes.clear();
                FOR_NEIGHBOR(node_central, node_leaf, it){
                    branch_nodes.push_back((*it)->node);
                }
                
                induced_part_tree_branch_1[i] = branch_nodes[0];
                induced_part_tree_branch_2[i] = branch_nodes[1];
                
                if(neiC1->link_neighbors[i]->node == branch_nodes[1] or nei1C->link_neighbors[i]->node == branch_nodes[1]){
                    induced_part_tree_branch_1[i] = branch_nodes[1];
                    induced_part_tree_branch_2[i] = branch_nodes[0];
                }
                
                /*if(neiC1->link_neighbors[i]->id == neiC2->link_neighbors[i]->id){
                    induced_part_tree_branch_1[i] = neiC1->link_neighbors[i]->node;
                    induced_part_tree_branch_2[i] = ((TerraceNeighbor*)node_1->findNeighbor(node_central_main))->link_neighbors[i]->node;
                    
                    if(induced_part_tree_branch_1[i] == node_leaf or induced_part_tree_branch_2[i] == node_leaf){
                        branch_nodes.clear();
                        FOR_NEIGHBOR(node_central, node_leaf, it){
                            branch_nodes.push_back((*it)->node);
                        }
                        induced_part_tree_branch_1[i] = branch_nodes[0];
                        induced_part_tree_branch_2[i] = branch_nodes[1];
                    }
                }else{
                    induced_part_tree_branch_1[i] = neiC1->link_neighbors[i]->node;
                    induced_part_tree_branch_2[i] = neiC2->link_neighbors[i]->node;
                }*/
                
                //cout<<"induced_part_tree_branch_1[i] = "<<induced_part_tree_branch_1[i]->id<<"\n";
                //cout<<"induced_part_tree_branch_2[i] = "<<induced_part_tree_branch_2[i]->id<<"\n";
                //cout<<"node_leaf="<<node_leaf->id<<"\n";
                //cout<<"node_central="<<node_central->id<<"\n";
                
                FOR_NEIGHBOR(node_central, nullptr, it){
                    nei1=(TerraceNeighbor*)(*it);
                    nei2=(TerraceNeighbor*)(*it)->node->findNeighbor(node_central);
                    
                    if(!((TerraceNeighbor*)(*it))->link_neighbors.empty()){
                        //cout<<"\n"<<"DEBUGGING: clearing forward nei's maps"<<"\n";
                        // clearing forward maps from agile (parent) tree to common partition subtree
                        for(j=0; j<nei1->link_neighbors.size(); j++){
                            //cout<<"Neibour "<<j<<":"<<"\n";
                            //cout<<"- nei1->link_neighbors[j]="<<nei1->link_neighbors[j]->node->id<<"\n";
                            //cout<<"- nei2->link_neighbors[j]="<<nei2->link_neighbors[j]->node->id<<"\n";
                            ((TerraceNeighbor*)nei1->link_neighbors[j])->link_neighbors[i] = nullptr;
                            ((TerraceNeighbor*)nei2->link_neighbors[j])->link_neighbors[i] = nullptr;
                            
                            
                            // INFO: you do not want any of the 6 neighbours that will be modified to appear in the list of neis to be updated, because these neis won't exist after removal and you'll get SEG_FAULT
                            appear=false;
                            FOR_NEIGHBOR_DECLARE(node_central_main, nullptr, it1){
                                if(nei1->link_neighbors[j] == (*it1)){
                                    appear = true;
                                    break;
                                }
                                if(nei1->link_neighbors[j] == (*it1)->node->findNeighbor(node_central_main)){
                                    appear = true;
                                    break;
                                }
                            }
                            if(!appear){
                                update_nei_aux.push_back(nei1->link_neighbors[j]);
                                //cout<<"+ added nei1->"<<nei1->link_neighbors[j]->node->id<<" to update list"<<"\n";
                            }
                            
                            appear=false;
                            FOR_NEIGHBOR(node_central_main, nullptr, it1){
                                if(nei2->link_neighbors[j] == (*it1)){
                                    appear = true;
                                    break;
                                }
                                if(nei2->link_neighbors[j] == (*it1)->node->findNeighbor(node_central_main)){
                                    appear = true;
                                    break;
                                }
                            }
                            if(!appear){
                                update_nei_aux.push_back(nei2->link_neighbors[j]);
                                //cout<<"+ added nei2->"<<nei2->link_neighbors[j]->node->id<<" to update list"<<"\n";
                            }
                        }
                        
                        // clearing backward maps from common partition subtree to agile (parent) tree
                        nei1->link_neighbors.clear();
                        nei2->link_neighbors.clear();
                    }
                
                    if(!((TerraceNeighbor*)(*it))->link_neighbors_lowtop_back.empty()){
                        
                        //cout<<"clearing forward maps from induced partition tree (parent) to common partition subtree"<<"\n";
                        for(j=0; j<((TerraceNeighbor*)(*it))->link_neighbors_lowtop_back.size(); j++){
                            ((TerraceNeighbor*)nei1->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                            ((TerraceNeighbor*)nei2->link_neighbors_lowtop_back[j])->link_neighbors[0] = nullptr;
                            
                            // INFO: note, that for top partition trees we do not need to avoid adding 6 neis to the list, because they will still exist after taxon deletion. However, when you'll update, set all neis to nei_low_12 and nei_low_21 correctly to nei_low_21 and not to nei_low_12
                            update_nei_top_aux.push_back(nei1->link_neighbors_lowtop_back[j]);
                            update_nei_top_aux.push_back(nei2->link_neighbors_lowtop_back[j]);
                            //cout<<"update_nei_top_aux.size()="<<update_nei_top_aux.size()<<"\n";
                        }
                        
                        // clearing backward maps from common partition subtree to induced partition tree (parent)
                        nei1->link_neighbors_lowtop_back.clear();
                        nei2->link_neighbors_lowtop_back.clear();
                    }
                }
                
                
                // Get your "to be updated list of neighbours" in place for agile tree
                update_nei[i] = update_nei_aux;
                
                // Get your "to be updated list of neighbours" in place for top induced partition tree
                update_nei_top[i] = update_nei_top_aux;
                
            } //else {
              //  cout<<"- with less than 2 leaves:"<<induced_trees[i]->leafNum<<"\n";
            //}
        } else {
            //cout<<"- without Leaf"<<"\n";
            if(induced_trees[i]->leafNum>2){
                //cout<<"- with more than 2 leaves:"<<induced_trees[i]->leafNum<<"\n";
                
                induced_part_tree_branch_1[i] = neiC1->link_neighbors[i]->node;
                induced_part_tree_branch_2[i] = nei1C->link_neighbors[i]->node;
                
                //cout<<"induced_part_tree_branch_1[i] = "<<induced_part_tree_branch_1[i]->id<<"\n";
                //cout<<"induced_part_tree_branch_2[i] = "<<induced_part_tree_branch_2[i]->id<<"\n";
                
                /*
                 * INFO: For partitions, which do not have corresponding taxon, remove 6 possible backward neighbors corresponding to three branches, which will be removed/modified
                 */
                nei1 = (TerraceNeighbor*)neiC1->link_neighbors[i];
                nei2 = (TerraceNeighbor*)nei1C->link_neighbors[i];
                
                for(j=nei1->link_neighbors.size()-1; j!=-1;j--){
                    for(h=0; h<node_central_main->neighbors.size();h++){
                        nei_aux = (TerraceNeighbor*)(node_central_main->neighbors[h]->node->findNeighbor(node_central_main));
                        if(nei1->link_neighbors[j] == node_central_main->neighbors[h]){
                            nei1->link_neighbors.erase(nei1->link_neighbors.begin()+j);
                            assert(nei2->link_neighbors[j] == nei_aux);
                            nei2->link_neighbors.erase(nei2->link_neighbors.begin()+j);
                        } else if(nei2->link_neighbors[j] == node_central_main->neighbors[h]){
                            nei2->link_neighbors.erase(nei2->link_neighbors.begin()+j);
                            assert(nei1->link_neighbors[j] == nei_aux);
                            nei1->link_neighbors.erase(nei1->link_neighbors.begin()+j);
                        }
                    }
                }
            }
            //else{
            //    cout<<"- with less than 2 leaves:"<<induced_trees[i]->leafNum<<"\n";
            //}
        }
    }
    
    //cleanAllLinkNeighboursAndTaxa(true);

    //cout<<"Removing the taxon from the low-level induced partition tree...."<<"\n";
    for(i=0; i<part_num; i++){
        if(induced_trees[i]->leafNum>1){
            //if(induced_trees[i]->findLeafName(taxon_name)){
            //if(part_tree_pairs[i]->matrix->findTaxonID(taxon_name)!=-1){
            if(part_tree_pairs[i]->leafNodes.find(taxon_name)!=part_tree_pairs[i]->leafNodes.end()){
                induced_trees[i]->remove_taxon(taxon_name,true);
            }
        } else if(induced_trees[i]->root){
            if(induced_trees[i]->root->name == taxon_name){
                induced_trees[i]->remove_taxon(taxon_name,true);
            }
        }
    }

    //cout<<"Removing the taxon from agile tree...."<<"\n";
    remove_taxon(taxon_name,true,true);
    
    ((TerraceNeighbor*)node_1->findNeighbor(node_2))->link_neighbors.resize(part_num,nullptr);
    ((TerraceNeighbor*)node_2->findNeighbor(node_1))->link_neighbors.resize(part_num,nullptr);
    
    //matrix->remove_taxon(taxon_name);
    taxa_num -= 1;
    
    //cout<<"Updating maps for partition trees...."<<"\n";
    // RE-MAP with new UPDATE using vectors of "to be modified" and not the update_map function
    for(i=0; i<part_num; i++){
        //cout<<"Updating maps for partition "<<i<<":"<<"\n";
        if(induced_trees[i]->leafNum>2){
            nei1 = (TerraceNeighbor*) induced_part_tree_branch_1[i]->findNeighbor(induced_part_tree_branch_2[i]);
            nei2 = (TerraceNeighbor*) induced_part_tree_branch_2[i]->findNeighbor(induced_part_tree_branch_1[i]);
            
            // AGILE TREE
            if(!update_nei[i].empty()){
                //cout<<"There are "<<update_nei[i].size()/2<<" branches ("<<update_nei[i].size()<<" neighbors) to be updated"<<"\n";
                for(j=0;j<update_nei[i].size()/2; j++){
                    ((TerraceNeighbor*)update_nei[i][2*j])->link_neighbors[i] = nei1;
                    nei1->link_neighbors.push_back(update_nei[i][2*j]);
                    ((TerraceNeighbor*)update_nei[i][2*j+1])->link_neighbors[i] = nei2;
                    nei2->link_neighbors.push_back(update_nei[i][2*j+1]);
                }
            }
            ((TerraceNeighbor*)node_1->findNeighbor(node_2))->link_neighbors[i] = nei1;
            nei1->link_neighbors.push_back(node_1->findNeighbor(node_2));
            ((TerraceNeighbor*)node_2->findNeighbor(node_1))->link_neighbors[i] = nei2;
            nei2->link_neighbors.push_back(node_2->findNeighbor(node_1));
            
            // TOP induced partition trees
            if(!update_nei_top[i].empty()){
                //cout<<"TOP PART TREE: There are "<<update_nei_top[i].size()/2<<" branches ("<<update_nei_top[i].size()<<" neighbors) to be updated"<<"\n";
                for(j=0;j<update_nei_top[i].size()/2; j++){
                    ((TerraceNeighbor*)update_nei_top[i][2*j])->link_neighbors[0] = nei1;
                    nei1->link_neighbors_lowtop_back.push_back(update_nei_top[i][2*j]);
                    ((TerraceNeighbor*)update_nei_top[i][2*j+1])->link_neighbors[0] = nei2;
                    nei2->link_neighbors_lowtop_back.push_back(update_nei_top[i][2*j+1]);
                }
            }
        }
    }
}

void Terrace::remove_one_taxon_naive(string taxon_name, vector<Terrace*> part_tree_pairs){
    //cout<<"-----------------------------------"<<"\n"<<"REMOVING TAXON: "<<taxon_name<<"\n"<<"-----------------------------------"<<"\n";
    
    bool clean_induced_part_maps = true;
    
    for(int i=0; i<part_num; i++){
        part_tree_pairs[i]->cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
        //if(part_tree_pairs[i]->findLeafName(taxon_name)){
        int id = part_tree_pairs[i]->matrix->findTaxonID(taxon_name);
        if(id!=-1){
            part_tree_pairs[i]->matrix->pr_ab_matrix[id][0] = 0;
        }
    }
    
    cleanAllLinkNeighboursAndTaxa(clean_induced_part_maps);
    
    for(int i=0; i<part_num; i++){
        //if(induced_trees[i]->findLeafName(taxon_name)){
        //if(part_tree_pairs[i]->matrix->findTaxonID(taxon_name)!=-1){
        if(part_tree_pairs[i]->leafNodes.find(taxon_name)!=part_tree_pairs[i]->leafNodes.end()){
            induced_trees[i]->remove_taxon(taxon_name,true);
        }
    }
    
    remove_taxon(taxon_name,true);
    matrix->remove_taxon(taxon_name);
    taxa_num -= 1;
    
    relinkALL(part_tree_pairs);
}

void Terrace::relinkALL(vector<Terrace*> part_tree_pairs){
    
    // re-link
    linkTrees(true, false);
    
    for(int i=0; i<part_tree_pairs.size(); i++){
        part_tree_pairs[i]->linkTrees(false, true);
    }
}

void Terrace::clearEmptyBranchAndTaxaINFO(TerraceNode *node, TerraceNode *dad){
    
    if (!node) {
        if (!root->isLeaf())
            node = (TerraceNode*) root;
        else
            node = (TerraceNode*) root->neighbors[0]->node;
        ASSERT(node);
    }
    
    if(node->empty_br_node_nei.size()>0){
        node->empty_br_node_nei.clear();
        node->empty_br_dad_nei.clear();
    }
    
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        clearEmptyBranchAndTaxaINFO((TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
}

void Terrace::cleanAllLinkNeighboursAndTaxa(bool clean_induced_part_maps){
    cleanAllLinkINFO(clean_induced_part_maps);
}

void Terrace::print_ALL_DATA(vector<Terrace*> part_tree_pairs){
    int i;
    
    cout<<"\n"<<"================ BEGIN: PRINTING all INFO ================"<<"\n";
    cout<<"Initial PARENT tree"<<"\n";
    print_terrace_tree();
    
    matrix->print_pr_ab_matrix();
    
    cout<<"================ INDUCED PARTITION PAIRS =========="<<"\n";
    for(i=0; i<part_num; i++){
        cout<<"------------------> TOP PART TREE "<<i<<":"<<"\n";
        part_tree_pairs[i]->print_terrace_tree();
        part_tree_pairs[i]->matrix->print_pr_ab_matrix();
        
        cout<<"\n"<<"------------------> LOW part tree "<<i<<":"<<"\n";
        induced_trees[i]->print_terrace_tree();
        cout<<"\n";
        cout<<"\n"<<"========================================================="<<"\n";
    }
    
    printMapInfo();
    printBackMapInfo();
    cout<<"================ END: PRINTING all INFO ================="<<"\n"<<"\n";

}

void Terrace::renameTaxa(){
    
    getTaxaName(taxa_names_orgn);
    
    NodeVector taxa;
    getTaxa(taxa);
    
    for(int i=0; i<taxa_num; i++){
        stringstream ss;
        ss << i;
        string name = "t" + ss.str();
        for(NodeVector::iterator it = taxa.begin(); it!=taxa.end(); it++){
            if((*it)->name == taxa_names_orgn[i]){
                (*it)->name = name;
                break;
            }
        }
        for(int j=0; j<taxa_num; j++){
            if(matrix->taxa_names[j] == taxa_names_orgn[i]){
                matrix->taxa_names[j] = name;
                break;
            }
        }
    }
}

bool Terrace::check_two_trees(MTree* query_tree){
    TerraceTree tree;
    tree.copyTree_byTaxonNames(query_tree, matrix->taxa_names);
    Terrace *query_terrace = new Terrace(tree, matrix);
    
    vector<string> trees, taxa_names;
    
    for(int part=0; part<part_num; part++){
        
        trees.clear();
        trees.push_back(getTreeTopologyString(induced_trees[part]));
        //trees.push_back(getTreeTopologyString(induced_trees_query[part]));
        
        taxa_names.clear();
        induced_trees[part]->getTaxaName(taxa_names);
        
        MTreeSet tree_set;
        tree_set.init(trees,taxa_names,induced_trees[part]->rooted);
        //tree_set[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
        
        MTreeSet tree_set_1;
        trees.clear();
        trees.push_back(getTreeTopologyString(query_terrace->induced_trees[part]));
        tree_set_1.init(trees,taxa_names,query_terrace->induced_trees[part]->rooted);
        //tree_set_1[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
        
        double *rfdist;
        int n=1;
        rfdist = new double [n*n];
        memset(rfdist, 0, n*n* sizeof(double));
        
        tree_set.computeRFDist(rfdist,&tree_set_1,true);
        //cout<<"Partition "<<part+1<<": RF-distance between induced trees is "<<rfdist[0]<<"\n";
        
        if(rfdist[0]!=0){
            //cout<<"---------------------------------------------"<<"\n";
            //cout<<"A query tree does not belong to the terrace:"<<"\n";
            //tree->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            
            //cout<<"The first encountered pair of induced partition trees that do not match. Partition "<<part+1<<"\n";
            //tree_set[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            //tree_set_1[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            
            delete [] rfdist;
            tree_set.clear();
            return false;
        }
        delete [] rfdist;
        tree_set.clear();
        
    }
    //cout<<"---------------------------------------------"<<"\n";
    //cout<<"A query tree is ON the terrace"<<"\n";
    
    return true;
}

void Terrace::write_terrace_trees_to_file(){
 
    cout<<"Wall-clock time used so far: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")"<<"\n";
    cout<<"CPU time used on generation of terrace trees: "
    << getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")" << "\n";
    cout<<"---------------------------------------------------------"<<"\n";
    
    cout<<"Printing "<<terrace_trees_num<<" terrace trees to file "<<"\n"<<out_file<<"..."<<"\n";
        
    if(trees_out_lim==0 or trees_out_lim > terrace_trees_num){
        trees_out_lim = terrace_trees_num;
    } else {
        cout<<"WARNING: The number of generated trees from the terrace ("<<terrace_trees_num<<") is larger than the output treshold ("<<trees_out_lim<<" trees). Only "<<trees_out_lim<<" trees will be written to the file."<<"\n";
    }
    
    ofstream out;
    out.exceptions(ios::failbit | ios::badbit);
    out.open(out_file,std::ios_base::app);
        
    for(int i=0; i<trees_out_lim; i++){
        out<<terrace_trees[i]<<"\n";
    }
    out.close();

}

void Terrace::write_summary_generation(){

    //Params::getInstance().run_time = (getCPUTime() - Params::getInstance().startCPUTime);
    
    if(rm_leaves > 0){
        cout<<"---------------------------------------------------------"<<"\n";
        cout<<"WARNING: Note, that the terrace trees were generated using BACKWARD approach.\nIt DOES NOT guarantee to generate all terrace trees and is only provided for exploratory purposes.\nIn current setting "<<rm_leaves<<" were removed from the input tree and used to generate terrace trees. \nYou can increase this number via option -t_rm_leaves <num> to explore terrace further.\n";
    }
    
    cout<<"---------------------------------------------------------"<<"\n";
    cout<<"SUMMARY:"<<"\n";
    cout<<"Number of trees on terrace: "<<terrace_trees_num<<"\n";
    cout<<"Number of intermediated trees visited: "<<intermediated_trees_num - terrace_trees_num<<"\n";
    cout<<"Number of dead ends encountered: "<<dead_ends_num<<"\n";
    cout<<"---------------------------------------------------------"<<"\n";
    
    // If there is an input tree (root != nullptr), try BACKWARD approach
    if(master_terrace->root){
        if(terrace_trees_num==0 && rm_leaves==-1){
            
            if(terrace_out){
                out.close();
            }
            cout<<"Current wall-clock time used: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")"<<"\n";
            cout<<"Current CPU time used: "
            << getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")" << "\n";
            cout<<"\n";
            
            cout<<"---------------------------------------------------------"<<"\n";
            cout<<"WARNING: Due to complexity of the input data, the FORWARD approach could not generate terrace trees.\nThe BACKWARD approach is triggered for exploratory purposes.\nThis approach DOES NOT guarantee generating all terrace trees.\n";
            cout<<"---------------------------------------------------------"<<"\n";
            
            Terrace *master = master_terrace;
            const int m = 10; // Maybe % of missing taxa?
            delete this;
            run_generate_trees(master, Params::getInstance(),m);
            
            exit(0);
        }
    
    }
    
    if(terrace_out){
        out.close();
        //write_terrace_trees_to_file();
        //cout<<"---------------------------------------------------------"<<"\n";
    }
    cout<<"Total wall-clock time used: "<<getRealTime()-Params::getInstance().start_real_time<<" seconds ("<<convert_time(getRealTime()-Params::getInstance().start_real_time)<<")"<<"\n";
    cout<<"Total CPU time used: "
    << getCPUTime()-Params::getInstance().startCPUTime << " seconds (" << convert_time(getCPUTime()-Params::getInstance().startCPUTime) << ")" << "\n";
    cout<<"\n";
    exit(0);
}

void Terrace::write_warning_stop(int type){
    
    cout<<"\n"<<"=========================================="<<"\n";
    cout<<"WARNING: stopping condition is active!"<<"\n";
    cout<<"The total number of trees on the terrace is NOT yet computed!"<<"\n";
    cout<<"Check summary at the current step. You can change the stopping rule via corresponding option (see below)."<<"\n";
    cout<<"=========================================="<<"\n";
    
    switch (type) {
        case 1:
            cout<<"Type of stopping rule: number of visited intermediate trees"<<"\n";
            cout<<"Current setting: stop if the number of visited intermediate trees reached "<<intermediate_max_trees<<"\n";
            cout<<"To change the value use -t_stop_i <number_of_intermediate_trees_to_stop>"<<"\n";
            cout<<"------------------------------------------"<<"\n";
            cout<<"The number of considered intermediated trees already reached "<<intermediate_max_trees<<". Exiting generation process..."<<"\n";
            break;
        case 2:
            cout<<"Type of stopping rule: terrace size"<<"\n";
            cout<<"Current setting: stop if the number of trees on the terrace reached "<<terrace_max_trees<<"\n";
            cout<<"To change the value use -t_stop_t <number_of_terrace_trees_to_stop>"<<"\n";
            cout<<"------------------------------------------"<<"\n";
            cout<<"Considered terrace already contains "<<terrace_max_trees<<" trees."<<"\n"<<"Exiting generation process..."<<"\n";
            break;
            
        case 3:
            cout<<"Type of stopping rule: CPU time used"<<"\n";
            cout<<"Current setting: stop if the CPU time used is larger than "<<seconds_max<<"\n";
            cout<<"To change the value use -t_stop_h <number_of_hours_to_stop>"<<"\n";
            cout<<"------------------------------------------"<<"\n";
            cout<<"Total CPU time is already larger than "<<seconds_max<<" trees."<<"\n"<<"Exiting generation process..."<<"\n";
            break;
            
        default:
            break;
    }
    
    cout<<"Printing summary at current step..."<<"\n";
    write_summary_generation();
    time_t end_time;
    time(&end_time);
    cout << "Date and Time: " << ctime(&end_time);
    exit(0);
}

string Terrace::getNextTaxon(vector<Terrace*> &part_tree_pairs, vector<string> *ordered_taxa_to_insert,NodeVector &node1_vec_main, NodeVector &node2_vec_main){
    
    string taxon_name;
    int len = 2*taxa_num-3;
    vector<string>::iterator it_NEO;
    
    //NodeVector branch_end_1, branch_end_2;
    //this->getBranches(branch_end_1, branch_end_2);
    
    int bond = 0, diff = ordered_taxa_to_insert->size()-matrix->uniq_taxa_num;
    //cout<<ordered_taxa_to_insert->size()<<"|"<<matrix->uniq_taxa_num<<"|"<<diff<<"|"<<ordered_taxa_to_insert->size()-matrix->uniq_taxa_num<<"\n";
    if(diff > 0){
        bond=matrix->uniq_taxa_num;
    }
    
    for(auto it=ordered_taxa_to_insert->begin(); it!=ordered_taxa_to_insert->end()-bond;it++){
        NodeVector node1_vec_branch, node2_vec_branch;
        //getAllowedBranches((*it), part_tree_pairs, &node1_vec_branch, &node2_vec_branch, &branch_end_1, &branch_end_2);
        getAllowedBranches((*it), part_tree_pairs, &node1_vec_branch, &node2_vec_branch);
        if(node1_vec_branch.size()==0){
            node1_vec_main.clear();
            node2_vec_main.clear();
            taxon_name = (*it);
            ordered_taxa_to_insert->erase(it);
            return taxon_name;
        }else if(node1_vec_branch.size()<len){
            len=node1_vec_branch.size();
            taxon_name = (*it);
            it_NEO = it;
            
            node1_vec_main.clear();
            node2_vec_main.clear();
            node1_vec_main=move(node1_vec_branch);
            node2_vec_main=move(node2_vec_branch);
        }else if(len==2*taxa_num-3 && node1_vec_branch.size()==len){
            taxon_name = (*it);
            it_NEO = it;
            
            node1_vec_main.clear();
            node2_vec_main.clear();
            node1_vec_main=move(node1_vec_branch);
            node2_vec_main=move(node2_vec_branch);
        }
    }

    ordered_taxa_to_insert->erase(it_NEO);
    return taxon_name;
};
