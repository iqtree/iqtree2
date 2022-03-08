/*
 *  terracetree.cpp
 *
 *  Created on: Sep 10, 2020
 *      Author: Olga
 */

#include "terracetree.hpp"
#include "terracenode.hpp"

TerraceTree::TerraceTree(): MTree(){
    branchNum = 0;
    leafNum = 0;
    nodeNum = 0;
    root = nullptr;
};
TerraceTree::~TerraceTree(){
   
    if(brNodes.size()>0){
        for(int id=0; id<branchNum; id++){
            brNodes[id]={nullptr,nullptr};
            brNodes.erase(id);
        }
    }
    if(leafNodes.size()>0){
        for(auto t=leafNodes.begin(); t!=leafNodes.end();t++){
            leafNodes.erase(t);
        }
    }

    if (root != NULL)
        freeNode();
    root = NULL;
};


void TerraceTree::readTree(const char *infile, bool &is_rooted){
    MTree::readTree(infile, is_rooted);
};

void TerraceTree::readTree(istream &in, bool &is_rooted) {
    MTree::readTree(in, is_rooted);
}

Node* TerraceTree::newNode(int node_id, const char* node_name) {
    Node* node = (Node*) (new (std::nothrow) TerraceNode(node_id, node_name));
    assert(node && "ERROR: cannot allocate memory... Exiting...");
    return node;
}

Node* TerraceTree::newNode(int node_id, int node_name) {
    Node* node = (Node*) (new (std::nothrow) TerraceNode(node_id, node_name));
    assert(node && "ERROR: cannot allocate memory... Exiting...");
    return node;
}


void TerraceTree::copyTree_byTaxonNames(MTree *tree, vector<string> taxon_names){
    
    //cout<<"Copying a tree using a vector of taxon names..."<<"\n";
    
    int i,j, sum = 0;
    int taxa_num = taxon_names.size();
    //bool check;
    
    string taxa_set = "";
    NodeVector taxa_nodes;
    
    NodeVector::iterator it2;
    vector<uint32_t> check_int;
    check_int.resize(tree->leafNum,0);
    
    tree->getTaxa(taxa_nodes);
    
    for(j=0; j<taxa_nodes.size(); j++){
        //check = false;
        for(i=0; i<taxa_num; i++){
            if(taxa_nodes[j]->name == taxon_names[i]){
                check_int[taxa_nodes[j]->id]=1;
                sum+=1;
                //check = true;
                break;
            }
        }
        //cout<<"Taxon["<<j<<"] = "<<check_int[j]<<"| main_tree:"<<taxa_nodes[j]->name<<"("<<taxa_nodes[j]->id<<") "<<"-> subtree: "<<((check) ? taxon_names[i] : "")<<"\n";
    }
    assert(sum == taxa_num && "Not all of the taxa appear in the complete tree!");
    taxa_set.clear();
    taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
    copyTree(tree,taxa_set);

    //printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
}

void TerraceTree::cleanAllLinkINFO(bool clean_induced_part_maps, TerraceNode *node, TerraceNode *dad){

    if(!node){
        if(root->isLeaf()){
            node = (TerraceNode*) root->neighbors[0]->node;
        }else{
            node = (TerraceNode*) root;
        }
        
        ASSERT(node);
    }
    
    int part;
    if(dad){
        TerraceNeighbor *nei = (TerraceNeighbor*)node->findNeighbor(dad);
        TerraceNeighbor *dad_nei = (TerraceNeighbor*)dad->findNeighbor(node);
        if(nei->link_neighbors.size()>0){
            //cout<<"| IF link_neighbours_exist -> clear them: size "<<nei->link_neighbors.size()<<"\n";
            
            // Clearing backward map from induced partition trees...
            if(clean_induced_part_maps){
                for(part=0; part<nei->link_neighbors.size(); part++){
                    // INFO: since for the trees with less than 3 taxa you do not do any maps, first check if there is a link_neighbor for the neighbor on the parent tree
                    if((TerraceNeighbor*)nei->link_neighbors[part]){
                        if(((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.size()>0){
                            ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors.clear();
                            ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors.clear();
                        }
                        if(((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors_lowtop_back.size()>0){
                            ((TerraceNeighbor*)nei->link_neighbors[part])->link_neighbors_lowtop_back.clear();
                            ((TerraceNeighbor*)dad_nei->link_neighbors[part])->link_neighbors_lowtop_back.clear();
                        }
                    
                    }
                }
            }
            
            // Clearing forward map from parent tree to induced partition trees...
            nei->link_neighbors.clear();
            dad_nei->link_neighbors.clear();
            
        }

        node->empty_br_dad_nei.clear();
        node->empty_br_node_nei.clear();
        
        dad->empty_br_dad_nei.clear();
        dad->empty_br_node_nei.clear();
    }
    
    FOR_NEIGHBOR_DECLARE(node, dad, it){
        cleanAllLinkINFO(clean_induced_part_maps, (TerraceNode*) (*it)->node, (TerraceNode*) node);
    }
    
}

TerraceNode* TerraceTree::insertNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, bool update_leafNode,bool update_brNodes){
    
    //cout<<"INSERTING: "<<node_name;
    
    TerraceNode *node_1, *node_2;
    
    node_1 = (TerraceNode*) newNode(nodeNum, node_name.c_str());
    leafNum += 1;
    nodeNum += 1;
    //assert(node_1->name == node_name);
    
    if(node_1_branch && node_2_branch){
        //cout<<"| case_3: br1-br2\n";
        node_2 = (TerraceNode*)(newNode(nodeNum));
        nodeNum += 1;
        
        int br_id = branchNum;
        branchNum += 2;
        
        node_1->addNeighbor(node_2, 0.0, br_id);
        node_2->addNeighbor(node_1, 0.0, br_id);
        
        if(update_brNodes){
            brNodes[br_id]={node_1,node_2};
        }
        
        br_id = node_1_branch->findNeighbor(node_2_branch)->id;
        node_1_branch->updateNeighbor(node_2_branch, node_2, 0.0);
        node_1_branch->findNeighbor(node_2)->id = br_id;
        node_2->addNeighbor(node_1_branch, 0.0, br_id);
        if(update_brNodes){
            brNodes[br_id]={node_1_branch,node_2};
        }
        
        br_id = branchNum-1;
        node_2_branch->updateNeighbor(node_1_branch, node_2, 0.0);
        node_2_branch->findNeighbor(node_2)->id = br_id;
        node_2->addNeighbor(node_2_branch, 0.0, br_id);
        if(update_brNodes){
            brNodes[br_id]={node_2_branch,node_2};
        }
    
    } else if(leafNum == 2){
        //cout<<"| case_2: root-leaf\n";
        root->addNeighbor(node_1, 0.0,0);
        node_1->addNeighbor(root, 0.0,0);
        branchNum = 1;
        if(update_brNodes){
            brNodes[0]={root,node_1};
        }
        //assert(nodeNum==2);
    } else {
        //cout<<"| case_1: root\n";
        root = node_1;
        //assert(nodeNum==1);
        //assert(branchNum=0);
    }
    
    //initializeTree();
    
    //assert(node_1->name == node_name);
    
    //drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    
    if(update_leafNode){
        leafNodes[node_name]=node_1;
    }
    return node_1;
    
}

void TerraceTree::remove_taxon(string taxon_name,bool update_leafNode,bool update_brNodes){

    //StrVector taxa;
    //taxa.push_back(taxon_name);
    
    //cout<<"REMOVING: "<<taxon_name;
    IntVector ids_brNodes;

    if(leafNum>=2){
        TerraceNode *node;
        if(update_leafNode){
            node = (TerraceNode*) leafNodes[taxon_name];
        }else{
            node = (TerraceNode*) findLeafName(taxon_name);
        }
        TerraceNeighbor* nei = (TerraceNeighbor*)node->neighbors[0];
        
        // nei is a central node (or a second leaf in a 2-taxon tree), the branch will be joined, therefore, the link_neis should be deleted (only the pointers, not objects).
        
        FOR_NEIGHBOR_DECLARE(nei->node,NULL, it){
            ((TerraceNeighbor*)(*it))->delete_ptr_members();
            ((TerraceNeighbor*)(*it)->node->findNeighbor(nei->node))->delete_ptr_members();
        }
        
        if(leafNum>2){
            //cout<<"| case_3: leafNUM>2\n";
            
            if (node == root)
            {    // find another root
                root = findFirstTaxon(root);
            }
            
            //    removeTaxa(taxa);
            //    initializeTree();
            
            Node *othernodes[2] = { nullptr, nullptr };
            int br_id=nei->id;
        
            FOR_NEIGHBOR(nei->node,node, it){
                if (othernodes[0] == nullptr){
                    othernodes[0] = (*it)->node;
                    if(update_brNodes){
                        ids_brNodes.push_back((*it)->id);
                    }
                }
                else if (othernodes[1] == nullptr){
                    othernodes[1] = (*it)->node;
                    if(update_brNodes){
                        if((*it)->id>ids_brNodes[0]){
                            ids_brNodes[0]=(*it)->id;
                        }
                    }
                }
                
                if((*it)->id<br_id){
                    br_id=(*it)->id;
                }
            }
        
            if(update_brNodes){
                ids_brNodes.push_back(nei->id);
                for(const auto& id: ids_brNodes){
                    brNodes[id]={nullptr,nullptr};
                    brNodes.erase(id);
                }
            }
            
            othernodes[0]->updateNeighbor(nei->node, othernodes[1], 0.0);
            othernodes[1]->updateNeighbor(nei->node, othernodes[0], 0.0);
            
            othernodes[0]->findNeighbor(othernodes[1])->id = br_id;
            othernodes[1]->findNeighbor(othernodes[0])->id = br_id;
            
            brNodes[br_id]={othernodes[0],othernodes[1]};
            
            delete nei->node;
            delete node;
            
            branchNum -=2;
            leafNum--;
            nodeNum -=2;
            
            //int branchNum_ = branchNum;
            //int leafNum_ = leafNum;
            //int nodeNum_ = nodeNum;
            
            
            //initializeTree(); // ???? why fail without this initialization???
            
            /*if(branchNum_ != branchNum or leafNum_ != leafNum or nodeNum_ != nodeNum){
                cout<<"-------------------------------------------------------\n";
                cout<<"BEFORE:"<<leafNum_<<"|"<<nodeNum_<<"|"<<branchNum_<<"\n";
                cout<<"AFTER:"<<leafNum<<"|"<<nodeNum<<"|"<<branchNum<<"\n";
                cout<<"-------------------------------------------------------\n";
            }*/
            
        } else {
            //cout<<"| case_2: leafNUM=2\n";
            //if(leafNum==2){
            //cout<<"two-taxon tree, remove one taxon"<<"\n";

            if(root->name == taxon_name){
                TerraceNode * new_root = (TerraceNode*) root->neighbors[0]->node;
                delete root;
                
                // Free pointers for the remaining node:
                new_root->deleteNode();
                
                // Set new root
                root = new_root;
                
            }else{
                delete root->neighbors[0]->node;
                
                // Free pointers for the remaining node:
                root->deleteNode();
            }
            
            leafNum = 1;
            nodeNum = 1;
            branchNum = 0;
            root->id = 0;
            
            if(update_brNodes){
                brNodes[0]={nullptr,nullptr};
                brNodes.erase(0);
            }
        }

    } else if(leafNum==1){
            //cout<<"| case_1: delete root\n";
            delete root;
            root = nullptr;
            leafNum = 0;
            nodeNum = 0;
            branchNum = 0;
    }
    
    if(update_leafNode){
        leafNodes[taxon_name]=nullptr;
        leafNodes.erase(taxon_name);
    }

}

void TerraceTree::print_terrace_tree(bool draw,ostream &out){
    
    if(leafNum>2 && draw){
        drawTree(out, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    }else if(leafNum==2 or (!draw && leafNum>1)){
        printTree(out, WT_BR_SCALE | WT_NEWLINE);
    }else if(leafNum==1){
        out<<"("<<root->name<<");"<<"\n";
    }else{
        out<<"();"<<"\n";
    }
    
}

string getTreeTopologyString(MTree* tree){
    stringstream tree_stream;
    tree->printTree(tree_stream, WT_BR_LEN_ROUNDING + WT_SORT_TAXA);
    return tree_stream.str();
}

void TerraceTree::fillLeafNodes(){
    
    NodeVector taxa_nodes;
    getTaxa(taxa_nodes);
    
    for(const auto& it: taxa_nodes){
        leafNodes[(*it).name]=it;
    }
    
    if(false){
        cout<<"-----------------------"<<"\n";
        for(const auto& entry: leafNodes){
            cout<<entry.first<<"->"<<entry.second->name<<"\n";
        }
    }
};

void TerraceTree::fillbrNodes(){
    //cout<<"Filling out branch nodes:\n";
    NodeVector br_1, br_2;
    getBranches(br_1, br_2);
    int id;
    for(int i=0; i<branchNum; i++){
        id=br_1[i]->findNeighbor(br_2[i])->id;
        brNodes[id]={br_1[i],br_2[i]};;
        /*for(const auto& b: brNodes[id]){
            cout<<"|Node_"<<b->id;
        }
        cout<<"|\n";*/
    }
};
