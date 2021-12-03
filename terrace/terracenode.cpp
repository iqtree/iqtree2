/*
 *  terracenode.cpp
 *
 *  Created on: Sep 9, 2020
 *      Author: Olga
 */

#include "terracenode.hpp"

TerraceNode::TerraceNode()
: Node()
{
    init();
}


TerraceNode::TerraceNode(int aid) : Node(aid)
{
    init();
}

TerraceNode::TerraceNode(int aid, int aname) : Node (aid, aname) {
    init();
}


TerraceNode::TerraceNode(int aid, const char *aname) : Node(aid, aname) {
    init();
}

void TerraceNode::init() {

}


void TerraceNode::addNeighbor(Node *node, double length, int id) {
    neighbors.push_back(new TerraceNeighbor(node, length, id));
}

void TerraceNode::deleteNode(){
    
    NeighborVec::reverse_iterator it;
    for (it = empty_br_dad_nei.rbegin(); it != empty_br_dad_nei.rend(); it++)
        delete (*it);
    empty_br_dad_nei.clear();
    
    for (it = empty_br_node_nei.rbegin(); it != empty_br_node_nei.rend(); it++)
        delete (*it);
    empty_br_node_nei.clear();
    
    for (it = neighbors.rbegin(); it != neighbors.rend(); it++)
        delete (*it);
    neighbors.clear();

}

TerraceNode::~TerraceNode(){
    deleteNode();
}

TerraceNeighbor::~TerraceNeighbor(){
    delete_ptr_members();
}

void TerraceNeighbor::delete_ptr_members(){
    
    //for(Neighbor* nei: link_neighbors)
    //    delete nei;
    link_neighbors.clear();
    
    //for(Neighbor* nei: link_neighbors_lowtop_back)
    //    delete nei;
    link_neighbors_lowtop_back.clear();
    
    //for(Node* node: taxa_to_insert)
    //    delete node;
    taxa_to_insert.clear();
    
}

void TerraceNeighbor::printInfo(Node *dad){
    
    cout<<"TerraceNeighbour information for node "<<this->node->id<<", a neighbour of node "<<dad->id<<":"<<"\n";
    
    NeighborVec::iterator it1;
    vector<Node*>::iterator it2;
    
    int link_neighbors_size = link_neighbors.size();
    int taxa_to_insert_size = taxa_to_insert.size();
    
    int i=0;
    cout<<"There are "<<link_neighbors_size<<" link neighbors."<<"\n";
    if(link_neighbors_size>0){
        for(it1=link_neighbors.begin(); it1<link_neighbors.end(); it1++){
            i++;
            cout<<"link_neighbor["<<i<<"]="<<(*it1)->node->id<<"\n";
        }
        cout<<"\n";
    }
    i=0;
    cout<<"There are "<<taxa_to_insert_size<<" taxa to insert."<<"\n";
    if(taxa_to_insert_size>0){
        for(it2=taxa_to_insert.begin(); it2<taxa_to_insert.end(); it2++){
            i++;
            cout<<"taxa_to_insert["<<i<<"]="<<(*it2)->name<<"\n";
        }
        cout<<"\n";
    }
    
}
