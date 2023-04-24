#ifndef MUTANNOTATEDTREE_H
#define MUTANNOTATEDTREE_H

#include "iqtree.h"
#include "proto/parsimony.pb.h"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

/**
 Mutation Annotated Tree class
 */
class MutAnnotatedTree{
    /**
     *  id of the reference node
     */
    int ref_node_id;
    
    /**
     *  name of the reference node
     */
    string ref_name;
    
    /**
     *  reference sequence
     */
   vector<vector<short int>> ref_seq_chunks;
    
    /**
     *  root sequence
     */
   vector<vector<short int>> root_seq_chunks;
    
    /**
    *  iteratively data at each node to MAT
    */
    void addNodeData2MAT(Parsimony::data& data, Node *node, Node *dad);
    
    /**
     *  get the ref state
     */
    int getRefState(int thread_id, int position);
    
    /**
        get the reference node (the first leaf in a pre-order traversal)
    */
    Node* getRefNode(Node* node, Node* dad = NULL);

public:
    /**
     *  init mat
     */
    void initMAT(Node* root, const int num_threads);
    
    /**
     *  get the id of the ref node
     */
    int getRefNodeID();
    
    /**
     *  get ref name
     */
    string getRefName();
    
    /**
     *  get chunks of reference sequence
     */
    vector<vector<short int>>& getRefSeqChunks();
    
    /**
     *  get chunks of root sequence
     */
    vector<vector<short int>>& getRootSeqChunks();
    
    /**
     *  save MAT to a file
     */
    void saveMAT(IQTree* tree, const std::string& ref_seq, const std::string& filename);

    /**
     *  read MAT from a file
     */
    void readMAT(const std::string& filename);
    
    /**
    *  traverse the tree update the ref state for all mutations
    */
    void updateRefState(int thread_id, int segment_start, Node *node, Node *dad);
};

#endif
