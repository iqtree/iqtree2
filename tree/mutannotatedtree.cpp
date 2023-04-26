#include "mutannotatedtree.h"

const string MutAnnotatedTree::CONDENSED_NAME = "MAT_condensed";
const string MutAnnotatedTree::SEPARATOR = "@@@";

void MutAnnotatedTree::initMAT(Node* root, const int num_threads)
{
    // seek the ref node
    Node* ref_node = getRefNode(root);
    
    ref_node_id = ref_node->id;
    ref_name = ref_node->name;
    ref_seq_chunks.resize(num_threads);
    
    // clone the root sequence
    root_seq_chunks = root->sequence->sequence_chunks;
}

int MutAnnotatedTree::getRefNodeID()
{
    return ref_node_id;
}

string MutAnnotatedTree::getRefName()
{
    return ref_name;
}

vector<vector<short int>>& MutAnnotatedTree::getRefSeqChunks()
{
    return ref_seq_chunks;
}

vector<vector<short int>>& MutAnnotatedTree::getRootSeqChunks()
{
    return root_seq_chunks;
}

void MutAnnotatedTree::addNodeData2MAT(Parsimony::data& data, Node *node, Node *dad)
{
    // ignore condensed leave (except the first one)
    if (node->isLeaf() && dad->name.find(CONDENSED_NAME) != std::string::npos)
    {
        // split node_name into token regarding SEPARATOR
        std::vector<std::string> tokens = splitString(dad->name, SEPARATOR);
        
        // CONDENSED_NAME@@@<internal_real_name>@@@<condense_count>@@@<first_ident_leaf>@@@<other_ident_leave>
        if (std::find(tokens.begin() + 4, tokens.end(), node->name) != tokens.end())
            return;
    }
        
    // add an empty metadata
    auto meta = data.add_metadata();
    
    // debug: MAT
    // std::cout<<node->name<<std::endl;
    
    Parsimony::mutation_list* node_mutations = data.add_node_mutations();
    for (auto mutation_list:node->sequence->node_mutations_vec)
        node_mutations->MergeFrom(mutation_list);
    
    // release memory
    vector<Parsimony::mutation_list>().swap(node->sequence->node_mutations_vec);
    
    /*for (auto mutation:node->sequence->node_mutations_vec[0].mutation())
    {
        std::string mut_nuc = "";
        for (auto nuc:mutation.mut_nuc())
            mut_nuc += std::to_string(nuc) + " ";

        std::cout << mutation.position() << " " << mutation.ref_nuc()<< " " << mutation.par_nuc()<< " " << mut_nuc << mutation.chromosome() << std::endl;
    }
    
    for (auto mutation:data.node_mutations(data.node_mutations_size() - 1).mutation())
    {
        std::string mut_nuc = "";
        for (auto nuc:mutation.mut_nuc())
            mut_nuc += std::to_string(nuc) + " ";
        std::cout << mutation.position() << " " << mutation.ref_nuc()<< " " << mutation.par_nuc()<< " " << mut_nuc << mutation.chromosome() << std::endl;
    }*/

    // Add condensed nodes
    /*for (auto cn: tree.condensed_nodes) {
        auto cn_ptr = data.add_condensed_nodes();
        cn_ptr->set_node_name(cn.first);
        for (auto lid: cn.second) {
            cn_ptr->add_condensed_leaves(lid);
        }
    }*/
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        addNodeData2MAT(data, (*it)->node, node);
    }
}

void MutAnnotatedTree::saveMAT(IQTree* tree, const std::string& ref_seq, const std::string& filename) {
    // output the ref sequence
    std::ofstream reffile(filename + ".ref.fa", std::ios::out);
    reffile << ">" << ref_name << "\n";
    reffile << ref_seq << "\n";
    reffile.close();
    
    Parsimony::data data;
    
    // condense identical leave and set the newick string
    condenseIdentLeave(tree, data);
    
    Node* start_node = tree->root;
    // handle the special case when we root an unrooted tree
    if (start_node->neighbors.size() == 1 && start_node->neighbors[0]->node->neighbors.size() > 2)
    {
        // move the mutations from the root the its child
        start_node->neighbors[0]->node->sequence->node_mutations_vec = std::move(start_node->sequence->node_mutations_vec);
        
        // start traversing the tree from the child
        start_node = start_node->neighbors[0]->node;
    }
    addNodeData2MAT(data, start_node, tree->root);

    /*auto dfs = tree.depth_first_expansion();
    
    std::cout << data.newick() << std::endl;
    std::cout << "metadata_size: " << data.metadata_size() << std::endl;
    for (auto metadata:data.metadata())
        std::cout << metadata.clade_annotations_size() << std::endl;
    std::cout << "condensed_nodes_size: " << data.condensed_nodes_size() << std::endl;
    for (auto condensed_node:data.condensed_nodes())
    {
        std::cout << condensed_node.node_name() << std::endl;
    }
    std::cout << "node_mutations_size: " << data.node_mutations_size() << std::endl;
    for (auto mutations:data.node_mutations())
    {
        std::cout << "----" << std::endl;
        for (auto mutation:mutations.mutation())
        {
            std::string mut_nuc = "";
            for (auto nuc:mutation.mut_nuc())
                mut_nuc += std::to_string(nuc) + " ";

            std::cout << mutation.position() << " " << mutation.ref_nuc()<< " " << mutation.par_nuc()<< " " << mut_nuc << mutation.chromosome() << std::endl;
        }
    }*/

    // write data to file
    std::ofstream outfile(filename + ".pd", std::ios::out | std::ios::binary);
    data.SerializeToOstream(&outfile);
    outfile.close();
    
    // show the output file name
    cout << "A mutation-annotated tree has been exported to " << filename + ".pd"  << endl;
    cout << "The reference sequence has been exported to " << filename + ".ref.fa"  << endl;
}

void MutAnnotatedTree::readMAT(const std::string& filename)
{
    Parsimony::data data;
    
    // read from file
    std::ifstream instream(filename, std::ios::in | std::ios::binary);
    google::protobuf::io::IstreamInputStream stream(&instream);
    google::protobuf::io::CodedInputStream input(&stream);
    //input.SetTotalBytesLimit(BIG_SIZE, BIG_SIZE);
    data.ParseFromCodedStream(&input);
    
    std::cout << "- NEWICK TREE:" << std::endl;
    std::cout << "  " <<  data.newick() << std::endl;
    std::cout << "- METADATA: " << data.metadata_size() << " element(s)" << std::endl;
    for (auto metadata:data.metadata())
        std::cout << "  " << metadata.clade_annotations_size() << std::endl;
    std::cout << "- CONDENSED NODES: " << data.condensed_nodes_size() << " element(s)" << std::endl;
    for (auto condensed_node:data.condensed_nodes())
    {
        std::cout << "  + " << condensed_node.node_name() << std::endl;
        for (auto condensed_leave:condensed_node.condensed_leaves())
            std::cout << "    ++ "<< condensed_leave << std::endl;
    }
    std::cout << "- NODE MUTATIONS: " << data.node_mutations_size() << " element(s)" << std::endl;
    std::cout << "  <position> <ref_nuc> <par_nuc> <mut_nuc> <chromosome>" << std::endl;
    for (auto mutations:data.node_mutations())
    {
        std::cout << "  ----------------" << std::endl;
        for (auto mutation:mutations.mutation())
        {
            std::string mut_nuc = "";
            for (auto nuc:mutation.mut_nuc())
                mut_nuc += std::to_string(nuc) + " ";

            std::cout << "  " << mutation.position() << " " << mutation.ref_nuc()<< " " << mutation.par_nuc()<< " " << mut_nuc << mutation.chromosome() << std::endl;
        }
    }
}

int MutAnnotatedTree::getRefState(int thread_id, int position)
{
    return ref_seq_chunks[thread_id][position];
}

void MutAnnotatedTree::updateRefState(int thread_id, int segment_start, Node *node, Node *dad)
{
    // update ref state for the mutations of the current node
    Parsimony::mutation_list& mutations = node->sequence->node_mutations_vec[thread_id];
    for (auto i = 0; i < mutations.mutation_size(); ++i)
    {
        Parsimony::mut* mutation = mutations.mutable_mutation(i);
        mutation->set_ref_nuc(getRefState(thread_id, mutation->position() - segment_start));
        mutation->set_chromosome(ref_name);
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        updateRefState(thread_id, segment_start, (*it)->node, node);
    }
}

Node* MutAnnotatedTree::getRefNode(Node* node, Node* dad)
{    
    // traverse the tree
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // return the first leaf
        if ((*it)->node->isLeaf())
            return (*it)->node;
        
        // browse 1-step deeper to the neighbor node
        return getRefNode((*it)->node, node);
    }
}

void MutAnnotatedTree::condenseIdentLeave(IQTree* tree, Parsimony::data& data)
{
    // traverse the tree to condense identical leave
    int condense_count = 0;
    condenseIdentLeaveAtInternal(data, condense_count, tree->root);
    
    // set the newick string
    stringstream nwk;
    // if condense_count > 0 -> clone the current tree, change and output the newick string
    if (condense_count > 0)
    {
        // clone the tree
        IQTree new_tree;
        new_tree.copyTree(tree);
        
        // update the tree -> collapse identical leave
        collapseIdentLeave(new_tree.root);
        
        // output the newick string
        new_tree.printTree(nwk);
    }
    // otherwise, output the newick string from the current tree
    else
    {
        tree->printTree(nwk);
    }
    data.set_newick(nwk.str());
}

void MutAnnotatedTree::condenseIdentLeaveAtInternal(Parsimony::data& data, int& condense_count, Node* node, Node* dad)
{
    std::vector<std::string> condensed_leaves;
    
    // browse all children of this node
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // child is a leaf -> check if no mutation occurs
        if ((*it)->node->isLeaf())
        {
            bool no_mutation = true;
            std::vector<Parsimony::mutation_list>& node_mutations_vec = (*it)->node->sequence->node_mutations_vec;
            // browse all chunks to check whether there is any mutations occurs at any chunk
            for (auto i = 0; i < node_mutations_vec.size(); ++i)
            {
                // stop checking if mutations occurs at any chunk
                if (node_mutations_vec[i].mutation_size() > 0)
                {
                    no_mutation = false;
                    break;
                }
            }
            
            // if no mutation occurs -> record this leave
            if (no_mutation)
                condensed_leaves.push_back((*it)->node->name);
        }
        // child is an internal node -> // browse 1-step deeper to the neighbor node
        else
            condenseIdentLeaveAtInternal(data, condense_count, (*it)->node, node);
    }
    
    // there are at least two leave in the list 'condensed_leaves' -> condense these leave
    if (condensed_leaves.size() > 1)
    {
        // add a new condense node
        Parsimony::condensed_node* condensed_node = data.add_condensed_nodes();
        ++condense_count;
        string condense_count_str = convertIntToString(condense_count);
        condensed_node->set_node_name("node_" + condense_count_str + "_condensed_" + convertIntToString(condensed_leaves.size()) + "_leaves");
        // CONDENSED_NAME@@@<internal_real_name>@@@<condense_count>@@@<first_ident_leaf>@@@<other_ident_leave>
        string tmp_name = CONDENSED_NAME + SEPARATOR + node->name + SEPARATOR + condense_count_str;
        for (auto i = 0; i < condensed_leaves.size(); ++i)
        {
            tmp_name += SEPARATOR + condensed_leaves[i];
            condensed_node-> add_condensed_leaves(std::move(condensed_leaves[i]));
        }
        // temporarily change the name of the internal node
        node->name = tmp_name;
    }
}

void MutAnnotatedTree::collapseIdentLeave(Node* node, Node* dad)
{
    // if node is an iternal node and contains idential leaves -> collapse them
    if (!node->isLeaf() && node->name.find(CONDENSED_NAME) != std::string::npos)
    {
        // split node_name into token regarding SEPARATOR
        string node_name = node->name;
        std::vector<std::string> tokens = splitString(node_name, SEPARATOR);
        
        int num_condensed_leave = tokens.size() - 3; // not count CONDENSED_NAME, <real_internal_node_name>, and condense_count
        ASSERT(num_condensed_leave <= node->neighbors.size() - 1);
        // recover the real name of the internal node
        node->name = tokens[1];
        // collapse some of children of this internal node
        node_name = "node_" + tokens[2] + "_condensed_" + convertIntToString(num_condensed_leave) + "_leaves";
        // delete some of identical children (except the first one)
        NeighborVec::iterator it;
        FOR_NEIGHBOR(node, dad, it) {
            // re-use the first identical child as the condensed node
            if ((*it)->node->name == tokens[3])
                (*it)->node->name = node_name;
            else
            {
                // delete other identical children
                if (std::find(tokens.begin() + 4, tokens.end(), (*it)->node->name) != tokens.end())
                {
                    delete (*it)->node;
                    node->neighbors.erase(it);
                    --it;
                }
            }
        }
    }
    
    // traverse the tree
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        collapseIdentLeave((*it)->node, node);
    }
}
