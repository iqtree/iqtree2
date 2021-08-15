#include "phylotree.h"

template <class Node> class WBTNode {
public:
    Node*  parent;
    Node*  leftChild;
    Node*  rightChild;
    double weight;
    WBTNode(): parent(nullptr), leftChild(nullptr), 
               rightChild(nullptr), weight(1.0) {
    }
    ~WBTNode() = default;
    inline Node* otherChild(Node* child) {
        if (leftChild==child) {
            return rightChild;
        } else if (rightChild==child) {
            return leftChild;
        }
        return nullptr;
    }
    inline void replaceChild(Node* out, Node* in) {
        leftChild  = (leftChild==out)  ? in : leftChild;
        rightChild = (rightChild==out) ? in : rightChild;
        in->parent = dynamic_cast<Node*>(this);
    }
    virtual void recalculate() {
        weight = (leftChild==nullptr ? 0 : leftChild->weight) 
               + (rightChild==nullptr ? 0 : rightChild->weight);
        weight = (weight==0.0 ? 1.0 : weight);
    }
};

template <class Node=WBTNode<Node> > class WeightBalancedTree {
protected:
    Node* root;
public:
    //Initialize and allocate
    WeightBalancedTree(): root(nullptr) {
    }
    WeightBalancedTree(const WeightBalancedTree& rhs) = delete;

    virtual Node* newNode(Node* left, Node* right) {
        Node* newNode = new Node();
        newNode->leftChild  = left;
        newNode->rightChild = right;
        newNode->recalculate();
        return newNode;
    }
protected:
    //Deallocate and destroy
    void deleteBelow(Node* top) {
        std::vector<Node*> stack;
        stack.push_back(top);
        do {
            top = stack.back();
            stack.pop_back();
            if (top->leftChild!=nullptr) {
                stack.push_back(top->leftChild);
            }
            if (top->rightChild!=nullptr) {
                stack.push_back(top->rightChild);
            }
            delete top;
        } 
        while (!stack.empty());
    }
public:
    ~WeightBalancedTree() {
        deleteBelow(root);
        root = nullptr;
    }
protected:
    Node* insertNextTo(Node* old, Node* added) {
        Node* parent      = old->parent;
        Node* replacement = newNode(old, added);
        old->parent       = replacement;
        added->parent     = replacement;
        if (parent!=nullptr) {
            //std::cout << "SIMA parent " << parent->taxon_id
            //          << " split link to " << old->taxon_id
            //          << " to accomodate " << added->taxon_id
            //          << "." << std::endl;
            parent->replaceChild(old, replacement);
            afterInsert(replacement);
        } else {
            //std::cout << "New root, above " << old->taxon_id
            //          << " and " << added->taxon_id
            //          << "." << std::endl;
            root = replacement;
        }
        return replacement;
    }
    void afterInsert(Node* node) {
        node->recalculate();
        if (false) {
            //Just recalculate, don't rebalance
            while (node->parent!=nullptr) {
                node = node->parent;
                node->recalculate();
            }
            return;
        }
        Node* parent = node->parent;
        if (parent==nullptr) {
            return;
        }
        Node* grandparent = parent->parent;
        if (grandparent==nullptr) {
            return;
        }
        do {
            if (grandparent->weight < node->weight + parent->weight) {
                //Rebalance
                Node* aunt   = grandparent->otherChild(parent);
                //Node* sister = parent->otherChild(node);
                //
                //FROM                    TO
                //====                    ====
                //     grandparent          grandparent
                //      /       \            /       \    ...
                //   parent     aunt      parent    node
                //  /     \                /  \  ...
                // node  sister        aunt  sister       
                //
                grandparent->replaceChild(aunt,  node);
                parent     ->replaceChild(node,  aunt);
                parent     ->recalculate();
                grandparent->recalculate();
                //
                //Note: If relative order had to be preserved, this would
                //      sometimes be more complicated (imagine sister on
                //      the left of node, in the above diagram). Like so:
                //
                //               granpdarent
                //                /        \        ...
                //            parent       node
                //             /   \      /    \    ... 
                //        sister  left right   aunt
                //                    
                //      where left and right were the subtrees of node.
                //                    
            } else {
                node = parent;
                node->recalculate();
            }
            parent      = grandparent;
            grandparent = grandparent->parent;
        }
        while (grandparent!=nullptr);
    }
};

class SimilarityNode: public WBTNode<SimilarityNode> {
public:
    typedef WBTNode<SimilarityNode> super;
    int          taxon_id;
    DoubleVector counts;

    virtual ~SimilarityNode() = default;

    virtual void recalculate() override {
        super::recalculate();
        if (leftChild==nullptr) {
            if (rightChild==nullptr) {
                return;
            } else {
                counts = rightChild->counts;
                return;
            }
        } else if (rightChild==nullptr) {
            counts = leftChild->counts;
            return;
        }
        const double* left   = leftChild->counts.data();
        const double* right  = rightChild->counts.data();
        intptr_t      states = leftChild->counts.size();

        counts.resize(states);
        double*       sum   = counts.data();
        for (intptr_t i=0; i<states; ++i) {
            sum[i] = left[i] + right[i];
        }
    }
    SimilarityNode* findPreferredNeighbor(SimilarityNode* top) {
        intptr_t states = counts.size();
        const double* here = counts.data();
        while (top->leftChild!=nullptr && 
               top->rightChild!=nullptr) {
            const double* left  = top->leftChild->counts.data();
            const double* right = top->rightChild->counts.data();
            double leftSum  = 0.0;
            double rightSum = 0.0;
            for (intptr_t i=0; i<states; ++i) {
                leftSum  += left[i]  * here[i];
                rightSum += right[i] * here[i];
            }
            double leftScore  = leftSum  / top->leftChild->weight  / (double)states;
            double rightScore = rightSum / top->rightChild->weight / (double)states;

            //std::cout << "SIMA leftScore=" << leftScore << "," 
            //          << " rightScore="    << rightScore << ","
            //          << " taxon ids "     << top->leftChild->taxon_id << ","
            //          << " and "           << top->rightChild->taxon_id
            //          << std::endl;

            if (leftScore<rightScore) {
                top = top->rightChild;
            } else if (rightScore<leftScore) {
                top = top->leftChild;
            } else if (top->leftChild->weight < top->rightChild->weight ) {
                top = top->rightChild;
            } else {
                top = top->leftChild;
            }
        }
        //std::cout << "SIMA found preferred neighbor." << std::endl;
        return top;
    }
    SimilarityNode* getLargerChild() const {
        if (leftChild==nullptr) {
            return rightChild;
        } else if (rightChild==nullptr) {
            return leftChild;
        }
        return (leftChild->weight < rightChild->weight) 
               ? rightChild : leftChild;
    }

};

class SimilarityTree: public WeightBalancedTree<SimilarityNode> {
public:

    SimilarityNode* insert(SimilarityNode* addMe, int parent_taxon_id) {
        if (root==nullptr) {
            root=addMe;
            return addMe;
        }
        SimilarityNode* neighbor  = addMe->findPreferredNeighbor(root);
        //std::cout << "SIMA preferred neighbor of " << addMe->taxon_id  
        //          << " has id " << neighbor->taxon_id << std::endl;

        SimilarityNode* parent = insertNextTo(neighbor, addMe);
        parent->taxon_id       = parent_taxon_id;
        while (root->parent!=nullptr) {
            root=root->parent;
        }
        return parent;
    }
    virtual SimilarityNode* newNode(SimilarityNode* left, SimilarityNode* right) override {
        SimilarityNode* node  = new SimilarityNode();
        node->taxon_id   = -1;
        node->leftChild  = left;
        node->rightChild = right;
        node->recalculate();
        return node;
    }
    SimilarityNode* getRoot() const {
        return root;
    }
};

void PhyloTree::calculateSimilarityWeights(std::vector<CharVector>& matrix,
                                           DoubleVector&            weights) {
    matrix.resize(aln->getNSeq());
    for (auto it = aln->begin(); it != aln->end(); ++it) {
        if (it->isInformative()) {
            auto     data = it->data();
            intptr_t size = it->size();
            for (int i=0; i<it->frequency; ++i) {
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for (intptr_t j=0; j<size; ++j) {
                    matrix[j].push_back(static_cast<char>(data[j]));
                }
                weights.push_back(1.0);
            }
            //std::cout << "SIMA weight " << it->frequency << std::endl;
        }
    }
}

void PhyloTree::constructSimilarityLeafNodes
        (SimilarityTree& sim_tree, const std::vector<CharVector>& matrix,
         const DoubleVector&   weights,
         std::vector<SimilarityNode*>& sim_nodes) {
    int      num_states = static_cast<char>(aln->num_states);
    intptr_t leaf_count = aln->getNSeq();
    sim_nodes.resize(leaf_count*2, nullptr);
    for (int taxon_id=0; taxon_id<leaf_count; ++taxon_id) {
        SimilarityNode* node = sim_tree.newNode(nullptr, nullptr);
        sim_nodes[taxon_id] = node;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int taxon_id=0; taxon_id<leaf_count; ++taxon_id) {
        SimilarityNode* node = sim_nodes[taxon_id];
        node->taxon_id       = taxon_id;
        const char* seq_data = matrix[taxon_id].data();
        intptr_t    seq_len  = matrix[taxon_id].size();
        DoubleVector& counts = node->counts;
        for (intptr_t pat=0; pat<seq_len; ++pat) {
            double w = weights[pat];
            int state_here = static_cast<int>(seq_data[pat]);
            for (int state = 0; state<num_states; ++state) {
                counts.push_back( (state==state_here) ? w : 0);
            }
        }
    }
}

void PhyloTree::constructSimilarityTree
        (std::vector<SimilarityNode*>& sim_nodes,
         SimilarityTree& sim_tree) {
    intptr_t leaf_count = aln->getNSeq();
    for (intptr_t taxon_id=0; taxon_id<leaf_count; ++taxon_id) {
        intptr_t        interior_node_id = leaf_count+taxon_id;
        SimilarityNode* taxon_node       = sim_nodes[taxon_id];
        SimilarityNode* interior_node    = sim_tree.insert(taxon_node, static_cast<int>(interior_node_id));
        sim_nodes[interior_node_id]      = interior_node;
    }
}

void PhyloTree::mirrorSimilarityTreeStructure
        (const std::vector<SimilarityNode*>& sim_nodes,
         std::vector<PhyloNode*>&      tree_nodes) {
    intptr_t leaf_count = aln->getNSeq();
    for (intptr_t node_id=0; node_id<leaf_count; ++node_id) {
        tree_nodes.push_back(newNode(static_cast<int>(node_id), aln->getSeqName(node_id).c_str()));
    }
    tree_nodes.push_back(nullptr);
    //sim_nodes[leaf_count] is *not* an interior node, it's an exterior.
    for (intptr_t node_id=leaf_count+1; 
         node_id<leaf_count+leaf_count; ++node_id) {
        tree_nodes.push_back(newNode(static_cast<int>(node_id)));
    }
    for (intptr_t node_id=leaf_count+1; 
         node_id<leaf_count+leaf_count; ++node_id) {
        PhyloNode*      dest   = tree_nodes[node_id];
        SimilarityNode* source = sim_nodes [node_id];
        PhyloNode*      left   = tree_nodes[source->leftChild->taxon_id];
        PhyloNode*      right  = tree_nodes[source->rightChild->taxon_id];

        ASSERT(dest!=left);
        ASSERT(dest!=right);
        ASSERT(left!=right);

        dest ->addNeighbor(left,  -1);
        dest ->addNeighbor(right, -1);
        left ->addNeighbor(dest,  -1);
        right->addNeighbor(dest,  -1);
        //std::cout << "SIMA linked children " << left->id << " and " << right->id
        //          << " to their parent " << dest->id << std::endl;
    }
}

void PhyloTree::correctSimilarityTreeStructure
        (const SimilarityTree& sim_tree,
         const std::vector<PhyloNode*>& tree_nodes) {
    //5. Need to correct tree structure.
    //
    //From this           To this
    //=========           =======
    //         centre        b
    //          /  \         |
    //         a    b      centre
    //        / \           /  \
    //      a1   a2       a1    a2
    //
    SimilarityNode* simroot      = sim_tree.getRoot();
    int             root_node_id = simroot->taxon_id;
    PhyloNode*      centre       = tree_nodes[root_node_id];
    int             a_node_id    = simroot->getLargerChild()->taxon_id;
    PhyloNode*      a            = tree_nodes[a_node_id];

    std::vector<PhyloNode*> aChildren;
    FOR_EACH_ADJACENT_PHYLO_NODE(a, centre, it, aChild) {
        aChildren.push_back(aChild);
    }
    ASSERT(centre->unlinkNeighbor(a));
    ASSERT(a->unlinkNeighbor(centre));
    for (PhyloNode* aChild : aChildren) {
        centre->addNeighbor(aChild, -1);
        aChild->updateNeighbor(a, centre);
    }
    delete a;
    root = tree_nodes[0];
}

int PhyloTree::createSimilarityMatchingTree() {

    deleteAllPartialParsimony();

    SimilarityTree               sim_tree;
    std::vector<SimilarityNode*> sim_nodes;
    {
        std::vector<CharVector>  matrix;
        DoubleVector             weights;
        calculateSimilarityWeights  (matrix, weights);
        constructSimilarityLeafNodes(sim_tree, matrix, weights, sim_nodes);
    }
    constructSimilarityTree         (sim_nodes,sim_tree);

    {
        std::vector<PhyloNode*> tree_nodes;
        mirrorSimilarityTreeStructure (sim_nodes, tree_nodes);
        correctSimilarityTreeStructure(sim_tree, tree_nodes);
        sim_nodes.clear();
    }

    initializeTree();
    initializeAllPartialPars();

    /* how long does this take?! */
    double parsimony_score = computeParsimony("Computing WBTree parsimony",
                                              true, true);
    setAllBranchLengthsFromParsimony(false, parsimony_score);
    deleteAllPartialParsimony();
    
    return static_cast<int>(floor(parsimony_score));
}
