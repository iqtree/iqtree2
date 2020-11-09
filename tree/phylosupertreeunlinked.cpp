//
//  phylosupertreeunlinked.cpp
//  tree
//
//  Created by Minh Bui on 2/5/18.
//

#include "phylosupertreeunlinked.h"
#include "utils/MPIHelper.h"
#include "utils/timeutil.h"

extern ostream cmust;

PhyloSuperTreeUnlinked::PhyloSuperTreeUnlinked(SuperAlignment *alignment)
: PhyloSuperTree(alignment, true)
{
}


void PhyloSuperTreeUnlinked::readTree(istream &in, bool &is_rooted) {
    for (iterator it = begin(); it != end(); it++) {
        (*it)->rooted = Params::getInstance().is_rooted;
        (*it)->readTree(in, (*it)->rooted);
        is_rooted |= (*it)->rooted;
    }
}

void PhyloSuperTreeUnlinked::setAlignment(Alignment *alignment) {
    ASSERT(alignment->isSuperAlignment());
    SuperAlignment *saln = (SuperAlignment*)alignment;
    ASSERT(saln->partitions.size() == size());
    for (int i = 0; i < size(); i++)
        at(i)->setAlignment(saln->partitions[i]);
}

/**
 * setup all necessary parameters  (declared as virtual needed for phylosupertree)
 */
void PhyloSuperTreeUnlinked::initSettings(Params& params) {
    PhyloSuperTree::initSettings(params);
    for (auto it = begin(); it != end(); it++)
        ((IQTree*)(*it))->initSettings(params);
}

void PhyloSuperTreeUnlinked::mapTrees() {
    // do nothing here as partition trees are unlinked
}

int PhyloSuperTreeUnlinked::computeParsimonyTree(const char *out_prefix, Alignment *alignment, int *rand_stream) {
    SuperAlignment *saln = (SuperAlignment*)alignment;
    int score = 0;
    int i;
    ASSERT(saln->partitions.size() == size());
    for (i = 0; i < size(); i++) {
        score += at(i)->computeParsimonyTree(NULL, saln->partitions[i], rand_stream);
    }
    if (out_prefix) {
        string file_name = out_prefix;
        file_name += ".parstree";
        try {
            ofstream out;
            out.open(file_name.c_str());
            for (i = 0; i < size(); i++) {
                at(i)->printTree(out, WT_NEWLINE);
            }
            out.close();
        } catch (...) {
            outError("Cannot write to file ", file_name);
        }
    }
    return score;
}

int PhyloSuperTreeUnlinked::wrapperFixNegativeBranch(bool force_change) {
    // Initialize branch lengths for the parsimony tree
    int numFixed = 0;
    for (auto tree = begin(); tree != end(); tree++) {
        numFixed += (*tree)->fixNegativeBranch(force_change);
        (*tree)->resetCurScore();
    }
    return numFixed;
}

bool PhyloSuperTreeUnlinked::isBifurcating(Node *node, Node *dad) {
    for (auto it = begin(); it != end(); it++)
        if (!(*it)->isBifurcating())
            return false;
    return true;
}

string PhyloSuperTreeUnlinked::getTreeString() {
    stringstream tree_stream;
    for (iterator it = begin(); it != end(); it++)
        (*it)->printTree(tree_stream, WT_TAXON_ID + WT_BR_LEN + WT_SORT_TAXA);
    return tree_stream.str();
}

void PhyloSuperTreeUnlinked::readTreeString(const string &tree_string) {
    stringstream str;
    str << tree_string;
    str.seekg(0, ios::beg);
    for (iterator it = begin(); it != end(); it++) {
        (*it)->freeNode();
        (*it)->readTree(str, rooted);
        (*it)->assignLeafNames();
        (*it)->resetCurScore();
    }
}

void PhyloSuperTreeUnlinked::saveCheckpoint() {
    for (iterator it = begin(); it != end(); it++) {
        checkpoint->startStruct((*it)->aln->name);
        (*it)->saveCheckpoint();
        checkpoint->endStruct();
    }
}

void PhyloSuperTreeUnlinked::restoreCheckpoint() {
    for (iterator it = begin(); it != end(); it++) {
        checkpoint->startStruct((*it)->aln->name);
        (*it)->restoreCheckpoint();
        checkpoint->endStruct();
    }
}

/**
 * save branch lengths into a vector
 */
void PhyloSuperTreeUnlinked::saveBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad) {
    int totalBranchNum = 0;
    for (auto it = begin(); it != end(); it++) {
        totalBranchNum += (*it)->branchNum * (*it)->getMixlen();
    }
    lenvec.resize(startid + totalBranchNum);
    
    for (auto  it = begin(); it != end(); it++) {
        (*it)->saveBranchLengths(lenvec, startid);
        startid += (*it)->branchNum * (*it)->getMixlen();
    }
}
/**
 * restore branch lengths from a vector previously called with saveBranchLengths
 */
void PhyloSuperTreeUnlinked::restoreBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad) {
    for (iterator it = begin(); it != end(); it++) {
        (*it)->restoreBranchLengths(lenvec, startid);
        startid += (*it)->branchNum * (*it)->getMixlen();
    }
}

void PhyloSuperTreeUnlinked::setRootNode(const char *my_root, bool multi_taxa) {
    // DOES NOTHING
}

void PhyloSuperTreeUnlinked::computeBranchLengths() {
    // DOES NOTHING
}

void PhyloSuperTreeUnlinked::printTree(ostream &out, int brtype) {
    for (iterator tree = begin(); tree != end(); tree++)
        (*tree)->printTree(out, brtype);
}

void PhyloSuperTreeUnlinked::printResultTree(string suffix) {
    if (MPIHelper::getInstance().isWorker()) {
        return;
    }
    if (params->suppress_output_flags & OUT_TREEFILE)
        return;
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    if (suffix.compare("") != 0) {
        tree_file_name += "." + suffix;
    }
    ofstream out;
    out.open(tree_file_name.c_str());
    for (iterator tree = begin(); tree != end(); tree++)
        (*tree)->printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    out.close();
    if (verbose_mode >= VB_MED)
        cout << "Best tree printed to " << tree_file_name << endl;
}

double PhyloSuperTreeUnlinked::treeLength(Node *node, Node *dad) {
    double len = 0.0;
    for (iterator tree = begin(); tree != end(); tree++)
        len += (*tree)->treeLength();
    return len;
}

double PhyloSuperTreeUnlinked::treeLengthInternal( double epsilon, Node *node, Node *dad) {
    double len = 0.0;
    for (iterator tree = begin(); tree != end(); tree++)
        len += (*tree)->treeLengthInternal(epsilon);
    return len;
}

pair<int, int> PhyloSuperTreeUnlinked::doNNISearch(bool write_info, const char* context) {
    int NNIs = 0, NNI_steps = 0;
    double score = 0.0;
#pragma omp parallel for schedule(dynamic) num_threads(num_threads) if (num_threads > 1) reduction(+: NNIs, NNI_steps, score)
    for (int i = 0; i < size(); i++) {
        IQTree *part_tree = (IQTree*)at(part_order[i]);
        Checkpoint *ckp = new Checkpoint;
        getCheckpoint()->getSubCheckpoint(ckp, part_tree->aln->name);
        part_tree->setCheckpoint(ckp);
        auto num_NNIs = part_tree->doNNISearch(false, context);
        NNIs += num_NNIs.first;
        NNI_steps += num_NNIs.second;
        score += part_tree->getCurScore();
#pragma omp critical
        {
        getCheckpoint()->putSubCheckpoint(ckp, part_tree->aln->name);
        getCheckpoint()->dump();
        }
        delete ckp;
        part_tree->setCheckpoint(getCheckpoint());
    }

    setCurScore(score);
    cout << "Log-likelihood: " << score << endl;
    return std::make_pair(NNIs, NNI_steps);
}

double PhyloSuperTreeUnlinked::doTreeSearch() {
    double tree_lh = 0.0;
    string bestTree;
    
    cout << "--------------------------------------------------------------------" << endl;
    cout << "|                SEPARATE TREE SEARCH FOR PARTITIONS               |" << endl;
    cout << "--------------------------------------------------------------------" << endl;

    if (part_order.empty())
        computePartitionOrder();

    int saved_flag = params->suppress_output_flags;
    params->suppress_output_flags |= OUT_TREEFILE + OUT_LOG;
    VerboseMode saved_mode = verbose_mode;
    verbose_mode = VB_QUIET;
    bool saved_print_ufboot_trees = params->print_ufboot_trees;
    params->print_ufboot_trees = false;

#pragma omp parallel for schedule(dynamic) num_threads(num_threads) if (num_threads > 1) reduction(+: tree_lh)
    for (int i = 0; i < size(); i++) {
        IQTree *part_tree = (IQTree*)at(part_order[i]);
        Checkpoint *ckp = new Checkpoint;
        getCheckpoint()->getSubCheckpoint(ckp, part_tree->aln->name);
        part_tree->setCheckpoint(ckp);
        double score = part_tree->doTreeSearch();
        tree_lh += score;
#pragma omp critical
        {
            getCheckpoint()->putSubCheckpoint(ckp, part_tree->aln->name);
            getCheckpoint()->dump();
            cmust << "Partition " << part_tree->aln->name
                 << " / Iterations: " << part_tree->stop_rule.getCurIt()
                 << " / LogL: " << score
                 << " / Time: " << convert_time(getRealTime() - params->start_real_time)
                 << endl;
        }
        delete ckp;
        part_tree->setCheckpoint(getCheckpoint());
    }

    verbose_mode = saved_mode;
    params->suppress_output_flags= saved_flag;
    params->print_ufboot_trees = saved_print_ufboot_trees;

    if (tree_lh < curScore)
        cout << "BETTER TREE FOUND: " << tree_lh << endl;
    curScore = tree_lh;

    bestTree = getTreeString();
    addTreeToCandidateSet(bestTree, getCurScore(), false, MPIHelper::getInstance().getProcessID());
    printResultTree();
    intermediateTrees.update(bestTree, getCurScore());
    candidateTrees.saveCheckpoint();
    return curScore;
}

void PhyloSuperTreeUnlinked::summarizeBootstrap(Params &params) {
    for (auto tree = begin(); tree != end(); tree++)
        ((IQTree*)*tree)->summarizeBootstrap(params);
}

void PhyloSuperTreeUnlinked::writeUFBootTrees(Params &params) {
    //    IntVector tree_weights;
    int i, j;
    string filename = params.out_prefix;
    filename += ".ufboot";
    ofstream out(filename.c_str());
    
    for (auto tree = begin(); tree != end(); tree++) {
        MTreeSet trees;

        trees.init(((IQTree*)*tree)->boot_trees, (*tree)->rooted);
        for (i = 0; i < trees.size(); i++) {
            NodeVector taxa;
            // change the taxa name from ID to real name
            trees[i]->getOrderedTaxa(taxa);
            for (j = 0; j < taxa.size(); j++)
                taxa[j]->name = aln->getSeqName(taxa[j]->id);
            if (removed_seqs.size() > 0) {
                // reinsert removed seqs into each tree
                trees[i]->insertTaxa(removed_seqs, twin_seqs);
            }
            // now print to file
            for (j = 0; j < trees.tree_weights[i]; j++)
                if (params.print_ufboot_trees == 1)
                    trees[i]->printTree(out, WT_NEWLINE);
                else
                    trees[i]->printTree(out, WT_NEWLINE + WT_BR_LEN);
        }
    }
    cout << "UFBoot trees printed to " << filename << endl;
    out.close();
}

/**
 Test all branches of the tree with aLRT SH-like interpretation
 */
int PhyloSuperTreeUnlinked::testAllBranches(int threshold, double best_score, double *pattern_lh,
                            int reps, int lbp_reps, bool aLRT_test, bool aBayes_test,
                            PhyloNode *node, PhyloNode *dad)
{
    int num_low_support = 0;
    double *ptn_lh[size()];
    ptn_lh[0] = pattern_lh;
    for (int id = 1; id < size(); id++) {
        ptn_lh[id] = ptn_lh[id-1] + at(id-1)->getAlnNPattern();
    }
#ifdef _OPENMP
#pragma omp parallel for reduction(+: num_low_support)
#endif
    for (int id = 0; id < size(); id++) {
        num_low_support += at(id)->testAllBranches(threshold, at(id)->getCurScore(), ptn_lh[id],
                            reps, lbp_reps, aLRT_test, aBayes_test);
    }
    return num_low_support;
}

int PhyloSuperTreeUnlinked::testNumThreads() {
#ifdef _OPENMP
    // unlinked partitions scales well with many cores
    int bestProc = min(countPhysicalCPUCores(), params->num_threads_max);
    bestProc = min(bestProc, (int)size());
    cout << "BEST NUMBER OF THREADS: " << bestProc << endl << endl;
    setNumThreads(bestProc);
    return bestProc;
#else
    return 1;
#endif
}
