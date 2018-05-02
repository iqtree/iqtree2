//
//  phylosupertreeunlinked.cpp
//  tree
//
//  Created by Minh Bui on 2/5/18.
//

#include "phylosupertreeunlinked.h"
#include "utils/MPIHelper.h"

PhyloSuperTreeUnlinked::PhyloSuperTreeUnlinked(SuperAlignment *alignment)
: PhyloSuperTree(alignment)
{
}


void PhyloSuperTreeUnlinked::readTree(istream &in, bool &is_rooted) {
    for (iterator it = begin(); it != end(); it++) {
        (*it)->readTree(in, is_rooted);
    }
}

void PhyloSuperTreeUnlinked::setAlignment(Alignment *alignment) {
    ASSERT(alignment->isSuperAlignment());
    SuperAlignment *saln = (SuperAlignment*)alignment;
    ASSERT(saln->partitions.size() == size());
    for (int i = 0; i < size(); i++)
        at(i)->setAlignment(saln->partitions[i]);
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
    iterator it;
    for (it = begin(); it != end(); it++) {
        totalBranchNum += (*it)->branchNum * (*it)->getMixlen();
    }
    lenvec.resize(startid + totalBranchNum);
    
    for (iterator it = begin(); it != end(); it++) {
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
