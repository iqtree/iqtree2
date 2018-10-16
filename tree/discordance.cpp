//
//  discordance.cpp
//  tree
//
//  Created by Minh Bui on 24/9/18.
//

#include "phylotree.h"

void PhyloTree::computeSiteConcordance() {
    BranchVector branches;
    getInnerBranches(branches);
    for (auto it = branches.begin(); it != branches.end(); it++) {
        double num_sites;
        double sup = computeSiteConcordance((*it), num_sites);
        string sup_str = convertDoubleToString(round(sup*1000)/10);
        Node *node = it->second;

        Neighbor *nei = it->second->findNeighbor(it->first);
        
        nei->putAttr("sCF", round(sup*1000)/10);
        nei->putAttr("sN", num_sites);
        
        if (Params::getInstance().newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']') {
                node->name += "[&sCF=" + sup_str + "]";
            } else
                node->name = node->name.substr(0, node->name.length()-1) + ",!sCF=" + sup_str + "]";
        } else {
            if (!node->name.empty())
                node->name += "/";
            node->name += sup_str;
        }
    }
}

double PhyloTree::computeSiteConcordance(Branch &branch, double &num_sites) {
    vector<IntVector> taxa;
    taxa.resize(4);

    // extract the taxa from the two left subtrees
    int id = 0;
    FOR_NEIGHBOR_DECLARE(branch.first, branch.second, it) {
        getTaxaID(taxa[id], (*it)->node, branch.first);
        id++;
        if (id > 2)
            outError(__func__, " only work with bifurcating tree");
    }

    // extract the taxa from the two right subtrees
    FOR_NEIGHBOR(branch.second, branch.first, it) {
        getTaxaID(taxa[id], (*it)->node, branch.second);
        id++;
        if (id > 4)
            outError(__func__, " only work with bifurcating tree");
    }
    
    double sum_support = 0.0;
    int sum_sites = 0;
    for (int i = 0; i < Params::getInstance().site_concordance; i++) {
        int j;
        // get a random quartet
        IntVector quartet;
        quartet.resize(taxa.size());
        for (j = 0; j < taxa.size(); j++) {
            quartet[j] = taxa[j][random_int(taxa[j].size())];
        }

        int support[3] = {0, 0, 0};
        for (auto pat = aln->begin(); pat != aln->end(); pat++) {
            if (!pat->isInformative()) continue;
            bool informative = true;
            for (j = 0; j < quartet.size(); j++)
                if (pat->at(quartet[j]) >= aln->num_states) {
                    informative = false;
                    break;
                }
            if (!informative) continue;
            if (pat->at(quartet[0]) == pat->at(quartet[1]) && pat->at(quartet[2]) == pat->at(quartet[3]) && pat->at(quartet[0]) != pat->at(quartet[2]))
                support[0] += pat->frequency;
            if (pat->at(quartet[0]) == pat->at(quartet[2]) && pat->at(quartet[1]) == pat->at(quartet[3]) && pat->at(quartet[0]) != pat->at(quartet[1]))
                support[1] += pat->frequency;
            if (pat->at(quartet[0]) == pat->at(quartet[3]) && pat->at(quartet[1]) == pat->at(quartet[2]) && pat->at(quartet[0]) != pat->at(quartet[1]))
                support[2] += pat->frequency;
        }
        int sum = support[0] + support[1] + support[2];
        sum_sites += sum;
        if (sum > 0)
            sum_support += ((double)support[0]) / sum;
    }
    num_sites = (double)sum_sites / Params::getInstance().site_concordance;
    return sum_support / Params::getInstance().site_concordance;
}

/**
 assign branch supports to a target tree
 */
void PhyloTree::computeGeneConcordance(MTreeSet &trees) {
    SplitGraph mysg;
    NodeVector mynodes, nodes1, nodes2;
    convertSplits(mysg, &mynodes, root->neighbors[0]->node);
    getBranches(nodes1, nodes2, root->neighbors[0]->node, NULL, true);
    vector<Split*> subtrees;
    extractQuadSubtrees(subtrees, root->neighbors[0]->node);
    IntVector decisive_counts;
    decisive_counts.resize(mynodes.size(), 0);
    StrVector occurence_trees; // list of tree IDs where each split occurs
    if (verbose_mode >= VB_MED)
        occurence_trees.resize(mynodes.size());
    SplitGraph::iterator sit;
    for (sit = mysg.begin(); sit != mysg.end(); sit++)
        (*sit)->setWeight(0.0);
    int treeid, taxid;
    for (treeid = 0; treeid < trees.size(); treeid++) {
        MTree *tree = trees[treeid];
        StrVector taxname;
        tree->getTaxaName(taxname);
        // create the map from taxa between 2 trees
        Split taxa_mask(leafNum);
        for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
            taxid = mysg.findLeafName(*it);
            if (taxid < 0)
                outError("Taxon not found in full tree: ", *it);
            taxa_mask.addTaxon(taxid);
        }
        // make the taxa ordering right before converting to split system
        taxname.clear();
        int smallid;
        for (taxid = 0, smallid = 0; taxid < leafNum; taxid++)
            if (taxa_mask.containTaxon(taxid)) {
                taxname.push_back(mysg.getTaxa()->GetTaxonLabel(taxid));
                string name = (string)mysg.getTaxa()->GetTaxonLabel(taxid);
                tree->findLeafName(name)->id = smallid++;
            }
        ASSERT(taxname.size() == tree->leafNum);
        
        SplitGraph sg;
        //NodeVector nodes;
        tree->convertSplits(sg);
        SplitIntMap hash_ss;
        for (sit = sg.begin(); sit != sg.end(); sit++)
            hash_ss.insertSplit((*sit), 1);
        
        // now scan through all splits in current tree
        int id, qid;
        for (sit = mysg.begin(), id = 0, qid = 0; sit != mysg.end(); sit++, id++)
            if ((*sit)->trivial() < 0) // it is an internal split
            {
                
                bool decisive = true;
                for (int i = 0; i < 4; i++) {
                    if (!taxa_mask.overlap(*subtrees[qid+i])) {
                        decisive = false;
                        break;
                    }
                }
                qid += 4;
                if (!decisive) continue;
                
                decisive_counts[id]++;
                Split *subsp = (*sit)->extractSubSplit(taxa_mask);
                if (subsp->shouldInvert())
                    subsp->invert();
                Split *sp = hash_ss.findSplit(subsp);
                if (sp && sp->trivial() < 0) {
                    (*sit)->setWeight((*sit)->getWeight()+1.0);
                    if (verbose_mode >= VB_MED)
                        occurence_trees[id] += convertIntToString(treeid+1) + " ";
                    if (verbose_mode >= VB_MAX) {
                        for (taxid = 0; taxid < (*sit)->getNTaxa(); taxid++)
                            if ((*sit)->containTaxon(taxid))
                                cout << " " << mysg.getTaxa()->GetTaxonLabel(taxid);
                        cout << " --> ";
                        for (taxid = 0; taxid < sp->getNTaxa(); taxid++)
                            if (sp->containTaxon(taxid))
                                cout << " " << taxname[taxid];
                        cout << endl;
                    }
                }
                delete subsp;
            }
        
    }
    
    ASSERT(nodes1.size() == mynodes.size());
    for (int i = 0; i < mysg.size(); i++)
        if (!mynodes[i]->isLeaf())
        {
            ASSERT(nodes1[i] == mynodes[i] || nodes2[i] == mynodes[i]);
            Neighbor *nei;
            if (nodes1[i] == mynodes[i])
                nei = nodes1[i]->findNeighbor(nodes2[i]);
            else
                nei = nodes2[i]->findNeighbor(nodes1[i]);
            nei->putAttr("gCF", round((mysg[i]->getWeight()/decisive_counts[i])*1000)/10);
            nei->putAttr("gN", decisive_counts[i]);
            
            stringstream tmp;
            if (mysg[i]->getWeight() == 0.0)
                tmp << "0";
            else
                tmp << round((mysg[i]->getWeight()/decisive_counts[i])*1000)/10;
            if (verbose_mode >= VB_MED)
                tmp << "%" << decisive_counts[i];
            
            if (Params::getInstance().newick_extended_format) {
                if (mynodes[i]->name.empty() || mynodes[i]->name.back() != ']')
                    mynodes[i]->name += "[&CF=" + tmp.str() + "]";
                else
                    mynodes[i]->name = mynodes[i]->name.substr(0, mynodes[i]->name.length()-1) + ",!CF=" + tmp.str() + "]";
            } else {
                if (!mynodes[i]->name.empty())
                    mynodes[i]->name.append("/");
                mynodes[i]->name.append(tmp.str());
            }
            if (verbose_mode >= VB_MED) {
                cout << mynodes[i]->name << " " << occurence_trees[i] << endl;
            }
        }
    for (vector<Split*>::reverse_iterator it = subtrees.rbegin(); it != subtrees.rend(); it++)
        delete (*it);
}

/**
 compute quartet internode certainty, similar to Zhou et al (biorxiv)
 */
void PhyloTree::computeQuartetConcordance(MTreeSet &trees) {
    
    outError("Not working yet, need consent from Zhou et al.");
    BranchVector branches;
    getInnerBranches(branches);
    
    for (auto treeit = trees.begin(); treeit != trees.end(); treeit++) {
        
    }
    
    for (auto it = branches.begin(); it != branches.end(); it++) {
        Node *node = it->second;
        double sup = computeQuartetConcordance(*it, trees);
        string sup_str = convertDoubleToString(sup);
        
        if (Params::getInstance().newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']') {
                node->name += "[&qCF=" + sup_str + "]";
            } else
                node->name = node->name.substr(0, node->name.length()-1) + ",!sCF=" + sup_str + "]";
        } else {
            if (!node->name.empty())
                node->name += "/";
            node->name += sup_str;
        }
    }

}

double PhyloTree::computeQuartetConcordance(Branch &branch, MTreeSet &trees) {
    vector<IntVector> taxa;
    taxa.resize(4);
    if (branch.first->degree() != 3 || branch.second->degree() != 3)
        outError(__func__, " only work with bifurcating tree");

    // extract the taxa from the two left subtrees
    int id = 0;
    FOR_NEIGHBOR_DECLARE(branch.first, branch.second, it) {
        getTaxaID(taxa[id], (*it)->node, branch.first);
        id++;
    }
    
    // extract the taxa from the two right subtrees
    FOR_NEIGHBOR(branch.second, branch.first, it) {
        getTaxaID(taxa[id], (*it)->node, branch.second);
        id++;
    }
    
    double sum_support = 0.0;
    int num_quartets = Params::getInstance().site_concordance; // TODO: change name
    for (int i = 0; i < num_quartets; i++) {
        int j;
        // get a random quartet
        IntVector quartet;
        quartet.resize(taxa.size());
        for (j = 0; j < taxa.size(); j++) {
            quartet[j] = taxa[j][random_int(taxa[j].size())];
        }
        
        int quartetCF[3] = {0, 0, 0};
        for (auto tree = trees.begin(); tree != trees.end(); tree++) {
        }
        int sum = quartetCF[0] + quartetCF[1] + quartetCF[2];
        if (sum > 0)
            sum_support += ((double)quartetCF[0]) / sum;
    }
    return sum_support / num_quartets;
}

