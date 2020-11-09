//
//  discordance.cpp
//  tree
//
//  Created by Minh Bui on 24/9/18.
//

#include "phylosupertree.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define PUT_MEANING(value, description) meanings.insert({#value, description})

void PhyloTree::computeSiteConcordance(map<string,string> &meanings) {
    BranchVector branches;
    getInnerBranches(branches);
#ifdef _OPENMP
#pragma omp parallel
    {
        int *rstream;
        init_random(params->ran_seed + omp_get_thread_num(), false, &rstream);
#pragma omp for
#else
        int *rstream = randstream;
#endif
    for (auto ii = 0; ii < branches.size(); ii++) {
        BranchVector::iterator it = branches.begin()+ii;
        computeSiteConcordance((*it), params->site_concordance, rstream);
        Neighbor *nei = it->second->findNeighbor(it->first);
        double sCF = 0.0;
        if (!GET_ATTR(nei, sCF))
            continue;

        stringstream tmp;
        tmp.precision(3);
        tmp << sCF;
        string sup_str = tmp.str();
        Node *node = it->second;
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
#ifdef _OPENMP
        finish_random(rstream);
    }
#endif

    PUT_MEANING(sCF, "Site concordance factor averaged over " + convertIntToString(params->site_concordance) +  " quartets (=sCF_N/sN %)");
    PUT_MEANING(sN, "Number of informative sites averaged over " + convertIntToString(params->site_concordance) +  " quartets");
    PUT_MEANING(sDF1, "Site discordance factor for alternative quartet 1 (=sDF1_N/sN %)");
    PUT_MEANING(sDF2, "Site discordance factor for alternative quartet 2 (=sDF2_N/sN %)");
    PUT_MEANING(sCF_N, "sCF in absolute number of sites");
    PUT_MEANING(sDF1_N, "sDF1 in absolute number of sites");
    PUT_MEANING(sDF2_N, "sDF2 in absolute number of sites");
}

void Alignment::computeQuartetSupports(IntVector &quartet, vector<int64_t> &support) {
    // sanity check e.g. when having rooted tree
    for (auto q = quartet.begin(); q != quartet.end(); q++)
        ASSERT(*q < getNSeq());
        
    for (auto pat = begin(); pat != end(); pat++) {
        if (!pat->isInformative()) continue;
        bool informative = true;
        for (int j = 0; j < quartet.size(); j++)
            if (pat->at(quartet[j]) >= num_states) {
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
}

void SuperAlignment::computeQuartetSupports(IntVector &quartet, vector<int64_t> &support) {
    for (int part = 0; part < partitions.size(); part++) {
        IntVector part_quartet;
        for (auto i = quartet.begin(); i != quartet.end(); i++) {
            if (taxa_index[*i][part] >= 0)
                part_quartet.push_back(taxa_index[*i][part]);
            else
                break;
        }
        if (part_quartet.size() != quartet.size())
            continue;
        if (Params::getInstance().site_concordance_partition) {
            vector<int64_t> part_support;
            part_support.resize(3, 0);
            partitions[part]->computeQuartetSupports(part_quartet, part_support);
            for (int j = 0; j < 3; j++) if (part_support[j] > 0) {
                ASSERT(support[part*3+3+j] >= 0);
                support[part*3+3+j] += part_support[j];
                support[j] += part_support[j];
            }
        } else
            partitions[part]->computeQuartetSupports(part_quartet, support);
    }
}

void PhyloTree::computeSiteConcordance(Branch &branch, int nquartets, int *rstream) {
    vector<IntVector> left_taxa, right_taxa;
    //taxa.resize(4);

    // extract the taxa from the two left subtrees
    left_taxa.resize(branch.first->neighbors.size()-1);
    int id = 0;
    FOR_NEIGHBOR_DECLARE(branch.first, branch.second, itNeighbor) {
        // 2018-12-11: do not consider internal branch at the root
        if (rooted && (*itNeighbor)->node == root)
            return;
        getTaxaID(left_taxa[id], (*itNeighbor)->node, branch.first);
        id++;
//        if (id > 2)
//            outError(__func__, " only work with bifurcating tree");
    }
    ASSERT(id == left_taxa.size());

    // extract the taxa from the two right subtrees
    right_taxa.resize(branch.second->neighbors.size()-1);
    id = 0;
    FOR_NEIGHBOR(branch.second, branch.first, itNeighbor) {
        // 2018-12-11: do not consider internal branch at the root
        if (rooted && (*itNeighbor)->node == root)
            return;
        getTaxaID(right_taxa[id], (*itNeighbor)->node, branch.second);
        id++;
//        if (id > 4)
//            outError(__func__, " only work with bifurcating tree");
    }
    
    ASSERT(id == right_taxa.size());
    
    // 2018-12-11: remove root taxon from taxa for rooted tree
    if (rooted) {
        for (auto it = left_taxa.begin(); it != left_taxa.end(); it++)
            for (auto it2 = it->begin(); it2 != it->end(); it2++)
                if (*it2 == leafNum-1) {
                    it->erase(it2);
                    break;
                }
        for (auto it = right_taxa.begin(); it != right_taxa.end(); it++)
            for (auto it2 = it->begin(); it2 != it->end(); it2++)
                if (*it2 == leafNum-1) {
                    it->erase(it2);
                    break;
                }
    }
    
    double sCF = 0.0; // concordance factor
    double sDF1 = 0.0;
    double sDF2 = 0.0;
    double sN = 0.0;
    size_t sum_sites = 0;
    double sCF_N = 0, sDF1_N = 0, sDF2_N = 0;
    vector<int64_t> support;
    support.resize(3, 0);
    // reserve size for partition-wise concordant/discordant sites
    if (Params::getInstance().site_concordance_partition && aln->isSuperAlignment()) {
        SuperAlignment *saln = (SuperAlignment*)aln;
        support.resize(saln->partitions.size()*3+3, 0);
        // check for gene trees not decisive for this branch
        StrVector taxname;
        getTaxaName(taxname);
        size_t part = 0;
        for (auto part_aln = saln->partitions.begin(); part_aln != saln->partitions.end(); part_aln++, part++) {
            // get the taxa names of the partition tree
            StringIntMap name_map;
            for (size_t i = 0; i < (*part_aln)->getNSeq(); i++)
                name_map[(*part_aln)->getSeqName(i)] = i;
            
            // check that at least one taxon from each subtree is present in partition tree
            int left_count = 0, right_count = 0;
            for (auto it = left_taxa.begin(); it != left_taxa.end(); ++it) {
                for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
                    if (name_map.find(taxname[*it2]) != name_map.end()) {
                        left_count++;
                        break;
                    }
                }
            }
            for (auto it = right_taxa.begin(); it != right_taxa.end(); ++it) {
                for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
                    if (name_map.find(taxname[*it2]) != name_map.end()) {
                        right_count++;
                        break;
                    }
                }
            }
            if (left_count < 2 || right_count < 2) {
                // not decisive
                support[part*3+3] = support[part*3+4] = support[part*3+5] = -1;
                break;
            }
        }
    }
    Neighbor *nei = branch.second->findNeighbor(branch.first);
    for (size_t i = 0; i < nquartets; ++i) {
        // get a random quartet
        IntVector quartet;
        quartet.resize(4);
        int left_id0 = 0, left_id1 = 1, right_id0 = 0, right_id1 = 1;
        if (left_taxa.size() > 2) {
            left_id0 = random_int(left_taxa.size(), rstream);
            do {
                left_id1 = random_int(left_taxa.size(), rstream);
            } while (left_id0 == left_id1);
        }
        if (right_taxa.size() > 2) {
            right_id0 = random_int(right_taxa.size(), rstream);
            do {
                right_id1 = random_int(right_taxa.size(), rstream);
            } while (right_id0 == right_id1);
        }
        quartet[0] = left_taxa[left_id0][random_int(left_taxa[left_id0].size(), rstream)];
        quartet[1] = left_taxa[left_id1][random_int(left_taxa[left_id1].size(), rstream)];
        quartet[2] = right_taxa[right_id0][random_int(right_taxa[right_id0].size(), rstream)];
        quartet[3] = right_taxa[right_id1][random_int(right_taxa[right_id1].size(), rstream)];

        support[0] = support[1] = support[2] = 0;
        aln->computeQuartetSupports(quartet, support);
        size_t sum = support[0] + support[1] + support[2];
        sum_sites += sum;
        if (sum > 0) {
            sCF += ((double)support[0]) / sum;
            sDF1 += ((double)support[1]) / sum;
            sDF2 += ((double)support[2]) / sum;
            sCF_N += support[0];
            sDF1_N += support[1];
            sDF2_N += support[2];
        }
        if (params->print_cf_quartets) {
            // print sCF for each quartet
            stringstream ss;
            ss << quartet[0]+1 << '\t' << quartet[1]+1 << '\t' << quartet[2]+1 << '\t' << quartet[3]+1
               << '\t' << (((double)support[0]) / sum) << '\t' << support[0]
               << '\t' << (((double)support[1]) / sum) << '\t' << support[1]
               << '\t' << (((double)support[2]) / sum) <<  '\t' << support[2]
               << '\t' << sum;
            nei->putAttr("q" + convertIntToString(i), ss.str());
        }
    }
    sN = (double)sum_sites / nquartets;
    // rounding
    sCF = round(sCF / nquartets * 10000)/100;
    sDF1 = round(sDF1 / nquartets * 10000)/100;
    sDF2 = round(sDF2 / nquartets * 10000)/100;
    sCF_N = round(sCF_N / nquartets * 100)/100;
    sDF1_N = round(sDF1_N / nquartets * 100)/100;
    sDF2_N = round(sDF2_N / nquartets * 100)/100;
    PUT_ATTR(nei, sCF);
    PUT_ATTR(nei, sN);
    PUT_ATTR(nei, sDF1);
    PUT_ATTR(nei, sDF2);
    PUT_ATTR(nei, sCF_N);
    PUT_ATTR(nei, sDF1_N);
    PUT_ATTR(nei, sDF2_N);
    stringstream s_factors;
    s_factors << sCF << "/" << sDF1 << "/" << sDF2;
    nei->putAttr("sCF/sDF1/sDF2", s_factors.str());
    stringstream s_factors_N;
    s_factors_N << sCF_N << "/" << sDF1_N << "/" << sDF2_N;
    nei->putAttr("sCF_N/sDF1_N/sDF2_N", s_factors_N.str());
    // insert key-value for partition-wise con/discordant sites
    string keys[] = {"sC", "sD1", "sD2"};
    for (size_t i = 3; i < support.size(); ++i) {
        if (support[i] >= 0)
            nei->putAttr(keys[i%3] + convertIntToString(i/3), (double)support[i]/nquartets);
        else
            nei->putAttr(keys[i%3] + convertIntToString(i/3), "NA");
    }
}

/**
 assign branch supports to a target tree
 */
void PhyloTree::computeGeneConcordance(MTreeSet &trees, map<string,string> &meanings) {
    StrVector names;
    getTaxaName(names);
    StringIntMap name_map;
    for (auto stri = names.begin(); stri != names.end(); stri++)
        name_map[*stri] = stri - names.begin();
    BranchVector branches;
    vector<Split*> subtrees;
    extractQuadSubtrees(subtrees, branches, root->neighbors[0]->node);
    IntVector decisive_counts; // number of decisive trees
    decisive_counts.resize(branches.size(), 0);
    IntVector supports[3]; // number of trees supporting 3 alternative splits
    supports[0].resize(branches.size(), 0);
    supports[1].resize(branches.size(), 0);
    supports[2].resize(branches.size(), 0);
    string prefix[3] = {"gC", "gD1", "gD2"};
    int treeid, taxid;
    for (treeid = 0; treeid < trees.size(); treeid++) {
        MTree *tree = trees[treeid];
        StrVector taxname;
        tree->getTaxaName(taxname);
        // create the map from taxa between 2 trees
        Split taxa_mask(leafNum);
        for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
            if (name_map.find(*it) == name_map.end())
                outError("Taxon not found in full tree: ", *it);
            taxa_mask.addTaxon(name_map[*it]);
        }
        // make the taxa ordering right before converting to split system
        taxname.clear();
        int smallid;
        for (taxid = 0, smallid = 0; taxid < leafNum; taxid++)
            if (taxa_mask.containTaxon(taxid)) {
                taxname.push_back(names[taxid]);
                tree->findLeafName(names[taxid])->id = smallid++;
            }
        ASSERT(taxname.size() == tree->leafNum);
        
        SplitGraph sg;
        //NodeVector nodes;
        tree->convertSplits(sg);
        SplitIntMap hash_ss;
        for (auto sit = sg.begin(); sit != sg.end(); sit++)
            hash_ss.insertSplit((*sit), 1);
        
        // now scan through all splits in current tree
        int id, qid;
        for (id = 0, qid = 0; qid < subtrees.size(); id++, qid += 4)
        {
            Neighbor *nei = branches[id].second->findNeighbor(branches[id].first);

            bool decisive = true;
            int i;
            for (i = 0; i < 4; i++) {
                if (!taxa_mask.overlap(*subtrees[qid+i])) {
                    decisive = false;
                    break;
                }
            }
            if (!decisive && params->site_concordance_partition) {
                for (i = 0; i < 3; i++)
                    nei->putAttr(prefix[i] + convertIntToString(treeid+1), "NA");
            }

            if (!decisive) continue;
            
            decisive_counts[id]++;
            for (i = 0; i < 3; i++) {
                Split this_split = *subtrees[qid]; // current split
                this_split += *subtrees[qid+i+1];
                Split *subsp = this_split.extractSubSplit(taxa_mask);
                if (subsp->shouldInvert())
                    subsp->invert();
                int concordant = 0;
                if (hash_ss.findSplit(subsp)) {
                    supports[i][id]++;
                    concordant = 1;
                }
                if (params->site_concordance_partition) {
                    nei->putAttr(prefix[i] + convertIntToString(treeid+1), concordant);
                }
                delete subsp;
            }
        }
        
    }
    
    for (int i = 0; i < branches.size(); i++) {
        if (decisive_counts[i] == 0)
            continue;
        Neighbor *nei = branches[i].second->findNeighbor(branches[i].first);
        int gN = decisive_counts[i];
        int gCF_N = supports[0][i];
        int gDF1_N = supports[1][i];
        int gDF2_N = supports[2][i];
        int gDFP_N = gN - gCF_N - gDF1_N - gDF2_N;
        double gCF = round((double)gCF_N/gN * 10000)/100;
        double gDF1 = round((double)gDF1_N/gN * 10000)/100;
        double gDF2 = round((double)gDF2_N/gN * 10000)/100;
        double gDFP = round((double)gDFP_N/gN * 10000)/100;
        PUT_ATTR(nei, gCF);
        PUT_ATTR(nei, gDF1);
        PUT_ATTR(nei, gDF2);
        PUT_ATTR(nei, gDFP);
        PUT_ATTR(nei, gN);
        PUT_ATTR(nei, gCF_N);
        PUT_ATTR(nei, gDF1_N);
        PUT_ATTR(nei, gDF2_N);
        PUT_ATTR(nei, gDFP_N);
        stringstream g_factors;
        g_factors << gCF << "/" << gDF1 << "/" << gDF2 << "/" << gDFP;
        nei->putAttr("gCF/gDF1/gDF2/gDFP", g_factors.str());
        stringstream g_factors_N;
        g_factors_N << gCF_N << "/" << gDF1_N << "/" << gDF2_N << "/" << gDFP_N;
        nei->putAttr("gCF_N/gDF1_N/gDF2_N/gDFP_N", g_factors_N.str());

        stringstream tmp;
        tmp.precision(3);
        tmp << (double)supports[0][i]/decisive_counts[i]*100;
        if (verbose_mode >= VB_MED)
            tmp << "%" << decisive_counts[i];
        
        Node *node = branches[i].second;
        if (Params::getInstance().newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']')
                node->name += "[&CF=" + tmp.str() + "]";
            else
                node->name = node->name.substr(0, node->name.length()-1) + ",!CF=" + tmp.str() + "]";
        } else {
            if (!node->name.empty())
                node->name.append("/");
            node->name.append(tmp.str());
        }
    }
    for (vector<Split*>::reverse_iterator it = subtrees.rbegin(); it != subtrees.rend(); it++)
        delete (*it);

    PUT_MEANING(gCF, "Gene concordance factor (=gCF_N/gN %)");
    PUT_MEANING(gDF1, "Gene discordance factor for NNI-1 branch (=gDF1_N/gN %)");
    PUT_MEANING(gDF2, "Gene discordance factor for NNI-2 branch (=gDF2_N/gN %)");
    PUT_MEANING(gDFP, "Gene discordance factor due to polyphyly (=gDFP_N/gN %)");
    PUT_MEANING(gN, "Number of trees decisive for the branch");
    PUT_MEANING(gCF_N, "Number of trees concordant with the branch");
    PUT_MEANING(gDF1_N, "Number of trees concordant with NNI-1 branch");
    PUT_MEANING(gDF2_N, "Number of trees concordant with NNI-2 branch");
    PUT_MEANING(gDFP_N, "Number of trees decisive but discordant due to polyphyly");
    meanings.insert({"*NOTE*", "(gCF+gDF1+gDF2+gDFP) = 100% and (gCF_N+gDF1_N+gDF2_N+gDFP_N) = gN"});
    if (Params::getInstance().print_df1_trees) {
        meanings.insert({"treeDF1", "Newick tree for gDF1"});
        meanings.insert({"treeDF2", "Newick tree for gDF2"});
    }
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

