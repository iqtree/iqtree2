/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "splitgraph.h"
#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include "tree/node.h"
#include "ncl/ncl.h"
#include "nclextra/myreader.h"
#include "tree/mtree.h"
#include "tree/mtreeset.h"


bool compareSplit(Split* sp1, Split* sp2) {
    if (sp1->countTaxa() != sp2->countTaxa())
        return sp1->countTaxa() < sp2->countTaxa();
    else
        return sp1->firstTaxon() < sp2->firstTaxon();
}

//#define MY_DEBUG
/********************************************************
    Defining SplitGraph methods
********************************************************/
SplitGraph::SplitGraph()
        : vector<Split*>()
{
    pda = NULL;
    taxa = NULL;
    splits = NULL;
    sets = NULL;
    trees = NULL;
    mtrees = NULL;
    areas_boundary = NULL;
}

SplitGraph::SplitGraph(Params &params) : vector<Split*>() {
    init(params);
}

void SplitGraph::createBlocks() {
    taxa = new NxsTaxaBlock();
    splits = new MSplitsBlock(this);
    pda = new MPdaBlock(this);
    sets = new MSetsBlock();
    trees = new TreesBlock(taxa);
    //mtrees = NULL;
}


void SplitGraph::init(Params &params)
{
    mtrees = NULL;
    if (params.intype == IN_NEWICK) {
        // read the input file, can contain more than 1 tree
        mtrees = new MTreeSet(params.user_file, params.is_rooted, params.tree_burnin, params.tree_max_count);
        //mtree = new MTree(params.user_file, params.is_rooted);

        if (params.is_rooted) {
            params.sub_size++;
            params.min_size++;
        }
        if (mtrees->isRooted() && params.root != NULL)
            outError(ERR_CONFLICT_ROOT);
        //SplitIntMap hash_ss;
        mtrees->convertSplits(*this, params.split_threshold, params.split_weight_summary, params.split_weight_threshold);

        if (verbose_mode >= VB_DEBUG)
            saveFileStarDot(cout);
    } else {
        createBlocks();
//        if (params.is_rooted)
//            outError(ERR_ROOT_NET);
    
         cout << "Reading input file " << params.user_file << "..." << endl;

        MyReader nexus(params.user_file);
    
        nexus.Add(taxa);
        nexus.Add(splits);
        nexus.Add(pda);
        nexus.Add(sets);
        nexus.Add(trees);

        MyToken token(nexus.inf);
        nexus.Execute(token);
        if (trees->GetNumTrees() > 0) { 
            if (getNSplits() > 0) 
                outError("Ambiguous input file, pls only specify either SPLITS block or TREES block");
            convertFromTreesBlock(params.tree_burnin, params.tree_max_count, params.split_threshold, 
                params.split_weight_summary, params.split_weight_threshold, params.tree_weight_file);
        }

    }
    
    if (verbose_mode >= VB_DEBUG)
        taxa->Report(cout);
    //splits->Report(cout);
    //reportConflict(cout);
    if (params.pdtaxa_file != NULL) {
        if (sets->getNSets() > 0)
            outError("Taxa sets were already specified in the input file");
        cout << "Reading taxa sets in file " << params.pdtaxa_file << "..." << endl;
    
        bool nexus_formated = (detectInputFile(params.pdtaxa_file) == IN_NEXUS);
        if (nexus_formated) {
            MyReader nexus(params.pdtaxa_file);
            nexus.Add(sets);
            MyToken token(nexus.inf);
            nexus.Execute(token);
        } else {
            readTaxaSets(params.pdtaxa_file, sets);
        }
        if (sets->getNSets() == 0)
            outError("No taxa sets found");
    }

    areas_boundary = NULL;
    if (params.areas_boundary_file) {
        if (sets->getNSets() == 0) outError("No taxon sets defined yet");
        areas_boundary = new double [sets->getNSets() * sets->getNSets()];
        cout << "Reading sets relation file " << params.areas_boundary_file << "..." << endl;
        readAreasBoundary(params.areas_boundary_file, sets, areas_boundary);
    }

    if (verbose_mode >= VB_DEBUG && sets->getNSets() > 0)
        sets->Report(cout);

    if (sets->getNSets() > 0 && taxa->GetNumTaxonLabels() == 0) {
        AddTaxaFromSets();
    }
    if (taxa->GetNumTaxonLabels() == 0)
        outError("No taxa found");
    if (getNSplits() == 0) {
        //outError(ERR_NO_SPLITS);
        createStarTree();
    }
    cout << getNTaxa()-params.is_rooted <<
        " taxa and " << getNSplits()-params.is_rooted << " splits." << endl;

}


void SplitGraph::saveCheckpoint() {
    if (empty()) return;
    int ntax = getNTaxa();
//    checkpoint->startStruct("S");
    CKP_SAVE(ntax);
    int nsplits = size();
    CKP_SAVE(nsplits);
    checkpoint->startList(size());
    for (iterator it = begin(); it != end(); it++) {
        checkpoint->addListElement();
        stringstream ss;
        ss << (*it)->getWeight();
        for (int i = 0; i < ntax; i++)
            if ((*it)->containTaxon(i))
                ss << " " << i;
        checkpoint->put("", ss.str());
    }
    checkpoint->endList();
//    checkpoint->endStruct();
    CheckpointFactory::saveCheckpoint();
}

void SplitGraph::restoreCheckpoint() {
    int ntax, nsplits;
    CheckpointFactory::restoreCheckpoint();
//    checkpoint->startStruct("S");

    if (!CKP_RESTORE(ntax)) return;
    CKP_RESTORE(nsplits);
    checkpoint->startList(nsplits);
    for (int split = 0; split < nsplits; split++) {
        checkpoint->addListElement();
        string str;
        bool found = checkpoint->getString("", str);
        ASSERT(found);
        stringstream ss(str);
        double weight;
        ss >> weight;
        Split *sp = new Split(ntax, weight);
        for (int i = 0; i < ntax; i++) {
            int tax;
            if (ss >> tax) {
                sp->addTaxon(tax);
            } else
                break;
        }
        push_back(sp);
    }
    checkpoint->endList();
//    checkpoint->endStruct();
}

int SplitGraph::getNTrivialSplits() {
    int count = 0;
    for (iterator it = begin(); it != end(); it++)
        if ((*it)->trivial() >= 0)
            count++;
    return count;
}


void SplitGraph::createStarTree() {
    cout << "No splits found, creating a star tree with branch length of 1..." << endl;
    int ntaxa = taxa->GetNumTaxonLabels();
    for (int i = 0; i < ntaxa; i++) {
        Split *sp = new Split(ntaxa, 1.0);
        sp->addTaxon(i);
        push_back(sp);
    }
    cout << "NOTE: subsequent PD will correspond to species richness." << endl;
}


void SplitGraph::AddTaxaFromSets() {
    cout << "Taking taxa from SETS block..." << endl;
    for (int i = 0; i < sets->getNSets(); i++)
        for(vector<string>::iterator it = sets->getSet(i)->taxlist.begin(); 
            it != sets->getSet(i)->taxlist.end(); it++) 
            if (!taxa->IsAlreadyDefined(NxsString(it->c_str()))) {
                taxa->AddTaxonLabel(NxsString(it->c_str()));
            }    
}

void SplitGraph::freeMem() {
    for (reverse_iterator it = rbegin(); it != rend(); it++) {
        //(*it)->report(cout);
        delete *it;
    }
    clear();
    if (areas_boundary) delete areas_boundary;
    if (trees) delete trees;
    if (sets) delete sets;
    if (pda) delete pda;
    if (splits) delete splits;
    if (taxa) delete taxa;
    if (mtrees) delete mtrees;
}

SplitGraph::~SplitGraph()
{
    freeMem();
}


void SplitGraph::convertFromTreesBlock(int burnin, int max_count, double split_threshold, 
    int split_weight_summary, double weight_threshold, const char *tree_weight_file) {
    cout << trees->GetNumTrees() << " tree(s) loaded" << endl;
    if (burnin >= trees->GetNumTrees())
        outError("Burnin value is too large");
    if (burnin > 0)
    cout << burnin << " beginning tree(s) discarded" << endl;
    mtrees = new MTreeSet();
    
    for (int i = burnin; i < trees->GetNumTrees() && (i < burnin+max_count); i++) {
        stringstream strs(trees->GetTranslatedTreeDescription(i), ios::in | ios::out | ios::app);
        strs << ";";
        MTree *tree = mtrees->newTree();
        bool myrooted = trees->IsRootedTree(i);
        tree->readTree(strs, myrooted);
        mtrees->push_back(tree);
        mtrees->tree_weights.push_back(1);
    }
    mtrees->checkConsistency();
    //SplitIntMap hash_ss;
    
    if (tree_weight_file) 
        readIntVector(tree_weight_file, burnin, max_count, mtrees->tree_weights);
/*    else if (!weights)
        tree_weights.resize(size(), 1);*/

    if (mtrees->size() != mtrees->tree_weights.size())
        outError("Tree file and tree weight file have different number of entries");    
    mtrees->convertSplits(*this, split_threshold, split_weight_summary, weight_threshold);
}



void SplitGraph::report(ostream &out)
{

    out << endl;
    out << "Split network contains ";

    if (size() == 0)
    {
        out << "no split" << endl;
    }
    else if (size() == 1)
        out << "one split" << endl;
    else
        out << size() << " splits" << endl;

    if (size() == 0)
        return;

    sort(begin(), end(), compareSplit);
    int k = 0;
    for (iterator it = begin(); it != end(); it++,k++)
    {
        out << '\t' << (k+1) << '\t';
        (*it)->report(out);
    }
}

void SplitGraph::reportConflict(ostream &out)
{
    int k = 0;
    out << "Compatible splits: " << endl;
    for (iterator i = begin(); i != end(); i++, k++)
    {
        out << (k+1) << '\t';
        int k2 = 1;
        for (iterator j = begin(); j != end(); j++, k2++)
            if ( j != i && (*i)->compatible(*(*j)))
            {
                out << k2 << " ";
            }
        out << endl;
    }
}

/**
    calculate sum of weights of preserved splits in the taxa_set
    @param taxa_set a set of taxa
*/
double SplitGraph::calcWeight(Split &taxa_set)
{
    double sum = 0.0;
    for (iterator it = begin(); it != end(); it++)
        if ((*it)->preserved(taxa_set))
            sum += (*it)->getWeight();
    return sum;
}

int SplitGraph::countSplits(Split &taxa_set)
{
    int cnt = 0;
    for (iterator it = begin(); it != end(); it++)
        if ((*it)->preserved(taxa_set))
            cnt++;
    return cnt;
}

int SplitGraph::countInternalSplits(Split &taxa_set)
{
    int cnt = 0;
    for (iterator it = begin(); it != end(); it++)
        if ((*it)->trivial() < 0 && (*it)->preserved(taxa_set))
            cnt++;
    return cnt;
}

/**
    calculate sum of weights of all splits
*/
double SplitGraph::calcWeight() {
    double sum = 0.0;
    for (iterator it = begin(); it != end(); it++)
        sum += (*it)->weight;
    return sum;
}

double SplitGraph::calcTrivialWeight() {
    double sum = 0.0;
    for (iterator it = begin(); it != end(); it++)
        if ((*it)->trivial() >= 0)
            sum += (*it)->weight;
    return sum;
}

double SplitGraph::maxWeight() {
    double m = -1e6;
    for (iterator it = begin(); it != end(); it++)
        if (m < (*it)->weight) m = (*it)->weight;
    return m;
}

void SplitGraph::generateTaxaSet(char *filename, int size, int overlap, int times) {
    ofstream out(filename);
    if (!out.is_open())
        outError(ERR_WRITE_OUTPUT, filename);
    ASSERT(overlap <= size);
    int total = 2*size - overlap;
    int ntaxa = getNTaxa();
    for (int cnt = 0; cnt < times; cnt++) {
        // generate random taxon index 
        IntVector ranvec;
        BoolVector occur(ntaxa, false);
        int i;
        for (i = 0; i < total; i++) {
            int rnum;
            do { rnum = random_int(ntaxa); } while (occur[rnum]);
            ranvec.push_back(rnum);
            occur[rnum] = true;
        }
        // now write the first set
        out << size << endl;
        for (i = 0; i < size; i++) 
            out << taxa->GetTaxonLabel(ranvec[i]) << endl;
        out << endl;
        // now write the second set
        out << size << endl;
        for (i = size-overlap; i < total; i++) 
            out << taxa->GetTaxonLabel(ranvec[i]) << endl;
        out << endl;
    }
    out.close();
}

void SplitGraph::calcDistance(char *filename) {
    ofstream out(filename);
    if (!out.is_open())
        outError(ERR_WRITE_OUTPUT, filename);
    mmatrix(double) dist;
    int i, j;    
    calcDistance(dist);

    int ntaxa = getNTaxa();

    // now write the distances in phylip .dist format
    out << ntaxa << endl;
    
    for (i = 0; i < ntaxa; i++) {
        out << taxa->GetTaxonLabel(i) << "   ";
        for (j = 0; j < ntaxa; j++) {
            out << dist[i][j] << "  ";
        }
        out << endl;
    }
    out.close();
}

void SplitGraph::calcDistance(mmatrix(double) &dist) {
    int ntaxa = getNTaxa();
    iterator it;
    vector<int> vi, vj;
    vector<int>::iterator i, j;

    dist.resize(ntaxa);
    for (mmatrix(double)::iterator di = dist.begin(); di != dist.end(); di++)
        (*di).resize(ntaxa, 0);

    for (it = begin(); it != end(); it++) {
        (*it)->getTaxaList(vi, vj);
        for (i = vi.begin(); i != vi.end(); i++)
            for (j = vj.begin(); j < vj.end(); j++) {
                dist[*i][*j] += (*it)->weight;
                dist[*j][*i] += (*it)->weight;
            }
    }

}


void SplitGraph::calcDistance(mmatrix(double) &dist, vector<int> &taxa_order) {
    int ntaxa = getNTaxa();
    int i, j;

    mmatrix(double) my_dist;
    calcDistance(my_dist);
    dist.resize(ntaxa);
    for (i = 0; i < ntaxa; i++) {
        dist[i].resize(ntaxa);
        for (j = 0; j < ntaxa; j++)
            dist[i][j] = my_dist[taxa_order[i]][taxa_order[j]];
    }
}

bool SplitGraph::checkCircular(mmatrix(double) &dist) {
    return true;
    int ntaxa = getNTaxa();
    Split taxa_set(ntaxa, 0.0);
    for (int i = 0; i < ntaxa-2; i++)
        for (int j = i+1; j < ntaxa-1; j++)
            for (int k = j+1; k < ntaxa; k++) {
                taxa_set.addTaxon(i);
                taxa_set.addTaxon(j);
                taxa_set.addTaxon(k);
                taxa_set.weight = calcWeight(taxa_set);
                if (fabs(2 * taxa_set.weight - (dist[i][j] + dist[i][k] + dist[j][k])) > 0.0001) {
                    cout << "Taxa " << i << " " << j << " " << k;
                    cout << " do not satisfy circular equation!" << endl;
                    cout << "Weight = " << taxa_set.weight << endl;
                    cout << "Sum dist/2 = " << (dist[i][j] + dist[i][k] + dist[j][k]) / 2.0 << endl;
                    cout << "dist = " << dist[i][j] << " " << dist[i][k] << " "
                         << dist[j][k] << endl;
                    return false;
                }
                taxa_set.removeTaxon(i);
                taxa_set.removeTaxon(j);
                taxa_set.removeTaxon(k);
            }
    return true;
}

void SplitGraph::generateCircular(Params &params) {
    int i, j;
    int ntaxa = params.sub_size;
    int num_splits = (params.num_splits > 0) ? params.num_splits : 3 * ntaxa;
    if (num_splits < ntaxa) 
        outError(ERR_FEW_SPLITS); 

    taxa = new NxsTaxaBlock();
    splits = new MSplitsBlock(this);

    double threshold = (ntaxa > 3) ? (double)(num_splits - ntaxa) / (ntaxa*(ntaxa-3)/2) : 0.0;

    // insert all trivial splits
    for (i = 0; i < ntaxa; i++) {
        double weight = randomLen(params);
        Split *sp = new Split(ntaxa, weight);
        sp->addTaxon(i);
        push_back(sp);
        ostringstream str;
        str << "T" << (i+1);
        taxa->AddTaxonLabel(NxsString(str.str().c_str()));
        splits->cycle.push_back(i);
    }

    // randomly insert internal splits
    for (i = 0; i < ntaxa-2 && getNSplits() < num_splits; i++)
        for (j = i+1; j < ntaxa && j < ntaxa-3+i; j++) {
            double choice = random_double();
            if (choice > threshold) continue;
            double weight = randomLen(params);
            Split *sp = new Split(ntaxa, weight);
            for (int k = i; k <= j; k++)
                sp->addTaxon(k);
            push_back(sp);
            if (getNSplits() >= num_splits) break;
        }

    ofstream out(params.user_file);
    if (!out.is_open()) {
        outError(ERR_WRITE_OUTPUT, params.user_file);
    }

    saveFileNexus(out);
    out.close();
} 

void SplitGraph::saveFileNexus(ostream &out, bool omit_trivial) {
    int ntaxa = getNTaxa();
    int i;
    out << "#nexus" << endl << endl;
    out << "BEGIN Taxa;" << endl;
    out << "DIMENSIONS ntax=" << ntaxa << ";" << endl;
    out << "TAXLABELS" << endl;
    for (i = 0; i < ntaxa; i++)
        out << "[" << i+1 << "] '" << taxa->GetTaxonLabel(i) << "'" << endl;
    out << ";" << endl << "END; [Taxa]" << endl << endl;
    out << "BEGIN Splits;" << endl;
    out << "DIMENSIONS ntax=" << ntaxa << " nsplits=" << ((omit_trivial) ? getNSplits() - getNTrivialSplits() : getNSplits()) << ";" << endl;
    out << "FORMAT labels=no weights=yes confidences=no intervals=no;" << endl;
    if (isCircular()) {
        out << "CYCLE";
        for (i = 0; i < ntaxa; i++) 
            out << " " << splits->cycle[i] + 1;
        out << ";" << endl;
    }
    out << "MATRIX" << endl;
    int near_zeros = 0;
    int zeros = 0;
    for (iterator it = begin(); it != end(); it++) {
        if (omit_trivial && (*it)->trivial() >= 0) continue;
        if ((*it)->weight == 0.0) zeros ++;
        if ((*it)->weight <= 1e-6) near_zeros ++;
        out << "\t" << (*it)->weight << "\t";
        for (i = 0; i < ntaxa; i++)
            if ((*it)->containTaxon(i))
                out << " " << i+1;
        out << "," << endl;
    }
    out << ";" << endl << "END; [Splits]" << endl << endl;
    if (near_zeros) {
        //outWarning("Some nearly-zero split weights observed");
        //cout << zeros << " zero-weights and " << near_zeros << " near zero weights!" << endl;
    }
}

void SplitGraph::saveFileStarDot(ostream &out, bool omit_trivial) {
    int ntaxa = getNTaxa();
    int i;
    for (iterator it = begin(); it != end(); it++) {
        if (omit_trivial && (*it)->trivial() >= 0) continue;
        bool swap_code = !(*it)->containTaxon(0);
        if (swap_code) {
            for (i = 0; i < ntaxa; i++)
                out << (((*it)->containTaxon(i)) ? '.' : '*');
        } else {
            for (i = 0; i < ntaxa; i++)
                out << (((*it)->containTaxon(i)) ? '*' : '.');
        }
        out << "\t" << (*it)->weight << endl;
    }
}

void SplitGraph::saveFile(const char* out_file, InputType file_format, bool omit_trivial) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(out_file);
        if (file_format == IN_NEXUS) 
            saveFileNexus(out, omit_trivial);
        else
            saveFileStarDot(out, omit_trivial);
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, out_file);
    }
}

void SplitGraph::scaleWeight(double norm, bool make_int, int precision) {
    for (iterator itg = begin(); itg != end(); itg ++ )
        if (make_int)
            (*itg)->setWeight( round((*itg)->getWeight()*norm) );
        else if (precision < 0)
            (*itg)->setWeight( (*itg)->getWeight()*norm);
        else 
            (*itg)->setWeight( round((*itg)->getWeight()*norm*pow((double)10.0,precision))/pow((double)10.0,precision));
}
// TODO Implement a more efficient function using Hash Table
bool SplitGraph::containSplit(Split &sp) {
    Split invert_sp(sp);
    invert_sp.invert();
    for (iterator it = begin(); it != end(); it++)
        if ((*(*it)) == sp || (*(*it)) == invert_sp)
            return true;
    return false;
}

double SplitGraph::computeBoundary(Split &area) {
    if (!areas_boundary) return 0.0;
    int nareas = sets->getNSets();
    double boundary = 0.0;
    for (int i = 0; i < nareas; i++) 
    if (area.containTaxon(i)) {
        boundary += areas_boundary[i*nareas+i];
        for (int j = i+1; j < nareas; j++) 
            if (area.containTaxon(j))
                boundary -= 2.0 * areas_boundary[i*nareas+j];
    }
    return boundary;
}

bool SplitGraph::compatible(Split *sp) {
    for (iterator it = begin(); it != end(); it++)
        if (!(*it)->compatible(*sp))
            return false;
    return true;
}

void SplitGraph::findMaxCompatibleSplits(SplitGraph &maxsg) {

    // maximum number of compatible splits = 2n-3!
    int max_splits = getNTaxa() * 2 - 3;

    // myset will be sorted by weight in descending order
    SplitSet myset;
    myset.insert(myset.end(), begin(), end());
    sort(myset.begin(), myset.end(), splitweightcmp);

    // now build the spset
    if (!maxsg.taxa)
        maxsg.taxa = new NxsTaxaBlock();
    if (!maxsg.splits)
        maxsg.splits = new MSplitsBlock(&maxsg);
    if (!maxsg.pda)
        maxsg.pda = new MPdaBlock(&maxsg);

    for (int i = 0; i < getNTaxa(); i++)
        maxsg.taxa->AddTaxonLabel(taxa->GetTaxonLabel(i));
    
    // make the cycle
    maxsg.splits->cycle = splits->cycle;
    // make the splits

    for (SplitSet::iterator it = myset.begin(); it != myset.end(); it++) 
        if (maxsg.compatible(*it)){
            maxsg.push_back(new Split(*(*it)));
            //(*it)->report(cout);
            if (maxsg.size() >= max_splits)
                break;
        }
    myset.clear();
}

bool SplitGraph::isWeaklyCompatible() {
    if (getNSplits() < 3) return true;
    for (iterator it1 = begin(); it1+2 != end(); it1++)
        for (iterator it2 = it1+1; it2+1 != end(); it2++)
            for (iterator it3 = it2+1; it3 != end(); it3++) {
                Split sp1(*(*it1));
                Split sp2(*(*it2));
                Split sp3(*(*it3));
                Split sp(sp1);
                sp *= sp2;
                sp *= sp3;
                if (sp.isEmpty()) continue;
                sp1.invert();
                sp2.invert();
                sp = sp1;
                sp *= sp2;
                sp *= sp3;
                if (sp.isEmpty()) continue;
                sp2.invert();
                sp3.invert();
                sp = sp1;
                sp *= sp2;
                sp *= sp3;
                if (sp.isEmpty()) continue;
                sp1.invert();
                sp2.invert();
                sp = sp1;
                sp *= sp2;
                sp *= sp3;
                if (sp.isEmpty()) continue;
                return false;
            }
    return true;
}


void SplitGraph::getTaxaName(vector<string> &taxname) {
    taxname.clear();
    for (int i = 0; i < getNTaxa(); i++)
        taxname.push_back(taxa->GetTaxonLabel(i));
}

int SplitGraph::findLeafName(string &name) {
    for (int i = 0; i < getNTaxa(); i++)
        if (taxa->GetTaxonLabel(i) == name)
            return i;
    return -1;
}

int SplitGraph::removeTrivialSplits() {
    int removed = 0;
    for (iterator itg = begin(); itg != end(); )  {
        if ((*itg)->trivial() >= 0) {
            removed++;
            delete (*itg);
            (*itg) = back();
            pop_back(); 
        } else itg++;
    }
    return removed;
}
