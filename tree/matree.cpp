/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
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
#include "matree.h"

void MaTree::printBrInfo(ostream& out) {
    //to store internal branch lengths
    DoubleVector inner;
    //to store external branch lengths
    DoubleVector outer;
    //to store all branch lengths
    DoubleVector all;
    //convert the tree into split graph (vector of split*)
    SplitGraph mySg;
    convertSplits(mySg);
    //get information about the branch length based on this SplitGraph
    for ( SplitGraph::iterator it = mySg.begin(); it != mySg.end(); it++)
    {
        (*it)->report(cout);
        //the split is an external branch
        if ( (*it)->countTaxa() == 1 )
            outer.push_back((*it)->getWeight());
        else //the split is an internal branch
            inner.push_back((*it)->getWeight());
        //a branch
        all.push_back((*it)->getWeight());
    }
    //sort the three vectors of branch lengths
    sort(inner.begin(),inner.end());
    sort(outer.begin(),outer.end());
    sort(all.begin(),all.end());
    //for the statistics
    int noInner = inner.size();
    int noOuter = outer.size();
    int noBr = all.size();
    double aveInner = 0;
    double aveOuter = 0;
    double treeLen = 0;

    for ( int i = 0; i < noInner; i++ )
        aveInner += inner[i];
    for ( int i = 0; i < noOuter; i++ )
        aveOuter += outer[i];
    for ( int i = 0; i < noBr; i++ )
        treeLen += all[i];
    aveInner /= (double)noInner;
    aveOuter /= (double)noOuter;
    out << "minInter maxInter aveInter minExter maxExter aveExter minBr maxBr treeLen noBr" << endl;
    out << inner[0] << " " << inner[noInner-1] << " " << aveInner << " " << outer[0] << " " << outer[noOuter-1] << " " << aveOuter << " " << all[0] << " " << all[noBr-1] << " " << treeLen << " " << noBr << endl;
}

void MaTree::comparedTo (MTreeSet &trees, DoubleMatrix &brLenMatrix, IntVector &RFs, DoubleVector &BSDs) {
    //for consistency reason
    NodeVector taxa;
    getTaxa(taxa);
    sort(taxa.begin(), taxa.end(), nodenamecmp);
    int i;
    NodeVector::iterator it;
    for (it = taxa.begin(), i = 0; it != taxa.end(); it++, i++)
        (*it)->id = i;

    //convert the tree into SplitIntMap
    SplitIntMap sim;
    Split *sp = new Split(leafNum);
    convertSplitIntMap(sim, sp, 0);
    //output to test
    /*	for ( SplitIntMap::iterator it = sim.begin(); it != sim.end(); it++ ){
    		cout << (*it).second << "\t";
    		(*it).first->report(cout);
    	}*/

    // get the taxa name
    vector<string> taxname;
    taxname.resize(leafNum);
    getTaxaName(taxname);

    int noTree = trees.size();
    if (noTree == 0 ) return;
    RFs.resize(noTree);
    BSDs.resize(noTree);

    //now check if it is consistent (rooting, same leaf set) with the input trees
    MTree *tree = trees.front();
//	if (tree->rooted != rooted)
//		outError("Rooted and unrooted trees are mixed up");
    if (tree->leafNum != leafNum)
        outError("Tree has different number of taxa!");
    vector<string> taxname1;
    taxname1.resize(leafNum);
    tree->getTaxaName(taxname1);

    vector<string>::iterator strit;
    for (strit = taxname1.begin(), i = 0; strit != taxname1.end(); strit++, i++) {
        if ((*strit) != taxname[i])
            outError("Tree has different taxa names!");
    }

    MTreeSet::iterator tit;
    for ( tit = trees.begin(), i=0; tit != trees.end(); tit++, i++ )
    {
        DoubleVector brVec(nodeNum,-2);
        SplitGraph *sg = new SplitGraph;
        SplitIntMap *hs = new SplitIntMap;
        (*tit)->convertSplits(taxname,*sg);
        // make sure that taxon 0 is included
        for (SplitGraph::iterator sit = sg->begin(); sit != sg->end(); sit++) {
            if (!(*sit)->containTaxon(0)) (*sit)->invert();
            //	(*sit)->report(cout);
            hs->insertSplit((*sit), 1);
        }

        int rf = 0;
        double bsd = 0;
        //go through each split in this tree (not the compared tree)
        for ( SplitIntMap::iterator tsit = sim.begin(); tsit != sim.end(); tsit++ )
        {
            Split* fSplit = hs->findSplit(tsit->first); // check whether the compared tree contains this split
            if (fSplit) { //yes
                brVec[tsit->second] = fSplit->getWeight(); //update brVec
                bsd += (fSplit->getWeight() - tsit->first->getWeight()) *  (fSplit->getWeight() - tsit->first->getWeight());       //update bsd
            }
            else {
                brVec[tsit->second] = -1;
                rf++;
                bsd += tsit->first->getWeight() * tsit->first->getWeight();
            }
        }
        //go through each split in the compared tree
        for ( SplitIntMap::iterator fsit = hs->begin(); fsit != hs->end(); fsit++ )
        {
            Split* fSplit = sim.findSplit(fsit->first);
            if (!fSplit) {
                rf++;
                bsd += fsit->first->getWeight() * fsit->first->getWeight();
            }
        }
        //insert the result
        RFs[i] = rf;
        BSDs[i] = bsd;
        brLenMatrix.push_back(brVec);
        delete sg;
        delete hs;
    }
}

//void MaTree::convertSplitIntMap(SplitIntMap &sim){}

void MaTree::convertSplitIntMap(SplitIntMap &sim, Split *resp, const int taxonID, Node *node, Node *dad) {
    if (!node) node = root;
    ASSERT(resp->getNTaxa() == leafNum);
    ASSERT (taxonID >= 0 && taxonID < leafNum);
    bool has_child = false;
    FOR_NEIGHBOR_IT(node, dad, it) {
        //vector<int> taxa;
        //getTaxaID((*it)->node, node, taxa);

        Split *sp = new Split(leafNum, (*it)->length);
        convertSplitIntMap(sim, sp, taxonID,(*it)->node, node);
        *resp += *sp;
        if (!sp->containTaxon(taxonID))
            sp->invert();
        //sg.push_back(sp);
        if ( node == root)
            sim.insertSplit(sp,node->id);
        else
            sim.insertSplit(sp,(*it)->node->id);
        has_child = true;
    }
    if (!has_child)
        resp->addTaxon(node->id);

}




