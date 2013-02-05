/* 
 * DPPDiv version 1.0b source code (git: 9c0ac3d2258f89827cfe9ba2b5038f0f656b82c1)
 * Copyright 2009-2011
 * Tracy Heath(1,2,3) (NSF postdoctoral fellowship in biological informatics DBI-0805631)
 * Mark Holder(1)
 * John Huelsenbeck(2)
 *
 * (1) Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS 66045
 * (2) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (3) email: tracyh@berkeley.edu
 *
 * DPPDiv is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * Some of this code is from publicly available source by John Huelsenbeck
 */

#include "Alignment.h"
#include "Calibration.h"
#include "MbMath.h"
#include "MbRandom.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_expcalib.h"
#include "Parameter_rate.h"
#include "Parameter_speciaton.h"
#include "Parameter_tree.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

Node::Node(void) {

	lft = NULL;
	rht = NULL;
	anc = NULL;
	idx = 0;
	name = "";
	activeCl = 0;
	activeTi = 0;
	isClDirty = true;
	isTiDirty = true;
	nodeDepth = 0.0;
	isLeaf = false;
	isCalib = false;
	youngt = -1.0;
	oldt = -1.0;
	branchTime = 0.0;
	rateGVal = 0.0;
	rateGrpIdx = -1;
	userBL = 0.0;
	nodeCalibPrD = 1;
	nodeExpCalRate = -1.0;
	taintFossil = false;
	redFlag = 0;
}

Tree::Tree(MbRandom *rp, Model *mp, Alignment *ap, string ts, bool ubl, bool allnm, 
		   bool rndNods, vector<Calibration *> clb, double rth, bool sb, bool exhpc,
		   ExpCalib *ec) : Parameter(rp, mp) {

	alignmentPtr = ap;
	numTaxa = 0;
	numNodes = 0;
	nodes = NULL;
	root = NULL;
	downPassSequence = NULL;
	useInputBLs = ubl;
	moveAllNodes = allnm;
	calibNds = clb;
	treeScale = rth;
	treeTimePrior = modelPtr->getTreeTimePriorNum();
	name = "TR";
	randShufNdMv = rndNods;
	softBounds = sb;
	isCalibTree = false;
	expHyperPrCal = exhpc;
	buildTreeFromNewickDescription(ts); 
	if(calibNds.size() > 0){
		isCalibTree = true;
		setNodeCalibrationPriors(ec);
		initializeCalibratedNodeDepths();
		while(checkTreeForCalibrationCompatibility() > 0){
			zeroNodeRedFlags();
			initializeCalibratedNodeDepths();
		}
	}
	else{
		isCalibTree = false;
		if(useInputBLs)
			initializeNodeDepthsFromUserBL(); 
		else
			initializeNodeDepths();
	}
	setAllNodeBranchTimes();
	
}

Tree::~Tree(void) {

	delete [] nodes;
	delete [] downPassSequence;
}

void Tree::buildTreeFromNewickDescription(string ts) {

	numTaxa = alignmentPtr->getNumTaxa();
	numNodes = 2 * numTaxa - 1;
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];
	for (int i=0; i<numNodes; i++)
		nodes[i].setIdx(i);

	bool readingBrlen = false;
	Node *p = NULL, *q = NULL;
	int nextInteriorNode = numTaxa;
	for (unsigned i=0; i<ts.size(); i++)
		{
		char c = ts[i];
		if ( c == '(' )
			{
			q = &nodes[nextInteriorNode++];
			if (p == NULL)
				{
				p = q;
				root = p;
				}
			else
				{
				q->setAnc(p);
				if (p->getLft() == NULL)
					p->setLft(q);
				else if (p->getRht() == NULL)
					p->setRht(q);
				else
					{
					cerr << "ERROR: Problem adding interior node to tree" << endl;
					exit(1);
					}
				}
			p = q;
			readingBrlen = false;
			}
		else if ( c == ')' )
			{
			if (p->getAnc() == NULL)
				{
				cerr << "ERROR: Problem moving down the tree" << endl;
				exit(1);
				}
			else
				p = p->getAnc();
			readingBrlen = false;
			}
		else if ( c == ',' )
			{
			if (p->getAnc() == NULL)
				{
				cerr << "ERROR: Problem moving down the tree" << endl;
				exit(1);
				}
			else
				p = p->getAnc();
			readingBrlen = false;
			}
		else if ( c == ':' )
			{
			readingBrlen = true;
			}
		else if ( c == ';' )
			{
			// we are finished with the tree...check
			}
		else
			{
			string s = "";
			while ( isValidChar(ts[i]) )
				s += ts[i++];
			i--;
			if (readingBrlen == false)
				{
				if ( alignmentPtr->isTaxonPresent(s) == false )
					{
					cerr << "ERROR: Cannot find taxon in alignment" << endl;
					exit(1);
					}
				int indexForTaxon = alignmentPtr->getIndexForTaxonNamed(s);
				q = &nodes[indexForTaxon];
				if (p == NULL)
					{
					cerr << "ERROR: Problem adding a tip to the tree" << endl;
					exit(1);
					}
				else
					{
					q->setAnc(p);
					if (p->getLft() == NULL)
						p->setLft(q);
					else if (p->getRht() == NULL)
						p->setRht(q);
					else
						{
						cerr << "ERROR: Problem adding interior node to tree" << endl;
						exit(1);
						}
					}
				p = q;
				p->setName(s);
				p->setIsLeaf(true);
				}
			else
				{
				double x;
				istringstream buf(s);
				buf >> x;
				p->setUerBL(x);
				}
			}
		}
	getDownPassSequence();
	root->setRtGrpVal(0.5); 
}

Tree& Tree::operator=(const Tree &t) {

	if (this != &t)
		clone(t);
	return *this;
}

void Tree::clone(const Tree &t) {

	if (numNodes != t.numNodes || numTaxa != t.numTaxa)
		{
		cerr << "ERROR: Attempting to clone trees of unequal size" << endl;
		exit(1);
		}
		
	for (int i=0; i<numNodes; i++){
		Node *pTo   = &nodes[i];
		Node *pFrom = &t.nodes[i];
				
		pTo->setName( pFrom->getName() );
		pTo->setActiveCl( pFrom->getActiveCl() );
		pTo->setActiveTi( pFrom->getActiveTi() );
		pTo->setIsClDirty( pFrom->getIsClDirty() );
		pTo->setIsTiDirty( pFrom->getIsTiDirty() );
		pTo->setNodeDepth( pFrom->getNodeDepth() );
		pTo->setIsLeaf( pFrom->getIsLeaf() );
		pTo->setRtGrpVal( pFrom->getRateGVal() );
		pTo->setBranchTime( pFrom->getBranchTime() ); 
		pTo->setIsCalibratedDepth( pFrom->getIsCalibratedDepth() );
		pTo->setNodeYngTime( pFrom->getNodeYngTime() );
		pTo->setNodeOldTime( pFrom->getNodeOldTime() );
		pTo->setRtGrpIdx( pFrom->getRateGrpIdx() );
		pTo->setNodeCalibPrDist( pFrom->getNodeCalibPrDist() );
		pTo->setNumDecendantTax( pFrom->getNumDecendantTax() );
		pTo->setNodeExpCalRate( pFrom->getNodeExpCalRate() );
		
		
		if (pFrom->getLft() == NULL)
			pTo->setLft(NULL);
		else
			pTo->setLft( &nodes[pFrom->getLft()->getIdx()] );
		
		if (pFrom->getRht() == NULL)
			pTo->setRht(NULL);
		else
			pTo->setRht( &nodes[pFrom->getRht()->getIdx()] );

		if (pFrom->getAnc() == NULL)
			pTo->setAnc(NULL);
		else
			pTo->setAnc( &nodes[pFrom->getAnc()->getIdx()] );
	}
		
	root = &nodes[t.root->getIdx()];
	treeScale = t.treeScale;
	treeTimePrior = t.treeTimePrior;
		
	for (int i=0; i<numNodes; i++)
		downPassSequence[i] = &nodes[ t.downPassSequence[i]->getIdx() ];
}

int Tree::dex(const Node *p) {
	return (p == NULL ? -1 : p->getIdx());
}

void Tree::getDownPassSequence(void) {

	int x = 0;
	passDown(root, &x);
}

void Tree::initializeNodeDepthsFromUserBL(void) {
	
	double inDepth = 0.0;
	int k = 0;
	Node *p = &nodes[k];
	while(!p->getIsLeaf())
		p = &nodes[k++];
	while(p != root){
		inDepth += p->getUerBL();
		p = p->getAnc();
	}
	double scaler = 1.0 / inDepth;
	Node *q = NULL;
	double nodeDepth = 0.0;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsLeaf())
			p->setNodeDepth(0.0);
		else if(p == root) 
			p->setNodeDepth(1.0);
		else {
			nodeDepth = p->getUerBL();
			q = p->getAnc();
			while(q != root){
				nodeDepth += q->getUerBL();
				q = q->getAnc();
			}
			p->setNodeDepth(1.0 - (nodeDepth * scaler));
		}
	}
	
}

void Tree::initializeNodeDepths(void) {

	vector<Node*> potentialSplits;
	vector<double> nodeTimes;
	
	for (int i=0; i<numTaxa-2; i++)
		nodeTimes.push_back( ranPtr->uniformRv() );
	sort( nodeTimes.begin(), nodeTimes.end(), greater<double>() );

	for (unsigned i=0; i<nodeTimes.size(); i++)
		cout << nodeTimes[i] << " ";
	cout << endl;

	root->setNodeDepth(1.0);
	if (root->getLft()->getIsLeaf() == false)
		potentialSplits.push_back(root->getLft());
	if (root->getRht()->getIsLeaf() == false)
		potentialSplits.push_back(root->getRht());
	
	int nextTime = 0;
	while (potentialSplits.size() > 0)
		{
		Node *p = potentialSplits[(int)(ranPtr->uniformRv()*potentialSplits.size())];
		p->setNodeDepth(nodeTimes[nextTime++]);
		for (vector<Node *>::iterator n=potentialSplits.begin(); n != potentialSplits.end(); n++)
			{
			if ( (*n) == p )
				{
				potentialSplits.erase( n );
				break;
				}
			}
		if (p->getLft()->getIsLeaf() == false)
			potentialSplits.push_back(p->getLft());
		if (p->getRht()->getIsLeaf() == false)
			potentialSplits.push_back(p->getRht());
		}
}

void Tree::initializeCalibratedNodeDepths(void) {
	root->setNodeDepth(1.0);
	bool goodinit = false;
	while(!goodinit){
		for(vector<Calibration *>::iterator v = calibNds.begin(); v != calibNds.end(); v++){
			Node *p = &nodes[(*v)->getNodeIndex()];
			if(p != root){
				Node *anc, *ldes, *rdes;
				anc = p->getAnc();
				ldes = p->getLft();
				rdes = p->getRht();
				double ltime = ldes->getNodeDepth();
				double rtime = rdes->getNodeDepth();
				double atime = anc->getNodeDepth();
				double oldbound = 0.0;
				double yngbound = p->getNodeYngTime() / treeScale;
				if(p->getNodeCalibPrDist() == 1)
					oldbound = p->getNodeOldTime() / treeScale;
				else
					oldbound = getTemporaryNodeMaxBound(p) / treeScale;
				if(p->getNodeOldTime() == p->getNodeYngTime()){
					double newNodeDepth = yngbound;
					p->setNodeDepth(newNodeDepth);
					goodinit = true;
				}
				else{
					if(ltime > 0.0 || rtime > 0.0){
						if(rtime > ltime){
							if(rtime > yngbound)
								yngbound = rtime;
						}
						else{
							if(ltime > yngbound)
								yngbound = ltime;
						}
					}
					if(atime > 0.0){
						if(atime < oldbound)
							oldbound = atime;
					}
					if(oldbound > root->getNodeDepth())
						oldbound = root->getNodeDepth();
					if(yngbound > oldbound){
						goodinit = false;
						break;
					}
					else{
						double newNodeDepth = yngbound + ranPtr->uniformRv()*(oldbound-yngbound);
						p->setNodeDepth(newNodeDepth);
						goodinit = true;
					}
				}
			}
			else
				goodinit = true;
		}
	}
	setNodesNumberDecendantTaxa(root);
	int nnc = 0;
	recursiveNodeDepthInitialization(root, nnc, 1.0);
}

double Tree::getTemporaryNodeMaxBound(Node *p){
	
	double tempmax = treeScale;
	Node *q = p->getAnc();
	while(q != root && tempmax == treeScale){
		if(q->getIsCalibratedDepth())
			tempmax = q->getNodeYngTime();
		else
			q = q->getAnc();
	}
	return tempmax;
}


vector<double> Tree::recursiveNodeDepthInitialization(Node *p, int &nCont, double maxD) {
		
	vector<double> ndTimesC;
	if(p->getIsLeaf()){
		for(int i=0; i< nCont; i++)
			ndTimesC.push_back( ranPtr->uniformRv()*(maxD) );
		if(ndTimesC.size() > 1)
			sort( ndTimesC.begin(), ndTimesC.end() );
		return ndTimesC;
	}
	else{
		bool getDFrmVec = false;
		int nAncN = nCont;
		double myDepth = maxD;
		if(p->getIsCalibratedDepth()){
			myDepth = p->getNodeDepth();
			nCont = 0;
			getDFrmVec = false;
		}
		else{
			nCont += 1;
			getDFrmVec = true;
		}
	
		Node *d1 = p->getLft();
		Node *d2 = p->getRht();
		if(p->getLft()->getNumDecendantTax() < p->getRht()->getNumDecendantTax()){
			if(p->getRht()->getIsCalibratedDepth() == false && p->getLft()->getIsCalibratedDepth() == true){
				d1 = p->getLft();
				d2 = p->getRht();
			}
			else{
				d2 = p->getLft();
				d1 = p->getRht();
			}
		}
		else{
			if(p->getRht()->getIsCalibratedDepth() == true && p->getLft()->getIsCalibratedDepth() == false){
				d2 = p->getLft();
				d1 = p->getRht();
			}
			else{
				d1 = p->getLft();
				d2 = p->getRht();
			}
		}
		ndTimesC = recursiveNodeDepthInitialization(d1, nCont, myDepth);
		if(getDFrmVec){
			if(p->getNodeDepth() > 0.0){
				myDepth = p->getNodeDepth();
				if(ndTimesC.size() > 0)
					ndTimesC.erase(ndTimesC.begin()+0);
			}
			else{
				if(ndTimesC.size() == 0){
					int myID = p->getIdx();
					cerr << "ERROR: node vector is empty at Node " << myID << endl;
					exit(1);
				}
				myDepth = ndTimesC[0];
				p->setNodeDepth(myDepth);
				ndTimesC.erase(ndTimesC.begin()+0);
			}
		}
		nCont = 0;
		recursiveNodeDepthInitialization(d2, nCont, myDepth);
		
		if(getDFrmVec)
			return ndTimesC;
		else{
			if(ndTimesC.size() > 0){
				int myID = p->getIdx();
				cout << "Calibrated " << myID << "  (MAX = " << maxD*treeScale; 
				cout << ")  (MIN = " << myDepth*treeScale << ") ";
				for(vector<double>::iterator v = ndTimesC.begin(); v != ndTimesC.end(); v++){
					cout << "  " << (*v)*treeScale << "  ";
				}
				cout << endl;
				cerr << "ERROR: node vector has stuff in it" << endl;
				exit(1);
			}
			if(myDepth > maxD){

				double badDep = maxD;
				vector<int> idNFix;
				Node *ancs = p->getAnc();
				while(badDep < myDepth){
					badDep = ancs->getNodeDepth();
					if(badDep < myDepth)
						idNFix.push_back(ancs->getIdx());
					ancs = ancs->getAnc();
				}
				for(int i=0; i<idNFix.size(); i++)
					ndTimesC.push_back( myDepth + ranPtr->uniformRv()*(badDep - myDepth) );
				sort( ndTimesC.begin(), ndTimesC.end() );
				for(int i=0; i<idNFix.size(); i++){
					nodes[idNFix[i]].setNodeDepth(ndTimesC[i]);
				}
				ndTimesC.clear();
				return ndTimesC;

			}
			for(int i=0; i< nAncN; i++)
				ndTimesC.push_back( myDepth + ranPtr->uniformRv()*(maxD - myDepth) );
			sort( ndTimesC.begin(), ndTimesC.end() );
			return ndTimesC;
		}
	}
	return ndTimesC;
}

void Tree::adjustNodesCompatibleWCalabrations(void) {
	
	
}

int Tree::setNodesNumberDecendantTaxa(Node *p) {
	
	if(p->getIsLeaf()){
		p->setNumDecendantTax(0);
		return 1;
	}
	else{
		int nds = setNodesNumberDecendantTaxa(p->getLft());
		nds += setNodesNumberDecendantTaxa(p->getRht());
		p->setNumDecendantTax(nds);
		return nds;
		
	}
	return 0;
}



bool Tree::isValidChar(char c) {

	if ( c == '(' || c == ')' || c == ',' || c == ':' || c == ';' )
		return false;
	return true;
}

void Tree::passDown(Node *p, int *x) {

	if ( p != NULL )
		{
		passDown(p->getLft(), x);
		passDown(p->getRht(), x);
		downPassSequence[(*x)++] = p;
		}
}

void Tree::print(std::ostream & o) const {

	o << "Tree:\n";
	showNodes(root, 3, o);
}

double Tree::update(double &oldLnL) {

	double lppr = 0.0;
	if(moveAllNodes){ 
		if(randShufNdMv)
			lppr = updateAllNodesRnd(oldLnL);
		else lppr = updateAllNodes(oldLnL);
	}	
	else lppr = updateOneNode();
	setAllNodeBranchTimes();
	return lppr;
}

double Tree::updateOneNode() {

	Node *p = NULL;
	do {
		p = downPassSequence[(int)(ranPtr->uniformRv()*numNodes)];
	} while (p->getLft() == NULL || p->getRht() == NULL || p->getAnc() == NULL);
		
	double largestTime = p->getAnc()->getNodeDepth();
	double smallestTime = p->getLft()->getNodeDepth();
	if (p->getRht()->getNodeDepth() > smallestTime)
		smallestTime = p->getRht()->getNodeDepth();
	if(p->getIsCalibratedDepth() && softBounds == false){
		double ycal = p->getNodeYngTime() / treeScale;
		double ocal =  largestTime;
		if(p->getNodeCalibPrDist() == 1)
			ocal = p->getNodeOldTime() / treeScale;
		if(ycal > smallestTime)
			smallestTime = ycal;
		if(ocal < largestTime)
			largestTime = ocal;
	}
	double currDepth = p->getNodeDepth();
	double newNodeDepth = currDepth;
	double lnPrRatio = 0.0;
	if(largestTime > smallestTime){
		newNodeDepth = smallestTime + ranPtr->uniformRv()*(largestTime-smallestTime);
		
		p->setNodeDepth(newNodeDepth);

		flipToRootClsTis(p);
		updateToRootClsTis(p);
		modelPtr->setTiProb();
		if(p->getIsCalibratedDepth()){
			if(softBounds && p->getNodeCalibPrDist() == 1){
				double ycal = p->getNodeYngTime();
				double ocal = p->getNodeOldTime();
				lnPrRatio += lnCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, ycal, ocal);
			}
			if(p->getNodeCalibPrDist() == 2){
				double offst = p->getNodeYngTime();
				double nodeExCR = p->getNodeExpCalRate();
				lnPrRatio += lnExpCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, offst, nodeExCR);
			}
		}
		
	}
	lnPrRatio += lnPriorRatio(newNodeDepth, currDepth);
	return lnPrRatio;
}

double Tree::updateAllNodes(double &oldLnL) {
	
	
	upDateAllCls(); 
	upDateAllTis();
	double oldLike = oldLnL;
	Node *p = NULL;
	for(int i = 0; i<numNodes; i++){
		p = downPassSequence[i];
		if(p != root && !p->getIsLeaf()){
			double currDepth = p->getNodeDepth();
			double largestTime = p->getAnc()->getNodeDepth();
			double smallestTime = p->getLft()->getNodeDepth();
			if (p->getRht()->getNodeDepth() > smallestTime)
				smallestTime = p->getRht()->getNodeDepth();
			
			if(p->getIsCalibratedDepth() && softBounds == false){
				double ycal = p->getNodeYngTime() / treeScale;
				double ocal =  largestTime;
				if(p->getNodeCalibPrDist() == 1)
					ocal = p->getNodeOldTime() / treeScale;
				if(ycal > smallestTime)
					smallestTime = ycal;
				if(ocal < largestTime)
					largestTime = ocal;
			}
			if (largestTime > smallestTime){
				double newNodeDepth = smallestTime + ranPtr->uniformRv()*(largestTime-smallestTime);
				p->setNodeDepth(newNodeDepth);
				flipToRootClsTis(p);
				updateToRootClsTis(p);
				modelPtr->setTiProb();
				double newLnl = modelPtr->lnLikelihood();
				double lnLRatio = newLnl - oldLike;
				double lnPrRatio = lnPriorRatio(newNodeDepth, currDepth);
				if(p->getIsCalibratedDepth()){
					if(softBounds && p->getNodeCalibPrDist() == 1){
						double ycal = p->getNodeYngTime();
						double ocal = p->getNodeOldTime();
						lnPrRatio += lnCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, ycal, ocal);
					}
					if(p->getNodeCalibPrDist() == 2){
						double offst = p->getNodeYngTime();
						double nodeExCR = p->getNodeExpCalRate();
						lnPrRatio += lnExpCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, offst, nodeExCR);
					}
				}
				double lnR = lnPrRatio + lnLRatio + 0.0;
				double r = modelPtr->safeExponentiation(lnR);
				
				if(ranPtr->uniformRv() < r)
					oldLike = newLnl;
				else{
					p->setNodeDepth(currDepth);
					flipToRootClsTis(p);
					updateToRootClsTis(p);
				}
			}
		}
	}
	oldLnL = oldLike;
	return 0.0;
}

double Tree::updateAllNodesRnd(double &oldLnL) {
	
	upDateAllCls();
	upDateAllTis();
	double oldLike = oldLnL;
	Node *p = NULL;
	vector<int> rndNodeIDs;
	for(int i=0; i<numNodes; i++)
		rndNodeIDs.push_back(i);
	random_shuffle(rndNodeIDs.begin(), rndNodeIDs.end());
	for(vector<int>::iterator it=rndNodeIDs.begin(); it!=rndNodeIDs.end(); it++){
		p = downPassSequence[(*it)];
		if(p != root && !p->getIsLeaf()){
			double currDepth = p->getNodeDepth();
			double largestTime = p->getAnc()->getNodeDepth();
			double smallestTime = p->getLft()->getNodeDepth();
			if (p->getRht()->getNodeDepth() > smallestTime)
				smallestTime = p->getRht()->getNodeDepth();
			if(p->getIsCalibratedDepth() && softBounds == false){
				double ycal = p->getNodeYngTime() / treeScale;
				double ocal =  largestTime;
				if(p->getNodeCalibPrDist() == 1)
					ocal = p->getNodeOldTime() / treeScale;
				if(ycal > smallestTime)
					smallestTime = ycal;
				if(ocal < largestTime)
					largestTime = ocal;
			}
			if (largestTime > smallestTime){
				double newNodeDepth = smallestTime + ranPtr->uniformRv()*(largestTime-smallestTime);
				p->setNodeDepth(newNodeDepth);
				flipToRootClsTis(p);
				updateToRootClsTis(p);
				modelPtr->setTiProb();
				double newLnl = modelPtr->lnLikelihood();
				double lnLRatio = newLnl - oldLike;
				double lnPrRatio = lnPriorRatio(newNodeDepth, currDepth);
				if(p->getIsCalibratedDepth()){
					if(softBounds && p->getNodeCalibPrDist() == 1){
						double ycal = p->getNodeYngTime();
						double ocal = p->getNodeOldTime();
						lnPrRatio += lnCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, ycal, ocal);
					}
					if(p->getNodeCalibPrDist() == 2){
						double offst = p->getNodeYngTime();
						double nodeExCR = p->getNodeExpCalRate();
						lnPrRatio += lnExpCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, offst, nodeExCR);
					}
				}
				
				double lnR = lnPrRatio + lnLRatio + 0.0;
				double r = modelPtr->safeExponentiation(lnR);
				
				if(ranPtr->uniformRv() < r){
					oldLike = newLnl;
				}
				else{
					p->setNodeDepth(currDepth);
					flipToRootClsTis(p);
					updateToRootClsTis(p);
				}
			}
		}
	}
	oldLnL = oldLike;	
	return 0.0;
}


void Tree::setAllNodeBranchTimes(void) {
	
	for (int n=0; n<numNodes; n++){
		Node *p = &nodes[n];
		if(p != root){
			double branchtime = (p->getAnc()->getNodeDepth() - p->getNodeDepth()) * treeScale;
			if(branchtime < 0){
				int myID = p->getIdx();
				cerr << "ERROR: The tree has a negative branch length! " << branchtime << " At Node " << myID << endl;
				exit(1);
			}
			else
				p->setBranchTime(branchtime);
		}
	}
}

double Tree::lnPrior() {
	
	return 0.0;
}

double Tree::lnPriorRatio(double nh, double oh) {
		
	nh = nh*treeScale;
	oh = oh*treeScale;
	
	if(treeTimePrior == 1) 
		return 0.0;
	else if(treeTimePrior == 2){
		double diff = modelPtr->getActiveSpeciation()->getNetDiversification();
		double nator = (-(diff)*nh);
		double dator = (-(diff)*oh);		
		return nator - dator;
	}
	else{
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	// br-dr
		double rel = s->getRelativeDeath();			// dr / br
		
		double zn = log(1 - (rel) * exp(-(diff)*nh));
		double nator = -2 * zn + (-(diff)*nh);
		
		double zd = log(1 - (rel) * exp(-(diff)*oh));
		double dator = -2 * zd + (-(diff)*oh);
		
		return nator - dator;
	}
	return 0.0;
}


double Tree::lnCalibPriorRatio(double nh, double oh, double lb, double ub) {
	
	/*
		Prior ratio for a uniform dist on calibrated nodes with soft bounds
		From Yang and Rannala, MBE (2006)
	*/
	double numr = 0.0;
	double dnom = 0.0;
	double th1 = (0.95 * lb) / (0.025 * (ub - lb));
	double th2 = 0.95 / (0.025 * (ub - lb));
	
	if(nh <= lb)
		numr += log(0.025) + log(th1 / lb) + ((th1 - 1) * log(nh / lb));
	else if(nh > ub)
		numr += log(0.025) + log(th2) + (-th2 * (nh - ub));
	else
		numr += log(0.95) - log(ub - lb);
	
	if(oh <= lb)
		dnom += log(0.025) + log(th1 / lb) + ((th1 - 1) * log(oh / lb));
	else if(oh > ub)
		dnom += log(0.025) + log(th2) + (-th2 * (oh - ub));
	else
		dnom += log(0.95) - log(ub - lb);
	
	return numr - dnom;
}

double Tree::lnExpCalibPriorRatio(double nh, double oh, double offSt, double expRate) {
	
	double numr = 0.0;
	double dnom = 0.0;
	numr = ranPtr->lnExponentialPdf(expRate, nh - offSt);
	dnom = ranPtr->lnExponentialPdf(expRate, oh - offSt);
	return numr - dnom;
}



void Tree::flipAllCls(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].flipCl();
}

void Tree::flipAllTis(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].flipTi();
}

void Tree::upDateAllCls(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].setIsClDirty(true);
}

void Tree::upDateAllTis(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].setIsTiDirty(true);
}

void Tree::flipToRootClsTis(Node *p){
	
	Node *q = p;
	if(!q->getIsLeaf()){
		q->getLft()->flipCl();
		q->getRht()->flipCl();
		q->getLft()->flipTi();
		q->getRht()->flipTi();
	}
	while(q != NULL){
		q->flipCl();
		q->flipTi();
		q = q->getAnc();
	}
}

void Tree::updateToRootClsTis(Node *p){
	
	Node *q = p;
	if(!q->getIsLeaf()){
		q->getLft()->setIsClDirty(true);
		q->getRht()->setIsClDirty(true);
		q->getLft()->setIsTiDirty(true);
		q->getRht()->setIsTiDirty(true);
	}
	while(q != NULL){
		q->setIsClDirty(true);
		q->setIsTiDirty(true);
		q = q->getAnc();
	}
}

void Tree::flipToRootClsTis(int ndID){

	Node *q = &nodes[ndID];
	while(q != NULL){
		q->flipCl();
		q->flipTi();
		q = q->getAnc();
	}
}

void Tree::updateToRootClsTis(int ndID){

	Node *q = &nodes[ndID];
	while(q != NULL){
		q->setIsClDirty(true);
		q->setIsTiDirty(true);
		q = q->getAnc();
	}
}


string Tree::getTreeDescription(void){ 
	
	stringstream ss;
	writeTree(root, ss);
	ss << ";";
	string treDesc = ss.str();
	return treDesc;
}

void Tree::writeTree(Node *p, stringstream &ss){
	
	if (p != NULL){
		if(p->getLft() == NULL){ 
			ss << p->getName();
		}
		else{
			ss << "(";
			writeTree(p->getLft(), ss);
			ss << ":" << p->getLft()->getBranchTime();
			ss << ",";
			writeTree(p->getRht(), ss);
			ss << ":" << p->getRht()->getBranchTime();
			ss << ")";
		}
	}
}

string Tree::getFigTreeDescription(void){ 
	
	stringstream ss;
	writeFigTree(root, ss);
	ss << ";";
	string treDesc = ss.str();
	return treDesc;
}

void Tree::writeFigTree(Node *p, stringstream &ss){
	
	if (p != NULL){
		if(p->getLft() == NULL){ 
			ss << p->getName();
		}
		else{
			ss << "(";
			writeFigTree(p->getLft(), ss);
			ss << "[&rate=" << p->getLft()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getLft()->getRateGrpIdx() << "]";
			ss << ":" << p->getLft()->getBranchTime();
			ss << ",";
			writeFigTree(p->getRht(), ss);
			ss << "[&rate=" << p->getRht()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getRht()->getRateGrpIdx() << "]";
			ss << ":" << p->getRht()->getBranchTime();
			ss << ")";
		}
	}
}

string Tree::getCalibInitialTree(void){ 
	
	stringstream ss;
	writeCalibrationFigTree(root, ss);
	ss << "[&rate=0.0,";
	ss << "rate_cat=0";
	if(root->getIsCalibratedDepth())
		ss << ",cal_range={" << root->getNodeYngTime() << "," << root->getNodeOldTime() << "}]";
	else
		ss << ",cal_range={" << treeScale << "," << treeScale << "}]";
	
	ss << ";";
	string treDesc = ss.str();
	return treDesc;
}

void Tree::writeCalibrationFigTree(Node *p, stringstream &ss){
	
	if (p != NULL){
		if(p->getLft() == NULL){ 
			ss << p->getName();
		}
		else{
			ss << "(";
			writeCalibrationFigTree(p->getLft(), ss);
			ss << "[&rate=" << p->getLft()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getLft()->getRateGrpIdx();
			if(p->getLft()->getIsCalibratedDepth())
				ss << ",cal_range={" << p->getLft()->getNodeYngTime() << "," << p->getLft()->getNodeOldTime() << "}]";
			else
				ss << ",cal_range={" << p->getLft()->getNodeDepth() * treeScale << "," << p->getLft()->getNodeDepth() * treeScale << "}]";
			ss << ":" << p->getLft()->getBranchTime();
			ss << ",";
			writeCalibrationFigTree(p->getRht(), ss);
			ss << "[&rate=" << p->getRht()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getRht()->getRateGrpIdx();
			if(p->getRht()->getIsCalibratedDepth())
				ss << ",cal_range={" << p->getRht()->getNodeYngTime() << "," << p->getRht()->getNodeOldTime() << "}]";
			else
				ss << ",cal_range={" << p->getRht()->getNodeDepth() * treeScale << "," << p->getRht()->getNodeDepth() * treeScale << "}]";
			ss << ":" << p->getRht()->getBranchTime();
			ss << ")";
		}
	}
}


string Tree::writeParam(void){
	
	stringstream ss;
	ss << "Tree:" << endl;
	showNodes(root, 3, ss);
	string outp = ss.str();
	return outp;
}

void Tree::showNodes(Node *p, int indent, std::ostream &ss) const {
	
	if (p != NULL){
		for (int i=0; i<indent; i++)
			ss << " ";
		ss << dex(p) << " (" << dex(p->getLft()) << ", " << dex(p->getRht()) << ", " << dex(p->getAnc()) 
		   << ") -- " << fixed << setprecision(5) << (p->getNodeDepth() * treeScale) << " -- ";
		if (p->getLft() == NULL && p->getRht() == NULL )
			ss << " (" << p->getName() << ") ";
		if (p == root)
			ss << " <- Root";
		ss << endl;
		showNodes (p->getLft(),  indent + 2, ss);
		showNodes (p->getRht(), indent + 2, ss);
	}
}

string Tree::getNodeInfoNames(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsLeaf() == false){
			if(p == root)
				ss << "\tRootDepth.Time(N" << p->getIdx() << ")";
			else
				ss << "\tTime(N" << p->getIdx() << ")";
		}
	}
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p != root)
			ss << "\tSRate(N" << p->getIdx() << ")";
	}
	string ni = ss.str();
	return ni;
}


string Tree::getNodeInfoList(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsLeaf() == false)
			ss << "\t" << p->getNodeDepth() * treeScale;
	}
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p != root)
			ss << "\t" << p->getRateGVal();
	}
	string ni = ss.str();
	return ni;
}

string Tree::getDownPNodeInfoNames(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
		if(p->getIsLeaf() == false){
			if(p == root)
				ss << "\tRootDepth.Time(DP" << i << "|N" << p->getIdx() << ")";
			else
				ss << "\tTime(DP" << i << "|N" << p->getIdx() << ")";
		}
	}
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
		if(p != root)
			ss << "\tSRate(DP" << i << "|N" << p->getIdx() << ")";
	}
	string ni = ss.str();
	return ni;
}


string Tree::getDownPNodeInfoList(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
		if(p->getIsLeaf() == false)
			ss << "\t" << p->getNodeDepth() * treeScale;
	}
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
		if(p != root)
			ss << "\t" << p->getRateGVal();
	}
	string ni = ss.str();
	return ni;
}

string Tree::getCalNodeInfoNames(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsCalibratedDepth()){
			if(p == root)
				ss << "\tRoot.expHP(N" << p->getIdx() << ")";
			else
				ss << "\texpHP(N" << p->getIdx() << ")";
		}
	}
	string ni = ss.str();
	return ni;
}


string Tree::getCalNodeInfoList(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsCalibratedDepth())
			ss << "\t" << p->getNodeExpCalRate();
	}
	string ni = ss.str();
	return ni;
}


void Tree::setNodeCalibrationPriors(ExpCalib *ec) {
	

	for(vector<Calibration *>::iterator v = calibNds.begin(); v != calibNds.end(); v++){
		int calbNo = -1;
		calbNo = findCalibNode((*v)->getTxN1(), (*v)->getTxN2());
		Node *p = &nodes[calbNo];
		int calDistrib = (*v)->getPriorDistributionType();
		p->setNodeCalibPrDist(calDistrib);
		
		if(calDistrib == 1){
			p->setNodeYngTime((*v)->getYngTime());
			p->setNodeOldTime((*v)->getOldTime());
		}
		else if(calDistrib == 2){
			if(expHyperPrCal){
				double lv = ec->getLambdaForNode();
				p->setNodeYngTime((*v)->getYngTime());
				p->setNodeOldTime(treeScale * 1.1);
				p->setNodeExpCalRate(lv);
				p->setIsContaminatedFossil(ec->getIsLambdaContaminationClass(lv));
			}
			else{
				p->setNodeYngTime((*v)->getYngTime());
				p->setNodeOldTime(treeScale * 1.1);
				p->setNodeExpCalRate((*v)->getCalExponRate());
			}
		}
		(*v)->setNodeIndex(calbNo);
	}

	for(vector<Calibration *>::iterator v = calibNds.begin(); v != calibNds.end(); v++){
		Node *p = &nodes[(*v)->getNodeIndex()];
		int myNdPrD = p->getNodeCalibPrDist();
		if(p != root){
			if(p->getLft()->getIsCalibratedDepth()){
				int yourNdPrD = p->getLft()->getNodeCalibPrDist();
				double pOldT = p->getNodeOldTime(), lfOldT = p->getLft()->getNodeOldTime();
				if(lfOldT > pOldT && myNdPrD == 1 && yourNdPrD == 1){
					cerr << "ERROR: problem with node calibration descendant old bound > anc old bound" << endl;
					cerr << "Node " << p->getIdx() << " conflicts with it's left descendant, node " 
						 << p->getLft()->getIdx() << endl;
					exit(1);
				}
			}
			if(p->getRht()->getIsCalibratedDepth()){
				int yourNdPrD = p->getRht()->getNodeCalibPrDist();
				double pOldT = p->getNodeOldTime(), rtOldT = p->getRht()->getNodeOldTime();
				if(rtOldT > pOldT && myNdPrD == 1 && yourNdPrD == 1){
					cerr << "ERROR: problem with node calibration descendant old bound > anc old bound" << endl;
					cerr << "Node " << p->getIdx() << " conflicts with it's right descendant, node " 
						 << p->getRht()->getIdx() << endl;
					exit(1);
				}
			}
			if(p->getAnc()->getIsCalibratedDepth()){
				double pYngT = p->getNodeYngTime(), anYngT = p->getAnc()->getNodeYngTime();
				if(anYngT < pYngT && myNdPrD == 1){
					cerr << "ERROR: problem with node calibration descendant young bound > anc young bound" << endl;
					cerr << "Node " << p->getIdx() << " conflicts with it's ancestor, node " 
						 << p->getAnc()->getIdx() << endl;
					exit(1);
				}
			}
		}
	}
}

int Tree::findCalibNode(std::string t1, std::string t2) {
	
	if(t1 == "root"){
		root->setIsCalibratedDepth(true);
		return root->getIdx();
	}
	int nidn = -1;
	int chk;
	chk = findTaxLftRht(root, t1, t2, nidn);
	if(chk != 3){
		cerr << "ERROR: problem finding calibration node (findTaxLftRht)" << endl;
		cout << t1 << "\t" << t2 << endl;
		exit(1);
	}
	if(nidn < 0){
		cerr << "ERROR: problem finding calibration node (setting nidn)" << endl;
		exit(1);
	}
	return nidn;
}

int Tree::findTaxLftRht(Node *p, std::string t1, std::string t2, int &setNd) {
	
	if(p->getIsLeaf()){
		string tn = p->getName();
		if(tn == t1)
			return 1;
		else if(tn == t2)
			return 2;
		else
			return 0;
	}
	else{
		int ld = 0;
		int rd = 0;
		ld = findTaxLftRht(p->getLft(), t1, t2, setNd);
		rd = findTaxLftRht(p->getRht(), t1, t2, setNd);
		if(ld > 0 && rd > 0){
			setNd = p->getIdx();
			if(p->getIsCalibratedDepth() == true){
				cerr << "ERROR: This node (ID = " << setNd << ") has already been calibrated. WTF?" << endl;
				exit(1);
			}
			else
				p->setIsCalibratedDepth(true);
		}
		return ld + rd;
	}
	return 0;
}


void Tree::verifyTreeDebug(int iter, string pn) {
	
	getTreeDotFormat(iter, pn);
}

void Tree::getTreeDotFormat(int ngen, string pn) {
	
	ofstream out;
	out.open("testdot.dbg.dot");
	out << "digraph T {\n   ";
	out << "g [shape=record, label=\"{" << ngen << "|" << pn << "}\"]\n";
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		out << "   n" << p->getIdx();
		out << " [shape=";
		stringstream rl;
		if(p->getIsLeaf()){
			rl << "{";
			if(p->getAnc() == root)
				rl << "anc = ROOT";
			else
				rl << "anc = N" << p->getAnc()->getIdx();
			rl << "|R = " << setprecision(4) << p->getRateGVal() << "|";
			rl << p->getName() << "}";
			out << "record, label=\"" << rl.str() << "\", fillcolor=lightblue, style=filled]\n";
		}
		else{
			if(p == root){
				rl << "{ROOT|H = " << setprecision(4) << treeScale << "}";
				out << "record, style=filled, label=\"" << rl.str() << "\"]\n";
			}
			else{
				rl << "{";
				if(p->getAnc() == root)
					rl << "anc = ROOT";
				else
					rl << "anc = N" << p->getAnc()->getIdx();
				rl << "|{H = " << setprecision(5) << (p->getNodeDepth() * treeScale);
				rl << "|R = " << setprecision(4) << p->getRateGVal() << "}|";
				rl << "N" << p->getIdx() << "}";
				out << "Mrecord, label=\"" << rl.str() << "\", fillcolor=pink, style=filled]\n";
			}
		}
	}
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false){
			out << "   n" << p->getIdx() << " -> n" << p->getLft()->getIdx();
			out << " [label=\"" << setprecision(4) << p->getLft()->getBranchTime()*treeScale << "\"]\n";
			out << "   n" << p->getIdx() << " -> n" << p->getRht()->getIdx();
			out << " [label=\"" << setprecision(4) << p->getRht()->getBranchTime()*treeScale << "\"]\n";
		}
	}
	out << "   {rank=same;";
	for(int i=0; i<numTaxa; i++)
		out << " n" << i << ";";
	out << "}\n";
	out << "}\n";
	out.close();
}

double Tree::getTreeCBDNodePriorProb() {
	
	double nprb = 0.0;
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	// br-dr
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf() == false){
				double nh = p->getNodeDepth() * treeScale;
				double l = (-(diff)*nh); 
				if(p == root){
					l += (-(diff)*nh);
				}
				nprb += l;
			}
		}
		return nprb;
	}
	else{
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	// br-dr
		double rel = s->getRelativeDeath();			// dr / br		
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf() == false){
				double nh = p->getNodeDepth() * treeScale;
				double zn = log(1 - (rel) * exp(-(diff)*nh));
				double l = -2 * zn + (-(diff)*nh); 
				if(p == root){
					l += (-(diff)*nh) - zn;
				}
				nprb += l;
			}
		}
		return nprb;
	}
	return 0.0;
}

double Tree::getTreeCBDNodePriorProb(double netDiv, double relDeath) {
	
	double nprb = 0.0;
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		double diff = netDiv;	// br-dr
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf() == false){
				double nh = p->getNodeDepth() * treeScale;
				double l = (-(diff)*nh); 
				if(p == root){
					l += (-(diff)*nh);
				}
				nprb += l;
			}
		}
		return nprb;
	}
	else{
		double diff = netDiv;	// br-dr
		double rel = relDeath;			// dr / br		
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf() == false){
				double nh = p->getNodeDepth() * treeScale;
				double zn = log(1 - ((rel) * exp(-(diff)*nh)));
				double l = (-2 * zn) + (-(diff)*nh); 
				if(p == root){
					l += (-(diff)*nh) - zn;
				}
				nprb += l;
			}
		}
		return nprb;
	}
	return 0.0;
}

double Tree::getTreeSpeciationProbability() {
	
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	
		double c1 = (numTaxa - 1) * log(diff); 
		double nps = getTreeCBDNodePriorProb(diff, 0.0);
		return c1 + nps;
	}
	else{
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	// br-dr
		double rel = s->getRelativeDeath();
		double lnC = 0.0; 
		double c1 = ((numTaxa - 1) * log(diff)) + (numTaxa * log(1 - rel));
		double nps = getTreeCBDNodePriorProb(diff, rel);
		return lnC + c1 + nps;
	}
	return 0.0;
}


double Tree::getSumOfNodeHeights() {
	
	double sumh = 0.0;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false){
			double nh = p->getNodeDepth() * treeScale;
			sumh += nh;
		}
	}
	return sumh;
}

void Tree::setNodeRateValues() {
	
	NodeRate *nr = modelPtr->getActiveNodeRate();
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p != root){
			int nID = p->getIdx();
			p->setRtGrpVal(nr->getRateForNodeIndexed(nID));
			p->setRtGrpIdx(nr->getTableNumForNodeIndexed(nID));			
		}
		else{
			p->setRtGrpVal(0.0);
			p->setRtGrpIdx(0.0);
		}
	}
}

void Tree::assignMixLambdaHyperPrToNode(Node *p){
	
	ExpCalib *ec = modelPtr->getActiveExpCalib();
	p->setNodeExpCalRate(ec->getLambdaForNode());
	
}

double Tree::getAMixLambdaHyperPrToNode(){
	
	ExpCalib *ec = modelPtr->getActiveExpCalib();
	return ec->getLambdaForNode();
	
}

vector<Node *> Tree::getListOfCalibratedNodes(){
	
	vector<Node *> calibList;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedDepth() && p->getNodeCalibPrDist() == 2)
			calibList.push_back(p);
	}
	return calibList;
}

void Tree::checkNodeCalibrationCompatibility(){

	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedDepth()){
			double lb = p->getNodeYngTime();
			double ub = p->getNodeOldTime();
			double myDepth = p->getNodeDepth() * treeScale;
			if(myDepth > ub || myDepth < lb){
				cerr << "Calibrated node: " << i << " has the depth: " << myDepth << " and is NOT in the range {" << lb << ", " << ub << "}" << endl;
				exit(1);
			}
			else
				cout << "Calibrated node: " << i << " has the depth: " << myDepth << " and is in the range {" << lb << ", " << ub << "}" << endl;
		}
	}
}

void Tree::zeroNodeRedFlags(){
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		p->setRedFlag(0);
		if(p->getIsLeaf() == false && p != root)
			p->setNodeDepth(0.0);
	}
}

int Tree::checkTreeForCalibrationCompatibility(){
	
	int numIncomp = 0;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p != root){
			double branchtime = (p->getAnc()->getNodeDepth() - p->getNodeDepth());
			if(branchtime < 0){
				p->setRedFlag(2);
				numIncomp++;
			}
			if(p->getIsCalibratedDepth() && p->getNodeCalibPrDist() == 1){
				double lb = p->getNodeYngTime();
				double ub = p->getNodeOldTime();
				double myDepth = p->getNodeDepth() * treeScale;
				
				if(myDepth > ub || myDepth < lb){
					p->setRedFlag(1);
					numIncomp++;
				}
			}
		}
	}
	cout << numIncomp << endl;
	return numIncomp;
}
