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

#ifndef PARAMETER_TREE_H
#define PARAMETER_TREE_H

#include <string>
#include <sstream>
#include <vector>

class Node {

	public:
						Node(void);
		void			flipCl(void) { (activeCl == 0 ? activeCl = 1 : activeCl = 0); }
		void			flipTi(void) { (activeTi == 0 ? activeTi = 1 : activeTi = 0); }
		Node*			getLft(void) { return lft; }
		Node*			getRht(void) { return rht; }
		Node*			getAnc(void) { return anc; }
		int				getIdx(void) const { return idx; }
		std::string		getName(void) { return name; }
		int				getActiveCl(void) { return activeCl; }
		int				getActiveTi(void) { return activeTi; }
		bool			getIsClDirty(void) { return isClDirty; }
		bool			getIsTiDirty(void) { return isTiDirty; }
		double			getNodeDepth(void) { return nodeDepth; }
		bool			getIsLeaf(void) { return isLeaf; }
		bool			getIsCalibratedDepth(void) { return isCalib; }
		void			setLft(Node *p) { lft = p; }
		void			setRht(Node *p) { rht = p; }
		void			setAnc(Node *p) { anc = p; }
		void			setIdx(int x) { idx = x; }
		void			setName(std::string s) { name = s; }
		void			setActiveCl(int x) { activeCl = x; }
		void			setActiveTi(int x) { activeTi = x; }
		void			setIsClDirty(bool x) { isClDirty = x; }
		void			setIsTiDirty(bool x) { isTiDirty = x; }
		void			setNodeDepth(double x) { nodeDepth = x; }
		void			setIsLeaf(bool x ) { isLeaf = x; }
		void			setIsCalibratedDepth(bool x ) { isCalib = x; }
		double			getBranchTime() { return branchTime; }
		double			getRateGVal() { return rateGVal; }
		double			getUerBL() { return userBL; }
		int				getRateGrpIdx() { return rateGrpIdx; }
		void			setBranchTime(double v) { branchTime = v; }
		void			setRtGrpVal(double v) { rateGVal = v; }
		void			setRtGrpIdx(int i) { rateGrpIdx = i; }
		void			setUerBL(double v) { userBL = v; }
		double			getNodeYngTime() { return youngt; }
		double			getNodeOldTime() { return oldt; }
		void			setNodeYngTime(double v) { youngt = v; }
		void			setNodeOldTime(double v) { oldt = v; }
		int				getNumDecendantTax() { return numDecTax; }
		void			setNumDecendantTax(int i) { numDecTax = i; }
		int				getNodeCalibPrDist() { return nodeCalibPrD; }
		void			setNodeCalibPrDist(int i) { nodeCalibPrD = i; }
		double			getNodeExpCalRate() { return nodeExpCalRate; }
		void			setNodeExpCalRate(double v) { nodeExpCalRate = v; }
		bool			getIsContaminatedFossil() { return taintFossil; }
		void			setIsContaminatedFossil(bool b) { taintFossil = b; }
		int				getRedFlag() { return redFlag; }
		void			setRedFlag(int i) { redFlag = i; }

	private:
		Node			*lft;
		Node			*rht;
		Node			*anc;
		int				idx;
		std::string		name;
		int				activeCl;
		int				activeTi;
		bool			isClDirty;
		bool			isTiDirty;
		double			nodeDepth;
		bool			isLeaf;
		bool			isCalib;
		double			youngt;
		double			oldt;
		double			branchTime;
		double			rateGVal;
		int				rateGrpIdx;
		double			userBL;
		int				numDecTax;
		int				nodeCalibPrD;
		double			nodeExpCalRate;
		bool			taintFossil;
		int				redFlag;
};

class Alignment;
class Calibration;
class MbRandom;
class Model;
class ExpCalib;
class Tree : public Parameter {

	public:
										Tree(MbRandom *rp, Model *mp, Alignment *ap, std::string ts, 
											 bool ubl, bool allnm, bool rndNods, std::vector<Calibration *> clb, 
											 double rth, bool sb, bool exhpc, ExpCalib *ec);
										~Tree(void); 
		Tree							&operator=(const Tree &t);
		void							clone(const Tree &t);
		void							getDownPassSequence(void);
		Node*							getDownPassNode(int i) { return downPassSequence[i]; }
		Node*							getNodeByIndex(int i) { return &nodes[i]; }
		int								getNumNodes(void) { return numNodes; }
		int								getNumTaxa(void) { return numTaxa; }
		Node*							getRoot(void) { return root; }
		double							update(double &oldLnL);
		double							updateOneNode();
		double							updateAllNodes(double &oldLnL);
		double							updateAllNodesRnd(double &oldLnL);
		double							lnPrior();
		double							lnPriorRatio(double nh, double oh);
		double							lnCalibPriorRatio(double nh, double oh, double lb, double ub);
		double							lnExpCalibPriorRatio(double nh, double oh, double offSt, double expRate);
		void							print(std::ostream & o) const;
		void							flipAllCls(void);
		void							flipAllTis(void);
		void							upDateAllCls(void);
		void							upDateAllTis(void);
		void							flipToRootClsTis(Node *p);
		void							updateToRootClsTis(Node *p);
		void							flipToRootClsTis(int ndID);
		void							updateToRootClsTis(int ndID);
		std::string						getTreeDescription(void);
		std::string						getFigTreeDescription(void);
		std::string						getCalibInitialTree(void);
		std::string						writeParam(void);
		std::string						getNodeInfoNames(void);
		std::string						getNodeInfoList(void);
		std::string						getDownPNodeInfoNames(void);
		std::string						getDownPNodeInfoList(void);
		std::string						getCalNodeInfoNames(void);
		std::string						getCalNodeInfoList(void);
		void							setRootRateValue(double v) { root->setRtGrpVal(v); }
		void							setAllNodeBranchTimes(void);
		void							setRndShufNdMv(bool b) { randShufNdMv = b; }
		double							getTreeScale() { return treeScale; }
		void							setTreeScale(double s) { treeScale = s; }
		void							verifyTreeDebug(int iter, std::string pn);
		bool							getIsCalibratedTree() { return isCalibTree; }
		double							getTreeCBDNodePriorProb();
		double							getTreeCBDNodePriorProb(double netDiv, double relDeath);
		double							getTreeSpeciationProbability();
		double							getSumOfNodeHeights();
		void							setNodeRateValues();
		std::vector<Node *>				getListOfCalibratedNodes();
		double							getRootCalibExpRate() { return root->getNodeExpCalRate(); }
		void							checkNodeCalibrationCompatibility();
		int								checkTreeForCalibrationCompatibility();
		void							zeroNodeRedFlags();
							
	private:
		void							buildTreeFromNewickDescription(std::string ts);
		void							initializeNodeDepthsFromUserBL(void);
		void							initializeNodeDepths(void);
		void							initializeCalibratedNodeDepths(void);
		std::vector<double>				recursiveNodeDepthInitialization(Node *p, int &nCont, double maxD);
		void							adjustNodesCompatibleWCalabrations(void);
		int								setNodesNumberDecendantTaxa(Node *p);
		static int						dex(const Node *p);
		bool							isValidChar(char c);
		void							passDown(Node *p, int *x);
		void							showNodes(Node *p, int indent, std::ostream &ss) const;
		void							writeTree(Node *p, std::stringstream &ss);
		void							writeFigTree(Node *p, std::stringstream &ss);
		void							writeCalibrationFigTree(Node *p, std::stringstream &ss);
		void							setNodeCalibrationPriors(ExpCalib *ec);
		int								findCalibNode(std::string t1, std::string t2);
		int								findTaxLftRht(Node *p, std::string t1, std::string t2, int &setNd);
		void							getTreeDotFormat(int ngen, std::string pn);
		double							getTemporaryNodeMaxBound(Node *p);
		void							assignMixLambdaHyperPrToNode(Node *p);
		double							getAMixLambdaHyperPrToNode();
		

		Alignment						*alignmentPtr;
		int								numTaxa;
		int								numNodes;
		Node							*nodes;
		Node							*root;
		Node							**downPassSequence;
		bool							useInputBLs;
		bool							moveAllNodes;
		bool							randShufNdMv;
		std::vector<Calibration*>		calibNds;
		double							treeScale;
		bool							isCalibTree;
		int								treeTimePrior;
		bool							softBounds;
		bool							expHyperPrCal;
};

#endif
