/*
 * upperbounds.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: olga
 */
#include "upperbounds.h"
#include "phylonode.h"
#include <string.h>
#include "timeutil.h"

void UpperBounds(Params *params, Alignment* alignment, IQTree* tree){

// Output details --------------------------------------------------
// UpperBounds File
	string out_file = params->out_prefix;
	//out_file += ".ub";
	out_file = "results.trueSplits.ub";
	ofstream out;
	out.exceptions(ios::failbit | ios::badbit);
	out.open((char*)out_file.c_str(),std::ofstream::out | std::ofstream::app);

// Details on Split: A|B
	string out_file_split = params->out_prefix;
	//out_file_split += ".split.ub";
	out_file_split = "results.trueSplits.ub.splits";
	ofstream out_split;
	out_split.exceptions(ios::failbit | ios::badbit);
	out_split.open((char*)out_file_split.c_str(),std::ofstream::out | std::ofstream::app);

// Within Family Info: A|B
	string out_file_within = params->out_prefix;
	//out_file_within += ".within.ub";
	out_file_within = "results.within.ub";
	ofstream out_within;
	out_within.exceptions(ios::failbit | ios::badbit);
	out_within.open((char*)out_file_within.c_str(),std::ofstream::out | std::ofstream::app);

// Between Families Info: A|B
	string out_file_between = params->out_prefix;
	//out_file_between += ".between.ub";
	out_file_between = "results.between.ub";
	ofstream out_between;
	out_between.exceptions(ios::failbit | ios::badbit);
	out_between.open((char*)out_file_between.c_str(),std::ofstream::out | std::ofstream::app);

	/* ------------------------------------------------------------------------------------------------------
	 * All output files:
	 * 	out 		-> "results.trueSplits.ub"			-> upper bounds for all splits from an input tree
	 * 	out_split   -> "results.trueSplits.ub.splits"	-> list of 			all splits from an input tree
	 * 	out_within	-> "results.within.ub"				-> comparison of upper bounds within  Split Families
	 * 	out_between	-> "results.between.ub"				-> comparison of upper bounds between Split Families
	 * ------------------------------------------------------------------------------------------------------
	 * FORMAT:
	 * general info (first columns of every file)...........................................................................
	 * 		leafNum		getNSite()	min(taxaA,taxaB)	brLen
	 *
	 * out .................................................................................................................
	 * 		L(A|B)		L(A)L(B)	cN*L(A)L(B)		L(A|B)/L(A)L(B)		L(A|B)/cN*L(A)L(B)	coef
	 *
	 *[4],[5] - the difference between likelihood and UB normalized by likelihood value
	 * if >1, the inequality is true.
	 * if <1, false.
	 *
	 * out_within ..........................................................................................................
	 * 		UB_true		[1..N] UB_random_AB/UB_true (how smaller is the bound for random tree)
	 *
	 * out_between..........................................................................................................
	 * 		UB_true		[1..N] UB_random_CD/UB_true (how smaller is the bound for random tree)
	 *
	 * ------------------------------------------------------------------------------------------------------ */

	int i=0;//, h=0;

	// Printing info about the TreeLogL changes during the tree search
/*	cout<<"mlInitial  = "<<tree->mlInitial<<endl;
	cout<<"mlFirstOpt = "<<tree->mlFirstOpt<<endl;
	cout<<"mlBestTree = "<<tree->getBestScore()<<endl;
	cout<<"mlUnConstr = "<<alignment->computeUnconstrainedLogL()<<endl;*/

	//double mlQuestionary = tree->mlInitial; //or tree->mlFirstOpt for example

	/* ------------------------------------------------------------------------------------------------------
	 * Main PART
	 * ------------------------------------------------------------------------------------------------------ */
	cout<<"Starting Upper Bounds analysis.."<<endl;

	NodeVector branch1, branch2;
	tree->getBranches(branch1, branch2);
	int allSplits = 0;
//	int R=10; // R is the number of random trees we will generate

// A loop over all A|B present on tree T
	for(i = 0; i != branch1.size(); i++){
		vector<int> taxaA, taxaB;
		vector<string> taxaAname, taxaBname;
		tree->getTaxaID(taxaA,branch1[i],branch2[i]);
		tree->getTaxaID(taxaB,branch2[i],branch1[i]);

		/* ------------------------------------------------------------------------------------------------------------
		 * TEST 1: This is the part for tests on [ai/(ai+bi)] and [bi/(ai+bi)] fractions
		 */
		int test1 = 1;
		if(test1 == 1){
			if(taxaA.size() > 3 && taxaB.size() > 3){ // IQTree does not compute lh of tree with less than 4 taxa.
				allSplits++;
				sumFraction(((PhyloNode*) branch1[i]), ((PhyloNode*) branch2[i]), tree);
		}
		}

		/* ------------------------------------------------------------------------------------------------------------
		 * TEST 2: This is the part for tests on random trees and evaluation of Upper Bounds for each split on the input tree
		 */
		int test2 = 0;
		if(test2 == 1){
		if(taxaA.size() > 3 && taxaB.size() > 3){ // IQTree does not compute lh of tree with less than 4 taxa.
			allSplits++;

			// Dealing with subtrees T_A and T_B
			PhyloTree *treeA, *treeB;
			treeA = extractSubtreeUB(taxaA,tree,params,1);
			treeB = extractSubtreeUB(taxaB,tree,params,1);

			// Upper Bound for a given split from the input tree
			double brLen = branch1[i]->findNeighbor(branch2[i])->length;
			double coef  = tree->aln->getNSite()*(log(1+3*exp(-brLen)) - log(1-exp(-brLen)));
			double coef2  = tree->aln->getNSite()*log(1+3*exp(-brLen));
			double UB_true  = coef + treeA->getCurScore() + treeB->getCurScore();
			double UB_true2 = coef2 + treeA->getCurScore() + treeB->getCurScore();

			//cout<<"UB_true = "<<UB_true<<endl;
			out<<tree->leafNum<<"\t"<<tree->aln->getNSite()<<"\t"<<min(taxaA.size(),taxaB.size())<<"\t"<<brLen<<"\t"
					<<tree->getCurScore()<<"\t"<<treeA->getCurScore() + treeB->getCurScore()<<"\t"<<UB_true<<"\t"<<UB_true2<<"\t"
					<<tree->getCurScore()/(treeA->getCurScore() + treeB->getCurScore())<<"\t"
					<<tree->getCurScore()/UB_true<<"\t"<<tree->getCurScore()/UB_true2<<"\t"<<coef<<"\t"<<coef2<<endl;

/*
			// Comparison of Upper Bounds within Split Family ----------------------------------------
			out_within<<tree->leafNum<<"\t"<<tree->aln->getNSite()<<"\t"<<min(taxaA.size(),taxaB.size())<<"\t"<<brLen<<"\t"<<UB_true<<"\t";

			cout<<"comparison within family...."<<endl;
			double UB_random_AB;
			for(j=0; j<30; j++){
				//cout<<"generating "<<j<<" random_AB tree..."<<endl;
				UB_random_AB = RandomTreeAB(tree, treeA, treeB, taxaA, taxaB, params,brLen);
				//cout<<"The upper bound for random tree: "<<UB_random_AB<<endl;
				out_within<<UB_random_AB/UB_true<<"\t";
			}
			out_within<<endl;


			// --------------------------------------------------------------------------------------

			// Comparison of Upper Bounds between Split Families ------------------------------------
			out_between<<tree->leafNum<<"\t"<<tree->aln->getNSite()<<"\t"<<min(taxaA.size(),taxaB.size())<<"\t"<<brLen<<"\t"<<UB_true<<"\t";

			cout<<"comparison between families...."<<endl;
			// creating split C|D which conflicts with A|B
			int n=0;
			n=int(min(taxaA.size(),taxaB.size())/2.);
			//cout<<"taxaA.size() = "<<taxaA.size()<<", taxaB.size() = "<<taxaB.size()<<", n = "<<n<<endl;

			vector<int> taxaC, taxaD;
			PhyloTree *treeC, *treeD;

			// ContraSplit1: changing 1/2 of taxa
			taxaC = taxaA;
			taxaD = taxaB;
			for(h=0; h<n; h++){
				taxaC[h]=taxaB[h];
				taxaD[h]=taxaA[h];
			}
			treeC = extractSubtreeUB(taxaC,tree,params);
			treeD = extractSubtreeUB(taxaD,tree,params);
			double UB_random_CD;
			for(j=0; j<R; j++){
				//cout<<"generating "<<j<<" random_CD1 tree..."<<endl;
				UB_random_CD = RandomTreeAB(tree, treeC, treeD, taxaC, taxaD, params);
				out_between<<UB_random_CD/UB_true<<"\t";
			}

			// ContraSplit 2: changing 1/4 of taxa
			taxaC = taxaA;
			taxaD = taxaB;
			for(h=0; h<int(n/2.); h++){
				taxaC[h]=taxaB[h];
				taxaD[h]=taxaA[h];
			}
			treeC = extractSubtreeUB(taxaC,tree,params);
			treeD = extractSubtreeUB(taxaD,tree,params);
			for(j=0; j<R; j++){
				//cout<<"generating "<<j<<" random_CD2 tree..."<<endl;
				UB_random_CD = RandomTreeAB(tree, treeC, treeD, taxaC, taxaD, params);
				out_between<<UB_random_CD/UB_true<<"\t";
			}
			out_between<<endl;
*/

/*

	// Printing Tree and its subtrees. This was just for check.
			cout<<"Tree T(A|B)"<<endl;
			tree->printTree(cout,2);
			cout<<endl<<"Tree T(A)"<<endl;
			treeA->printTree(cout,2);
			cout<<endl<<"Tree T(B)"<<endl;
			treeB->printTree(cout,2);
			cout<<endl;

			cout<<"Tree T(A|B)"<<endl;
			printTreeUB(tree);
			cout<<endl<<endl<<"Tree T(A)"<<endl;
			printTreeUB(treeA);
			cout<<endl<<"Tree T(B)"<<endl;
			printTreeUB(treeB);

	// Printing out the results ----------------------------------------------------------
			// Split A|B ------------------------------------------
			out_split<<min(taxaA.size(),taxaB.size())<<"|"<<((double) max(taxaA.size(),taxaB.size()))<<"\t";
			if(min(taxaA.size(),taxaB.size()) == taxaA.size()){
				for(int f = 0; f < taxaA.size()-1; f++)
					out_split<<taxaA[f]<<",";
				out_split<<taxaA[taxaA.size()-1]<<"\t|\t";
				for(int f = 0; f < taxaB.size()-1; f++)
					out_split<<taxaB[f]<<",";
				out_split<<taxaB[taxaB.size()-1];
			} else {
				for(int f = 0; f < taxaB.size()-1; f++)
					out_split<<taxaB[f]<<",";
				out_split<<taxaB[taxaB.size()-1]<<"\t|\t";
				for(int f = 0; f < taxaA.size()-1; f++)
					out_split<<taxaA[f]<<",";
				out_split<<taxaA[taxaA.size()-1];
			}
			out_split<<endl;
			// ----------------------------------------------------

			//out<<min(taxaA.size(),taxaB.size())<<"|"<<((double) max(taxaA.size(),taxaB.size()))<<"\t"<<br_len<<"\t"<<tree->curScore<<"\t";
			out<<params->aln_file<<"\t";
			if(tree->mlInitial == 0)
				out<<"0"<<"\t";
			else
				out<<"1"<<"\t";
			out<<tree->leafNum<<"\t"<<tree->aln->getNSite()<<"\t"<<min(taxaA.size(),taxaB.size())<<"\t"<<br_len<<"\t"<<tree->curScore<<"\t";


			if(min(taxaA.size(),taxaB.size()) == taxaA.size()){
				out<<treeA->curScore<<"\t"<<treeB->curScore<<"\t";
			}
			else{
				out<<treeB->curScore<<"\t"<<treeA->curScore<<"\t";
			}

			out<<tree->aln->size()*(log(1+3*exp(-br_len)) - log(1-exp(-br_len)))<<"\t"<<diff_1<<"\t"<<diff_2;

			if(diff_1>0){
				out<<"\t"<<"FALSE\t0";
				BadSplits1++;
			}else{
				out<<"\t"<<"TRUE\t1";
			}

			if(diff_2>0){
				out<<"\t"<<"FALSE\t0";
				BadSplits2++;
			}else{
				out<<"\t"<<"TRUE\t1";
			}
			out<<endl;

	// END: Printing out the results -----------------------------------------------------
*/
		}} // END: if taxaA.size() and taxaB.size() >3

	}
	// END: the loop over all A|B present on tree T

	out_within.close();
	out_between.close();
	out.close();
	out_split.close();
}

PhyloTree* extractSubtreeUB(IntVector &ids, MTree* tree, Params *params, int sw) {
	string taxa_set;
	int i;
	for(i = 0; i < tree->leafNum; i++)
		taxa_set.push_back(0);
	for (i = 0; i < ids.size(); i++)
		taxa_set[ids[i]]=1;

	PhyloTree *treeCopy = new PhyloTree(); // this will be a new subtree
	Alignment *alignment = new Alignment();
	alignment->extractSubAlignment(((PhyloTree*)tree)->aln,ids,0);

	treeCopy->copyTree(tree, taxa_set);
	treeCopy->setAlignment(alignment);
	if(sw == 1){
		treeCopy->setModel(((PhyloTree*)tree)->getModel());
		treeCopy->setRate(((PhyloTree*)tree)->getRate());
		treeCopy->setModelFactory(((PhyloTree*)tree)->getModelFactory());
		treeCopy->initializeAllPartialLh();
		treeCopy->setCurScore(treeCopy->computeLikelihood());
	}

	return treeCopy;
}

void printTreeUB(MTree *tree){
	int i=0, j=0;

	NodeVector nodeLeaves;
	tree->getTaxa(nodeLeaves);
	cout<<"Taxa nodes:"<<endl;
	for(i=0; i<nodeLeaves.size(); i++){
		cout<<nodeLeaves[i]->name<<":"<<nodeLeaves[i]->id<<"(";
		for(j=0; j<nodeLeaves[i]->neighbors.size(); j++)
			cout<<nodeLeaves[i]->neighbors[j]->node->name<<":"<<nodeLeaves[i]->neighbors[j]->node->id<<"["<<nodeLeaves[i]->neighbors[j]->length<<"]"<<",";
		cout<<")"<<endl;
	}

	NodeVector nodeInternal;
	tree->getInternalNodes(nodeInternal);
	cout<<"Internal nodes:"<<endl;
	if(nodeInternal.size() == 0)
		cout<<"no internal nodes"<<endl;
	else{
		for(i=0; i<nodeInternal.size(); i++){
			cout<<nodeInternal[i]->name<<":"<<nodeInternal[i]->id<<"(";
			for(j=0; j<nodeInternal[i]->neighbors.size(); j++)
				cout<<nodeInternal[i]->neighbors[j]->node->name<<":"<<nodeInternal[i]->neighbors[j]->node->id<<"["<<nodeInternal[i]->neighbors[j]->length<<"]"<<",";
			cout<<")"<<endl;
		}
	}
}

MTree* generateRandomYH_UB(Params &params, PhyloTree *tree){
	MExtTree* treeR = new MExtTree();
	bool binary = TRUE;

	int size = tree->leafNum;
	if (size < 3)
		outError(ERR_FEW_TAXA);

	treeR->root = treeR->newNode();
	int i;
	NodeVector myleaves;
	NodeVector innodes;
	Node *node;
	double len;

	innodes.push_back(treeR->root);
	// create initial tree with 3 leaves
	for (i = 0; i < 3; i++) {
		node = treeR->newNode();
		len = randomLen(params);
		treeR->root->addNeighbor(node, len);
		node->addNeighbor(treeR->root, len);
		myleaves.push_back(node);
	}

	// additionally add a leaf
	for (i = 3; i < size; i++)
	{
		int index;
		if (binary) {
			index = random_int(i);
		} else {
 			index = random_int(i + innodes.size());
		}
		if (index < i) {
			node = myleaves[index];
			innodes.push_back(node);
			// add the first leaf
			Node *newleaf = treeR->newNode();
			len = randomLen(params);
			node->addNeighbor(newleaf, len);
			newleaf->addNeighbor(node, len);
			myleaves[index] = newleaf;

			// add the second leaf
			newleaf = treeR->newNode();
			len = randomLen(params);
			node->addNeighbor(newleaf, len);
			newleaf->addNeighbor(node, len);
			myleaves.push_back(newleaf);
		}
		else {
			node = innodes[index-i];
			// add only 1 new leaf
			Node *newleaf = treeR->newNode();
			len = randomLen(params);
			node->addNeighbor(newleaf, len);
			newleaf->addNeighbor(node, len);
			myleaves.push_back(newleaf);
		}
	}

	treeR->root = myleaves[0];
	// indexing the leaves
	treeR->setLeavesName(myleaves);
	treeR->leafNum = myleaves.size();
	treeR->nodeNum = treeR->leafNum;
	treeR->initializeTree();

	NodeVector taxa;
	treeR->getTaxa(taxa);
	ASSERT(taxa.size() == size);
	for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++)
		(*it)->name = tree->aln->getSeqName((*it)->id);

	return (MTree*)treeR;
}

double RandomTreeAB(PhyloTree* treeORGN, PhyloTree* treeAorgn, PhyloTree* treeBorgn, IntVector &taxaA, IntVector &taxaB, Params* params, double brLen){
	PhyloTree *tree  = new PhyloTree();
	MTree *treeA = new MTree();
	MTree *treeB = new MTree();

	treeA = generateRandomYH_UB(*params,treeAorgn);
	treeB = generateRandomYH_UB(*params,treeBorgn);

/*
	// PrintTree ---------------
	cout<<"TreeA.root:"<<treeA->root->name<<treeA->root->id<<endl;
	cout<<"TreeB.root:"<<treeB->root->name<<treeB->root->id<<endl;
	cout<<"TreeA:"<<endl;
	treeA->printTree(cout);
	cout<<endl<<"TreeB:"<<endl;
	treeB->printTree(cout);
	cout<<endl;
	// -------------------------
*/

	extendingTree(treeA,params);
	extendingTree(treeB,params);

/*
	// PrintTree ---------------
	cout<<"TreeA.root:"<<treeA->root->name<<treeA->root->id<<endl;
	cout<<"TreeB.root:"<<treeB->root->name<<treeB->root->id<<endl;

	cout<<"extended TreeA:"<<endl;
	treeA->printTree(cout);
	cout<<endl<<"extended TreeB:"<<endl;
	treeB->printTree(cout);
	cout<<endl;
	// -------------------------
*/

	treeA->root->name = "NewNodeA";
	treeB->root->name = "NewNodeB";
	treeA->root->addNeighbor(treeB->root,0.0,tree->branchNum);
	treeB->root->addNeighbor(treeA->root,0.0,tree->branchNum);

	tree->copyTree(treeA);
/*	cout<<"Leaves number = "<<tree->leafNum<<endl;
	cout<<"Nodes  number = "<<tree->nodeNum<<endl;
	cout<<"Branch number:"<<tree->branchNum<<endl;
	*/
	//tree->printTree(cout);
	//cout<<endl;

	NodeVector brID;
	//brID= getBranchABid(brLen, tree);
	brID.push_back(tree->findNodeName(treeA->root->name));
	brID.push_back(tree->findNodeName(treeB->root->name));

	if(brLen == 0){
		brLen = randomLen(*params);
	}
	//tree->findNodeName(treeA->root->name)->findNeighbor(treeB->root)->length = brLen;
	//tree->findNodeName(treeB->root->name)->findNeighbor(treeA->root)->length = brLen;


	tree->findNodeID(brID[0]->id)->findNeighbor(brID[1])->length = brLen;
	tree->findNodeID(brID[1]->id)->findNeighbor(brID[0])->length = brLen;


	tree->setAlignment(treeORGN->aln);
	tree->setModel(((PhyloTree*)treeORGN)->getModel());
	tree->setRate(((PhyloTree*)treeORGN)->getRate());
	tree->setModelFactory(((PhyloTree*)treeORGN)->getModelFactory());
	tree->initializeAllPartialLh();

	tree->setCurScore(tree->computeLikelihood());
	//cout<<"LogLh score before optimization: "<<tree->curScore<<endl;
	tree->params = params;
	//tree->curScore = tree->optimizeAllBranches(50);
	//cout<<"LogLh score after  optimization: "<<tree->curScore<<endl;

	//double len = tree->findNodeName(treeA->root->name)->findNeighbor(treeB->root)->length;
	double len = tree->findNodeID(brID[0]->id)->findNeighbor(brID[1])->length;
	//cout<<"The length of corresponding branch after optimization: "<<len<<endl;
	//cout<<"before it was equal to "<<brLen<<endl;

	string out_file = "results.branches.ub";
	ofstream out;
	out.exceptions(ios::failbit | ios::badbit);
	out.open((char*)out_file.c_str(),std::ofstream::out | std::ofstream::app);

	//len = 1;

	double coef = tree->aln->getNSite()*(log(1+3*exp(-len)) - log(1-exp(-len)));
	double U = coef + UpperBoundAB(taxaA, taxaB, tree, params);

	// leafNum		alnLen		brLen (before opt)		brLen (after opt)		coef 		UB
	//out<<treeORGN->leafNum<<"\t"<<treeORGN->aln->getNSite()<<"\t"<<brLen<<"\t"<<len<<"\t"<<coef<<"\t"<<U<<endl;

	out.close();
	return U;
}

double UpperBoundAB(IntVector &taxaA, IntVector &taxaB, PhyloTree* tree, Params *params){
	double U = 0.0;

	PhyloTree *treeA, *treeB;
	treeA = extractSubtreeUB(taxaA,tree,params,1);
	treeB = extractSubtreeUB(taxaB,tree,params,1);

	U = treeA->getCurScore() + treeB->getCurScore();

	return U;
}

NodeVector getBranchABid(double brLen, PhyloTree* tree){
	NodeVector branch1, branch2;
	NodeVector branch;
	tree->getBranches(branch1, branch2);
	for(int i = 0; i != branch1.size(); i++){
		if(branch1[i]->findNeighbor(branch2[i])->length == 0.0){
			branch.push_back(branch1[i]);
			branch.push_back(branch2[i]);
			return branch;
		}
	}
	outError("UpperBounds: did not find matching branch:(");
	return branch;
}

void extendingTree(MTree *tree, Params* params){

	// Choose random internal node
	int maxR = tree->nodeNum-1;
	int randomNodeID = rand() % maxR;
	//cout<<"randomNodeID = "<<randomNodeID<<endl;
	if(randomNodeID<tree->leafNum){
		if(randomNodeID+tree->leafNum > tree->nodeNum-1){
			randomNodeID += tree->nodeNum - tree->leafNum;
			//cout<<"adding nodeNum-leafNum"<<endl;
		}
		else{
			randomNodeID+=tree->leafNum;
			//cout<<"adding leafNum"<<endl;
		}
	}

	//cout<<"leafNum = "<<tree->leafNum-1<<" < random = "<<randomNodeID<<" < nodeNum = "<<tree->nodeNum-1<<endl;

	ASSERT(randomNodeID < tree->nodeNum && randomNodeID > tree->leafNum-1);

	// Choose random neighbor
	int randomNeiID = rand() % 2;
	//cout<<"randomNeiID = "<<randomNeiID<<endl;


	Node *randomNode;
	randomNode = tree->findNodeID(randomNodeID);
	Node *randomNeiNode = tree->findNodeID(randomNodeID)->neighbors[randomNeiID]->node;

	// Create new node ----------------------------
	string str;
	str = "NewNode";
	const char *ch = str.c_str();
	Node* newNode1 = tree->newNode(tree->nodeNum,ch);
	tree->nodeNum++;

	// Add new node as a neighbor to randomNei->node
	double len = tree->findNodeID(randomNodeID)->neighbors[randomNeiID]->length;
	int    id  = tree->findNodeID(randomNodeID)->neighbors[randomNeiID]->id;

	randomNeiNode->findNeighbor(randomNode)->node = newNode1;
	newNode1->addNeighbor(randomNeiNode,len,id);


	//Change randomNei with this new node for randomNode. Create new branch.
	randomNode->neighbors[randomNeiID]->node = newNode1;
	randomNode->neighbors[randomNeiID]->id = tree->branchNum;
	tree->branchNum++;
	randomNode->neighbors[randomNeiID]->length = randomLen(*params);

	newNode1->addNeighbor(randomNode,randomNode->neighbors[randomNeiID]->length,randomNode->neighbors[randomNeiID]->id);

	tree->root = newNode1;

}

NNIMove getBestNNIForBranUB(PhyloNode *node1, PhyloNode *node2, PhyloTree *tree){

	NNIMove nniMoves[2];

    // Initialize node1 and node2 in nniMoves
	nniMoves[0].node1 = nniMoves[1].node1 = node1;
	nniMoves[0].node2 = nniMoves[1].node2 = node2;

	// Initialize two NNIs
	int cnt;
	double t[4];
    FOR_NEIGHBOR_IT(node1, node2, node1_it) {
			cnt = 0;
			t[cnt]=(*node1_it)->length;
			FOR_NEIGHBOR_IT(node2, node1, node2_it) {
				//   Initialize the 2 NNI moves
				nniMoves[cnt].node1Nei_it = node1_it; // for both cnt = 0,1 this is the same neighbor of node1,
													  // which will be swapped with nei1 and nei2 of node2
				nniMoves[cnt].node2Nei_it = node2_it;
				t[cnt+2] = (*node2_it)->length;
				cnt++;
			}
			break;
    }

    NeighborVec::iterator node1Nei2_it;

    FOR_NEIGHBOR_IT(node1, node2, node1_it){
    	if ((*node1_it)->node != (*nniMoves[0].node1Nei_it)->node){
    		t[cnt]=(*node1_it)->length;
    		node1Nei2_it = node1_it;
    		break;
    	}
    }

    /*
     * Correspondence:
     *
     * Nodes, incident to node1 with corresponding branches:
     * nniMoves[0].node1Nei_it	| t[0]
     * node1Nei2_it				| t[1]
     *
     * Nodes, incident to node2 with corresponding branches:
     * nniMoves[0].node2Nei_it	| t[2]
     * nniMoves[1].node2Nei_it	| t[3]
     *
     * NNIs:
     * nniMoves[0] -> swapping (nniMoves[0].node1Nei_it	| t[0]) with (nniMoves[0].node2Nei_it	| t[2])
     * corresponding coef: q1
     *
     * nniMoves[1] -> swapping (nniMoves[1].node1Nei_it	| t[0]) with (nniMoves[1].node2Nei_it	| t[3])
     * corresponding coef: q2
     */

    double L[4]; // likelihoods of 4 subtrees
    double score[4];
    L[0] = L[1] = L[2] = L[3] = 0.0;
    score[0] = score[1] = score[2] = score[3] = 0.0;

    double UB = 0.0; // in log terms
    int nsite = tree->aln->getNSite();
    UB = nsite*logC(node1->findNeighbor(node2)->length,tree); // coefficient c

    //int ncat = tree->site_rate->getNDiscreteRate();
    int nptn = tree->aln->getNPattern();
    int nstates = tree->aln->num_states;
    int i,x;
    //int cat;
    IntVector ptnFreq;
    tree->aln->getPatternFreq(ptnFreq);

    int clear_pl_lh[4]; // if equals to 1, partial likelihoods were computed, don't clear.
    clear_pl_lh[0] = clear_pl_lh[1] = clear_pl_lh[2] = clear_pl_lh[3] = 1;

    double* T1_partial_lh;
    if(((PhyloNeighbor*) (*nniMoves[0].node1Nei_it))->get_partial_lh_computed() == 0){
    	tree->computePartialLikelihood((PhyloNeighbor*) (*nniMoves[0].node1Nei_it), node1);
    	clear_pl_lh[0] = 0;
    }
    T1_partial_lh = ((PhyloNeighbor*) (*nniMoves[0].node1Nei_it))->get_partial_lh();

    double* T2_partial_lh;
    if(((PhyloNeighbor*) (*node1Nei2_it))->get_partial_lh_computed() == 0){
    	tree->computePartialLikelihood(((PhyloNeighbor*) (*node1Nei2_it)), node1);
    	clear_pl_lh[1] = 0;
    }
    T2_partial_lh = ((PhyloNeighbor*) (*node1Nei2_it))->get_partial_lh();

    double* T3_partial_lh;
    if(((PhyloNeighbor*) (*nniMoves[0].node2Nei_it))->get_partial_lh_computed() == 0){
    	tree->computePartialLikelihood(((PhyloNeighbor*) (*nniMoves[0].node2Nei_it)), node1);
    	clear_pl_lh[2] = 0;
    }
    T3_partial_lh = ((PhyloNeighbor*) (*nniMoves[0].node2Nei_it))->get_partial_lh();

    double* T4_partial_lh;
    if(((PhyloNeighbor*) (*nniMoves[1].node2Nei_it))->get_partial_lh_computed() == 0){
    	tree->computePartialLikelihood(((PhyloNeighbor*) (*nniMoves[1].node2Nei_it)), node1);
    	clear_pl_lh[3] = 0;
    }
    T4_partial_lh = ((PhyloNeighbor*) (*nniMoves[1].node2Nei_it))->get_partial_lh();

    for(i = 0; i<nptn; i++){
    	score[0] = score[1] = score[2] = score[3] = 0.0;
    	// Sum over Gamma categories and over states
    	//for(cat = 0; cat < ncat; cat++){
    		for(x = 0; x < nstates; x++){
    		// First  subtree --------------------------
    			score[0] += tree->getModel()->state_freq[x]*T1_partial_lh[i*nstates+x];
    		// Second subtree --------------------------
    			score[1] += tree->getModel()->state_freq[x]*T2_partial_lh[i*nstates+x];
    	   	// Third  subtree --------------------------
    			score[2] += tree->getModel()->state_freq[x]*T3_partial_lh[i*nstates+x];
    	   	// Fourth subtree --------------------------
    			score[3] += tree->getModel()->state_freq[x]*T4_partial_lh[i*nstates+x];
    		}
   	//}
    	L[0] += log(score[0])*ptnFreq[i];
    	L[1] += log(score[1])*ptnFreq[i];
    	L[2] += log(score[2])*ptnFreq[i];
    	L[3] += log(score[3])*ptnFreq[i];

        ASSERT(isnormal(L[0] + L[1] + L[2] + L[3]));

    }

/*
    if(clear_pl_lh[0] == 0){
    	((PhyloNeighbor*) (*nniMoves[0].node1Nei_it))->clearPartialLh();
    }
    if(clear_pl_lh[1] == 0){
    	((PhyloNeighbor*) (*node1Nei2_it))->clearPartialLh();
    }
    if(clear_pl_lh[2] == 0){
    	((PhyloNeighbor*) (*nniMoves[0].node2Nei_it))->clearPartialLh();
    }
    if(clear_pl_lh[3] == 0){
    	((PhyloNeighbor*) (*nniMoves[1].node2Nei_it))->clearPartialLh();
    }*/
    //cout<<"Clear_pl_lh:"<<clear_pl_lh[0]<<" "<<clear_pl_lh[1]<<" "<<clear_pl_lh[2]<<" "<<clear_pl_lh[3]<<endl;

    //double logNcat = log(((double)ncat));
    L[0] = L[0] + ((PhyloNeighbor*) (*nniMoves[0].node1Nei_it))->get_lh_scale_factor();
    L[1] = L[1] + ((PhyloNeighbor*) (*node1Nei2_it))->get_lh_scale_factor();
    L[2] = L[2] + ((PhyloNeighbor*) (*nniMoves[0].node2Nei_it))->get_lh_scale_factor();
    L[3] = L[3] + ((PhyloNeighbor*) (*nniMoves[1].node2Nei_it))->get_lh_scale_factor();

/*   // Print some info:
    cout<<"The log likelihood  of the parent tree T:"<<tree->computeLikelihood()<<endl;
    cout<<"The log likelihoods of the four subtrees:"<<endl;
    cout<<"Node"<<(*nniMoves[0].node1Nei_it)->node->id<<": L[0] = "<<L[0]<<endl;
    cout<<"Node"<<(*node1Nei2_it)->node->id<<": L[1] = "<<L[1]<<endl;
    cout<<"Node"<<(*nniMoves[0].node2Nei_it)->node->id<<": L[2] = "<<L[2]<<endl;
    cout<<"Node"<<(*nniMoves[1].node2Nei_it)->node->id<<": L[3] = "<<L[3]<<endl;*/

    UB += L[0] + L[1] + L[2] + L[3];

    double q1 = logC(t[0]+t[3],tree) + logC(t[1]+t[2],tree);
    double q2 = logC(t[0]+t[2],tree) + logC(t[1]+t[3],tree);
    //cout<<"Coefficients q1 and q2:"<<endl<<q1<<endl<<q2<<endl;

    double UBq1 = UB + nsite*q1;
    double UBq2 = UB + nsite*q2;

	string out_file_UB = tree->params->out_prefix;
	out_file_UB += ".UB.NNI.upperBounds";
	ofstream out_UB;
	out_UB.exceptions(ios::failbit | ios::badbit);
	out_UB.open((char*)out_file_UB.c_str(),std::ofstream::out | std::ofstream::app);

	out_UB << tree->getCurScore() << "\t" << UBq1 << "\t" << UBq2 << "\t" << tree->getCurScore() - UBq1 << "\t" << tree->getCurScore() - UBq2 << endl;

	out_UB.close();

    if(UBq1 < tree->getCurScore()){
    	tree->skippedNNIub += 1;
  /*  	tree->meanUB += UBq1;
    	if(UBq1 < tree->minUB){
    		tree->minUB = UBq1;
    	} else if(UBq1 > tree->maxUB){
    		tree->maxUB = UBq1;
    	}*/
    	//cout<<"----------------- UBq1 < L !!!"<<endl;
    }
    if(UBq2 < tree->getCurScore()){
    	tree->skippedNNIub += 1;
    	/*tree->meanUB += UBq2;
    	if(UBq2 < tree->minUB){
    		tree->minUB = UBq2;
    	} else if(UBq2 > tree->maxUB){
    		tree->maxUB = UBq2;
    	}*/
    	//cout<<"----------------- UBq2 < L !!!"<<endl;
    }

	// Decide which NNI has a larger UB (we base our decision on q coefficients)
	if(q1 > q2){
		// NNI 1:
		//nniMoves[0].newLen[0] = NULL;
		nniMoves[0].newloglh = UBq1;
		//cout<<"q1 and NNI1 is chosen with UB "<<UBq1<<endl;
		return nniMoves[0];
	} else {
		// NNI 2:
		//nniMoves[1].newLen[0] = NULL;
		nniMoves[1].newloglh = UBq2;
		//cout<<"q2 and NNI2 is chosen with UB "<<UBq2<<endl;
		return nniMoves[1];
	}
}

double logC(double t, PhyloTree* tree){
	//double c = log((1+3*exp(-t)))-log(1-exp(-t));

	int i, m = tree->aln->num_states*tree->aln->num_states, n = tree->aln->num_states;
	double* TransMatrix = new double[m];
	tree->getModelFactory()->computeTransMatrix(t,TransMatrix);
	double maxTransProb = 0.0;
	for(i = 0; i < m; i++)
		if(TransMatrix[i]>maxTransProb)
			maxTransProb = TransMatrix[i];
	//maxTransProb = 0.25*(1+3*exp(-3*t/4));
	//maxTransProb = 1;

	if(tree->minStateFreq == 0.0){
		tree->minStateFreq = tree->getModel()->state_freq[0];
		for(i = 1; i < n; i++){
			if(tree->minStateFreq > tree->getModel()->state_freq[i])
				tree->minStateFreq = tree->getModel()->state_freq[i];
		}

	}
	//cout<<tree->minStateFreq<<endl;
	//tree->minStateFreq = 0.25;
	//assert(isnormal(log(maxTransProb/tree->minStateFreq)));
	return log(maxTransProb/tree->minStateFreq);
}

void sumFraction(PhyloNode *node1, PhyloNode *node2, PhyloTree *tree){
	PhyloNeighbor* nei1 = node1->findNeighbor(node2);
	PhyloNeighbor* nei2 = node2->findNeighbor(node1);

//	int loglh = tree->computeLikelihood();

    double* T1_partial_lh;
    if(nei1->get_partial_lh_computed() == 0){
    	tree->computePartialLikelihood(nei1, node1);
    }
    T1_partial_lh = nei1->get_partial_lh();

    double* T2_partial_lh;
    if(nei2->get_partial_lh_computed() == 0){
    	tree->computePartialLikelihood(nei2, node2);
    }
    T2_partial_lh = nei2->get_partial_lh();

    double score[3];
    score[0] = score[1] = score[2] = 0.0;

    double plh[3];
    plh[0]=plh[1]=plh[2]=0.0;

    int nptn = tree->aln->getNPattern();
    int nstates = tree->aln->num_states;
    int j,i,x,y;

    double plhx[nstates];
    double plhy[nstates];

    double *eigen = tree->getModel()->getEigenvectors();

    for(i = 0; i<nptn; i++){
    	score[0] = score[1] = score[2] = 0.0;

    	// computing partial likelihoods
		for(x = 0; x < nstates; x++){
			plhx[x] = 0.0;
			plhy[x] = 0.0;
			for(j = 0; j<nstates; j++){
				plhx[x]+= T1_partial_lh[i*nstates+j]*eigen[x*nstates+j];
				plhy[x]+= T2_partial_lh[i*nstates+j]*eigen[x*nstates+j];
			}
		}

		for(x = 0; x < nstates; x++){
			for(y = 0; y < nstates; y++){
				if(x == y){
				// Term for a pair of matching nucleotides --------------------------
					score[0] += plhx[x]*plhy[y];
				} else {
				// Term for a pair of non-matching nucleotides ----------------------
					score[1] += plhx[x]*plhy[y];
				}
			// Full sum ---------------------------------------------------------
			score[2] += plhx[x]*plhy[y];
			}
		}

        ASSERT(isnormal(score[0] + score[1] + score[2]));

		cout<<"BranchLEN |"<< nei1->length
			<<"| FRACTION of ai (sum over matching pairs) |"<<score[0]/score[2]
		    <<"| FRACTION of bi (sum over non-matching pairs) |"<<score[1]/score[2]
		    <<"| likelihood |"<<score[2]
		    <<endl;

    }


}
