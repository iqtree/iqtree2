/*
 * upperbounds.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: olga
 */
#include "upperbounds.h"

void UpperBounds(Params *params, Alignment* alignment, IQTree* tree){

	// Output details --------------------------------------------------
	// UpperBounds File
	string out_file = params->out_prefix;
	out_file += ".ub";
	out_file = "results.trueSplits.ub";
	ofstream out;
	out.exceptions(ios::failbit | ios::badbit);
	out.open((char*)out_file.c_str(),std::ofstream::out | std::ofstream::app);
	// Details on Split: A|B
	string out_file_split = params->out_prefix;
	out_file_split += ".ub_split";
	out_file_split = "results.trueSplits.ub.splits";
	ofstream out_split;
	out_split.exceptions(ios::failbit | ios::badbit);
	out_split.open((char*)out_file_split.c_str(),std::ofstream::out | std::ofstream::app);

	// Within Family Info: A|B
		string out_file_within = params->out_prefix;
		out_file_within += ".ub_split";
		out_file_within = "results.within.ub";
		ofstream out_within;
		out_within.exceptions(ios::failbit | ios::badbit);
		out_within.open((char*)out_file_within.c_str(),std::ofstream::out | std::ofstream::app);

	// Auxiliary variables ---------------------------------------------
	int i=0, j=0, h=0;

	// STARTing the Analysis -------------------------------------------

	cout<<"Starting Upper Bounds analysis.."<<endl;

	// Printing info about the TreeLogL changes during the tree search
/*	cout<<"mlInitial  = "<<tree->mlInitial<<endl;
	cout<<"mlFirstOpt = "<<tree->mlFirstOpt<<endl;
	cout<<"mlBestTree = "<<tree->getBestScore()<<endl;
	cout<<"mlUnConstr = "<<alignment->computeUnconstrainedLogL()<<endl;*/

	//double mlQuestionary = tree->mlInitial; //or tree->mlFirstOpt for example

	NodeVector branch1, branch2;
	tree->getBranches(branch1, branch2);
	int BadSplits1 = 0, BadSplits2 = 0, allSplits = 0;

// HERE must be a loop over all A|B present on tree T
	for(i = 0; i != branch1.size(); i++){
		vector<int> taxaA, taxaB;
		vector<string> taxaAname, taxaBname;
		tree->getTaxaID(taxaA,branch1[i],branch2[i]);
		tree->getTaxaID(taxaB,branch2[i],branch1[i]);

		if(taxaA.size() > 3 and taxaB.size() > 3){ // IQTree does not compute lh of tree with less than 4 taxa.
			allSplits++;

	// Getting taxa_names for taxon subsets A and B
			for(h = 0; h < taxaA.size(); h++)
				taxaAname.push_back(tree->findNodeID(taxaA[h])->name);
			for(h = 0; h < taxaB.size(); h++)
				taxaBname.push_back(tree->findNodeID(taxaB[h])->name);

			/* ---------------------------------------------------------------------------------------
			 * Collective Stuff for artificial Root. Possibly not needed anymore
			 * --------------------------------------------------------------------------------------- */
	/*
			PhyloTree *treeR; // tree with artificial root
			string ch;
			ch = "ATGC"; // If one is interested in adding different nucleotides in sequence for artificial root
			string art_root_name = "art_root";
			int treeW = 0; // specifies the subtree with artificial root; 0 for treeA, 1 for treeB.
			if(treeW == 0){
				treeR= extractSubtreeUB(taxaA,tree,params,1,ch[1]);
			} else {
				treeR= extractSubtreeUB(taxaB,tree,params,1,ch[1]);
			}
			Node *art_root_node = treeR->findNodeName(art_root_name);
			art_root_node->neighbors[0]->length = tree->findNodeID(branch1[i]->id,branch1[i],branch2[i])->findNeighbor(branch2[i])->length;
			treeR->findNodeID(art_root_node->neighbors[0]->node->id)->findNeighbor(art_root_node)->length = art_root_node->neighbors[0]->length;
	*/
			/* ---------------------------------------------------------------------------------------
			 * END of artificial root.
			 * --------------------------------------------------------------------------------------- */

	// Dealing with subtrees T_A and T_B
			PhyloTree *treeA, *treeB;
			treeA = extractSubtreeUB(taxaA,tree,params);
			treeB = extractSubtreeUB(taxaB,tree,params);


			//Random trees
			out_within<<min(taxaA.size(),taxaB.size())<<"\t"<<treeA->curScore+treeB->curScore<<"\t";
			for(j=0; j<50; j++){
				cout<<"generating "<<j<<" random tree..."<<endl;
				double llh_randomA = 0.0, llh_randomB = 0.0;
				llh_randomA = generateRandomYH_UB(*params, treeA);
				llh_randomB = generateRandomYH_UB(*params, treeB);
				out_within<<llh_randomA+llh_randomB<<"\t";
				//cout<<"TreeA_true  : llh = "<<treeA->curScore<<endl;
				//cout<<"TreeA_random: llh = "<<llh_randomA<<endl;
			}
			out_within<<endl;

	/*		//Between Families comparison: true split vs contradicting splits
			int n=0;
			n=int(min(taxaA.size(),taxaB.size())/2.);
			cout<<"taxaA.size() = "<<taxaA.size()<<", taxaB.size() = "<<taxaB.size()<<", n = "<<n<<endl;

			vector<int> taxaAnew, taxaBnew;
			PhyloTree *treeAnew, *treeBnew;

			// ContraSplit1: changing 1/2 of taxa
			taxaAnew = taxaA;
			taxaBnew = taxaB;

			for(h=0; h<n; h++){
				taxaAnew[h]=taxaB[h];
				taxaBnew[h]=taxaA[h];
			}

			treeAnew = extractSubtreeUB(taxaAnew,tree,params);
			treeBnew = extractSubtreeUB(taxaBnew,tree,params);
			out_within<<min(taxaA.size(),taxaB.size())<<"\t"<<treeA->curScore+treeB->curScore<<"\t";

			for(j=0; j<10; j++){
				double randomA = 0.0, randomB = 0.0;
				randomA = generateRandomYH_UB(*params, treeAnew);
				randomB = generateRandomYH_UB(*params, treeBnew);

				out_within<<randomA+randomB<<"\t";

				cout<<"trueSplit="<<treeA->curScore+treeB->curScore<<", contraSplit="<<randomA+randomB<<endl;
			}

			// ContraSplit 2: changing 1/4 of taxa
			taxaAnew = taxaA;
			taxaBnew = taxaB;

			for(h=0; h<int(n/2.); h++){
				taxaAnew[h]=taxaB[h];
				taxaBnew[h]=taxaA[h];
			}

			treeAnew = extractSubtreeUB(taxaAnew,tree,params);
			treeBnew = extractSubtreeUB(taxaBnew,tree,params);

			for(j=0; j<10; j++){
				double randomA = 0.0, randomB = 0.0;
				randomA = generateRandomYH_UB(*params, treeAnew);
				randomB = generateRandomYH_UB(*params, treeBnew);

				out_within<<randomA+randomB<<"\t";

				cout<<"trueSplit="<<treeA->curScore+treeB->curScore<<", contraSplit="<<randomA+randomB<<endl;
			}

			// end of contra split for one split*/
			out_within<<endl;

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

	// Getting all interesting info:)


			 *  LogL T_A|B - does not change for A|B
			 *  double llhT = tree->curScore;
			 *   Interesting quantities:
			 *  	diff_1 = LogL(T_A|B) - LogL(T_A) - LogL(T_B)
			 *  	diff_2 = LogL(T_A|B) - LogL(T_A) - LogL(T_B) - N*log(c(A|B))
			 *  	if diff_# < 0 everything is ok, if > 0, the corresponding inequality is violated
			 *
			 *  LogL(T_B'_fixed_root_nucleotide) < LogL(T_B) < LogL(T_B'_not_fixed_root_nucleotide)
			 *


			double diff_1, diff_2;
			double br_len;
			br_len = branch1[i]->findNeighbor(tree->findNodeID(branch2[i]->id))->length;
			diff_2 = tree->curScore - treeA->curScore - treeB->curScore - tree->aln->size()*(log(1+3*exp(-br_len)) - log(1-exp(-br_len)));
			diff_1 = tree->curScore - treeA->curScore - treeB->curScore;

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

		} // END: if taxaA.size() and taxaB.size() >3.

	}
	// END: the loop over all A|B present on tree T

	out.close();
	out_split.close();
	out_within.close();
/*	cout<<"BadSplits1 = "<<BadSplits1<<endl;
	cout<<"BadSplits2 = "<<BadSplits2<<endl;
	cout<<"All tested = "<<allSplits<<endl;*/

	// Some more questions:
	//		- is the UB too far from the TreeLogL?
	//		- dependency on the length of t(A|B)
	//		- evaluate these quantities also for the splits on the initial tree
	//		- can we exclude some splits based on their UB?
}

PhyloTree* extractSubtreeUB(IntVector &ids, MTree* tree, Params *params,int type, char ch) {
	string taxa_set;
	int i;
	for(i = 0; i < tree->leafNum; i++)
		taxa_set.push_back(0);
	for (i = 0; i < ids.size(); i++)
		taxa_set[ids[i]]=1;

	PhyloTree *treeCopy = new PhyloTree(); // this will be a new subtree
	Alignment *alignment = new Alignment();
	alignment->extractSubAlignment(((PhyloTree*)tree)->aln,ids,0);

	if(type == 0){
		//cout<<"Copying tree. Normal procedure."<<endl;
		treeCopy->copyTree(tree, taxa_set);
		treeCopy->setAlignment(alignment);
	} else {
		/*
		 * This part was written for an artificial root. There are might be some problems with it.
		 * In fact we do not need it anymore. Delete.
		 */
		copyTreeUB(tree, treeCopy, taxa_set);
		reindexTaxonIDs(treeCopy);
		string outfile;
		outfile = "aux.alignment";
		alignment->printFasta(outfile.c_str(), false,
				params->aln_site_list, params->aln_nogaps, params->aln_no_const_sites,
				params->ref_seq_name);
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(outfile.c_str(),std::ofstream::out | std::ofstream::app);
		out<<">art_root"<<endl;
		for(i=0; i<alignment->getNSite(); i++)
			out<<"?";
			//out<<ch;
		out<<endl;
		out.close();
		delete(alignment);
		char *cstr = new char[outfile.length() + 1];
		strcpy(cstr, outfile.c_str());
		alignment = new Alignment(cstr, params->sequence_type, params->intype);
		treeCopy->setAlignment(alignment);
		treeCopy->aln->printFasta(outfile.c_str(), false,
						params->aln_site_list, params->aln_nogaps, params->aln_no_const_sites,
						params->ref_seq_name);
		remove (outfile.c_str());
		delete [] cstr;
	}

	treeCopy->setModel(((PhyloTree*)tree)->getModel());
	treeCopy->setRate(((PhyloTree*)tree)->getRate());
	treeCopy->setModelFactory(((PhyloTree*)tree)->getModelFactory());
	treeCopy->initializeAllPartialLh();

	treeCopy->curScore = treeCopy->computeLikelihood();
	//cout<<"SubtreeLikelihood : lh = "<<treeCopy->computeLikelihood()<<endl;
	//printTreeUB(treeCopy);

	return treeCopy;
}

void copyTreeUB(MTree *tree, MTree *treeCopy, string &taxa_set) {
    if (tree->leafNum != taxa_set.length()) outError("#leaves and taxa_set do not match!");
    treeCopy->leafNum = treeCopy->nodeNum = treeCopy->branchNum = 0;
    for (string::iterator it = taxa_set.begin(); it != taxa_set.end(); it++)
    	treeCopy->nodeNum += (*it);
    double new_len;
    if (treeCopy->root) treeCopy->freeNode();
    treeCopy->root = NULL;
    treeCopy->root = copyTreeUBnode(tree,treeCopy,taxa_set, new_len, NULL, NULL);
}

Node* copyTreeUBnode(MTree *tree, MTree *treeCopy, string &taxa_set, double &len, Node *node, Node *dad) {
    if (!node) {
        if (taxa_set[tree->root->id]) {
            node = tree->root;
            //cout<<"The tree root is in subtree."<<endl;
        } else {
            for (int i = 0; i < tree->leafNum; i++)
                if (taxa_set[i]) {
                    node = tree->findNodeID(i);
                    break;
                }
        }
        //cout<<"First node name:id "<<node->name<<":"<<node->id<<endl;
    }
    Node *new_node = NULL;
    NodeVector new_nodes;
    DoubleVector new_lens;
    if (node->isLeaf()) {
        len = 0.0;
        if (taxa_set[node->id]) {
            new_node = treeCopy->newNode(treeCopy->leafNum++, node->name.c_str());
            //nodesB_OLD.push_back(node);
            //nodesB_NEW.push_back(new_node);
        }
        if (dad) return new_node;
    }
    if (new_node) {
        new_nodes.push_back(new_node);
        new_lens.push_back(len);
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        double new_len;
        new_node = copyTreeUBnode(tree, treeCopy, taxa_set, new_len, (*it)->node, node);
        if (new_node) {
            new_nodes.push_back(new_node);
            //new_lens.push_back((*it)->length + new_len); //change here!!!
            new_lens.push_back((*it)->length);
        }
    }
    if (new_nodes.empty()) return NULL;

    // This was changed to add an artificial root to a subtree

    /*
     * !"%$&$%ZFVFGRT Tut kakoe-to portachivo vishlo:(((( ili s reindexing. Koroche gde-to zdes ein Feller!!!:(
     */

    if (new_nodes.size() == 1) {
    	new_node = treeCopy->newNode(treeCopy->leafNum++, dad->name.c_str());

        new_node->addNeighbor(new_nodes[0], new_lens[0], treeCopy->branchNum);
        new_nodes[0]->addNeighbor(new_node, new_lens[0], treeCopy->branchNum);
        treeCopy->branchNum++;

        // Artificial Root
        string root_name = "art_root";
        Node *art_root = new Node(treeCopy->leafNum++, root_name.c_str());
        new_node->addNeighbor(art_root,0.0,treeCopy->branchNum);
        art_root->addNeighbor(new_node,0.0,treeCopy->branchNum);
        treeCopy->branchNum++;

        //len = new_lens[0];
        return new_node;
    }

    // Here change (for cherry case): add an artificial ROOT
    if (!dad && new_nodes.size() == 2) {
        double sum_len = new_lens[0] + new_lens[1]; //change here!!!
        new_nodes[0]->addNeighbor(new_nodes[1], sum_len, treeCopy->branchNum);
        new_nodes[1]->addNeighbor(new_nodes[0], sum_len, treeCopy->branchNum);
        treeCopy->branchNum++;
        return new_nodes[0];
    }

    Node* int_node = treeCopy->newNode(treeCopy->nodeNum++, node->name.c_str());
    len = 0.0;
    for (int i = 0; i < new_nodes.size(); i++) {
        int_node->addNeighbor(new_nodes[i], new_lens[i], treeCopy->branchNum);
        new_nodes[i]->addNeighbor(int_node, new_lens[i], treeCopy->branchNum);
        treeCopy->branchNum++;
    }
    return int_node;
}

void reindexTaxonIDs(MTree *tree){
	int i=0, j=0;
	NodeVector nodesTaxa, nodesInternal;
	tree->getTaxa(nodesTaxa);
	for(i=0; i<nodesTaxa.size(); i++)
		nodesTaxa[i]->id = i;
	tree->getInternalNodes(nodesInternal);
	for(j=0; j<nodesInternal.size(); j++)
		nodesInternal[j]->id = i++;
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

double generateRandomYH_UB(Params &params, PhyloTree *tree){
	MExtTree* treeR = new MExtTree();
	bool binary = FALSE;

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
	assert(taxa.size() == size);
	for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++)
		(*it)->name = tree->aln->getSeqName((*it)->id);

	PhyloTree* treeRphylo = new PhyloTree();
	treeRphylo->copyTree((MTree*)treeR);
	treeRphylo->setAlignment(tree->aln);
	treeRphylo->setModel(((PhyloTree*)tree)->getModel());
	treeRphylo->setRate(((PhyloTree*)tree)->getRate());
	treeRphylo->setModelFactory(((PhyloTree*)tree)->getModelFactory());
	treeRphylo->initializeAllPartialLh();

	treeRphylo->curScore = treeRphylo->computeLikelihood();

/*	cout<<endl;
	tree->printTree(cout);
	cout<<endl;
	treeRphylo->printTree(cout);
	cout<<endl;*/

	double score = treeRphylo->curScore;
	//delete treeRphylo;
	delete treeR;
	return score;
}
