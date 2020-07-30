/*
 * ecopd.cpp
 *
 *  Created on: Oct 30, 2013
 *      Author: Olga
 */

#include <sstream>
#include "ecopd.h"
#include "split.h"
#include "pdnetwork.h"
#include "tree/mtreeset.h"
#include "ecopdmtreeset.h"
#include "graph.h"



ECOpd::ECOpd(const char *userTreeFile, bool &is_rooted) : MTree(userTreeFile,is_rooted) {}

ECOpd::ECOpd() :MTree(){}

ECOpd::~ECOpd(){}

void ECOpd::initializeEcoPD(Params &params){
}

void ECOpd::readInitialTaxa(const char *infile){
	ifstream in;
	//cout<<"Reading taxa to be included in the final optimal subset from file: "<<infile<< endl;
	try {
		in.exceptions(ios::failbit | ios::badbit);
		in.open(infile);
		in.exceptions(ios::badbit);
		readInitialTaxa(in);
		in.close();
	} catch (const char* str) {
		outError(str);
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, infile);
	}
}
void ECOpd::readInitialTaxa(istream &in){
	string name;
	while(!in.eof()){
		in>>name;
		initialTaxa.push_back(name);
	}
	if(initialTaxa.size() != 0){
		initialTaxa.erase(initialTaxa.end());
	}
}

bool ECOpd::OUT_tree(int i){
	bool check=false;
	for(int j=0;j<OUTtreeTaxa.size();j++)
		if(OUTtreeTaxa[j]==i){
			check=true;
			//cout<<"Taxon "<<i<<" is not present in the tree."<<endl;
		}
	return(check);
}

void ECOpd::readDAG(const char* infile) {
	ifstream in;
	//cout<<endl<<"-----------------------------------------------------"<<endl;
	if(weighted)
		cout<<"Reading Diet Composition matrix from file: "<<infile<<endl;
	else
		cout<<"Reading Food Web matrix from file: "<<infile<<endl;
	try {
		in.exceptions(ios::failbit | ios::badbit);
		in.open(infile);
		in.exceptions(ios::badbit);
		readDAG(in);
		in.close();
	} catch (const char* str) {
		outError(str);
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, infile);
	}
}

void ECOpd::readDAG(istream &in) {
	int i=0,j=0;

/* ---------------------------------------------------------------------------------------------------------
 * Reading Diet Composition matrix from the input file
 * ---------------------------------------------------------------------------------------------------------*/
	if(!(in>>SpeciesNUM)) throw "The first line must contain the number of species in this Food Web!!";
	string str_rest, speciesName;
	getline(in, str_rest);

	if(rooted)
		SpeciesNUM++;

	vector<double*> MM;
	for(i=0;i<SpeciesNUM;i++){
		MM.push_back(new double [SpeciesNUM]);
	}

	nvar = (TaxaNUM > SpeciesNUM) ? TaxaNUM : SpeciesNUM;
	for(i=0;i<nvar;i++){
		DAG.push_back(new double [nvar]);
		for(j=0; j<nvar; j++){
			DAG[i][j] = 0.0;
		}
	}
	i = 0;
	j = 0;
	if(rooted){
		while(i != SpeciesNUM-1){
			if(!(in >> speciesName)) throw "Each line should start with a species name!";
			dagNames.push_back(speciesName);
			j = 0;
			while(j != SpeciesNUM-1){
				if(!(in >> MM[i][j])) throw "Could not read matrix entry! For each species make sure there are as many entries as the number of species specified in the file. Only square matrices are accepted.";
				if(MM[i][j] < 0) throw "The Food Web matrix should not contain negative values.Use either 0, 1 or a positive number to indicate the portion of diet.";
				j++;
			}
			MM[i][SpeciesNUM-1] = 0;
			i++;
		}
		for(j=0; j<SpeciesNUM; j++)
			MM[SpeciesNUM-1][j] = 0;
		dagNames.push_back("_root");
	} else {
		while(i != SpeciesNUM){
				if(!(in >> speciesName)) throw "Each line should start with a species name!";
				dagNames.push_back(speciesName);
				j = 0;
				while(j != SpeciesNUM){
					if(!(in >> MM[i][j])) throw "Could not read matrix entry! For each species make sure there are as many entries as the number of species specified in the file. Only square matrices are accepted.";
					if(MM[i][j] < 0) throw "The Food Web matrix should not contain negative values.Use either 0, 1 or a positive number to indicate the portion of diet.";
					j++;
				}
				i++;
			}
	}

/* ---------------------------------------------------------------------------------------------------------
 * Input data
 * ---------------------------------------------------------------------------------------------------------*/
	if(verbose_mode == VB_MAX){
		cout<<endl<<"Food web is defined by the following matrix"<<endl;
		for(i=0;i<SpeciesNUM;i++) {
			cout<<dagNames[i]<<"\t";
			for(j=0;j<SpeciesNUM;j++)
				cout<<MM[i][j]<<"\t";
			cout<<endl;
		}
		// Species in the food web and their ids
		for(i=0; i<SpeciesNUM;i++)
			cout<<"["<<i<<"] "<<dagNames[i]<<endl;
	}
/* ---------------------------------------------------------------------------------------------------------
 * Processing the input data
 * ---------------------------------------------------------------------------------------------------------*/
//Ignoring cannibalism -------------------------------------------------------------------------
	int cannibals=0;
	for(i=0;i<SpeciesNUM;i++)
		if(MM[i][i]!=0){
			cannibals++;
			if(weighted){
				if(cannibals == 1){
					cout<<"------------------------------------------"<<endl;
					cout<<"    Cannibal species         link weight  "<<endl;
					cout<<"------------------------------------------"<<endl;
				}
				cout.width(30);
				cout<<left<<dagNames[i];
				cout<<" | "<<MM[i][i]<<endl;
			}else{
				if(cannibals == 1){
					cout<<"-----------------------------"<<endl;
					cout<<"       Cannibal species      "<<endl;
					cout<<"-----------------------------"<<endl;
				}
				cout<<dagNames[i]<<endl;
			}
			MM[i][i]=0;
		}
	if(cannibals!=0){
		cout<<endl<<"Deleted "<<cannibals;
		if(cannibals == 1)
			cout<<" cannibalistic link."<<endl;
		else
			cout<<" cannibalistic links."<<endl;
	}

//Check whether the graph is acyclic or not
	Graph g(SpeciesNUM);
	for(i=0; i<SpeciesNUM; i++)
		for(j=0; j<SpeciesNUM; j++)
			if(MM[i][j]>0)
				g.addEdge(i,j);
	if(g.isCyclic()){
		if(cannibals != 0)
			cout<<endl<<"ERROR: Even after deleting cannibalistic links, there are still some cycles present."<<endl;
		else
			cout<<endl<<"ERROR: ";
		cout<<"Cyclic food webs are not supported. Delete the links which cause cycles and run again."<<endl;
		cout<<"SOLUTION:"<<endl;
		cout<<"Detect species in the cycle and choose one link to be deleted in order to break the cycle."<<endl;
		cout<<"One possibility is to delete the link with least weight. This can be done by setting the corresponding value in the matrix to 0."<<endl;
		exit(0);
	}

// The number of links -------------------------------------------------------------------------
	linksNUM = 0;
	for(i = 0; i<SpeciesNUM; i++)
		for(j = 0; j<SpeciesNUM; j++)
			if(MM[i][j]>0)
				linksNUM++;

//Rescaling the diet if necessary --------------------------------------------------------------
	if(weighted){
		int dietReScaled = 0;
		vector<double> colsum;
		//cout<<"Food web is weighted."<<endl;
		for(j=0;j<SpeciesNUM;j++){

			colsum.push_back(0);
			for(i=0;i<SpeciesNUM;i++)
				colsum[j]=colsum[j]+MM[i][j];
			if(colsum[j]!=1 && colsum[j]!=0){
				dietReScaled++;
				//cout<<"    WARNING: rescaled diet composition of species "<<j<<". Column sum = "<<colsum[j]<<endl;
				for(i=0;i<SpeciesNUM;i++)
					MM[i][j]=MM[i][j]/colsum[j];
			}
			colsum[j]=0;
			//for(i=0;i<SpeciesNUM;i++)
			//	colsum[j]=colsum[j]+MM[i][j];
			//cout<<j<<"  Column sum = "<<colsum[j]<<endl;
		}
		cout<<"Rescaled diet composition of "<<dietReScaled<<" species."<<endl;
	}else{
		for(i=0; i<SpeciesNUM; i++)
			for(j=0; j<SpeciesNUM; j++)
				if( MM[i][j] > 0)
					MM[i][j] = 1;
		//cout<<"Since the -eco option was chosen, the entries of Food Web matrix will be converted to 0/1 [not prey / prey]. You can use -ecoW option to account for the Diet Composition."<<endl;
	}

// Technical: in case of rooted trees, we check which species are basal ones, i.e. for which check = 0, and set them to "feed on" root M[i,j] = 1
	if(rooted){
		vector<double> check;
		for(j=0;j<SpeciesNUM-1;j++){
			check.push_back(0);
			for(i=0;i<SpeciesNUM-1;i++)
				check[j]=check[j]+MM[i][j];
			if(check[j]==0)
				MM[SpeciesNUM-1][j]=1;
		}
	}

//Detecting which species are not present in either FoodWeb or Tree/SplitNetwork-----------------
	detectMissingSpecies();

//Check whether all the species from initialTaxa set are actually present on Tree/SplitSys or in Food Web
	checkInitialTaxa();

// Synchronization of species in Tree/SplitSys and species in FoodWeb ---------------------------
	synchronizeSpecies();

	for(i=0; i<SpeciesNUM; i++){
		for(j=0; j<SpeciesNUM; j++){
			DAG[phylo_order[i]][phylo_order[j]]=MM[i][j];
		}
	}

	for(i=SpeciesNUM-1;i>=0;i--)
		delete[] MM[i];

	if(verbose_mode == VB_MAX){
	// Print info about synchronization
		cout<<endl<<"Synchronization:"<<endl;
		cout<<"PhyloInfo id | FoodWeb id, name"<<endl;
		for(i=0; i<SpeciesNUM; i++){
			cout<<"["<<phylo_order[i]<<"] | ["<<i<<"] "<<dagNames[i]<<endl;
		}
		cout<<"PhyloInfo id | name"<<endl;
		for(i=0; i<TaxaNUM;i++){
			cout<<"["<<i<<"] "<<findNodeID(i)->name<<endl;
		}

	 // Input data after processing: cannibalism, rescaling, reordering
		cout<<endl<<"Food web is defined by the following matrix"<<endl;
		for(i=0;i<nvar;i++) {
			if(findFoodWebID(i) != -1)
				cout<<dagNames[findFoodWebID(i)]<<"\t";
			else
				cout<<"\t\t";
			for(j=0;j<nvar;j++)
				cout<<DAG[i][j]<<"\t";
			cout<<endl;
		}
	}
/* ---------------------------------------------------------------------------------------------------------
 * Filling out taxaDAG vector: node corresponds to taxa, neighbors to preys, length (node-neighbor) to weight
 * ---------------------------------------------------------------------------------------------------------*/
 	vector<int> vec2;//the value of vec[j] is the height of the species in the DAG
	taxaDAG.resize(nvar,NULL);
	for(j=0;j<nvar;j++){
		taxaDAG[j] = newNode(j,j);
		//cout<<"taxonDAG[j="<<j+1<<"]->id="<<taxaDAG[j]->id<<endl;
	}

	for(j=0;j<nvar;j++){
		for(i=0;i<nvar;i++)
			if(DAG[i][j]>0){
				//cout<<"cheking matrix"<<i<<j<<endl;
				taxaDAG[j]->addNeighbor(taxaDAG[i], DAG[i][j], taxaDAG[i]->id);
				//cout<<"neighbors[i="<<taxaDAG[j]->degree()-1<<"]->id="<<taxaDAG[j]->neighbors[taxaDAG[j]->degree()-1]->node->id<<endl;
			}
		//cout<<endl;
	}

/* ---------------------------------------------------------------------------------------------------------
 * Defining levels in the Food Web based on the longest food chain of predators
 * ---------------------------------------------------------------------------------------------------------*/
	for(j=0;j<nvar;j++){
		levelDAG.push_back(0);
		if(taxaDAG[j]->degree()>0)
			vec2.push_back(1);
		else
			vec2.push_back(0);
//  		if(taxaDAG[j]->degree()>0){
//  			cout<<"Children of taxonDAG[j="<<j<<"]->id="<<taxaDAG[j]->id<<":"<<endl;
// 			for(i=0;i<taxaDAG[j]->degree();i++)
// 				cout<<"taxaDAG["<<j<<"]->neighbors["<<i<<"]->node->id "<<taxaDAG[j]->neighbors[i]->node->id<<endl;
// 				//cout<<"id of the child "<<i<<" node id "<<taxaDAG[j]->neighbors[i]->node->id+1<<" "<<endl;
// 				//cout<<"           neighbors[i="<<i<<"]->id="<<taxaDAG[j]->neighbors[i]->node->id<<endl;
// 			cout<<endl;
//
//  		}
	}
//	for(j=0;j<nvar;j++)
//		cout<<j<<" "<<levelDAG[j]<<" "<<vec2[j]<<endl;

	int eq=0,step=0;
	//cout<<"Starting while..."<<endl;
	while(eq!=1){
		eq=1;
		step++;
// 		if(step==1 or step==2 or step==3)
//		cout<<"-------STEP "<<step<<"-------"<<endl<<"j v1 v2"<<endl;
		for(j=0;j<nvar;j++){
			if(levelDAG[j]!=vec2[j])
				eq=0;
// 			if(step==1 or step==2 or step==3)
//			cout<<j<<" "<<levelDAG[j]<<" "<<vec2[j]<<endl;
			levelDAG[j]=vec2[j];
		}
		for(j=0;j<nvar;j++){
			if(taxaDAG[j]->degree()>0){
			//cout<<"taxaDAG["<<j<<"]->neighbors[0]->node->id "<<taxaDAG[j]->neighbors[0]->node->id<<endl;
			vec2[j]=vec2[taxaDAG[j]->neighbors[0]->node->id]+1;
			for(i=1;i<taxaDAG[j]->degree();i++)
				if(vec2[taxaDAG[j]->neighbors[i]->node->id]>=vec2[j])
					vec2[j]=vec2[taxaDAG[j]->neighbors[i]->node->id]+1;
			}
		}
	}

// For each predator the level corresponds to its longest food chain----------------------------
	if(verbose_mode == VB_MAX){
		cout<<"For each species its longest chain according to a food web"<<endl;
		for(j=0;j<nvar;j++)
			//if(findFoodWebID(j) != -1)
			//	cout<<dagNames[findFoodWebID(j)]<<"\t| "<<levelDAG[j]<<endl;
			//else
				cout<<*names[j]<<"\t| "<<levelDAG[j]<<endl;
	}
 	//cout<<"Species - level"<<endl;
	//ofstream levelF;
	//levelF.open("Level",ios::app);
 	//for(j=0;j<SpeciesNUM;j++)
 	//	levelF<<j+1<<" "<<levelDAG[j]<<endl;
	// 	for(i=0;i<tree.leafNum;i++)
	// 	myfile<<"taxon id: "<<taxaTree[i]->id<<" | taxon name: "<<taxaTree[i]->name<<endl;
	// 	myfile<<"root  id: "<<root->id<<" | root  name: "<<root->name<<endl;
	// // myfile.close();

// The maximum level is the longest food chain of the food web ---------------------------------
//	 int maxlevel;
//	 maxlevel=0;
//	 for(i=0;i<SpeciesNUM;i++)
//		if(maxlevel<levelDAG[i])
//			maxlevel=levelDAG[i];

// Decrease SpeciesNUM since you do not need to include the root to the Species anymore---------
	if(rooted)
		SpeciesNUM--;
}

/* =========================================================================================================
 *	ROOTED TREES
 * =========================================================================================================*/
void ECOpd::printECOlpRooted(const char* fileOUT,ECOpd &tree){
	ofstream out;
	out.exceptions(ios::failbit | ios::badbit);
	out.open(fileOUT);
	int m,i,j;
	//int step=0,step_all=0;
	//cout<<"# of species to conserve:"<<nspecies<<endl;
	int nspecies=k;
	nspecies++; //you have to include also one place for the root
//----------------------------------------------------------------------------------------------------------------
// Dealing with d levels
//----------------------------------------------------------------------------------------------------------------
// max d level-------------------------------------------------
	int maxlevel;
	maxlevel=levelDAG[0]; //i=0;
// 	cout<<"DAG levels:"<<endl;
// 	cout<<"LevelDAG[0]:"<<levelDAG[0]<<endl;
 	for(i=1;i<nvar;i++){
//  		//cout<<"LevelDAG["<<i<<"]:"<<levelDAG[i]<<endl;
 		if(maxlevel<levelDAG[i])
 			maxlevel=levelDAG[i];
	}
	//cout<<"max DAG level:"<<maxlevel+1<<endl;
//# of species at each d level---------------------------------
// 	int hit=0;
// 	hvec.resize(maxlevel+1,0);
// 	while(hit<=maxlevel){
// 		for(i=0;i<nvar;i++)
// 			if(levelDAG[i]==hit) hvec[hit]++;
// 		hit++;
// 	}
	/*cout<<"# of species at each level"<<endl;
	for(i=0;i<hvec.size();i++)
		cout<<hvec[i]<<" ";
	cout<<endl;
	*/

 /****************************************************************************************************************
  * Integer Programming formulation
  ****************************************************************************************************************/
//----------------------------------------------------------------------------------------------------------------
// Printing objective function
//----------------------------------------------------------------------------------------------------------------
	out<<"Maximize"<<endl;
	tree.getBranchOrdered(nodes1,nodes2);
	for(i=0;i<tree.branchNum;i++){
		nodes1[i]->findNeighbor(nodes2[i])->id=i;
		nodes2[i]->findNeighbor(nodes1[i])->id=i;
		if(i<tree.branchNum-1)
			out<<nodes1[i]->findNeighbor(nodes2[i])->length<<" "<<"y"<<i<<" + ";
		else
			out<<nodes1[i]->findNeighbor(nodes2[i])->length<<" "<<"y"<<i<<endl;
		}
//----------------------------------------------------------------------------------------------------------------
// Printing constraints
//----------------------------------------------------------------------------------------------------------------
	out<<"Subject To"<<endl;
//----------------------------------------------------------------------------------------------------------------
// 1. constraint: species present in the set
	if(initialTaxa.size()!=0)
		for(m=0;m<initialTaxa.size();m++)
			out<<"x"<<findSpeciesIDname(&initialTaxa[m])<<" = 1"<<endl;
//----------------------------------------------------------------------------------------------------------------
// 2. constraint: the sum of all species is <= k
	for(i=0;i<nvar-1;i++)
		out<<"x"<<i<<" + ";
	out<<"x"<<nvar-1<<" <= "<<nspecies<<endl;

//----------------------------------------------------------------------------------------------------------------
// 3. constraint: the sum of leaves in the DAG is >= to 1
// 	SpeciesNUM++;
// 	int nleafDAG=0,nleaf=0;
// 	for(i=0;i<SpeciesNUM;i++)
// 		if(levelDAG[i]==0)
// 			nleafDAG++;
// 		for(j=0;j<SpeciesNUM;j++){
// 			if(taxaDAG[j]->degree()==0){
// 				nleaf++;
// 			if(nleaf<nleafDAG)
// 				out<<"x"<<taxaDAG[j]->id<<" + ";
// 			else
// 				out<<"x"<<taxaDAG[j]->id<<" >= 1"<<endl;
// 			}
// 		}

//----------------------------------------------------------------------------------------------------------------
// 4. constraints: SURVIVAL CONSTRAINT
	if(weighted){
		//weighted food web: sum of weights is greater than a given threshold--------------------------------
		for(j=0;j<nvar;j++)
			if(taxaDAG[j]->degree()>0){//the ones that have children in the DAG
				for(i=0;i<taxaDAG[j]->degree();i++){
					if(i<taxaDAG[j]->degree()-1){
						out<<taxaDAG[j]->neighbors[i]->length<<" x"<<taxaDAG[j]->neighbors[i]->node->id<<" + ";
					} else {
						out<<taxaDAG[j]->neighbors[i]->length<<" x"<<taxaDAG[j]->neighbors[i]->node->id<<" - "<<T<<" x"<<taxaDAG[j]->id<<" >= 0"<<endl;
					}
				}
			}
	}else{
		//for each predator the sum of children in the DAG is >= to its value-------------------------------
		for(j=0;j<nvar;j++)
			if(taxaDAG[j]->degree()>0){//the ones that have children in the DAG
				for(i=0;i<taxaDAG[j]->degree();i++){
					if(i<taxaDAG[j]->degree()-1){
						out<<"x"<<taxaDAG[j]->neighbors[i]->node->id<<" + ";
					} else {
						out<<"x"<<taxaDAG[j]->neighbors[i]->node->id<<" - x"<<taxaDAG[j]->id<<" >= 0"<<endl;
					}
				}
			}
	}
//----------------------------------------------------------------------------------------------------------------
// 5. constraints for edges in the PhyloTree
	//cout<<"root "<<tree.root->id<<endl;
	vector<int> taxaBelow;
	for(i=0;i<tree.branchNum;i++)
		//constraints: SUM{Xv in T(e)}(Xv)>=Ye -----------------------------------------------
		if((nodes1[i]->isLeaf()) && (nodes1[i]!=root))
			out<<"x"<<nodes1[i]->id<<" - y"<<nodes1[i]->findNeighbor(nodes2[i])->id<<" >= 0"<<endl;
		else {
			tree.getTaxaID(taxaBelow,nodes2[i],nodes1[i]);
			for(j=0;j<taxaBelow.size();j++)
				if(j<taxaBelow.size()-1)
					out<<"x"<<taxaBelow[j]<<" + ";
				else
					out<<"x"<<taxaBelow[j];
			taxaBelow.clear();
			out<<" - y"<<nodes1[i]->findNeighbor(nodes2[i])->id<<" >= 0"<<endl;
		}
//----------------------------------------------------------------------------------------------------------------
// Printing bounds for variables
//----------------------------------------------------------------------------------------------------------------
	out<<"Bounds"<<endl;
	for(j=0;j<nvar;j++)
			out<<"0 <= x"<<taxaDAG[j]->id<<" <= 1"<<endl;
	for(i=0;i<tree.branchNum;i++)
		out<<"0 <= y"<<i<<" <= 1"<<endl;
//----------------------------------------------------------------------------------------------------------------
// Printing variables (For IP model)
//----------------------------------------------------------------------------------------------------------------
	out<<"Generals"<<endl;
	for(j=0;j<nvar;j++)
		out<<"x"<<taxaDAG[j]->id<<" ";
	for(i=0;i<tree.branchNum;i++)
		out<<"y"<<i<<" ";
//----------------------------------------------------------------------------------------------------------------
	out<<endl<<"End"<<endl;
	out.close();
}

/* =========================================================================================================
 *	UNROOTED TREES and d-levels
 * =========================================================================================================*/
void ECOpd::printECOlpUnrooted(const char* fileOUT,ECOpd &tree){
	ofstream myfile;
	string str_out = fileOUT;
	string str_out_1,str_out_2;
	//myfile.open(fileOUT);

	int i,m,j,step=0,step_all=0;
	int nspecies=k;
//---------------------------------------------Dealing with d levels---------------------------------------------
	{
//--------------------------------------------------max d level--------------------------------------------------
	int maxlevel;
	maxlevel=levelDAG[0];
// 	cout<<"DAG levels:"<<endl;
// 	cout<<"LevelDAG[0]:"<<levelDAG[0]<<endl;
	for(i=1;i<TaxaNUM;i++){
// 		cout<<"LevelDAG["<<i<<"]:"<<levelDAG[i]<<endl;
		if(maxlevel<levelDAG[i])
			maxlevel=levelDAG[i];
	}
// 	cout<<"max DAG level:"<<maxlevel+1<<endl;
//-----------------------------------------------------end-------------------------------------------------------


//---------------------------------------generating first vector of d levels-------------------------------------
	generateFirstMultinorm(dvec, nspecies-1, maxlevel+1); //nspecies-1 - we are saving at least 1 place at d[0] level
							      //maxlevel+1 - as we start counting levels from 0 one
// 	cout<<"vector of d levels:"<<endl;
// 	dvec[0]++;
// 	for(i=0;i<dvec.size();i++)
// 		cout<<dvec[i]<<" ";
// 	cout<<endl;
// 	dvec[0]--;
//-----------------------------------------------------end-------------------------------------------------------


//------------------------------------------# of species at each d level-----------------------------------------
	int hit=0;
	hvec.resize(maxlevel+1,0);
	while(hit<=maxlevel){
		for(i=0;i<TaxaNUM;i++)
			if(levelDAG[i]==hit) hvec[hit]++;
		hit++;
	}
	/*cout<<"# of species at each level"<<endl;
	for(i=0;i<hvec.size();i++)
		cout<<hvec[i]<<" ";
	cout<<endl;
	*/
//-----------------------------------------------------end-------------------------------------------------------
	}
//-----------------------------------------END:Dealing with d levels---------------------------------------------

	int DAGlevels=1,check_print=1;

//--------------------------------------Printing all cases for different d levels--------------------------------
	while(DAGlevels==1) {
		step_all++;
		//cout<<endl<<"STEP ALL:"<<step_all<<endl;
// 		cout<<"vector of d levels:"<<endl;
 		dvec[0]++;
// 		for(i=0;i<dvec.size();i++)
// 			cout<<dvec[i]<<" ";
// 		cout<<endl;

//--------------------------------------CHECKPOINT:is vector d good or we should ignore it?----------------------
{
//Print only if d[i] is <= than the # of species on this level, otherwise it's a waste of places for conservation
		check_print=1;
		for(i=0;i<hvec.size();i++){
			//cout<<"dvec["<<i<<"]="<<dvec[i]<<"   hvec["<<i<<"]="<<hvec[i]<<endl;
			if(dvec[i]>hvec[i])
				check_print=0;
		}
		//cout<<"CHECKPOINT="<<check_print<<endl<<endl;

		check_print=1;//IGNORE checkpoint: when only one model for each run is needed, used together with DAGlevels=0; (below)

		if(check_print==1){
			step++;
			//cout<<"Vector d is SUITABLE, step="<<step<<endl;
			str_out_1 = convertIntToString(step);
			str_out_2 = str_out  + str_out_1 + ".lp";
			//myfile.open(str_out_2.c_str());

			//str_out_2 = str_out + "lp";		//IGNORED d levels: only one model for each run
			str_out_2 = str_out;
			myfile.open(str_out_2.c_str());

/**----------------------------------------------Printing objective function---------------------------------------*/
	{

			myfile<<"Maximize"<<endl;
			tree.getBranchOrdered(nodes1,nodes2);
	{
			for(i=0;i<tree.branchNum;i++){
				nodes1[i]->findNeighbor(nodes2[i])->id=i;
				nodes2[i]->findNeighbor(nodes1[i])->id=i;
				if(i<tree.branchNum-1)
					myfile<<nodes1[i]->findNeighbor(nodes2[i])->length<<" "<<"y"<<i<<" + ";
				else
					myfile<<nodes1[i]->findNeighbor(nodes2[i])->length<<" "<<"y"<<i<<endl;
			}
	}
			//IDEA: objective**********************************************************************
	{
// 			double lambda_sum=0;
// 			for(i=0;i<tree.branchNum;i++){
// 				nodes1[i]->findNeighbor(nodes2[i])->id=i;
// 				nodes2[i]->findNeighbor(nodes1[i])->id=i;
// 				lambda_sum = lambda_sum + nodes1[i]->findNeighbor(nodes2[i])->length;
// 				if(i<tree.branchNum-1)
// 					myfile<<nodes1[i]->findNeighbor(nodes2[i])->length<<" "<<"y"<<i<<" + ";
// 				else
// 					myfile<<nodes1[i]->findNeighbor(nodes2[i])->length<<" "<<"y"<<i;
// 			}
//
// 			for(j=0;j<tree.leafNum;j++)
// 				if(taxaDAG[j]->degree()>0)
// 					myfile<<" - "<<lambda_sum<<" z"<<j;
// 			myfile<<endl;
	}
			//*************************************************************************************
}

/**--------------------------------------------------Printing constraints------------------------------------------*/
	{
			myfile<<"Subject To"<<endl;
			int c=0;
/**species present in the set-----------------------------------------------*/
	if(initialTaxa.size()!=0)
		for(m=0;m<initialTaxa.size();m++)
			myfile<<"x"<<findSpeciesIDname(&initialTaxa[m])<<" = 1"<<endl;
/**the sum of all species is <= k-------------------------------------------------------------------*/
	{
			NodeVector taxaTree;
			tree.getTaxa(taxaTree);
			myfile<<"c"<<c<<": ";
			c++;

			for(i=0;i<nvar-1;i++)
				myfile<<"x"<<i<<" + ";
			myfile<<"x"<<nvar-1<<" <= "<<nspecies<<endl;
	}
/**the sum of leaves in the DAG is >= to 1----------------------------------------------------------*/
	{ //it might be incorrect, so check it out before using...
// 			myfile<<"c"<<c<<": ";
// 			c++;
// 			int nleafDAG=0,nleaf=0;
// 			for(i=0;i<TaxaNUM;i++)
// 				if(levelDAG[i]==0)
// 					nleafDAG++;
// 			for(j=0;j<TaxaNUM;j++){
// 				if(taxaDAG[j]->degree()==0){
// 				nleaf++;
// 				if(nleaf<nleafDAG)
// 					myfile<<"x"<<taxaDAG[j]->id<<" + ";
// 				else
// 					myfile<<"x"<<taxaDAG[j]->id<<" >= 1"<<endl;
// 				}
// 			}
	}

//constraints: SURVIVAL CONSTRAINT

if(weighted){//weighted food web: sum of weights is greater than a given threshold--------------------------------

	for(j=0;j<nvar;j++){
		if(taxaDAG[j]->degree()>0){//the ones that have children in the DAG
			for(i=0;i<taxaDAG[j]->degree();i++)
				if(i<taxaDAG[j]->degree()-1)
					myfile<<taxaDAG[j]->neighbors[i]->length<<" x"<<taxaDAG[j]->neighbors[i]->node->id<<" + ";
				else
					myfile<<taxaDAG[j]->neighbors[i]->length<<" x"<<taxaDAG[j]->neighbors[i]->node->id<<" - "<<T<<" x"<<taxaDAG[j]->id<<" >= 0"<<endl;
			}

	}
}
else {//for each predator the sum of children in the DAG is >= to its value-----------------------------
			for(j=0;j<TaxaNUM;j++){
				if(taxaDAG[j]->degree()>0){//the ones that have children in the DAG
					myfile<<"c"<<c<<": ";
					c++;
					for(i=0;i<taxaDAG[j]->degree();i++)
						if(i<taxaDAG[j]->degree()-1)
							myfile<<"x"<<taxaDAG[j]->neighbors[i]->node->id<<" + ";
						else
							myfile<<"x"<<taxaDAG[j]->neighbors[i]->node->id<<" - x"<<taxaDAG[j]->id<<" >= 0"<<endl;
				}
			}


}
			//IDEA:new variables Z'tas*************************************************************
	{

// 			for(j=0;j<tree.leafNum;j++){
// 				if(taxaDAG[j]->degree()>0){
// 					myfile<<"c"<<c<<": ";
// 					c++;
// 					myfile<<"z"<<j<<" - x"<<j<<" >= 0"<<endl;
// 					for(i=0;i<taxaDAG[j]->degree();i++){
// 						myfile<<"c"<<c<<": ";
// 						c++;
// 						myfile<<"z"<<j<<" - x"<<taxaDAG[j]->neighbors[i]->node->id<<" >= 0"<<endl;
// 					}
// 				}
// 			}

}
			//*************************************************************************************

//d levels----------------------------------------------------------------------------------------
	{


// 			int hit=0;
// 			int h=0;
// 			int maxlevel=dvec.size()-1;
// 			while(hit<=maxlevel){
// 				myfile<<"c"<<c<<": ";
// 				c++;
// 				h=0;
// 				for(i=0;i<SpeciesNUM;i++)
// 					if(levelDAG[i]==hit){
// 						h++;
// 						if(h<hvec[hit])
// 							myfile<<"x"<<i<<" + ";
// 						else {
// 							//myfile<<"x"<<i<<" = "<<dvec[hit];
// 							myfile<<"x"<<i<<" - d"<<hit<<" = 0";
//  							if(step==1)
//  								cout<<"dvec["<<hit<<"]="<<dvec[hit]<<endl;
// 						}
// 					}
// 				myfile<<endl;
// 				hit++;
// 			}
// 			myfile<<"c"<<c<<": ";
// 			c++;
// 			for(i=0;i<=maxlevel;i++){
// 				if(i<maxlevel)
// 					myfile<<"d"<<i<<" + ";
// 				else
// 					myfile<<"d"<<i<<" = "<<nspecies<<endl;
// 			}


}

//for edges in the PhyloTree--------------------------------------------------------------------------
	{

			vector<int> taxaBelow;
			for(i=0;i<tree.branchNum;i++){

		//constraints: SUM{Xv in T(e)}(Xv)>=Ye -----------------------------------------------
				myfile<<"c"<<c<<": ";
				c++;
				tree.getTaxaID(taxaBelow,nodes2[i],nodes1[i]);
				for(j=0;j<taxaBelow.size();j++)
					if(j<taxaBelow.size()-1)
						myfile<<"x"<<taxaBelow[j]<<" + ";
					else
						myfile<<"x"<<taxaBelow[j];
				taxaBelow.clear();
				myfile<<" - y"<<nodes1[i]->findNeighbor(nodes2[i])->id<<" >= 0"<<endl;
		//constraints: SUM{Xv not in T(e)}(Xv)>=Ye -------------------------------------------
				myfile<<"c"<<c<<": ";
				c++;
				tree.getTaxaID(taxaBelow,nodes1[i],nodes2[i]);
				for(j=0;j<taxaBelow.size();j++)
					if(j<taxaBelow.size()-1)
						myfile<<"x"<<taxaBelow[j]<<" + ";
					else
						myfile<<"x"<<taxaBelow[j];
				taxaBelow.clear();
				myfile<<" - y"<<nodes1[i]->findNeighbor(nodes2[i])->id<<" >= 0"<<endl;
			}

	}//end printing constraints for edges (PhyloTree)



}//-------------------------------------END printing constraints ALL--------------------------------------
	//int maxlevel=dvec.size()-1;


//---------------------------------------------Printing bounds for variables--------------------------------------
	{
			myfile<<"Bounds"<<endl;
			for(j=0;j<nvar;j++)
				myfile<<"0 <= x"<<taxaDAG[j]->id<<" <= 1"<<endl;
			for(i=0;i<tree.branchNum;i++)
				myfile<<"0 <= y"<<i<<" <= 1"<<endl;
	}
//------------------------------------------Printing variables (For IP model)---------------------------------
	{
 			myfile<<"Generals"<<endl;
			for(j=0;j<nvar;j++)
				myfile<<"x"<<taxaDAG[j]->id<<" ";
			for(i=0;i<tree.branchNum;i++)
				myfile<<"y"<<i<<" ";

//			myfile<<"Generals"<<endl;
// 			int maxlevel=dvec.size()-1;
// 			for(j=0;j<=maxlevel;j++)
// 				myfile<<"d"<<j<<" ";
// 			if(fractVAR.size()!=0)
// 				for(i=0;i<fractVAR.size();i++)
// 					myfile<<fractVAR[i]<<" ";
			myfile<<endl;
	}
 			myfile<<"End"<<endl;
			myfile.close();
		}//checkpoint "if" ends here
}
//----------------------------------------------CHECKPOINT ENDs here----------------------------------------------


//---------------------------------------generating next vector of d levels--------------------------------------
		dvec[0]--;
		if(generateNextMultinorm(dvec))
			DAGlevels=1;
		else
			DAGlevels=0;
//-----------------------------------------------------end-------------------------------------------------------

		DAGlevels=0;//IGNORE d levels: only one model to solve for each run
	}
//----------------------------------END:Printing all cases for different d levels--------------------------------

//cout<<"ALL STEPS:"<<step_all<<endl;
//cout<<"STEPs:"<<step<<endl;

}




/* =========================================================================================================
 * SPLIT systems
 * =========================================================================================================*/
void ECOpd::printInfDAG (const char* fileOUT,PDNetwork &splitsys, Params &params) {
	ofstream out;
	out.exceptions(ios::failbit | ios::badbit);
	out.open(fileOUT,ios::app);
	int i,j,nspecies=k;
	int maxlevel;
	maxlevel=levelDAG[0];
	for(i=1;i<TaxaNUM;i++)
		if(maxlevel<levelDAG[i])
			maxlevel=levelDAG[i];
	int hit=0;
	hvec.resize(maxlevel+1,0);
	while(hit<=maxlevel){
		for(i=0;i<TaxaNUM;i++)
			if(levelDAG[i]==hit) hvec[hit]++;
		hit++;
	}
//Constraints----------------------------------------------------------------------
	//species present in the set-----------------------------------------------
	if (initialTaxa.size() != 0) {
		for (i = 0; i < initialTaxa.size(); i++) {
			out << "x" << findSpeciesIDname(&initialTaxa[i]) << " = 1" << endl;
		}
	}
	//the sum of all species is <= k-------------------------------------------
	for (i = 0; i < nvar - 1; i++) {
		out << "x" << i << " + ";
	}
	out<<"x"<<nvar-1<<" <= "<<nspecies<<endl;
	//the sum of leaves in the DAG is >= to 1----------------------------------
	int nleafDAG=0,nleaf=0;
	for (i = 0; i < nvar; i++) {
		if (levelDAG[i] == 0) {
			nleafDAG++;
		}
	}
	for(j=0;j<nvar;j++){
		if (taxaDAG[j]->degree() == 0) {
			nleaf++;
			if (nleaf < nleafDAG) {
				out << "x" << taxaDAG[j]->id << " + ";
			}
			else {
				out << "x" << taxaDAG[j]->id << " >= 1" << endl;
			}
		}
	}
	//SURVIVAL CONSTRAINT
	if(weighted){
		//constraint: Weighted food web. Sum of weights is greater than a given threshold--------------------------------
		for(j=0;j<nvar;j++){
			if(taxaDAG[j]->degree()>0){//the ones that have children in the DAG
				for(i=0;i<taxaDAG[j]->degree();i++)
					if(i<taxaDAG[j]->degree()-1)
						out<<taxaDAG[j]->neighbors[i]->length<<" x"<<taxaDAG[j]->neighbors[i]->node->id<<" + ";
					else
						out<<taxaDAG[j]->neighbors[i]->length<<" x"<<taxaDAG[j]->neighbors[i]->node->id<<" - "<<T<<" x"<<taxaDAG[j]->id<<" >= 0"<<endl;
				}
		}
	} else {
		//for each predator the sum of children in the DAG is >= to its value-------
		for(j=0;j<nvar;j++)
			if(taxaDAG[j]->degree()>0){//the ones that have children in the DAG
				for(i=0;i<taxaDAG[j]->degree();i++){
					if(i<taxaDAG[j]->degree()-1){
						out<<"x"<<taxaDAG[j]->neighbors[i]->node->id<<" + ";
					} else {
						out<<"x"<<taxaDAG[j]->neighbors[i]->node->id<<" - x"<<taxaDAG[j]->id<<" >= 0"<<endl;
					}
				}
			}
	}

//Bounds-----------------------------------------------------------------------------
	Split included_tax(TaxaNUM);
	IntVector::iterator it2;
// 	for (it2 = splitsys.initialset.begin(); it2 != splitsys.initialset.end(); it2++)
// 		included_tax.addTaxon(*it2);
	vector<int> y_value;
	y_value.resize(splitsys.getNSplits(), -1);
	//splitsys.checkYValue(nspecies, y_value);
	splitsys.lpVariableBound(out, params, included_tax, y_value);

//Generals for IP or MIP--------------------------------------------------------------
			out<<"Generals"<<endl;
 			for(i=0;i<nvar;i++)
 				out<<"x"<<i<<" ";
 			for(i=0;i<splitsys.getNSplits();i++)
 				out<<"y"<<i<<" ";
// 			for(j=0;j<=maxlevel;j++)
// 				out<<"d"<<j<<" ";
/*			if(fractVAR.size()!=0)
				for(i=0;i<fractVAR.size();i++)
					out<<fractVAR[i]<<" ";*/
			out<<endl;
			out<<"End"<<endl;

	out.close();

//	ofstream out1;
//	out1.exceptions(ios::failbit | ios::badbit);
//	out1.open("variablesNUM.data",ios::app);
//	out1<<this->k<<" "<<nvar<<" "<<splitsys.getNSplits()<<endl;//" "<<maxlevel+1<<endl;
//	out1.close();
}

//Fractional stuff-----------------------------------------------------------------
void ECOpd::readREC(const char* infile) {
	ifstream in;
	cout<<endl<<"-----------------------------------------------------"<<endl;
	cout<<"Reading file with fractional variables from "<<infile<<endl;
	try {
		in.exceptions(ios::failbit | ios::badbit);
		in.open(infile);
		in.exceptions(ios::badbit);
		readREC(in);
		in.close();
	} catch (const char* str) {
		outError(str);
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, infile);
	}
}

void ECOpd::readREC(istream &in) {
	int i;
	string str,name;
	while (getline(in, str)) {
		stringstream ss(str);
		getline(ss,name,':');
		fractVAR.push_back(name);
	}
	for(i=0; i<fractVAR.size();i++)
	cout<<fractVAR[i]<<endl;
}

//Generating all d vectors ---------------------------------------------------------
void ECOpd::generateFirstMultinorm(vector<int> &x, int n, int k) {
     x.resize(k, 0);
     x.back() = n;
}

bool ECOpd::generateNextMultinorm(vector<int> &x) {
     if (x.size() < 2) return false;
     int id = x.size()-1;
     while (id >= 0 && x[id] == 0) id--;
     if (id <= 0) return false;
     x[id-1]++;
     x.back() = x[id]-1;
     if (id < x.size()-1) x[id] = 0;
     return true;
}

void ECOpd::getBranchOrdered(NodeVector &nodes, NodeVector &nodes2, Node *node, Node *dad){
	if(!node) node = root;
	FOR_NEIGHBOR_IT(node, dad, it){
		nodes.push_back(node);
		nodes2.push_back((*it)->node);
		getBranchOrdered(nodes,nodes2,(*it)->node,node);
	}
}

void ECOpd::synchTreeDAG(ECOpd &tree){
	if(rooted)
		tree.root->id=SpeciesNUM;

	for(int i=0; i<SpeciesNUM;i++){
		if(tree.findLeafName(dagNames[i]))
			tree.findLeafName(dagNames[i])->id = i;
	}
}

int ECOpd::findPhyloID(string name){
	for(int i=0; i<TaxaNUM; i++)
		if((phyloNames[i]).compare(name) == 0) return(i);
	return(-1);
}


int ECOpd::findFoodWebID(int id){
	for(int i=0; i<phylo_order.size();i++){
		if(phylo_order[i] == id) return i;
	}
	return(-1);
}

void ECOpd::randomBranLenTrees(Params &params){
	ECOpd tree = *this;
	//Trees with random branch length---------------------------------------------------------------------
	NodeVector nodes_1,nodes_2;
	tree.getBranchOrdered(nodes_1,nodes_2);
	for(int i=0;i<tree.branchNum;i++){
		if(nodes_1[i]!=tree.root && nodes_2[i]!=tree.root){
			nodes_1[i]->findNeighbor(nodes_2[i])->id=i;
			nodes_2[i]->findNeighbor(nodes_1[i])->id=i;
			nodes_1[i]->findNeighbor(nodes_2[i])->length=randomLen(params);
			nodes_2[i]->findNeighbor(nodes_1[i])->length=nodes_1[i]->findNeighbor(nodes_2[i])->length;
			//cout<<"Branch: y"<<i<<"-> length = "<<nodes_1[i]->findNeighbor(nodes_2[i])->length<<endl;
		} else {
			nodes_1[i]->findNeighbor(nodes_2[i])->id=i;
			nodes_2[i]->findNeighbor(nodes_1[i])->id=i;
			nodes_1[i]->findNeighbor(nodes_2[i])->length=1;
			nodes_2[i]->findNeighbor(nodes_1[i])->length=nodes_1[i]->findNeighbor(nodes_2[i])->length;
			//cout<<"Branch: y"<<i<<"-> length = "<<nodes_1[i]->findNeighbor(nodes_2[i])->length<<"              y"<<i<<" - ROOT edge"<<endl;
		}
	}
	string str_out = params.user_file;
	string str_out_1,str_out_2, str1, str2;

	str_out_1 = convertIntToString(params.eco_run);
	if(params.k_percent)
		str1 = convertIntToString(params.k_percent);
	else
		str1 = convertIntToString(params.sub_size);
	str2 = convertIntToString(params.diet_max);
	str_out_2 = str_out  + "." + str1 + "." + str2 + "." + str_out_1;
	const char *outfile=str_out_2.c_str();
	tree.printTree(outfile);
}

void ECOpd::detectMissingSpecies(){
	int i;
	for(i=0; i<TaxaNUM; i++)
		if(!findTaxaDAG(i)){
			missInDAG.push_back(phyloNames[i]);
			OUTdagTaxa.push_back(i+1);
		}
	if(missInDAG.size() != 0){
		cout<<endl<<"-------------------------------------------------------------------"<<endl;
		cout<<" There are "<<missInDAG.size()<<" species missing in the Food Web: "<<endl;
		cout<<"-------------------------------------------------------------------"<<endl;
		for(i=0; i<missInDAG.size(); i++)
			cout<<missInDAG[i]<<endl;
		cout<<endl;
	}


	for(i=0; i<SpeciesNUM; i++)
		if(!findSpeciesPhylo(i)){
			missInPhylo.push_back(dagNames[i]);
			OUTtreeTaxa.push_back(i+1);
		}
	if(missInPhylo.size() != 0){
		cout<<endl<<"-----------------------------------------------------------------------------"<<endl;
		cout<<" There are "<<missInPhylo.size()<<" species missing on the Tree/SplitSystem: "<<endl;
		cout<<"-----------------------------------------------------------------------------"<<endl;
		for(i=0; i<missInPhylo.size(); i++){
			cout<<missInPhylo[i]<<endl;
			names.push_back(&missInPhylo[i]);
		}
		cout<<endl;
	}
}

bool ECOpd::findTaxaDAG(int i){
	for(int j=0; j<SpeciesNUM; j++)
		if(phyloNames[i].compare(dagNames[j]) == 0)
			return true;
	return false;
}
bool ECOpd::findSpeciesPhylo(int i){
	for(int j=0; j<TaxaNUM; j++)
		if(dagNames[i].compare(phyloNames[j]) == 0)
			return true;
	return false;
}

void ECOpd::synchronizeSpecies(){
	int i;
	//cout<<"Synchronization starts..."<<endl;
	if(rooted)
		SpeciesNUM--;
	int num = 0;
	for(i=0; i<SpeciesNUM; i++){
		if(!OUT_tree(i+1)){
			if(phyloType == "t"){
				if(findLeafName(dagNames[i])){
					phylo_order.push_back(findLeafName(dagNames[i])->id);
				}
			}else{
				if(findPhyloID(dagNames[i]) != -1){
					phylo_order.push_back(findPhyloID(dagNames[i]));
				}
			}
		}else{
			phylo_order.push_back(TaxaNUM + num);
			num++;
		}
	}
	if(rooted){
		phylo_order.push_back(TaxaNUM-1);
		SpeciesNUM++;
	}
	// Filling out OUTtreeTaxa vector with new ids, based on phylo_order
	//cout<<"OUTtreeTaxa after reordering:"<<endl;
	for(i=0; i<OUTtreeTaxa.size(); i++){
		OUTtreeTaxa[i] = phylo_order[OUTtreeTaxa[i]-1];
		//cout<<OUTtreeTaxa[i]<<" ";
	}
	//cout<<endl;
}

int ECOpd::findSpeciesIDname(string *name){
	for(int i=0; i<nvar; i++){
		if(name->compare(*names[i]) == 0)
			return i;
	}
	return -1;
}

void ECOpd::defineK(Params &params){
	cout<<"Defining the subset size, k..."<<endl;
	if(rooted)
		nvar--;

	if(params.k_percent)
		k = params.k_percent*0.01*nvar;
	else if(params.sub_size)
		k = params.sub_size;

	if(k<2){
		cout<<"k = "<<k<<endl;
		cout<<"ERROR: Wrong value of parameter k. The subset size must be larger than 1."<<endl;
		exit(0);
	}else if(k>nvar){
		cout<<"k = "<<k<<endl;
		cout<<"Total number of species in the analysis | "<<nvar<<endl;
		cout<<"ERROR: Wrong value of parameter k. The subset size must be less or equal to the number of all species in the analysis."<<endl;
		exit(0);
	}
	cout<<"k = "<<k<<endl;
	if(initialTaxa.size() > k){
		cout<<endl<<"Initial set "<<initialTaxa.size()<<" taxa | Subset size k = "<<k<<endl;
		cout<<"ERROR: the initial set is already larger than the specified subset size! Increase k or reduce the initial set."<<endl;
		exit(0);
	}

	if(rooted)
		nvar++;
	if(T != 0){
		cout<<"Defining the minimum diet, d..."<<endl;
		cout<<"d = "<<(int) (T*100)<<endl;
	}

}

void ECOpd::checkInitialTaxa(){
	int i = 0, j = 0;
	vector<int> eraseSET;
	if(initialTaxa.size() != 0){
		cout<<"Reading taxa to be included in the final optimal subset.."<<endl;
		//for(i=0; i<initialTaxa.size(); i++)
		//	cout<<initialTaxa[i]<<endl;
		for(i=0; i<initialTaxa.size(); i++){
			if(findSpeciesIDname(&initialTaxa[i]) == -1){
				j++;
				if(j == 1){
					cout<<"---------------------------------------------------------------------------------------------------------"<<endl;
					cout<<"The following species are not present on Tree/SplitSystem nor in the Food Web, therefore will be ignored:"<<endl;
					cout<<"---------------------------------------------------------------------------------------------------------"<<endl;
				}
				cout<<initialTaxa[i]<<endl;
				eraseSET.push_back(i);
			}
		}
		cout<<endl;
		if(eraseSET.size()!=0)
			for(i = eraseSET.size()-1; i>=0; i--){
				initialTaxa.erase(initialTaxa.begin() + eraseSET[i]);
			}
		cout<<"------------------------------------------"<<endl;
		cout<<"Taxa to be included in the optimal subset:"<<endl;
		cout<<"------------------------------------------"<<endl;
		for(i=0; i<initialTaxa.size(); i++)
			cout<<initialTaxa[i]<<endl;
		//cout<<"The initial subset size is"<<initialTaxa.size()<<endl;
		cout<<endl;
	}
}

void ECOpd::printSubFoodWeb(char* fileOUT, double* variables){
	ofstream out;
	out.exceptions(ios::failbit | ios::badbit);
	out.open(fileOUT);
	int i,j;
	out<<k<<endl;
	for(i=0; i<nvar; i++){
		if(variables[i] == 1){
			out<<*(names[i])<<" ";
			for(j=0; j<nvar; j++)
				if(variables[j] == 1)
					out<<DAG[i][j]<<" ";
			out<<endl;
		}
	}
	out.close();
}

void ECOpd::dietConserved(double *variables){
	int i,j;
	double c;
	for(i=0; i<nvar; i++){
		c = 0;
		if(variables[i] == 1){
			for(j=0; j<nvar; j++){
				if(variables[j] == 1)
					c += DAG[j][i];
			}
		}
		c *= 100;
		dietVAL.push_back(c);
	}
}

void ECOpd::printResults(char* fileOUT,double* variables, double score,Params &params){
	cout<<endl<<"Results of the analysis are printed to "<<fileOUT<<endl<<endl;
	ofstream out;
	out.exceptions(ios::failbit | ios::badbit);
	out.open(fileOUT);
	int i;

	if(phyloType == "t")
		summarizeHeader(out, params, false, IN_NEWICK);
	else
		summarizeHeader(out, params, false, IN_NEXUS);

	// Analyze the results of IP and print the information
	out<<endl<<"------Results of biodiversity analysis--------------------------------------------------"<<endl;
	out<<endl;
	out<<"Food Web  | # of species "<<SpeciesNUM<<"\t| # of links    "<<linksNUM;
	//printf("Food Web  | # of species  %3i | # of links    %7i ", SpeciesNUM,linksNUM);
	if(weighted)
		out<<"\t| weighted \t|"<<endl;
	else
		out<<"\t| non weighted \t|"<<endl;

	if(phyloType == "t"){
		if(rooted){
			TaxaNUM--;
		}
		//printf("PhyloTree | # of species  %3i | # of branches %7i ", TaxaNUM,branchNum);
		out<<"PhyloTree | # of species "<<TaxaNUM<<"\t| # of branches "<<branchNum;
		if(rooted)
			out<<"\t| rooted \t|";
		else
			out<<"\t| unrooted \t|";
		out<<" total PD "<<treeLength()<<endl;
	}else{
		//printf("SplitSys  | # of species  %3i | # of splits %7i | total SD %f \n", ecoInfDAG.TaxaNUM,splitSYS.getNSplits(),splitSYS.calcWeight());
		out<<"SplitSys  | # of species "<<TaxaNUM<<"\t| # of splits   "<<splitsNUM<<"\t|\t\t| total SD "<<totalSD<<endl;
	}

	out<<endl;
	out<<"SubsetSize| k = "<<k<<endl;
	if(T!= 0)
		out<<"Constraint| "<<(int) (T*100)<<"%-viability"<<endl;
	else
		out<<"Constraint| naive viability"<<endl;
	out<<endl;

	if(phyloType == "t"){
		out<<"PD of the optimal subset: "<<score<<" (constitutes "<< score / treeLength() * 100<<"% of the total PD)"<<endl;
	} else {
		out<<"SD of the optimal subset: "<<score<<" (constitutes "<< score / totalSD * 100<<"% of the total SD)"<<endl;
	}

	if(weighted){
		out<<"--------------------------------------------------"<<endl;
		out<<" Optimal subset of species  (% of diet conserved) "<<endl;
		out<<"--------------------------------------------------"<<endl;
		if(rooted){
			for(i=0; i<nvar; i++)
				if(variables[i] == 1 && i!=root->id){
					if(dietVAL[i]!=0){
						out<<" ";
						out.width(30);
						out<<left<<*(names[i]);
						out<<"\t("<<dietVAL[i]<<"%)"<<endl;
					}else
						out<<" "<<*(names[i])<<endl;
				}
		}else{
			for(i=0; i<nvar; i++)
				if(variables[i] == 1){
					if(dietVAL[i]!=0){
						out<<" ";
						out.width(30);
						out<<left<<*(names[i]);
						out<<"\t("<<dietVAL[i]<<"%)"<<endl;
					}else
						out<<" "<<*(names[i])<<endl;
				}
		}
	}else{
		out<<"-----------------------------"<<endl;
		out<<" Optimal subset of species  "<<endl;
		out<<"-----------------------------"<<endl;
		if(rooted){
			for(i=0; i<nvar; i++)
				if(variables[i] == 1 && i!=root->id)
					out<<" "<<*(names[i])<<endl;
		}else{
			for(i=0; i<nvar; i++)
				if(variables[i] == 1)
					out<<" "<<*(names[i])<<endl;
		}
	}
	//out<<"----------------------------------------------------------------------------------------"<<endl;
	summarizeFooter(out, params);
	out.close();
}
