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

/*
	collection of classes for Next-generation sequencing 
*/

#include "ngs.h"
#include "modeltest_wrapper.h"

/****************************************************************************
        NGSAlignment
 ****************************************************************************/

NGSAlignment::NGSAlignment(const char *filename) : AlignmentPairwise() {
	readFritzFile(filename);
}

NGSAlignment::NGSAlignment(int nstate, int ncat, int *freq) : AlignmentPairwise() {
	num_states = nstate;
	ncategory = ncat;
	int total_size = ncategory*num_states*num_states;
	pair_freq = new int[total_size];
	memcpy(pair_freq, freq, total_size * sizeof(int));
}

NGSAlignment::NGSAlignment(int nstate, string &seq1, string &seq2) {
	num_states = nstate;
	ncategory = 1;
	pair_freq = new int[nstate*nstate];
	memset(pair_freq, 0, sizeof(int)*nstate*nstate);
	assert(seq1.length() == seq2.length());
	int len = seq1.length();
	int i;
	for (i = 0; i < len; i++) {
		int state1 = convertState(seq1[i], SEQ_DNA);
		int state2 = convertState(seq2[i], SEQ_DNA);
		if (state1 < num_states && state2 < num_states) 
			pair_freq[state1*num_states+state2] ++;
	}
}

void NGSAlignment::readFritzFile(const char *filename) {
	cout << "Reading Fritz file " << filename << " ..." << endl;
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(filename);
		in.clear();
		int i, total_size;
		string tmp;
		in >> tmp;
		ncategory = convert_int(tmp.c_str());
		if (ncategory < 1) throw "Wrong number of positions";
		in >> tmp;
		num_states = convert_int(tmp.c_str());
		total_size = ncategory*num_states*num_states;
		if (num_states < 1) throw "Wrong number of states";
		pair_freq = new int[total_size];
		for (i=0; i < total_size; i++) {
			in >> tmp;
			int count = convert_int(tmp.c_str());
			if (count < 0) throw "Wrong count";
			pair_freq[i] = count;
		}
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (const char *str) {
		outError(str);
	} catch (string str) {
		outError(str);
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	}

	cout << ncategory << " matrices of size " << num_states << endl;
}

void NGSAlignment::computeStateFreq (double *stateFrqArr) {
	int cat, i, j, id = 0;
	int state_count[num_states];
	memset(state_count, 0, sizeof(int)*num_states);
	for (cat = 0, id = 0; cat < ncategory; cat++) {
		for (i = 0; i < num_states; i++)
			for (j = 0; j < num_states; j++, id++) {
				state_count[i] += pair_freq[id];
				state_count[j] += pair_freq[id];
			}
	}

	int sum_count = 0;
	for (i = 0; i < num_states; i++) sum_count += state_count[i];
	if (sum_count == 0) throw "Empty data observed";
	for (i = 0; i < num_states; i++) stateFrqArr[i] = double(state_count[i]) / sum_count;
	if (verbose_mode >= VB_MIN) {
		cout << "Empirical state frequencies: ";
		for (i = 0; i < num_states; i++) 
			cout << stateFrqArr[i] << " ";
		cout << endl;
	}
}

void NGSAlignment::computeSumPairFreq (int *sum_pair_freq) {
	int cat, id, i, j;
	memset(sum_pair_freq, 0, sizeof(int)*num_states*num_states);
	for (cat = 0, id = 0; cat < ncategory; cat++) {
		for (i = 0; i < num_states; i++)
			for (j = 0; j < num_states; j++, id++) {
				sum_pair_freq[i*num_states+j] += pair_freq[id];
			}
	}
}

void NGSAlignment::computeEmpiricalRate (double *rates) {
	int i, j, k, cat, id;
	assert(rates);
	double **pair_rates = (double**) new double[num_states];
	for (i = 0; i < num_states; i++) {
		pair_rates[i] = new double[num_states];
		memset(pair_rates[i], 0, sizeof(double)*num_states);
	}

	for (cat = 0, id = 0; cat < ncategory; cat++) {
		for (i = 0; i < num_states; i++)
			for (j = 0; j < num_states; j++, id++) {
				pair_rates[i][j] += pair_freq[id];
			}
	}

	k = 0;
	double last_rate = pair_rates[num_states-2][num_states-1] + pair_rates[num_states-1][num_states-2];
	if (last_rate == 0.0) throw "Last rate entry is ZERO";
	for (i = 0; i < num_states-1; i++)
		for (j = i+1; j < num_states; j++)
			rates[k++] = (pair_rates[i][j] + pair_rates[j][i]) / last_rate;
	if (verbose_mode >= VB_MIN) {
		cout << "Empirical rates: ";
		for (k = 0; k < num_states*(num_states-1)/2; k++)
			cout << rates[k] << " ";
		cout << endl;
	}
	for (i = num_states-1; i >= 0; i--) {
		delete [] pair_rates[i];
	}
	delete [] pair_rates;
}

double NGSAlignment::computeEmpiricalDist(int cat) {
	int i;
	int trans_size = num_states*num_states;
	int *pair_pos = pair_freq + (cat*trans_size);
	int match_pos = 0, total_pos = 0;
	for (i = 0; i < num_states; i++) 
		match_pos += pair_pos[i*num_states+i];
	for (i = 0; i < trans_size; i++) 
		total_pos += pair_pos[i];
	if (total_pos == 0) total_pos = 1;
	return (double)(total_pos - match_pos) / total_pos;
}


double NGSAlignment::computeFunctionCat(int cat, double value) {
	int trans_size = num_states*num_states;
	double lh = 0.0;
	double trans_mat[trans_size];
	int i;

	tree->getModelFactory()->computeTransMatrix(value, trans_mat);
	int *pair_pos = pair_freq + cat*trans_size;

	for (i = 0; i < trans_size; i++) if (pair_pos[i]) {
		if (trans_mat[i] <= 0) throw "Negative transition probability";
		lh -= pair_pos[i] * log(trans_mat[i]);
	}
	return lh;
}


double NGSAlignment::computeFuncDervCat(int cat, double value, double &df, double &ddf) {
	int trans_size = num_states*num_states;
	double lh = 0.0;
	df = 0.0; ddf = 0.0;
	int i;
	double derv1 = 0.0, derv2 = 0.0;
	double trans_mat[trans_size], trans_derv1[trans_size], trans_derv2[trans_size];
	

	tree->getModelFactory()->computeTransDerv(value, trans_mat, trans_derv1, trans_derv2);
	int *pair_pos = pair_freq + cat*trans_size;
	for (i = 0; i < trans_size; i++) if (pair_pos[i] > 0) {
		if (trans_mat[i] <= 0) throw "Negative transition probability";
		double d1 = trans_derv1[i] / trans_mat[i];
		derv1 += pair_pos[i] * d1;
		derv2 += pair_pos[i] * (trans_derv2[i]/trans_mat[i] - d1 * d1);
		lh -= pair_pos[i] * log(trans_mat[i]);
	}
	//df -= derv1 * rate_val;
	//ddf -= derv2 * rate_val * rate_val;
	df -= derv1;
	ddf -= derv2;
	return lh;
}

/****************************************************************************
        NGSRate
 ****************************************************************************/
NGSRate::NGSRate(PhyloTree *tree) {
	phylo_tree = tree;
	ncategory = ((NGSAlignment*)tree->aln)->ncategory;
	rates = new double[ncategory];
	int i;
	for (i = 0; i < ncategory; i++) {
		rates[i] = ((NGSAlignment*)tree->aln)->computeEmpiricalDist(i);
		if (rates[i] < 1e-6) rates[i] = 1e-6;
	}
	
	name = "+F";
	name += convertIntToString(ncategory);
	full_name = name;
	is_categorized = true;

}

double NGSRate::optimizeParameters() {
	int cat;
	double negative_lh;
	for (cat = 0; cat < ncategory; cat++) {
		optimizing_cat = cat;
		if (phylo_tree->optimize_by_newton) 
			rates[cat] = minimizeNewtonSafeMode(1e-6, rates[cat], 10.0, 1e-6, negative_lh);
		else
			rates[cat] = minimizeOneDimenSafeMode(1e-6, rates[cat], 10.0, 1e-6, &negative_lh);
	}
	return phylo_tree->computeLikelihood();
}


double NGSRate::computeFunction(double value) {
	return ((NGSAlignment*)phylo_tree->aln)->computeFunctionCat(optimizing_cat, value);
}
double NGSRate::computeFuncDerv(double value, double &df, double &ddf) {
	return ((NGSAlignment*)phylo_tree->aln)->computeFuncDervCat(optimizing_cat, value, df, ddf);
}

void NGSRate::writeInfo(ostream &out) {
}

/****************************************************************************
        NGSTree
 ****************************************************************************/

NGSTree::NGSTree(Params &params, NGSAlignment *alignment) {
	aln = alignment;
    model = NULL;
    site_rate = NULL;
    model_factory = NULL;
    optimize_by_newton = params.optimize_by_newton;
    //tree.sse = params.SSE;
    sse = false;
}

double NGSTree::computeLikelihood(double *pattern_lh) {
	return -((NGSAlignment*)aln)->computeFunction(1.0);
}

double NGSTree::optimizeAllBranches(int iterations) {
	return computeLikelihood();
}


/****************************************************************************
        NGSRead
 ****************************************************************************/


void NGSRead::init() {
	scaff.clear();
	read.clear();
	id = -2;
	match_pos= -2;
	flag=true;
	identity=-2;
	times=1.0;
	direction=true;
}


double NGSRead::computeFunction(double value) {

	RateHeterogeneity *site_rate = tree->getRate();
	int i, rate_id;
	int nptn = scaff.length();
	double lh = 0.0;

	// site-specific rates
	for (i = 0, rate_id = 0; i < nptn; i++) {
		int state1 = scaff[i];
		int state2 = read[i];
		if (state1 >= num_states || state2 >= num_states) continue;
		double rate_val = site_rate->getRate(rate_id);
		if (homo_rate > 0)
			rate_val = homo_rate;
		double trans = tree->getModelFactory()->computeTrans(value * rate_val, state1, state2);
		lh -= log(trans);
		rate_id++;
	}
	return lh;
}

double NGSRead::computeFuncDerv(double value, double &df, double &ddf) {
	RateHeterogeneity *site_rate = tree->getRate();
	int i, rate_id;
	int nptn = scaff.length();
	double lh = 0.0;
	df = 0.0; ddf = 0.0;

	for (i = 0, rate_id = 0; i < nptn; i++) {
		int state1 = scaff[i];
		int state2 = read[i];
		if (state1 >= num_states || state2 >= num_states) continue;
		double rate_val = site_rate->getRate(rate_id);
		if (homo_rate > 0)
			rate_val = homo_rate;
		double rate_sqr = rate_val * rate_val;
		double derv1, derv2;
		double trans = tree->getModelFactory()->computeTrans(value * rate_val, state1, state2, derv1, derv2);
		lh -= log(trans);
		double d1 = derv1 / trans;
		df -= rate_val * d1;
		ddf -= rate_sqr * (derv2/trans - d1*d1);
		rate_id++;	
	}
	return lh;
}

/****************************************************************************
        NGSReadSet
 ****************************************************************************/

void reverseComplement(string &str) {
	string out;
	out.resize(str.length(), ' ');
	string::reverse_iterator it;
	string::iterator oit;
	for (it = str.rbegin(), oit = out.begin(); it != str.rend(); it++, oit++) {
		char c = toupper(*it);
		switch (c) {
			case 'A': *oit = 'T'; break;
			case 'T': *oit = 'A'; break;
			case 'G': *oit = 'C'; break;
			case 'C': *oit = 'G'; break;
			default: *oit = *it; break;
		}
	}
	//cout << str << endl << out << endl;
	str = out;
}

//("File","total",0.8,-1)
void NGSReadSet::parsNextGen(string filename, string ref_ID,double ident,int mismatches) 
{
//	cout<<"start"<<endl;
	string a= "total";
	size_t buffer_size = 1200;
	ifstream myfile; //test2
	myfile.open(filename.c_str(),ifstream::in);
	if(!myfile.good()){
		cout<<"No such file "<<filename.c_str()<<endl;
		exit(0);
	}
	char* line = new char[buffer_size];
//	cout<<"start"<<endl;

	NGSRead tempread;
	tempread.init();
	tempread.tree = tree;
	tempread.num_states = tree->aln->num_states;

	ReadInfo read_info;

	myfile.getline(line,buffer_size);
	string ref;
	for (; !myfile.eof(); ) {
		if(line[0]=='S'&& line[1]=='e'){
			for(size_t i=0;i<buffer_size;i++){
				if(line[i]=='\0' ||line[i]=='\n' ){
					break;
				}
				if(tempread.id ==-2 && strncmp(&line[i],"ID: ",4)==0){
					tempread.id = atoi(&line[i+4]);
				}else if(tempread.id !=-2 && strncmp(&line[i],"ID: ",4)==0){
					int id = atoi(&line[i+4]);

					if(id==0){
						tempread.flag=true;
					}else{
						tempread.flag=false;
					}

				}else if(strncmp(&line[i],"forward",7)==0){
					tempread.direction=true;
//					cout<<i<<endl;
				}else if(strncmp(&line[i],"backward",8)==0){
					tempread.direction=false;
				}

				if(strncmp(&line[i],"me: ",4)==0){
					i=i+4;
					while(i<buffer_size&&line[i]!=' '){
						tempread.name+=line[i];
						i++;
					}
				}
				if(strncmp(&line[i],"re: ",4)==0){
					tempread.score= atoi(&line[i+4]);
					break;
				}

				if(strncmp(&line[i],"at: ",4)==0){
					tempread.match_pos= atoi(&line[i+4])+1;
				}
				if(strncmp(&line[i],"ld: ",4)==0){
					tempread.chr.clear();
					size_t t=i+4;
					while(t<buffer_size && line[t]!='\n' &&  line[t]!='\0'){
						//tempread.chr.size()>3 &&
						if( line[t]==' '){
							break;
						}
						tempread.chr+=line[t];
						t++;
					}
				}
			}

			if( (strcmp(tempread.chr.c_str(),ref_ID.c_str())==0 || strcmp(a.c_str(),ref_ID.c_str())==0 )){

				myfile.getline(line,buffer_size);
				for(size_t i=0;i<buffer_size;i++){
					if(line[i]=='\0' ||line[i]=='\n' ){
						break;
					}
					if(strncmp(&line[i],"es: ",4)==0){
						tempread.times= atof(&line[i+4]);
					}
					if(strncmp(&line[i],"ty: ",4)==0){
						tempread.identity=atof(&line[i+4]);
						break;
					}
				}

				if(tempread.identity>=ident){
					string scaff;
					string read;
					myfile.getline(line,buffer_size);
					size_t i=0;
					while(i<buffer_size &&line[i]!=' '  &&line[i]!='\t'&&line[i]!='\0'&&line[i]!='\n'){
						scaff+=line[i];
						i++;
					}

					myfile.getline(line,buffer_size);
					i=0;
					int count=0;
					while(i<buffer_size && line[i]!=' ' &&line[i]!='\t'&&line[i]!='\0'&&line[i]!='\n'){
						read+=line[i];
						if(line[i]!='-' && scaff[i]!='-' && scaff[i]!=line[i]){
							count++;
						}
						i++;
					}

					tempread.scaff=scaff;
					tempread.read=read;
					if (!tempread.direction) {
						reverseComplement(tempread.scaff);
						reverseComplement(tempread.read);
					}
					tempread.convertState(tempread.scaff, SEQ_DNA);
					tempread.convertState(tempread.read, SEQ_DNA);
					assert(tempread.scaff.length() == tempread.read.length());

					if(count==mismatches || mismatches < 0){
						tempread.homo_rate = homo_rate;
						read_info.homo_distance = tempread.optimizeDist(1.0-tempread.identity);
						read_info.homo_logl = -tempread.computeFunction(read_info.homo_distance);
						tempread.homo_rate = 0.0;
						read_info.distance = tempread.optimizeDist(read_info.homo_distance);
						read_info.logl = -tempread.computeFunction(read_info.distance);
						read_info.id = tempread.id;
						read_info.identity = tempread.identity;
						push_back(read_info);
					}
					scaff.clear();
					read.clear();
				}
			}
			tempread.chr.clear();
			tempread.init();
		}

		myfile.getline(line,buffer_size);
		if (size()>0 && size() % 10000 == 0) cout << size() << " reads processed" << endl;

	}

	cout << size() << " total reads processed" << endl;

	myfile.close();
	delete [] line;
}

void NGSReadSet::writeInfo() {
	cout << size() << " reads process in total" << endl;
	return;
}

/****************************************************************************
        main function
 ****************************************************************************/

void reportNGSAnalysis(const char *file_name, Params &params, NGSAlignment &aln, NGSTree &tree, 
	DoubleMatrix &rate_info, StrVector &rate_name) {
	ofstream out(file_name);
	out.setf(ios::fixed,ios::floatfield);

	int i, j, k;


	double rate_param[aln.num_states * aln.num_states];
	double rate_matrix[aln.num_states * aln.num_states];

	out << "Input file: " << params.ngs_file << endl;
	out << "Model of evolution: " << tree.getModel()->name << endl << endl;

	out << "Substitution process assuming one homogeneous model among all positions:" << endl;

	out << "Rate parameters: " << endl;

	tree.getModel()->getRateMatrix(rate_param);

	if (tree.getModel()->name == "UNREST") {
		for (i = 0, k=0; i < aln.num_states; i++)
			for (j = 0; j < aln.num_states; j++)
				if (i != j)
					rate_matrix[i*aln.num_states+j] = rate_param[k++];
	} else {
		for (i = 0, k=0; i < aln.num_states-1; i++)
			for (j = i+1; j < aln.num_states; j++, k++)
				rate_matrix[i*aln.num_states+j] = rate_matrix[j*aln.num_states+i] = rate_param[k];
	}

	for (i = 0; i < aln.num_states; i++) {
		for (j = 0; j < aln.num_states; j++) {
			if (j > 0) out << " \t";
			if (j != i) out << rate_matrix[i*aln.num_states+j]; else out << "-";
		}
		out << endl;
	}
	out << endl;
	out << "State frequencies: " << endl;

	double state_freq[aln.num_states];
	tree.getModel()->getStateFrequency(state_freq);

	for (i = 0; i < aln.num_states; i++) out << state_freq[i] << " \t";
	out << endl << endl;

	out << "Q matrix can be obtained by multiplying rate parameters with state frequencies" << endl << endl;

	out << "Log-likelihood: " << tree.computeLikelihood() << endl << endl;

	out << "Inferred posisiton-specific rates under one model or position-specific model: " << endl;

	out << "Position\tOne_model_rate";
	for (StrVector::iterator it = rate_name.begin(); it != rate_name.end(); it++)
		out << "\t" << (*it);
	out << endl;
	for (i = 0; i < aln.ncategory; i++) {
		out << i+1 << '\t' << tree.getRate()->getRate(i);
		DoubleVector *rate_vec = &rate_info[i];
		for (DoubleVector::iterator dit = rate_vec->begin(); dit != rate_vec->end(); dit++)
			out << "\t" << *dit;
		out << endl;
	}
	out.close();
	cout << endl << "Results written to: " << file_name << endl << endl;
}

bool checkFreq(int *pair_freq, int n) {
	int i, count = 0;
	for (i=0; i < n*n; i++)
		if (pair_freq[i] != 0) count++;
	if (count <= n) return false;
	return true;
}

void testSingleRateModel(Params &params, NGSAlignment &aln, NGSTree &tree, string model, 
	int *freq, DoubleVector &rate_info, StrVector &rate_name, bool write_info) {

	char model_name[20];
	NGSAlignment sum_aln(aln.num_states, 1, freq);

	NGSTree sum_tree(params, &sum_aln);
	sum_aln.tree = &sum_tree;

	if (model == "") 
		sprintf(model_name, "GTR+F1");
	else
		sprintf(model_name, "%s+F1", model.c_str());
	try {
		params.model_name = model_name;
		sum_tree.setModelFactory(new ModelFactory(params, &sum_tree));
		sum_tree.setModel(sum_tree.getModelFactory()->model);
		sum_tree.setRate(sum_tree.getModelFactory()->site_rate);
    	double bestTreeScore = sum_tree.getModelFactory()->optimizeParameters(false, write_info);
		cout << "LogL: " << bestTreeScore;
		cout << " / Rate: " << sum_tree.getRate()->getRate(0) << endl;
    } catch (...) {
		cout << "Skipped due to sparse matrix" << endl;
		//rate_info.push_back(MIN_SITE_RATE);
		rate_info.insert(rate_info.end(), rate_name.size(), MIN_SITE_RATE);
		return;
    }
    //return sum_tree.getRate()->getRate(0);
	rate_info.push_back(sum_tree.getRate()->getRate(0));

    double rate_mat[aln.num_states*aln.num_states];
    memset(rate_mat, 0, aln.num_states*aln.num_states*sizeof(double));
    sum_tree.getModel()->getRateMatrix(rate_mat);
    rate_info.insert(rate_info.end(), rate_mat, rate_mat+sum_tree.getModel()->getNumRateEntries());

	if (tree.getModel()->isReversible()) {
		sum_tree.getModel()->getStateFrequency(rate_mat);
		rate_info.insert(rate_info.end(), rate_mat, rate_mat+aln.num_states);
    }
}


/*

void testSingleRateModel(Params &params, NGSAlignment &aln, NGSTree &tree, string model, int *sum_freq) {
	char model_name[20];

	NGSAlignment sum_aln(aln.num_states, 1, sum_freq);

	NGSTree sum_tree(params, &sum_aln);
	sum_aln.tree = &sum_tree;

	if (model == "") 
		sprintf(model_name, "GTR+F1");
	else
		sprintf(model_name, "%s+F1", model.c_str());
	params.model_name = model_name;
    sum_tree.setModelFactory(new ModelFactory(params, &sum_tree));
    sum_tree.setModel(sum_tree.getModelFactory()->model);
    sum_tree.setRate(sum_tree.getModelFactory()->site_rate);

    double bestTreeScore = sum_tree.getModelFactory()->optimizeParameters(false, false);
    cout << "Log-likelihood of null model: " << bestTreeScore << endl;
    cout << "Rate (or distance) of null model: " << sum_tree.getRate()->getRate(0) << endl;
    double lh_diff = 2*(tree.computeLikelihood() - bestTreeScore);
    cout << "2(lnL1 - lnL0) = " << lh_diff << endl;
    cout << "p-value (chi-square test, df = " << aln.ncategory-1 << "): " << computePValueChiSquare(lh_diff, aln.ncategory-1) << endl;

	string out_file = params.out_prefix;
	out_file += ".ngs_e";
	DoubleVector tmp;
	reportNGSAnalysis(out_file.c_str(), params, sum_aln, sum_tree, tmp);

}*/

void reportNGSReads(const char *file_name, Params &params, NGSReadSet &ngs_reads)
{
	ofstream out(file_name);
	out.setf(ios::fixed,ios::floatfield);
	out << "Read\tHamm_dist\tHomo_dist\tHete_dist\tHomo_logl\tHete_logl" << endl;
	for (int i = 0; i < ngs_reads.size(); i++)
		out << ngs_reads[i].id << '\t' << 1.0 - ngs_reads[i].identity << 
		'\t' << ngs_reads[i].homo_distance << '\t' << ngs_reads[i].distance << 
		'\t' << ngs_reads[i].homo_logl << '\t' << ngs_reads[i].logl << endl;
	out.close();
	cout << endl << "Results written to: " << file_name << endl << endl;

}

void runNGSAnalysis(Params &params) {
	char model_name[20];
	// read input file, initialize NGSAlignment
	NGSAlignment aln(params.ngs_file);
	cout.setf(ios::fixed,ios::floatfield);

	params.freq_type = FREQ_ESTIMATE;

	// initialize NGSTree
	NGSTree tree(params, &aln);
	aln.tree = &tree;

	// initialize Model 
	string original_model = params.model_name;
	if (params.model_name == "") 
		sprintf(model_name, "GTR+F%d", aln.ncategory);
	else
		sprintf(model_name, "%s+F%d", params.model_name.c_str(), aln.ncategory);
	params.model_name = model_name;
    tree.setModelFactory(new ModelFactory(params, &tree));
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);

    int model_df = tree.getModel()->getNDim() + tree.getRate()->getNDim();
    cout << endl;
    cout << "Model of evolution: " << tree.getModelName() << " (" << model_df << " free parameters)" << endl;
    cout << endl;

	// optimize model parameters and rate scaling factors
    cout << "Optimizing model parameters" << endl;
    double bestTreeScore = tree.getModelFactory()->optimizeParameters(false, true);
    cout << "Log-likelihood: " << bestTreeScore << endl;


	DoubleMatrix part_rate(aln.ncategory);
	StrVector rate_name;


	int i, j;

	rate_name.push_back("Varying_model_rate");

	if (tree.getModel()->isReversible()) {
		for (i = 0; i < aln.num_states-1; i++) 
			for (j = i+1; j < aln.num_states; j++) {
				stringstream x;
				x << aln.convertStateBack(i) << "->" << aln.convertStateBack(j);
				rate_name.push_back(x.str());
			}
		for (i = 0; i < aln.num_states; i++) {
			stringstream x;
			x << aln.convertStateBack(i);
			rate_name.push_back(x.str());
		}
	} else {
		for (i = 0; i < aln.num_states; i++) 
			for (j = 0; j < aln.num_states; j++) if (j != i) {
				stringstream x;
				x << aln.convertStateBack(i) << "->" << aln.convertStateBack(j);
				rate_name.push_back(x.str());
			}
	}


	VerboseMode vb_saved = verbose_mode;
	verbose_mode = VB_QUIET;

	cout << endl << "--> INFERING RATE ASSUMING POSITION-SPECIFIC MODEL..." << endl << endl;
	for (int pos = 0; pos < aln.ncategory; pos++) {
		cout << "Position " << pos+1 << " / ";
		int *pair_pos = aln.pair_freq + (pos*aln.num_states*aln.num_states);
		testSingleRateModel(params, aln, tree, original_model, pair_pos, part_rate[pos], rate_name, false);
	}


	verbose_mode = vb_saved;

	int sum_freq[aln.num_states*aln.num_states];
	cout << endl << "-->INFERING RATE UNDER EQUAL-RATE NULL MODEL..." << endl << endl;
	aln.computeSumPairFreq(sum_freq);
	DoubleVector null_rate;
	testSingleRateModel(params, aln, tree, original_model, sum_freq, null_rate, rate_name, true);

	// report running results
	string out_file = params.out_prefix;
	out_file += ".ngs";
	reportNGSAnalysis(out_file.c_str(), params, aln, tree, part_rate, rate_name);

	if (!params.ngs_mapped_reads) return;
	
	NGSReadSet ngs_reads;
	ngs_reads.tree = &tree;
	ngs_reads.homo_rate = null_rate[0];
	cout << "Homogeneous rate: " << ngs_reads.homo_rate << endl;
	cout << "Reading mapped reads file " << params.ngs_mapped_reads << " ... " << endl;
	ngs_reads.parsNextGen(params.ngs_mapped_reads);
	ngs_reads.writeInfo();

	out_file = params.ngs_mapped_reads;
	out_file += ".dist";
	reportNGSReads(out_file.c_str(), params, ngs_reads);
}
