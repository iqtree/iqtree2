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
#include "model/modelmarkov.h"
//#include "modeltest_wrapper.h"

/****************************************************************************
        NGSAlignment
 ****************************************************************************/

NGSAlignment::NGSAlignment(PhyloTree *atree) : AlignmentPairwise() {
    tree = atree;
}

NGSAlignment::NGSAlignment(const char *filename) : AlignmentPairwise() {
    readFritzFile(filename);
}

NGSAlignment::NGSAlignment(int nstate, int ncat, double *freq) : AlignmentPairwise() {
    num_states = nstate;
    ncategory = ncat;
    int total_size = ncategory*num_states*num_states;
    pair_freq = new double[total_size];
    memcpy(pair_freq, freq, total_size * sizeof(double));
}

NGSAlignment::NGSAlignment(int nstate, string &seq1, string &seq2) {
    num_states = nstate;
    ncategory = 1;
    pair_freq = new double[nstate*nstate];
    memset(pair_freq, 0, sizeof(double)*nstate*nstate);
    ASSERT(seq1.length() == seq2.length());
    int len = seq1.length();
    int i;
    for (i = 0; i < len; i++) {
        int state1 = convertState(seq1[i], SEQ_DNA);
        int state2 = convertState(seq2[i], SEQ_DNA);
        if (state1 < num_states && state2 < num_states)
            pair_freq[state1*num_states+state2] += 1;
    }
}

char NGSAlignment::convertState(char state, SeqType seq_type) {
    char c = Alignment::convertState(state, SEQ_DNA);
    if (c == STATE_UNKNOWN) return 4;
    if (c >= 4) return 5;
    return c;
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
        pair_freq = new double[total_size];
        for (i=0; i < total_size; i++) {
            in >> tmp;
            double count = convert_double(tmp.c_str());
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
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    }

    cout << ncategory << " matrices of size " << num_states << endl;
}

void NGSAlignment::computeStateFreq (double *stateFrqArr, size_t num_unknown_states) {
    int cat, i, j, id = 0;
    double *state_count = new double[num_states];
    memset(state_count, 0, sizeof(double)*num_states);
    for (cat = 0, id = 0; cat < ncategory; cat++) {
        for (i = 0; i < num_states; i++)
            for (j = 0; j < num_states; j++, id++) {
                state_count[i] += pair_freq[id];
                state_count[j] += pair_freq[id];
            }
    }

    double sum_count = 0;
    for (i = 0; i < num_states; i++) sum_count += state_count[i];
    if (sum_count == 0) throw "Empty data observed";
    for (i = 0; i < num_states; i++) stateFrqArr[i] = double(state_count[i]) / sum_count;
    /*if (verbose_mode >= VB_MIN)*/ {
        cout << "Empirical state frequencies: ";
        for (i = 0; i < num_states; i++)
            cout << stateFrqArr[i] << " ";
        cout << endl;
    }
    delete [] state_count;
}

void NGSAlignment::computeSumPairFreq (double *sum_pair_freq) {
    int cat, id, i, j;
    memset(sum_pair_freq, 0, sizeof(double)*num_states*num_states);
    for (cat = 0, id = 0; cat < ncategory; cat++) {
        for (i = 0; i < num_states; i++)
            for (j = 0; j < num_states; j++, id++) {
                sum_pair_freq[i*num_states+j] += pair_freq[id];
            }
    }
}

void NGSAlignment::computeDivergenceMatrix (double *rates) {
    int i, j, k, cat, id;
    ASSERT(rates);
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
    /*if (verbose_mode >= VB_MIN)*/ {
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
    double *pair_pos = pair_freq + (cat*trans_size);
    double match_pos = 0, total_pos = 0;
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
    double *trans_mat = new double[trans_size];
    int i;

    tree->getModelFactory()->computeTransMatrix(value, trans_mat);
    double *pair_pos = pair_freq + cat*trans_size;

    for (i = 0; i < trans_size; i++) if (pair_pos[i] > 1e-6) {
            if (trans_mat[i] <= 0) throw "Negative transition probability";
            lh -= pair_pos[i] * log(trans_mat[i]);
        }
    delete [] trans_mat;
    return lh;
}


void NGSAlignment::computeFuncDervCat(int cat, double value, double &df, double &ddf) {
    int trans_size = num_states*num_states;
//    double lh = 0.0;
    df = 0.0;
    ddf = 0.0;
    int i;
    double derv1 = 0.0, derv2 = 0.0;
    double *trans_mat = new double[trans_size];
    double *trans_derv1 = new double[trans_size];
    double *trans_derv2 = new double[trans_size];


    tree->getModelFactory()->computeTransDerv(value, trans_mat, trans_derv1, trans_derv2);
    double *pair_pos = pair_freq + cat*trans_size;
    for (i = 0; i < trans_size; i++) if (pair_pos[i] > 1e-6) {
            if (trans_mat[i] <= 0) throw "Negative transition probability";
            double d1 = trans_derv1[i] / trans_mat[i];
            derv1 += pair_pos[i] * d1;
            derv2 += pair_pos[i] * (trans_derv2[i]/trans_mat[i] - d1 * d1);
//            lh -= pair_pos[i] * log(trans_mat[i]);
        }
    //df -= derv1 * rate_val;
    //ddf -= derv2 * rate_val * rate_val;
    df -= derv1;
    ddf -= derv2;
	delete [] trans_derv2;
	delete [] trans_derv1;
	delete [] trans_mat;
//    return lh;
    return;
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

double NGSRate::optimizeParameters(double epsilon) {
    int cat;
    double negative_lh;
    for (cat = 0; cat < ncategory; cat++) {
        optimizing_cat = cat;
        if (phylo_tree->optimize_by_newton)
            rates[cat] = minimizeNewtonSafeMode(1e-6, rates[cat], 10.0, max(epsilon,1e-6), negative_lh);
        else
            rates[cat] = minimizeOneDimenSafeMode(1e-6, rates[cat], 10.0, max(epsilon, 1e-6), &negative_lh);
    }
    return phylo_tree->computeLikelihood();
}


double NGSRate::computeFunction(double value) {
    return ((NGSAlignment*)phylo_tree->aln)->computeFunctionCat(optimizing_cat, value);
}
void NGSRate::computeFuncDerv(double value, double &df, double &ddf) {
    ((NGSAlignment*)phylo_tree->aln)->computeFuncDervCat(optimizing_cat, value, df, ddf);
}

void NGSRate::writeInfo(ostream &out) {
}

/****************************************************************************
        NGSRateCat
 ****************************************************************************/
NGSRateCat::NGSRateCat(PhyloTree *tree, int ncat) {
    phylo_tree = tree;
    ncategory = ncat;
    rates = new double[ncategory];
    proportion = new double[ncategory];
    int i;
    for (i = 0; i < ncategory; i++) {
        rates[i] = random_double();
        proportion[i] = 1.0/ncategory;
		
    }

    sum_pair_freq = new double[tree->aln->num_states * tree->aln->num_states];
    ((NGSAlignment*)tree->aln)->computeSumPairFreq(sum_pair_freq);

    name = "+FC";
    name += convertIntToString(ncategory);
    full_name = name;
    is_categorized = true;
}


/**
	return the number of dimensions
*/
int NGSRateCat::getNDim() {
    return 2*ncategory-1;
}

void NGSRateCat::setVariables(double *variables) {
    memcpy(variables+1, rates, ncategory * sizeof(double));
    memcpy(variables+ncategory+1, proportion, (ncategory-1)*sizeof(double));
}

bool NGSRateCat::getVariables(double *variables) {
    bool changed = memcmpcpy(rates, variables+1, ncategory * sizeof(double));
    changed |= memcmpcpy(proportion, variables+ncategory+1, (ncategory-1)*sizeof(double));
    double sum = 0.0;
    for (int i = 0; i < ncategory-1; i++)
        sum += proportion[i];
    proportion[ncategory-1] = 1.0 - sum;
    return changed;
}


/**
	the target function which needs to be optimized
	@param x the input vector x
	@return the function value at x
*/
double NGSRateCat::targetFunk(double x[]) {
    getVariables(x);
    if (proportion[ncategory-1] <= 1e-6) return 1e9;
    return -phylo_tree->computeLikelihood();
}


double NGSRateCat::optimizeParameters(double epsilon) {
    int ndim = getNDim();

    // return if nothing to be optimized
    if (ndim == 0) return 0.0;

    cout << "Optimizing " << name << " model parameters..." << endl;


    double *variables = new double[ndim+1];
    double *upper_bound = new double[ndim+1];
    double *lower_bound = new double[ndim+1];
    bool *bound_check = new bool[ndim+1];
    int i;
    double score;

    // by BFGS algorithm
    setVariables(variables);
    for (i = 1; i <= ndim; i++) {
        //cout << variables[i] << endl;
        lower_bound[i] = 1e-4;
        upper_bound[i] = 100.0;
        bound_check[i] = false;
    }
    for (i = ndim-ncategory+2; i <= ndim; i++)
        upper_bound[i] = 1.0;
    //packData(variables, lower_bound, upper_bound, bound_check);
    score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(epsilon, 1e-6));

    getVariables(variables);

    delete [] bound_check;
    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;

    return score;
}


void NGSRateCat::writeInfo(ostream &out) {
    int i;
    double avg = 0.0;
    out << "Rates: ";
    for (i = 0; i < ncategory; i++)
        out << " " << rates[i];
    out << endl << "Proportion: ";
    for (i = 0; i < ncategory; i++)
        out << " " << proportion[i];
    out << endl;
    for (i = 0; i < ncategory; i++)
        avg += rates[i]*proportion[i];
    cout << avg << endl;
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
    setLikelihoodKernel(LK_386, params.num_threads);
}

double NGSTree::computeLikelihood(double *pattern_lh) {
    return -((NGSAlignment*)aln)->computeFunction(1.0);
}

double NGSTree::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
    return computeLikelihood();
}


/****************************************************************************
        NGSTreeCat
 ****************************************************************************/

NGSTreeCat::NGSTreeCat(Params &params, NGSAlignment *alignment) : NGSTree(params, alignment) {
}

double NGSTreeCat::computeLikelihood(double *pattern_lh) {
    int num_states = getModel()->num_states;
    int trans_size = num_states*num_states;
    double *sum_trans_mat = new double[trans_size];
    double *trans_mat = new double[trans_size];
    int i, cat;
    NGSRateCat *site_rate = (NGSRateCat*)getRate();

    memset(sum_trans_mat, 0, trans_size * sizeof(double));
    for (cat = 0; cat < site_rate->getNDiscreteRate(); cat++) {
        getModel()->computeTransMatrix(site_rate->getRate(cat), trans_mat);
        for (i = 0; i < trans_size; i++)
            sum_trans_mat[i] += site_rate->proportion[cat]*trans_mat[i];
    }
    double lh = 0.0;
    for (i = 0; i < trans_size; i++)
        lh += ((NGSAlignment*)aln)->pair_freq[i] * log(sum_trans_mat[i]);
    delete [] trans_mat;
    delete [] sum_trans_mat;
    return lh;
}


/****************************************************************************
        NGSRead
 ****************************************************************************/

NGSRead::NGSRead(PhyloTree *atree) : NGSAlignment(atree) {
    init();
    if (tree) {
        num_states = tree->aln->num_states;
    } else num_states = 4;
    pair_freq = new double[(num_states+1) * (num_states+1)];
}

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

void NGSRead::computePairFreq() {
    int len = scaff.length();
    ASSERT(len == read.length());
    memset(pair_freq, 0, sizeof(double)*num_states*num_states);
    for (int i = 0; i < len; i++)
        if (scaff[i] < num_states && read[i] < num_states)
            pair_freq[scaff[i]*num_states+read[i]] += 1;
}


double NGSRead::computeFunction(double value) {

    RateHeterogeneity *site_rate = tree->getRate();
    int i, rate_id;
    int nptn = scaff.length();
    double lh = 0.0;
    if (homo_rate > 0.0) {
        int trans_size = num_states*num_states;
        double *trans_mat = new double[trans_size];
        tree->getModelFactory()->computeTransMatrix(value * homo_rate, trans_mat);
        for (i = 0; i < trans_size; i++) if (pair_freq[i] > 1e-6) {
                lh -= pair_freq[i] * log(trans_mat[i]);
            }
        delete [] trans_mat;
        return lh;
    }
    // site-specific rates
    for (i = 0, rate_id = 0; i < nptn; i++) {
        int state1 = scaff[i];
        int state2 = read[i];
        if (state1 >= num_states || state2 >= num_states) continue;
        double trans;
        double rate_val = site_rate->getRate(rate_id);
        trans = tree->getModelFactory()->computeTrans(value * rate_val, state1, state2);
        lh -= log(trans);
        rate_id++;
    }
    return lh;
}

void NGSRead::computeFuncDerv(double value, double &df, double &ddf) {
    RateHeterogeneity *site_rate = tree->getRate();
    int i, rate_id;
    int nptn = scaff.length();
//    double lh = 0.0;
    df = 0.0;
    ddf = 0.0;

    if (homo_rate > 0.0) { // homogeneous rate
        int trans_size = num_states*num_states;
        double *trans_mat = new double[trans_size];
        double *trans_derv1 = new double[trans_size];
        double *trans_derv2 = new double[trans_size];
        tree->getModelFactory()->computeTransDerv(value * homo_rate, trans_mat, trans_derv1, trans_derv2);
        for (i = 0; i < trans_size; i++) if (pair_freq[i] > 1e-6) {
//                lh -= pair_freq[i] * log(trans_mat[i]);
                double d1 = trans_derv1[i] / trans_mat[i];
                df -=  pair_freq[i] * d1;
                ddf -= pair_freq[i] * (trans_derv2[i]/trans_mat[i] - d1*d1);
            }
        df *= homo_rate;
        ddf *= homo_rate * homo_rate;
        delete [] trans_derv2;
        delete [] trans_derv1;
        delete [] trans_mat;
//        return lh;
        return;
    }

    for (i = 0, rate_id = 0; i < nptn; i++) {
        int state1 = scaff[i];
        int state2 = read[i];
        if (state1 >= num_states || state2 >= num_states) continue;
        double trans, derv1, derv2;
        double rate_val = site_rate->getRate(rate_id);
        double rate_sqr = rate_val * rate_val;
        trans = tree->getModelFactory()->computeTrans(value * rate_val, state1, state2, derv1, derv2);
//        lh -= log(trans);
        double d1 = derv1 / trans;
        df -= rate_val * d1;
        ddf -= rate_sqr * (derv2/trans - d1*d1);
        rate_id++;
    }

//    return lh;
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
        //*oit = c;
        switch (c) {
        case 'A':
            *oit = 'T';
            break;
        case 'T':
            *oit = 'A';
            break;
        case 'G':
            *oit = 'C';
            break;
        case 'C':
            *oit = 'G';
            break;
        default:
            *oit = c;
            break;
        }
    }
    //cout << str << endl << out << endl;
    str = out;
}

//("File","total",0.8,-1)
void NGSReadSet::parseNextGen(string filename, string ref_ID,double ident,int mismatches)
{
//	cout<<"start"<<endl;
    string a= "total";
    size_t buffer_size = 1200;
    ifstream myfile; //test2
    myfile.open(filename.c_str(),ifstream::in);
    if (!myfile.good()) {
        cout<<"No such file "<<filename.c_str()<<endl;
        exit(0);
    }
    char* line = new char[buffer_size];
//	cout<<"start"<<endl;

    NGSRead tempread(tree);

    myfile.getline(line,buffer_size);
    string ref;
    for (; !myfile.eof(); ) {
        if (line[0]=='S'&& line[1]=='e') {
            for (size_t i=0;i<buffer_size;i++) {
                if (line[i]=='\0' ||line[i]=='\n' ) {
                    break;
                }
                if (tempread.id ==-2 && strncmp(&line[i],"ID: ",4)==0) {
                    tempread.id = atoi(&line[i+4]);
                } else if (tempread.id !=-2 && strncmp(&line[i],"ID: ",4)==0) {
                    int id = atoi(&line[i+4]);

                    if (id==0) {
                        tempread.flag=true;
                    } else {
                        tempread.flag=false;
                    }

                } else if (strncmp(&line[i],"forward",7)==0) {
                    tempread.direction=true;
//					cout<<i<<endl;
                } else if (strncmp(&line[i],"backward",8)==0) {
                    tempread.direction=false;
                }

                if (strncmp(&line[i],"me: ",4)==0) {
                    i=i+4;
                    while (i<buffer_size&&line[i]!=' ') {
                        tempread.name+=line[i];
                        i++;
                    }
                }
                if (strncmp(&line[i],"re: ",4)==0) {
                    tempread.score= atoi(&line[i+4]);
                    break;
                }

                if (strncmp(&line[i],"at: ",4)==0) {
                    tempread.match_pos= atoi(&line[i+4])+1;
                }
                if (strncmp(&line[i],"ld: ",4)==0) {
                    tempread.chr.clear();
                    size_t t=i+4;
                    while (t<buffer_size && line[t]!='\n' &&  line[t]!='\0') {
                        //tempread.chr.size()>3 &&
                        if ( line[t]==' ') {
                            break;
                        }
                        tempread.chr+=line[t];
                        t++;
                    }
                }
            }

            if ( (strcmp(tempread.chr.c_str(),ref_ID.c_str())==0 || strcmp(a.c_str(),ref_ID.c_str())==0 )) {

                myfile.getline(line,buffer_size);
                for (size_t i=0;i<buffer_size;i++) {
                    if (line[i]=='\0' ||line[i]=='\n' ) {
                        break;
                    }
                    if (strncmp(&line[i],"es: ",4)==0) {
                        tempread.times= atof(&line[i+4]);
                    }
                    if (strncmp(&line[i],"ty: ",4)==0) {
                        tempread.identity=atof(&line[i+4]);
                        break;
                    }
                }

                if (tempread.identity>=ident) {
                    string scaff;
                    string read;
                    myfile.getline(line,buffer_size);
                    size_t i=0;
                    while (i<buffer_size &&line[i]!=' '  &&line[i]!='\t'&&line[i]!='\0'&&line[i]!='\n') {
                        scaff+=line[i];
                        i++;
                    }

                    myfile.getline(line,buffer_size);
                    i=0;
                    int count=0;
                    while (i<buffer_size && line[i]!=' ' &&line[i]!='\t'&&line[i]!='\0'&&line[i]!='\n') {
                        read+=line[i];
                        if (line[i]!='-' && scaff[i]!='-' && scaff[i]!=line[i]) {
                            count++;
                        }
                        i++;
                    }

                    tempread.scaff=scaff;
                    tempread.read=read;

                    if (count==mismatches || mismatches < 0) {
                        processReadWhileParsing(tempread);
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

    cout << size() << " reads processed in total" << endl;

    myfile.close();
    delete [] line;
}

void NGSReadSet::processReadWhileParsing(NGSRead &tempread) {

    //if (!tempread.flag) return;
    int i, id;

    if (!tempread.direction) {
        reverseComplement(tempread.scaff);
        reverseComplement(tempread.read);
    }
    tempread.convertStateStr(tempread.scaff, SEQ_DNA);
    tempread.convertStateStr(tempread.read, SEQ_DNA);
    ASSERT(tempread.scaff.length() == tempread.read.length());

    int nstates = 4 + (!ngs_ignore_gaps);

    for (i = 0, id = 0; i < tempread.scaff.length(); i++) {
        int state1 = tempread.scaff[i];
        int state2 = tempread.read[i];
        if (state1 >= nstates || state2 >= nstates) continue;
        double *pair_pos, *state_pos;
        while (id >= state_freq.size()) {
            state_pos = new double[nstates];
            memset(state_pos, 0, sizeof(double)*(nstates));
            state_freq.push_back(state_pos);
        }
        state_pos = state_freq[id];
        state_pos[state2] += 1.0/tempread.times;
        while (id >= pair_freq.size()) {
            pair_pos = new double[(nstates) * (nstates)];
            memset(pair_pos, 0, sizeof(double)*(nstates) * (nstates));
            pair_freq.push_back(pair_pos);
        }
        pair_pos = pair_freq[id];
        pair_pos[state1*(nstates) + state2] += 1.0/tempread.times;
        id++;
    }

    if (tree) {
        ReadInfo read_info;
        tempread.homo_rate = homo_rate;
        tempread.computePairFreq();
        read_info.homo_distance = tempread.optimizeDist(1.0-tempread.identity);
        read_info.homo_logl = -tempread.computeFunction(read_info.homo_distance);
        tempread.homo_rate = 0.0;
        read_info.distance = tempread.optimizeDist(read_info.homo_distance);
        read_info.logl = -tempread.computeFunction(read_info.distance);
        read_info.id = tempread.id;
        read_info.identity = tempread.identity;
        push_back(read_info);
    }


}

void NGSReadSet::writeInfo() {
    //cout << size() << " reads process in total" << endl;
    return;
}

void NGSReadSet::writeFreqMatrix(ostream &out) {
    int num_states = 4 + (!ngs_ignore_gaps);
    out << pair_freq.size() << " " << num_states << endl;
    vector<double*>::iterator it;
    vector<double*>::iterator pit;

    for (it = pair_freq.begin(), pit = state_freq.begin(); it != pair_freq.end(); it++, pit++) {
        for (int i = 0; i < num_states; i++) {
            for (int j = 0; j < num_states; j++) {
                if (!ngs_ignore_gaps && i == num_states-1 && j == num_states-1)
                    out << int(round((*pit)[i]*((*pit)[i]-1)/2));
                else out << int(round((*it)[i*num_states+j])) << ((j<num_states-1) ? "\t" : "");
            }
            out << endl;
        }
        out << endl;
    }
}

/****************************************************************************
        main function
 ****************************************************************************/

void reportNGSAnalysis(const char *file_name, Params &params, NGSAlignment &aln, NGSTree &tree,
                       DoubleMatrix &rate_info, StrVector &rate_name) {
    ofstream out(file_name);
    out.setf(ios::fixed,ios::floatfield);

    int i, j, k;


    double *rate_param = new double[aln.num_states * aln.num_states];
    double *rate_matrix = new double[aln.num_states * aln.num_states];

    out << "Input file: " << params.ngs_file << endl;
    out << "Model of evolution: " << tree.getModel()->name << endl << endl;

    out << "Substitution process assuming one homogeneous model among all positions:" << endl;

    out << "Rate parameters: " << endl;

    tree.getModel()->getRateMatrix(rate_param);

    /*
     * This isn't a good way of doing it. Rather somewhere high up in the model heirarchy
     * define a bool isSymmetric() method, which is true for time reversible models and
     * not true for nonTR models. ! ModelGTR has 'half_matrix' member, which should do the job.
     */
    if (ModelMarkov::validModelName(tree.getModel()->name)) {
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
            if (j != i) out << rate_matrix[i*aln.num_states+j];
            else out << "-";
        }
        out << endl;
    }
    out << endl;
    out << "State frequencies: ";
    switch (tree.getModel()->getFreqType()) {
    case FREQ_EMPIRICAL:
        out << "(empirical counts from alignment)" << endl;
        break;
    case FREQ_ESTIMATE:
        out << "(estimated with maximum likelihood)" << endl;
        break;
    case FREQ_USER_DEFINED:
        out << "(user-defined)" << endl;
        break;
    case FREQ_EQUAL:
        out << "(equal frequencies)" << endl;
        break;
    default:
        break;
    }

    double *state_freq = new double[aln.num_states];
    tree.getModel()->getStateFrequency(state_freq);

    for (i = 0; i < aln.num_states; i++) out << state_freq[i] << " \t";
    out << endl << endl;

    out << "Q matrix can be obtained by multiplying rate parameters with state frequencies" << endl << endl;

    double *q_mat = new double[tree.aln->num_states * tree.aln->num_states];
    tree.getModel()->getQMatrix(q_mat);

    for (i = 0, k = 0; i < tree.aln->num_states; i++) {
        for (j = 0; j < tree.aln->num_states; j++, k++)
            out << "  " << q_mat[k];
        out << endl;
    }

    delete [] q_mat;

    out << endl;

    out << "Log-likelihood: " << tree.computeLikelihood() << endl << endl;

    out << "Inferred posisiton-specific rates under one model or position-specific model: " << endl;

    out << "Position\tSeq_error";
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
    delete [] state_freq;
    delete [] rate_matrix;
    delete [] rate_param;
}

/*
bool checkFreq(int *pair_freq, int n) {
	int i, count = 0;
	for (i=0; i < n*n; i++)
		if (pair_freq[i] != 0) count++;
	if (count <= n) return false;
	return true;
}*/

void testSingleRateModel(Params &params, NGSAlignment &aln, NGSTree &tree, string model,
                         double *freq, DoubleVector &rate_info, StrVector &rate_name,
                         bool write_info, const char *report_file)
{
    char model_name[20];
    NGSAlignment sum_aln(aln.num_states, 1, freq);
    ModelsBlock *models_block = new ModelsBlock;

    NGSTree sum_tree(params, &sum_aln);
    sum_aln.tree = &sum_tree;

    if (model == "")
        sprintf(model_name, "GTR+F1");
    else
        sprintf(model_name, "%s+F1", model.c_str());
    try {
        params.model_name = model_name;
        sum_tree.setModelFactory(new ModelFactory(params, &sum_tree, models_block));
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

    double *rate_mat = new double[aln.num_states*aln.num_states];
    memset(rate_mat, 0, aln.num_states*aln.num_states*sizeof(double));
    sum_tree.getModel()->getRateMatrix(rate_mat);
    rate_info.insert(rate_info.end(), rate_mat, rate_mat+sum_tree.getModel()->getNumRateEntries());

    if (tree.getModel()->isReversible()) {
        sum_tree.getModel()->getStateFrequency(rate_mat);
        rate_info.insert(rate_info.end(), rate_mat, rate_mat+aln.num_states);
    }
	delete [] rate_mat;
	delete models_block;

    if (report_file) {
        DoubleMatrix tmp(1);
        tmp[0] = rate_info;
        reportNGSAnalysis(report_file, params, sum_aln, sum_tree, tmp, rate_name);
    }
}

void testTwoRateModel(Params &params, NGSAlignment &aln, NGSTree &tree, string model,
                      double *freq, DoubleVector &rate_info, StrVector &rate_name,
                      bool write_info, const char *report_file)
{
    char model_name[20];
    NGSAlignment sum_aln(aln.num_states, 1, freq);


    NGSTreeCat sum_tree(params, &sum_aln);
    sum_aln.tree = &sum_tree;

    ModelsBlock *models_block = new ModelsBlock;

    if (model == "")
        sprintf(model_name, "GTR+FC2");
    else
        sprintf(model_name, "%s+FC2", model.c_str());
    try {
        params.model_name = model_name;
        sum_tree.setModelFactory(new ModelFactory(params, &sum_tree, models_block));
        sum_tree.setModel(sum_tree.getModelFactory()->model);
        sum_tree.setRate(sum_tree.getModelFactory()->site_rate);
        double bestTreeScore = sum_tree.getModelFactory()->optimizeParameters(false, write_info);
        cout << "LogL: " << bestTreeScore;
        cout << " / Rate: " << sum_tree.getRate()->getRate(0) << endl;
    } catch (const char*) {
        cout << "Skipped due to sparse matrix" << endl;
        //rate_info.insert(rate_info.end(), rate_name.size(), MIN_SITE_RATE);
        return;
    } catch (string &str) {
        cout << str;
        return;
    }
    delete models_block;
    //return sum_tree.getRate()->getRate(0);

    /*
    	rate_info.push_back(sum_tree.getRate()->getRate(0));

        double rate_mat[aln.num_states*aln.num_states];
        memset(rate_mat, 0, aln.num_states*aln.num_states*sizeof(double));
        sum_tree.getModel()->getRateMatrix(rate_mat);
        rate_info.insert(rate_info.end(), rate_mat, rate_mat+sum_tree.getModel()->getNumRateEntries());

    	if (tree.getModel()->isReversible()) {
    		sum_tree.getModel()->getStateFrequency(rate_mat);
    		rate_info.insert(rate_info.end(), rate_mat, rate_mat+aln.num_states);
        }

    	if (report_file) {
    		DoubleMatrix tmp(1);
    		tmp[0] = rate_info;
    		reportNGSAnalysis(report_file, params, sum_aln, sum_tree, tmp, rate_name);
    	}*/
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
    cout << endl << "Read distances to the reference written to: " << file_name << endl << endl;

    string count_file = params.ngs_mapped_reads;
    count_file += ".freq";
    out.open(count_file.c_str());
    ngs_reads.writeFreqMatrix(out);
    out.close();
    cout << "Position-specific pair counts written to: " << count_file << endl << endl;

}

void computePairCount(Params &params, NGSTree *tree, double homo_rate) {
    NGSReadSet ngs_reads;
    ngs_reads.tree = tree;
    ngs_reads.homo_rate = homo_rate;
    ngs_reads.ngs_ignore_gaps = params.ngs_ignore_gaps;
    //cout << "Homogeneous rate: " << ngs_reads.homo_rate << endl;
    cout << "Computing read distances to reference from file " << params.ngs_mapped_reads << " ... " << endl;
    ngs_reads.parseNextGen(params.ngs_mapped_reads);
    ngs_reads.writeInfo();

    string out_file = params.ngs_mapped_reads;
    out_file += ".dist";
    reportNGSReads(out_file.c_str(), params, ngs_reads);
}



void runNGSAnalysis(Params &params) {

    time_t begin_time;
    time(&begin_time);

    char model_name[20];

    if (!params.ngs_file) {
        computePairCount(params, NULL, 0.0);
        return;
    }

    // read input file, initialize NGSAlignment
    NGSAlignment aln(params.ngs_file);
    cout.setf(ios::fixed,ios::floatfield);

    //params.freq_type = FREQ_ESTIMATE;

    // initialize NGSTree
    NGSTree tree(params, &aln);
    aln.tree = &tree;
    ModelsBlock *models_block = new ModelsBlock;

    // initialize Model
    string original_model = params.model_name;
    if (params.model_name == "") {
        sprintf(model_name, "GTR+F%d", aln.ncategory);
        params.freq_type = FREQ_ESTIMATE;
    }
    else
        sprintf(model_name, "%s+F%d", params.model_name.c_str(), aln.ncategory);
    params.model_name = model_name;
    tree.setModelFactory(new ModelFactory(params, &tree, models_block));
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);

    delete models_block;

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

    rate_name.push_back("Hete_error");

    if (tree.getModel()->isReversible()) {
        for (i = 0; i < aln.num_states-1; i++)
            for (j = i+1; j < aln.num_states; j++) {
                stringstream x;
                x << aln.convertStateBackStr(i) << "<->" << aln.convertStateBackStr(j);
                rate_name.push_back(x.str());
            }
        for (i = 0; i < aln.num_states; i++) {
            stringstream x;
            x << aln.convertStateBackStr(i);
            rate_name.push_back(x.str());
        }
    } else {
        for (i = 0; i < aln.num_states; i++)
            for (j = 0; j < aln.num_states; j++) if (j != i) {
                    stringstream x;
                    x << aln.convertStateBackStr(i) << "->" << aln.convertStateBackStr(j);
                    rate_name.push_back(x.str());
                }
    }


    VerboseMode vb_saved = verbose_mode;
    verbose_mode = VB_QUIET;

    cout << endl << "--> INFERING RATE ASSUMING POSITION-SPECIFIC MODEL..." << endl << endl;
    for (int pos = 0; pos < aln.ncategory; pos++) {
        cout << "Position " << pos+1 << " / ";
        double *pair_pos = aln.pair_freq + (pos*aln.num_states*aln.num_states);
        testSingleRateModel(params, aln, tree, original_model, pair_pos, part_rate[pos], rate_name, false, NULL);
    }


    verbose_mode = vb_saved;

    double *sum_freq = new double[aln.num_states*aln.num_states];
    cout << endl << "-->INFERING RATE UNDER EQUAL-RATE NULL MODEL..." << endl << endl;
    aln.computeSumPairFreq(sum_freq);
    DoubleVector null_rate;
    string out_file = params.out_prefix;
    out_file += ".ngs_e";
    for (i = 0; i < aln.num_states*aln.num_states; i++)
        cout << sum_freq[i] << " ";
    cout << endl;
    testSingleRateModel(params, aln, tree, original_model, sum_freq, null_rate, rate_name, true, out_file.c_str());

    DoubleVector two_rate;

    cout << endl << "-->INFERING RATE UNDER TWO-RATE MODEL..." << endl << endl;
    testTwoRateModel(params, aln, tree, original_model, sum_freq, two_rate, rate_name, true, NULL);


    // report running results
    out_file = params.out_prefix;
    out_file += ".ngs";
    reportNGSAnalysis(out_file.c_str(), params, aln, tree, part_rate, rate_name);

    if (params.ngs_mapped_reads) {
        computePairCount(params, &tree, null_rate[0]);
    }


    time_t end_time;
    time(&end_time);

    cout << "Total run time: " << difftime(end_time, begin_time) << " seconds" << endl << endl;
	delete [] sum_freq;
}
