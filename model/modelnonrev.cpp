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
#include "modelnonrev.h"
#include "modelunrest.h"
#include "modelliemarkov.h"
//#include "whtest/eigen.h"

ModelNonRev::ModelNonRev(PhyloTree *tree)
        : ModelGTR(tree, false)
{
	// model_parameters must be initialized by subclass
//    ModelGTR::freeMem();

	fixed_parameters = false;
    model_parameters = NULL;
    ignore_state_freq = true;

    int num_rates = getNumRateEntries();
    
    // reallocate the mem spaces
    delete [] rates;
    rates = new double [num_rates];
    memset(rates, 0, sizeof(double) * (num_rates));

    rate_matrix = aligned_alloc<double>(num_states*num_states);
    temp_space =  aligned_alloc<double>(num_states*num_states);
    eigenvalues_imag = aligned_alloc<double>(num_states);
    
    ceval = aligned_alloc<complex<double> >(num_states);
    cevec = aligned_alloc<complex<double> >(num_states*num_states);
    cinv_evec = aligned_alloc<complex<double> >(num_states*num_states);
    
    if (!tree->rooted) {
        cout << "Converting unrooted to rooted tree..." << endl;
        tree->convertToRooted();
    }
}

ModelNonRev::~ModelNonRev() {
    // this is not necessary and causes double free memory, because freeMem is a virtual function
    // thus will be called by the destructor of ModelGTR anyway.
//	freeMem();
}

void ModelNonRev::freeMem() {
    ModelGTR::freeMem();
    aligned_free(cinv_evec);
    aligned_free(cevec);
    aligned_free(ceval);
    aligned_free(eigenvalues_imag);
    aligned_free(temp_space);
    aligned_free(rate_matrix);
    delete [] model_parameters;
}

/* static */ ModelNonRev* ModelNonRev::getModelByName(string model_name, PhyloTree *tree, string model_params, bool count_rates) {
	if (ModelUnrest::validModelName(model_name)) {
		return((ModelNonRev*)new ModelUnrest(tree, model_params, count_rates));
	} else if (ModelLieMarkov::validModelName(model_name)) {
	        return((ModelNonRev*)new ModelLieMarkov(model_name, tree, model_params, count_rates));
	} else {
		cerr << "Unrecognized model name " << model_name << endl;
		return((ModelNonRev*)NULL);
	}
}

/* static */ bool ModelNonRev::validModelName(string model_name) {
	return ModelUnrest::validModelName(model_name) 
	  || ModelLieMarkov::validModelName(model_name);
}

void ModelNonRev::getQMatrix(double *rate_mat) {
    memmove(rate_mat, rate_matrix, num_states*num_states*sizeof(double));
}

/* BQM: Ziheng Yang code which fixed old matinv function */
int matinv (double x[], int n, int m, double space[])
{
    /* x[n*m]  ... m>=n
       space[n].  This puts the fabs(|x|) into space[0].  Check and calculate |x|.
       Det may have the wrong sign.  Check and fix.
    */
    int i,j,k;
    int *irow=(int*) space;
    double ee=1e-100, t,t1,xmax, det=1;

    for (i=0; i<n; i++) irow[i]=i;

    for (i=0; i<n; i++)  {
        xmax = fabs(x[i*m+i]);
        for (j=i+1; j<n; j++)
            if (xmax<fabs(x[j*m+i]))
            {
                xmax = fabs(x[j*m+i]);
                irow[i]=j;
            }
        det *= x[irow[i]*m+i];
        if (xmax < ee)   {
            cout << endl << "xmax = " << xmax << " close to zero at " << i+1 << "!\t" << endl;
            exit(-1);
        }
        if (irow[i] != i) {
            for (j=0; j < m; j++) {
                t = x[i*m+j];
                x[i*m+j] = x[irow[i]*m+j];
                x[irow[i]*m+j] = t;
            }
        }
        t = 1./x[i*m+i];
        for (j=0; j < n; j++) {
            if (j == i) continue;
            t1 = t*x[j*m+i];
            for (k=0; k<m; k++)  x[j*m+k] -= t1*x[i*m+k];
            x[j*m+i] = -t1;
        }
        for (j=0; j < m; j++)   x[i*m+j] *= t;
        x[i*m+i] = t;
    }                            /* for(i) */
    for (i=n-1; i>=0; i--) {
        if (irow[i] == i) continue;
        for (j=0; j < n; j++)  {
            t = x[j*m+i];
            x[j*m+i] = x[j*m + irow[i]];
            x[j*m + irow[i]] = t;
        }
    }
    space[0]=det;
    return(0);
}

int computeStateFreqFromQMatrix (double Q[], double pi[], int n, double space[])
{
    /* from rate matrix Q[] to pi, the stationary frequencies:
       Q' * pi = 0     pi * 1 = 1
       space[] is of size n*(n+1).
    */
    int i,j;
    double *T = space;      /* T[n*(n+1)]  */

    for (i=0;i<n+1;i++) T[i]=1;
    for (i=1;i<n;i++) {
        for (j=0;j<n;j++)
            T[i*(n+1)+j] =  Q[j*n+i];     /* transpose */
        T[i*(n+1)+n] = 0.;
    }
    matinv(T, n, n+1, pi);
    for (i=0;i<n;i++)
        pi[i] = T[i*(n+1)+n];
    return (0);
}
/* End of Ziheng Yang code */

void ModelNonRev::decomposeRateMatrix() {
    int i, j, k;
    double sum;
    //double m[num_states];
    double *space = new double[num_states*(num_states+1)];

    for (i = 0; i < num_states; i++)
        state_freq[i] = 1.0/num_states;

    for (i = 0, k = 0; i < num_states; i++) {
        rate_matrix[i*num_states+i] = 0.0;
        double row_sum = 0.0;
        for (j = 0; j < num_states; j++)
            if (j != i) {
                row_sum += (rate_matrix[i*num_states+j] = rates[k++]);
            }
        rate_matrix[i*num_states+i] = -row_sum;
    }
    computeStateFreqFromQMatrix(rate_matrix, state_freq, num_states, space);


    for (i = 0, sum = 0.0; i < num_states; i++) {
        sum -= rate_matrix[i*num_states+i] * state_freq[i]; /* exp. rate */
    }

    if (sum == 0.0) throw "Empty Q matrix";

    double delta = total_num_subst / sum; /* 0.01 subst. per unit time */

    for (i = 0; i < num_states; i++) {
        for (j = 0; j < num_states; j++) {
            rate_matrix[i*num_states+j] *= delta;
        }
    }
    delete [] space;

    if (phylo_tree->params->matrix_exp_technique == MET_EIGEN_DECOMPOSITION) {
        eigensystem_nonrev(rate_matrix, state_freq, eigenvalues, eigenvalues_imag, eigenvectors, inv_eigenvectors, num_states);
    }
}


void ModelNonRev::writeInfo(ostream &out) {
	int i, j, k;
    out << "Model parameters: ";
    if (num_params>0) out << model_parameters[0];
    for (i=1; i < num_params; i++) out << "," << model_parameters[i];
    out << endl;

    if (num_states != 4) return;
    out << "Rate parameters:" << endl;
    for (i = 0, k = 0; i < num_states; i++) {
        switch (i) {
        case 0:
            out << "A";
            break;
        case 1:
            out << "C";
            break;
        case 2:
            out << "G";
            break;
        case 3:
            out << "T";
            break;
        }
        for (j = 0; j < num_states; j++)
            if (j != i)
                out << '\t' << rates[k++];
            else out << '\t' << "-";
        out << endl;
    }
}


int matby (double a[], double b[], double c[], int n,int m,int k)
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b
*/
{
    int i,j,i1;
    double t;
    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++) {
            for (i1=0,t=0; i1<m; i1++) t+=a[i*m+i1]*b[i1*k+j];
            c[i*k+j] = t;
        }
    return (0);
}

int matexp (double Q[], double t, int n, int TimeSquare, double space[])
{
    /* This calculates the matrix exponential P(t) = exp(t*Q).
       Input: Q[] has the rate matrix, and t is the time or branch length.
              TimeSquare is the number of times the matrix is squared and should
              be from 5 to 31.
       Output: Q[] has the transition probability matrix, that is P(Qt).
       space[n*n]: required working space.

          P(t) = (I + Qt/m + (Qt/m)^2/2)^m, with m = 2^TimeSquare.

       T[it=0] is the current matrix, and T[it=1] is the squared result matrix,
       used to avoid copying matrices.
       Use an even TimeSquare to avoid one round of matrix copying.
    */
    int it, i;
    double *T[2];

    if (TimeSquare<2 || TimeSquare>31) cout << "TimeSquare not good" << endl;
    T[0]=Q;
    T[1]=space;
    for (i=0; i<n*n; i++)  T[0][i] = ldexp( Q[i]*t, -TimeSquare );

    matby (T[0], T[0], T[1], n, n, n);
    for (i=0; i<n*n; i++)  T[0][i] += T[1][i]/2;
    for (i=0; i<n; i++)  T[0][i*n+i] ++;

    for (i=0,it=0; i<TimeSquare; i++) {
        it = !it;
        matby (T[1-it], T[1-it], T[it], n, n, n);
    }
    if (it==1)
        for (i=0;i<n*n;i++) Q[i]=T[1][i];
    return(0);
}

const int TimeSquare = 10;


void ModelNonRev::computeTransMatrixEigen(double time, double *trans_matrix) {
	/* compute P(t) */
	double evol_time = time / total_num_subst;
    int nstates_2 = num_states*num_states;
	double *exptime = new double[nstates_2];
	int i, j, k;

    memset(exptime, 0, sizeof(double)*nstates_2);
	for (i = 0; i < num_states; i++)
        if (eigenvalues_imag[i] == 0.0) {
            exptime[i*num_states+i] = exp(evol_time * eigenvalues[i]);
        } else {
            assert(i < num_states-1 && eigenvalues_imag[i+1] != 0.0 && eigenvalues_imag[i] > 0.0);
            complex<double> exp_eval(eigenvalues[i] * evol_time, eigenvalues_imag[i] * evol_time);
            exp_eval = exp(exp_eval);
            exptime[i*num_states+i] = exp_eval.real();
            exptime[i*num_states+i+1] = exp_eval.imag();
            i++;
            exptime[i*num_states+i] = exp_eval.real();
            exptime[i*num_states+i-1] = -exp_eval.imag();
        }


    // compute V * exp(L t)
    for (i = 0; i < num_states; i++)
        for (j = 0; j < num_states; j++) {
            double val = 0;
            for (k = 0; k < num_states; k++)
                val += eigenvectors[i*num_states+k] * exptime[k*num_states+j];
            trans_matrix[i*num_states+j] = val;
        }

    memcpy(exptime, trans_matrix, sizeof(double)*nstates_2);

    // then compute V * exp(L t) * V^{-1}
    for (i = 0; i < num_states; i++) {
        double row_sum = 0.0;
        for (j = 0; j < num_states; j++) {
            double val = 0;
            for (k = 0; k < num_states; k++)
                val += exptime[i*num_states+k] * inv_eigenvectors[k*num_states+j];
            // make sure that trans_matrix are non-negative
            assert(val >= -0.001);
            val = fabs(val);
            trans_matrix[i*num_states+j] = val;
            row_sum += val;
        }
        assert(fabs(row_sum-1.0) < 1e-4);
    }

    delete [] exptime;
}

void ModelNonRev::computeTransMatrix(double time, double *trans_matrix) {
    if (phylo_tree->params->matrix_exp_technique == MET_EIGEN_DECOMPOSITION) {
        computeTransMatrixEigen(time, trans_matrix);
    } else {
        // scaling and squaring technique
        int statesqr = num_states*num_states;
        memcpy(trans_matrix, rate_matrix, statesqr*sizeof(double));
        matexp(trans_matrix, time, num_states, TimeSquare, temp_space);
    }
    // 2016-04-05: 2nd version
//    for (int i = 0; i < statesqr; i++) 
//        trans_matrix[i] *= time;
//    double space[NCODE*NCODE*3] = {0};
//    matexp2(trans_matrix, num_states, 7, 5, space);
}

double ModelNonRev::computeTrans(double time, int state1, int state2) {
    double *trans_matrix = new double[num_states*num_states];
//    memcpy(trans_matrix, rate_matrix, num_states*num_states*sizeof(double));
//    matexp(trans_matrix, time, num_states, TimeSquare, temp_space);
    computeTransMatrix(time, trans_matrix);
    double trans = trans_matrix[state1*num_states+state2];
    delete [] trans_matrix;
    return trans;
}

void ModelNonRev::computeTransDerv(double time, double *trans_matrix, double *trans_derv1, double *trans_derv2) {
    computeTransMatrix(time, trans_matrix);
    int i, j, k;
    // First derivative = Q * e^(Qt)
    for (i = 0; i < num_states; i++)
        for (j = 0; j < num_states; j++) {
            double val = 0.0;
            for (k = 0; k < num_states; k++)
                val += rate_matrix[i*num_states+k] * trans_matrix[k*num_states+j];
            trans_derv1[i*num_states+j] = val;
        }
        
    // Second derivative = Q * Q * e^(Qt)
    for (i = 0; i < num_states; i++)
        for (j = 0; j < num_states; j++) {
            double val = 0.0;
            for (k = 0; k < num_states; k++)
                val += rate_matrix[i*num_states+k] * trans_derv1[k*num_states+j];
            trans_derv2[i*num_states+j] = val;
        }
}

int ModelNonRev::getNDim() {
	if (fixed_parameters)
		return 0;
	return(num_params);
}

void ModelNonRev::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	// I don't know the proper C++ way to handle this: got error if I didn't define something here.
	cerr << "setBounds should only be called on subclass of ModelNonRev\n";
}

void ModelNonRev::setVariables(double *variables) {
	int nrate = getNDim();
	if (nrate > 0)
		memcpy(variables+1, model_parameters, nrate*sizeof(double));
}

bool ModelNonRev::getVariables(double *variables) {
	int nrate = getNDim();
	int i;
	bool changed = false;
	for (i = 0; i < nrate && !changed; i++) changed = (model_parameters[i] != variables[i+1]);
	if (changed) {
		memcpy(model_parameters, variables+1, nrate * sizeof(double));
		this->setRates();
	}
	return changed;
}

void ModelNonRev::setRates() {
	// I don't know the proper C++ way to handle this: got error if I didn't define something here.
	cerr << "setRates should only be called on subclass of ModelNonRev\n";
}

double ModelNonRev::targetFunk(double x[]) {
	bool changed = getVariables(x);
//	if (state_freq[num_states-1] < 1e-4) return 1.0e+12;
	if (changed) {
		decomposeRateMatrix();
		assert(phylo_tree);
		phylo_tree->clearAllPartialLH();
	}
    for (int i = 0; i < num_states; i++)
        if (state_freq[i] < 1e-4)
            return 1e12;
	return -phylo_tree->computeLikelihood();
}

double ModelNonRev::optimizeParameters(double gradient_epsilon) {
	int ndim = getNDim();
	
	// return if nothing to be optimized
	if (ndim == 0) return 0.0;

	if (verbose_mode >= VB_MAX)
		cout << "Optimizing " << name << " model parameters..." << endl;

	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	double score;

	// by BFGS algorithm
	setVariables(variables);
	setBounds(lower_bound, upper_bound, bound_check);
//    if (phylo_tree->params->optimize_alg.find("BFGS-B") == string::npos)
        score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_RATE));
//    else
//        score = -L_BFGS_B(ndim, variables+1, lower_bound+1, upper_bound+1, max(gradient_epsilon, TOL_RATE));

	bool changed = getVariables(variables);
    if (changed) {
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
        score = phylo_tree->computeLikelihood();
    }
	
	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}


void ModelNonRev::saveCheckpoint() {
    checkpoint->startStruct("ModelNonRev");
    CKP_ARRAY_SAVE(num_params, model_parameters);
    checkpoint->endStruct();
    ModelSubst::saveCheckpoint();
}

void ModelNonRev::restoreCheckpoint() {
    ModelSubst::restoreCheckpoint();
    checkpoint->startStruct("ModelNonRev");
    CKP_ARRAY_RESTORE(num_params, model_parameters);
    checkpoint->endStruct();
    this->setRates();
}
