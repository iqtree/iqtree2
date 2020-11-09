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

#include "tree/phylotree.h"
#include "ratemeyerdiscrete.h"
//#include "kmeans/KMeans.h"
//#include "modeltest_wrapper.h"

/************************************************
	Huy's k-means dynamic programming algorithm
************************************************/

void quicksort(double arr[], int weight[], int index[], int left, int right) {
      int i = left, j = right, tmp2;
      double tmp;
      double pivot = arr[(left + right) / 2];
 
      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  tmp2 = index[i];
                  index[i] = index[j];
                  index[j] = tmp2;
                  tmp2 = weight[i];
                  weight[i] = weight[j];
                  weight[j] = tmp2;
                  i++;
                  j--;
            }
      };
 
      /* recursion */
      if (left < j)
            quicksort(arr, weight, index, left, j);
      if (i < right)
            quicksort(arr, weight, index, i, right);
}

double mean_sum(int l, int r, double *sumA, double *sumAsquare, int *sumW) {
/*	double mean = (sumA[r]-sumA[l-1])/(r-l+1);
	return sumAsquare[r]-sumAsquare[l-1]- 2*(sumA[r]- sumA[l-1])*mean + mean*mean*(r-l+1);*/

	double sum = (sumA[r]- sumA[l-1]);
	return sumAsquare[r]-sumAsquare[l-1] - sum*sum/(sumW[r] - sumW[l-1]);

/*
	double mean = (sumA[r]-sumA[l-1]);
	return sumAsquare[r]-sumAsquare[l-1]- 2*(sumA[r]- sumA[l-1])*mean + mean*mean*(r-l+1);*/
}



// Runs k-means on the given set of points.
//   - n: The number of points in the data set
//   - k: The number of clusters to look for
//   - d: The number of dimensions that the data set lives in
//   - points: An array of size n*d where points[d*i + j] gives coordinate j of poi
//   - attempts: The number of times to independently run k-means with different starting centers.
//               The best result is always returned (as measured by the cost function).
//   - centers: This can either be null or an array of size k*d. In the latter case, it will be
//              filled with the locations of all final cluster centers. Specifically
//              centers[d*i + j] will give coordinate j of center i. If the cluster is unused, it
//              will contain NaN instead.
//   - assignments: This can either be null or an array of size n. In the latter case, it will be
//                  filled with the cluster that each pois assigned to (an integer between 0
//                  and k-1 inclusive).
// The final cost of the clustering is also returned.

double RunKMeans1D(int n, int k, double *points_orig, int *weights, double *centers, int *assignments) {
	double *sumA;
	double *sumAsquare;
	int *sumW;
	double **Cost;
	int **trace;
	
	sumA = new double[n+1];
	sumAsquare = new double[n+1];
	sumW = new int[n+1];
	Cost = new double*[n+1];
	for (int i=0; i<=n; i++) Cost[i] = new double[k+1];
	trace = new int*[n+1];
	for (int i=0; i<=n; i++) trace[i] = new int[k+1];

	int *index = new int[n+1];
	for (int i=0; i<n; i++) index[i] = i;
    double *points = new double[n];
    memcpy(points, points_orig, n*sizeof(double));

	//for (int i=1; i<=n; i++) cout <<index[i] <<"\t" <<points[i] <<endl;

	quicksort(points, weights, index, 0, n-1);
	
	//for (int i=n; i>0; i--) {points[i] = points[i-1]; index[i] = index[i-1];}
	//for (int i=1; i<=n; i++) cout <<index[i] <<"\t" <<points[i] <<endl;
	
	//exit(1);
	
	sumA[0] = 0; sumAsquare[0] =0; sumW[0] = 0;
	for (int i=1; i<=n; i++) {
		/*sumA[i] = sumA[i-1] + points[i-1];
		sumAsquare[i] = sumAsquare[i-1] + points[i-1]*points[i-1];*/
		sumA[i] = sumA[i-1] + points[i-1] * weights[i-1];
		sumAsquare[i] = sumAsquare[i-1] + points[i-1]*points[i-1] * weights[i-1];
		sumW[i] = sumW[i-1] + weights[i-1];
	}
	
	Cost[0][0] = 0;
	for (int i=1; i<=n; i++) {
		Cost[i][1] = mean_sum(1, i, sumA, sumAsquare, sumW);
		trace[i][1] = 0;
		for (int j=2; j<=(i<=k?i:k); j++) {
			Cost[i][j] = Cost[j-1][j-1]+ mean_sum(j, i, sumA, sumAsquare, sumW);
			trace[i][j] = j-1;
			for (int _k=j; _k<=i-1; _k++) {
				double temp = mean_sum(_k+1, i, sumA, sumAsquare, sumW);
				if (Cost[i][j] >= Cost[_k][j-1]+ temp) {
					Cost[i][j] = Cost[_k][j-1]+ temp;			
					trace[i][j] = _k;
				}
			}
		}
	}
	
	double min_cost = Cost[n][k];
	
    for (int i=n, j=k; i>0; --j) {
        int t= trace[i][j];
        centers[j-1] = (sumA[i]-sumA[t])/(sumW[i]-sumW[t]);
        //cout << "category " <<k-j<<endl;
        for (int _i=t+1; _i<=i; _i++) {
            //cout <<index[_i] << "\t" <<points[_i] <<endl;
            assignments[index[_i-1]] = j-1; //points[_i] \in category k-j
        }
        i=t;
    }
	
    delete [] points;
	for (int i=n; i>=0; i--) delete [] trace[i];
	delete [] trace;
	for (int i=n; i>=0; i--) delete [] Cost[i];
	delete [] Cost;

	delete [] sumW;	
	delete [] sumAsquare;
	delete [] sumA;

	return min_cost;
}
     

/************************************************
	RateMeyerDiscrete
************************************************/
RateMeyerDiscrete::RateMeyerDiscrete(int ncat, int cat_type, char *file_name, PhyloTree *tree, bool rate_type)
 : RateMeyerHaeseler(file_name, tree, rate_type)
{
	ncategory = ncat;
	rates = NULL;
	ptn_cat = NULL;
	is_categorized = false;
	mcat_type = cat_type;
	if (ncat > 0) {
		rates = new double[ncategory];
		memset(rates, 0, sizeof(double) * ncategory);
	}
	name += convertIntToString(ncategory);
	if (ncategory > 0)
		full_name += " with " + convertIntToString(ncategory) + " categories";
	else
		full_name += "auto-detect #categories";
}

RateMeyerDiscrete::RateMeyerDiscrete() {
	ncategory = 0;
	rates = NULL;
	ptn_cat = NULL;
	is_categorized = false;
	mcat_type = 0;
	rates = NULL;
	name = full_name = "";
	rate_mh = true;
}


RateMeyerDiscrete::~RateMeyerDiscrete()
{
	delete [] rates;
}

bool RateMeyerDiscrete::isSiteSpecificRate() { 
	return !is_categorized; 
}

int RateMeyerDiscrete::getNDiscreteRate() { 
	if (!is_categorized) return RateMeyerHaeseler::getNDiscreteRate();
	ASSERT(ncategory > 0);
	return ncategory; 
}

double RateMeyerDiscrete::getRate(int category) {
	if (!is_categorized) return RateMeyerHaeseler::getRate(category);
	ASSERT(category < ncategory); 
	return rates[category]; 
}

int RateMeyerDiscrete::getPtnCat(int ptn) {
	if (!is_categorized) return RateMeyerHaeseler::getPtnCat(ptn);
	ASSERT(ptn_cat);
	return ptn_cat[ptn];
}

double RateMeyerDiscrete::getPtnRate(int ptn) {
	if (!is_categorized) return RateMeyerHaeseler::getPtnRate(ptn);
	ASSERT(ptn_cat && rates);
	return rates[ptn_cat[ptn]];
}

int RateMeyerDiscrete::computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat) {
	pattern_rates.insert(pattern_rates.begin(), begin(), end());
	pattern_cat.insert(pattern_cat.begin(), ptn_cat, ptn_cat + size());
    return ncategory;
}

/*double RateMeyerDiscrete::optimizeParameters() {
	if (is_categorized) {
		is_categorized = false;
		phylo_tree->clearAllPartialLh();
		return phylo_tree->computeLikelihood();
	}
	double tree_lh = RateMeyerHaeseler::optimizeParameters();
	return tree_lh;
}*/

double RateMeyerDiscrete::optimizeParameters(double epsilon) {
	if (!is_categorized) return RateMeyerHaeseler::optimizeParameters(epsilon);
	phylo_tree->calcDist(dist_mat);
	for (int i = 0; i < ncategory; i++)
		optimizeCatRate(i);
	normalizeRates();
	phylo_tree->clearAllPartialLH();
	return phylo_tree->computeLikelihood();
	//return phylo_tree->optimizeAllBranches(2);
}


double RateMeyerDiscrete::computeFunction(double value) {
	if (!is_categorized) return RateMeyerHaeseler::computeFunction(value);
	if (!rate_mh) {
		if (value != cur_scale) {
			ptn_tree->scaleLength(value/cur_scale);
			cur_scale = value;
			ptn_tree->clearAllPartialLH();
		}
		return -ptn_tree->computeLikelihood();
	}

	double lh = 0.0;
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;
	ModelSubst *model = phylo_tree->getModel();
    int trans_size = nstate * nstate;
	double *trans_mat = new double[trans_size];
	int *pair_freq = new int[trans_size];
    
    auto frequencies = phylo_tree->getConvertedSequenceFrequencies();
    for (size_t i = 0; i < nseq-1; i++) {
        auto eyeSequence = phylo_tree->getConvertedSequenceByNumber(i);
        for (size_t j = i + 1; j < nseq; j++) {
            auto jaySequence = phylo_tree->getConvertedSequenceByNumber(j);
            memset(pair_freq, 0, trans_size * sizeof(int));
            if (jaySequence!=nullptr) { 
                for (size_t k = 0; k < size(); k++) {
                    if (ptn_cat[k] != optimizing_cat) {
                        continue;
                    }
                    int state1 = eyeSequence[k];
                    auto pairRow = pair_freq + state1*nstate;
                    if (nstate<=state1) {
                        continue;
                    }
                    int state2 = jaySequence[k];
                    if ( state2 < nstate) {
                        pairRow[state2] += frequencies[k];
                    }
                }
            } else {
                for (size_t k = 0; k < size(); k++) {
                    if (ptn_cat[k] != optimizing_cat) {
                        continue;
                    }
                    Pattern *pat = & phylo_tree->aln->at(k);
                    int state1 = pat->at(i);
                    int state2 = pat->at(j);
                    if ( state1 < nstate && state2 < nstate) {
                        pair_freq[state1*nstate + state2] += pat->frequency;
                    }
                }
            }
            model->computeTransMatrix(value * dist_mat[i*nseq + j], trans_mat);
            for (size_t k = 0; k < trans_size; k++)
            {
                lh -= pair_freq[k] * log(trans_mat[k]);
            }
        }
    }
delete [] pair_freq;
delete [] trans_mat;
return lh;
}

void RateMeyerDiscrete::computeFuncDerv(double value, double &df, double &ddf) {
	if (!is_categorized) {
		RateMeyerHaeseler::computeFuncDerv(value, df, ddf);
		return;
	}

//	double lh = 0.0;
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;
	ModelSubst *model = phylo_tree->getModel();
    int trans_size = nstate * nstate;
	double *trans_mat = new double[trans_size];
	double *trans_derv1 = new double[trans_size];
	double *trans_derv2 = new double[trans_size];
	df = ddf = 0.0;

    int *pair_freq = new int[trans_size];
    auto frequencies = phylo_tree->getConvertedSequenceFrequencies();
    for (size_t i = 0; i + 1 < nseq; ++i) {
        auto eyeSequence = phylo_tree->getConvertedSequenceByNumber(i);
        for (size_t j = i+1; j < nseq; ++j) {
            auto jaySequence = phylo_tree->getConvertedSequenceByNumber(j);
            memset(pair_freq, 0, trans_size * sizeof(int));
            if (frequencies!=nullptr && eyeSequence!=nullptr && jaySequence!=nullptr) {
                for (size_t k = 0; k < size(); ++k) {
                    if (ptn_cat[k] != optimizing_cat) {
                        continue;
                    }
                    int state1 = eyeSequence[k];
                    if (nstate<=state1) {
                        continue;
                    }
                    auto pairRow = pair_freq + state1*nstate;
                    int state2 = jaySequence[k];
                    if (nstate<=state2) {
                        continue;
                    }
                    pairRow[state2] += frequencies[k];
                }
            } else {
                for (size_t k = 0; k < size(); ++k) {
                    if (ptn_cat[k] != optimizing_cat) {
                        continue;
                    }
                    Pattern *pat = & phylo_tree->aln->at(k);
                    int state1 = pat->at(i);
                    int state2 = pat->at(j);
                    if (state1 < nstate && state2 < nstate)
                        pair_freq[state1*nstate + state2] += pat->frequency;
                }
            }
            double dist = dist_mat[i*nseq + j];
            double derv1 = 0.0, derv2 = 0.0;
            model->computeTransDerv(value * dist, trans_mat, trans_derv1, trans_derv2);
            for (size_t k = 0; k < trans_size; ++k) {
                if (pair_freq[k]) {
                    double t1 = trans_derv1[k] / trans_mat[k];
                    double t2 = trans_derv2[k] / trans_mat[k];
                    trans_derv1[k] = t1;
                    trans_derv2[k] = (t2 - t1*t1);
                    //lh -= log(trans_mat[k]) * pair_freq[k];
                    derv1 += trans_derv1[k] * pair_freq[k];
                    derv2 += trans_derv2[k] * pair_freq[k];
                }
            }
            df -= derv1 * dist;
            ddf -= derv2 * dist * dist;
        }
    }
    delete [] pair_freq;
    delete [] trans_derv2;
    delete [] trans_derv1;
    delete [] trans_mat;
}


double RateMeyerDiscrete::optimizeCatRate(int cat) {
	optimizing_cat = cat;
	double negative_lh;
	double current_rate = rates[cat];
	double ferror, optx;

	if (!rate_mh) {
		IntVector ptn_id;
		for (int i = 0; i < size(); i++)
			if (ptn_cat[i] == optimizing_cat)
				ptn_id.push_back(i);
		prepareRateML(ptn_id);
	}

    if (phylo_tree->optimize_by_newton && rate_mh) // Newton-Raphson method 
	{
    	optx = minimizeNewtonSafeMode(MIN_SITE_RATE, current_rate, MAX_SITE_RATE, TOL_SITE_RATE, negative_lh);
    }
    else {
		optx = minimizeOneDimen(MIN_SITE_RATE, current_rate, MAX_SITE_RATE, TOL_SITE_RATE, &negative_lh, &ferror);
		double fnew;
		if ((optx < MAX_SITE_RATE) && (fnew = computeFunction(MAX_SITE_RATE)) <= negative_lh+TOL_SITE_RATE) {
			optx = MAX_SITE_RATE;
			negative_lh = fnew;
		}
		if ((optx > MIN_SITE_RATE) && (fnew = computeFunction(MIN_SITE_RATE)) <= negative_lh+TOL_SITE_RATE) {
			optx = MIN_SITE_RATE;
			negative_lh = fnew;
		}
	}
	//negative_lh = brent(MIN_SITE_RATE, current_rate, max_rate, 1e-3, &optx);
	if (optx > MAX_SITE_RATE*0.99) optx = MAX_SITE_RATE;
	if (optx < MIN_SITE_RATE*2) optx = MIN_SITE_RATE;
	rates[cat] = optx;
//#ifndef NDEBUG		
//#endif

	if (!rate_mh) completeRateML();
	return optx;	
}

void RateMeyerDiscrete::normalizeRates() {
    double sum = 0.0;
    double ok  = 0.0;
    size_t nptn = size();
    
    auto frequencies = phylo_tree->getConvertedSequenceFrequencies();
    if (frequencies!=nullptr) {
        for (size_t i = 0; i < nptn; i++) {
            //at(i) = rates[ptn_cat[i]];
            if (getPtnRate(i) < MAX_SITE_RATE) {
                double freq = frequencies[i];
                sum += getPtnRate(i) * freq;
                ok  += freq;
            }
        }
    } else {
        for (size_t i = 0; i < nptn; i++) {
            //at(i) = rates[ptn_cat[i]];
            if (getPtnRate(i) < MAX_SITE_RATE) {
                double freq = phylo_tree->aln->at(i).frequency;
                sum += getPtnRate(i) * freq;
                ok  += freq;
            }
        }
    }
    if (fabs(sum - ok) > 1e-3) {
        //cout << "Normalizing rates " << sum << " / " << ok << endl;
        double scale_f = ok / sum;
        for (int i = 0; i < ncategory; i++) {
            if (rates[i] > 2*MIN_SITE_RATE && rates[i] < MAX_SITE_RATE) {
                rates[i] *= scale_f;
            }
        }
    }
}

double RateMeyerDiscrete::classifyRatesKMeans() {

	ASSERT(ncategory > 0);
	int nptn = size();

	// clustering the rates with k-means
	//AddKMeansLogging(&cout, false);
	double *points = new double[nptn];
	int *weights = new int[nptn];
	int i;
	if (!ptn_cat) ptn_cat = new int[nptn];
	for (i = 0; i < nptn; i++) {
		points[i] = at(i);
		if (mcat_type & MCAT_LOG) points[i] = log(points[i]);
		weights[i] = 1;
		if (!(mcat_type & MCAT_PATTERN)) 
			weights[i] = phylo_tree->aln->at(i).frequency;
	}
	memset(rates, 0, sizeof(double) * ncategory);

	//double cost = RunKMeansPlusPlus(nptn, ncategory, 1, points, sqrt(nptn), rates, ptn_cat);
	double cost = RunKMeans1D(nptn, ncategory, points, weights, rates, ptn_cat);
	// assign the categorized rates
	if  (mcat_type & MCAT_LOG) 
		for (i = 0; i < ncategory; i++) rates[i] = exp(rates[i]);
	if (rates[0] < MIN_SITE_RATE) rates[0] = MIN_SITE_RATE;
	if (rates[ncategory-1] > MAX_SITE_RATE - 1e-6) rates[ncategory-1] = MAX_SITE_RATE;
	if (verbose_mode >= VB_MED) {
		cout << "K-means cost: " << cost << endl;
		for (i = 0; i < ncategory; i++) cout << rates[i] << " ";
		cout << endl;
	}

	normalizeRates();
	phylo_tree->clearAllPartialLH();
	double cur_lh = phylo_tree->computeLikelihood();

	delete [] weights;
	delete [] points;
	
	if (mcat_type & MCAT_MEAN)
		return cur_lh;

	return phylo_tree->getModelFactory()->optimizeParameters(false,false, TOL_LIKELIHOOD);

	// optimize category rates again by ML
/*	for (int k = 0; k < 100; k++) {
		phylo_tree->calcDist(dist_mat);
		for (i = 0; i < ncategory; i++)
			optimizeCatRate(i);
		normalizeRates();
		phylo_tree->clearAllPartialLh();
		double new_lh = phylo_tree->optimizeAllBranches(k+2);
		if (new_lh > cur_lh + 1e-2) {
			cur_lh = new_lh; 
			cout << "Current log-likelihood: " << cur_lh << endl;
		} else {
			cur_lh = new_lh;
			break;
		}
	}
	*/
	return cur_lh;
}


double RateMeyerDiscrete::classifyRates(double tree_lh) {
	if (is_categorized) return tree_lh;

	double new_tree_lh;
	is_categorized = true;
	if (ncategory > 0) {
		cout << endl << "Classifying rates into " << ncategory << " categories..." << endl;
		return classifyRatesKMeans();
	}

	// identifying proper number of categories
	int nptn = phylo_tree->aln->getNPattern();
	rates = new double[nptn];

	for (ncategory = 2; ; ncategory++) {
		cout << endl << "Classifying rates into " << ncategory << " categories..." << endl;
		new_tree_lh = classifyRatesKMeans();
		new_tree_lh = phylo_tree->optimizeAllBranches();
		cout << "For " << ncategory << " categories, LogL = " << new_tree_lh;
		double lh_diff = 2*(tree_lh - new_tree_lh);
		int df = (nptn - ncategory);
		double pval = computePValueChiSquare(lh_diff, df);
		cout << ", p-value = " << pval;
		cout << endl;
		//if (new_tree_lh > tree_lh - 3.0) break;
		if (pval > 0.05) break;
	}

	cout << endl << "Number of categories is set to " << ncategory << endl;
	return new_tree_lh;
}




void RateMeyerDiscrete::writeInfo(ostream &out) {
	//out << "Number of categories: " << ncategory << endl;
}

