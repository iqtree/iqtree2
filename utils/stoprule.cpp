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
#include "stoprule.h"
#include "timeutil.h"
#include "MPIHelper.h"

StopRule::StopRule() : CheckpointFactory()
{
//	nTime_ = 0;
	predicted_iteration = 0;

	stop_condition = SC_FIXED_ITERATION;
	confidence_value = 0.95;
	min_iteration = 0;
	max_iteration = 0;
	unsuccess_iteration = 100;
	min_correlation = 0.99;
	step_iteration = 100;
	start_real_time = getRealTime();
	max_run_time = -1.0;
	curIteration = 0;
    should_stop = false;
}

void StopRule::initialize(Params &params) {
	stop_condition = params.stop_condition;
	confidence_value = params.stop_confidence;
	min_iteration = params.min_iterations;
	max_iteration = params.max_iterations;
	unsuccess_iteration = params.unsuccess_iteration;
	min_correlation = params.min_correlation;
	step_iteration = params.step_iterations;
	start_real_time = getRealTime();
	max_run_time = params.maxtime * 60; // maxtime is in minutes
}

void StopRule::getUFBootCountCheck(int &ufboot_count, int &ufboot_count_check) {
    int step = step_iteration;
    while (step*2 < MPIHelper::getInstance().getNumProcesses())
        step *= 2;
    ufboot_count = (curIteration/(step/2)+1)*(step/2);
    ufboot_count_check = (curIteration/step+1)*step;
}

StopRule::~StopRule()
{
}

void StopRule::saveCheckpoint() {
    checkpoint->startStruct("StopRule");
    CKP_SAVE(curIteration);
    CKP_SAVE(start_real_time);
    CKP_VECTOR_SAVE(time_vec);
    checkpoint->endStruct();
    CheckpointFactory::saveCheckpoint();
}

void StopRule::restoreCheckpoint() {
    CheckpointFactory::restoreCheckpoint();
    checkpoint->startStruct("StopRule");
    CKP_RESTORE(curIteration);
    CKP_RESTORE(start_real_time);
    CKP_VECTOR_RESTORE(time_vec);
    checkpoint->endStruct();
}


//
//int StopRule::getNumIterations() {
//	if (stop_condition == SC_FIXED_ITERATION || predicted_iteration == 0)
//		return min_iteration;
//	return predicted_iteration;
//}

//int StopRule::getPredictedIteration(int cur_iteration) {
//	double realtime_secs = getRealTime() - start_real_time;
//
//	switch (stop_condition) {
//	case SC_FIXED_ITERATION:
//		return min_iteration;
//	case SC_WEIBULL:
//		if (predicted_iteration == 0)
//			return min_iteration;
//		else
//			return predicted_iteration;
//	case SC_UNSUCCESS_ITERATION:
//		return getLastImprovedIteration() + unsuccess_iteration;
//	case SC_BOOTSTRAP_CORRELATION:
//		return ((cur_iteration+step_iteration-1)/step_iteration)*step_iteration;
//	case SC_REAL_TIME:
////		return ((max_run_time - realtime_secs)/max_run_time);
//		assert(0);
//		return 0;
//	}
//}

bool StopRule::meetStopCondition(int cur_iteration, double cur_correlation) {
	if (should_stop) {
		return true;
	}
	switch (stop_condition) {
		case SC_FIXED_ITERATION:
			return cur_iteration >= min_iteration;
		case SC_WEIBULL:
			if (predicted_iteration == 0)
				return cur_iteration > min_iteration;
			else
				return cur_iteration > predicted_iteration;
		case SC_UNSUCCESS_ITERATION:
			return cur_iteration > getLastImprovedIteration() + unsuccess_iteration;
		case SC_BOOTSTRAP_CORRELATION:
			return ((cur_correlation >= min_correlation) && (cur_iteration > getLastImprovedIteration() + unsuccess_iteration))
				   || cur_iteration > max_iteration;
		case SC_REAL_TIME:
			return (getRealTime() - start_real_time >= max_run_time);
	}
	return false;
}

double StopRule::getRemainingTime(int cur_iteration) {
	double realtime_secs = getRealTime() - start_real_time;
	int niterations;
	switch (stop_condition) {
	case SC_REAL_TIME:
		return max_run_time - realtime_secs;
	case SC_FIXED_ITERATION:
		niterations = min_iteration;
		break;
	case SC_WEIBULL:
		niterations = (predicted_iteration == 0) ? min_iteration : predicted_iteration;
		break;
	case SC_UNSUCCESS_ITERATION:
		niterations = getLastImprovedIteration() + unsuccess_iteration;
		break;
	case SC_BOOTSTRAP_CORRELATION:
		niterations = max(((cur_iteration+step_iteration-1)/step_iteration)*step_iteration, getLastImprovedIteration() + unsuccess_iteration);
//		if (cur_correlation >= min_correlation)
//			niterations = getLastImprovedIteration() + unsuccess_iteration;
		break;
	}
	return (niterations - cur_iteration) * realtime_secs / (cur_iteration - 1);
}

//void StopRule::setStopCondition(STOP_CONDITION sc) {
//	stop_condition = sc;
//}
//
//void StopRule::setIterationNum(int min_it, int max_it) {
//	min_iteration = min_it;
//	max_iteration = max_it;
//}
//
//void StopRule::setConfidenceValue(double confidence_val)
//{
//	confidence_value = confidence_val;
//	assert(confidence_value > 0 && confidence_value < 1);
//}
//
//void StopRule::setUnsuccessIteration(int unsuccess_iteration) {
//	this->unsuccess_iteration = unsuccess_iteration;
//}
//
//void StopRule::setMinCorrelation(double min_correlation, int step_iteration) {
//	this->min_correlation = min_correlation;
//	this->step_iteration = step_iteration;
//}
//
//void StopRule::setRealTime(double start_real_time, double max_un_time) {
//	this->start_real_time = start_real_time;
//	this->max_run_time = max_run_time;
//}


double StopRule::predict (double &upperTime) {
	if (time_vec.size() < 4) return 0;
	//readVector(time_vec);
	double predictedTime_ = cmpExtinctTime (time_vec.size());
	upperTime = cmpUpperTime (time_vec.size(), 1.0 - confidence_value);
	return predictedTime_;
}

void StopRule::addImprovedIteration(int iteration) {
	time_vec.insert(time_vec.begin(), iteration);
//	nTime_++;
	if (stop_condition != SC_WEIBULL) return;
	double upperTime;
	if (predict(upperTime) == 0) return;
	predicted_iteration = upperTime;
	if (stop_condition == SC_WEIBULL && predicted_iteration > max_iteration)
		predicted_iteration = max_iteration;
	if (predicted_iteration < min_iteration)
			predicted_iteration = min_iteration;
	//cout << "Stopping rule suggests " << predicted_iteration << " iterations ("
	//	<< (predicted_iteration - iteration) << " more iterations)" << endl;
}

int StopRule::getLastImprovedIteration() {
	if (time_vec.empty())
		return 0;
	return time_vec[0];
}

void StopRule::cmpInvMat (DoubleMatrix &oriMat, DoubleMatrix &invMat, int size) {
	//invMat.setLimit (size, size);
	double eps = 1.0e-20; /* ! */
	int i, j, k, l, maxi=0, idx, ix, jx;
	double sum, tmp, maxb, aw;

	invMat.resize(size);
	for (i = 0; i < size; i++) invMat[i].resize(size);

	IntVector index (size);
	double *wk;
	DoubleMatrix omtrx (size);
	for (i = 0; i < size; i++) omtrx[i].resize(size);



	/* copy oriMat to omtrx */
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			omtrx[i][j] = oriMat[i][j];

	wk = (double *) calloc((size_t)size, sizeof(double));
	aw = 1.0;
	for (i = 0; i < size; i++) {
		maxb = 0.0;
		for (j = 0; j < size; j++) {
			if (fabs(omtrx[i][j]) > maxb)
				maxb = fabs(omtrx[i][j]);
		}
		if (maxb == 0.0) {
			/* Singular matrix */
			cout << "\n\n\nHALT: PLEASE REPORT ERROR D TO DEVELOPERS\n\n\n";
			//OutStream::write(oriMat, cout);
			exit(1);
		}
		wk[i] = 1.0 / maxb;
	}
	for (j = 0; j < size; j++) {
		for (i = 0; i < j; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < i; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
		}
		maxb = 0.0;
		for (i = j; i < size; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < j; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
			tmp = wk[i] * fabs(sum);
			if (tmp >= maxb) {
				maxb = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < size; k++) {
				tmp = omtrx[maxi][k];
				omtrx[maxi][k] = omtrx[j][k];
				omtrx[j][k] = tmp;
			}
			aw = -aw;
			wk[maxi] = wk[j];
		}
		index[j] = maxi;
		if (omtrx[j][j] == 0.0)
			omtrx[j][j] = eps;
		if (j != size - 1) {
			tmp = 1.0 / omtrx[j][j];
			for (i = j + 1; i < size; i++)
				omtrx[i][j] *= tmp;
		}
	}
	for (jx = 0; jx < size; jx++) {
		for (ix = 0; ix < size; ix++)
			wk[ix] = 0.0;
		wk[jx] = 1.0;
		l = -1;
		for (i = 0; i < size; i++) {
			idx = index[i];
			sum = wk[idx];
			wk[idx] = wk[i];
			if (l != -1) {
				for (j = l; j < i; j++)
					sum -= omtrx[i][j] * wk[j];
			} else if (sum != 0.0)
				l = i;
			wk[i] = sum;
		}
		for (i = size - 1; i >= 0; i--) {
			sum = wk[i];
			for (j = i + 1; j < size; j++)
				sum -= omtrx[i][j] * wk[j];
			wk[i] = sum / omtrx[i][i];
		}
		for (ix = 0; ix < size; ix++)
			invMat[ix][jx] = wk[ix];
	}
	free((char *)wk);
	wk = NULL;
} /* luinverse */

void StopRule::readMat (char *fileName, DoubleMatrix &oriMat, int &size) {
	std::ifstream inFile_;
	inFile_.open(fileName);
	inFile_ >> size;
	oriMat.resize(size);
	for (int i = 0; i < size; i++) oriMat[i].resize(size);
	for (int row_ = 0; row_ < size; row_ ++)
		for (int col_ = 0; col_ < size; col_ ++)
			inFile_ >> oriMat[row_][col_];
	inFile_.close ();

}

void StopRule::multiple (DoubleMatrix &mat1, DoubleMatrix &mat2, DoubleMatrix &proMat) {
	int row_, col_;
	//proMat.setLimit (mat1.getNRow (), mat2.getNCol ());
	proMat.resize(mat1.size());
	int nrow_ = proMat.size();
	int ncol_ = mat2[0].size();
	for (row_ = 0; row_ < proMat.size(); row_++)   proMat[row_].resize(ncol_);
	for (row_ = 0; row_ < nrow_; row_ ++)
		for (col_ = 0; col_ < ncol_; col_ ++) {
			proMat[row_][col_] = 0.0;
			for (int count_ = 0; count_ < mat1[0].size(); count_ ++) {
				proMat[row_][col_] += mat1[row_][count_] * mat2[count_][col_];
				//         std::cout << mat1[row_][count_] << " --> " << mat2[count_][col_] << endl;
			}
		}
}

void StopRule::multiple (DoubleMatrix &mat1, DoubleVector &vec2, DoubleVector &proVec) {
	int row_, col_;
	proVec.resize(mat1.size());

	for (row_ = 0; row_ < mat1.size (); row_ ++) {
		proVec[row_] = 0.0;
		for (col_ = 0; col_ < mat1[0].size(); col_ ++)
			proVec[row_] += mat1[row_][col_] * vec2[col_];
	}
}

void StopRule::multiple (DoubleVector &vec1, DoubleMatrix &mat2, DoubleVector &proVec) {
	int row_, col_;
	proVec.resize(mat2[0].size());
	for (col_ = 0; col_ < mat2[0].size(); col_ ++) {
		proVec[col_] = 0.0;
		for (row_ = 0; row_ < mat2.size(); row_ ++)
			proVec[col_] += vec1[row_] * mat2[row_][col_];
	}
}

void StopRule::multiple (DoubleVector &vec1, DoubleVector &vec2, DoubleMatrix &proMat) {
	int row_, col_;
	proMat.resize(vec1.size());
	for (row_ = 0; row_ < vec1.size(); row_++)
		proMat[row_].resize(vec2.size());

	for (row_ = 0; row_ < vec1.size(); row_ ++)
		for (col_ = 0; col_ < vec2.size(); col_ ++)
			proMat[row_][col_] = vec1[row_] * vec2[col_];
}

double StopRule::multiple (DoubleVector &vec1, DoubleVector &vec2) {
	double sum_ = 0.0;
	for (int count_ = 0; count_ < vec1.size(); count_ ++)
		sum_ += vec1[count_] * vec2[count_];
	return sum_;
}

/* THE FOLLOWING CODE COMES FROM tools.c in Yang's PAML package */
//----------------------------------------------------------------------------------------
double StopRule::cmpLnGamma (double alpha) {
	/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
	   Stirling's formula is used for the central polynomial part of the procedure.
	   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
	   Communications of the Association for Computing Machinery, 9:684
	*/
	double x=alpha, f=0, z;

	if (x<7) {
		f=1;  z=x-1;
		while (++z<7)  f*=z;
		x=z;   f=-log(f);
	}
	z = 1/(x*x);
	return  f + (x-0.5)*log(x) - x + .918938533204673
	        + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	           +.083333333333333)/x;
} //end of function cmpLnGamma


void StopRule::readVector(DoubleVector &tmpTimeVec_) {
//	nTime_ = tmpTimeVec_.size();
	time_vec.resize(tmpTimeVec_.size());
	for (int count_ = 0; count_ < tmpTimeVec_.size(); count_ ++)
		time_vec[count_] = tmpTimeVec_[tmpTimeVec_.size() - count_ - 1];
}

void StopRule::readFile (const char *fileName) {
	std::ifstream inFile_;
	inFile_.open (fileName);

//	int nTime_ = 0;


	DoubleVector tmpTimeVec_;// (MAX_ITERATION, MAX_ITERATION);

	double old_time = -1.0;
	while (inFile_.eof () == 0) {
		double tmpTime_ = -1.0;
		inFile_ >> tmpTime_;
		if (tmpTime_ > old_time) {
			tmpTimeVec_.push_back(tmpTime_);
//			nTime_ ++;
			old_time = tmpTime_;
		}
	}
	inFile_.close ();

	time_vec.resize(tmpTimeVec_.size());
	for (int count_ = 0; count_ < tmpTimeVec_.size(); count_ ++)
		time_vec[count_] = tmpTimeVec_[tmpTimeVec_.size() - count_ - 1];
}

double StopRule::cmpMuy (int k) {
	double sum_ = 0.0;

	for (int i = 0; i < k - 2; i ++)
		sum_ += log ( (time_vec[0] - time_vec[ k - 1]) / (time_vec[0] - time_vec[i + 1]) );

	double lamda_;
	lamda_ = (1.0 / (k - 1.0) ) * sum_;
	return lamda_;
}


void StopRule::cmpLamdaMat (int k, DoubleMatrix &lamdaMat) {
	int i, j;
	lamdaMat.resize(k);
	for (i = 0; i < k; i ++)
		lamdaMat[i].resize(k);
	double muy_ = cmpMuy (k);
	for (i = 0; i < k; i ++)
		for (j = 0; j <= i; j ++) {
			/*
			lamdaMat[i][j] = (cmpGamma (2*muy_ + i + 1) * cmpGamma (muy_ + j + 1) ) /
			                 ( cmpGamma (muy_ + i + 1) * cmpGamma (j + 1) );*/

			// to fix divide by zero PROBLEM!
			lamdaMat[i][j] = cmpLnGamma (2*muy_ + i + 1) + cmpLnGamma (muy_ + j + 1) -
			                 cmpLnGamma (muy_ + i + 1) - cmpLnGamma (j + 1);

			//if (i == 98 && j == 97)
			 //     std::cout << i << "," << j << " -> " << lamdaMat[i][j] << endl;
			lamdaMat[i][j] = exp(lamdaMat[i][j]);
			lamdaMat[j][i] = lamdaMat[i][j];
		}
}

void StopRule::cmpVecA (int k, DoubleVector &aVec) {
	DoubleVector eVec_ (k, k);
	int count_;
	for (count_ = 0; count_ < k; count_ ++)
		eVec_[count_] = 1.0;

	DoubleMatrix lamdaMat_;
	cmpLamdaMat (k, lamdaMat_);
	 // OutStream::write (lamdaMat_, std::cout);

	DoubleMatrix invLamdaMat_;
	cmpInvMat (lamdaMat_, invLamdaMat_, k);

	//  OutStream::write (invLamdaMat_, std::cout);

	DoubleMatrix proMat_;
	multiple (lamdaMat_, invLamdaMat_, proMat_);
	//OutStream::write (proMat_, std::cout);


	DoubleVector tmp1Vec_;
	multiple (eVec_, invLamdaMat_, tmp1Vec_);
	//OutStream::write (tmp1Vec_, std::cout);
	double tmp2_ = multiple (tmp1Vec_, eVec_);
	double invTmp2_ = 1.0 / tmp2_;

	for (int row_ = 0; row_ < k; row_ ++)
		for (int col_ = 0; col_ < k; col_ ++)
			invLamdaMat_[row_][col_] *= invTmp2_;


	DoubleVector tmp3Vec_;
	multiple (invLamdaMat_, eVec_, aVec);
}

double StopRule::cmpExtinctTime (int k) {
	DoubleVector a;
	cmpVecA (k, a);
	double extinctTime_ = 0.0;
	for (int count_ = 0; count_ < k; count_ ++)
		extinctTime_ += a[count_] * time_vec[count_];
	return extinctTime_;
}



double StopRule::cmpUpperTime (int k, double alpha) {
	double muy_ = cmpMuy (k);
	double priSu_ = -log (alpha) / k ;
	double su_ = pow (priSu_, -muy_);
	return time_vec[0] + (time_vec[0] - time_vec[k - 1]) / (su_ - 1.0);
}

