/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
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

#include "tools.h"

VerboseMode verbose_mode;

//FIXME Where should I put this switches ?
bool phyml_opt;
bool nni_lh;
int speedUpFromIter;
double cmdLambda;
//bool ils;
//int perLevel;

/*
	WIN32 does not define gettimeofday() function.
	Here declare it extra for WIN32 only.
*/
#if defined(WIN32)
#include <sys/timeb.h>
#include <sys/types.h>
#include <winsock.h>
struct timezone {};
void gettimeofday(struct timeval* t,void* timezone)
{       struct _timeb timebuffer;
        _ftime( &timebuffer );
        t->tv_sec=timebuffer.time;
        t->tv_usec=1000*timebuffer.millitm;
}
#else
#include <sys/time.h>
#endif


/********************************************************
	Defining DoubleMatrix methods
********************************************************/

/*DoubleMatrix::DoubleMatrix(int arows, int acols) {
	rows = arows;
	cols = acols;
	size =  rows * cols;
	value = new double[size];
}

void DoubleMatrix::setZero() {
	memset(value, 0, size * sizeof(double));
}


DoubleMatrix::~DoubleMatrix() {
	if (value) delete value;
	value = NULL;
}
*/

/********************************************************
	Miscellaneous
********************************************************/

/**
	Output an error to screen, then exit program
	@param error error message
*/
/*
void outError(char *error)
{
	cerr << "ERROR: " << error << endl;
	exit(2);
}
*/

/**
	Output an error to screen, then exit program
	@param error error message
*/
void outError(const char *error)
{
	cerr << "ERROR: " << error << endl;
	exit(2);
}

/**
	Output an error to screen, then exit program
	@param error error message
*/
void outError(string error)
{
	outError(error.c_str());
}

void outError(const char *error, const char *msg) {
	string str = error;
	str += msg;
	outError(str);
}

void outError(const char *error, string msg) {
	string str = error;
	str += msg;
	outError(str);
}

/**
	Output a warning message to screen
	@param error warning message
*/
void outWarning(const char *warn)
{
	cerr << "WARNING: " << warn << endl;
}

void outWarning(string warn)
{
	outWarning(warn.c_str());
}


double randomLen(Params &params) {
	double ran = static_cast<double> (rand () % 999 + 1) / 1000;
	double len = - params.mean_len * log (ran);

    if (len < params.min_len) {
      int fac = rand () % 1000;
      double delta = static_cast<double> (fac)  / 1000.0; //delta < 1.0
      len = params.min_len + delta / 1000.0;
    }

    if (len > params.max_len) {
      int fac = rand () % 1000;
      double delta = static_cast<double> (fac)  / 1000.0; //delta < 1.0
      len = params.max_len - delta / 1000.0;
    }
	return len;
}

//From Tung
string convertIntToString(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

//From Tung
bool copyFile (const char SRC[], const char DEST[])
{
    std::ifstream src; // the source file
    std::ofstream dest; // the destination file

    src.open (SRC, std::ios::binary); // open in binary to prevent jargon at the end of the buffer
    dest.open (DEST, std::ios::binary); // same again, binary
    if (!src.is_open() || !dest.is_open())
        return false; // could not be copied

    dest << src.rdbuf (); // copy the content
    dest.close (); // close destination file
    src.close (); // close source file

    return true; // file copied successfully
}

bool fileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  return(blnReturn);
}

int convert_int(const char *str) throw (string) {
	char *endptr;
	long i = strtol(str, &endptr, 10);

	if ((i == 0 && ((long) endptr - (const long) str) == 0) || abs(i) == HUGE_VALL || *endptr != 0) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}

	return i;
}

double convert_double(const char *str) throw (string) {
	char *endptr;
	double d = strtod(str, &endptr);
	if ((d == 0.0 && ((long) endptr - (const long) str) == 0) || fabs(d) == HUGE_VALF || *endptr != 0) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	return d;
}

string convert_time(const double sec) {
    int sec_int = (int) floor(sec);
    int secs = sec_int % 60;
    int mins = (sec_int % 3600) / 60;
    int hours = sec_int / 3600;
    stringstream ss;
    ss << hours << ":" << mins << ":" << secs;
    return ss.str();
}

void convert_range(const char *str, int &lower, int &upper, int &step_size) throw (string) {
	char *endptr;
	char *beginptr = (char*) str;

	// parse the lower bound of the range
	int d = strtol(str, &endptr, 10);
	if ((d == 0 && ((long) endptr - (const long) str) == 0) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	//lower = d;
	int d_save = d;
	upper = d;
	if (*endptr == 0) return;


	// parse the upper bound of the range
	str = endptr+1;
	d = strtol(str, &endptr, 10);
	if ((d == 0 && ((long) endptr - (const long) str) == 0) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}

	lower = d_save;
	upper = d;
	if (*endptr == 0) return;

	// parse the step size of the range
	str = endptr+1;
	d = strtol(str, &endptr, 10);
	if ((d == 0 && ((long) endptr - (const long) str) == 0) || abs(d) == HUGE_VALL || *endptr != 0) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}

	step_size = d;
	str = beginptr;

}


void convert_range(const char *str, double &lower, double &upper, double &step_size) throw (string) {
	char *endptr;
	char *beginptr = (char*) str;

	// parse the lower bound of the range
	double d = strtod(str, &endptr);
	if ((d == 0.0 && ((long) endptr - (const long) str) == 0) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	//lower = d;
	double d_save = d;
	upper = d;
	if (*endptr == 0) return;


	// parse the upper bound of the range
	str = endptr+1;
	d = strtod(str, &endptr);
	if ((d == 0.0 && ((long) endptr - (const long) str) == 0) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}

	lower = d_save;
	upper = d;
	if (*endptr == 0) return;

	// parse the step size of the range
	str = endptr+1;
	d = strtod(str, &endptr);
	if ((d == 0.0 && ((long) endptr - (const long) str) == 0) || fabs(d) == HUGE_VALF || *endptr != 0) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}

	step_size = d;
	str = beginptr;

}

void readWeightFile(Params &params, int ntaxa, double &scale, StrVector &tax_name, DoubleVector &tax_weight) {
	cout << "Reading scale factor and taxa weights file " << params.param_file << " ..." << endl;
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(params.param_file);
		string name, tmp;

		in >> tmp;
		scale = convert_double(tmp.c_str());

		for (; !in.eof() && ntaxa > 0; ntaxa--) {
			// remove the failbit
			in.exceptions(ios::badbit);
			if(!(in >> name)) break;
			// set the failbit again
			in.exceptions(ios::failbit | ios::badbit);

			tax_name.push_back(name);
			// read the sequence weight
			in >> tmp;
			tax_weight.push_back(convert_double(tmp.c_str()));
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (string str) {
		outError(str);
	}
}


void readStringFile(const char* filename, int max_num, StrVector &strv) {
	try {
		ifstream in;
		// set the failbit and badbit
		in.exceptions(ios::failbit | ios::badbit);
		in.open(filename);
		string name;

		// remove the failbit
		in.exceptions(ios::badbit);
		for (; !in.eof() && max_num > 0; max_num--) {
			if (!(in >> name)) break;
			strv.push_back(name);
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT);
	}
}

void readInitTaxaFile(Params &params, int ntaxa, StrVector &tax_name) {
	cout << "Reading initial taxa set file " << params.initial_file << " ..." << endl;
	readStringFile(params.initial_file, ntaxa, tax_name);
}

void readInitAreaFile(Params &params, int nareas, StrVector &area_name) {
	cout << "Reading initial area file " << params.initial_area_file << " ..." << endl;
	readStringFile(params.initial_area_file, nareas, area_name);
}

void readAreasBoundary(char *file_name, MSetsBlock *areas, double *areas_boundary) {
	
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(file_name);

		int nset;
		in >> nset;
		if (nset != areas->getNSets())
			throw "File has different number of areas";
		int pos = 0, seq1, seq2;
		for (seq1 = 0; seq1 < nset; seq1 ++)  {
			string seq_name;
			in >> seq_name;
			if (seq_name != areas->getSet(seq1)->name)
				throw "Area name " + seq_name + " is different from " + areas->getSet(seq1)->name;
			for (seq2 = 0; seq2 < nset; seq2 ++) {
				in >> areas_boundary[pos++];
			}	
		}
		// check for symmetric matrix
		for (seq1 = 0; seq1 < nset-1; seq1++) {
			if (areas_boundary[seq1*nset+seq1] <= 1e-6)
				throw "Diagonal elements of distance matrix should represent the boundary of single areas";
			for (seq2 = seq1+1; seq2 < nset; seq2++)
				if (areas_boundary[seq1*nset+seq2] != areas_boundary[seq2*nset+seq1])
					throw "Shared boundary between " + areas->getSet(seq1)->name + " and " + areas->getSet(seq2)->name + " is not symmetric";
		}


		in.close();
		cout << "Areas relation matrix was read from " << file_name << endl;
	} catch (const char *str) {
		outError(str);
	} catch (string str) {
		outError(str);
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, file_name);
	}
	
}

void readTaxaSets(char *filename, MSetsBlock *sets) {
	TaxaSetNameVector *allsets = sets->getSets();
	try {
		int count = 0;
		ifstream in;
		// set the failbit and badbit
		in.exceptions(ios::failbit | ios::badbit);
		in.open(filename);
		string name;

		// remove the failbit
		in.exceptions(ios::badbit);
		while (!in.eof()) {
			int ntaxa = 0;
			string str;
			if (!(in >> str)) break;
			ntaxa = convert_int(str.c_str());
			if (ntaxa <= 0) throw "Number of taxa must be > 0";
			count++;
			//allsets->resize(allsets->size()+1);
			TaxaSetName *myset = new TaxaSetName;
			allsets->push_back(myset);
			myset->name = "";
			myset->name += count;
			for (; ntaxa > 0; ntaxa--) {
				string str;
				if (!(in >> str)) throw "Cannot read in taxon name";
				if ((ntaxa > 1) && in.eof()) throw "Unexpected end of file while reading taxon names";
				myset->taxlist.push_back(str);
			}
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
		if (count == 0) throw "No set found, you must specify at least 1 set";
	} catch (ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (const char *str) {
		outError(str);
	} catch (string str) {
		outError(str);
	}
}



void parseArg(int argc, char *argv[], Params &params) {
	int cnt;
	verbose_mode = VB_MIN;
	params.tree_gen = NONE;
	params.user_file = NULL;
	params.out_prefix = NULL;
	params.out_file = NULL;
	params.sub_size = 0;
	params.pd_proportion = 0.0;
	params.min_proportion = 0.0;
	params.step_proportion = 0.01;
	params.min_size = 0;
	params.step_size = 1;
	params.find_all = false;
	params.run_mode = DETECTED;
	params.detected_mode = DETECTED;
	params.param_file = NULL;
	params.initial_file = NULL;
	params.initial_area_file = NULL;
	params.pdtaxa_file = NULL;
	params.areas_boundary_file = NULL;
	params.boundary_modifier = 1.0;
	params.dist_file = NULL;
	params.compute_ml_dist = false;
	params.compute_ml_tree = true;
	params.budget_file = NULL;
	params.overlap = 0;
	params.is_rooted = false;
	params.sample_size = -1;

	params.repeated_time = 1;
	//params.nr_output = 10000;
	params.nr_output = 0;
	//params.smode = EXHAUSTIVE;
	params.intype = IN_OTHER;
	params.budget = -1;
	params.min_budget = -1;
	params.step_budget = 1;
	params.root = NULL;
	params.num_splits = 0;
	params.min_len = 0.001;
	params.mean_len = 0.1;
	params.max_len = 0.999;
	params.pd_limit = 100;
	params.calc_pdgain = false;
	params.multi_tree = false;
	params.second_tree = NULL;
	params.consensus_type = CT_NONE;
	params.find_pd_min = false;
	params.branch_cluster = 0;
	params.taxa_order_file = NULL;
	params.endemic_pd = false;
	params.exclusive_pd = false;
	params.complement_area = NULL;
	params.scaling_factor = 0.0;
	params.binary_programming = false;
	params.quad_programming = false;
	params.test_input = TEST_NONE;
	params.tree_burnin = 0;
	params.split_threshold = 0.0;
	params.split_weight_summary = SW_SUM;
	params.gurobi_format = false;
	params.gurobi_threads = 1;
	params.num_bootstrap_samples = 0;


	params.aln_file = NULL;
	params.sequence_type = NULL;
	params.aln_output = NULL;
	params.parsimony = false;
	params.tree_spr = false;
	params.nexus_output = false;
	params.k_representative = 4;
	params.p_delete = 0.0;
	params.min_iterations = -1;
	params.max_iterations = 1;
	params.stop_condition = SC_FIXED_ITERATION;
	params.stop_confidence = 0.95;
	params.model_name = "HKY";
	params.store_trans_matrix = false;
	params.freq_type = FREQ_EMPIRICAL;
	//params.freq_type = FREQ_UNKNOWN;
	params.num_rate_cats = 4;
	params.optimize_by_newton = true;
	params.fixed_branch_length = false;
	params.iqp_assess_quartet = IQP_DISTANCE;
	params.write_intermediate_trees = false;
	params.rf_dist_mode = 0;
	params.mvh_site_rate = false;
	params.aLRT_threshold = 101;
	params.aLRT_replicates = 1000;
	params.SSE = true;
	params.print_site_lh = false;			
	params.nni_lh = false;
	params.cmdLambda = 0.75;
        params.speed_conf = 0.75;
	
	struct timeval tv;
	struct timezone tz;
	// initialize random seed based on current time
	gettimeofday(&tv, &tz);
	//params.ran_seed = (unsigned) (tv.tv_sec+tv.tv_usec);
	params.ran_seed = (unsigned) (tv.tv_usec);

	for (cnt = 1; cnt < argc; cnt++) {
		try {

			if (strcmp(argv[cnt],"-h") == 0 || strcmp(argv[cnt],"--help") == 0) {
				usage(argv, false);
			} else if (strcmp(argv[cnt],"-ho") == 0 || strcmp(argv[cnt],"-?") == 0) {
				usage_iqtree(argv, false);
			} else if (strcmp(argv[cnt],"-hh") == 0 || strcmp(argv[cnt],"-hhh") == 0) {
				usage(argv, true);
			} else if (strcmp(argv[cnt],"-v") == 0) {
				verbose_mode = VB_MED;
			} else if (strcmp(argv[cnt],"-vv") == 0) {
				verbose_mode = VB_MAX;
			} else if (strcmp(argv[cnt],"-vvv") == 0) {
				verbose_mode = VB_DEBUG;
			} else if (strcmp(argv[cnt],"-k") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -k <num_taxa>";
				convert_range(argv[cnt], params.min_size, params.sub_size, params.step_size);
				params.k_representative = params.min_size;
			} else if (strcmp(argv[cnt],"-pre") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pre <output_prefix>";
				params.out_prefix = argv[cnt];
			} else if (strcmp(argv[cnt],"-pp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pp <pd_proportion>";
				convert_range(argv[cnt], params.min_proportion, params.pd_proportion, params.step_proportion);
				if (params.pd_proportion < 0 || params.pd_proportion > 1)
					throw "PD proportion must be between 0 and 1";
			} else if (strcmp(argv[cnt],"-mk") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mk <min_taxa>";
				params.min_size = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-bud") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bud <budget>";
				convert_range(argv[cnt], params.min_budget, params.budget, params.step_budget);
			} else if (strcmp(argv[cnt],"-mb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mb <min_budget>";
				params.min_budget = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-o") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -o <taxon>";
				params.root = argv[cnt];
			} else if (strcmp(argv[cnt],"-root") == 0) {
				params.is_rooted = true;
			} else if (strcmp(argv[cnt],"-a") == 0 || strcmp(argv[cnt],"-all") == 0) {
				params.find_all = true;
			} else if (strcmp(argv[cnt],"-g") == 0 || strcmp(argv[cnt],"--greedy") == 0) {
				params.run_mode = GREEDY;
			} else if (strcmp(argv[cnt],"-pr") == 0 || strcmp(argv[cnt],"--pruning") == 0) {
				params.run_mode = PRUNING;
			//} else if (strcmp(argv[cnt],"--both") == 0) {
				//params.run_mode = BOTH_ALG;
			} else if (strcmp(argv[cnt],"-e") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -e <file>";
				params.param_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-i") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -i <file>";
				params.initial_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-ia") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ia <file>";
				params.initial_area_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-u") == 0) {
				// file containing budget information
				cnt++;
				if (cnt >= argc)
					throw "Use -u <file>";
				params.budget_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-dd") == 0) {
				// compute distribution of PD score on random sets
				cnt++;
				if (cnt >= argc)
					throw "Use -dd <sample_size>";
				params.run_mode = PD_DISTRIBUTION;
				params.sample_size = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-ts") == 0) {
				// calculate PD score a taxa set listed in the file
				cnt++;
				//params.run_mode = PD_USER_SET;
				if (cnt >= argc)
					throw "Use -ts <taxa_file>";
				params.pdtaxa_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-bound") == 0) {
				// boundary length of areas
				cnt++;
				if (cnt >= argc)
					throw "Use -bound <file>";
				params.areas_boundary_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-blm") == 0) {
				// boundary length modifier
				cnt++;
				if (cnt >= argc)
					throw "Use -blm <boundary_modifier>";
				params.boundary_modifier = convert_double(argv[cnt]);
			} else if (strcmp(argv[cnt],"-dist") == 0 || strcmp(argv[cnt],"-d") == 0) {
				// calculate distance matrix from the tree
				params.run_mode = CALC_DIST;
				cnt++;
				if (cnt >= argc)
					throw "Use -dist <distance_file>";
				params.dist_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-djc") == 0) {
				params.compute_ml_dist = false;
			} else if (strcmp(argv[cnt],"-r") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -r <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = YULE_HARDING;
			} else if (strcmp(argv[cnt],"-ru") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ru <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = UNIFORM;
			} else if (strcmp(argv[cnt],"-rcat") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcat <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CATERPILLAR;
			} else if (strcmp(argv[cnt],"-rbal") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rbal <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = BALANCED;
			} else if (strcmp(argv[cnt],"-rcsg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcsg <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CIRCULAR_SPLIT_GRAPH;
			} else if (strcmp(argv[cnt],"-rpam") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rpam <num_splits>";
				params.num_splits = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-rlen") == 0) {
				cnt++;
				if (cnt >= argc-2)
					throw "Use -rlen <min_len> <mean_len> <max_len>";
				params.min_len = convert_double(argv[cnt]);
				params.mean_len = convert_double(argv[cnt+1]);
				params.max_len = convert_double(argv[cnt+2]);
				cnt += 2;

			} else if (strcmp(argv[cnt],"-rset") == 0) {
				cnt++;
				if (cnt >= argc-1)
					throw "Use -rset <overlap> <outfile>";
				params.overlap = convert_int(argv[cnt]);
				cnt++;
				params.pdtaxa_file = argv[cnt];
				params.tree_gen = TAXA_SET;
			} else if (strcmp(argv[cnt],"-rep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rep <repeated_times>";
				params.repeated_time = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-lim") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lim <pd_limit>";
				params.pd_limit = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-noout") == 0) {
				params.nr_output = 0;
			} else if (strcmp(argv[cnt],"-1out") == 0) {
				params.nr_output = 1;
			} else if (strcmp(argv[cnt],"-oldout") == 0) {
				params.nr_output = 100;
			} else if (strcmp(argv[cnt],"-nexout") == 0) {
				params.nexus_output = true;
			} else if (strcmp(argv[cnt],"-exhaust") == 0) {
				params.run_mode = EXHAUSTIVE;
			} else if (strcmp(argv[cnt],"-seed") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -seed <random_seed>";
				params.ran_seed = (unsigned)convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-pdgain") == 0) {
				params.calc_pdgain = true;
			} else if (strcmp(argv[cnt],"-sup") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sup <target_tree_file>";
				params.second_tree = argv[cnt];
				params.consensus_type = CT_ASSIGN_SUPPORT;
			} else if (strcmp(argv[cnt],"-con") == 0) {
				params.consensus_type = CT_CONSENSUS_TREE;
			} else if (strcmp(argv[cnt],"-net") == 0) {
				params.consensus_type = CT_CONSENSUS_NETWORK;
			} 
			/**MINH ANH: to serve some statistics*/
			else if (strcmp(argv[cnt],"-comp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -comp <treefile>";
				params.consensus_type = COMPARE;
				params.second_tree = argv[cnt];
			}else if (strcmp(argv[cnt],"-stats") == 0) {
				params.run_mode = STATS;
				
			} // MA
			else if (strcmp(argv[cnt],"-min") == 0) {
				params.find_pd_min = true;
			} else if (strcmp(argv[cnt],"-excl") == 0) {
				params.exclusive_pd = true;
			} else if (strcmp(argv[cnt],"-endem") == 0) {
				params.endemic_pd = true;
			} else if (strcmp(argv[cnt],"-compl") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -compl <area_name>";
				params.complement_area = argv[cnt];
			} else if (strcmp(argv[cnt],"-cluster") == 0) {
				params.branch_cluster = 4;
				cnt++;
				if (cnt >= argc)
					throw "Use -cluster <taxa_order_file>";
				params.taxa_order_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-taxa") == 0) {
				params.run_mode = PRINT_TAXA;
			} else if (strcmp(argv[cnt],"-area") == 0) {
				params.run_mode = PRINT_AREA;
			} else if (strcmp(argv[cnt],"-scalebranch") == 0) {
				params.run_mode = SCALE_BRANCH_LEN;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalebranch <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
			} else if (strcmp(argv[cnt],"-scalenode") == 0) {
				params.run_mode = SCALE_NODE_NAME;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalenode <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
			} else if (strcmp(argv[cnt],"-lp") == 0) {
				params.run_mode = LINEAR_PROGRAMMING;
			} else if (strcmp(argv[cnt],"-lpbin") == 0) {
				params.run_mode = LINEAR_PROGRAMMING;
				params.binary_programming = true;
			} else if (strcmp(argv[cnt],"-qp") == 0) {
				params.gurobi_format = true;
				params.quad_programming = true;
			} else if (strcmp(argv[cnt],"-mult") == 0) {
				params.multi_tree = true;
			} else if (strcmp(argv[cnt],"-bi") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bi <burnin_value>";
				params.tree_burnin = convert_int(argv[cnt]);
				if (params.tree_burnin < 0)
					throw "Burnin value must not be negative";
			} else if (strcmp(argv[cnt],"-t") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -t <split_threshold>";
				params.split_threshold = convert_double(argv[cnt]);
				if (params.split_threshold < 0 || params.split_threshold > 1)
					throw "Split threshold must be between 0 and 1";
			} else if (strcmp(argv[cnt],"-swc") == 0) {
				params.split_weight_summary = SW_COUNT;			
			} else if (strcmp(argv[cnt],"-swa") == 0) {
				params.split_weight_summary = SW_AVG_ALL;
			} else if (strcmp(argv[cnt],"-swp") == 0) {
				params.split_weight_summary = SW_AVG_PRESENT;
			} else if (strcmp(argv[cnt],"-iwc") == 0) {
				params.test_input = TEST_WEAKLY_COMPATIBLE;
			} else if (strcmp(argv[cnt],"-aln") == 0 || strcmp(argv[cnt],"-s") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -aln, -s <alignment_file>";
				params.aln_file = argv[cnt];
			} else if (strcmp(argv[cnt],"-st") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -st <BIN|DNA|AA>";
				params.sequence_type = argv[cnt];
			} else if (strcmp(argv[cnt],"-alnout") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -alnout <alignment_file>";
				params.aln_output = argv[cnt];
			} else if (strcmp(argv[cnt],"-pars") == 0) {
				// maximum parsimony
				params.parsimony = true;
			} else if (strcmp(argv[cnt],"-spr") == 0) {
				// subtree pruning and regrafting
				params.tree_spr = true;
			} else if (strcmp(argv[cnt],"-krep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -krep <num_k>";
				params.k_representative = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-pdel") == 0 || strcmp(argv[cnt],"-p") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pdel <probability>";
				params.p_delete = convert_double(argv[cnt]);
				if (params.p_delete < 0.0 || params.p_delete > 1.0)
					throw "Probability of deleting a leaf must be between 0 and 1";
			} else if (strcmp(argv[cnt],"-n") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -n <#iterations>";
				params.min_iterations = convert_int(argv[cnt]);
			} else if (strcmp(argv[cnt],"-mod") == 0 || strcmp(argv[cnt],"-m") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mod <model_name>";
				params.model_name = argv[cnt];
			} else if (strcmp(argv[cnt],"-mvh") == 0) {
				params.mvh_site_rate = true;
			} else if (strcmp(argv[cnt],"-mstore") == 0) {
				params.store_trans_matrix = true;
			} else if (strcmp(argv[cnt], "-phyml_opt") == 0) {
				phyml_opt = true;
			} else if (strcmp(argv[cnt], "-nni_lh") == 0) {
				params.nni_lh = true;
			} else if (strcmp(argv[cnt], "-lambda") == 0) {
				cnt++;
				params.cmdLambda = convert_double(argv[cnt]);
                                if (params.cmdLambda > 1.0)
                                    throw "Lambda must be in (0,1]";
			}
                        else if (strcmp(argv[cnt], "-spc") == 0) {
				cnt++;
                                if (cnt >= argc)
                                    throw "Please specify the speed-up confidence level";
                                params.speed_conf = convert_double(argv[cnt]);
                                if (  params.speed_conf <= 0.5 || params.speed_conf >1 )
                                    throw "Speed up confidence level must be bigger than 0.5 and smaller or equal 1";
			} else if (strcmp(argv[cnt], "-nosse") == 0) {
				params.SSE = false;
                        }
			else if (strcmp(argv[cnt],"-f") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -f <EQ | EM | ES | UD>";
				if (strcmp(argv[cnt],"EQ") == 0)
					params.freq_type = FREQ_EQUAL;
				else if (strcmp(argv[cnt],"EM") == 0)
					params.freq_type = FREQ_EMPIRICAL;
				else if (strcmp(argv[cnt],"ES") == 0)
					params.freq_type = FREQ_ESTIMATE;
				else if (strcmp(argv[cnt],"UD") == 0)
					params.freq_type = FREQ_USER_DEFINED;
				else
					throw "Use -f <EQ | EM | ES | UD>";

			} else if (strcmp(argv[cnt],"-c") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -c <#rate_category>";
				params.num_rate_cats = convert_int(argv[cnt]);
				if (params.num_rate_cats < 1) throw "Wrong number of rate categories";
			} else if (strcmp(argv[cnt],"-brent") == 0) {
				params.optimize_by_newton = false;
			} else if (strcmp(argv[cnt],"-fixbr") == 0) {
				params.fixed_branch_length = true;
			} else if (strcmp(argv[cnt],"-sr") == 0) {
				params.stop_condition = SC_STOP_PREDICT;
				cnt++;
				if (cnt >= argc)
					throw "Use -sr <#max_iteration>";
				params.max_iterations = convert_int(argv[cnt]);
				if (params.max_iterations <= params.min_iterations)
					throw "Specified max iteration must be greater than min iteration";
			} else if (strcmp(argv[cnt],"-sc") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sc <stop_confidence_value>";
				params.stop_confidence = convert_double(argv[cnt]);
				if (params.stop_confidence <= 0.5 || params.stop_confidence >= 1)
					throw "Stop confidence value must be in range (0.5,1)";
			} else if (strcmp(argv[cnt],"-gurobi") == 0) {
				params.gurobi_format = true;
			} else if (strcmp(argv[cnt],"-gthreads") == 0) {
				params.gurobi_format = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -gthreads <gurobi_threads>";
				params.gurobi_threads = convert_int(argv[cnt]);
				if (params.gurobi_threads < 1)
					throw "Wrong number of threads";
			} else if (strcmp(argv[cnt],"-b") == 0) {
				params.multi_tree = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -b <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1) 
					throw "Wrong number of bootstrap samples";
			} else if (strcmp(argv[cnt],"-bo") == 0) {
				params.multi_tree = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -b <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1) 
					throw "Wrong number of bootstrap samples";
				params.compute_ml_tree = false;
			} else if (strcmp(argv[cnt],"-iqppars") == 0) {
				params.iqp_assess_quartet = IQP_PARSIMONY;
			} else if (strcmp(argv[cnt],"-wt") == 0) {
				params.write_intermediate_trees = true;
			} else if (strcmp(argv[cnt],"-rf_all") == 0) {
				params.rf_dist_mode = RF_ALL_PAIR;
			} else if (strcmp(argv[cnt],"-rf_adj") == 0) {
				params.rf_dist_mode = RF_ADJACENT_PAIR;
			} else if (strcmp(argv[cnt],"-aLRT") == 0) {
				cnt++;
				if (cnt+1 >= argc)
					throw "Use -aLRT <threshold%> <#replicates>";
				params.aLRT_threshold = convert_int(argv[cnt]);
				if (params.aLRT_threshold < 85 || params.aLRT_threshold > 101) 
					throw "aLRT thresold must be between 85 and 100";
				cnt++;
				params.aLRT_replicates = convert_int(argv[cnt]);
				if (params.aLRT_replicates < 1000 && params.aLRT_replicates != 0) 
					throw "aLRT replicates must be at least 1000";
			} else if (strcmp(argv[cnt],"-wsl") == 0) {
				params.print_site_lh = true;
			} else if (argv[cnt][0] == '-') {
				string err = "Invalid \"";
				err += argv[cnt];
				err += "\" option.";
				throw err;
			} else {
				if (params.user_file == NULL)
					params.user_file = argv[cnt];
				else
					params.out_file = argv[cnt];
			}
			if (params.root != NULL && params.is_rooted)
				throw "Not allowed to specify both -o <taxon> and -root";

		}
		 
		// try
		catch (const char *str) {
			outError(str);
		//} catch (char *str) {
			//outError(str);
		} catch (string str) {
			outError(str);
		} catch (...) {
			string err = "Unknown argument \"";
			err += argv[cnt];
			err += "\"";
			outError(err);
		}

	} // for
	if (params.user_file == NULL && params.aln_file == NULL)
		usage(argv, false);
	if (!params.out_prefix) {
		if (params.aln_file) 
			params.out_prefix = params.aln_file;
		else 
			params.out_prefix = params.user_file;
	}
}

void usage(char* argv[], bool full_command) {
	cout << "Usage: " << argv[0] << " [OPTIONS] <file_name> [<output_file>]" << endl;
	cout << "GENERAL OPTIONS:" << endl;
	cout << "  -h                Print this help dialog. Use -hh to display all options." << endl;
	cout << "  -?                Print help options for phylogenetic inference." << endl;
	cout << "  <file_name>       User tree in NEWICK format or split network in NEXUS format." << endl;
	cout << "  <output_file>     Output file to store results, default is '<file_name>.pda'." << endl;
	cout << "  -k <num_taxa>     Find optimal PD set of size <num_taxa>." << endl;
	cout << "  -k <min>:<max>    Find optimal PD sets of size from <min> to <max>." << endl;
	cout << "  -k <min>:<max>:<step>" << endl;
	cout << "    Find optimal PD sets of size <min>, <min>+<step>, <min>+2*<step>,..." << endl;
	cout << "  -o <taxon>        Root name to compute rooted PD, default is unrooted. " << endl;
	cout << "  -i <file>         File containing taxa to be included into PD set." << endl;
	cout << "  -e <file>         File containing branch/split scale and taxa weights." << endl;
	cout << "  -a, -all          Identify multiple optimal PD sets." << endl;
	cout << "  -lim <max_limit>  The maximum number of PD sets for each k if -a is specified." << endl;
	cout << "  -min              Compute minimal PD sets, default is maximal PD sets." << endl;
	cout << "  -1out             Also print taxa sets and scores to separate files." << endl;
	cout << "  -oldout           Also print output compatible with version 0.3." << endl;
	cout << "  -v                Verbose mode." << endl;
	cout << endl;
	cout << "OPTIONS FOR TREE:" << endl;
	cout << "  -root             Make the tree ROOTED, default is unrooted." << endl;
	cout << "    NOTE: this option and -o <taxon> cannot be both specified." << endl;
	cout << "  -g, --greedy      Run greedy algorithm only." << endl;
	cout << "  -pr, --pruning    Run pruning algorithm only." << endl;
	cout << "    NOTE: by default, the program automatically chooses suitable algorithm." << endl;
	cout << endl;
	cout << "OPTIONS FOR SPLIT-NETWORK:" << endl;
	cout << "  -exhaust          Force to use exhaustive search." << endl;
	cout << "    NOTE: by default, the program applies dynamic programming algorithm" << endl;
	cout << "          on circular networks and exhaustive search on general networks." << endl;
	cout << endl;
	cout << "OPTIONS FOR BUDGET-CONSTRAINT:" << endl;
	cout << "  -u <file>         File containing total budget and taxa preservation costs." << endl;
	cout << "  -b <budget>       Total budget to conserve taxa." << endl;
	cout << "  -b <min>:<max>    Find all PD sets with budget from <min> to <max>." << endl;
	cout << "  -b <min>:<max>:<step>" << endl;
	cout << "    Find optimal PD sets with budget <min>, <min>+<step>, <min>+2*<step>,..." << endl;
	cout << endl;
	cout << "OPTIONS FOR AREA ANALYSIS:" << endl;
	cout << "  -ts <taxa_file>   Compute PD of areas (user-defined sets) in <taxa_file>." << endl;
	cout << "  -excl             Compute area exclusive PD." << endl;
	cout << "  -endem            Compute area endemic PD." << endl;
	cout << "  -compl <areas>    Compute area PD-complementarity given the listed <areas>." << endl;
	cout << endl;

	if (!full_command) exit(0);

	cout << "GENERATING RANDOM TREES:" << endl;
	cout << "  -r <num_taxa>     Create a random tree under Yule-Harding model." << endl;
	cout << "  -ru <num_taxa>    Create a random tree under Uniform model." << endl;
	cout << "  -rcat <num_taxa>  Create a random caterpillar tree." << endl;
	cout << "  -rbal <num_taxa>  Create a random balanced tree." << endl;
	cout << "  -rcsg <num_taxa>  Create a random circular split network." << endl;
	cout << "  -rlen <min_len> <mean_len> <max_len>  " << endl;
	cout << "      minimum, mean, and maximum branch lengths of the random trees." << endl;
	cout << endl;

	cout << "MISCELLANEOUS:" << endl;
	cout << "  -dd <sample_size> Compute PD distribution of random sets of size k." << endl;
	cout << "  -d <outfile>      Calculate the distance matrix inferred from tree." << endl;
	cout << "  -seed <number>    Set the seed for random number generator." << endl;
	cout << "  -stats <outfile>  Output some statistics about branch lengths on the tree." << endl;
	cout << "  -comp <treefile> Compare the tree with each in the input trees." << endl; 


//	cout << "  -rep <times>        Repeat algorithm a number of times." << endl;
//	cout << "  -noout              Print no output file." << endl;
	cout << endl;
	//cout << "HIDDEN OPTIONS: see the source code file pda.cpp::parseArg()" << endl;

	exit(0);
}

void usage_iqtree(char* argv[], bool full_command) {
	cout << "Usage: " << argv[0] << " -s <alignment> [OPTIONS] [<tree_file>] " << endl << endl;
	cout << "GENERAL OPTIONS:" << endl
		 << "  -?                   Print this help dialog" << endl
		 << "  -s <alignment>       Input alignment (REQUIRED) in PHYLIP or NEXUS format"  << endl
		 << "  -st <BIN|DNA|AA>     Binary, DNA, or Protein sequences (default: auto-detect)"  << endl
		 << "  -b <#replicates>     Non-parametric bootstrap (default: none)" << endl
		 << "  <tree_file>          Initial tree for tree reconstruction (default: BIONJ)" << endl
		 << "                       Or set of trees for consensus reconstruction (see below)" << endl
		 << "  -o <outgroup_taxon>  Outgroup taxon name, used when writing .treefile" << endl 
		 << "  -pre <PREFIX>        Use <PREFIX> for all output files (default: alignment" << endl
		 << "                       file name)" << endl
	<< endl << "SUBSTITUTION MODEL OPTIONS:" << endl
		 << "  -m <substitution_model_name>" << endl
		 << "                  DNA: JC (default), F81, K2P, HKY, K3P, K81uf, TN/TrN, TNef," << endl
		 << "                       TIM, TIMef, HKY, TVM, TVMef, SYM, GTR, or 6-letter model" << endl
		 << "                       specification, e.g., '010010' is equivalent to HKY" << endl
		 << "              Protein: JC (default), WAG, cpREV, mtREV, PAM, mtMAM, JTT, LG," << endl
		 << "                       mtART, mtZOA, VT, or rtREV" << endl
		 << "               Binary: JC (default)" << endl
		 << "   Rate heterogeneity: Append '+I', '+G[n]', or '+I+G[n]' to the model name for" << endl
		 << "                       Invar, Gamma, or Invar+Gamma rates. 'n' is the number of" << endl
		 << "                       categories for Gamma rates (default: n=4)" << endl
		 << "            Modeltest: TEST or TESTONLY if you have no idea which model to use." << endl
		 << "                       With 'TEST', IQTREE will evaluate all models and select" << endl
		 << "                       the most suitable one. With 'TESTONLY' IQTREE will stop" << endl
		 << "                       after finishing Modeltest" << endl
		 << "            Otherwise: Name of the file containing user model parameters" << endl
		 << "                       incl. rate matrix and state frequencies" << endl
		 //<< "  -c <#categories>     Number of Gamma rate categories (default: 4)" << endl
		 << "  -f <EQ|EM|ES|UD>     EQual, EMpirical, EStimated, or User-Defined state" << endl 
		 << "                       frequency (default: detected from the model name)" << endl
	<< endl << "TREE INFERENCE OPTIONS:" << endl
		 << "  -p <probability>     IQP: Probability of deleting a leaf (default: 0.1)" << endl
		 << "  -k <#representative> IQP: Size of representative leaf set (default: 5)" << endl
		 << "  -n <#iterations>     Number of iterations  (default: 1, equiv. to PHYML)" << endl
		 << "  -sr <#iterations>    Turn on stopping rule, not exceeding #iterations" << endl
		 << "  -sc <confidence>     Confidence value for stopping rule (default: 0.95)" << endl
		 << "  -wt                  Write all intermediate trees into .treels file" << endl
		 << "  -d <file>            Read the genetic distances from the file (default: " << endl
		 << "                       Juke-Cantor distances)" << endl
		 << "  -fixbr               Fix branch lengths of <tree_file>, only estimate" << endl
		 << "                       model parameters" << endl
		 << "  -seed <number>       Random seed number, normally used for debugging purpose" << endl
		 << "  -v, -vv, -vvv        Verbose mode, print more messages to the screen" << endl
	<< endl << "CONSENSUS RECONSTRUCTION OPTIONS:" << endl
		 << "  -t <threshold>       Number between 0 and 1 for minimum split support." << endl
		 << "                       Use '-t 0.5 -con' for majority-rule consensus'" << endl
		 << "                       (default: 0, i.e. extended consensus)" << endl
		 << "  -bi <burnin>         Discard <burnin> trees at the beginning of <tree_file>" << endl
		 << "  -con                 Compute consensus tree of the trees in <tree_file>." << endl
		 << "                       Output to <tree_file>.contree or <PREFIX>.contree" << endl
		 << "  -net                 Compute consensus network of the trees in <trees_file>." << endl
		 << "                       Output to <tree_file>.nex or <PREFIX>.nex" << endl
		 << "  -sup <target_tree>   Assign support values for branches in <target_tree>" << endl
		 << "                       based on the trees from <tree_file>. Output to" << endl
		 << "                       <target_tree>.suptree or <PREFIX>.suptree" << endl
	<< endl << "ROBINSON-FOULDS DISTANCE OPTIONS:" << endl
		 << "  -rf_all              Compute all-to-all Robinson-Foulds distances of the" << endl
		 << "                       trees in <tree_file>. Output to .rfdist file" << endl
		 << "  -rf_adj              Compute Robinson-Foulds distances between adjacent pairs" << endl
		 << "                       of the trees in <tree_file>. Output to .rfdist file" << endl
		 << endl;

	if (full_command) {
		//TODO Print other options here (to be added)
	}

	exit(0);
}
InputType detectInputFile(char *input_file) {

	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(input_file);

		unsigned char ch;
		int count=0;
		do {
			in >> ch;
		} while (ch <= 32 && !in.eof() && count++ < 20);
		in.close();
		return (ch == '#') ? IN_NEXUS : ( (ch == '(' || ch == '[') ? IN_NEWICK : IN_OTHER );
	} catch (ios::failure) {
		outError("Cannot read file ", input_file);
	}
	return IN_OTHER;
}


bool overwriteFile(char *filename) {
	ifstream infile(filename);
	if (infile.is_open())
	{
		cout << "Overwrite " << filename << " (y/n)? ";
		char ch;
		cin >> ch;
		if (ch != 'Y' && ch != 'y') {
			infile.close();
			return false;
		}
	}
	infile.close();
	return true;
}

void parseAreaName(char *area_names, set<string> &areas) {
	string all = area_names;
	int pos;
	while (!all.empty()) {
		pos = all.find(',');
		if (pos < 0) pos = all.length();
		areas.insert(all.substr(0, pos));
		if (pos >= all.length())
			all = "";
		else
			all = all.substr(pos+1);
	}
}
