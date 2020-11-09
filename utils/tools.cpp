/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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
#include "starttree.h" //for START_TREE_RECOGNIZED macro.
#include "timeutil.h"
#include "MPIHelper.h"
#ifndef CLANG_UNDER_VS
    #include <dirent.h>
#else
     //James B. Workaround for Windows builds where these macros might not be defined
    #ifndef S_ISDIR
    #define S_ISDIR(mode) (((mode) & S_IFMT) == S_IFDIR)
    #endif
    #ifndef S_ISREG
    #define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
    #endif
#endif
#include <thread>


#if defined(Backtrace_FOUND)
#include <execinfo.h>
#include <cxxabi.h>
#endif

#include "tools.h"
#include "timeutil.h"
#include "progress.h"
#include "gzstream.h"
#include "MPIHelper.h"
#include "alignment/alignment.h"

VerboseMode verbose_mode;
extern void printCopyright(ostream &out);

#if defined(WIN32)
#include <sstream>
#endif

/********************************************************
        Miscellaneous
 ********************************************************/

/**
        Output an error to screen, then exit program
        @param error error message
 */
void outError(const char *error, bool quit) {
	if (error == ERR_NO_MEMORY) {
        print_stacktrace(cerr);
	}
	cerr << error << endl;
    if (quit)
    	exit(2);
}

/**
        Output an error to screen, then exit program
        @param error error message
 */
void outError(string error, bool quit) {
    outError(error.c_str(), quit);
}

void outError(const char *error, const char *msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

void outError(const char *error, string msg, bool quit) {
    string str = error;
    str += msg;
    outError(str, quit);
}

/**
        Output a warning message to screen
        @param error warning message
 */
void outWarning(const char *warn) {
    cout << "WARNING: " << warn << endl;
}

void outWarning(string warn) {
    outWarning(warn.c_str());
}

double randomLen(Params &params) {
    double ran = static_cast<double> (random_int(999) + 1) / 1000;
    double len = -params.mean_len * log(ran);

    if (len < params.min_len) {
        int fac = random_int(1000);
        double delta = static_cast<double> (fac) / 1000.0; //delta < 1.0
        len = params.min_len + delta / 1000.0;
    }

    if (len > params.max_len) {
        int fac = random_int(1000);
        double delta = static_cast<double> (fac) / 1000.0; //delta < 1.0
        len = params.max_len - delta / 1000.0;
    }
    return len;
}


//From Tung

string convertIntToString(int number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertInt64ToString(int64_t number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

string convertDoubleToString(double number) {
    stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

bool iEquals(const string a, const string b)
{
    unsigned int sz = a.size();
    if (b.size() != sz)
        return false;
    for (unsigned int i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}

//From Tung

bool copyFile(const char SRC[], const char DEST[]) {
    std::ifstream src; // the source file
    std::ofstream dest; // the destination file

    src.open(SRC, std::ios::binary); // open in binary to prevent jargon at the end of the buffer
    dest.open(DEST, std::ios::binary); // same again, binary
    if (!src.is_open() || !dest.is_open())
        return false; // could not be copied

    dest << src.rdbuf(); // copy the content
    dest.close(); // close destination file
    src.close(); // close source file

    return true; // file copied successfully
}

bool fileExists(string strFilename) {
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0) {
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
    return (blnReturn);
}

int isDirectory(const char *path) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0)
        return 0;
    return S_ISDIR(statbuf.st_mode);
}

int isFile(const char *path) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0)
        return 0;
    return S_ISREG(statbuf.st_mode);
}

size_t getFilesInDir(const char *path, StrVector &filenames)
{
    if (!isDirectory(path)) {
        return 0;
    }
    string path_name = path;
    if (path_name.back() != '/') {
        path_name.append("/");
    }
    size_t oldCount = filenames.size();
#ifndef CLANG_UNDER_VS
    DIR* dp = opendir (path);
    if (dp == nullptr) {
        return 0;
    }
    struct dirent* ep;
    while ((ep = readdir (dp)) != NULL) {
        if (isFile((path_name + ep->d_name).c_str()))
            filenames.push_back(ep->d_name);
    }
    (void) closedir (dp);
    
#else
    path_name += "*";
    WIN32_FIND_DATA find_data;
    HANDLE search_handle = FindFirstFile(path_name.c_str(), &find_data);
    if (search_handle == INVALID_HANDLE_VALUE) {
        return 0;
    }
    do {
        if ((find_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0) {
            filenames.emplace_back(find_data.cFileName);
        }
    } while (FindNextFile(search_handle, &find_data));
    FindClose(search_handle);
#endif
    return filenames.size() - oldCount;
}
int convert_int(const char *str) {
    char *endptr;
    int i = strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    return i;
}

int convert_int_nothrow(const char* str, int defaultValue) throw() {
    char *endptr;
    int i = strtol(str, &endptr, 10);
    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        return defaultValue;
    }
    return i;
}


int convert_int(const char *str, int &end_pos) {
	char *endptr;
	int i = strtol(str, &endptr, 10);

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return i;
}

void convert_int_vec(const char *str, IntVector &vec) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		int i = strtol(beginptr, &endptr, 10);

		if ((i == 0 && endptr == beginptr) || abs(i) == HUGE_VALL) {
			string err = "Expecting integer, but found \"";
			err += beginptr;
			err += "\" instead";
			throw err;
		}
		vec.push_back(i);
		if (*endptr == ',') endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}


int64_t convert_int64(const char *str) {
    char *endptr;
    int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting large integer , but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    return i;
}

int64_t convert_int64(const char *str, int &end_pos) {
	char *endptr;
	int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

	if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL) {
		string err = "Expecting large integer, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return i;
}


double convert_double(const char *str) {
    char *endptr;
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    return d;
}

double convert_double_nothrow(const char *str, double defaultValue) throw() {
    char *endptr;
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        return defaultValue;
    }
    return d;
}


double convert_double(const char *str, int &end_pos) {
	char *endptr;
	double d = strtod(str, &endptr);
	if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF) {
		string err = "Expecting floating-point number, but found \"";
		err += str;
		err += "\" instead";
		throw err;
	}
	end_pos = endptr - str;
	return d;
}

void convert_double_vec(const char *str, DoubleVector &vec, char separator) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    do {
		double d = strtod(beginptr, &endptr);

		if ((d == 0.0 && endptr == beginptr) || fabs(d) == HUGE_VALF) {
			string err = "Expecting floating-point number, but found \"";
			err += beginptr;
			err += "\" instead";
			throw err;
		}
		vec.push_back(d);
		if (*endptr == separator) endptr++;
		beginptr = endptr;
    } while (*endptr != 0);
}

string convert_time(const double sec) {
    int sec_int = (int) floor(sec);
    int secs = sec_int % 60;
    int mins = (sec_int % 3600) / 60;
    int hours = sec_int / 3600;
    stringstream ss;
    ss << hours << "h:" << mins << "m:" << secs << "s";
    return ss.str();
}

void convert_range(const char *str, int &lower, int &upper, int &step_size) {
    char *endptr;

    // parse the lower bound of the range
    int d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
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
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || *endptr != 0) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    step_size = d;
}

void convert_range(const char *str, double &lower, double &upper, double &step_size) {
    char *endptr;

    // parse the lower bound of the range
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
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
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || (*endptr != 0 && *endptr != ':')) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    lower = d_save;
    upper = d;
    if (*endptr == 0) return;

    // parse the step size of the range
    str = endptr + 1;
    d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    step_size = d;
}

void convert_string_vec(const char *str, StrVector &vec, char separator) {
    char *beginptr = (char*)str, *endptr;
    vec.clear();
    string elem;
    do {
    	endptr = strchr(beginptr, separator);
    	if (!endptr) {
    		elem.assign(beginptr);
    		vec.push_back(elem);
    		return;
    	}
    	elem.assign(beginptr, endptr-beginptr);
    	vec.push_back(elem);
		beginptr = endptr+1;
    } while (*endptr != 0);

}

bool renameString(string &name) {
    bool renamed = false;
    for (string::iterator i = name.begin(); i != name.end(); i++) {
        if ((*i) == '/') {
            // PLL does not accept '/' in names, turn it off
            if (Params::getInstance().start_tree == STT_PLL_PARSIMONY)
                Params::getInstance().start_tree = STT_PARSIMONY;
        }
        if (!isalnum(*i) && (*i) != '_' && (*i) != '-' && (*i) != '.' && (*i) != '|' && (*i) != '/') {
            (*i) = '_';
            renamed = true;
        }
    }
    return renamed;
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
            if (!(in >> name)) break;
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
    } catch (ios::failure) {
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

void printString2File(string myString, string filename) {
    ofstream myfile(filename.c_str());
    if (myfile.is_open()) {
        myfile << myString;
        myfile.close();
    } else {
        cout << "Unable to open file " << filename << endl;
    }
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
        for (seq1 = 0; seq1 < nset; seq1++) {
            string seq_name;
            in >> seq_name;
            if (seq_name != areas->getSet(seq1)->name)
                throw "Area name " + seq_name + " is different from " + areas->getSet(seq1)->name;
            for (seq2 = 0; seq2 < nset; seq2++) {
                in >> areas_boundary[pos++];
            }
        }
        // check for symmetric matrix
        for (seq1 = 0; seq1 < nset - 1; seq1++) {
            if (areas_boundary[seq1 * nset + seq1] <= 1e-6)
                throw "Diagonal elements of distance matrix should represent the boundary of single areas";
            for (seq2 = seq1 + 1; seq2 < nset; seq2++)
                if (areas_boundary[seq1 * nset + seq2] != areas_boundary[seq2 * nset + seq1])
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
            string str_taxa;
            if (!(in >> str_taxa)) break;
            ntaxa = convert_int(str_taxa.c_str());
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

void get2RandNumb(const int size, int &first, int &second) {
    // pick a random element
    first = random_int(size);
    // pick a random element from what's left (there is one fewer to choose from)...
    second = random_int(size - 1);
    // ...and adjust second choice to take into account the first choice
    if (second >= first) {
        ++second;
    }
}

void quickStartGuide();

namespace {
    std::string string_to_lower(const char* input) {
        std::string answer = input;
        std::transform(answer.begin(), answer.end(), answer.begin(),
                       []( char c){ return std::tolower(c); });
        return answer;
    }

    std::string string_to_upper(const char* input) {
        std::string answer = input;
        std::transform(answer.begin(), answer.end(), answer.begin(),
                       []( char c){ return std::toupper(c); });
        return answer;
    }

    std::string next_argument(int argc, char* argv[], const char* desc, int& cnt ) {
        cnt++;
        if (cnt >= argc) {
            std::string problem = std::string("Use ") + argv[cnt-1] + " <" + desc + ">";
            throw problem;
        }
        return argv[cnt];
    }


    template <class V, class S> void throw_if_not_in_set
    ( const char* name, const V& value, S set, size_t setCount ) {
        for ( size_t i=0; i<setCount; ++i ) {
            if (value == set[i]) {
                return;
            }
        }
        std::string complaint = std::string(name) + " was " + value + " but must be one of ";
        for ( size_t i=0; i<setCount; ++i ) {
            complaint += (0<i) ? " , " : "";
            complaint += set[i];
        }
        throw complaint;
    }

    int strip_number_suffix(std::string &stripMe, int defaultValue) {
        size_t c = stripMe.length();
        while (0<c && '0'<=stripMe[c-1] && stripMe[c-1]<='9') {
            --c;
        }
        if (c==stripMe.length()) {
            return defaultValue;
        }
        int rv = convert_int( stripMe.c_str() + c );
        stripMe = stripMe.substr(0, c);
        return rv;
    }
};

void parseArg(int argc, char *argv[], Params &params) {
    int cnt;
    progress_display::setProgressDisplay(false);
    verbose_mode = VB_MIN;
    params.tree_gen = NONE;
    params.user_file = NULL;
    params.constraint_tree_file = NULL;
    params.opt_gammai = true;
    params.opt_gammai_fast = false;
    params.opt_gammai_keep_bran = false;
    params.testAlphaEpsAdaptive = false;
    params.randomAlpha = false;
    params.testAlphaEps = 0.1;
    params.exh_ai = false;
    params.alpha_invar_file = NULL;
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
    params.dist_file = nullptr;
    params.dist_format = "square";
    params.incremental = false;
    params.dist_compression_level = 1;
    params.compute_obs_dist = false;
    params.count_unknown_as_different = false;
    params.compute_jc_dist = true;
    params.experimental = true;
    params.compute_ml_dist = true;
    params.compute_ml_tree = true;
    params.compute_ml_tree_only = false;
    params.budget_file = NULL;
    params.overlap = 0;
    params.is_rooted = false;
    params.root_move_dist = 2;
    params.root_find = false;
    params.root_test = false;
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
    params.num_zero_len = 0;
    params.pd_limit = 100;
    params.calc_pdgain = false;
    params.multi_tree = false;
    params.second_tree = NULL;
    params.support_tag = NULL;
    params.site_concordance = 0;
    params.site_concordance_partition = false;
    params.print_cf_quartets = false;
    params.print_df1_trees = false;
    params.internode_certainty = 0;
    params.tree_weight_file = NULL;
    params.consensus_type = CT_NONE;
    params.find_pd_min = false;
    params.branch_cluster = 0;
    params.taxa_order_file = NULL;
    params.endemic_pd = false;
    params.exclusive_pd = false;
    params.complement_area = NULL;
    params.scaling_factor = -1;
    params.numeric_precision = -1;
    params.binary_programming = false;
    params.quad_programming = false;
    params.test_input = TEST_NONE;
    params.tree_burnin = 0;
    params.tree_max_count = 1000000;
    params.split_threshold = 0.0;
    params.split_threshold_str = NULL;
    params.split_weight_threshold = -1000;
    params.collapse_zero_branch = false;
    params.split_weight_summary = SW_SUM;
    params.gurobi_format = true;
    params.gurobi_threads = 1;
    params.num_bootstrap_samples = 0;
    params.bootstrap_spec = NULL;
    params.transfer_bootstrap = 0;

    params.aln_file = NULL;
    params.phylip_sequential_format = false;
    params.symtest = SYMTEST_NONE;
    params.symtest_only = false;
    params.symtest_remove = 0;
    params.symtest_keep_zero = false;
    params.symtest_type = 0;
    params.symtest_pcutoff = 0.05;
    params.symtest_stat = false;
    params.symtest_shuffle = 1;
    //params.treeset_file = NULL;
    params.topotest_replicates = 0;
    params.topotest_optimize_model = false;
    params.do_weighted_test = false;
    params.do_au_test = false;
    params.siteLL_file = NULL; //added by MA
    params.partition_file = NULL;
    params.partition_type = BRLEN_OPTIMIZE;
    params.partfinder_rcluster = 100;
    params.partfinder_rcluster_max = 0;
    params.partition_merge = MERGE_NONE;
    params.merge_models = "1";
    params.merge_rates = "1";
    params.partfinder_log_rate = true;
    params.remove_empty_seq = true;
    params.terrace_aware = true;
#ifdef IQTREE_TERRAPHAST
    params.terrace_analysis = false;
#else
    params.terrace_analysis = false;
#endif
    params.sequence_type = NULL;
    params.aln_output = NULL;
    params.aln_site_list = NULL;
    params.aln_output_format = IN_PHYLIP;
    params.output_format = FORMAT_NORMAL;
    params.newick_extended_format = false;
    params.gap_masked_aln = NULL;
    params.concatenate_aln = NULL;
    params.aln_nogaps = false;
    params.aln_no_const_sites = false;
    params.print_aln_info = false;
//    params.parsimony = false;
//    params.parsimony_tree = false;
    params.tree_spr = false;
    params.nexus_output = false;
    params.k_representative = 4;
    params.loglh_epsilon = 0.001;
    params.numSmoothTree = 1;
    params.nni5 = true;
    params.nni5_num_eval = 1;
    params.brlen_num_traversal = 1;
    params.leastSquareBranch = false;
    params.pars_branch_length = false;
    params.bayes_branch_length = false;
    params.manuel_analytic_approx = false;
    params.leastSquareNNI = false;
    params.ls_var_type = OLS;
    params.maxCandidates = 20;
    params.popSize = 5;
    params.p_delete = -1;
    params.min_iterations = -1;
    params.max_iterations = 1000;
    params.num_param_iterations = 100;
    params.stop_condition = SC_UNSUCCESS_ITERATION;
    params.stop_confidence = 0.95;
    params.num_runs = 1;
    params.model_name = "";
    params.model_name_init = NULL;
    params.model_opt_steps = 10;
    params.model_set = "ALL";
    params.model_extra_set = NULL;
    params.model_subset = NULL;
    params.state_freq_set = NULL;
    params.ratehet_set = "AUTO";
    params.score_diff_thres = 10.0;
    params.model_def_file = NULL;
    params.modelomatic = false;
    params.model_test_again = false;
    params.model_test_and_tree = 0;
    params.model_test_separate_rate = false;
    params.optimize_mixmodel_weight = false;
    params.optimize_rate_matrix = false;
    params.store_trans_matrix = false;
    //params.freq_type = FREQ_EMPIRICAL;
    params.freq_type = FREQ_UNKNOWN;
    params.keep_zero_freq = true;
    params.min_state_freq = MIN_FREQUENCY;
    params.min_rate_cats = 2;
    params.num_rate_cats = 4;
    params.max_rate_cats = 10;
    params.gamma_shape = -1.0;
    params.min_gamma_shape = MIN_GAMMA_SHAPE;
    params.gamma_median = false;
    params.p_invar_sites = -1.0;
    params.optimize_model_rate_joint = false;
    params.optimize_by_newton = true;
    params.optimize_alg_freerate = "2-BFGS,EM";
    params.optimize_alg_mixlen = "EM";
    params.optimize_alg_gammai = "EM";
    params.optimize_from_given_params = false;
    params.fixed_branch_length = BRLEN_OPTIMIZE;
    params.min_branch_length = 0.0; // this is now adjusted later based on alignment length
    // TODO DS: This seems inappropriate for PoMo.  It is handled in
    // phyloanalysis::2908.
    params.max_branch_length = 10.0; // Nov 22 2016: reduce from 100 to 10!
    params.iqp_assess_quartet = IQP_DISTANCE;
    params.iqp = false;
    params.write_intermediate_trees = 0;
//    params.avoid_duplicated_trees = false;
    params.writeDistImdTrees = false;
    params.rf_dist_mode = 0;
    params.rf_same_pair = false;
    params.normalize_tree_dist = false;
    params.mvh_site_rate = false;
    params.rate_mh_type = true;
    params.discard_saturated_site = false;
    params.mean_rate = 1.0;
    params.aLRT_threshold = 101;
    params.aLRT_replicates = 0;
    params.aLRT_test = false;
    params.aBayes_test = false;
    params.localbp_replicates = 0;
#ifdef __AVX512KNL
    params.SSE = LK_AVX512;
#else
    params.SSE = LK_AVX_FMA;
#endif
    params.lk_safe_scaling = false;
    params.numseq_safe_scaling = 2000;
    params.ignore_any_errors = false;
    params.kernel_nonrev = false;
    params.print_site_lh = WSL_NONE;
    params.print_partition_lh = false;
    params.print_site_prob = WSL_NONE;
    params.print_site_state_freq = WSF_NONE;
    params.print_site_rate = 0;
    params.print_trees_site_posterior = 0;
    params.print_ancestral_sequence = AST_NONE;
    params.min_ancestral_prob = 0.0;
    params.print_tree_lh = false;
    params.lambda = 1;
    params.speed_conf = 1.0;
    params.whtest_simulations = 1000;
    params.mcat_type = MCAT_LOG + MCAT_PATTERN;
    params.rate_file = NULL;
    params.ngs_file = NULL;
    params.ngs_mapped_reads = NULL;
    params.ngs_ignore_gaps = true;
    params.do_pars_multistate = false;
    params.gene_pvalue_file = NULL;
    params.gene_scale_factor = -1;
    params.gene_pvalue_loga = false;
    params.second_align = NULL;
    params.ncbi_taxid = 0;
    params.ncbi_taxon_level = NULL;
    params.ncbi_names_file = NULL;
    params.ncbi_ignore_level = NULL;

	params.eco_dag_file  = NULL;
	params.eco_type = NULL;
	params.eco_detail_file = NULL;
	params.k_percent = 0;
	params.diet_min = 0;
	params.diet_max = 0;
	params.diet_step = 0;
	params.eco_weighted = false;
	params.eco_run = 0;

	params.upper_bound = false;
	params.upper_bound_NNI = false;
	params.upper_bound_frac = 0.0;

    params.gbo_replicates = 0;
	params.ufboot_epsilon = 0.5;
    params.check_gbo_sample_size = 0;
    params.use_rell_method = true;
    params.use_elw_method = false;
    params.use_weighted_bootstrap = false;
    params.use_max_tree_per_bootstrap = true;
    params.max_candidate_trees = 0;
    params.distinct_trees = false;
    params.online_bootstrap = true;
    params.min_correlation = 0.99;
    params.step_iterations = 100;
//    params.store_candidate_trees = false;
	params.print_ufboot_trees = 0;
    params.jackknife_prop = 0.0;
    params.robust_phy_keep = 1.0;
    params.robust_median = false;
    //const double INF_NNI_CUTOFF = -1000000.0;
    params.nni_cutoff = -1000000.0;
    params.estimate_nni_cutoff = false;
    params.nni_sort = false;
    //params.nni_opt_5branches = false;
    params.testNNI = false;
    params.approximate_nni = false;
    params.do_compression = false;

    params.new_heuristic = true;
    params.iteration_multiple = 1;
    params.initPS = 0.5;
#ifdef USING_PLL
    params.pll = true;
#else
    params.pll = false;
#endif
    params.modelEps = 0.01;
    params.modelfinder_eps = 0.1;
    params.parbran = false;
    params.binary_aln_file = NULL;
    params.maxtime = 1000000;
    params.reinsert_par = false;
    params.bestStart = true;
    params.snni = true; // turn on sNNI default now
//    params.autostop = true; // turn on auto stopping rule by default now
    params.unsuccess_iteration = 100;
    params.speednni = true; // turn on reduced hill-climbing NNI by default now
    params.numInitTrees = 100;
    params.fixStableSplits = false;
    params.stableSplitThreshold = 0.9;
    params.five_plus_five = false;
    params.memCheck = false;
    params.tabu = false;
    params.adaptPertubation = false;
    params.numSupportTrees = 20;
//    params.sprDist = 20;
    params.sprDist = 6;
    params.sankoff_cost_file = NULL;
    params.numNNITrees = 20;
    params.avh_test = 0;
    params.bootlh_test = 0;
    params.bootlh_partitions = NULL;
    params.site_freq_file = NULL;
    params.tree_freq_file = NULL;
    params.num_threads = 1;
    params.num_threads_max = 10000;
    params.openmp_by_model = false;
    params.model_test_criterion = MTC_BIC;
//    params.model_test_stop_rule = MTC_ALL;
    params.model_test_sample_size = 0;
    params.root_state = NULL;
    params.print_bootaln = false;
    params.print_boot_site_freq = false;
	params.print_subaln = false;
	params.print_partition_info = false;
	params.print_conaln = false;
	params.count_trees = false;
    params.pomo = false;
    params.pomo_random_sampling = false;
	// params.pomo_counts_file_flag = false;
	params.pomo_pop_size = 9;
	params.print_branch_lengths = false;
	params.lh_mem_save = LM_PER_NODE; // auto detect
    params.buffer_mem_save = false;
	params.start_tree = STT_PLL_PARSIMONY;
    params.start_tree_subtype_name = StartTree::Factory::getNameOfDefaultTreeBuilder();

    params.modelfinder_ml_tree = true;
    params.final_model_opt = true;
	params.print_splits_file = false;
    params.print_splits_nex_file = true;
    params.ignore_identical_seqs = true;
    params.write_init_tree = false;
    params.write_candidate_trees = false;
    params.write_branches = false;
    params.freq_const_patterns = NULL;
    params.no_rescale_gamma_invar = false;
    params.compute_seq_identity_along_tree = false;
    params.compute_seq_composition = true;
    params.lmap_num_quartets = -1;
    params.lmap_cluster_file = NULL;
    params.print_lmap_quartet_lh = false;
    params.num_mixlen = 1;
    params.link_alpha = false;
    params.link_model = false;
    params.model_joint = NULL;
    params.ignore_checkpoint = false;
    params.checkpoint_dump_interval = 60;
    params.force_unfinished = false;
    params.suppress_output_flags = 0;
    params.ufboot2corr = false;
    params.u2c_nni5 = false;
    params.date_with_outgroup = true;
    params.date_debug = false;
    params.date_replicates = 0;
    params.clock_stddev = -1.0;
    params.date_outlier = -1.0;
    
    params.matrix_exp_technique = MET_EIGEN3LIB_DECOMPOSITION;

	if (params.nni5) {
	    params.nni_type = NNI5;
	} else {
	    params.nni_type = NNI1;
	}

    struct timeval tv;
    struct timezone tz;
    // initialize random seed based on current time
    gettimeofday(&tv, &tz);
    //params.ran_seed = (unsigned) (tv.tv_sec+tv.tv_usec);
    params.ran_seed = (tv.tv_usec);
    params.subsampling_seed = params.ran_seed;
    params.subsampling = 0;
    
    params.suppress_list_of_sequences = false;
    params.suppress_zero_distance_warnings = false;
    params.suppress_duplicate_sequence_warnings = false;

    for (cnt = 1; cnt < argc; cnt++) {
        try {

            if (strcmp(argv[cnt], "-h") == 0 || strcmp(argv[cnt], "--help") == 0) {
#ifdef IQ_TREE
                usage_iqtree(argv, strcmp(argv[cnt], "--help") == 0);
#else
                usage(argv, false);
#endif
                continue;
            }
            if (strcmp(argv[cnt], "-V") == 0 || strcmp(argv[cnt], "-version") == 0 || strcmp(argv[cnt], "--version") == 0) {
                printCopyright(cout);
                exit(EXIT_SUCCESS);
            }
			if (strcmp(argv[cnt], "-ho") == 0 || strcmp(argv[cnt], "-?") == 0) {
				usage_iqtree(argv, false);
				continue;
			}
			if (strcmp(argv[cnt], "-hh") == 0 || strcmp(argv[cnt], "-hhh") == 0) {
#ifdef IQ_TREE
                usage_iqtree(argv, true);
#else
				usage(argv);
#endif
				continue;
			}
			if (strcmp(argv[cnt], "-v0") == 0) {
				verbose_mode = VB_QUIET;
				continue;
			}
			if (strcmp(argv[cnt], "-v") == 0 || strcmp(argv[cnt], "--verbose") == 0) {
				verbose_mode = VB_MED;
				continue;
			}
			if (strcmp(argv[cnt], "-vv") == 0
					|| strcmp(argv[cnt], "-v2") == 0) {
				verbose_mode = VB_MAX;
				continue;
			}
			if (strcmp(argv[cnt], "-vvv") == 0
					|| strcmp(argv[cnt], "-v3") == 0) {
				verbose_mode = VB_DEBUG;
				continue;
			}
			if (strcmp(argv[cnt], "-k") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -k <num_taxa>";
				convert_range(argv[cnt], params.min_size, params.sub_size,
						params.step_size);
				params.k_representative = params.min_size;
				continue;
			}
			if (strcmp(argv[cnt], "-pre") == 0 || strcmp(argv[cnt], "--prefix") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pre <output_prefix>";
				params.out_prefix = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-pp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pp <pd_proportion>";
				convert_range(argv[cnt], params.min_proportion,
						params.pd_proportion, params.step_proportion);
				if (params.pd_proportion < 0 || params.pd_proportion > 1)
					throw "PD proportion must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-mk") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mk <min_taxa>";
				params.min_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bud") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bud <budget>";
				convert_range(argv[cnt], params.min_budget, params.budget,
						params.step_budget);
				continue;
			}
			if (strcmp(argv[cnt], "-mb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mb <min_budget>";
				params.min_budget = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-o") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -o <taxon>";
				params.root = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-optalg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -opt_alg <1-BFGS|2-BFGS|EM>";
				params.optimize_alg_freerate = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-optlen") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -optlen <BFGS|EM>";
				params.optimize_alg_mixlen = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "-optalg_gammai") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -optalg_gammai <Brent|BFGS|EM>";
                params.optimize_alg_gammai = argv[cnt];
                continue;
            }
			if (strcmp(argv[cnt], "-root") == 0 || strcmp(argv[cnt], "-rooted") == 0) {
				params.is_rooted = true;
				continue;
			}
            
            if (strcmp(argv[cnt], "--root-dist") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --root-dist <maximum-root-move-distance>";
                params.root_move_dist = convert_int(argv[cnt]);
                continue;
            }

            if (strcmp(argv[cnt], "--root-find") == 0) {
                params.root_find = true;
                continue;
            }

            if (strcmp(argv[cnt], "--root-test") == 0) {
                params.root_test = true;
                continue;
            }
            
			if (strcmp(argv[cnt], "-all") == 0) {
				params.find_all = true;
				continue;
			}
			if (strcmp(argv[cnt], "--greedy") == 0) {
				params.run_mode = GREEDY;
				continue;
			}
			if (strcmp(argv[cnt], "-pr") == 0
					|| strcmp(argv[cnt], "--pruning") == 0) {
				params.run_mode = PRUNING;
				//continue; } if (strcmp(argv[cnt],"--both") == 0) {
				//params.run_mode = BOTH_ALG;
				continue;
			}
			if (strcmp(argv[cnt], "-e") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -e <file>";
				params.param_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-if") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -if <file>";
				params.initial_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nni_nr_step") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nni_nr_step <newton_raphson_steps>";
				NNI_MAX_NR_STEP = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-ia") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ia <file>";
				params.initial_area_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-u") == 0) {
				// file containing budget information
				cnt++;
				if (cnt >= argc)
					throw "Use -u <file>";
				params.budget_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dd") == 0) {
				// compute distribution of PD score on random sets
				cnt++;
				if (cnt >= argc)
					throw "Use -dd <sample_size>";
				params.run_mode = PD_DISTRIBUTION;
				params.sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-ts") == 0) {
				// calculate PD score a taxa set listed in the file
				cnt++;
				//params.run_mode = PD_USER_SET;
				if (cnt >= argc)
					throw "Use -ts <taxa_file>";
				params.pdtaxa_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-bound") == 0) {
				// boundary length of areas
				cnt++;
				if (cnt >= argc)
					throw "Use -bound <file>";
				params.areas_boundary_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-blm") == 0) {
				// boundary length modifier
				cnt++;
				if (cnt >= argc)
					throw "Use -blm <boundary_modifier>";
				params.boundary_modifier = convert_double(argv[cnt]);
				continue;
			}
            std::string arg = argv[cnt];
            //Todo; move this up, use == rather than strcmp elsewhere, too.
            if (arg=="-dist-format") {
                std::string nextArg = next_argument(argc, argv, "<distance_file_format>", cnt);
                params.dist_format = string_to_lower(nextArg.c_str());
                params.dist_compression_level
                    = strip_number_suffix(params.dist_format,
                                          params.dist_compression_level);
                if (params.dist_compression_level<0) {
                    params.dist_compression_level=0;
                } else if (9<params.dist_compression_level) {
                    params.dist_compression_level=9;
                }
                const char* allowed[] = {
                    "square", "lower", "upper"
                    , "square.gz", "lower.gz", "upper.gz"
                };
                throw_if_not_in_set ( "dist-format", params.dist_format
                                    , allowed, sizeof(allowed)/sizeof(allowed[0]));
                continue;
            }
            if (arg=="-update") {
                params.incremental = true;
                std::string method = next_argument(argc, argv, "incremental method", cnt);
                params.incremental_method = string_to_upper(method.c_str());
                continue;
            }
            if (arg=="-merge") {
                std::string alignment_file = next_argument(argc, argv, "alignment file (to merge)", cnt);
                params.additional_alignment_files.emplace_back(alignment_file);            
                continue;
            }
            
            if (strcmp(argv[cnt], "-dist") == 0
                || strcmp(argv[cnt], "-d") == 0) {
                // calculate distance matrix from the tree
                params.run_mode = CALC_DIST;
                cnt++;
                if (cnt >= argc)
                    throw "Use -dist <distance_file>";
                params.dist_file = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "-djc") == 0) {
                params.compute_ml_dist = false;
                continue;
            }
            if (arg=="-mlnj-only" || arg=="--mlnj-only") {
                params.compute_ml_tree_only = true;
                continue;
            }
            if (strcmp(argv[cnt], "-dobs") == 0) {
                params.compute_obs_dist = true;
                continue;
            }
            if (strcmp(argv[cnt], "-cud") == 0) {
                params.count_unknown_as_different = true;
                continue;
            }
            if (strcmp(argv[cnt], "-experimental") == 0 || strcmp(argv[cnt], "--experimental") == 0) {
                params.experimental = true;
                continue;
            }
            if (strcmp(argv[cnt], "--no-experimental") == 0) {
                params.experimental = false;
                continue;
            }
			if (strcmp(argv[cnt], "-r") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -r <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = YULE_HARDING;
				continue;
			}
			if (strcmp(argv[cnt], "-rs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rs <alignment_file>";
				params.tree_gen = YULE_HARDING;
				params.aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-rstar") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rstar <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = STAR_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-ru") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ru <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = UNIFORM;
				continue;
			}
			if (strcmp(argv[cnt], "-rcat") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcat <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CATERPILLAR;
				continue;
			}
			if (strcmp(argv[cnt], "-rbal") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rbal <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = BALANCED;
				continue;
			}
            if (strcmp(argv[cnt], "--rand") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -rand UNI | CAT | BAL | NET";
                if (strcmp(argv[cnt], "UNI") == 0)
                    params.tree_gen = UNIFORM;
                else if (strcmp(argv[cnt], "CAT") == 0)
                    params.tree_gen = CATERPILLAR;
                else if (strcmp(argv[cnt], "BAL") == 0)
                    params.tree_gen = BALANCED;
                else if (strcmp(argv[cnt], "NET") == 0)
                    params.tree_gen = CIRCULAR_SPLIT_GRAPH;
                else
                    throw "wrong --rand option";
                continue;
            }
            
            if (strcmp(argv[cnt], "--keep-ident") == 0 || strcmp(argv[cnt], "-keep-ident") == 0) {
                params.ignore_identical_seqs = false;
                continue;
            }
			if (strcmp(argv[cnt], "-rcsg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcsg <num_taxa>";
				params.sub_size = convert_int(argv[cnt]);
				params.tree_gen = CIRCULAR_SPLIT_GRAPH;
				continue;
			}
			if (strcmp(argv[cnt], "-rpam") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rpam <num_splits>";
				params.num_splits = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-rlen") == 0 || strcmp(argv[cnt], "--rlen") == 0) {
				cnt++;
				if (cnt >= argc - 2)
					throw "Use -rlen <min_len> <mean_len> <max_len>";
				params.min_len = convert_double(argv[cnt]);
				params.mean_len = convert_double(argv[cnt + 1]);
				params.max_len = convert_double(argv[cnt + 2]);
				cnt += 2;
				continue;
			}
			if (strcmp(argv[cnt], "-rzero") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rzero <num_zero_branch>";
				params.num_zero_len = convert_int(argv[cnt]);
				if (params.num_zero_len < 0)
					throw "num_zero_len must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-rset") == 0) {
				cnt++;
				if (cnt >= argc - 1)
					throw "Use -rset <overlap> <outfile>";
				params.overlap = convert_int(argv[cnt]);
				cnt++;
				params.pdtaxa_file = argv[cnt];
				params.tree_gen = TAXA_SET;
				continue;
			}
			if (strcmp(argv[cnt], "-rep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rep <repeated_times>";
				params.repeated_time = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lim") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lim <pd_limit>";
				params.pd_limit = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-noout") == 0) {
				params.nr_output = 0;
				continue;
			}
			if (strcmp(argv[cnt], "-1out") == 0) {
				params.nr_output = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-oldout") == 0) {
				params.nr_output = 100;
				continue;
			}
			if (strcmp(argv[cnt], "-nexout") == 0) {
				params.nexus_output = true;
				continue;
			}
			if (strcmp(argv[cnt], "-exhaust") == 0) {
				params.run_mode = EXHAUSTIVE;
				continue;
			}
			if (strcmp(argv[cnt], "-seed") == 0 || strcmp(argv[cnt], "--seed") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -seed <random_seed>";
				params.ran_seed = abs(convert_int(argv[cnt]));
				continue;
			}
			if (strcmp(argv[cnt], "-pdgain") == 0) {
				params.calc_pdgain = true;
				continue;
			}
			if (strcmp(argv[cnt], "-sup") == 0 || strcmp(argv[cnt], "--support") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sup <target_tree_file>";
				params.second_tree = argv[cnt];
				params.consensus_type = CT_ASSIGN_SUPPORT;
				continue;
			}
			if (strcmp(argv[cnt], "-suptag") == 0 || strcmp(argv[cnt], "--suptag") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -suptag <tagname or ALL>";
				params.support_tag = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "-sup2") == 0) {
                outError("Deprecated -sup2 option, please use --gcf --tree FILE");
            }
            
            if (strcmp(argv[cnt], "--gcf") == 0) {
				params.consensus_type = CT_ASSIGN_SUPPORT_EXTENDED;
                cnt++;
                if (cnt >= argc)
                    throw "Use --gcf <user_trees_file>";
                params.treeset_file = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "--scf") == 0) {
                params.consensus_type = CT_ASSIGN_SUPPORT_EXTENDED;
                cnt++;
                if (cnt >= argc)
                    throw "Use --scf NUM_QUARTETS";
                params.site_concordance = convert_int(argv[cnt]);
                if (params.site_concordance < 1)
                    throw "Positive --scf please";
                continue;
            }
            if (strcmp(argv[cnt], "--scf-part") == 0 || strcmp(argv[cnt], "--cf-verbose") == 0) {
                params.site_concordance_partition = true;
                continue;
            }
            if (strcmp(argv[cnt], "--cf-quartet") == 0) {
                params.print_cf_quartets = true;
                continue;
            }
            if (strcmp(argv[cnt], "--df-tree") == 0) {
                params.print_df1_trees = true;
                continue;
            }
            if (strcmp(argv[cnt], "--qic") == 0) {
                params.internode_certainty = 1;
                continue;
            }
			if (strcmp(argv[cnt], "-treew") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -treew <tree_weight_file>";
				params.tree_weight_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-con") == 0 || strcmp(argv[cnt], "--contree") == 0) {
				params.consensus_type = CT_CONSENSUS_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-net") == 0 || strcmp(argv[cnt], "--connet") == 0) {
				params.consensus_type = CT_CONSENSUS_NETWORK;
                continue;
			}
            
            /**MINH ANH: to serve some statistics on tree*/
			if (strcmp(argv[cnt], "-comp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -comp <treefile>";
				params.consensus_type = COMPARE;
				params.second_tree = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-stats") == 0) {
				params.run_mode = STATS;
				continue;
			}
			if (strcmp(argv[cnt], "-gbo") == 0) { //guided bootstrap
				cnt++;
				if (cnt >= argc)
					throw "Use -gbo <site likelihod file>";
				params.siteLL_file = argv[cnt];
				//params.run_mode = GBO;
                continue;
			} // MA
            
			if (strcmp(argv[cnt], "-mprob") == 0) { //compute multinomial distribution probability
				cnt++;
				if (cnt >= argc)
					throw "Use -mprob <ref_alignment>";
				params.second_align = argv[cnt];
				//params.run_mode = MPRO;
                continue;
			} // MA
            
			if (strcmp(argv[cnt], "-min") == 0) {
				params.find_pd_min = true;
				continue;
			}
			if (strcmp(argv[cnt], "-excl") == 0) {
				params.exclusive_pd = true;
				continue;
			}
			if (strcmp(argv[cnt], "-endem") == 0) {
				params.endemic_pd = true;
				continue;
			}
			if (strcmp(argv[cnt], "-compl") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -compl <area_name>";
				params.complement_area = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-cluster") == 0) {
				params.branch_cluster = 4;
				cnt++;
				if (cnt >= argc)
					throw "Use -cluster <taxa_order_file>";
				params.taxa_order_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-taxa") == 0) {
				params.run_mode = PRINT_TAXA;
				continue;
			}
			if (strcmp(argv[cnt], "-area") == 0) {
				params.run_mode = PRINT_AREA;
				continue;
			}
			if (strcmp(argv[cnt], "-scale") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -scale <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scaleg") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -scaleg <gene_scale_factor>";
				params.gene_scale_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scalebranch") == 0) {
				params.run_mode = SCALE_BRANCH_LEN;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalebranch <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-scalenode") == 0) {
				params.run_mode = SCALE_NODE_NAME;
				cnt++;
				if (cnt >= argc)
					throw "Use -scalenode <scaling_factor>";
				params.scaling_factor = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-prec") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -prec <numeric_precision>";
				params.numeric_precision = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lp") == 0) {
				params.run_mode = LINEAR_PROGRAMMING;
				continue;
			}
			if (strcmp(argv[cnt], "-lpbin") == 0) {
				params.run_mode = LINEAR_PROGRAMMING;
				params.binary_programming = true;
				continue;
			}
			if (strcmp(argv[cnt], "-qp") == 0) {
				params.gurobi_format = true;
				params.quad_programming = true;
				continue;
			}
			if (strcmp(argv[cnt], "-quiet") == 0 || strcmp(argv[cnt], "--quiet") == 0) {
				verbose_mode = VB_QUIET;
				continue;
			}
			if (strcmp(argv[cnt], "-mult") == 0) {
				params.multi_tree = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bi") == 0 || strcmp(argv[cnt], "--burnin") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bi <burnin_value>";
				params.tree_burnin = convert_int(argv[cnt]);
				if (params.tree_burnin < 0)
					throw "Burnin value must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-tm") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -tm <tree_max_count>";
				params.tree_max_count = convert_int(argv[cnt]);
				if (params.tree_max_count < 0)
					throw "tree_max_count must not be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-minsup") == 0 || strcmp(argv[cnt], "--sup-min") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -minsup <split_threshold>";
				params.split_threshold = convert_double(argv[cnt]);
				if (params.split_threshold < 0 || params.split_threshold > 1)
					throw "Split threshold must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-minsupnew") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -minsupnew <split_threshold_1/.../split_threshold_k>";
				params.split_threshold_str = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-tw") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -tw <split_weight_threshold>";
				params.split_weight_threshold = convert_double(argv[cnt]);
				if (params.split_weight_threshold < 0)
					throw "Split weight threshold is negative";
				continue;
			}

            if (strcmp(argv[cnt], "-czb") == 0 || strcmp(argv[cnt], "--polytomy") == 0) {
                params.collapse_zero_branch = true;
                continue;
            }

			if (strcmp(argv[cnt], "-swc") == 0) {
				params.split_weight_summary = SW_COUNT;
				continue;
			}
			if (strcmp(argv[cnt], "-swa") == 0) {
				params.split_weight_summary = SW_AVG_ALL;
				continue;
			}
			if (strcmp(argv[cnt], "-swp") == 0) {
				params.split_weight_summary = SW_AVG_PRESENT;
				continue;
			}
			if (strcmp(argv[cnt], "-iwc") == 0) {
				params.test_input = TEST_WEAKLY_COMPATIBLE;
				continue;
			}
			if (strcmp(argv[cnt], "--aln") == 0 || strcmp(argv[cnt], "--msa") == 0 || strcmp(argv[cnt], "-s") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use --aln, -s <alignment_file>";
				params.aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "--sequential") == 0) {
                params.phylip_sequential_format = true;
                continue;
            }
            if (strcmp(argv[cnt], "--symtest") == 0) {
                params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--bisymtest") == 0) {
                params.symtest = SYMTEST_BINOM;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-only") == 0) {
                params.symtest_only = true;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-remove-bad") == 0) {
                params.symtest_remove = 1;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-remove-good") == 0) {
                params.symtest_remove = 2;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-keep-zero") == 0) {
                params.symtest_keep_zero = true;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-type") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --symtest-type SYM|MAR|INT";
                if (strcmp(argv[cnt], "SYM") == 0)
                    params.symtest_type = 0;
                else if (strcmp(argv[cnt], "MAR") == 0)
                    params.symtest_type = 1;
                else if (strcmp(argv[cnt], "INT") == 0)
                    params.symtest_type = 2;
                else
                    throw "Use --symtest-type SYM|MAR|INT";
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-pval") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --symtest-pval PVALUE_CUTOFF";
                params.symtest_pcutoff = convert_double(argv[cnt]);
                if (params.symtest_pcutoff <= 0 || params.symtest_pcutoff >= 1)
                    throw "--symtest-pval must be between 0 and 1";
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }
            
            if (strcmp(argv[cnt], "--symstat") == 0) {
                params.symtest_stat = true;
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "--symtest-perm") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --symtest-perm INT";
                params.symtest_shuffle = convert_int(argv[cnt]);
                if (params.symtest_shuffle <= 0)
                    throw "--symtest-perm must be positive";
                if (params.symtest == SYMTEST_NONE)
                    params.symtest = SYMTEST_MAXDIV;
                continue;
            }

            if (strcmp(argv[cnt], "-z") == 0 || strcmp(argv[cnt], "--trees") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -z <user_trees_file>";
				params.treeset_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-zb") == 0 || strcmp(argv[cnt], "--test") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -zb <#replicates>";
				params.topotest_replicates = convert_int(argv[cnt]);
				if (params.topotest_replicates < 1000)
					throw "Please specify at least 1000 replicates";
				continue;
			}
            if (strcmp(argv[cnt], "--estimate-model") == 0) {
                params.topotest_optimize_model = true;
                continue;
            }
			if (strcmp(argv[cnt], "-zw") == 0 || strcmp(argv[cnt], "--test-weight") == 0) {
				params.do_weighted_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-au") == 0 || strcmp(argv[cnt], "--test-au") == 0) {
				params.do_au_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-sp") == 0 || strcmp(argv[cnt], "-Q") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sp <partition_file>";
				params.partition_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-spp") == 0 || strcmp(argv[cnt], "-p") == 0 || strcmp(argv[cnt], "--partition") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -p <partition_file>";
				params.partition_file = argv[cnt];
				params.partition_type = BRLEN_SCALE;
                params.opt_gammai = false;
				continue;
			}
			if (strcmp(argv[cnt], "-spj") == 0 || strcmp(argv[cnt], "-q") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -q <partition_file>";
				params.partition_file = argv[cnt];
				params.partition_type = BRLEN_FIX;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
				continue;
			}
			if (strcmp(argv[cnt], "-M") == 0) {
                params.partition_type = BRLEN_OPTIMIZE;
                continue;
            }

            if (strcmp(argv[cnt], "-spu") == 0 || strcmp(argv[cnt], "-S") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -spu <partition_file>";
                params.partition_file = argv[cnt];
                params.partition_type = TOPO_UNLINKED;
                params.ignore_identical_seqs = false;
                params.buffer_mem_save = true;
                params.print_splits_nex_file = false;
                continue;
            }
            
            if (strcmp(argv[cnt], "--edge") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --edge equal|scale|unlink";
                if (strcmp(argv[cnt], "equal") == 0)
                    params.partition_type = BRLEN_FIX;
                else if (strcmp(argv[cnt], "scale") == 0)
                    params.partition_type = BRLEN_SCALE;
                else if (strcmp(argv[cnt], "unlink") == 0)
                    params.partition_type = BRLEN_OPTIMIZE;
                else
                    throw "Use --edge equal|scale|unlink";
            }
            
            if (strcmp(argv[cnt], "-rcluster") == 0 || strcmp(argv[cnt], "--rcluster") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcluster <percent>";
                params.partfinder_rcluster = convert_double(argv[cnt]);
                if (params.partfinder_rcluster < 0 || params.partfinder_rcluster > 100)
                    throw "rcluster percentage must be between 0 and 100";
                params.partition_merge = MERGE_RCLUSTER;
				continue;
            }
            if (strcmp(argv[cnt], "-rclusterf") == 0 || strcmp(argv[cnt], "--rclusterf") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rclusterf <percent>";
                params.partfinder_rcluster = convert_double(argv[cnt]);
                if (params.partfinder_rcluster < 0 || params.partfinder_rcluster > 100)
                    throw "rcluster percentage must be between 0 and 100";
                params.partition_merge = MERGE_RCLUSTERF;
				continue;
            }

            if (strcmp(argv[cnt], "-rcluster-max") == 0 || strcmp(argv[cnt], "--rcluster-max") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -rcluster-max <num>";
                params.partfinder_rcluster_max = convert_int(argv[cnt]);
                if (params.partfinder_rcluster_max <= 0)
                    throw "rcluster-max must be between > 0";
                if (params.partfinder_rcluster == 100)
                    params.partfinder_rcluster = 99.9999;
                if (params.partition_merge != MERGE_RCLUSTER && params.partition_merge != MERGE_RCLUSTERF)
                    params.partition_merge = MERGE_RCLUSTERF;
				continue;
            }

            if (strcmp(argv[cnt], "--merge") == 0) {
                if (cnt >= argc-1 || argv[cnt+1][0] == '-') {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                    continue;
                }
                cnt++;
                if (cnt >= argc)
                    throw "Use --merge [none|greedy|rcluster|rclusterf|kmeans]";
                if (strcmp(argv[cnt], "none") == 0)
                    params.partition_merge = MERGE_NONE;
                else if (strcmp(argv[cnt], "greedy") == 0)
                    params.partition_merge = MERGE_GREEDY;
                else if (strcmp(argv[cnt], "rcluster") == 0) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTER;
                } else if (strcmp(argv[cnt], "rclusterf") == 0) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                } else if (strcmp(argv[cnt], "rcluster") == 0)
                    params.partition_merge = MERGE_KMEANS;
                else
                    throw "Use --merge [none|greedy|rcluster|rclusterf|kmeans]";
                continue;
            }

            if (strcmp(argv[cnt], "--merge-model") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --merge-model 1|4|ALL|model1,...,modelK";
                params.merge_models = argv[cnt];
                if (params.partition_merge == MERGE_NONE) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                    continue;
                }
                continue;
            }

            if (strcmp(argv[cnt], "--merge-rate") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --merge-rate rate1,...,rateK";
                params.merge_rates = argv[cnt];
                if (params.partition_merge == MERGE_NONE) {
                    if (params.partfinder_rcluster == 100)
                        params.partfinder_rcluster = 99.9999;
                    params.partition_merge = MERGE_RCLUSTERF;
                    continue;
                }
                continue;
            }

            if (strcmp(argv[cnt], "--merge-log-rate") == 0) {
                params.partfinder_log_rate = true;
                continue;
            }

            if (strcmp(argv[cnt], "--merge-normal-rate") == 0) {
                params.partfinder_log_rate = false;
                continue;
            }

			if (strcmp(argv[cnt], "-keep_empty_seq") == 0) {
				params.remove_empty_seq = false;
				continue;
			}
			if (strcmp(argv[cnt], "-no_terrace") == 0) {
				params.terrace_aware = false;
                params.terrace_analysis = false;
				continue;
			}
            if (strcmp(argv[cnt], "--terrace") == 0) {
#ifdef IQTREE_TERRAPHAST
                params.terrace_analysis = true;
#else
                    throw "Unsupported command: --terrace.\n"
                        "Please build IQ-TREE with the USE_TERRAPHAST flag.";
#endif
                continue;
            }

            if (strcmp(argv[cnt], "--no-terrace") == 0) {
                params.terrace_analysis = false;
                continue;
            }
            
            if (strcmp(argv[cnt], "-sf") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sf <ngs_file>";
				params.ngs_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-sm") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sm <ngs_mapped_read_file>";
				params.ngs_mapped_reads = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-ngs_gap") == 0) {
				params.ngs_ignore_gaps = false;
				continue;
			}
			if (strcmp(argv[cnt], "-st") == 0 || strcmp(argv[cnt], "--seqtype") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -st BIN or -st DNA or -st AA or -st CODON or -st MORPH or -st CRXX or -st CFxx.";
                params.sequence_type = argv[cnt];
                // if (arg.substr(0,2) == "CR") params.pomo_random_sampling = true;
                // if (arg.substr(0,2) == "CF" || arg.substr(0,2) == "CR") {
                //     outWarning("Setting the sampling method and population size with this flag is deprecated.");
                //     outWarning("Please use the model string instead (see `iqtree --help`).");
                //     if (arg.length() > 2) {
                //         int ps = convert_int(arg.substr(2).c_str());
                //         params.pomo_pop_size = ps;
                //         if (((ps != 10) && (ps != 2) && (ps % 2 == 0)) || (ps < 2) || (ps > 19)) {
                //             std::cout << "Please give a correct PoMo sequence type parameter; e.g., `-st CF09`." << std::endl;
                //             outError("Custom virtual population size of PoMo not 2, 10 or any other odd number between 3 and 19.");   
                //         }
                //     }
                // }
				continue;
			}
            
			if (strcmp(argv[cnt], "-starttree") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -starttree BIONJ|PARS|PLLPARS|PJ";
                else if (strcmp(argv[cnt], "PARS") == 0)
					params.start_tree = STT_PARSIMONY;
                else if (strcmp(argv[cnt], "PJ") == 0)
                    params.start_tree = STT_PARSIMONY_JOINING;
				else if (strcmp(argv[cnt], "PLLPARS") == 0)
					params.start_tree = STT_PLL_PARSIMONY;
                else if (START_TREE_RECOGNIZED(argv[cnt])) {
                    params.start_tree_subtype_name = argv[cnt];
                    params.start_tree = STT_BIONJ;
                }
                else
					throw "Invalid option, please use -starttree with BIONJ or PARS or PLLPARS";
				continue;
			}

			if (strcmp(argv[cnt], "-ao") == 0 || strcmp(argv[cnt], "--out-alignment") == 0 || strcmp(argv[cnt], "--out-aln") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ao <alignment_file>";
				params.aln_output = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-as") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -as <aln_site_list>";
				params.aln_site_list = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-an") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -an <ref_seq_name>";
				params.ref_seq_name = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-af") == 0 || strcmp(argv[cnt], "--out-format") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -af phy|fasta";
				if (strcmp(argv[cnt], "phy") == 0)
					params.aln_output_format = IN_PHYLIP;
				else if (strcmp(argv[cnt], "fasta") == 0)
					params.aln_output_format = IN_FASTA;
                else if (strcmp(argv[cnt], "nexus") == 0)
                    params.aln_output_format = IN_NEXUS;
				else
					throw "Unknown output format";
				continue;
			}


            if (strcmp(argv[cnt], "--out-csv") == 0) {
                params.output_format = FORMAT_CSV;
                continue;
            }
            
            if (strcmp(argv[cnt], "--out-tsv") == 0) {
                params.output_format = FORMAT_TSV;
                continue;
            }            

            if (strcmp(argv[cnt], "--figtree") == 0) {
                params.newick_extended_format = true;
                continue;
            }

            if (strcmp(argv[cnt], "-am") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -am <gap_masked_aln>";
				params.gap_masked_aln = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-ac") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ac <concatenate_aln>";
				params.concatenate_aln = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nogap") == 0) {
				params.aln_nogaps = true;
				continue;
			}
			if (strcmp(argv[cnt], "-noconst") == 0) {
				params.aln_no_const_sites = true;
				continue;
			}
			if (strcmp(argv[cnt], "-alninfo") == 0) {
				params.print_aln_info = true;
				continue;
			}
//			if (strcmp(argv[cnt], "-parstree") == 0) {
				// maximum parsimony
//				params.parsimony_tree = true;
//            continue; } if (strcmp(argv[cnt], "-pars") == 0) {
//                // maximum parsimony
//                params.parsimony = true;
//				continue;
//			}
			if (strcmp(argv[cnt], "-spr") == 0) {
				// subtree pruning and regrafting
				params.tree_spr = true;
				continue;
			}
			if (strcmp(argv[cnt], "-krep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -krep <num_k>";
				params.k_representative = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pdel") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pdel <probability>";
				params.p_delete = convert_double(argv[cnt]);
				if (params.p_delete < 0.0 || params.p_delete > 1.0)
					throw "Probability of deleting a leaf must be between 0 and 1";
				continue;
			}
			if (strcmp(argv[cnt], "-pers") == 0 || strcmp(argv[cnt], "--perturb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pers <perturbation_strength>";
				params.initPS = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-n") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -n <#iterations>";
                if (params.gbo_replicates != 0) {
                    throw("Ultrafast bootstrap does not work with -n option");
                }
				params.min_iterations = convert_int(argv[cnt]);
				params.stop_condition = SC_FIXED_ITERATION;
//                params.autostop = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nparam") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nparam <#iterations>";
				params.num_param_iterations = convert_int(argv[cnt]);
				if (params.num_param_iterations < 0)
					throw "Number of parameter optimization iterations (-nparam) must be non negative";
				continue;
			}

			if (strcmp(argv[cnt], "-nb") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nb <#bootstrap_replicates>";
				params.min_iterations = convert_int(argv[cnt]);
				params.iqp_assess_quartet = IQP_BOOTSTRAP;
//				params.avoid_duplicated_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "--model") == 0 || strcmp(argv[cnt], "-m") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use --model <model_name>";
				params.model_name = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "--init-model") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --init-model FILE";
                params.model_name_init = argv[cnt];
                continue;
            }
            if (strcmp(argv[cnt], "--loop-model") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --loop-model NUM";
                params.model_opt_steps = convert_int(argv[cnt]);
                continue;
            }
			if (strcmp(argv[cnt], "-mset") == 0 || strcmp(argv[cnt], "--mset") == 0 || strcmp(argv[cnt], "--models") == 0 ) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mset <model_set>";
				params.model_set = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-madd") == 0 || strcmp(argv[cnt], "--madd") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -madd <extra_model_set>";
				params.model_extra_set = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-msub") == 0 || strcmp(argv[cnt], "--msub") == 0 || strcmp(argv[cnt], "--model-sub") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -msub <model_subset>";
				params.model_subset = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-mfreq") == 0 || strcmp(argv[cnt], "--freqs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mfreq <state_freq_set>";
				params.state_freq_set = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-mrate") == 0 || strcmp(argv[cnt], "--mrate") == 0 || strcmp(argv[cnt], "--rates") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mrate <rate_set>";
				params.ratehet_set = argv[cnt];
				continue;
			}
            
            if (strcmp(argv[cnt], "--score-diff") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --score-diff <score>";
                if (iEquals(argv[cnt], "all"))
                    params.score_diff_thres = -1.0;
                else
                    params.score_diff_thres = convert_double(argv[cnt]);
                continue;
            }
            
			if (strcmp(argv[cnt], "-mdef") == 0 || strcmp(argv[cnt], "--mdef") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mdef <model_definition_file>";
				params.model_def_file = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "--modelomatic") == 0) {
                params.modelomatic = true;
                continue;
            }
			if (strcmp(argv[cnt], "-mredo") == 0 || strcmp(argv[cnt], "--mredo") == 0 || strcmp(argv[cnt], "--model-redo") == 0) {
				params.model_test_again = true;
				continue;
			}
			if (strcmp(argv[cnt], "-mtree") == 0 || strcmp(argv[cnt], "--mtree") == 0) {
				params.model_test_and_tree = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-mretree") == 0) {
				params.model_test_and_tree = 2;
				continue;
			}
			if (strcmp(argv[cnt], "-msep") == 0) {
				params.model_test_separate_rate = true;
				continue;
			}
			if (strcmp(argv[cnt], "-mwopt") == 0 || strcmp(argv[cnt], "--mix-opt") == 0) {
				params.optimize_mixmodel_weight = true;
				continue;
			}
			if (strcmp(argv[cnt], "--opt-rate-mat") == 0) {
				params.optimize_rate_matrix = true;
				continue;
			}
//			if (strcmp(argv[cnt], "-mh") == 0) {
//				params.mvh_site_rate = true;
//				params.discard_saturated_site = false;
//				params.SSE = LK_NORMAL;
//				continue;
//			}
//			if (strcmp(argv[cnt], "-mhs") == 0) {
//				params.mvh_site_rate = true;
//				params.discard_saturated_site = true;
//				params.SSE = LK_NORMAL;
//				continue;
//			}
			if (strcmp(argv[cnt], "-rl") == 0) {
				params.rate_mh_type = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nr") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nr <mean_rate>";
				params.mean_rate = convert_double(argv[cnt]);
				if (params.mean_rate < 0)
					throw "Wrong mean rate for MH model";
				continue;
			}
			if (strcmp(argv[cnt], "-mstore") == 0) {
				params.store_trans_matrix = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nni_lh") == 0) {
				params.nni_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lmd") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lmd <lambda>";
                params.lambda = convert_double(argv[cnt]);
				if (params.lambda > 1.0)
					throw "Lambda must be in (0,1]";
				continue;
			}

            if (strcmp(argv[cnt], "-lk") == 0) {
				cnt++;
				if (cnt >= argc)
                    throw "-lk x86|SSE|AVX|FMA|AVX512";
                if (strcmp(argv[cnt], "x86") == 0)
                    params.SSE = LK_386;
                else if (strcmp(argv[cnt], "SSE") == 0)
                    params.SSE = LK_SSE2;
                else if (strcmp(argv[cnt], "AVX") == 0)
                    params.SSE = LK_AVX;
                else if (strcmp(argv[cnt], "FMA") == 0)
                    params.SSE = LK_AVX_FMA;
                else if (strcmp(argv[cnt], "AVX512") == 0)
                    params.SSE = LK_AVX512;
                else
                    throw "Incorrect -lk likelihood kernel option";
				continue;
			}

			if (strcmp(argv[cnt], "-safe") == 0 || strcmp(argv[cnt], "--safe") == 0) {
				params.lk_safe_scaling = true;
				continue;
			}

			if (strcmp(argv[cnt], "-safe-seq") == 0) {
				cnt++;
				if (cnt >= argc)
                    throw "-safe-seq <number of sequences>";
				params.numseq_safe_scaling = convert_int(argv[cnt]);
                if (params.numseq_safe_scaling < 10)
                    throw "Too small -safe-seq";
				continue;
			}
            
            if (strcmp(argv[cnt], "--ignore-errors")==0) {
                params.ignore_any_errors = true;
                continue;
            }

            if (strcmp(argv[cnt], "--kernel-nonrev") == 0) {
                params.kernel_nonrev = true;
                continue;
            }

			if (strcmp(argv[cnt], "-f") == 0) {
				cnt++;
				if (cnt >= argc)
				        throw "Use -f <c | o | u | q | ry | ws | mk | <digits>>";
				if (strcmp(argv[cnt], "q") == 0 || strcmp(argv[cnt], "EQ") == 0)
					params.freq_type = FREQ_EQUAL;
				else if (strcmp(argv[cnt], "c") == 0
						|| strcmp(argv[cnt], "EM") == 0)
					params.freq_type = FREQ_EMPIRICAL;
				else if (strcmp(argv[cnt], "o") == 0
						|| strcmp(argv[cnt], "ES") == 0)
					params.freq_type = FREQ_ESTIMATE;
				else if (strcmp(argv[cnt], "u") == 0
						|| strcmp(argv[cnt], "UD") == 0)
					params.freq_type = FREQ_USER_DEFINED;
				else if (strcmp(argv[cnt], "ry") == 0
						|| strcmp(argv[cnt], "RY") == 0)
					params.freq_type = FREQ_DNA_RY;
				else if (strcmp(argv[cnt], "ws") == 0
						|| strcmp(argv[cnt], "WS") == 0)
					params.freq_type = FREQ_DNA_WS;
				else if (strcmp(argv[cnt], "mk") == 0
						|| strcmp(argv[cnt], "MK") == 0)
					params.freq_type = FREQ_DNA_MK;
				else
				        // throws error message if can't parse
				        params.freq_type = parseStateFreqDigits(argv[cnt]);
				continue;
			}

            if (strcmp(argv[cnt], "--keep-zero-freq") == 0) {
                params.keep_zero_freq = true;
                continue;
            }

            if (strcmp(argv[cnt], "--min-freq") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --min-freq NUM";
                params.min_state_freq = convert_double(argv[cnt]);
                if (params.min_state_freq <= 0)
                    throw "--min-freq must be positive";
                if (params.min_state_freq >= 1.0)
                    throw "--min-freq must be < 1.0";
                continue;
            }

            if (strcmp(argv[cnt], "--inc-zero-freq") == 0) {
                params.keep_zero_freq = false;
                continue;
            }

			if (strcmp(argv[cnt], "-fs") == 0 || strcmp(argv[cnt], "--site-freq") == 0) {
                if (params.tree_freq_file)
                    throw "Specifying both -fs and -ft not allowed";
				cnt++;
				if (cnt >= argc)
					throw "Use -fs <site_freq_file>";
				params.site_freq_file = argv[cnt];
//				params.SSE = LK_EIGEN;
				continue;
			}
			if (strcmp(argv[cnt], "-ft") == 0 || strcmp(argv[cnt], "--tree-freq") == 0) {
                if (params.site_freq_file)
                    throw "Specifying both -fs and -ft not allowed";
                cnt++;
				if (cnt >= argc)
					throw "Use -ft <treefile_to_infer_site_frequency_model>";
                params.tree_freq_file = argv[cnt];
                if (params.print_site_state_freq == WSF_NONE)
                    params.print_site_state_freq = WSF_POSTERIOR_MEAN;
                continue;
            }

			if (strcmp(argv[cnt], "-fconst") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -fconst <const_pattern_frequencies>";
				params.freq_const_patterns = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "--nrate") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -c <#rate_category>";
				params.num_rate_cats = convert_int(argv[cnt]);
				if (params.num_rate_cats < 1)
					throw "Wrong number of rate categories";
				continue;
			}
			if (strcmp(argv[cnt], "-cmin") == 0 || strcmp(argv[cnt], "--cmin") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -cmin <#min_rate_category>";
				params.min_rate_cats = convert_int(argv[cnt]);
				if (params.min_rate_cats < 2)
					throw "Wrong number of rate categories for -cmin";
				continue;
			}
			if (strcmp(argv[cnt], "-cmax") == 0 || strcmp(argv[cnt], "--cmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -cmax <#max_rate_category>";
				params.max_rate_cats = convert_int(argv[cnt]);
				if (params.max_rate_cats < 2)
					throw "Wrong number of rate categories for -cmax";
				continue;
			}
			if (strcmp(argv[cnt], "-a") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -a <gamma_shape>";
				params.gamma_shape = convert_double(argv[cnt]);
				if (params.gamma_shape <= 0)
					throw "Wrong gamma shape parameter (alpha)";
				continue;
			}

			if (strcmp(argv[cnt], "-amin") == 0 || strcmp(argv[cnt], "--alpha-min") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -amin <min_gamma_shape>";
				params.min_gamma_shape = convert_double(argv[cnt]);
				if (params.min_gamma_shape <= 0)
					throw "Wrong minimum gamma shape parameter (alpha)";
				continue;
			}

			if (strcmp(argv[cnt], "-gmean") == 0 || strcmp(argv[cnt], "--gamma-mean") == 0) {
				params.gamma_median = false;
				continue;
			}
			if (strcmp(argv[cnt], "-gmedian") == 0 || strcmp(argv[cnt], "--gamma-median") == 0) {
				params.gamma_median = true;
				continue;
			}
			if (strcmp(argv[cnt], "-i") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -i <p_invar_sites>";
				params.p_invar_sites = convert_double(argv[cnt]);
				if (params.p_invar_sites < 0)
					throw "Wrong number of proportion of invariable sites";
				continue;
			}
			if (strcmp(argv[cnt], "-optfromgiven") == 0) {
				params.optimize_from_given_params = true;
				continue;
			}
			if (strcmp(argv[cnt], "-brent") == 0) {
				params.optimize_by_newton = false;
				continue;
			}
			if (strcmp(argv[cnt], "-jointopt") == 0) {
				params.optimize_model_rate_joint = true;
				continue;
			}
			if (strcmp(argv[cnt], "-brent_ginvar") == 0) {
				params.optimize_model_rate_joint = false;
				continue;
			}
			if (strcmp(argv[cnt], "-fixbr") == 0 || strcmp(argv[cnt], "-blfix") == 0) {
				params.fixed_branch_length = BRLEN_FIX;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
				continue;
			}
			if (strcmp(argv[cnt], "-blscale") == 0) {
				params.fixed_branch_length = BRLEN_SCALE;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
				continue;
			}
			if (strcmp(argv[cnt], "-blmin") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -blmin <min_branch_length>";
				params.min_branch_length = convert_double(argv[cnt]);
				if (params.min_branch_length < 0.0)
					throw("Negative -blmin not allowed!");
				if (params.min_branch_length == 0.0)
					throw("Zero -blmin is not allowed due to numerical problems");
				if (params.min_branch_length > 0.1)
					throw("-blmin must be < 0.1");

				continue;
			}
			if (strcmp(argv[cnt], "-blmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -blmax <max_branch_length>";
				params.max_branch_length = convert_double(argv[cnt]);
				if (params.max_branch_length < 0.5)
					throw("-blmax smaller than 0.5 is not allowed");
				continue;
			}
            if (strcmp(argv[cnt], "--show-lh") == 0) {
                params.ignore_identical_seqs = false;
                params.fixed_branch_length = BRLEN_FIX;
                params.optimize_alg_gammai = "Brent";
                params.opt_gammai = false;
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
                verbose_mode = VB_DEBUG;
                params.ignore_checkpoint = true;
                continue;
            }
			if (strcmp(argv[cnt], "-sr") == 0) {
				params.stop_condition = SC_WEIBULL;
				cnt++;
				if (cnt >= argc)
					throw "Use -sr <#max_iteration>";
				params.max_iterations = convert_int(argv[cnt]);
				if (params.max_iterations <= params.min_iterations)
					throw "Specified max iteration must be greater than min iteration";
				continue;
			}
			if (strcmp(argv[cnt], "-nm") == 0 || strcmp(argv[cnt], "--nmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nm <#max_iteration>";
				params.max_iterations = convert_int(argv[cnt]);
				if (params.max_iterations <= params.min_iterations)
					throw "Specified max iteration must be greater than min iteration";
				continue;
			}
			if (strcmp(argv[cnt], "-sc") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sc <stop_confidence_value>";
				params.stop_confidence = convert_double(argv[cnt]);
				if (params.stop_confidence <= 0.5
						|| params.stop_confidence >= 1)
					throw "Stop confidence value must be in range (0.5,1)";
				continue;
			}
            if (strcmp(argv[cnt], "--runs") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --runs <number_of_runs>";
                params.num_runs = convert_int(argv[cnt]);
                if (params.num_runs < 1)
                    throw "Positive --runs please";
                continue;
            }
			if (strcmp(argv[cnt], "-gurobi") == 0) {
				params.gurobi_format = true;
				continue;
			}
			if (strcmp(argv[cnt], "-gthreads") == 0) {
				params.gurobi_format = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -gthreads <gurobi_threads>";
				params.gurobi_threads = convert_int(argv[cnt]);
				if (params.gurobi_threads < 1)
					throw "Wrong number of threads";
				continue;
			}
			if (strcmp(argv[cnt], "-b") == 0 || strcmp(argv[cnt], "--boot") == 0 ||
                strcmp(argv[cnt], "-j") == 0 || strcmp(argv[cnt], "--jack") == 0 ||
                strcmp(argv[cnt], "-bo") == 0 || strcmp(argv[cnt], "--bonly") == 0) {
				params.multi_tree = true;
				if (strcmp(argv[cnt], "-bo") == 0 || strcmp(argv[cnt], "--bonly") == 0)
					params.compute_ml_tree = false;
				else
					params.consensus_type = CT_CONSENSUS_TREE;
                if ((strcmp(argv[cnt], "-j") == 0 || strcmp(argv[cnt], "--jack") == 0) && params.jackknife_prop == 0.0)
                    params.jackknife_prop = 0.5;
				cnt++;
				if (cnt >= argc)
					throw "Use -b <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1)
					throw "Wrong number of bootstrap samples";
				if (params.num_bootstrap_samples == 1)
					params.compute_ml_tree = false;
				if (params.num_bootstrap_samples == 1)
					params.consensus_type = CT_NONE;
				continue;
			}
			if (strcmp(argv[cnt], "--bsam") == 0 || strcmp(argv[cnt], "-bsam") == 0 || strcmp(argv[cnt], "--sampling") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bsam <bootstrap_specification>";
				params.bootstrap_spec = argv[cnt];
                params.remove_empty_seq = false;
				continue;
			}
            
            if (strcmp(argv[cnt], "--subsample") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --subsample NUM";
                params.subsampling = convert_int(argv[cnt]);
                continue;
            }
            
            if (strcmp(argv[cnt], "--subsample-seed") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --subsample-seed <random_seed>";
                params.subsampling_seed = convert_int(argv[cnt]);
                continue;
            }

#ifdef USE_BOOSTER
            if (strcmp(argv[cnt], "--tbe") == 0) {
                params.transfer_bootstrap = 1;
                continue;
            }

            if (strcmp(argv[cnt], "--tbe-raw") == 0) {
                params.transfer_bootstrap = 2;
                continue;
            }
#endif

            if (strcmp(argv[cnt], "-bc") == 0 || strcmp(argv[cnt], "--bcon") == 0) {
				params.multi_tree = true;
				params.compute_ml_tree = false;
				cnt++;
				if (cnt >= argc)
					throw "Use -bc <num_bootstrap_samples>";
				params.num_bootstrap_samples = convert_int(argv[cnt]);
				if (params.num_bootstrap_samples < 1)
					throw "Wrong number of bootstrap samples";
				if (params.num_bootstrap_samples > 1)
					params.consensus_type = CT_CONSENSUS_TREE;
				continue;
			}
			if (strcmp(argv[cnt], "-iqppars") == 0) {
				params.iqp_assess_quartet = IQP_PARSIMONY;
				continue;
			}
			if (strcmp(argv[cnt], "-iqp") == 0) {
				params.iqp = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wct") == 0) {
				params.write_candidate_trees = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wt") == 0 || strcmp(argv[cnt], "--treels") == 0) {
				params.write_intermediate_trees = 1;
				continue;
			}

            if (strcmp(argv[cnt], "-wdt") == 0) {
                params.writeDistImdTrees = true;
                continue;
            }

            if (strcmp(argv[cnt], "-wtc") == 0) {
                params.write_intermediate_trees = 1;
                params.print_tree_lh = true;
                continue;
            }

			if (strcmp(argv[cnt], "-wt2") == 0) {
				params.write_intermediate_trees = 2;
//				params.avoid_duplicated_trees = true;
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wt3") == 0) {
				params.write_intermediate_trees = 3;
//				params.avoid_duplicated_trees = true;
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wbl") == 0) {
				params.print_branch_lengths = true;
				continue;
			}
            if (strcmp(argv[cnt], "-wit") == 0) {
                params.write_init_tree = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--write-branches") == 0) {
                params.write_branches = true;
                continue;
            }
            
//			if (strcmp(argv[cnt], "-nodup") == 0) {
//				params.avoid_duplicated_trees = true;
//				continue;
//			}
			if (strcmp(argv[cnt], "-rf_all") == 0 || strcmp(argv[cnt], "--tree-dist-all") == 0) {
				params.rf_dist_mode = RF_ALL_PAIR;
				continue;
			}
			if (strcmp(argv[cnt], "-rf_adj") == 0) {
				params.rf_dist_mode = RF_ADJACENT_PAIR;
				continue;
			}
			if (strcmp(argv[cnt], "-rf") == 0 || strcmp(argv[cnt], "--tree-dist") == 0) {
				params.rf_dist_mode = RF_TWO_TREE_SETS;
				cnt++;
				if (cnt >= argc)
					throw "Use -rf <second_tree>";
				params.second_tree = argv[cnt];
				continue;
			}
            if (strcmp(argv[cnt], "-rf1") == 0 || strcmp(argv[cnt], "--tree-dist1") == 0) {
                params.rf_dist_mode = RF_TWO_TREE_SETS;
                params.rf_same_pair = true;
                cnt++;
                if (cnt >= argc)
                    throw "Use --tree-dist1 <second_tree>";
                params.second_tree = argv[cnt];
                continue;
            }
			if (strcmp(argv[cnt], "-rf2") == 0 || strcmp(argv[cnt], "--tree-dist2") == 0) {
				params.rf_dist_mode = RF_TWO_TREE_SETS_EXTENDED;
				cnt++;
				if (cnt >= argc)
					throw "Use -rf2 <second_tree>";
				params.second_tree = argv[cnt];
				continue;
			}
            
            if (strcmp(argv[cnt], "--normalize-dist") == 0) {
                params.normalize_tree_dist = true;
                continue;
            }
            
			if (strcmp(argv[cnt], "-aLRT") == 0) {
				cnt++;
				if (cnt + 1 >= argc)
					throw "Use -aLRT <threshold%> <#replicates>";
				params.aLRT_threshold = convert_int(argv[cnt]);
				if (params.aLRT_threshold < 85 || params.aLRT_threshold > 101)
					throw "aLRT threshold must be between 85 and 100";
				cnt++;
				params.aLRT_replicates = convert_int(argv[cnt]);
				if (params.aLRT_replicates < 1000
						&& params.aLRT_replicates != 0)
					throw "aLRT replicates must be at least 1000";
				continue;
			}
			if (strcmp(argv[cnt], "-alrt") == 0 || strcmp(argv[cnt], "--alrt") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -alrt <#replicates | 0>";
                int reps = convert_int(argv[cnt]);
                if (reps == 0)
                    params.aLRT_test = true;
                else {
                    params.aLRT_replicates = reps;
                    if (params.aLRT_replicates < 1000)
                        throw "aLRT replicates must be at least 1000";
                }
				continue;
			}
			if (strcmp(argv[cnt], "-abayes") == 0 || strcmp(argv[cnt], "--abayes") == 0) {
				params.aBayes_test = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lbp") == 0 || strcmp(argv[cnt], "--lbp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lbp <#replicates>";
				params.localbp_replicates = convert_int(argv[cnt]);
				if (params.localbp_replicates < 1000
						&& params.localbp_replicates != 0)
					throw "Local bootstrap (LBP) replicates must be at least 1000";
				continue;
			}
			if (strcmp(argv[cnt], "-wsl") == 0 || strcmp(argv[cnt], "--sitelh") == 0) {
				params.print_site_lh = WSL_SITE;
				continue;
			}

			if (strcmp(argv[cnt], "-wpl") == 0 || strcmp(argv[cnt], "--partlh") == 0) {
				params.print_partition_lh = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wslg") == 0 || strcmp(argv[cnt], "-wslr") == 0) {
				params.print_site_lh = WSL_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-wslm") == 0) {
				params.print_site_lh = WSL_MIXTURE;
				continue;
			}
			if (strcmp(argv[cnt], "-wslmr") == 0 || strcmp(argv[cnt], "-wslrm") == 0) {
				params.print_site_lh = WSL_MIXTURE_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-wspr") == 0) {
				params.print_site_prob = WSL_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-wspm") == 0) {
				params.print_site_prob = WSL_MIXTURE;
				continue;
			}
			if (strcmp(argv[cnt], "-wspmr") == 0 || strcmp(argv[cnt], "-wsprm") == 0) {
				params.print_site_prob = WSL_MIXTURE_RATECAT;
				continue;
			}

			if (strcmp(argv[cnt], "-asr") == 0 || strcmp(argv[cnt], "--ancestral") == 0) {
				params.print_ancestral_sequence = AST_MARGINAL;
                params.ignore_identical_seqs = false;
				continue;
			}

			if (strcmp(argv[cnt], "-asr-min") == 0 || strcmp(argv[cnt], "--asr-min") == 0) {
                cnt++;
				if (cnt >= argc)
					throw "Use -asr-min <probability>";
                
                params.min_ancestral_prob = convert_double(argv[cnt]);
                if (params.min_ancestral_prob < 0 || params.min_ancestral_prob > 1)
                    throw "Minimum ancestral probability [-asr-min] must be between 0 and 1.0";
                continue;
            }

			if (strcmp(argv[cnt], "-asr-joint") == 0) {
				params.print_ancestral_sequence = AST_JOINT;
                params.ignore_identical_seqs = false;
				continue;
			}

			if (strcmp(argv[cnt], "-wsr") == 0 || strcmp(argv[cnt], "--rate") == 0) {
				params.print_site_rate |= 1;
				continue;
			}

            if (strcmp(argv[cnt], "--mlrate") == 0) {
                params.print_site_rate |= 2;
                continue;
            }

            if (strcmp(argv[cnt], "-wsptrees") == 0) {
				params.print_trees_site_posterior = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-wsf") == 0) {
				params.print_site_state_freq = WSF_POSTERIOR_MEAN;
				continue;
			}
			if (strcmp(argv[cnt], "--freq-max") == 0 || strcmp(argv[cnt], "-fmax") == 0) {
				params.print_site_state_freq = WSF_POSTERIOR_MAX;
				continue;
			}
			if (strcmp(argv[cnt], "-wba") == 0) {
				params.print_bootaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wbsf") == 0) {
				params.print_boot_site_freq = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wsa") == 0) {
				params.print_subaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wtl") == 0) {
				params.print_tree_lh = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wpi") == 0) {
				params.print_partition_info = true;
				params.print_conaln = true;
				continue;
			}
			if (strcmp(argv[cnt], "-wca") == 0) {
				params.print_conaln = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wsplits") == 0) {
				params.print_splits_file = true;
				continue;
			}
            if (strcmp(argv[cnt], "--no-splits.nex") == 0) {
                params.print_splits_nex_file = false;
                continue;
            }
			if (strcmp(argv[cnt], "-ns") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ns <num_simulations>";
				params.whtest_simulations = convert_int(argv[cnt]);
				if (params.whtest_simulations < 1)
					throw "Wrong number of simulations for WH-test";
				continue;
			}
			if (strcmp(argv[cnt], "-mr") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mr <rate_file>";
				params.rate_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-cat_mean") == 0) {
				params.mcat_type |= MCAT_MEAN;
				continue;
			}
			if (strcmp(argv[cnt], "-cat_nolog") == 0) {
				params.mcat_type &= (127 - MCAT_LOG);
				continue;
			}
			if (strcmp(argv[cnt], "-cat_site") == 0) {
				params.mcat_type &= (127 - MCAT_PATTERN);
				continue;
			}
			if (strcmp(argv[cnt], "-tina") == 0) {
				params.do_pars_multistate = true;
                params.ignore_checkpoint = true;
				continue;
			}
			if (strcmp(argv[cnt], "-pval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -pval <gene_pvalue_file>";
				params.gene_pvalue_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-nnitest") == 0) {
				params.testNNI = true;
				continue;
			}
			if (strcmp(argv[cnt], "-anni") == 0) {
				params.approximate_nni = true;
				continue;
			}
			if (strcmp(argv[cnt], "-nnicut") == 0) {
				params.estimate_nni_cutoff = true;
				//nni_cutoff = -5.41/2;
				continue;
			}
			if (strcmp(argv[cnt], "-nnichi2") == 0) {
				params.nni_cutoff = -5.41 / 2;
				continue;
			}
			if (strcmp(argv[cnt], "-nnicutval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nnicutval <log_diff_value>";
				params.nni_cutoff = convert_double(argv[cnt]);
				if (params.nni_cutoff >= 0)
					throw "cutoff value for -nnicutval must be negative";
				continue;
			}
			if (strcmp(argv[cnt], "-nnisort") == 0) {
				params.nni_sort = true;
				continue;
			}
			if (strcmp(argv[cnt], "-plog") == 0) {
				params.gene_pvalue_loga = true;
				continue;
			}
			if (strcmp(argv[cnt], "-dmp") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmp <ncbi_taxid>";
				params.ncbi_taxid = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-dmplevel") == 0
					|| strcmp(argv[cnt], "-dmprank") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmprank <ncbi_taxon_rank>";
				params.ncbi_taxon_level = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dmpignore") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmpignore <ncbi_ignore_level>";
				params.ncbi_ignore_level = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-dmpname") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -dmpname <ncbi_names_file>";
				params.ncbi_names_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-eco") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -eco <eco_dag_file>";
				params.eco_dag_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-k%") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -k% <k in %>";
				//convert_range(argv[cnt], params.k_percent, params.sub_size, params.step_size);
				params.k_percent = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-diet") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -diet <d in %>";
				convert_range(argv[cnt], params.diet_min, params.diet_max,
						params.diet_step);
				//params.diet = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-up") == 0) {
				params.upper_bound = true;
				continue;
			}
			if (strcmp(argv[cnt], "-upNNI") == 0) {
 				params.upper_bound_NNI = true;
			}
			if (strcmp(argv[cnt], "-upFrac") == 0) {
				cnt++;
				if (cnt >= argc)
				  throw "Use -upFrac <fraction>";
				params.upper_bound_frac = convert_double(argv[cnt]);
			}
			if (strcmp(argv[cnt], "-ecoR") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ecoR <run number>";
				params.eco_run = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bb") == 0 || strcmp(argv[cnt], "-B") == 0 || strcmp(argv[cnt], "--ufboot") == 0 ||
                strcmp(argv[cnt], "-J") == 0 || strcmp(argv[cnt], "--ufjack") == 0) {
                if ((strcmp(argv[cnt], "-J") == 0 || strcmp(argv[cnt], "--ufjack") == 0) && params.jackknife_prop == 0.0)
                    params.jackknife_prop = 0.5;
				cnt++;
				if (cnt >= argc)
					throw "Use -B <#replicates>";
                if (params.stop_condition == SC_FIXED_ITERATION) {
                    throw("Ultrafast bootstrap does not work with -fast, -te or -n option");
                }
				params.gbo_replicates = convert_int(argv[cnt]);
//				params.avoid_duplicated_trees = true;
				if (params.gbo_replicates < 1000)
					throw "#replicates must be >= 1000";
				params.consensus_type = CT_CONSENSUS_TREE;
				params.stop_condition = SC_BOOTSTRAP_CORRELATION;
				//params.nni5Branches = true;
				continue;
			}
			if (strcmp(argv[cnt], "-beps") == 0 || strcmp(argv[cnt], "--beps") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -beps <epsilon>";
				params.ufboot_epsilon = convert_double(argv[cnt]);
				if (params.ufboot_epsilon <= 0.0)
					throw "Epsilon must be positive";
				continue;
			}
			if (strcmp(argv[cnt], "-wbt") == 0 || strcmp(argv[cnt], "--wbt") == 0 || strcmp(argv[cnt], "--boot-trees") == 0) {
				params.print_ufboot_trees = 1;
				continue;
			}
			if (strcmp(argv[cnt], "-wbtl") == 0 || strcmp(argv[cnt], "--wbtl") == 0) {
                // print ufboot trees with branch lengths
				params.print_ufboot_trees = 2;
				continue;
			}
			if (strcmp(argv[cnt], "-bs") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bs <begin_sampling_size>";
				params.check_gbo_sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bmax") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bmax <max_candidate_trees>";
				params.max_candidate_trees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bcor") == 0 | strcmp(argv[cnt], "--bcor") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bcor <min_correlation>";
				params.min_correlation = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "--bnni") == 0 || strcmp(argv[cnt], "-bnni") == 0) {
				params.ufboot2corr = true;
                // print ufboot trees with branch lengths
//				params.print_ufboot_trees = 2; // Diep: relocate to be below this for loop
				continue;
			}
			if (strcmp(argv[cnt], "-u2c_nni5") == 0) {
				params.u2c_nni5 = true;
				continue;
			}

			if (strcmp(argv[cnt], "-nstep") == 0 || strcmp(argv[cnt], "--nstep") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nstep <step_iterations>";
				params.step_iterations = convert_int(argv[cnt]);
				if (params.step_iterations < 10
						|| params.step_iterations % 2 == 1)
					throw "At least step size of 10 and even number please";
				params.min_iterations = params.step_iterations;
				continue;
			}
			if (strcmp(argv[cnt], "-boff") == 0) {
				params.online_bootstrap = false;
				continue;
			}
//			if (strcmp(argv[cnt], "-nostore") == 0
//					|| strcmp(argv[cnt], "-memsave") == 0) {
//				params.store_candidate_trees = false;
//				continue;
//			}
            
            if (strcmp(argv[cnt], "--jack-prop") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --jack-prop jackknife_proportion";
                params.jackknife_prop = convert_double(argv[cnt]);
                if (params.jackknife_prop <= 0.0 || params.jackknife_prop >= 1.0)
                    throw "Jackknife proportion must be between 0.0 and 1.0";
                continue;
            }
            
            if (strcmp(argv[cnt], "--robust-phy") == 0) {
                if (params.robust_median)
                    throw "Can't couple --robust-phy with --robust-median";
                cnt++;
                if (cnt >= argc)
                    throw "Use --robust-phy proportion_of_best_sites_to_keep";
                params.robust_phy_keep = convert_double(argv[cnt]);
                if (params.robust_phy_keep <= 0.0 || params.robust_phy_keep > 1.0)
                    throw "--robust-phy parameter must be between 0 and 1";
                // TODO: use Brent (instead of Newton) optimisation of branch lengths
                params.optimize_by_newton = false;
                params.optimize_alg_gammai = "Brent";
                params.optimize_alg_freerate = "2-BFGS";
                continue;
            }

            if (strcmp(argv[cnt], "--robust-median") == 0) {
                if (params.robust_phy_keep < 1.0)
                    throw "Can't couple --robust-phy with --robust-median";
                params.robust_median = true;
                // TODO: use Brent (instead of Newton) optimisation of branch lengths
                params.optimize_by_newton = false;
                params.optimize_alg_gammai = "Brent";
                params.optimize_alg_freerate = "2-BFGS";
                continue;
            }

			if (strcmp(argv[cnt], "-mem") == 0 || strcmp(argv[cnt], "--mem") == 0) {
				cnt++;
				if (cnt >= argc)
                    throw "Use -mem max_mem_size";
				params.lh_mem_save = LM_MEM_SAVE;
                int end_pos;
                double mem = convert_double(argv[cnt], end_pos);
                if (mem < 0)
                    throw "-mem must be non-negative";
                if (argv[cnt][end_pos] == 'G') {
                    params.max_mem_size = mem * 1073741824.0;
                } else if (argv[cnt][end_pos] == 'M') {
                    params.max_mem_size = mem * 1048576.0;
                } else if (argv[cnt][end_pos] == '%'){
                    params.max_mem_size = mem * 0.01;
                    if (params.max_mem_size > 1)
                        throw "-mem percentage must be between 0 and 100";
                } else {
                    if (mem > 1)
                        throw "Invalid -mem option. Example: -mem 200M, -mem 10G";
                    params.max_mem_size = mem;
                }
				continue;
			}
            if (strcmp(argv[cnt], "--save-mem-buffer") == 0) {
                params.buffer_mem_save = true;
                continue;
            }
            if (strcmp(argv[cnt], "--no-save-mem-buffer") == 0) {
                params.buffer_mem_save = false;
                continue;
            }
//			if (strcmp(argv[cnt], "-storetrees") == 0) {
//				params.store_candidate_trees = true;
//				continue;
//			}
			if (strcmp(argv[cnt], "-nodiff") == 0) {
				params.distinct_trees = false;
				continue;
			}
			if (strcmp(argv[cnt], "-treediff") == 0) {
				params.distinct_trees = true;
				continue;
			}
			if (strcmp(argv[cnt], "-norell") == 0) {
				params.use_rell_method = false;
				continue;
			}
			if (strcmp(argv[cnt], "-elw") == 0) {
				params.use_elw_method = true;
				continue;
			}
			if (strcmp(argv[cnt], "-noweight") == 0) {
				params.use_weighted_bootstrap = false;
				continue;
			}
			if (strcmp(argv[cnt], "-nomore") == 0) {
				params.use_max_tree_per_bootstrap = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bweight") == 0) {
				params.use_weighted_bootstrap = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bmore") == 0) {
				params.use_max_tree_per_bootstrap = false;
				continue;
			}
			if (strcmp(argv[cnt], "-gz") == 0) {
				params.do_compression = true;
				continue;
			}
			if (strcmp(argv[cnt], "-newheu") == 0) {
				params.new_heuristic = true;
				// Enable RAxML kernel
				continue;
			}
			if (strcmp(argv[cnt], "-maxtime") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -maxtime <time_in_minutes>";
				params.maxtime = convert_double(argv[cnt]);
				params.min_iterations = 1000000;
				params.stop_condition = SC_REAL_TIME;
				continue;
			}
			if (strcmp(argv[cnt], "--ninit") == 0 || strcmp(argv[cnt], "-ninit") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ninit <number_of_parsimony_trees>";
				params.numInitTrees = convert_int(argv[cnt]);
                if (params.numInitTrees < 0)
                    throw "-ninit must be non-negative";
				if (params.numInitTrees < params.numNNITrees)
					params.numNNITrees = params.numInitTrees;
				continue;
			}
			if (strcmp(argv[cnt], "-fast") == 0 || strcmp(argv[cnt], "--fast") == 0) {
                // fast search option to resemble FastTree
                if (params.gbo_replicates != 0) {
                    throw("Ultrafast bootstrap (-bb) does not work with -fast option");
                }
                params.numInitTrees = 2;
                if (params.min_iterations == -1)
                    params.min_iterations = 2;
				params.stop_condition = SC_FIXED_ITERATION;
                params.modelEps = 0.05;
                params.suppress_list_of_sequences = true;
                params.suppress_zero_distance_warnings = true;
                params.suppress_duplicate_sequence_warnings = true;
                params.optimize_alg_freerate = "1-BFGS";
                params.opt_gammai = false;
                continue;
            }
			if (strcmp(argv[cnt], "-fss") == 0) {
				params.fixStableSplits = true;
//				params.five_plus_five = true;
				continue;
			}
            if (strcmp(argv[cnt], "--stable-thres") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --stable-thres <support_value_threshold>";
                params.stableSplitThreshold = convert_double(argv[cnt]);
                continue;
            }
			if (strcmp(argv[cnt], "-ff") == 0) {
				params.five_plus_five = true;
				continue;
			}

			if (strcmp(argv[cnt], "-tabu") == 0) {
                params.fixStableSplits = true;
				params.tabu = true;
                params.maxCandidates = params.numSupportTrees;
				continue;
			}

            if (strcmp(argv[cnt], "--adt-pert") == 0) {
                if (params.tabu == true) {
                    throw("option -tabu and --adt-pert cannot be combined");
                }
                params.adaptPertubation = true;
                params.stableSplitThreshold = 1.0;
                continue;
            }

            if (strcmp(argv[cnt], "-memcheck") == 0) {
                params.memCheck = true;
                continue;
            }

			if (strcmp(argv[cnt], "--ntop") == 0 || strcmp(argv[cnt], "-ntop") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ntop <number_of_top_parsimony_trees>";
				params.numNNITrees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "--num-sup-trees") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use --num-sup-trees <number_of_support_trees>";
				params.numSupportTrees = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-fixai") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -fixai <alpha_invar_file>";
				params.alpha_invar_file = argv[cnt];
				continue;
			}

            if (strcmp(argv[cnt], "--opt-gamma-inv") == 0) {
                params.opt_gammai = true;
                continue;
            }

            if (strcmp(argv[cnt], "--no-opt-gamma-inv") == 0) {
                params.opt_gammai = false;
                continue;
            }

            if (strcmp(argv[cnt], "--opt-gammai-fast") == 0) {
                params.opt_gammai_fast = true;
                params.opt_gammai = true;
                continue;
            }

            if (strcmp(argv[cnt], "--opt-gammai-kb") == 0) {
                params.opt_gammai_keep_bran = true;
                params.opt_gammai = true;
                continue;
            }

            if (strcmp(argv[cnt], "--adaptive-eps") == 0) {
                params.testAlphaEpsAdaptive = true;
                continue;
            }
            if (strcmp(argv[cnt], "--rand-alpha") == 0) {
                params.randomAlpha = true;
                continue;
            }

            if (strcmp(argv[cnt], "-eai") == 0) {
                params.exh_ai = true;
                continue;
            }
			if (strcmp(argv[cnt], "-poplim") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -poplim <max_pop_size>";
				params.maxCandidates = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "--nbest") == 0 ||strcmp(argv[cnt], "-nbest") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nbest <number_of_candidate_trees>";
				params.popSize = convert_int(argv[cnt]);
				ASSERT(params.popSize < params.numInitTrees);
				continue;
			}
			if (strcmp(argv[cnt], "-beststart") == 0) {
				params.bestStart = true;
				cnt++;
				if (cnt >= argc)
					throw "Use -best_start <binary_alignment_file>";
				params.binary_aln_file = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-pll") == 0) {
                throw("-pll option is discontinued.");
				params.pll = true;
				continue;
			}
			if (strcmp(argv[cnt], "-me") == 0 || strcmp(argv[cnt], "--epsilon") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -me <model_epsilon>";
				params.modelEps = convert_double(argv[cnt]);
				if (params.modelEps <= 0.0)
					throw "Model epsilon must be positive";
				if (params.modelEps > 1.0)
					throw "Model epsilon must not be larger than 1.0";
				continue;
			}

            if (strcmp(argv[cnt], "--mf-epsilon") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --mf-epsilon <modelfinder_epsilon>";
                params.modelfinder_eps = convert_double(argv[cnt]);
                if (params.modelfinder_eps <= 0.0)
                    throw "ModelFinder epsilon must be positive";
                if (params.modelEps > 1.0)
                    throw "ModelFinder epsilon must not be larger than 1.0";
                continue;
            }

            if (strcmp(argv[cnt], "-pars_ins") == 0) {
				params.reinsert_par = true;
				continue;
			}
			if (strcmp(argv[cnt], "-allnni") == 0 || strcmp(argv[cnt], "--allnni") == 0) {
				params.speednni = false;
				continue;
			}
            
			if (strcmp(argv[cnt], "-snni") == 0) {
				params.snni = true;
				// dont need to turn this on here
				//params.autostop = true;
				//params.speednni = true;
				// Minh: why do you turn this on? it doubles curPerStrength at some point
				//params.adaptPert = true;
				continue;
			}
			if (strcmp(argv[cnt], "-iqpnni") == 0) {
				params.snni = false;
				params.start_tree = STT_BIONJ;
				params.numNNITrees = 1;
//            continue; } if (strcmp(argv[cnt], "-auto") == 0) {
//            	params.autostop = true;
				continue;
			}
			if (strcmp(argv[cnt], "--nstop") == 0 || strcmp(argv[cnt], "-nstop") == 0) {
				if (params.stop_condition != SC_BOOTSTRAP_CORRELATION)
					params.stop_condition = SC_UNSUCCESS_ITERATION;
				cnt++;
				if (cnt >= argc)
					throw "Use -nstop <#iterations>";
				params.unsuccess_iteration = convert_int(argv[cnt]);
                if (params.unsuccess_iteration <= 0)
                    throw "-nstop iterations must be positive";
                params.max_iterations = max(params.max_iterations, params.unsuccess_iteration*10);
				continue;
			}
			if (strcmp(argv[cnt], "-lsbran") == 0) {
				params.leastSquareBranch = true;
				continue;
			}
			if (strcmp(argv[cnt], "-manuel") == 0) {
				params.manuel_analytic_approx = true;
				continue;
			}
			if (strcmp(argv[cnt], "-parsbran") == 0) {
				params.pars_branch_length = true;
				continue;
			}
			if (strcmp(argv[cnt], "-bayesbran") == 0) {
				params.bayes_branch_length = true;
				continue;
			}
			if (strcmp(argv[cnt], "-fivebran") == 0
					|| strcmp(argv[cnt], "-nni5") == 0) {
				params.nni5 = true;
				params.nni_type = NNI5;
				continue;
			}
			if (strcmp(argv[cnt], "-onebran") == 0
					|| strcmp(argv[cnt], "-nni1") == 0) {
				params.nni_type = NNI1;
				params.nni5 = false;
				continue;
			}
            
            if (strcmp(argv[cnt], "-nni-eval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -nni-eval <num_evaluation>";
                params.nni5_num_eval = convert_int(argv[cnt]);
                if (params.nni5_num_eval < 1)
                    throw("Positive -nni-eval expected");
                continue;
            }

            if (strcmp(argv[cnt], "-bl-eval") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bl-eval <num_evaluation>";
                params.brlen_num_traversal = convert_int(argv[cnt]);
                if (params.brlen_num_traversal < 1)
                    throw("Positive -bl-eval expected");
                continue;
            }
            
			if (strcmp(argv[cnt], "-smooth") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -smooth <num_iterations>";
				params.numSmoothTree = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-lsnni") == 0) {
				params.leastSquareNNI = true;
				continue;
			}
			if (strcmp(argv[cnt], "-lsvar") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lsvar <o|ft|fm|st|p>";
				if (strcmp(argv[cnt], "o") == 0
						|| strcmp(argv[cnt], "ols") == 0) {
					params.ls_var_type = OLS;
					continue;
				}
				if (strcmp(argv[cnt], "ft") == 0
						|| strcmp(argv[cnt], "first_taylor") == 0) {
					params.ls_var_type = WLS_FIRST_TAYLOR;
					continue;
				}
				if (strcmp(argv[cnt], "fm") == 0
						|| strcmp(argv[cnt], "fitch_margoliash") == 0) {
					params.ls_var_type = WLS_FITCH_MARGOLIASH;
					continue;
				}
				if (strcmp(argv[cnt], "st") == 0
						|| strcmp(argv[cnt], "second_taylor") == 0) {
					params.ls_var_type = WLS_SECOND_TAYLOR;
					continue;
				}
				if (strcmp(argv[cnt], "p") == 0
						|| strcmp(argv[cnt], "pauplin") == 0) {
					params.ls_var_type = WLS_PAUPLIN;
				} else {
					throw "Use -lsvar <o|ft|fm|st|p>";
				}
				continue;
			}
			if (strcmp(argv[cnt], "-eps") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -eps <log-likelihood epsilon>";
				params.loglh_epsilon = convert_double(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-pb") == 0) { // Enable parsimony branch length estimation
				params.parbran = true;
				continue;
			}
			if (strcmp(argv[cnt], "-x") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -x <iteration_multiple>";
				params.iteration_multiple = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-sp_iter") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sp_iter <number_iteration>";
				params.speedup_iter = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-avh") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -avh <arndt_#bootstrap>";
				params.avh_test = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bootlh") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bootlh <#replicates>";
				params.bootlh_test = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-bootpart") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -bootpart <part1_length,part2_length,...>";
				params.bootlh_partitions = argv[cnt];
				continue;
			}
			if (strcmp(argv[cnt], "-AIC") == 0) {
				params.model_test_criterion = MTC_AIC;
				continue;
			}
			if (strcmp(argv[cnt], "-AICc") == 0 || strcmp(argv[cnt], "-AICC") == 0) {
				params.model_test_criterion = MTC_AICC;
				continue;
			}
			if (strcmp(argv[cnt], "-merit") == 0 || strcmp(argv[cnt], "--merit") == 0) {
                cnt++;
				if (cnt >= argc)
					throw "Use -merit AIC|AICC|BIC";
                if (strcmp(argv[cnt], "AIC") == 0)
                    params.model_test_criterion = MTC_AIC;
                else if (strcmp(argv[cnt], "AICc") == 0 || strcmp(argv[cnt], "AICC") == 0)
                    params.model_test_criterion = MTC_AICC;
                else if (strcmp(argv[cnt], "BIC") == 0)
                    params.model_test_criterion = MTC_BIC;
                else throw "Use -merit AIC|AICC|BIC";
				continue;
			}
			if (strcmp(argv[cnt], "-ms") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -ms <model_test_sample_size>";
				params.model_test_sample_size = convert_int(argv[cnt]);
				continue;
			}
			if (strcmp(argv[cnt], "-nt") == 0 || strcmp(argv[cnt], "-c") == 0 ||
                strcmp(argv[cnt], "-T") == 0  || strcmp(argv[cnt], "--threads") == 0) {
				cnt++;
				if (cnt >= argc)
				throw "Use -nt <num_threads|AUTO>";
                if (iEquals(argv[cnt], "AUTO"))
                    params.num_threads = 0;
                else {
                    params.num_threads = convert_int(argv[cnt]);
                    if (params.num_threads < 1)
                        throw "At least 1 thread please";
                }
				continue;
			}
            
            if (strcmp(argv[cnt], "-ntmax") == 0 || strcmp(argv[cnt], "--threads-max") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -ntmax <num_threads_max>";
                params.num_threads_max = convert_int(argv[cnt]);
                if (params.num_threads_max < 1)
                    throw "At least 1 thread please";
                continue;
            }
            
            if (strcmp(argv[cnt], "--thread-model") == 0) {
                params.openmp_by_model = true;
                continue;
            }

            if (strcmp(argv[cnt], "--thread-site") == 0) {
                params.openmp_by_model = false;
                continue;
            }

//			if (strcmp(argv[cnt], "-rootstate") == 0) {
//                cnt++;
//                if (cnt >= argc)
//                    throw "Use -rootstate <rootstate>";
//                params.root_state = argv[cnt];
//                params.SSE = LK_NORMAL;
//                continue;
//			}
			if (strcmp(argv[cnt], "-ct") == 0) {
            	params.count_trees = true;
            	continue;
			}
			if (strcmp(argv[cnt], "--sprrad") == 0 || strcmp(argv[cnt], "--radius") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -sprrad <SPR radius used in parsimony search>";
				params.sprDist = convert_int(argv[cnt]);
				continue;
			}
            
            if (strcmp(argv[cnt], "--mpcost") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --mpcost <parsimony_cost_file>";
                params.sankoff_cost_file = argv[cnt];
                continue;
            }
            
			if (strcmp(argv[cnt], "-no_rescale_gamma_invar") == 0) {
				params.no_rescale_gamma_invar = true;
				continue;
			}

			if (strcmp(argv[cnt], "-wsi") == 0) {
				params.compute_seq_identity_along_tree = true;
				continue;
			}
            
            if (strcmp(argv[cnt], "--no-seq-comp") == 0) {
                params.compute_seq_composition = false;
                continue;
            }
            
			if (strcmp(argv[cnt], "-t") == 0 || strcmp(argv[cnt], "-te") == 0 || strcmp(argv[cnt], "--tree") == 0) {
                if (strcmp(argv[cnt], "-te") == 0) {
                    if (params.gbo_replicates != 0) {
                        throw("Ultrafast bootstrap does not work with -te option");
                    }
                    params.min_iterations = 0;
                    params.stop_condition = SC_FIXED_ITERATION;
                }
				cnt++;
				if (cnt >= argc)
					throw "Use -t,-te <start_tree | BIONJ | PARS | PLLPARS | PJ | RANDOM>";
				else if (strcmp(argv[cnt], "PARS") == 0)
					params.start_tree = STT_PARSIMONY;
                else if (strcmp(argv[cnt], "PJ") == 0)
                    params.start_tree = STT_PARSIMONY_JOINING;
				else if (strcmp(argv[cnt], "PLLPARS") == 0)
					params.start_tree = STT_PLL_PARSIMONY;
                else if (strcmp(argv[cnt], "RANDOM") == 0 || strcmp(argv[cnt], "RAND") == 0)
					params.start_tree = STT_RANDOM_TREE;
                else if (START_TREE_RECOGNIZED(argv[cnt])) {
                    params.start_tree_subtype_name = argv[cnt];
                    params.start_tree = STT_BIONJ;
                } else {
                    params.user_file = argv[cnt];
                    if (params.min_iterations == 0)
                        params.start_tree = STT_USER_TREE;
                }
				continue;
			}
            
            if (strcmp(argv[cnt], "--no-ml-tree") == 0) {
                params.modelfinder_ml_tree = false;
                continue;
            }
            
            if (strcmp(argv[cnt], "--tree-fix") == 0) {
                if (params.gbo_replicates != 0) {
                    outError("Ultrafast bootstrap does not work with -te option");
                }
                params.min_iterations = 0;
                params.stop_condition = SC_FIXED_ITERATION;
            }
                
            if (strcmp(argv[cnt], "-g") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use -g <constraint_tree>";
                params.constraint_tree_file = argv[cnt];
                continue;
            }
            
			if (strcmp(argv[cnt], "-lmap") == 0 || strcmp(argv[cnt], "--lmap") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lmap <likelihood_mapping_num_quartets>";
                if (iEquals(argv[cnt], "ALL")) {
                    params.lmap_num_quartets = 0;
                } else {
                    params.lmap_num_quartets = convert_int64(argv[cnt]);
                    if (params.lmap_num_quartets < 0)
                        throw "Number of quartets must be >= 1";
                }
				continue;
			}

			if (strcmp(argv[cnt], "-lmclust") == 0 || strcmp(argv[cnt], "--lmclust") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -lmclust <likelihood_mapping_cluster_file>";
				params.lmap_cluster_file = argv[cnt];
				// '-keep_ident' is currently required to allow a 1-to-1 mapping of the 
				// user-given groups (HAS) - possibly obsolete in the future versions
				params.ignore_identical_seqs = false;
                if (params.lmap_num_quartets < 0)
                    params.lmap_num_quartets = 0;
				continue;
			}

			if (strcmp(argv[cnt], "-wql") == 0 || strcmp(argv[cnt], "--quartetlh") == 0) {
				params.print_lmap_quartet_lh = true;
				continue;
			}

			if (strcmp(argv[cnt], "-mixlen") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -mixlen <number of mixture branch lengths for heterotachy model>";
				params.num_mixlen = convert_int(argv[cnt]);
				if (params.num_mixlen < 1)
					throw("-mixlen must be >= 1");
				continue;
			}
            
			if (strcmp(argv[cnt], "--link-alpha") == 0) {
				params.link_alpha = true;
				continue;
			}

            if (strcmp(argv[cnt], "--link-model") == 0) {
                params.link_model = true;
                continue;
            }

            if (strcmp(argv[cnt], "--model-joint") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --model-joint MODEL_NAME";
                params.model_joint = argv[cnt];
                params.link_model = true;
                continue;
            }

            if (strcmp(argv[cnt], "--unlink-tree") == 0) {
                params.partition_type = TOPO_UNLINKED;
                params.ignore_identical_seqs = false;
                continue;
            }
            
			if (strcmp(argv[cnt], "-redo") == 0 || strcmp(argv[cnt], "--redo") == 0) {
				params.ignore_checkpoint = true;
                // 2020-04-27: SEMANTIC CHANGE: also redo ModelFinder
                params.model_test_again = true;
				continue;
			}

            if (strcmp(argv[cnt], "-tredo") == 0 || strcmp(argv[cnt], "--tredo") == 0 || strcmp(argv[cnt], "--redo-tree") == 0) {
                params.ignore_checkpoint = true;
                continue;
            }

			if (strcmp(argv[cnt], "-undo") == 0 || strcmp(argv[cnt], "--undo") == 0) {
				params.force_unfinished = true;
				continue;
			}

			if (strcmp(argv[cnt], "-cptime") == 0 || strcmp(argv[cnt], "--cptime") == 0) {
				cnt++;
				if (cnt >= argc)
					throw "Use -cptime <checkpoint_time_interval>";
				params.checkpoint_dump_interval = convert_int(argv[cnt]);
				continue;
			}
            
			if (strcmp(argv[cnt], "--no-log") == 0) {
				params.suppress_output_flags |= OUT_LOG;
				continue;
			}

			if (strcmp(argv[cnt], "--no-treefile") == 0) {
				params.suppress_output_flags |= OUT_TREEFILE;
				continue;
			}
			if (strcmp(argv[cnt], "--no-iqtree") == 0) {
				params.suppress_output_flags |= OUT_IQTREE;
				continue;
			}
			if (strcmp(argv[cnt], "--no-outfiles") == 0) {
				params.suppress_output_flags |= OUT_LOG + OUT_TREEFILE + OUT_IQTREE;
				continue;
			}

            // -- Mon Apr 17 21:18:23 BST 2017
            // DONE Minh: merged correctly.
            if (strcmp(argv[cnt], "--scaling-squaring") == 0) {
                params.matrix_exp_technique = MET_SCALING_SQUARING;
                continue;
            }
            if (strcmp(argv[cnt], "--eigenlib") == 0) {
                params.matrix_exp_technique = MET_EIGEN3LIB_DECOMPOSITION;
                continue;
            }
            if (strcmp(argv[cnt], "--eigen") == 0) {
                params.matrix_exp_technique = MET_EIGEN_DECOMPOSITION;
                continue;
            }
            if (strcmp(argv[cnt], "--lie-markov") == 0) {
                params.matrix_exp_technique = MET_LIE_MARKOV_DECOMPOSITION;
                continue;
            }            
			if (strcmp(argv[cnt], "--no-uniqueseq") == 0) {
				params.suppress_output_flags |= OUT_UNIQUESEQ;
				continue;
			}
            // --

            if (strcmp(argv[cnt], "--dating") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --dating LSD";
                params.dating_method = argv[cnt];
                if (params.dating_method != "LSD")
                    throw "Currently only LSD (least-square dating) method is supported";
                continue;
            }

            if (strcmp(argv[cnt], "--date") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date <date_file>|TAXNAME";
                if (params.dating_method == "")
                    params.dating_method = "LSD";
                params.date_file = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "--date-tip") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-tip <YYYY[-MM-DD]>";
                if (params.dating_method == "")
                    params.dating_method = "LSD";
                params.date_tip = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "--date-root") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-root <YYYY[-MM-DD]>";
                if (params.dating_method == "")
                    params.dating_method = "LSD";
                params.date_root = argv[cnt];
                continue;
            }

            if (strcmp(argv[cnt], "--date-no-outgroup") == 0) {
                params.date_with_outgroup = false;
                continue;
            }

            if (strcmp(argv[cnt], "--date-outgroup") == 0) {
                params.date_with_outgroup = true;
                continue;
            }

            if (strcmp(argv[cnt], "--date-ci") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-ci <number_of_replicates>";
                params.date_replicates = convert_int(argv[cnt]);
                if (params.date_replicates < 1)
                    throw "--date-ci must be positive";
                continue;
            }

            if (strcmp(argv[cnt], "--clock-sd") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --clock-sd <standard_dev_of_lognormal_relaxed_lock>";
                params.clock_stddev = convert_double(argv[cnt]);
                if (params.clock_stddev < 0)
                    throw "--clock-sd must be non-negative";
                continue;
            }

            if (strcmp(argv[cnt], "--date-outlier") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-outlier <z_score_for_removing_outlier_nodes>";
                params.date_outlier = convert_double(argv[cnt]);
                if (params.date_outlier < 0)
                    throw "--date-outlier must be non-negative";
                continue;
            }

            if (strcmp(argv[cnt], "--date-debug") == 0) {
                params.date_debug = true;
                continue;
            }
            
            if (strcmp(argv[cnt], "--suppress-list-of-sequences") == 0) {
                params.suppress_list_of_sequences = true;
                continue;
            }

            if (strcmp(argv[cnt], "--suppress-zero-distance") == 0) {
                params.suppress_zero_distance_warnings = true;
                continue;
            }

            if (strcmp(argv[cnt], "--suppress-duplicate-sequence") == 0) {
                params.suppress_duplicate_sequence_warnings = true;
                continue;
            }

            if (strcmp(argv[cnt], "--date-options") == 0) {
                cnt++;
                if (cnt >= argc)
                    throw "Use --date-options <extra_options_for_dating_method>";
                params.dating_options = argv[cnt];
                continue;
            }
            if (arg=="-progress-bar" || arg=="--progress-bar" || arg=="-bar") {
                progress_display::setProgressDisplay(true);
                continue;
            }

            if (argv[cnt][0] == '-') {
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
        }
        // try
        catch (const char *str) {
            if (MPIHelper::getInstance().isMaster())
                outError(str);
            else
                exit(EXIT_SUCCESS);
            //} catch (char *str) {
            //outError(str);
        } catch (string str) {
            if (MPIHelper::getInstance().isMaster())
                outError(str);
            else
                exit(EXIT_SUCCESS);
        } catch (...) {
            string err = "Unknown argument \"";
            err += argv[cnt];
            err += "\"";
            if (MPIHelper::getInstance().isMaster())
                outError(err);
            else
                exit(EXIT_SUCCESS);
        }

    } // for
    if (!params.user_file && !params.aln_file && !params.ngs_file && !params.ngs_mapped_reads && !params.partition_file) {
#ifdef IQ_TREE
        quickStartGuide();
//        usage_iqtree(argv, false);
#else
        usage(argv, false);
#endif
    }

//    if (params.do_au_test)
//        outError("The AU test is temporarily disabled due to numerical issue when bp-RELL=0");

    if (params.root != NULL && params.is_rooted)
        outError("Not allowed to specify both -o <taxon> and -root");
    
    if (params.model_test_and_tree && params.partition_type != BRLEN_OPTIMIZE)
        outError("-mtree not allowed with edge-linked partition model (-spp or -q)");
    
    if (params.do_au_test && params.topotest_replicates == 0)
        outError("For AU test please specify number of bootstrap replicates via -zb option");
    
    if (params.lh_mem_save == LM_MEM_SAVE && params.partition_file)
        outError("-mem option does not work with partition models yet");
    
    if (params.gbo_replicates && params.num_bootstrap_samples)
        outError("UFBoot (-bb) and standard bootstrap (-b) must not be specified together");
    
    if ((params.model_name.find("ONLY") != string::npos || (params.model_name.substr(0,2) == "MF" && params.model_name.substr(0,3) != "MFP")) && (params.gbo_replicates || params.num_bootstrap_samples))
        outError("ModelFinder only cannot be combined with bootstrap analysis");
    
    if (params.num_runs > 1 && !params.treeset_file.empty())
        outError("Can't combine --runs and -z options");
    
    if (params.num_runs > 1 && params.lmap_num_quartets >= 0)
        outError("Can't combine --runs and -lmap options");

    if (params.terrace_analysis && !params.partition_file)
        params.terrace_analysis = false;

    if (params.constraint_tree_file && params.partition_type == TOPO_UNLINKED)
        outError("-g constraint tree option does not work with -S yet.");

    if (params.num_bootstrap_samples && params.partition_type == TOPO_UNLINKED)
        outError("-b bootstrap option does not work with -S yet.");

    if (params.dating_method != "") {
    #ifndef USE_LSD2
        outError("IQ-TREE was not compiled with LSD2 library, rerun cmake with -DUSE_LSD2=ON option");
    #endif
    }

    if (params.date_file.empty()) {
        if (params.date_root.empty() ^ params.date_tip.empty())
            outError("Both --date-root and --date-tip must be provided when --date file is absent");
    }
    
	// Diep:
	if(params.ufboot2corr == true){
		if(params.gbo_replicates <= 0) params.ufboot2corr = false;
		else params.stop_condition = SC_UNSUCCESS_ITERATION;

		params.print_ufboot_trees = 2; // 2017-09-25: fix bug regarding the order of -bb 1000 -bnni -wbt
	}

    if (!params.out_prefix) {
    	if (params.eco_dag_file)
    		params.out_prefix = params.eco_dag_file;
        else if (params.user_file && params.consensus_type == CT_ASSIGN_SUPPORT_EXTENDED)
            params.out_prefix = params.user_file;
        else if (params.partition_file) {
            params.out_prefix = params.partition_file;
            if (params.out_prefix[strlen(params.out_prefix)-1] == '/' || params.out_prefix[strlen(params.out_prefix)-1] == '\\') {
                params.out_prefix[strlen(params.out_prefix)-1] = 0;
            }
        } else if (params.aln_file) {
            params.out_prefix = params.aln_file;
            if (params.out_prefix[strlen(params.out_prefix)-1] == '/' || params.out_prefix[strlen(params.out_prefix)-1] == '\\') {
                params.out_prefix[strlen(params.out_prefix)-1] = 0;
            }
        } else if (params.ngs_file)
            params.out_prefix = params.ngs_file;
        else if (params.ngs_mapped_reads)
            params.out_prefix = params.ngs_mapped_reads;
        else
            params.out_prefix = params.user_file;
    }

    if (params.model_name.find("LINK") != string::npos || params.model_name.find("MERGE") != string::npos)
        if (params.partition_merge == MERGE_NONE)
            params.partition_merge = MERGE_RCLUSTERF;

    //    if (MPIHelper::getInstance().isWorker()) {
    // BUG: setting out_prefix this way cause access to stack, which is cleaned up after returning from this function
//        string newPrefix = string(params.out_prefix) + "."  + NumberToString(MPIHelper::getInstance().getProcessID()) ;
//        params.out_prefix = (char *) newPrefix.c_str();
//    }
    
    if (!params.additional_alignment_files.empty()) {
        params.incremental_method = true;
    }
}

void usage(char* argv[]) {
    printCopyright(cout);
    cout << "Usage: " << argv[0] << " [OPTIONS] <file_name> [<output_file>]" << endl;
    cout << "GENERAL OPTIONS:" << endl;
    cout << "  -hh               Print this help dialog" << endl;
    cout << "  -h                Print help options for phylogenetic inference" << endl;
    cout << "  <file_name>       User tree in NEWICK format or split network in NEXUS format" << endl;
    cout << "  <output_file>     Output file to store results, default is '<file_name>.pda'" << endl;
    cout << "  -k <num_taxa>     Find optimal set of size <num_taxa>" << endl;
    cout << "  -k <min>:<max>    Find optimal sets of size from <min> to <max>" << endl;
    cout << "  -k <min>:<max>:<step>" << endl;
    cout << "                    Find optimal sets of size min, min+step, min+2*step,..." << endl;
    cout << "  -o <taxon>        Root name to compute rooted PD (default: unrooted)" << endl;
    cout << "  -if <file>        File containing taxa to be included into optimal sets" << endl;
    cout << "  -e <file>         File containing branch/split scale and taxa weights" << endl;
    cout << "  -all              Identify all multiple optimal sets" << endl;
    cout << "  -lim <max_limit>  The maximum number of optimal sets for each k if -a is specified" << endl;
    cout << "  -min              Compute minimal sets (default: maximal)" << endl;
    cout << "  -1out             Print taxa sets and scores to separate files" << endl;
    cout << "  -oldout           Print output compatible with version 0.3" << endl;
    cout << "  -v                Verbose mode" << endl;
    cout << endl;
    cout << "OPTIONS FOR PHYLOGENETIC DIVERSITY (PD):" << endl;
    cout << "  -root             Make the tree ROOTED, default is unrooted" << endl;
    cout << "    NOTE: this option and -o <taxon> cannot be both specified" << endl;
    cout << "  -g                Run greedy algorithm only (default: auto)" << endl;
    cout << "  -pr               Run pruning algorithm only (default: auto)" << endl;
    cout << endl;
    /*
    cout << "OPTIONS FOR SPLIT DIVERSITY:" << endl;
    cout << "  -exhaust          Force to use exhaustive search" << endl;
    cout << "    NOTE: by default, the program applies dynamic programming algorithm" << endl;
    cout << "          on circular networks and exhaustive search on general networks" << endl;
    cout << endl;*/
    cout << "OPTIONS FOR BUDGET CONSTRAINTS:" << endl;
    cout << "  -u <file>         File containing total budget and taxa preservation costs" << endl;
    cout << "  -b <budget>       Total budget to conserve taxa" << endl;
    cout << "  -b <min>:<max>    Find all sets with budget from <min> to <max>" << endl;
    cout << "  -b <min>:<max>:<step>" << endl;
    cout << "                    Find optimal sets with budget min, min+step, min+2*step,..." << endl;
    cout << endl;
    cout << "OPTIONS FOR AREA ANALYSIS:" << endl;
    cout << "  -ts <taxa_file>   Compute/maximize PD/SD of areas (combine with -k to maximize)" << endl;
    cout << "  -excl             Compute exclusive PD/SD" << endl;
    cout << "  -endem            Compute endemic PD/SD" << endl;
    cout << "  -compl <areas>    Compute complementary PD/SD given the listed <areas>" << endl;
    cout << endl;

    cout << "OPTIONS FOR VIABILITY CONSTRAINTS:" << endl;
    cout << "  -eco <food_web>   File containing food web matrix" << endl;
    cout << "  -k% <n>           Find optimal set of size relative the total number of taxa" << endl;
    cout << "  -diet <min_diet>  Minimum diet portion (%) to be preserved for each predator" << endl;
    cout << endl;
    //if (!full_command) exit(0);

    cout << "MISCELLANEOUS:" << endl;
    cout << "  -dd <sample_size> Compute PD distribution of random sets of size k" << endl;
    /*
    cout << "  -gbo <sitelh_file> Compute and output the alignment of (normalized)" << endl;
    cout << "                    expected frequencies given in site_ll_file" << endl;
	*/

    //	cout << "  -rep <times>        Repeat algorithm a number of times." << endl;
    //	cout << "  -noout              Print no output file." << endl;
    cout << endl;
    //cout << "HIDDEN OPTIONS: see the source code file pda.cpp::parseArg()" << endl;

    exit(0);
}

void usage_iqtree(char* argv[], bool full_command) {
    printCopyright(cout);
    cout << "Usage: iqtree [-s ALIGNMENT] [-p PARTITION] [-m MODEL] [-t TREE] ..." << endl << endl;
    cout << "GENERAL OPTIONS:" << endl
    << "  -h, --help           Print (more) help usages" << endl
    << "  -s FILE[,...,FILE]   PHYLIP/FASTA/NEXUS/CLUSTAL/MSF alignment file(s)" << endl
    << "  -s DIR               Directory of alignment files" << endl
    << "  --seqtype STRING     BIN, DNA, AA, NT2AA, CODON, MORPH (default: auto-detect)" << endl
    << "  -t FILE|PARS|RAND    Starting tree (default: 99 parsimony and BIONJ)" << endl
    << "  -o TAX[,...,TAX]     Outgroup taxon (list) for writing .treefile" << endl
    << "  --prefix STRING      Prefix for all output files (default: aln/partition)" << endl
    << "  --seed NUM           Random seed number, normally used for debugging purpose" << endl
    << "  --safe               Safe likelihood kernel to avoid numerical underflow" << endl
    << "  --mem NUM[G|M|%]     Maximal RAM usage in GB | MB | %" << endl
    << "  --runs NUM           Number of indepedent runs (default: 1)" << endl
    << "  -v, --verbose        Verbose mode, printing more messages to screen" << endl
    << "  -V, --version        Display version number" << endl
    << "  --quiet              Quiet mode, suppress printing to screen (stdout)" << endl
    << "  -fconst f1,...,fN    Add constant patterns into alignment (N=no. states)" << endl
    << "  --epsilon NUM        Likelihood epsilon for parameter estimate (default 0.01)" << endl
#ifdef _OPENMP
    << "  -T NUM|AUTO          No. cores/threads or AUTO-detect (default: 1)" << endl
    << "  --threads-max NUM    Max number of threads for -T AUTO (default: all cores)" << endl
#endif
    << endl << "CHECKPOINT:" << endl
    << "  --redo               Redo both ModelFinder and tree search" << endl
    << "  --redo-tree          Restore ModelFinder and only redo tree search" << endl
    << "  --undo               Revoke finished run, used when changing some options" << endl
    << "  --cptime NUM         Minimum checkpoint interval (default: 60 sec and adapt)" << endl
    << endl << "PARTITION MODEL:" << endl
    << "  -p FILE|DIR          NEXUS/RAxML partition file or directory with alignments" << endl
    << "                       Edge-linked proportional partition model" << endl
    << "  -q FILE|DIR          Like -p but edge-linked equal partition model " << endl
    << "  -Q FILE|DIR          Like -p but edge-unlinked partition model" << endl
    << "  -S FILE|DIR          Like -p but separate tree inference" << endl
    << "  --subsample NUM      Randomly sub-sample partitions (negative for complement)" << endl
    << "  --subsample-seed NUM Random number seed for --subsample" << endl
    << endl << "LIKELIHOOD/QUARTET MAPPING:" << endl
    << "  --lmap NUM           Number of quartets for likelihood mapping analysis" << endl
    << "  --lmclust FILE       NEXUS file containing clusters for likelihood mapping" << endl
    << "  --quartetlh          Print quartet log-likelihoods to .quartetlh file" << endl
    << endl << "TREE SEARCH ALGORITHM:" << endl
//            << "  -pll                 Use phylogenetic likelihood library (PLL) (default: off)" << endl
    << "  --ninit NUM          Number of initial parsimony trees (default: 100)" << endl
    << "  --ntop NUM           Number of top initial trees (default: 20)" << endl
    << "  --nbest NUM          Number of best trees retained during search (defaut: 5)" << endl
    << "  -n NUM               Fix number of iterations to stop (default: OFF)" << endl
    << "  --nstop NUM          Number of unsuccessful iterations to stop (default: 100)" << endl
    << "  --perturb NUM        Perturbation strength for randomized NNI (default: 0.5)" << endl
    << "  --radius NUM         Radius for parsimony SPR search (default: 6)" << endl
    << "  --allnni             Perform more thorough NNI search (default: OFF)" << endl
    << "  -g FILE              (Multifurcating) topological constraint tree file" << endl
    << "  --fast               Fast search to resemble FastTree" << endl
    << "  --polytomy           Collapse near-zero branches into polytomy" << endl
    << "  --tree-fix           Fix -t tree (no tree search performed)" << endl
    << "  --treels             Write locally optimal trees into .treels file" << endl
    << "  --show-lh            Compute tree likelihood without optimisation" << endl
#ifdef IQTREE_TERRAPHAST
    << "  --terrace            Check if the tree lies on a phylogenetic terrace" << endl
#endif
//            << "  -iqp                 Use the IQP tree perturbation (default: randomized NNI)" << endl
//            << "  -iqpnni              Switch back to the old IQPNNI tree search algorithm" << endl
    << endl << "ULTRAFAST BOOTSTRAP/JACKKNIFE:" << endl
    << "  -B, --ufboot NUM     Replicates for ultrafast bootstrap (>=1000)" << endl
    << "  -J, --ufjack NUM     Replicates for ultrafast jackknife (>=1000)" << endl
    << "  --jack-prop NUM      Subsampling proportion for jackknife (default: 0.5)" << endl
    << "  --sampling STRING    GENE|GENESITE resampling for partitions (default: SITE)" << endl
    << "  --boot-trees         Write bootstrap trees to .ufboot file (default: none)" << endl
    << "  --wbtl               Like --boot-trees but also writing branch lengths" << endl
//            << "  -n <#iterations>     Minimum number of iterations (default: 100)" << endl
    << "  --nmax NUM           Maximum number of iterations (default: 1000)" << endl
    << "  --nstep NUM          Iterations for UFBoot stopping rule (default: 100)" << endl
    << "  --bcor NUM           Minimum correlation coefficient (default: 0.99)" << endl
    << "  --beps NUM           RELL epsilon to break tie (default: 0.5)" << endl
    << "  --bnni               Optimize UFBoot trees by NNI on bootstrap alignment" << endl
    << endl << "NON-PARAMETRIC BOOTSTRAP/JACKKNIFE:" << endl
    << "  -b, --boot NUM       Replicates for bootstrap + ML tree + consensus tree" << endl
    << "  -j, --jack NUM       Replicates for jackknife + ML tree + consensus tree" << endl
    << "  --jack-prop NUM      Subsampling proportion for jackknife (default: 0.5)" << endl
    << "  --bcon NUM           Replicates for bootstrap + consensus tree" << endl
    << "  --bonly NUM          Replicates for bootstrap only" << endl
#ifdef USE_BOOSTER
    << "  --tbe                Transfer bootstrap expectation" << endl
#endif
//            << "  -t <threshold>       Minimum bootstrap support [0...1) for consensus tree" << endl
    << endl << "SINGLE BRANCH TEST:" << endl
    << "  --alrt NUM           Replicates for SH approximate likelihood ratio test" << endl
    << "  --alrt 0             Parametric aLRT test (Anisimova and Gascuel 2006)" << endl
    << "  --abayes             approximate Bayes test (Anisimova et al. 2011)" << endl
    << "  --lbp NUM            Replicates for fast local bootstrap probabilities" << endl
    << endl << "MODEL-FINDER:" << endl
    << "  -m TESTONLY          Standard model selection (like jModelTest, ProtTest)" << endl
    << "  -m TEST              Standard model selection followed by tree inference" << endl
    << "  -m MF                Extended model selection with FreeRate heterogeneity" << endl
    << "  -m MFP               Extended model selection followed by tree inference" << endl
    << "  -m ...+LM            Additionally test Lie Markov models" << endl
    << "  -m ...+LMRY          Additionally test Lie Markov models with RY symmetry" << endl
    << "  -m ...+LMWS          Additionally test Lie Markov models with WS symmetry" << endl
    << "  -m ...+LMMK          Additionally test Lie Markov models with MK symmetry" << endl
    << "  -m ...+LMSS          Additionally test strand-symmetric models" << endl
    << "  --mset STRING        Restrict search to models supported by other programs" << endl
    << "                       (raxml, phyml or mrbayes)" << endl
    << "  --mset STR,...       Comma-separated model list (e.g. -mset WAG,LG,JTT)" << endl
    << "  --msub STRING        Amino-acid model source" << endl
    << "                       (nuclear, mitochondrial, chloroplast or viral)" << endl
    << "  --mfreq STR,...      List of state frequencies" << endl
    << "  --mrate STR,...      List of rate heterogeneity among sites" << endl
    << "                       (e.g. -mrate E,I,G,I+G,R is used for -m MF)" << endl
    << "  --cmin NUM           Min categories for FreeRate model [+R] (default: 2)" << endl
    << "  --cmax NUM           Max categories for FreeRate model [+R] (default: 10)" << endl
    << "  --merit AIC|AICc|BIC  Akaike|Bayesian information criterion (default: BIC)" << endl
//            << "  -msep                Perform model selection and then rate selection" << endl
    << "  --mtree              Perform full tree search for every model" << endl
    << "  --madd STR,...       List of mixture models to consider" << endl
    << "  --mdef FILE          Model definition NEXUS file (see Manual)" << endl
    << "  --modelomatic        Find best codon/protein/DNA models (Whelan et al. 2015)" << endl

    << endl << "PARTITION-FINDER:" << endl
    << "  --merge              Merge partitions to increase model fit" << endl
    << "  --merge greedy|rcluster|rclusterf" << endl
    << "                       Set merging algorithm (default: rclusterf)" << endl
    << "  --merge-model 1|all  Use only 1 or all models for merging (default: 1)" << endl
    << "  --merge-model STR,..." << endl
    << "                       Comma-separated model list for merging" << endl
    << "  --merge-rate 1|all   Use only 1 or all rate heterogeneity (default: 1)" << endl
    << "  --merge-rate STR,..." << endl
    << "                       Comma-separated rate list for merging" << endl
    << "  --rcluster NUM       Percentage of partition pairs for rcluster algorithm" << endl
    << "  --rclusterf NUM      Percentage of partition pairs for rclusterf algorithm" << endl
    << "  --rcluster-max NUM   Max number of partition pairs (default: 10*partitions)" << endl

    << endl << "SUBSTITUTION MODEL:" << endl
    << "  -m STRING            Model name string (e.g. GTR+F+I+G)" << endl
    << "                 DNA:  HKY (default), JC, F81, K2P, K3P, K81uf, TN/TrN, TNef," << endl
    << "                       TIM, TIMef, TVM, TVMef, SYM, GTR, or 6-digit model" << endl
    << "                       specification (e.g., 010010 = HKY)" << endl
    << "             Protein:  LG (default), Poisson, cpREV, mtREV, Dayhoff, mtMAM," << endl
    << "                       JTT, WAG, mtART, mtZOA, VT, rtREV, DCMut, PMB, HIVb," << endl
    << "                       HIVw, JTTDCMut, FLU, Blosum62, GTR20, mtMet, mtVer, mtInv, FLAVI," << endl
    << "			Q.LG, Q.pfam, Q.pfam_gb, Q.bird, Q.mammal, Q.insect, Q.plant, Q.yeast" << endl
    << "     Protein mixture:  C10,...,C60, EX2, EX3, EHO, UL2, UL3, EX_EHO, LG4M, LG4X" << endl
    << "              Binary:  JC2 (default), GTR2" << endl
    << "     Empirical codon:  KOSI07, SCHN05" << endl
    << "   Mechanistic codon:  GY (default), MG, MGK, GY0K, GY1KTS, GY1KTV, GY2K," << endl
    << "                       MG1KTS, MG1KTV, MG2K" << endl
    << "Semi-empirical codon:  XX_YY where XX is empirical and YY is mechanistic model" << endl
    << "      Morphology/SNP:  MK (default), ORDERED, GTR" << endl
    << "      Lie Markov DNA:  1.1, 2.2b, 3.3a, 3.3b, 3.3c, 3.4, 4.4a, 4.4b, 4.5a," << endl
    << "                       4.5b, 5.6a, 5.6b, 5.7a, 5.7b, 5.7c, 5.11a, 5.11b, 5.11c," << endl
    << "                       5.16, 6.6, 6.7a, 6.7b, 6.8a, 6.8b, 6.17a, 6.17b, 8.8," << endl
    << "                       8.10a, 8.10b, 8.16, 8.17, 8.18, 9.20a, 9.20b, 10.12," << endl
    << "                       10.34, 12.12 (optionally prefixed by RY, WS or MK)" << endl
    << "      Non-reversible:  STRSYM (strand symmetric model, equiv. WS6.6)," << endl
    << "                       NONREV, UNREST (unrestricted model, equiv. 12.12)" << endl
    << "           Otherwise:  Name of file containing user-model parameters" << endl
    << endl << "STATE FREQUENCY:" << endl
    << "  -m ...+F             Empirically counted frequencies from alignment" << endl
    << "  -m ...+FO            Optimized frequencies by maximum-likelihood" << endl
    << "  -m ...+FQ            Equal frequencies" << endl
    << "  -m ...+FRY           For DNA, freq(A+G)=1/2=freq(C+T)" << endl
    << "  -m ...+FWS           For DNA, freq(A+T)=1/2=freq(C+G)" << endl
    << "  -m ...+FMK           For DNA, freq(A+C)=1/2=freq(G+T)" << endl
    << "  -m ...+Fabcd         4-digit constraint on ACGT frequency" << endl
    << "                       (e.g. +F1221 means f_A=f_T, f_C=f_G)" << endl
    << "  -m ...+FU            Amino-acid frequencies given protein matrix" << endl
    << "  -m ...+F1x4          Equal NT frequencies over three codon positions" << endl
    << "  -m ...+F3x4          Unequal NT frequencies over three codon positions" << endl

    << endl << "RATE HETEROGENEITY AMONG SITES:" << endl
    << "  -m ...+I             A proportion of invariable sites" << endl
    << "  -m ...+G[n]          Discrete Gamma model with n categories (default n=4)" << endl
    << "  -m ...*G[n]          Discrete Gamma model with unlinked model parameters" << endl
    << "  -m ...+I+G[n]        Invariable sites plus Gamma model with n categories" << endl
    << "  -m ...+R[n]          FreeRate model with n categories (default n=4)" << endl
    << "  -m ...*R[n]          FreeRate model with unlinked model parameters" << endl
    << "  -m ...+I+R[n]        Invariable sites plus FreeRate model with n categories" << endl
    << "  -m ...+Hn            Heterotachy model with n classes" << endl
    << "  -m ...*Hn            Heterotachy model with n classes and unlinked parameters" << endl
    << "  --alpha-min NUM      Min Gamma shape parameter for site rates (default: 0.02)" << endl
    << "  --gamma-median       Median approximation for +G site rates (default: mean)" << endl
    << "  --rate               Write empirical Bayesian site rates to .rate file" << endl
    << "  --mlrate             Write maximum likelihood site rates to .mlrate file" << endl
//            << "  --mhrate             Computing site-specific rates to .mhrate file using" << endl
//            << "                       Meyer & von Haeseler (2003) method" << endl

    << endl << "POLYMORPHISM AWARE MODELS (PoMo):"                                           << endl
    << "  -s FILE              Input counts file (see manual)"                               << endl
    << "  -m ...+P             DNA substitution model (see above) used with PoMo"            << endl
    << "  -m ...+N<POPSIZE>    Virtual population size (default: 9)"                         << endl
// TODO DS: Maybe change default to +WH.
    << "  -m ...+WB|WH|S]      Weighted binomial sampling"       << endl
    << "  -m ...+WH            Weighted hypergeometric sampling" << endl
    << "  -m ...+S             Sampled sampling"              << endl
    << "  -m ...+G[n]          Discrete Gamma rate with n categories (default n=4)"    << endl
// TODO DS: Maybe change default to +WH.

    << endl << "COMPLEX MODELS:" << endl
    << "  -m \"MIX{m1,...,mK}\"  Mixture model with K components" << endl
    << "  -m \"FMIX{f1,...fK}\"  Frequency mixture model with K components" << endl
    << "  --mix-opt            Optimize mixture weights (default: detect)" << endl
    << "  -m ...+ASC           Ascertainment bias correction" << endl
    << "  --tree-freq FILE     Input tree to infer site frequency model" << endl
    << "  --site-freq FILE     Input site frequency model file" << endl
    << "  --freq-max           Posterior maximum instead of mean approximation" << endl

    << endl << "TREE TOPOLOGY TEST:" << endl
    << "  --trees FILE         Set of trees to evaluate log-likelihoods" << endl
    << "  --test NUM           Replicates for topology test" << endl
    << "  --test-weight        Perform weighted KH and SH tests" << endl
    << "  --test-au            Approximately unbiased (AU) test (Shimodaira 2002)" << endl
    << "  --sitelh             Write site log-likelihoods to .sitelh file" << endl

    << endl << "ANCESTRAL STATE RECONSTRUCTION:" << endl
    << "  --ancestral          Ancestral state reconstruction by empirical Bayes" << endl
    << "  --asr-min NUM        Min probability of ancestral state (default: equil freq)" << endl

    << endl << "TEST OF SYMMETRY:" << endl
    << "  --symtest               Perform three tests of symmetry" << endl
    << "  --symtest-only          Do --symtest then exist" << endl
//    << "  --bisymtest             Perform three binomial tests of symmetry" << endl
//    << "  --symtest-perm NUM      Replicates for permutation tests of symmetry" << endl
    << "  --symtest-remove-bad    Do --symtest and remove bad partitions" << endl
    << "  --symtest-remove-good   Do --symtest and remove good partitions" << endl
    << "  --symtest-type MAR|INT  Use MARginal/INTernal test when removing partitions" << endl
    << "  --symtest-pval NUMER    P-value cutoff (default: 0.05)" << endl
    << "  --symtest-keep-zero     Keep NAs in the tests" << endl

    << endl << "CONCORDANCE FACTOR ANALYSIS:" << endl
    << "  -t FILE              Reference tree to assign concordance factor" << endl
    << "  --gcf FILE           Set of source trees for gene concordance factor (gCF)" << endl
    << "  --df-tree            Write discordant trees associated with gDF1" << endl
    << "  --scf NUM            Number of quartets for site concordance factor (sCF)" << endl
    << "  -s FILE              Sequence alignment for --scf" << endl
    << "  -p FILE|DIR          Partition file or directory for --scf" << endl
    << "  --cf-verbose         Write CF per tree/locus to cf.stat_tree/_loci" << endl
    << "  --cf-quartet         Write sCF for all resampled quartets to .cf.quartet" << endl

#ifdef USE_LSD2
    << endl << "TIME TREE RECONSTRUCTION:" << endl
    << "  --date FILE          File containing dates of tips or ancestral nodes" << endl
    << "  --date TAXNAME       Extract dates from taxon names after last '|'" << endl
    << "  --date-tip STRING    Tip dates as a real number or YYYY-MM-DD" << endl
    << "  --date-root STRING   Root date as a real number or YYYY-MM-DD" << endl
    << "  --date-ci NUM        Number of replicates to compute confidence interval" << endl
    << "  --clock-sd NUM       Std-dev for lognormal relaxed clock (default: 0.2)" << endl
    << "  --date-no-outgroup   Exclude outgroup from time tree" << endl
    << "  --date-outlier NUM   Z-score cutoff to remove outlier tips/nodes (e.g. 3)" << endl
    << "  --date-options \"..\"  Extra options passing directly to LSD2" << endl
    << "  --dating STRING      Dating method: LSD for least square dating (default)" << endl
#endif
    << endl;
    

//            << endl << "TEST OF MODEL HOMOGENEITY:" << endl
//            << "  -m WHTEST            Testing model (GTR+G) homogeneity assumption using" << endl
//            << "                       Weiss & von Haeseler (2003) method" << endl
//            << "  -ns <#simulations>   #Simulations to obtain null-distribution (default: 1000)" << endl

    if (full_command)
    cout
        << endl << "CONSENSUS RECONSTRUCTION:" << endl
        << "  -t FILE              Set of input trees for consensus reconstruction" << endl
        << "  --sup-min NUM        Min split support, 0.5 for majority-rule consensus" << endl
        << "                       (default: 0, extended consensus)" << endl
        << "  --burnin NUM         Burnin number of trees to ignore" << endl
        << "  --con-tree           Compute consensus tree to .contree file" << endl
        << "  --con-net            Computing consensus network to .nex file" << endl
        << "  --support FILE       Assign support values into this tree from -t trees" << endl
        //<< "  -sup2 FILE           Like -sup but -t trees can have unequal taxon sets" << endl
        << "  --suptag STRING      Node name (or ALL) to assign tree IDs where node occurs" << endl
        << endl << "TREE DISTANCE BY ROBINSON-FOULDS (RF) METRIC:" << endl
        << "  --tree-dist-all      Compute all-to-all RF distances for -t trees" << endl
        << "  --tree-dist FILE     Compute RF distances between -t trees and this set" << endl
        << "  --tree-dist2 FILE    Like -rf but trees can have unequal taxon sets" << endl
    //            << "  -rf_adj              Computing RF distances of adjacent trees in <treefile>" << endl
    //            << "  -wja                 Write ancestral sequences by joint reconstruction" << endl


        << endl

        << "GENERATING RANDOM TREES:" << endl
        << "  -r NUM               No. taxa for Yule-Harding random tree" << endl
        << "  --rand UNI|CAT|BAL   UNIform | CATerpillar | BALanced random tree" << endl
        //<< "  --rand NET           Random circular split network" << endl
        << "  --rlen NUM NUM NUM   min, mean, and max random branch lengths" << endl

        << endl << "MISCELLANEOUS:" << endl
        << "  --keep-ident         Keep identical sequences (default: remove & finally add)" << endl
        << "  -blfix               Fix branch lengths of user tree passed via -te" << endl
        << "  -blscale             Scale branch lengths of user tree passed via -t" << endl
        << "  -blmin               Min branch length for optimization (default 0.000001)" << endl
        << "  -blmax               Max branch length for optimization (default 100)" << endl
        << "  -wslr                Write site log-likelihoods per rate category" << endl
        << "  -wslm                Write site log-likelihoods per mixture class" << endl
        << "  -wslmr               Write site log-likelihoods per mixture+rate class" << endl
        << "  -wspr                Write site probabilities per rate category" << endl
        << "  -wspm                Write site probabilities per mixture class" << endl
        << "  -wspmr               Write site probabilities per mixture+rate class" << endl
        << "  --partlh             Write partition log-likelihoods to .partlh file" << endl
        << "  --no-outfiles        Suppress printing output files" << endl
        << "  --eigenlib           Use Eigen3 library" << endl
        << "  -alninfo             Print alignment sites statistics to .alninfo" << endl
    //            << "  -d <file>            Reading genetic distances from file (default: JC)" << endl
    //			<< "  -d <outfile>         Calculate the distance matrix inferred from tree" << endl
    //			<< "  -stats <outfile>     Output some statistics about branch lengths" << endl
    //			<< "  -comp <treefile>     Compare tree with each in the input trees" << endl;


        << endl;

    if (full_command) {
        //TODO Print other options here (to be added)
    }

    exit(0);
}

void quickStartGuide() {
    printCopyright(cout);
    cout << "Command-line examples (replace 'iqtree2 ...' by actual path to executable):" << endl << endl
         << "1. Infer maximum-likelihood tree from a sequence alignment (example.phy)" << endl
         << "   with the best-fit model automatically selected by ModelFinder:" << endl
         << "     iqtree2 -s example.phy" << endl << endl
         << "2. Perform ModelFinder without subsequent tree inference:" << endl
         << "     iqtree2 -s example.phy -m MF" << endl
         << "   (use '-m TEST' to resemble jModelTest/ProtTest)" << endl << endl
         << "3. Combine ModelFinder, tree search, ultrafast bootstrap and SH-aLRT test:" << endl
         << "     iqtree2 -s example.phy --alrt 1000 -B 1000" << endl << endl
         << "4. Perform edge-linked proportional partition model (example.nex):" << endl
         << "     iqtree2 -s example.phy -p example.nex" << endl
         << "   (replace '-p' by '-Q' for edge-unlinked model)" << endl << endl
         << "5. Find best partition scheme by possibly merging partitions:" << endl
         << "     iqtree2 -s example.phy -p example.nex -m MF+MERGE" << endl
         << "   (use '-m TESTMERGEONLY' to resemble PartitionFinder)" << endl << endl
         << "6. Find best partition scheme followed by tree inference and bootstrap:" << endl
         << "     iqtree2 -s example.phy -p example.nex -m MFP+MERGE -B 1000" << endl << endl
#ifdef _OPENMP
         << "7. Use 4 CPU cores to speed up computation: add '-T 4' option" << endl << endl
#endif
         << "8. Polymorphism-aware model with HKY nucleotide model and Gamma rate:" << endl
         << "     iqtree2 -s counts_file.cf -m HKY+P+G" << endl << endl
         << "9. PoMo mixture with virtual popsize 5 and weighted binomial sampling:" << endl
         << "     iqtree2 -s counts_file.cf -m \"MIX{HKY+P{EMP},JC+P}+N5+WB\"" << endl << endl
         << "To show all available options: run 'iqtree2 -h'" << endl << endl
         << "Have a look at the tutorial and manual for more information:" << endl
         << "     http://www.iqtree.org" << endl << endl;
    exit(0);
}

InputType detectInputFile(const char *input_file) {

    if (!fileExists(input_file))
        outError("File not found ", input_file);

    try {
        igzstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(input_file);

        unsigned char ch, ch2;
        int count = 0;
        do {
            in >> ch;
        } while (ch <= 32 && !in.eof() && count++ < 20);
        in >> ch2;
        in.close();
        switch (ch) {
            case '#': return IN_NEXUS;
            case '(': return IN_NEWICK;
            case '[': return IN_NEWICK;
            case '>': return IN_FASTA;
            case 'C': if (ch2 == 'L') return IN_CLUSTAL;
                      else if (ch2 == 'O') return IN_COUNTS;
                      else return IN_OTHER;
            case '!': if (ch2 == '!') return IN_MSF; else return IN_OTHER;
            default:
                if (isdigit(ch)) return IN_PHYLIP;
                return IN_OTHER;
        }
    } catch (ios::failure) {
        outError("Cannot read file ", input_file);
    } catch (...) {
        outError("Cannot read file ", input_file);
    }
    return IN_OTHER;
}

bool overwriteFile(char *filename) {
    ifstream infile(filename);
    if (infile.is_open()) {
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
        if (pos >= (signed int) all.length())
            all = "";
        else
            all = all.substr(pos + 1);
    }
}

double logFac(const int num) {
    if (num < 0) return -1.0;
    if (num == 0) return 0.0;
    double ret = 0;
    for (int i = 1; i <= num; i++)
        ret += log((double) i);
    return ret;
}

template <typename I>
I random_element(I begin, I end)
{
    const unsigned long n = std::distance(begin, end);
    const unsigned long divisor = (RAND_MAX + 1) / n;

    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);

    return std::advance(begin, k);
}

template <class T>
inline T quantile(const vector<T>& v, const double q) {
    unsigned int size = v.size();
    if (q <= 0) return *std::min_element(v.begin(), v.end());
    if (q >= 1) return *std::max_element(v.begin(), v.end());
    //double pos = (size - 1) * q;
    //unsigned int ind = (unsigned int)(pos);
    //double delta = pos - ind;
    vector<T> w(size);
    std::copy(v, v.begin() + size, w.begin());
}

#define RAN_STANDARD 1
#define RAN_SPRNG    2
#define RAN_RAND4    3

#define RAN_TYPE 2

#if RAN_TYPE == RAN_STANDARD

int init_random(int seed) {
    srand(seed);
    cout << "(Using rand() - Standard Random Number Generator)" << endl;
    return seed;
}

int finish_random() {
	return 0;
}


#elif RAN_TYPE == RAN_RAND4
/******************************************************************************/
/* random numbers generator  (Numerical recipes)                              */
/******************************************************************************/

/* variable */
long _idum;

/* definitions */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double randomunitintervall()
/* Long period (> 2e18) random number generator. Returns a uniform random
   deviate between 0.0 and 1.0 (exclusive of endpoint values).

   Source:
   Press et al., "Numerical recipes in C", Cambridge University Press, 1992
   (chapter 7 "Random numbers", ran2 random number generator) */ {
    int j;
    long k;
    static long _idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (_idum <= 0) {
        if (-(_idum) < 1)
            _idum = 1;
        else
            _idum = -(_idum);
        _idum2 = (_idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (_idum) / IQ1;
            _idum = IA1 * (_idum - k * IQ1) - k*IR1;
            if (_idum < 0)
                _idum += IM1;
            if (j < NTAB)
                iv[j] = _idum;
        }
        iy = iv[0];
    }
    k = (_idum) / IQ1;
    _idum = IA1 * (_idum - k * IQ1) - k*IR1;
    if (_idum < 0)
        _idum += IM1;
    k = _idum2 / IQ2;
    _idum2 = IA2 * (_idum2 - k * IQ2) - k*IR2;
    if (_idum2 < 0)
        _idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - _idum2;
    iv[j] = _idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
} /* randomunitintervall */

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

int init_random(int seed) /* RAND4 */ {
    //    srand((unsigned) time(NULL));
    //    if (seed < 0)
    // 	seed = rand();
    _idum = -(long) seed;
#ifndef PARALLEL
    cout << "(Using RAND4 Random Number Generator)" << endl;
#else /* PARALLEL */
    {
        int n;
        if (PP_IamMaster) {
            cout << "(Using RAND4 Random Number Generator with leapfrog method)" << endl;
        }
        for (n = 0; n < PP_Myid; n++)
            (void) randomunitintervall();
        if (verbose_mode >= VB_MED) {
            cout << "(" << PP_Myid << ") !!! random seed set to " << seed << ", " << n << " drawn !!!" << endl;
        }
    }
#endif
    return (seed);
} /* initrandom */

int finish_random() {
	return 0;
}
/******************/

#else /* SPRNG */

/******************/

int *randstream;

int init_random(int seed, bool write_info, int** rstream) {
    //    srand((unsigned) time(NULL));
    if (seed < 0)
        seed = make_sprng_seed();
#ifndef PARALLEL
    if (write_info)
    	cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    if (rstream) {
        *rstream = init_sprng(0, 1, seed, SPRNG_DEFAULT); /*init stream*/
    } else {
        randstream = init_sprng(0, 1, seed, SPRNG_DEFAULT); /*init stream*/
        if (verbose_mode >= VB_MED) {
            print_sprng(randstream);
        }
    }
#else /* PARALLEL */
    if (PP_IamMaster && write_info) {
        cout << "(Using SPRNG - Scalable Parallel Random Number Generator)" << endl;
    }
    /* MPI_Bcast(&seed, 1, MPI_UNSIGNED, PP_MyMaster, MPI_COMM_WORLD); */
    if (rstream) {
        *rstream = init_sprng(PP_Myid, PP_NumProcs, seed, SPRNG_DEFAULT); /*initialize stream*/
    } else {
        randstream = init_sprng(PP_Myid, PP_NumProcs, seed, SPRNG_DEFAULT); /*initialize stream*/
        if (verbose_mode >= VB_MED) {
            cout << "(" << PP_Myid << ") !!! random seed set to " << seed << " !!!" << endl;
            print_sprng(randstream);
        }
    }
#endif /* PARALLEL */
    return (seed);
} /* initrandom */

int finish_random(int *rstream) {
    if (rstream)
        return free_sprng(rstream);
    else
        return free_sprng(randstream);
}

#endif /* USE_SPRNG */

/******************/

/* returns a random integer in the range [0; n - 1] */
int random_int(int n, int *rstream) {
    return (int) floor(random_double(rstream) * n);
} /* randominteger */

/* returns a random integer in the range [a; b] */
int random_int(int a, int b) {
	ASSERT(b > a);
	//return a + (RAND_MAX * rand() + rand()) % (b + 1 - a);
	return a + random_int(b - a);
}

double random_double(int *rstream) {
#ifndef FIXEDINTRAND
#ifndef PARALLEL
#if RAN_TYPE == RAN_STANDARD
    return ((double) rand()) / ((double) RAND_MAX + 1);
#elif RAN_TYPE == RAN_SPRNG
    if (rstream)
        return sprng(rstream);
    else
        return sprng(randstream);
#else /* NO_SPRNG */
    return randomunitintervall();
#endif /* NO_SPRNG */
#else /* NOT PARALLEL */
#if RAN_TYPE == RAN_SPRNG
    if (rstream)
        return sprng(rstream);
    else
        return sprng(randstream);
#else /* NO_SPRNG */
    int m;
    for (m = 1; m < PP_NumProcs; m++)
        (void) randomunitintervall();
    PP_randn += (m - 1);
    PP_rand++;
    return randomunitintervall();
#endif /* NO_SPRNG */
#endif /* NOT PARALLEL */
#else /* FIXEDINTRAND */
    cerr << "!!! fixed \"random\" integers for testing purposes !!!" << endl;
    return 0.0;
#endif /* FIXEDINTRAND */

}

void random_resampling(int n, IntVector &sample, int *rstream) {
    sample.resize(n, 0);
    if (Params::getInstance().jackknife_prop == 0.0) {
        // boostrap resampling
        for (int i = 0; i < n; i++) {
            int j = random_int(n, rstream);
            sample[j]++;
        }
    } else {
        // jackknife resampling
        int total = floor((1.0 - Params::getInstance().jackknife_prop)*n);
        if (total <= 0)
            outError("Jackknife sample size is zero");
        // make sure jackknife samples have exacly the same size
        for (int num = 0; num < total; ) {
            for (int i = 0; i < n; i++) if (!sample[i]) {
                if (random_double(rstream) < Params::getInstance().jackknife_prop)
                    continue;
                sample[i] = 1;
                num++;
                if (num >= total)
                    break;
            }
        }
    }
}


/* Following part is taken from ModelTest software */
#define	BIGX            20.0                                 /* max value to represent exp (x) */
#define	LOG_SQRT_PI     0.5723649429247000870717135          /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795          /* 1 / sqrt (pi) */
#define	Z_MAX           6.0                                  /* maximum meaningful z value */
#define	ex(x)           (((x) < -BIGX) ? 0.0 : exp (x))

/************** Normalz: probability of normal z value *********************/

/*
ALGORITHM:	Adapted from a polynomial approximation in:
                        Ibbetson D, Algorithm 209
                        Collected Algorithms of the CACM 1963 p. 616
                Note:
                        This routine has six digit accuracy, so it is only useful for absolute
                        z values < 6.  For z values >= to 6.0, Normalz() returns 0.0.
 */

double Normalz(double z) /*VAR returns cumulative probability from -oo to z VAR normal z value */ {
    double y, x, w;

    if (z == 0.0)
        x = 0.0;
    else {
        y = 0.5 * fabs(z);
        if (y >= (Z_MAX * 0.5))
            x = 1.0;
        else if (y < 1.0) {
            w = y*y;
            x = ((((((((0.000124818987 * w
                    - 0.001075204047) * w + 0.005198775019) * w
                    - 0.019198292004) * w + 0.059054035642) * w
                    - 0.151968751364) * w + 0.319152932694) * w
                    - 0.531923007300) * w + 0.797884560593) * y * 2.0;
        } else {
            y -= 2.0;
            x = (((((((((((((-0.000045255659 * y
                    + 0.000152529290) * y - 0.000019538132) * y
                    - 0.000676904986) * y + 0.001390604284) * y
                    - 0.000794620820) * y - 0.002034254874) * y
                    + 0.006549791214) * y - 0.010557625006) * y
                    + 0.011630447319) * y - 0.009279453341) * y
                    + 0.005353579108) * y - 0.002141268741) * y
                    + 0.000535310849) * y + 0.999936657524;
        }
    }
    return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}


/**************  ChiSquare: probability of chi square value *************/

/*ALGORITHM Compute probability of chi square value.
Adapted from: 	Hill, I. D. and Pike, M. C.  Algorithm 299.Collected Algorithms for the CACM 1967 p. 243
Updated for rounding errors based on remark inACM TOMS June 1985, page 185. Found in Perlman.lib*/

double computePValueChiSquare(double x, int df) /* x: obtained chi-square value,  df: degrees of freedom */ {
    double a, y, s;
    double e, c, z;
    int even; /* true if df is an even number */

    if (x <= 0.0 || df < 1)
        return (1.0);

    y = 1;

    a = 0.5 * x;
    even = (2 * (df / 2)) == df;
    if (df > 1)
        y = ex(-a);
    s = (even ? y : (2.0 * Normalz(-sqrt(x))));
    if (df > 2) {
        x = 0.5 * (df - 1.0);
        z = (even ? 1.0 : 0.5);
        if (a > BIGX) {
            e = (even ? 0.0 : LOG_SQRT_PI);
            c = log(a);
            while (z <= x) {
                e = log(z) + e;
                s += ex(c * z - a - e);
                z += 1.0;
            }
            return (s);
        } else {
            e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
            c = 0.0;
            while (z <= x) {
                e = e * (a / z);
                c = c + e;
                z += 1.0;
            }
            return (c * y + s);
        }
    } else
        return (s);
}

void trimString(string &str) {
    str.erase(0, str.find_first_not_of(" \n\r\t"));
    str.erase(str.find_last_not_of(" \n\r\t")+1);
}



Params& Params::getInstance() {
    static Params instance;
    return instance;
}


int countPhysicalCPUCores() {
    #ifdef _OPENMP
    return omp_get_num_procs();
    #else
    return std::thread::hardware_concurrency();
    #endif
    /*
    uint32_t registers[4];
    unsigned logicalcpucount;
    unsigned physicalcpucount;
#if defined(_WIN32) || defined(WIN32)
    SYSTEM_INFO systeminfo;
    GetSystemInfo( &systeminfo );
    logicalcpucount = systeminfo.dwNumberOfProcessors;
#else
    logicalcpucount = sysconf( _SC_NPROCESSORS_ONLN );
#endif
    if (logicalcpucount < 1) logicalcpucount = 1;
    return logicalcpucount;
    
    if (logicalcpucount % 2 != 0)
        return logicalcpucount;
    __asm__ __volatile__ ("cpuid " :
                          "=a" (registers[0]),
                          "=b" (registers[1]),
                          "=c" (registers[2]),
                          "=d" (registers[3])
                          : "a" (1), "c" (0));

    unsigned CPUFeatureSet = registers[3];
    bool hyperthreading = CPUFeatureSet & (1 << 28);    
    if (hyperthreading){
        physicalcpucount = logicalcpucount / 2;
    } else {
        physicalcpucount = logicalcpucount;
    }
    if (physicalcpucount < 1) physicalcpucount = 1;
    return physicalcpucount;
     */
}

// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0

/** Print a demangled stack backtrace of the caller function to FILE* out. */

#if  !defined(Backtrace_FOUND)

// donothing for WIN32
void print_stacktrace(ostream &out, unsigned int max_frames) {}

#else

void print_stacktrace(ostream &out, unsigned int max_frames)
{
#ifdef _OPENMP
#pragma omp master
{
#endif
    out << "STACK TRACE FOR DEBUGGING:" << endl;

    // storage array for stack trace address data
    void* addrlist[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

//    if (addrlen == 0) {
//        out << "  <empty, possibly corrupt>" << endl;
//        return;
//    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
    {
	char *begin_name = 0, *begin_offset = 0;

	// find parentheses and +address offset surrounding the mangled name:
#ifdef __clang__
      // OSX style stack trace
      for ( char *p = symbollist[i]; *p; ++p )
      {
         if (( *p == '_' ) && ( *(p-1) == ' ' ))
            begin_name = p-1;
         else if ( *p == '+' )
            begin_offset = p-1;
      }

      if ( begin_name && begin_offset && ( begin_name < begin_offset ))
      {
         *begin_name++ = '\0';
         *begin_offset++ = '\0';

         // mangled name is now in [begin_name, begin_offset) and caller
         // offset in [begin_offset, end_offset). now apply
         // __cxa_demangle():
         int status;
         char* ret = abi::__cxa_demangle( begin_name, &funcname[0],
                                          &funcnamesize, &status );
         if ( status == 0 )
         {
            funcname = ret; // use possibly realloc()-ed string
//            out << "  " << symbollist[i] << " : " << funcname << "+"<< begin_offset << endl;
            out << i << "   "  << funcname << endl;
         } else {
            // demangling failed. Output function name as a C function with
            // no arguments.
//             out << "  " << symbollist[i] << " : " << begin_name << "()+"<< begin_offset << endl;
            out << i << "   " << begin_name << "()" << endl;
         }

#else // !DARWIN - but is posix
         // ./module(function+0x15c) [0x8048a6d]
    char *end_offset = 0;
	for (char *p = symbollist[i]; *p; ++p)
	{
	    if (*p == '(')
		begin_name = p;
	    else if (*p == '+')
		begin_offset = p;
	    else if (*p == ')' && begin_offset) {
		end_offset = p;
		break;
	    }
	}

	if (begin_name && begin_offset && end_offset
	    && begin_name < begin_offset)
	{
	    *begin_name++ = '\0';
	    *begin_offset++ = '\0';
	    *end_offset = '\0';

	    // mangled name is now in [begin_name, begin_offset) and caller
	    // offset in [begin_offset, end_offset). now apply
	    // __cxa_demangle():

	    int status;
	    char* ret = abi::__cxa_demangle(begin_name,
					    funcname, &funcnamesize, &status);
	    if (status == 0) {
            funcname = ret; // use possibly realloc()-ed string
//            out << "  " << symbollist[i] << " : " << funcname << "+"<< begin_offset << endl;
            out << i << "   " << funcname << endl;
	    }
	    else {
            // demangling failed. Output function name as a C function with
            // no arguments.
//            out << "  " << symbollist[i] << " : " << begin_name << "()+"<< begin_offset << endl;
            out << i << "   " << begin_name << "()" << endl;
	    }
#endif
	}
	else
	{
	    // couldn't parse the line? print the whole line.
//	    out << i << ". " << symbollist[i] << endl;
	}
    }

    free(funcname);
    free(symbollist);
#ifdef _OPENMP
}
#endif

}

#endif // Backtrace_FOUND

bool memcmpcpy(void * destination, const void * source, size_t num) {
    bool diff = (memcmp(destination, source, num) != 0);
    memcpy(destination, source, num);
    return diff;
}

// Pairing function: see https://en.wikipedia.org/wiki/Pairing_function
int pairInteger(int int1, int int2) {
    if (int1 <= int2) {
        return ((int1 + int2)*(int1 + int2 + 1)/2 + int2);
    } else {
        return ((int1 + int2)*(int1 + int2 + 1)/2 + int1);
    }
}

/*
 * Given a model name, look in it for "+F..." and 
 * determine the StateFreqType. Returns FREQ_EMPIRICAL if
 * unable to find a good +F... specifier
 */
StateFreqType parseStateFreqFromPlusF(string model_name) {
//    StateFreqType freq_type = FREQ_EMPIRICAL;

    // BQM 2017-05-02: change back to FREQ_UNKNOWN to resemble old behavior
    StateFreqType freq_type = FREQ_UNKNOWN;
    size_t plusFPos;
    if (model_name.find("+F1X4") != string::npos)
        freq_type = FREQ_CODON_1x4;
    else if (model_name.find("+F3X4C") != string::npos)
        freq_type = FREQ_CODON_3x4C;
    else if (model_name.find("+F3X4") != string::npos)
        freq_type = FREQ_CODON_3x4;
    else if (model_name.find("+FQ") != string::npos)
        freq_type = FREQ_EQUAL;
    else if (model_name.find("+FO") != string::npos)
        freq_type = FREQ_ESTIMATE;
    else if (model_name.find("+FU") != string::npos)
        freq_type = FREQ_USER_DEFINED;
    else if (model_name.find("+FRY") != string::npos)
        freq_type = FREQ_DNA_RY;
    else if (model_name.find("+FWS") != string::npos)
        freq_type = FREQ_DNA_WS;
    else if (model_name.find("+FMK") != string::npos)
        freq_type = FREQ_DNA_MK;
    else if ((plusFPos = model_name.find("+F")) != string::npos) {
        freq_type = FREQ_EMPIRICAL;
        // Now look for +F#### where #s are digits
        if (model_name.length() > plusFPos+2 && isdigit(model_name[plusFPos+2]))
        try {
            // throws if string is not 4 digits
            freq_type = parseStateFreqDigits(model_name.substr(plusFPos+2,4));
        } catch (const char *str) {
            // +F exists, but can't parse it as anything else
            outError(str);
        }
    }
    return(freq_type);
}

/*
 * Given a string of 4 digits, return a StateFreqType according to
 * equality constraints expressed by those digits.
 * E.g. "1233" constrains pi_G=pi_T (ACGT order, 3rd and 4th equal)
 * which results in FREQ_DNA_2311. "5288" would give the same result.
 */

StateFreqType parseStateFreqDigits(string digits) {
    bool good = true;
    if (digits.length()!=4) {
        good = false;
    } else {
        // Convert digits to canonical form, first occuring digit becomes 1 etc.
        int digit_order[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        int first_found = 0;
        for (int i=0; i<4; i++) {
            int digit = digits[i]-'0';
            if (digit<0 || digit>9) {
	good = false; // found a non-digit
	break;
            }
            if (digit_order[digit]==-1) {
	// haven't seen this digit before
	digit_order[digit]=++first_found;
            }
            // rewrite digit in canonical form
            digits[i] = '0'+digit_order[digit];
        }
        // e.g. if digits was "5288", digit_order will end up as {-1,-1,2,-1,-1,1,-1,-1,3,-1}
    }
    if (!good) throw "Use -f <c | o | u | q | ry | ws | mk | <digit><digit><digit><digit>>";
    StateFreqType freq_type = FREQ_UNKNOWN;
    // Now just exhaustively list all canonical digit possibilities
    if (digits.compare("1111")==0) {
        freq_type = FREQ_EQUAL;
    } else if (digits.compare("1112")==0) {
        freq_type = FREQ_DNA_1112;
    } else if (digits.compare("1121")==0) {
        freq_type = FREQ_DNA_1121;
    } else if (digits.compare("1211")==0) {
        freq_type = FREQ_DNA_1211;
    } else if (digits.compare("1222")==0) {
        freq_type = FREQ_DNA_2111;
    } else if (digits.compare("1122")==0) {
        freq_type = FREQ_DNA_1122;
    } else if (digits.compare("1212")==0) {
        freq_type = FREQ_DNA_1212;
    } else if (digits.compare("1221")==0) {
        freq_type = FREQ_DNA_1221;
    } else if (digits.compare("1123")==0) {
        freq_type = FREQ_DNA_1123;
    } else if (digits.compare("1213")==0) {
        freq_type = FREQ_DNA_1213;
    } else if (digits.compare("1231")==0) {
        freq_type = FREQ_DNA_1231;
    } else if (digits.compare("1223")==0) {
        freq_type = FREQ_DNA_2113;
    } else if (digits.compare("1232")==0) {
        freq_type = FREQ_DNA_2131;
    } else if (digits.compare("1233")==0) {
        freq_type = FREQ_DNA_2311;
    } else if (digits.compare("1234")==0) {
        freq_type = FREQ_ESTIMATE;
    } else
        throw ("Unrecognized canonical digits - Can't happen"); // paranoia is good.
    return freq_type;
}


/*
 * All params in range [0,1] 
 * returns true if base frequencies have changed as a result of this call
 */

bool freqsFromParams(double *freq_vec, double *params, StateFreqType freq_type) {

    // BQM 2017-05-02: Note that only freq for A, C, G are free parameters and stored
    // in params, whereas freq_T is not free and should be handled properly

    double pA, pC, pG, pT; // base freqs
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        return false;
    case FREQ_ESTIMATE:
    // Minh: in code review, please pay extra attention to ensure my treadment of FREQ_ESTIMATE is equivalent to old treatment.
    // BQM: DONE!
        pA=params[0];
        pC=params[1];
        pG=params[2];
        //pT=1-pA-pC-pG;
        pT=freq_vec[3];
        break;
    case FREQ_DNA_RY:
        pA = params[0]/2;
        pG = 0.5-pA;
        pC = params[1]/2;
        pT = 0.5-pC;
        break;
    case FREQ_DNA_WS:
        pA = params[0]/2;
        pT = 0.5-pA;
        pC = params[1]/2;
        pG = 0.5-pC;
        break;
    case FREQ_DNA_MK:
        pA = params[0]/2;
        pC = 0.5-pA;
        pG = params[1]/2;
        pT = 0.5-pG;
        break;
    case FREQ_DNA_1112:
        pA = pC = pG = params[0]/3;
        pT = 1-3*pA;
        break;
    case FREQ_DNA_1121:
        pA = pC = pT = params[0]/3;
        pG = 1-3*pA;
        break;
    case FREQ_DNA_1211:
        pA = pG = pT = params[0]/3;
        pC = 1-3*pA;
        break;
    case FREQ_DNA_2111:
        pC = pG = pT = params[0]/3;
        pA = 1-3*pC;
        break;
    case FREQ_DNA_1122:
        pA = params[0]/2;
        pC = pA;
        pG = 0.5-pA;
        pT = pG;
        break;
    case FREQ_DNA_1212:
        pA = params[0]/2;
        pG = pA;
        pC = 0.5-pA;
        pT = pC;
        break;
    case FREQ_DNA_1221:
        pA = params[0]/2;
        pT = pA;
        pC = 0.5-pA;
        pG = pC;
        break;
    case FREQ_DNA_1123:
        pA = params[0]/2;
        pC = pA;
        pG = params[1]*(1-2*pA);
        pT = 1-pG-2*pA;
        break;
    case FREQ_DNA_1213:
        pA = params[0]/2;
        pG = pA;
        pC = params[1]*(1-2*pA);
        pT = 1-pC-2*pA;
        break;
    case FREQ_DNA_1231:
        pA = params[0]/2;
        pT = pA;
        pC = params[1]*(1-2*pA);
        pG = 1-pC-2*pA;
        break;
    case FREQ_DNA_2113:
        pC = params[0]/2;
        pG = pC;
        pA = params[1]*(1-2*pC);
        pT = 1-pA-2*pC;
        break;
    case FREQ_DNA_2131:
        pC = params[0]/2;
        pT = pC;
        pA = params[1]*(1-2*pC);
        pG = 1-pA-2*pC;
        break;
    case FREQ_DNA_2311:
        pG = params[0]/2;
        pT = pG;
        pA = params[1]*(1-2*pG);
        pC = 1-pA-2*pG;
        break;
    default:
        throw("Unrecognized freq_type in freqsFromParams - can't happen");
    }

    // To MDW, 2017-05-02: please make sure that frequencies are positive!
    // Otherwise, numerical issues will occur.

    bool changed = freq_vec[0]!=pA || freq_vec[1]!=pC || freq_vec[2]!=pG || freq_vec[3]!=pT;
    if (changed) {
        freq_vec[0]=pA;
        freq_vec[1]=pC;
        freq_vec[2]=pG;
        freq_vec[3]=pT;
    }
    return(changed);
}

/*
 * For given freq_type, derives frequency parameters from freq_vec
 * All parameters are in range [0,1] (assuming freq_vec is valid)
 */

void paramsFromFreqs(double *params, double *freq_vec, StateFreqType freq_type) {
    double pA = freq_vec[0]; // These just improve code readability
    double pC = freq_vec[1];
    double pG = freq_vec[2];
//    double pT = freq_vec[3]; // pT is not used below
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        break; // freq_vec never changes
    case FREQ_ESTIMATE:
        params[0]=pA;
        params[1]=pC;
        params[2]=pG;
        break;
    case FREQ_DNA_RY:
        params[0]=2*pA;
        params[1]=2*pC;
        break;
    case FREQ_DNA_WS:
        params[0]=2*pA;
        params[1]=2*pC;
        break;
    case FREQ_DNA_MK:
        params[0]=2*pA;
        params[1]=2*pG;
        break;
    case FREQ_DNA_1112:
        params[0]=3*pA;
        break;
    case FREQ_DNA_1121:
        params[0]=3*pA;
        break;
    case FREQ_DNA_1211:
        params[0]=3*pA;
        break;
    case FREQ_DNA_2111:
        params[0]=3*pC;
        break;
    case FREQ_DNA_1122:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1212:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1221:
        params[0]=2*pA;
        break;
    case FREQ_DNA_1123:
        params[0]=2*pA;
        params[1]=pG/(1-params[0]);
        break;
    case FREQ_DNA_1213:
        params[0]=2*pA;
        params[1]=pC/(1-params[0]);
        break;
    case FREQ_DNA_1231:
        params[0]=2*pA;
        params[1]=pC/(1-params[0]);
        break;
    case FREQ_DNA_2113:
        params[0]=2*pC;
        params[1]=pA/(1-params[0]);
        break;
    case FREQ_DNA_2131:
        params[0]=2*pC;
        params[1]=pA/(1-params[0]);
        break;
    case FREQ_DNA_2311:
        params[0]=2*pG;
        params[1]=pA/(1-params[0]);
        break;
    default:
        throw("Unrecognized freq_type in paramsFromFreqs - can't happen");
    }
}

/* 
 * Given a DNA freq_type and a base frequency vector, alter the
 * base freq vector to conform with the constraints of freq_type
 */
void forceFreqsConform(double *base_freq, StateFreqType freq_type) {
    double pA = base_freq[0]; // These just improve code readability
    double pC = base_freq[1];
    double pG = base_freq[2];
    double pT = base_freq[3];
    double scale;
    switch (freq_type) {
    case FREQ_EQUAL:
        // this was already handled, thus not necessary to check here 
//        base_freq[0] = base_freq[1] = base_freq[2] = base_freq[3] = 0.25;
//        break;
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
    case FREQ_ESTIMATE:
        break; // any base_freq is legal
    case FREQ_DNA_RY:
        scale = 0.5/(pA+pG);
        base_freq[0] = pA*scale;
        base_freq[2] = pG*scale;
        scale = 0.5/(pC+pT);
        base_freq[1] = pC*scale;
        base_freq[3] = pT*scale;
        break;
    case FREQ_DNA_WS:
        scale = 0.5/(pA+pT);
        base_freq[0] = pA*scale;
        base_freq[3] = pT*scale;
        scale = 0.5/(pC+pG);
        base_freq[1] = pC*scale;
        base_freq[2] = pG*scale;
        break;
    case FREQ_DNA_MK:
        scale = 0.5/(pA+pC);
        base_freq[0] = pA*scale;
        base_freq[1] = pC*scale;
        scale = 0.5/(pG+pT);
        base_freq[2] = pG*scale;
        base_freq[3] = pT*scale;
        break;
    case FREQ_DNA_1112:
        base_freq[0]=base_freq[1]=base_freq[2]=(pA+pC+pG)/3;
        break;
    case FREQ_DNA_1121:
        base_freq[0]=base_freq[1]=base_freq[3]=(pA+pC+pT)/3;
        break;
    case FREQ_DNA_1211:
        base_freq[0]=base_freq[2]=base_freq[3]=(pA+pG+pT)/3;
        break;
    case FREQ_DNA_2111:
        base_freq[1]=base_freq[2]=base_freq[3]=(pC+pG+pT)/3;
        break;
    case FREQ_DNA_1122:
        base_freq[0]=base_freq[1]=(pA+pC)/2;
        base_freq[2]=base_freq[3]=(pG+pT)/2;
        break;
    case FREQ_DNA_1212:
        base_freq[0]=base_freq[2]=(pA+pG)/2;
        base_freq[1]=base_freq[3]=(pC+pT)/2;
        break;
    case FREQ_DNA_1221:
        base_freq[0]=base_freq[3]=(pA+pT)/2;
        base_freq[1]=base_freq[2]=(pC+pG)/2;
        break;
    case FREQ_DNA_1123:
        base_freq[0]=base_freq[1]=(pA+pC)/2;
        break;
    case FREQ_DNA_1213:
        base_freq[0]=base_freq[2]=(pA+pG)/2;
        break;
    case FREQ_DNA_1231:
        base_freq[0]=base_freq[3]=(pA+pT)/2;
        break;
    case FREQ_DNA_2113:
        base_freq[1]=base_freq[2]=(pC+pG)/2;
        break;
    case FREQ_DNA_2131:
        base_freq[1]=base_freq[3]=(pC+pT)/2;
        break;
    case FREQ_DNA_2311:
        base_freq[2]=base_freq[3]=(pG+pT)/2;
        break;
    default:
        throw("Unrecognized freq_type in forceFreqsConform - can't happen");
    }
    ASSERT(base_freq[0]>=0 && base_freq[1]>=0 && base_freq[2]>=0 && base_freq[3]>=0 && fabs(base_freq[0]+base_freq[1]+base_freq[2]+base_freq[3]-1)<1e-7);
}

/*
 * For given freq_type, how many parameters are needed to
 * determine frequenc vector?
 * Currently, this is for DNA StateFreqTypes only.
 */

int nFreqParams(StateFreqType freq_type) {
    switch (freq_type) {
    case FREQ_DNA_1112:
    case FREQ_DNA_1121:
    case FREQ_DNA_1211:
    case FREQ_DNA_2111:
    case FREQ_DNA_1122:
    case FREQ_DNA_1212:
    case FREQ_DNA_1221:
        return(1);
    case FREQ_DNA_RY:
    case FREQ_DNA_WS:
    case FREQ_DNA_MK:
    case FREQ_DNA_1123:
    case FREQ_DNA_1213:
    case FREQ_DNA_1231:
    case FREQ_DNA_2113:
    case FREQ_DNA_2131:
    case FREQ_DNA_2311:
        return(2);   
    default:
        return 0; // BQM: don't care about other cases
    }
}

/*
 * For freq_type, and given every base must have frequency >= min_freq, set upper
 * and lower bounds for parameters.
 */
 void setBoundsForFreqType(double *lower_bound, 
                           double *upper_bound, 
                           bool *bound_check, 
                           double min_freq, 
                           StateFreqType freq_type) {
    // Sanity check: if min_freq==0, lower_bound=0 and upper_bound=1 
    // (except FREQ_ESTIMATE, which follows legacy code way of doing things.)
    switch (freq_type) {
    case FREQ_EQUAL:
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
        break; // There are no frequency determining parameters
    case FREQ_DNA_1112:
    case FREQ_DNA_1121:
    case FREQ_DNA_1211:
    case FREQ_DNA_2111:
        // one frequency determining parameter
        lower_bound[0] = 3*min_freq;
        upper_bound[0] = 1-min_freq;
        bound_check[0] = true;
        break;
    case FREQ_DNA_1122:
    case FREQ_DNA_1212:
    case FREQ_DNA_1221:
        // one frequency determining parameter
        lower_bound[0] = 2*min_freq;
        upper_bound[0] = 1-2*min_freq;
        bound_check[0] = true;
        break;
    case FREQ_DNA_RY:
    case FREQ_DNA_WS:
    case FREQ_DNA_MK:
        // two frequency determining parameters
        lower_bound[0] = lower_bound[1] = 2*min_freq;
        upper_bound[0] = upper_bound[1] = 1-2*min_freq;
        bound_check[0] = bound_check[1] = true;
	break;
    case FREQ_DNA_1123:
    case FREQ_DNA_1213:
    case FREQ_DNA_1231:
    case FREQ_DNA_2113:
    case FREQ_DNA_2131:
    case FREQ_DNA_2311:
        // two frequency determining parameters
        lower_bound[0] = 2*min_freq;
        upper_bound[0] = 1-2*min_freq;
	lower_bound[1] = min_freq/(1-2*min_freq);
        upper_bound[1] = (1-3*min_freq)/(1-2*min_freq);
        bound_check[0] = bound_check[1] = true;
	break;
        /* NOTE:
	 * upper_bound[1] and lower_bound[1] are not perfect. Some in-bounds parameters
         * will give base freqs for '2' or '3' base below minimum. This is
         * the best that can be done without passing min_freq to freqsFromParams
         */
    case FREQ_ESTIMATE:
        lower_bound[0] = lower_bound[1] = lower_bound[2] = min_freq;
        upper_bound[0] = upper_bound[1] = upper_bound[2] = 1;
        bound_check[0] = bound_check[1] = bound_check[2] = false;
        break;
    default:
        throw("Unrecognized freq_type in setBoundsForFreqType - can't happen");
    }
}
 
double binomial_coefficient_log(unsigned int N, unsigned int n) {
  static DoubleVector logv;
  if (logv.size() <= 0) {
    logv.push_back(0.0);
    logv.push_back(0.0);
  }
  if (n < N-n)
    n = N-n;
  if (n==0)
    return 0.0;
  if (N >= logv.size()) {
    for (unsigned int i = logv.size(); i <= N; i++)
      logv.push_back(log((double) i));
  }
  double binom_log = 0.0;
  for (unsigned int i = n+1; i <= N; i++)
    binom_log += logv[i] - logv[i-n];
  return binom_log;
}

double binomial_dist(unsigned int k, unsigned int N, double p) {
  double binom_log = binomial_coefficient_log(N, k);
  double res_log = binom_log + log(p)*k + log(1-p)*(N-k);
  return exp(res_log);
}

double hypergeometric_dist(unsigned int k, unsigned int n, unsigned int K, unsigned int N) {
  if (n > N)
    outError("Invalid parameters for hypergeometric distribution.");
  if (k > K || (n-k) > (N-K))
    return 0.0;
  double num_successes_log = binomial_coefficient_log(K, k);
  double num_failures_log = binomial_coefficient_log(N-K, n-k);
  double num_total_log = binomial_coefficient_log(N,n);
  return exp(num_successes_log + num_failures_log - num_total_log);
}

// Calculate the Frobenius norm of an N x N matrix M (flattened, rows
// concatenated) and linearly scaled by SCALE.
 double frob_norm(double m[], int n, double scale) {
   double sum = 0;
   for (int i = 0; i < n; i++) {
     for (int j = 0; j < n; j++) {
       sum += m[i*n + j] * m[i*n + j] * scale * scale;
     }
   }
   return sqrt(sum);
 }
