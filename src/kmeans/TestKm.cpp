// BEWARE: BETA VERSION
// --------------------
//
// A framework for testing k-means and k-means++ against data sets loaded from a file or generated
// at random. The actual k-means code is all contained within KMeans.h. This file is devoted
// entirely to testing.
//
// Author: David Arthur (darthur@gmail.com), 2009

#include "KMeans.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
using namespace std;

// Logging
static vector<ostream*> gLogOutputs[3];
#define LOG(verbosity, text) {                                             \
  vector<ostream*> &outputs = gLogOutputs[verbosity];                      \
  if (outputs.size() > 0) {                                                \
    ostringstream string_stream;                                           \
    string_stream << text;                                                 \
    for (int i = 0; i < (int)outputs.size(); i++)                          \
      *(outputs[i]) << string_stream.str();                                \
  }                                                                        \
}
void ClearLogging() {
  for (int i = 0; i < 3; i++)
    gLogOutputs[i].clear();
  ClearKMeansLogging();
}
void AddLogging(std::ostream *out, int verbosity) {
  if (verbosity == 0)
    return;
  for (int i = 0; i < verbosity; i++)
    gLogOutputs[i].push_back(out);
  if (verbosity > 1)
    AddKMeansLogging(out, verbosity == 3);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Utilities for loading text files
// ================================

// Quits with an error message if condition == false, assuming we are running unit tests.
void UnitTestAssert(bool condition, const string &description) {
  if (condition)
    return;
  LOG(0, "ERROR: Unit test failed:" << endl);
  LOG(0, description << "!" << endl);
  LOG(0, "(Note: unit test errors could *possibly* just be precision problems.)" << endl << endl);
  exit(-1);
}

// Quits with an error message if condition == false, assuming we are loading a text file.
void FileAssert(bool condition, const string &filename, const string &error) {
  if (condition)
    return;
  LOG(0, "ERROR: Could not load file \"" << filename << "\":" << endl);
  LOG(0, error << "!" << endl << endl);
  exit(-1);
}

// Quits with an error message if we ran out of memory.
void MemoryAssertFail(const string &task) {
  LOG(0, "ERROR: Ran out of memory while " << task << "!" << endl << endl);
  exit(-1);
}

// Returns all tokens on the next non-empty line in the given input stream. Text after // comments
// is ignored, and tokens are converted to lower case. True is returned if we do not hit
// end-of-file.
bool TokenizeNextLine(ifstream &in, vector<string> *tokens) {
  string line;

  // Read until end-of-file
  tokens->clear();
  while (getline(in, line)) {
    istringstream line_in(line);
    string token;

    // Read all tokens on this line
    while (line_in >> token) {
      if (int(token.size()) >= 2 && token[0] == '/' && token[1] == '/')
        break;
      for (int i = 0; i < (int)token.size(); i++)
        if (token[i] >= 'A' && token[i] <= 'Z')
          token[i] += 'a' - 'A';
      tokens->push_back(token);
    }

    // If there were valid tokens, return
    if (int(tokens->size()) > 0)
      return true;
  }

  // Return end-of-file
  return false;
}

// If s begins with prefix, sets suffix to the rest of s, and returns true.
// Otherwise returns false.
bool GetSuffix(const string &s, const string &prefix, string *suffix) {
  if (s.size() >= prefix.size() && s.substr(0, prefix.size()) == prefix) {
    *suffix = s.substr(prefix.size());
    return true;
  } else {
    return false;
  }
}

// Returns whether a string appears to be representing a valid numerical type. Non-integers can
// optionally be allowed.
bool IsValidNumericType(const string &s, bool allow_floats) {
  // Handle the initial minus sign if it exists
  int start_i = 0;
  if (s[0] == '-')
    start_i++;

  // Handle leading zeroes
  if (s.size() == start_i)
    return false;
  if (s[start_i] == '0' && int(s.size()) > start_i+1 && s[start_i+1] != '.')
    return false;
  if (s == "-0")
    return false;

  // Check intermediate characters
  for (int i = start_i; i < (int)s.size(); i++) {
    if (s[i] == '.' && allow_floats) {
      allow_floats = false;
      continue;
    } else if (s[i] < '0' || s[i] > '9') {
      return false;
    }
  }
  if (s[int(s.size())-1] == '.')
    return false;

  // No mistakes found
  return true;
}

// Attempts to convert s into a boolean. If s is not valid, it quits with an error message.
bool FileParseBool(const string &s, const string &filename, const string &bool_name) {
  if (s == "true")
    return true;
  else if (s == "false")
    return false;
  LOG(0, "ERROR: Could not read " << bool_name << " while loading file \""
      << filename << "\":" << endl);
  LOG(0, "\"" << s << "\" is not a boolean value!" << endl << endl);
  exit(-1);
}

// Attempts to convert s into an integer in the range [lower_bound, upper_bound]. If s is not
// valid, it quits with an error message.
int FileParseInt(const string &s, int lower_bound, int upper_bound, const string &filename,
                 const string &int_name) {
  bool is_good = true;
  if (is_good) is_good = IsValidNumericType(s, false);

  // Check against lower and upper bound, accounting for overflow
  if (is_good) {
    int max_length = (s[0] == '-'? 11 : 10);
    is_good = (int(s.size()) <= max_length);
  }
  if (is_good) {
    int value = atoi(s.c_str());
    is_good = ((value%10) == (s[int(s.size())-1]-'0'));
    if (is_good && value >= lower_bound && value <= upper_bound)
      return value;
  }

  // Output the error
  LOG(0, "ERROR: Could not read " << int_name << " while loading file \""
      << filename << "\":" << endl);
  LOG(0, "\"" << s << "\" is not an integer between " << lower_bound
      << " and " << upper_bound << "!" << endl << endl);
  exit(-1);
}

// Attempts to convert s into a Scalar. If s is not valid, it quits with an error message.
Scalar FileParseScalar(const string &s, const string &filename, const string &scalar_name) {
  if (IsValidNumericType(s, true))
    return atof(s.c_str());
  LOG(0, "ERROR: Could not read " << scalar_name << " while loading file \""
      << filename << "\":" << endl);
  LOG(0, "\"" << s << "\" is not a real number!" << endl << endl);
  exit(-1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Load a RunSpec, which details what tests should we do
// =====================================================

struct RunSpec {
  RunSpec(): cout_log_level(2), run_km(true), run_kmpp(true), run_unit_tests(false),
             restarts(20) {}

  // Logging info
  int cout_log_level;
  vector<string> log_filename;
  vector<int> log_file_level;

  // Which algorithms should we run?
  bool run_km;
  bool run_kmpp;

  // Should we run a unit test on the kd-tree that is used by both k-means and k-means++?
  bool run_unit_tests;

  // What values of k should we use, and how many times should we restart k-means per data set?
  int restarts;
  vector<int> k_values;

  // What generated data sets should we test on?
  struct GeneratedDataSpec {
    int n, k, d;
    Scalar r, R;
  };
  vector<bool> is_data_generated;
  vector<GeneratedDataSpec> generated_data_specs;
  vector<string> loaded_data_files;
};

// Loads a run spec from the given file and returns the result.
RunSpec LoadRunSpec(const string &filename) {
  ifstream in(filename.c_str());
  FileAssert(!in.fail(), filename, "File not found");
  RunSpec result;

  vector<string> tokens;
  while (TokenizeNextLine(in, &tokens)) {
    for (int i = 0; i < int(tokens.size()); i++) {
      string suffix;

      // Logging
      if (GetSuffix(tokens[i], "coutloglevel=", &suffix)) {
        result.cout_log_level = FileParseInt(suffix, 0, 3, filename, "CoutLogLevel");
      } else if (tokens[i] == "addlog:") {
        FileAssert(i+1 < int(tokens.size()) && GetSuffix(tokens[i+1], "file=", &suffix) &&
                   suffix != "", filename,
                   "Expected token \"File=...\" in log spec");
        result.log_filename.push_back(suffix);
        FileAssert(i+2 < int(tokens.size()) && GetSuffix(tokens[i+2], "level=", &suffix) &&
                   suffix != "", filename,
                   "Expected token \"Level=...\" in log spec");
        result.log_file_level.push_back(FileParseInt(suffix, 1, 3, filename, "file log level"));
        i+=2;
      }

      // Algorithms
      else if (GetSuffix(tokens[i], "runkm=", &suffix)) {
        result.run_km = FileParseBool(suffix, filename, "RunKm");
      } else if (GetSuffix(tokens[i], "runkmpp=", &suffix)) {
        result.run_kmpp = FileParseBool(suffix, filename, "RunKmpp");
      }

      // Unit tests
      else if (GetSuffix(tokens[i], "unittests=", &suffix)) {
        result.run_unit_tests = FileParseBool(suffix, filename, "UnitTests");
      }

      // Restarts
      else if (GetSuffix(tokens[i], "restarts=", &suffix)) {
        result.restarts = FileParseInt(suffix, 1, 1000, filename, "Restarts");
      }

      // k
      else if (GetSuffix(tokens[i], "addk=", &suffix)) {
        int k = FileParseInt(suffix, 2, 1000, filename, "k");
        result.k_values.push_back(k);
      }

      // Generated data sets
      else if (tokens[i] == "addgenerateddata:") {
        RunSpec::GeneratedDataSpec spec;
        FileAssert(i+1 < int(tokens.size()) && GetSuffix(tokens[i+1], "n=", &suffix), filename,
                   "Expected token \"n=...\" in generated data spec");
        spec.n = FileParseInt(suffix, 1, 1000000, filename, "generated n");
        FileAssert(i+2 < int(tokens.size()) && GetSuffix(tokens[i+2], "k=", &suffix), filename,
                   "Expected token \"k=...\" in generated data spec");
        spec.k = FileParseInt(suffix, 1, 1000000, filename, "generated k");
        FileAssert(i+3 < int(tokens.size()) && GetSuffix(tokens[i+3], "d=", &suffix), filename,
                   "Expected token \"d=...\" in generated data spec");
        spec.d = FileParseInt(suffix, 1, 1000, filename, "generated d");
        FileAssert(i+4 < int(tokens.size()) && GetSuffix(tokens[i+4], "r1=", &suffix), filename,
                   "Expected token \"r=...\" in generated data spec");
        spec.r = FileParseScalar(suffix, filename, "generated r1");
        FileAssert(i+5 < int(tokens.size()) && GetSuffix(tokens[i+5], "r2=", &suffix), filename,
                   "Expected token \"R=...\" in generated data spec");
        spec.R = FileParseScalar(suffix, filename, "generated r2");
        result.is_data_generated.push_back(true);
        result.generated_data_specs.push_back(spec);
        i += 5;
      }

      // Loaded data sets
      else if (GetSuffix(tokens[i], "addloadeddata=", &suffix)) {
        FileAssert(suffix != "", filename, "\"" + suffix + "\" is not a filename");
        result.is_data_generated.push_back(false);
        result.loaded_data_files.push_back(suffix);
      }

      // Unrecognized token
      else {
        FileAssert(false, filename, "Unrecognized token: \"" + tokens[i] + "\"");
      }
    }
  }

  // Sanity check the spec
  FileAssert(result.cout_log_level > 0 || int(result.log_file_level.size()) > 0, filename,
             "No logging of any sort");
  FileAssert(int(result.k_values.size()) > 0, filename, "No k-values added");
  FileAssert(result.run_km || result.run_kmpp || result.run_unit_tests, filename,
             "No algorithms are requested");
  FileAssert(int(result.is_data_generated.size()) > 0, filename, "No data sets added");
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Load data sets from a file
// ==========================

// Attempts to load a set of data points from a file
void LoadPointsFromFile(const string &filename, int *n, int *d, Scalar **points) {
  *n = -1;
  *d = -1;
  *points = 0;
  int cur_point = 0;

  // Open the file
  LOG(1, "Loading data file: " << filename << "..." << endl);
  ifstream in(filename.c_str());
  FileAssert(!in.fail(), filename, "File not found");

  // Read each line
  vector<string> tokens;
  while (TokenizeNextLine(in, &tokens)) {
    // Read each word on the line
    int cur_dim = 0;
    for (int i = 0; i < int(tokens.size()); i++) {
      string suffix;

      // n
      if (GetSuffix(tokens[i], "n=", &suffix)) {
        FileAssert(*n == -1, filename, "Redefinition of n");
        *n = FileParseInt(suffix, 1, 1000000, filename, "n");
        cur_dim = -1;
      }

      // d
      else if (GetSuffix(tokens[i], "d=", &suffix)) {
        FileAssert(*d == -1, filename, "Redefinition of d");
        *d = FileParseInt(suffix, 1, 1000000, filename, "d");
        cur_dim = -1;
      }

      // Coordinates
      else if (IsValidNumericType(tokens[i], true)) {
        FileAssert(cur_dim != -1, filename, "Coordinates on the same line as n, d definitions");
        FileAssert(*n != -1 && *d != -1, filename, "Coordinates found before n, d definitions");
        FileAssert(cur_point != *n, filename, "More than n points");
        FileAssert(cur_dim != *d, filename, "Point has more than d dimensions");
        (*points)[cur_point * (*d) + cur_dim] = (Scalar)atof(tokens[i].c_str());
        cur_dim++;
      }

      // Unrecognized token
      else {
        FileAssert(false, filename, "Unrecognized token: \"" + tokens[i] + "\"");
      }
    }

    // If the line was a point, check the number of dimensions and update our counters
    if (cur_dim > 0) {
      FileAssert(cur_dim == *d, filename, "Point has fewer than d dimensions");
      cur_point++;
      if ((cur_point % 10000 == 0) || (cur_point == *n)) {
        LOG(1, cur_point << "/" << (*n) << " points loaded..." << endl);
      }
    }

    // If n and d are now given, allocate points
    if (*n >= 0 && *d >= 0 && *points == 0) {
      *points = (Scalar*)malloc((*n) * (*d) * sizeof(Scalar));
      if (points == 0)
        MemoryAssertFail("reading data points from \"" + filename + "\"");
    }
  }

  // final check
  FileAssert(*n != -1 && *d != -1, filename, "No definition for n, d found");
  FileAssert(cur_point == *n, filename, "Fewer than n points");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Create random data sets that admit a good clustering
// ====================================================

// Returns a number chosen from the 1-dimensional normal distribution with the given mean and
// standard deviation.
Scalar GetNormalDist(Scalar mean, Scalar std_dev) {
  // First choose a number from a distribution with mean 0, standard deviation 1, then translate
  // and scale it to get the right mean and standard deviation. The algorithm for the first part
  // uses the polar form of the Box-Muller transformation and is described in more detail here:
  // http://www.taygeta.com/random/gaussian.html
  while (1) {
    Scalar x = -1 + 2 * Scalar(rand()) / RAND_MAX;
    Scalar y = -1 + 2 * Scalar(rand()) / RAND_MAX;
    Scalar r_squared = x*x + y*y;
    if (r_squared <= 1 && r_squared > 0) {
      Scalar uniform_sample = sqrt(-2 * log(r_squared) / r_squared) * x;
      return mean + uniform_sample * std_dev;
    }
  }
}

// Generates n points in R^d as follows:
//   - First choose k centers according to a normal distribution with standard deviation R.
//   - For each data point, choose a center uniformly at random for it to be associated with.
//     Set its positioned to be that center, perturbed by a normal distribution with standard
//     deviation r.
// The resulting set of points is stored in points.
void ChooseClusterablePoints(int n, int k, int d, Scalar r, Scalar R, Scalar **points) {
  LOG(1, "Generating random data points..." << endl);

  // Allocate memory
  *points = (Scalar*)malloc(n * d * sizeof(Scalar));
  Scalar *centers = (Scalar*)malloc(k * d * sizeof(Scalar));
  if (points == 0 || centers == 0)
    MemoryAssertFail("generating random data points");

  // Choose the centers
  for (int i = 0; i < k * d; i++)
    centers[i] = GetNormalDist(0, R);

  // Choose the data points
  Scalar opt_ub = 0;
  for (int i = 0; i < n; i++) {
    int c = GetRandom(k);
    for (int j = 0; j < d; j++) {
      Scalar offset = GetNormalDist(0, r);
      opt_ub += offset*offset;
      (*points)[i*d + j] = centers[c*d + j] + offset;
    }
  }
  LOG(1, "Min cost for k=" << k << " is at most " << opt_ub << endl);

  // Clean up and return
  free(centers);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Test the implementation of KmTree.h against a simpler but unoptimized implementation to make
// sure it is correct. Unlike the rest of this file, this is not for evaluating k-means++ per se -
// it is a sanity check that the optimizations in KmTree.h are all correct.
// ===============================================================================================

#include "KmTree.h"

// Returns the amount of error we can expect from precision issues while comparing x1 and x2.
Scalar GetTestThreshold(Scalar x1, Scalar x2) {
  Scalar m = 1;
  if (x1 > m) m = x1;
  if (x2 > m) m = x2;
  return m * Scalar(1e-8);
}

// Returns whether x1 == x2 within a small tolerance, used to account for precision error.
bool TestScalarsEq(Scalar x1, Scalar x2) {
  Scalar diff = (Scalar)fabs(x1 - x2);
  return diff <= GetTestThreshold(x1, x2);
}

// Returns whether x1 >= x2 within a small tolerance, used to account for precision error.
bool TestScalarsGe(Scalar x1, Scalar x2) {
  const Scalar kEpsilon = Scalar(1e-5);
  return x1 >= x2 - GetTestThreshold(x1, x2);
}

// Executes a k-means step with a KmTree on the given point-set with the given set of centers, and
// tests if it is correct by comparing results with the naive O(nkd) implementation. The cost is
// returned, which is used to determine when k-means has finished.
//
// Generates an assertion failure if there is an error.
Scalar TestKMeansStep(const KmTree &tree, int n, int k, int d, Scalar *points, Scalar *centers) {
  // Allocate memory
  int *assignment = (int*)malloc(n * sizeof(int));
  Scalar *old_centers = (Scalar*)malloc(k * d * sizeof(Scalar));
  Scalar *new_sums = (Scalar*)calloc(k * d, sizeof(Scalar));
  int *new_counts = (int*)calloc(k, sizeof(int));
  Scalar *bad_center = PointAllocate(d);
  if (assignment == 0 || old_centers == 0 || new_sums == 0 || new_counts == 0 || bad_center == 0)
    MemoryAssertFail("running unit test");
  memset(bad_center, 0xff, d * sizeof(Scalar));
  memcpy(old_centers, centers, k * d * sizeof(Scalar));

  // Run fancy k-means
  Scalar fancy_cost = tree.DoKMeansStep(k, centers, assignment);

  // Test the assignments and build the correct aggregate data
  Scalar correct_cost = 0;
  for (int i = 0; i < n; i++) {
    Scalar fancy_dist_sq = PointDistSq(points + i*d, old_centers + assignment[i]*d, d);
    for (int j = 0; j < k; j++)
    if (memcmp(old_centers + j*d, bad_center, d*sizeof(Scalar)) != 0) {
      Scalar dist_sq = PointDistSq(points + i*d, old_centers + j*d, d);
      UnitTestAssert(TestScalarsGe(dist_sq, fancy_dist_sq),
                     "k-means assigned point to the wrong cluster");
    }

    // Note that the cost is measured from the OLD centers, not from the new centers
    correct_cost += fancy_dist_sq;
    PointAdd(new_sums + assignment[i]*d, points + i*d, d);
    new_counts[assignment[i]]++;
  }

  // Test the costs
  UnitTestAssert(TestScalarsEq(correct_cost, fancy_cost),
                 "k-means calculated the cost function incorrectly");

  // Test the centers
  for (int i = 0; i < k; i++) {
    bool fancy_is_void = (memcmp(centers + i*d, bad_center, d*sizeof(Scalar)) == 0);
    bool correct_is_void = (new_counts[i] == 0);
    UnitTestAssert(fancy_is_void == correct_is_void,
                   "k-means failed to correctly mark whether a center was being used");
    if (!fancy_is_void) {
      PointScale(new_sums + i*d, Scalar(1) / new_counts[i], d);
      for (int j = 0; j < d; j++) {
        UnitTestAssert(TestScalarsEq(new_sums[i*d + j], centers[i*d + j]),
                       "k-means failed to set a center correctly");
      }
    }
  }

  // Free memory
  PointFree(bad_center);
  free(new_counts);
  free(new_sums);
  free(old_centers);
  free(assignment);
  return fancy_cost;
}

// Executes k-means multiple times to completion with a KmTree on the given point-set with centers
// chosen at random, and tests if it is correct by comparing results with the naive O(nkd)
// implementation.
//
// Generates an assertion failure if there is an error.
void TestKMeans(int n, int k, int d, Scalar *points, int attempts) {
  const Scalar kEpsilon = Scalar(1e-8);  // Used to determine when to terminate k-means
  LOG(1, "Unit-testing k-means implementation..." << endl);

  // Initialization
  KmTree tree(n, d, points);
  Scalar *centers = (Scalar*)malloc(sizeof(Scalar)*k*d);
  int *unused_centers = (int*)malloc(sizeof(int)*n);
  if (centers == 0 || unused_centers == 0)
    MemoryAssertFail("running unit test");

  for (int attempt = 0; attempt < attempts; attempt++) {
    // Choose centers uniformly at random
    for (int i = 0; i < n; i++)
      unused_centers[i] = i;
    int num_unused_centers = n;
    for (int i = 0; i < k; i++) {
      int j = GetRandom(num_unused_centers--);
      memcpy(centers + i*d, points + unused_centers[j]*d, d*sizeof(Scalar));
      unused_centers[j] = unused_centers[num_unused_centers];
    }
    
    // Run k-means until it stabilizes
    Scalar old_cost;
    bool is_done = false;
    for (int iteration = 0; !is_done; iteration++) {
      Scalar new_cost = TestKMeansStep(tree, n, k, d, points, centers);
      is_done = (iteration > 0 && new_cost >= (1 - kEpsilon) * old_cost);
      old_cost = new_cost;
      LOG(2, "Unit-test completed iteration " << (iteration+1) << " of run " << (attempt+1) << "..." << endl);
    }
  }
  LOG(1, "Unit test passed!" << endl);

  // Clean up and return
  free(unused_centers);
  free(centers);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Main
// ====

// Executes k-means or k-means++ with the given parameters and returns the cost and time.
struct TestResult {
  Scalar cost;
  double time;
};
TestResult DoTest(bool kmpp, int n, int k, int d, Scalar *points, int restarts) {
  TestResult result;
  double start_time = double(clock()) / CLOCKS_PER_SEC;
  srand(0);
  if (kmpp)
    result.cost = RunKMeansPlusPlus(n, k, d, points, restarts, 0, 0);
  else
    result.cost = RunKMeans(n, k, d, points, restarts, 0, 0);
  result.time = double(clock()) / CLOCKS_PER_SEC - start_time;
  return result;
}

int main(int argnum, char **argv) {
  string spec_filename;
  vector<string> data_set_descriptions;

  // Get the input spec filename
  if (argnum == 1) {
    cout << "Enter the test spec filename: ";
    cin >> spec_filename;
  } else if (argnum == 2) {
    spec_filename = argv[1];
  } else {
    cout << "Usage is TestKm <file>, where <file> is an optional parameter that "
         << "specifies the test spec filename." << endl << endl;
    exit(-1);
  }

  // Load the spec
  AddLogging(&cout, 1);
  RunSpec run_spec = LoadRunSpec(spec_filename);
  ClearLogging();
  AddLogging(&cout, run_spec.cout_log_level);
  for (int i = 0; i < int(run_spec.log_filename.size()); i++)
    AddLogging(new ofstream(run_spec.log_filename[i].c_str()), run_spec.log_file_level[i]);

  // Allocate our results
  int num_algs = 0, km_alg_index = 0, kmpp_alg_index = 0;
  if (run_spec.run_km) {
    km_alg_index = num_algs;
    num_algs++;
  } if (run_spec.run_kmpp) {
    kmpp_alg_index = num_algs;
    num_algs++;
  }
  vector<vector<vector<TestResult> > > results(int(run_spec.is_data_generated.size()),
    vector<vector<TestResult> >(int(run_spec.k_values.size()),
      vector<TestResult>(num_algs)));

  // Handle each data set
  int generated_i = 0, loaded_i = 0;
  for (int i = 0; i < int(run_spec.is_data_generated.size()); i++) {
    int n, d;
    Scalar *points;

    // Get the data
    ostringstream data_set_description;
    if (run_spec.is_data_generated[i]) {
      n = run_spec.generated_data_specs[generated_i].n;
      d = run_spec.generated_data_specs[generated_i].d;
      srand(0);
      ChooseClusterablePoints(n, run_spec.generated_data_specs[generated_i].k, d,
                              run_spec.generated_data_specs[generated_i].r,
                              run_spec.generated_data_specs[generated_i].R, &points);
      data_set_description << "Generated (n=" << n << " d=" << d
        << " k=" << run_spec.generated_data_specs[generated_i].k
        << " r=" << run_spec.generated_data_specs[generated_i].r
        << " R=" << run_spec.generated_data_specs[generated_i].R << ")";
      generated_i++;
    } else {
      LoadPointsFromFile(run_spec.loaded_data_files[loaded_i], &n, &d, &points);
      data_set_description << "Loaded (n=" << n << " d=" << d
        << " file=" << run_spec.loaded_data_files[loaded_i] << ")";
      loaded_i++;
    }
    data_set_descriptions.push_back(data_set_description.str());

    // Loop through each k
    for (int ki = 0; ki < int(run_spec.k_values.size()); ki++) {
      int k = run_spec.k_values[ki];
      LOG(1, "Testing with n = " << n << ", d = " << d << ", k = " << k << "..." << endl);

      // Run the unit test
      if (run_spec.run_unit_tests)
        TestKMeans(n, k, d, points, run_spec.restarts);

      // Run k-means and k-means++
      if (run_spec.run_km)
        results[i][ki][km_alg_index] = DoTest(false, n, k, d, points, run_spec.restarts);
      if (run_spec.run_kmpp)
        results[i][ki][kmpp_alg_index] = DoTest(true, n, k, d, points, run_spec.restarts);
    }

    // Clear the memory
    free(points);
  }

  // Output the results header
  LOG(0, "Final results summary" << endl);
  LOG(0, "=====================" << endl << endl);

  LOG(0, "Data sets:" << endl);
  for (int i = 0; i < int(data_set_descriptions.size()); i++)
    LOG(0, (i+1) << ": " << data_set_descriptions[i] << endl);
  LOG(0, endl);

  // Output the column heads
  LOG(0, "Algorithm\tk\t|\t");
  for (int i = 0; i < int(data_set_descriptions.size()); i++)
    LOG(0, "Cost #" << (i+1) << "\t");
  LOG(0, "|\t");
  for (int i = 0; i < int(data_set_descriptions.size()); i++)
    LOG(0, "Time #" << (i+1) << "\t");
  LOG(0, endl);

  // Output the results
  for (int alg_i = 0; alg_i < num_algs; alg_i++) {
    for (int ki = 0; ki < int(run_spec.k_values.size()); ki++) {
      // Output the algorithm
      if (alg_i == km_alg_index)
        LOG(0, "km\t");
      if (alg_i == kmpp_alg_index)
        LOG(0, "km++\t");

      // Output k
      LOG(0, run_spec.k_values[ki] << "\t|\t");

      // Output costs
      for (int i = 0; i < int(data_set_descriptions.size()); i++)
        LOG(0, results[i][ki][alg_i].cost << "\t");
      LOG(0, "|\t");

      // Output times
      for (int i = 0; i < int(data_set_descriptions.size()); i++)
        LOG(0, results[i][ki][alg_i].time << "\t");
      LOG(0, endl);
    }
  }
  LOG(0, endl);

  return 0;
}
