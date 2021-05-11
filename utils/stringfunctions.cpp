//
//  stringfunctions.cpp
//  Created by James Barbetti on 05-Feb-2021.
//  (But most functions written by others, notably Minh Bui and Tung)
//  Unless tagged otherwise, functions here came from tools.cpp.
//
#include "stringfunctions.h"
#include <sstream>   //for std::stringstream
#include <math.h>    //for HUGE_VALL symbol and abs() function
#include <algorithm> //for std::transform (needed on Windows)
#include <cstring>   //for std::strchr

int convert_int(const char *str) {
    char *endptr;
    int i = strtol(str, &endptr, 10);

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        std::string err = "Expecting integer, but found \"";
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
        std::string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    end_pos = static_cast<int>(endptr - str);
    return i;
}

void convert_int_vec(const char *str, IntVector &vec) {
    const char *beginptr = str;
    char* endptr;
    vec.clear();
    do {
        int i = strtol(beginptr, &endptr, 10);

        if ((i == 0 && endptr == beginptr) || abs(i) == HUGE_VALL) {
            std::string err = "Expecting integer, but found \"";
            err += beginptr;
            err += "\" instead";
            throw err;
        }
        vec.push_back(i);
        if (*endptr == ',') {
            ++endptr;
        }
        beginptr = endptr;
    } while (*endptr != 0);
}


int64_t convert_int64(const char *str) {
    char *endptr;
    int64_t i = (int64_t)strtoll(str, &endptr, 10); // casted because 'long long' may be larger than int64_t

    if ((i == 0 && endptr == str) || abs(i) == HUGE_VALL || *endptr != 0) {
        std::string err = "Expecting large integer , but found \"";
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
        std::string err = "Expecting large integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    end_pos = static_cast<int>(endptr - str);
    return i;
}


double convert_double(const char *str) {
    char *endptr;
    double d = strtod(str, &endptr);
    if ((d == 0.0 && endptr == str) || fabs(d) == HUGE_VALF || *endptr != 0) {
        std::string err = "Expecting floating-point number, but found \"";
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
        std::string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    end_pos = static_cast<int>(endptr - str);
    return d;
}

void convert_double_vec(const char *str, DoubleVector &vec, char separator) {
    const char *beginptr = str;
    char* endptr;
    vec.clear();
    do {
        double d = strtod(beginptr, &endptr);

        if ((d == 0.0 && endptr == beginptr) || fabs(d) == HUGE_VALF) {
            std::string err = "Expecting floating-point number, but found \"";
            err += beginptr;
            err += "\" instead";
            throw err;
        }
        vec.push_back(d);
        if (*endptr == separator) {
            ++endptr;
        }
        beginptr = endptr;
    } while (*endptr != 0);
}

std::string convert_time(const double sec) {
    int sec_int = (int) floor(sec);
    int secs = sec_int % 60;
    int mins = (sec_int % 3600) / 60;
    int hours = sec_int / 3600;
    std::stringstream ss;
    ss << hours << "h:" << mins << "m:" << secs << "s";
    return ss.str();
}

void convert_range(const char *str, int &lower, int &upper, int &step_size) {
    char *endptr;

    // parse the lower bound of the range
    int d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL || (*endptr != 0 && *endptr != ':')) {
        std::string err = "Expecting integer, but found \"";
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
        std::string err = "Expecting integer, but found \"";
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
        std::string err = "Expecting integer, but found \"";
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
        std::string err = "Expecting floating-point number, but found \"";
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
        std::string err = "Expecting floating-point number, but found \"";
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
        std::string err = "Expecting floating-point number, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    step_size = d;
}

void convert_string_vec(const char *str, StrVector &vec, char separator) {
    const char *beginptr = str, *endptr;
    vec.clear();
    std::string elem;
    do {
        endptr = std::strchr(beginptr, separator);
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

//From Tung

std::string convertIntToString(int number) {
    std::stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

std::string convertInt64ToString(int64_t number) {
    std::stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

std::string convertDoubleToString(double number) {
    std::stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

bool iEquals(const std::string& a, const std::string& b)
{
    auto sz = a.size();
    if (b.size() != sz) {
        return false;
    }
    for (size_t i = 0; i < sz; ++i) {
        if (tolower(a[i]) != tolower(b[i])) {
            return false;
        }
    }
    return true;
}

std::string string_to_lower(const std::string& input_string) {
    std::string answer = input_string;
    std::transform(answer.begin(), answer.end(), answer.begin(),
                   []( char c){ return std::tolower(c); });
    return answer;
}

std::string string_to_lower(const char* input) {
    std::string answer = input;
    std::transform(answer.begin(), answer.end(), answer.begin(),
                   []( char c){ return std::tolower(c); });
    return answer;
}

std::string string_to_upper(const std::string& input_string) {
    std::string answer = input_string;
    std::transform(answer.begin(), answer.end(), answer.begin(),
                   []( char c){ return std::toupper(c); });
    return answer;
}

std::string string_to_upper(const char* input) {
    std::string answer = input;
    std::transform(answer.begin(), answer.end(), answer.begin(),
                   []( char c){ return std::toupper(c); });
    return answer;
}

std::string next_argument(int argc, char* argv[], const char* desc, int& cnt ) {
    ++cnt;
    if (cnt >= argc) {
        std::string problem = std::string("Use ") + argv[cnt-1] + " <" + desc + ">";
        throw problem;
    }
    return argv[cnt];
}

//These three functions moved here from model/modelinfofromyamlfile.cpp, 30-Apr-2021:

bool startsWith(const std::string& s, const char* front) {
    auto frontLen = strlen(front);
    return (s.substr(0, frontLen) == front);
}

bool endsWith(const std::string& s, const char* suffix) {
    auto suffixLen = strlen(suffix);
    if (s.length() < suffixLen) {
        return false;
    }
    return s.substr(s.length() - suffixLen, suffixLen) == suffix;
}

bool contains(const std::string &s, const char* pattern) {
    return s.find(pattern) != std::string::npos;
}

bool contains(const char* s, const char* pattern) {
    return strstr(s, pattern) != nullptr;
}

bool is_string_all_digits(const std::string& s) {
    for (auto ch: s) {
        if ( ch < '0' || '9' < ch ) {
            return false;
        }
    }
    return true;
}

bool is_string_all_digits(const char* s) {
    for (; *s!='\0'; ++s) {
        auto ch = *s;
        if ( ch < '0' || '9' < ch ) {
            return false;
        }
    }
    return true;
}

StrVector split_string(const std::string& splitme, const char* withme) {
    StrVector answer;
    if (splitme.empty()) {
        return answer;
    }
    size_t separator_len = strlen(withme);
    if (separator_len==0) {
        answer.push_back(splitme);
        return answer;
    }
    size_t last_cut = 0L;
    for ( size_t next_cut = splitme.find(withme, last_cut);
          next_cut != std::string::npos;
          last_cut = next_cut + separator_len,
          next_cut = splitme.find(withme, last_cut)) {
        std::string next_string = splitme.substr(last_cut, next_cut-last_cut);
        answer.push_back(next_string);
    }
    std::string last_string = splitme.substr(last_cut, splitme.length()-last_cut);
    answer.push_back( last_string );
    return answer;
}

StrVector split_string(const std::string& splitme, const std::string& withme) {
    return split_string(splitme, withme.c_str());
}

