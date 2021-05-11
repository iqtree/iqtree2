//
//  stringfunctions.hpp
//  Created by James Barbetti on 5/2/21.
//  but, unless otherwise tagged, functions here came from tools.h
//

#ifndef stringfunctions_h
#define stringfunctions_h

#include <string>
#include <sstream>
#include <algorithm>     //for std::reverse
#include "vectortypes.h" //for IntVector, DoubleVector, StrVector

/**
        convert string to int, with error checking
        @param str original string
        @return the number
 */
int convert_int(const char *str);

/**
       convert string to int, with error checking (but not throwing on error)
       @param str original string
       @param defaultValue value to return if the string isn't numeric
       @return the number
*/
int convert_int_nothrow(const char* str, int defaultValue) throw();

/**
    convert string to int64, with error checking
    @param str original string
    @return the number
 */
int64_t convert_int64(const char *str);

/**
        convert string to int, with error checking
        @param str original string
        @param end_pos end position
        @return the number
 */
int convert_int(const char *str, int &end_pos);

/**
        convert comma-separated string to integer vector, with error checking
        @param str original string with integers separated by comma
        @param vec (OUT) integer vector
 */
void convert_int_vec(const char *str, IntVector &vec);

/**
        convert string to int64_t, with error checking
        @param str original string
        @return the number
 */
int64_t convert_int64(const char *str);

/**
        convert string to int64_t, with error checking
        @param str original string
        @param end_pos end position
        @return the number
 */
int64_t convert_int64(const char *str, int &end_pos);

/**
        convert string to double, with error checking
        @param str original string
        @return the double
 */
double convert_double(const char *str);

/**
        convert string to double, with error checking
        @param str original string
        @param end_pos end position
        @return the double
 */
double convert_double(const char *str, int &end_pos);

/**
       convert string to double, with error checking (but not throwing on error)
       @param str original string
       @param defaultValue value to return if the string isn't numeric
       @return the number
*/
double convert_double_nothrow(const char* str, double defaultValue) throw();


/**
        convert comma-separated string to integer vector, with error checking
        @param str original string with integers separated by comma
        @param vec (OUT) integer vector
        @param separator char separating elements
 */
void convert_double_vec(const char *str, DoubleVector &vec, char separator = ',');

/**
 * Convert seconds to hour, minute, second
 * @param sec
 * @return string represent hour, minute, second
 */
std::string convert_time(const double sec);


/**
        convert a string to to range lower:upper:step_size with error checking
        @param str original string
        @param lower (OUT) lower bound of the range
        @param upper (OUT) upper bound of the range
        @param step_size (OUT) step size of the range
 */
void convert_range(const char *str, int &lower, int &upper, int &step_size);

/**
        convert a string to to range lower:upper:step_size with error checking
        @param str original string
        @param lower (OUT) lower bound of the range
        @param upper (OUT) upper bound of the range
        @param step_size (OUT) step size of the range
 */
void convert_range(const char *str, double &lower, double &upper, double &step_size);

void convert_string_vec(const char *str, StrVector &str_vec, char separator = ',');


/**
 * convert int to string
 * @param int
 * @return string
 */
std::string convertIntToString(int number);
std::string convertInt64ToString(int64_t number);

std::string convertDoubleToString(double number);

/**
 case-insensitive comparison between two strings
 @return true if two strings are equal.
 */
bool iEquals(const std::string& a, const std::string& b);
 
std::string string_to_lower(const std::string& input_string); //added 30-Mar-2021

std::string string_to_lower(const char* input); //was only in tools.cpp

std::string string_to_upper(const std::string& input_string); //added 18-Feb-2021

std::string string_to_upper(const char* input); //was only in tools.cpp

std::string next_argument(int argc, char* argv[],
                          const char* desc, int& cnt );
                            //was only in tools.cpp

template <typename T>
std::string NumberToString ( T Number )
{
    std::ostringstream ss;
    ss << Number;
    return ss.str();
}

template <typename T>
T StringToNumber ( const std::string &Text )
{
    std::istringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

//These three functions moved here from model/modelinfofromyamlfile.cpp, 30-Apr-2021:
bool startsWith(const std::string& s, const char* front);
bool endsWith  (const std::string& s, const char* suffix);
bool contains  (const std::string& s, const char* pattern);
bool contains  (const char*        s, const char* pattern);

//These functions added, 04-May-2021
bool is_string_all_digits(const std::string& s);
bool is_string_all_digits(const char* s);

//These functions added, 11-May-2021
StrVector split_string(const std::string& splitme, const std::string& withme);
StrVector split_string(const std::string& splitme, const char* withme);

#endif /* stringfunctions_hpp */
