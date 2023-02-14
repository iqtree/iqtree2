#pragma once
#ifndef utils_argument_h
#define utils_argument_h

#include <string>  //for std::string
#include <sstream> //for std::stringstream
#include <map>     //for std::map

/**
 * @brief A formal argument which might appear in a command-line
 * @note  This is an abstract class, with a virtual dtor, so that
 *        maps to pointers to instances of Argument can be used.  
 *        Every subclass must overwrite the accept member function.
 */
class Argument {
public:
    std::string name;
    explicit Argument(const char* arg_name);
    virtual ~Argument() = default;
    /**
     * @brief          Accept a recognized argument (found in the command-line)
     *                 (if the argument takes a parameter - or parameters - 
     *                 this entails doing something with it; parsing it, possibly
     *                 validating it.  The details vary between subclasses).
     * @param arg      The argument that matched the name
     * @param nextArg  The next argument
     * @param argv     The command-line argument array
     * @param argc     The number of command line arguments
     * @param argNum   Index into the command-line argument array
     * @param problems A reference to a stringstream, to which problems
     *                 are to be appended (followed by linefeeds).
     *                 If there is a problem (like, say, an 
     *                 argument that should be followed by a numeric 
     *                 string, that is not, a description of the
     *                 problem, and a line-feed should be appended.
     */
    virtual void accept(const std::string& arg, const std::string& nextArg, 
                        char* argv[], int argc, int &argNum,
                        std::stringstream& problems) = 0;
};

class DoubleArgument: public Argument {
    std::string description;
    double& dbl_var;
public: 
    typedef Argument super;
    /**
     * @brief An argument that should be followed by a string that represents
     *        a double precision number, to be assigned to a double variable.
     * @param arg_name - the name of the argument
     * @param desc     - a description of the argument (to use in error messages)
     * @param var      - a reference to a double (a variable), to be updated
     *                   with the double precision number represented by the
     *                   following argument, if an arg_name argument appears 
     *                   in the command-line.
     * @note  there isn't any check that the argument doesn't appear more than
     *        once.  if it does, the *later* setting will override the *earlier*.
     */
    DoubleArgument(const char* arg_name, const char* desc, double& var);
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, int &argNum, 
                std::stringstream& problems) override;
};

class IntArgument: public Argument {
    std::string description;
    int& int_var;
public:
    typedef Argument super;
    /**
     * @brief An argument that should be followed by a string that represents
     *        an integer, that is to be assigned to an int variable.
     * @param arg_name - the name of the argument
     * @param desc     - a description of the argument (to use in error messages)
     * @param var      - a reference to a int (a variable), to be updated
     *                   with the integer represented by the following argument,
     *                   if an arg_name argument appears in the command-line.
     * @note  there isn't any check that the argument doesn't appear more than
     *        once.  if it does, the *later* setting will override the *earlier*.
     */
    IntArgument(const char* arg_name, const char* desc, int& var);
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, int &argNum, 
                std::stringstream& problems) override;
};

class StringArgument: public Argument {
public:
    typedef Argument super;
    std::string  description;
    std::string& mapped_to;
    /**
     * @brief An argument that should be followed by a string, that is to 
     *        be assigned to std::string variable.
     * @param argument_name        - the name of the argument
     * @param argument_description - a description of the argument 
     *                               (to use in error messages)
     * @param variable - a reference to a std::string (a variable), to be 
     *                   updated with the value of the following argument,
     *                   if an arg_name argument appears in the command-line.
     * @note  there isn't any check that the argument doesn't appear more than
     *        once.  if it does, the *later* setting will override the *earlier*.
     */
    StringArgument(const char* argument_name, 
                   const char* argument_description, 
                   std::string& variable);
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, 
                int &argNum, std::stringstream& problems) override;
};

/**
 * @brief A switch argument (the classic switch argument is -v for verbose)
 *        If the argument is provided, a variable will be set (perhaps
 *        to true, perhaps to false).
 */
class SwitchArgument: public Argument {
private:
    bool& switch_var;
    bool  switch_setting;
public:
    typedef Argument super;
    /**
     * @brief Construct a Switch argument, that will, if a particular
     *        argument is passed in a command-line, set a boolean variable
     *        to true, or to false.
     * @param arg_name - the name of the argument
     * @param var      - the variable to set (perhaps to true, perhaps to false)
     *                   if the argument is found in the command-line
     * @param setting  - what to set var to if an arg_name argument appears
     */
    SwitchArgument(const char* arg_name, bool& var, bool setting);
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, int &argNum,
                std::stringstream& problems) override;
};

class ArgumentMap: std::map<std::string, Argument*> {
public:
    /**
     * @brief  Add an argument to the this map (return a reference to *this).
     * @param  arg - a pointer to the argument
     * @return ArgumentMap& a reference to *this
     * @note   Calling this member function hands ownership of arg to
     *         this ArgumentMap.  When this ArgumentMap goes out of scope,
     *         all the Argument instances it owns the pointers to will be
     *         deleted.
     */
    ArgumentMap& operator << (Argument* arg);
    ~ArgumentMap();
    /**
     * @brief  Given a name, return the pointer, to the Argument, with that
     *         name, owned by this ArgumentMap, if there is one, or null if
     *         there is not.
     * @param  name - the name to look for
     * @return Argument* either a pointer to the matching Argument, owned
     *         by this ArgumentMap, or nullptr (if there isn't one).
     */
    Argument* findByName(const std::string& name);
};

#endif
