#pragma once
#ifndef utils_argument_h
#define utils_argument_h

#include <string>  //for std::string
#include <sstream> //for std::stringstream
#include <map>     //for std::map

class Argument {
public:
    std::string name;
    explicit Argument(const char* arg_name);
    virtual ~Argument() = default;
    virtual void accept(const std::string& arg, const std::string& nextArg, 
                        char* argv[], int argc, int &argNum,
                        std::stringstream& problems) = 0;
};

class DoubleArgument: public Argument {
    std::string description;
    double& dbl_var;
public: 
    typedef Argument super;
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
    StringArgument(const char* argument_name, 
                   const char* argument_description, 
                   std::string& variable);
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, 
                int &argNum, std::stringstream& problems) override;
};

class SwitchArgument: public Argument {
private:
    bool& switch_var;
    bool  switch_setting;
public:
    typedef Argument super;
    SwitchArgument(const char* arg_name, bool& var, bool setting);
    void accept(const std::string& arg, const std::string& nextArg, 
                char* argv[], int argc, int &argNum,
                std::stringstream& problems) override;
};

class ArgumentMap: std::map<std::string, Argument*> {
public:
    ArgumentMap& operator << (Argument* arg);
    ~ArgumentMap();
    Argument* findByName(const std::string& name);
};

#endif
