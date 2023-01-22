//
//  argument.cpp - classes for parsing command-line parameters:
//
//      Argument
//      +--IntArgument
//      +--StringArgument
//      +--SwitchArgument
//
//      ArgumentMap
//
//  Copyright (C) 2021, James Barbetti.
//
//  LICENSE:
//* This program is free software; you can redistribute it and/or modify
//* it under the terms of the GNU General Public License as published by
//* the Free Software Foundation; either version 2 of the License, or
//* (at your option) any later version.
//*
//* This program is distributed in the hope that it will be useful,
//* but WITHOUT ANY WARRANTY; without even the implied warranty of
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//* GNU General Public License for more details.
//*
//* You should have received a copy of the GNU General Public License
//* along with this program; if not, write to the
//* Free Software Foundation, Inc.,
//* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//

#include "argument.h"

Argument::Argument(const char* arg_name): name(arg_name) {}

DoubleArgument::DoubleArgument(const char* arg_name, const char* desc, 
                               double& var) 
    : super(arg_name), description(desc), dbl_var(var) { }
void DoubleArgument::accept(const std::string& arg, 
                            const std::string& nextArg, 
                            char* argv[], int argc, int &argNum, 
                            std::stringstream& problems) {
    ++argNum;
    if (argNum==argc) {
        problems << name << " should be followed by " << description << ".\n";
    }
    dbl_var = atof(nextArg.c_str());
}

IntArgument::IntArgument(const char* arg_name, const char* desc, int& var) 
    : super(arg_name), description(desc), int_var(var) { }
void IntArgument::accept(const std::string& arg, const std::string& nextArg, 
            char* argv[], int argc, int &argNum, 
            std::stringstream& problems) {
    ++argNum;
    if (argNum==argc) {
        problems << name << " should be followed by " << description << ".\n";
    }
    int_var = atoi(nextArg.c_str());
}

StringArgument::StringArgument(const char* argument_name, 
                               const char* argument_description, 
                               std::string& variable) 
        : super(argument_name), description(argument_description),
          mapped_to(variable) { }
          
void StringArgument::accept(const std::string& arg, const std::string& nextArg, 
                            char* argv[], int argc, 
                            int &argNum, std::stringstream& problems) {
    ++argNum;
    if (argNum==argc) {
        problems << name << " should be followed by " << description << "\n";
    }
    mapped_to = nextArg;
}

SwitchArgument::SwitchArgument(const char* arg_name, bool& var, bool setting)
        : super(arg_name), switch_var(var), switch_setting(setting) { }
void SwitchArgument::accept(const std::string& arg, const std::string& nextArg, 
            char* argv[], int argc, int &argNum,
            std::stringstream& problems) {
    switch_var = switch_setting;
}

/**
 * @brief  Add an argument, indicated by a pointer, to
 *         this argument, and hand over ownership of it.
 * @param  arg - pointer to the argument
 * @return ArgumentMap& a reference to *this.
 * @note   If there is already an argument, x, in the map, with the same name as
 *         arg, delete it, and replace it.
 */
ArgumentMap& ArgumentMap::operator << (Argument* arg) {
    Argument* already_there = findByName(arg->name);
    if (already_there!=nullptr) {
        erase(find(arg->name)); //remove it from the map, and
        delete already_there;   //delete it (this ArgumentMap owns it)
    }
    insert(value_type(arg->name, arg));
    return *this;
}

/**
 * @brief Destroy the ArgumentMap
 * @note  Deletes all of the Argument instance that this ArgumentMap owns.
 */
ArgumentMap::~ArgumentMap() {
    for (auto it=begin(); it!=end(); ++it) {
        delete it->second;
    }
    clear();
}

/**
 * @brief  Find an argument, owned by the ArgumentMap, with the specified
 *         name, if there is one, and return a pointer to it.
 * @param  name - the name to look for
 * @return Argument* pointer to the argument, if this ArgumentMap 
 *         owns one, with that name, or nullptr, if not.
 */
Argument* ArgumentMap::findByName(const std::string& name) {
    auto it = find(name);
    if (it == end()) {
        return nullptr;
    }
    return it->second;
}
