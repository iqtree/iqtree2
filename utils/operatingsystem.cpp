//
//  operatingsystem.cpp
//  alignment
//
//  Created by James Barbetti on 29/7/20.
//

#include "operatingsystem.h"
#include <sstream>

std::string getOSName() {
    std::stringstream out;
    #if defined _WIN32 || defined WIN32 || defined WIN64
        out << "Windows";
    #elif defined __APPLE__ || defined __MACH__
        out << "Mac OS X";
    #elif defined __linux__
        out << "Linux";
    #elif defined __unix__ || defined __unix
        out << "Unix";
    #else
        out << "Unknown Platform";
    #endif
    out << " " << 8*sizeof(void*) << "-bit";
    return out.str();
}

