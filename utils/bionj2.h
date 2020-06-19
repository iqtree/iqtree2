//
//  bionj2.hpp
//  alignment
//
//  Created by James Barbetti on 18/6/20.
//

#ifndef bionj2_hpp
#define bionj2_hpp

#include <string>

class BIONJ2
{
public:
    void constructTree ( const std::string &distanceMatrixFilePath, const std::string & newickTreeFilePath);
};

#endif /* bionj2_hpp */
