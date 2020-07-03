//
//  bionj2.h
//  alignment
//
//  Created by James Barbetti on 18/6/20.
//

#ifndef bionj2_h
#define bionj2_h

#include <string>
#include <algorithm>
class BIONJ2
{
public:
    void constructTree
        ( const std::string &distanceMatrixFilePath
        , const std::string & newickTreeFilePath);
    void constructTreeRapid
        ( const std::string &distanceMatrixFilePath
        , const std::string & newickTreeFilePath);
};

#endif /* bionj2_hpp */
