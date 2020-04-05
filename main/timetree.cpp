/*
 * timetree.cpp
 * Interface to call dating method incl. LSD2
 *  Created on: Apr 4, 2020
 *      Author: minh
 */

#include "timetree.h"
#include "lsd2/src/lsd.h"

#ifdef USE_LSD2
void runLSD2(PhyloTree *tree) {
    cout << "Building time tree with least-square dating (LSD)...";
    cout << endl;
    string filename = (string)Params::getInstance().out_prefix + ".lsd";
    tree->printTree(filename.c_str());
    StrVector arg = {"lsd", "-i", filename, "-r", "a", "-s", convertIntToString(tree->getAlnNSite()), "-c"};
    int argc = arg.size();
    char *argv[argc];
    for (int i = 0; i < argc; i++)
        argv[i] = (char*)arg[i].c_str();
    std::copy(arg.begin(), arg.end(), std::ostream_iterator<string>(std::cout, " "));
    lsd_main(argc, argv);
    cout << "LSD results written to:" << endl;
    cout << "  LSD report:                  " << filename << ".result" << endl;
    cout << "  Time tree in nexus format:   " << filename << ".result.nexus" << endl;
    cout << "  Time tree with dates:        " << filename << ".result.date.nexus" << endl;
    cout << "  Time tree in newick format:  " << filename << ".result.nwk" << endl;
    cout << endl;
}
#endif

void doTimeTree(PhyloTree *tree) {
#ifdef USE_LSD2
    if (Params::getInstance().dating_method == "LSD") {
        runLSD2(tree);
        return;
    }
#endif
    // This line shouldn't be reached
    outError("Unsupported " + Params::getInstance().dating_method + " dating method");
}
