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
    cout << "Building time tree by least-square dating (LSD) with command:" << endl;
    string basename = (string)Params::getInstance().out_prefix + ".timetree";
    string treefile = basename + ".subst";
    tree->printTree(treefile.c_str());
    StrVector arg = {"lsd", "-i", treefile, "-s", convertIntToString(tree->getAlnNSite()), "-c", "-o", basename};
    if (Params::getInstance().root) {
        // print outgroup file for LSD
        StrVector outgroup_names;
        convert_string_vec(Params::getInstance().root, outgroup_names);
        string outgroup_file = basename + ".outgroup";
        ofstream out;
        out.open(outgroup_file);
        out << outgroup_names.size() << endl;
        for (auto outgroup : outgroup_names) {
            out << outgroup << endl;
        }
        out.close();
        arg.push_back("-g");
        arg.push_back(outgroup_file);
        arg.push_back("-k");
    } else {
        // search for all possible rootings
        arg.push_back("-r");
        arg.push_back("a");
    }

    if (Params::getInstance().dating_options != "") {
        // extra options for LSD
        StrVector options;
        convert_string_vec(Params::getInstance().dating_options.c_str(), options, ' ');
        for (auto opt : options)
            if (!opt.empty())
                arg.push_back(opt);
    }
    
    int argc = arg.size();
    char *argv[argc];
    for (int i = 0; i < argc; i++)
        argv[i] = (char*)arg[i].c_str();
    std::copy(arg.begin(), arg.end(), std::ostream_iterator<string>(std::cout, " "));
    cout << endl;
    lsd_main(argc, argv);
    cout << "LSD results written to:" << endl;
    cout << "  LSD report:                  " << basename << endl;
    cout << "  Time tree in nexus format:   " << basename << ".nexus" << endl;
    cout << "  Time tree with dates:        " << basename << ".date.nexus" << endl;
    cout << "  Time tree in newick format:  " << basename << ".nwk" << endl;
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
