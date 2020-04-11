/*
 * timetree.cpp
 * Interface to call dating method incl. LSD2
 *  Created on: Apr 4, 2020
 *      Author: minh
 */

#include "timetree.h"
#include "lsd2/src/lsd.h"

typedef vector<pair<string, string> > DateVector;
#define YEAR_SCALE 100000

/**
 read a date file. Each line has two strings: name and date
 */
void readDateFile(string date_file, DateVector &dates) {
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(date_file);
        int line_num;
        for (line_num = 1; !in.eof(); line_num++) {
            string line;
            if (!safeGetline(in, line))
                break;
            if (line.empty()) // ignore empty line
                continue;
            string name, date;
            istringstream line_in(line);
            if (!(line_in >> name >> date))
                throw "Line " + convertIntToString(line_num) + ": '" + line + "' does not contain name and date";
            dates.push_back(std::make_pair(name, date));
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, date_file);
    } catch (string str) {
        outError(str);
    } catch (...) {
        outError("Error reading date file " + date_file);
    }
}


void writeOutgroup(ostream &out, const char *outgroup) {
    StrVector outgroup_names;
    convert_string_vec(outgroup, outgroup_names);
    try {
        out << outgroup_names.size() << endl;
        for (auto outgroup : outgroup_names) {
            out << outgroup << endl;
        }
    } catch (...) {
        ASSERT(0 && "Error writing outgroup stream");
    }
}

void writeDate(ostream &out, PhyloTree *tree) {
    DateVector dates;
    cout << "Reading date file " << Params::getInstance().date_file << " ..." << endl;
    readDateFile(Params::getInstance().date_file, dates);
    // only retain taxon appearing in alignment
    DateVector retained_dates;
    cout << "ID\tTaxon\tDate" << endl;
    for (auto date: dates) {
        int taxon_id = tree->aln->getSeqID(date.first);
        if (taxon_id >= 0) {
            retained_dates.push_back(date);
            cout << taxon_id+1 << "\t" << date.first << "\t" << date.second << endl;
        }
    }
    cout << retained_dates.size() << " dates extracted" << endl;
    try {
        out << retained_dates.size() << endl;
        for (auto date : retained_dates)
            out << date.first << " " << date.second << endl;
    } catch (...) {
        ASSERT(0 && "Error writing date stream");
    }
}

#ifdef USE_LSD2
void runLSD2(PhyloTree *tree) {
    cout << "Building time tree by least-square dating (LSD) with command:" << endl;
    
    string basename = (string)Params::getInstance().out_prefix + ".timetree";
    string treefile = basename + ".subst";
    stringstream tree_stream, outgroup_stream, date_stream;
    tree->printTree(tree_stream);
    StrVector arg = {"lsd", "-i", treefile, "-s", convertIntToString(tree->getAlnNSite()), "-c", "-o", basename};
    
    if (Params::getInstance().root) {
        // print outgroup file for LSD
        writeOutgroup(outgroup_stream, Params::getInstance().root);
        string outgroup_file = basename + ".outgroup";
        arg.push_back("-g");
        arg.push_back(outgroup_file); // only fake file
        arg.push_back("-k");
    } else {
        // search for all possible rootings
        arg.push_back("-r");
        arg.push_back("a");
    }

    if (Params::getInstance().date_file != "") {
        // parse the date file
        writeDate(date_stream, tree);
        string date_file = basename + ".date";
        arg.push_back("-d");
        arg.push_back(date_file);  // only fake file
    }
    
    lsd::InputOutputStream io(tree_stream.str(), outgroup_stream.str(), date_stream.str());

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
    
    // main call to LSD!
    lsd::buildTimeTree(argc, argv, &io);
    
    // fetch the output
    string report_file = basename + ".lsd";
    //string tree1_file = basename + ".raw";
    string tree2_file = basename + ".nex";
    string tree3_file = basename + ".nwk";
    try {
        ofstream out;
        out.open(report_file);
        out << ((ostringstream*)io.outResult)->str();
        out.close();
//        out.open(tree1_file);
//        out << ((ostringstream*)io.outTree1)->str();
//        out.close();
        out.open(tree2_file);
        out << ((ostringstream*)io.outTree2)->str();
        out.close();
        out.open(tree3_file);
        out << ((stringstream*)io.outTree3)->str();
        out.close();
    } catch (...) {
        outError("Couldn't write LSD output files");
    }
    
    cout << "LSD results written to:" << endl;
    cout << "  LSD report:                  " << report_file << endl;
//    cout << "  Time tree in nexus format:   " << tree1_file << endl;
    cout << "  Time tree in nexus format:   " << tree2_file << endl;
    cout << "  Time tree in newick format:  " << tree3_file << endl;
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
