/*
 * timetree.cpp
 * Interface to call dating method incl. LSD2
 *  Created on: Apr 4, 2020
 *      Author: minh
 */

#include "timetree.h"
#include <utils/io.h> //for safeGetLine

#ifdef USE_LSD2
#include "lsd2/src/lsd.h"
#endif

/** map from taxon name to date */
typedef unordered_map<string, string> TaxonDateMap;
#define YEAR_SCALE 100000

/**
 @param[in] date date string
 @return converted date as a float or YYYY-MM[-DD] format
 */
string convertDate(string date) {
    // check for range in x:y format
    if (date.find(':') != string::npos) {
        StrVector vec;
        convert_string_vec(date.c_str(), vec, ':');
        if (vec.size() != 2)
            outError("Invalid date range " + date);
        if (vec[0].empty() || vec[0] == "NA")
            return "u(" + vec[1] + ")";
        if (vec[1].empty() || vec[1] == "NA")
            return "l(" + vec[0] + ")";

        return "b(" + vec[0] + "," + vec[1] + ")";
    }
    if (date.empty() || !isdigit(date[0]) || date[0] == '-')
        return date;
    DoubleVector vec;
    try {
        convert_double_vec(date.c_str(), vec, '-');
    } catch (...) {
        outError("Invalid date " + date);
    }
    // otherwise, return the original date string
    return date;
}

/**
 read a date file. Each line has two strings: name and date
 */
void readDateFile(string date_file, set<string> &node_names, TaxonDateMap &dates) {
    try {
        cout << "Reading date file " << date_file << " ..." << endl;
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(date_file);
        int line_num;
        for (line_num = 1; !in.eof(); line_num++) {
            string line_out = "Line " + convertIntToString(line_num) + ": ";
            string line;
            if (!safeGetLine(in, line))
                break;
            // ignore comment
            if (line.find('#') != string::npos)
                line = line.substr(0, line.find('#'));
            trimString(line);
            if (line.empty()) // ignore empty line
                continue;
            string name, date;
            istringstream line_in(line);
            if (!(line_in >> name >> date))
                throw line_out + "'" + line + "' does not contain name and date";
            // error checking, make sure that name appear in tree
            StrVector name_vec;
            convert_string_vec(name.c_str(), name_vec);
            for (auto s : name_vec)
                if (node_names.find(s) == node_names.end())
                    throw line_out + "'" + s + "' does not appear in tree";
            // error checking, make sure is date is valid
            if (date.empty())
                throw line_out + "date is empty";
            if (date.substr(0,2) != "NA")
            try {
                int end_pos;
                convert_double(date.c_str(), end_pos);
            } catch (string str) {
                throw line_out + str;
            }
            dates[name] = date;
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

/** read the date information from the alignment taxon names */
void readDateTaxName(set<string> &nodenames, TaxonDateMap &dates) {
    cout << "Extracting date from node names..." << endl;
    for (string name : nodenames) {
        // get the date in the taxon name after the '|' sign
        auto pos = name.rfind('|');
        if (pos == string::npos)
            continue;
        string date = name.substr(pos+1);
        try {
            // try to parse
            int end_pos;
            convert_double(date.c_str(), end_pos);
            // it works! so get the date
            dates[name] = date;
        } catch (...) {
            // does not work, ignore the taxon name
            continue;
        }
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

void writeDate(string date_file, ostream &out, set<string> &nodenames) {
    TaxonDateMap dates;
    if (date_file == "TAXNAME") {
        // read the dates from alignment taxon names
        readDateTaxName(nodenames, dates);
    } else {
        readDateFile(date_file, nodenames, dates);
    }
    // only retain taxon appearing in alignment
    TaxonDateMap retained_dates;
    set<string> outgroup_set;
    if (Params::getInstance().root) {
        StrVector outgroup_names;
        convert_string_vec(Params::getInstance().root, outgroup_names);
        for (auto name : outgroup_names)
            outgroup_set.insert(name);
    }
    if (verbose_mode >= VB_MED)
        cout << "Node\tDate" << endl;
    for (auto name: nodenames) {
        string date = "NA";
        if (dates.find(name) == dates.end()) {
            // taxon present in the dates
//            if (!Params::getInstance().date_tip.empty())
//                date = Params::getInstance().date_tip;
        } else if (outgroup_set.find(name) == outgroup_set.end() || Params::getInstance().date_with_outgroup) {
            // ignore the date of the outgroup
            date = dates[name];
        }
        if (date != "NA") {
            retained_dates[name] = date;
            dates.erase(name);
        }
        if (verbose_mode >= VB_MED)
            cout << name << "\t" << date << endl;
    }
    
    // add remaining ancestral dates
    for (auto date : dates) {
        if (date.first.substr(0,4) == "mrca" || date.first.substr(0,8) == "ancestor")
            retained_dates[date.first] = date.second;
        else if (date.first.find(',') != string::npos) {
            retained_dates["ancestor(" + date.first + ")"] = date.second;
        } else if (outgroup_set.find(date.first) == outgroup_set.end() || Params::getInstance().date_with_outgroup) {
            retained_dates[date.first] = date.second;
        }
    }

//    if (!Params::getInstance().date_root.empty()) {
//        retained_dates["root"] = Params::getInstance().date_root;
//    }
    
    cout << retained_dates.size() << " dates extracted" << endl;
    try {
        out << retained_dates.size() << endl;
        for (auto date : retained_dates) {
            out << date.first << " " << convertDate(date.second) << endl;
        }
    } catch (...) {
        ASSERT(0 && "Error writing date stream");
    }
}

#ifdef USE_LSD2
void runLSD2(PhyloTree *tree) {
    string basename = (string)Params::getInstance().out_prefix + ".timetree";
    string treefile = basename + ".subst";
    stringstream tree_stream, outgroup_stream, date_stream;
    tree->printTree(tree_stream);
    StrVector arg = {"lsd", "-i", treefile, "-s", convertIntToString(tree->getAlnNSite()), "-o", basename};
    if (Params::getInstance().date_debug) {
        ofstream out(treefile);
        out << tree_stream.str();
        out.close();
        cout << "Tree printed to " << treefile << endl;
    }
    
    if (Params::getInstance().date_replicates > 0) {
        arg.push_back("-f");
        arg.push_back(convertIntToString(Params::getInstance().date_replicates));
        if (Params::getInstance().clock_stddev >= 0) {
            arg.push_back("-q");
            arg.push_back(convertDoubleToString(Params::getInstance().clock_stddev));
        }
    }

    if (Params::getInstance().date_outlier >= 0) {
        arg.push_back("-e");
        arg.push_back(convertIntToString(Params::getInstance().date_outlier));
    }
    
    if (Params::getInstance().root) {
        // print outgroup file for LSD
        writeOutgroup(outgroup_stream, Params::getInstance().root);
        string outgroup_file = basename + ".outgroup";
        arg.push_back("-g");
        arg.push_back(outgroup_file); // only fake file
        if (!Params::getInstance().date_with_outgroup)
            arg.push_back("-G");
        if (Params::getInstance().date_debug) {
            ofstream out(outgroup_file);
            out << outgroup_stream.str();
            out.close();
            cout << "Outgroup printed to " << outgroup_file << endl;
        }
    } else {
        // search for all possible rootings
        arg.push_back("-r");
        arg.push_back("a");
    }

    if (Params::getInstance().date_file != "") {
        // parse the date file
        set<string> nodenames;
        tree->getNodeName(nodenames);
        writeDate(Params::getInstance().date_file, date_stream, nodenames);
        string date_file = basename + ".date";
        arg.push_back("-d");
        arg.push_back(date_file);  // only fake file
        if (Params::getInstance().date_debug) {
            ofstream out(date_file);
            out << date_stream.str();
            out.close();
            cout << "Date file printed to " << date_file << endl;
        }
    }
    // input tip and root date
    if (Params::getInstance().date_root != "") {
        arg.push_back("-a");
        arg.push_back(convertDate(Params::getInstance().date_root));
    }
    
    if (Params::getInstance().date_tip != "") {
        arg.push_back("-z");
        arg.push_back(convertDate(Params::getInstance().date_tip));
    }

    lsd::InputOutputStream io(tree_stream.str(), outgroup_stream.str(), date_stream.str(), "", "", "");

    if (Params::getInstance().dating_options != "") {
        // extra options for LSD
        StrVector options;
        convert_string_vec(Params::getInstance().dating_options.c_str(), options, ' ');
        for (auto opt : options)
            if (!opt.empty())
                arg.push_back(opt);
    }
    
    cout << "Building time tree by least-square dating (LSD) with command:" << endl;
    
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
    
    if (((stringstream*)io.outTree3)->str().empty()) {
        outError("Something went wrong, LSD could not date the tree");
    } else {
        cout << "LSD results written to:" << endl;
        cout << "  LSD report:                  " << report_file << endl;
    //    cout << "  Time tree in nexus format:   " << tree1_file << endl;
        cout << "  Time tree in nexus format:   " << tree2_file << endl;
        cout << "  Time tree in newick format:  " << tree3_file << endl;
        cout << endl;
    }
}
#endif

void doTimeTree(PhyloTree *tree) {

    cout << "--- Start phylogenetic dating ---" << endl;
    cout.unsetf(ios::fixed);

#ifdef USE_LSD2
    if (Params::getInstance().dating_method == "LSD") {
        runLSD2(tree);
        cout << "--- End phylogenetic dating ---" << endl;
        return;
    }
#endif
    // This line shouldn't be reached
    outError("Unsupported " + Params::getInstance().dating_method + " dating method");
}
