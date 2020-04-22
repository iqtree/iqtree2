/*
 * timetree.cpp
 * Interface to call dating method incl. LSD2
 *  Created on: Apr 4, 2020
 *      Author: minh
 */

#include "timetree.h"
#include "lsd2/src/lsd.h"

/** map from taxon name to date */
typedef unordered_map<string, string> TaxonDateMap;
#define YEAR_SCALE 100000

/**
 @param[in] date date string
 @param is_ISO true to force the date being in YYYY-MM-DD format
 @return converted date as a float or YYYY-MM[-DD] format
 */
string convertDate(string date, bool is_ISO) {
    if (date.empty() || !isdigit(date[0]))
        return date;
    DoubleVector vec;
    convert_double_vec(date.c_str(), vec, '-');
    IntVector month_days = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (vec.size() == 1 && fabs(vec[0] - floor(vec[0])) < 1e-6 && is_ISO) {
        // incomplete YYYY format, convert it to range YYYY-01-01 to YYYY-12-31
        return "b(" + date + "-01-01," + date + "-12-31)";
    }
    if (vec.size() == 2) {
        // incomplete YYYY-MM date string, convert it to range from 1st to last day of month
        if (vec[1] < 1 || vec[1] > month_days.size())
            outError("Invalid month in ", date);
        return "b(" + date + "-01," + date + "-" + convertIntToString(month_days[vec[1]-1]) + ")";
    }
    // otherwise, return the original date string
    return date;
}

/**
 check if any date is in YYYY-MM-DD format
 */
bool hasISODate(TaxonDateMap &dates) {
    for (auto date : dates) {
        DoubleVector vec;
        convert_double_vec(date.second.c_str(), vec, '-');
        if (vec.size() > 1)
            return true;
    }
    return false;
}
/**
 read a date file. Each line has two strings: name and date
 */
void readDateFile(string date_file, TaxonDateMap &dates) {
    try {
        cout << "Reading date file " << date_file << " ..." << endl;
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
void readDateAlignment(Alignment *aln, TaxonDateMap &dates) {
    cout << "Extracting date from taxa names..." << endl;
    for (int i = 0; i < aln->getNSeq(); i++) {
        string name = aln->getSeqName(i);
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

void writeDate(string date_file, ostream &out, Alignment *aln) {
    TaxonDateMap dates;
    if (date_file == "TAXNAME") {
        // read the dates from alignment taxon names
        readDateAlignment(aln, dates);
    } else {
        readDateFile(date_file, dates);
    }
    bool is_ISO = hasISODate(dates);
    // only retain taxon appearing in alignment
    TaxonDateMap retained_dates;
    cout << "ID\tTaxon\tDate" << endl;
    for (int i = 0; i < aln->getNSeq(); i++) {
        string name = aln->getSeqName(i);
        string date;
        if (dates.find(name) == dates.end())
            date = "NA";
        else {
            date = dates[name];
            retained_dates[name] = date;
        }
        cout << i+1 << "\t" << name << "\t" << date << endl;
    }
    cout << retained_dates.size() << " dates extracted" << endl;
    try {
        out << retained_dates.size() << endl;
        for (auto date : retained_dates) {
            out << date.first << " " << convertDate(date.second, is_ISO) << endl;
        }
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
    if (Params::getInstance().date_debug) {
        ofstream out(treefile);
        out << tree_stream.str();
        out.close();
        cout << "Tree printed to " << treefile << endl;
    }
    
    if (Params::getInstance().root) {
        // print outgroup file for LSD
        writeOutgroup(outgroup_stream, Params::getInstance().root);
        string outgroup_file = basename + ".outgroup";
        arg.push_back("-g");
        arg.push_back(outgroup_file); // only fake file
        if (Params::getInstance().date_with_outgroup)
            arg.push_back("-k");
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
        writeDate(Params::getInstance().date_file, date_stream, tree->aln);
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
    
    if (Params::getInstance().date_root != "") {
        arg.push_back("-a");
        arg.push_back(Params::getInstance().date_root);
    }
    
    if (Params::getInstance().date_tip != "") {
        arg.push_back("-z");
        arg.push_back(Params::getInstance().date_tip);
    }
    
    lsd::InputOutputStream io(tree_stream.str(), outgroup_stream.str(), date_stream.str(), "");

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
