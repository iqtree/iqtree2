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


/**
 convert a date string in YYYY-MM-DD format into integer in unit of year*million
 */
string convertDateForLSD(string date_str) {
    
    StrVector str_vec;
    convert_string_vec(date_str.c_str(), str_vec, '-');
    if (str_vec.size() > 3)
        outError("Invalid date string " + date_str + ": not in YYYY-MM-DD format");

    int64_t year = convert_int64(str_vec[0].c_str());
    int64_t month = 1;
    int64_t day = 1;
    if (str_vec.size() >= 2)
        month = convert_int64(str_vec[1].c_str());
    if (str_vec.size() == 3) {
        day = convert_int64(str_vec[2].c_str());
    }
    if (month < 1 || month > 12)
        outError("Invalid date string " + date_str + ": month is not between 1 and 12");
    int64_t days_per_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int64_t days_of_year = 365;
    if (year % 4 == 0) {
        // leap year
        days_per_month[1]++;
        days_of_year++;
    }
    if (day < 1 || day > days_per_month[month-1])
        outError("Invalid date string " + date_str + ": day is not between 1 and " +
                 convertIntToString(days_per_month[month-1]));

    int64_t elapsed_days = day-1;
    for (int i = 0; i < month-1; i++)
        elapsed_days += days_per_month[i];
    return convertInt64ToString(year*YEAR_SCALE + elapsed_days*YEAR_SCALE/days_of_year);
}


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
        try {
            ofstream out;
            out.open(outgroup_file);
            out << outgroup_names.size() << endl;
            for (auto outgroup : outgroup_names) {
                out << outgroup << endl;
            }
            out.close();
        } catch (...) {
            outError("Error writing file " + outgroup_file);
        }
        arg.push_back("-g");
        arg.push_back(outgroup_file);
        arg.push_back("-k");
    } else {
        // search for all possible rootings
        arg.push_back("-r");
        arg.push_back("a");
    }

    if (Params::getInstance().date_file != "") {
        // parse the date file
        DateVector dates;
        cout << "Reading date file " << Params::getInstance().date_file << " ..." << endl;
        readDateFile(Params::getInstance().date_file, dates);
        // only retain taxon appearing in alignment
        DateVector converted_dates;
        cout << "ID\tTaxon\tDate\tConverted" << endl;
        for (auto date: dates) {
            int taxon_id = tree->aln->getSeqID(date.first);
            if (taxon_id >= 0) {
                string converted_date = convertDateForLSD(date.second);
                converted_dates.push_back(std::make_pair(date.first, converted_date));
                cout << taxon_id+1 << "\t" << date.first << "\t" << date.second << "\t" << converted_date << endl;
            }
        }
        cout << converted_dates.size() << " dates extracted" << endl;
        string date_file = basename + ".date";
        try {
            ofstream out;
            out.open(date_file);
            out << converted_dates.size() << endl;
            for (auto date : converted_dates)
                out << date.first << " " << date.second << endl;
            out.close();
        } catch (...) {
            outError("Error writing file " + date_file);
        }
        arg.push_back("-d");
        arg.push_back(date_file);
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
