/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*
	Geneset selection (GSS) for Roland
*/

#include "gss.h"
#include "pda/lpwrapper.h"
#include "pda/gurobiwrapper.h"
#include "tree/mtreeset.h"


GSSNetwork::GSSNetwork(Params &params) : PDNetwork(params) {
    readGenePValues(params);
}

bool GSSNetwork::isPDArea() {
    return false;
}

void GSSNetwork::readGenePValues(Params &params) {

    //taxa->Report(cout);
    // first build the gene list
    TaxaSetNameVector *allsets = sets->getSets();
    TaxaSetNameVector::iterator i;
    for (i = allsets->begin(); i != allsets->end(); i++) {
        for (vector<string>::iterator it2 = (*i)->taxlist.begin(); it2 != (*i)->taxlist.end(); it2++) {
            if (gene_index.find(*it2) == gene_index.end()) {
                gene_index[*it2] = genes.size();
                genes.push_back(*it2);
            }
        }
    }
    int ntaxa = genes.size();

    // build the area_taxa structure
    if (allsets->size() != getNTaxa())
        outError("Number of gene sets do not match between tree file and set file");
    area_taxa.resize(getNTaxa(), NULL);
    for (i = allsets->begin(); i != allsets->end(); i++) {
        int id = -1;
        try {
            id = taxa->FindTaxon(NxsString((*i)->name.c_str()));
        } catch (NxsTaxaBlock::NxsX_NoSuchTaxon) {
            outError(ERR_NO_TAXON, (*i)->name);
        }
        if (area_taxa[id]) outError("Duplicated set name in set file", (*i)->name);
        Split *sp = new Split(ntaxa);
        for (vector<string>::iterator it2 = (*i)->taxlist.begin(); it2 != (*i)->taxlist.end(); it2++) {
            sp->addTaxon(gene_index[*it2]);
        }
        area_taxa[id] = sp;
        cout << id << "\t" << (*i)->name << endl;
    }
    cout << ntaxa << " genes and " << area_taxa.size() << " gene sets detected" << endl;

    cout << "Reading p-values file " << params.gene_pvalue_file << " ..." << endl;
    gene_pvalues.resize(ntaxa, -1);
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(params.gene_pvalue_file);
        string name, tmp;

        for (; !in.eof() && ntaxa > 0; ntaxa--) {
            // remove the failbit
            in.exceptions(ios::badbit);
            if (!(in >> name)) break;
            // set the failbit again
            in.exceptions(ios::failbit | ios::badbit);
            if (gene_index.find(name) == gene_index.end())
                outError("A gene not found in gene p-values file");
            // read the sequence weight
            in >> tmp;
            double pval = convert_double(tmp.c_str());
            if (pval < 0 || pval > 1) outError("Some pvalue is out of range [0, 1]");
            if (gene_pvalues[gene_index[name]] != -1) outError("Duplicated p-value entry");
            gene_pvalues[gene_index[name]] = pval;
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (string str) {
        outError(str);
    }

    if (params.gene_scale_factor < 0 || params.gene_scale_factor > 1)
        outError("gene_scale_factor must be in range [0,1]");
    cout << "Rescaling split weights with " << params.gene_scale_factor <<
         " and gene p-values with " << 1 - params.gene_scale_factor << endl;
    // incoporate into the split system
    for (iterator it = begin(); it != end(); it++) {
        // first, multiply split weight with the coefficient
        (*it)->setWeight((*it)->getWeight() * params.gene_scale_factor);
    }

    for (DoubleVector::iterator it2 = gene_pvalues.begin(); it2 != gene_pvalues.end(); it2++)
        if (params.gene_pvalue_loga)
            (*it2) = (-log(*it2)) * (1 - params.gene_scale_factor);
        else
            (*it2) = (1 - (*it2)) * (1 - params.gene_scale_factor);

}

void GSSNetwork::checkZValue(int total_size, vector<int> &z_value) {
    z_value.resize(genes.size(), -1);
    int i, j;
    for (i = 0; i < genes.size(); i++) {
        int genesetid = -1;
        for (j = 0; j < area_taxa.size(); j++)
            if (area_taxa[j]->containTaxon(i))  {
                if (genesetid < 0)
                    genesetid = j;
                else {
                    genesetid = -1;
                    break;
                }
            }
        if (genesetid >= 0) z_value[i] = genesetid+2;
    }
}


void GSSNetwork::lpObjectiveGSS(ostream &out, Params &params, IntVector &y_value, IntVector &z_value, int total_size) {
    //IntVector y_value, count1, count2;
    iterator spit;
    int i;
    // define the objective function
    if (params.gurobi_format)
        out << "Maximize" << endl;
    else
        out << "max: ";

    // first compute the coefficient for x variable
    DoubleVector xweights;
    xweights.resize(getNTaxa(), 0.0);
    for (spit = begin(),i=0; spit != end(); spit++,i++)	{
        if (y_value[i] >= 2)
            xweights[y_value[i] - 2] += (*spit)->getWeight();
    }
    for (i = 0; i < gene_pvalues.size(); i++)
        if (z_value[i] >= 2)
            xweights[z_value[i]-2] += gene_pvalues[i];

    // now write down the objective function

    for (i = 0; i < xweights.size(); i++)
        out << " +" << xweights[i] << " x" << i;


    for (spit = begin(),i=0; spit != end(); spit++,i++)	{
        if (y_value[i] < 0)
            out << " +" << (*spit)->getWeight() << " y" << i;
    }

    for (i = 0; i < gene_pvalues.size(); i++)
        if (z_value[i] < 0)
            out << " +" << gene_pvalues[i] << " z" << i;

    if (params.gurobi_format)
        out << endl << "Subject to" << endl;
    else
        out << ";" << endl;
}


void GSSNetwork::lpVariableBound(ostream &out, Params &params, Split &included_vars, IntVector &y_value, IntVector &z_value) {
    int i;
    PDNetwork::lpVariableBound(out, params, included_vars, y_value);

    for (i = 0; i < gene_pvalues.size(); i++) {
        if (z_value[i] >= 0) continue;
        if (params.gurobi_format)
            out << "0 <= ";
        out << "z" << i << " <= 1";
        if (params.gurobi_format)
            out << endl;
        else
            out << ";" << endl;
    }
}

void GSSNetwork::lpGeneConstraint(ostream &out, Params &params, IntVector &z_value) {
    int i, j;
    for (i = 0; i < genes.size(); i++) {
        if (z_value[i] >= 0) continue;
        out << "z" << i;
        for (j = 0; j < area_taxa.size(); j++)
            if (area_taxa[j]->containTaxon(i))
                out << " -x" << j;
        out << " <= 0";
        if (params.gurobi_format)
            out << endl;
        else
            out << ";" << endl;
    }
}

void GSSNetwork::transformLP_GSS(Params &params, const char *outfile, int total_size, bool make_bin) {
    Split included_tax(getNTaxa());
    IntVector::iterator it2;
    for (it2 = initialset.begin(); it2 != initialset.end(); it2++)
        included_tax.addTaxon(*it2);
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(outfile);
        vector<int> y_value;
        vector<int> z_value;
        checkYValue(total_size, y_value);
        checkZValue(total_size, z_value);

        lpObjectiveGSS(out, params, y_value, z_value, total_size);
        lpSplitConstraint_TS(out, params, y_value, total_size);
        lpK_BudgetConstraint(out, params, total_size);
        lpGeneConstraint(out, params, z_value);
        lpVariableBound(out, params, included_tax, y_value, z_value);
        if (make_bin)
            lpVariableBinary(out, params, included_tax);

        out.close();
        //cout << "Transformed LP problem printed to " << outfile << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, outfile);
    }
}

void GSSNetwork::findPD(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order) {
    // call the entering function
    if (isBudgetConstraint()) { // non-budget case
        cout << "Please specify k";
        return;
    }
    enterFindPD(params);
    if (params.find_all)
        outError("Current linear programming does not support multiple optimal sets!");

    string ofile = params.out_prefix;
    ofile += ".lp";
    double score;
    int lp_ret, i, ntaxa = getNTaxa();
    int k, min_k, max_k, step_k, index;

    double *variables = new double[ntaxa];

    if (isBudgetConstraint()) { // non-budget case
        min_k = params.min_budget;
        max_k = params.budget;
        step_k = params.step_budget;
    } else {
        min_k = params.min_size;
        max_k = params.sub_size;
        step_k = params.step_size;
    }
    taxa_set.resize((max_k - min_k)/step_k + 1);

    // now construction the optimal PD sets
    if (isBudgetConstraint())
        cout << "running budget = ";
    else
        cout << "running k = ";
    for (k = min_k; k <= max_k; k += step_k) {
        index = (k - min_k) / step_k;
        if (!params.binary_programming) {
            transformLP_GSS(params, ofile.c_str(), k, false);
            cout << " " << k;
            cout.flush();
            if (params.gurobi_format)
                lp_ret = gurobi_solve((char*)ofile.c_str(), ntaxa, &score, variables, verbose_mode, params.gurobi_threads);
            else
                lp_ret = lp_solve((char*)ofile.c_str(), ntaxa, &score, variables, verbose_mode);
        } else lp_ret = 7;
        if (lp_ret != 0 && lp_ret != 7)
            outError("Something went wrong with LP solver!");
        if (lp_ret == 7) { // fail with non-binary case, do again with strict binary
            if (params.binary_programming)
                transformLP_GSS(params, ofile.c_str(), k, true);
            else
                lpVariableBinary(ofile.c_str(), params, initialset);
            cout << " " << k << "(bin)";
            cout.flush();
            if (params.gurobi_format)
                lp_ret = gurobi_solve((char*)ofile.c_str(), ntaxa, &score, variables, verbose_mode, params.gurobi_threads);
            else
                lp_ret = lp_solve((char*)ofile.c_str(), ntaxa, &score, variables, verbose_mode);
            if (lp_ret != 0) // check error again without allowing non-binary
                outError("Something went wrong with LP solver!");
        }

        Split *pd_set = new Split(ntaxa, score);
        for (i = 0; i < ntaxa; i++)
            if (1.0 - variables[i] < tolerance) {
                //pd_set->addTaxon(taxa_order[i]);
                pd_set->addTaxon(i);
            }
        calcPD(*pd_set);
        taxa_set[index].push_back(pd_set);
    }
    cout << endl;
    delete [] variables;
    // call the leaving function
    leaveFindPD(taxa_set);
}

extern void summarizeSplit(Params &params, PDNetwork &sg, vector<SplitSet> &pd_set, PDRelatedMeasures &pd_more, bool full_report);

void runGSSAnalysis(Params &params) {
    cout << "Dedicated for Roland..." << endl;
    vector<SplitSet> taxa_set;
    IntVector taxa_order;
    StrVector genes;
    DoubleVector gene_pvalues;
    PDRelatedMeasures pd_more;

    params.intype = detectInputFile(params.user_file);

    GSSNetwork sg(params);

    sg.findPD(params, taxa_set, taxa_order);

    summarizeSplit(params, sg, taxa_set, pd_more, true);
}

