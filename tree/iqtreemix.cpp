//
//  iqtreemix.cpp
//  tree
//
//  Created by Thomas Wong on 14/12/20.
//

#include "iqtreemix.h"
const double MIN_PROP = 0.001;
const double MAX_PROP = 1000.0;
// const double MIN_LEN = 1e-3;
const double MAX_LEN = 1.0;
const double LIKE_THRES = 0.1; // 10% more in team of likelihood value
const double WEIGHT_EPSILON = 0.001;
const int OPTIMIZE_STEPS = 10000;

// Input formats for the tree-mixture model
// 1. linked models and site rates: GTR+G4+T
// 2. unlinked models and linked site rates: MIX{GTR,GTR}+G4+T
// 3. linked models and unlinked site rates: GTR+MIX{G4,E}+T
// 4. unlinked models and unlinked site rates: MIX{GTR+G4,GTR}+T
// The situation that a part of the model is linked while another part is unlinked is not allowed.
//    For example, MIX{GTR,GTR}+FO+T or GTR+MIX{FO,FO}+T is not be accepted
// Similarly, the situation that a part of the site rate is linked while another part is unlinked is also not allowed.
//    For example, GTR+MIX{I,I}+G4+T or GTR+I+MIX{G4+G4}+T is not be accepted
// The tree weights can be fixed by: T{0.2,0.3,0.5}, otherwise, optimization on tree weights will be performed.

bool isRHS(string m) {
    if (m.length() > 0 && (m[0]=='I' || m[0]=='G' || m[0]=='R')) {
        if (m.length() > 1) {
            return (isdigit(m[1]) || m[1]=='{');
        }
        return true;
    }
    return false;
}

// return how many char c inside the infile
int checkCharInFile(char* infile, char c) {
    ifstream fin;
    string aline;
    size_t i,k;
    k=0;
    fin.open(infile);
    while (getline(fin,aline)) {
        for (i=0; i<aline.length(); i++) {
            if (aline[i] == c)
                k++;
        }
    }
    fin.close();
    return k;
}

// get the number of trees for the tree-mixture model
int getTreeMixNum(Params& params) {
    int n = 0;
    size_t p = params.model_name.find("+T");
    string str_n;
    if (p != string::npos && p < params.model_name.length()-2) {
        str_n = params.model_name.substr(p+2);
        n = atoi(str_n.c_str());
    }
    if (n == 1) {
        outError("The number after +T has to be greater than 1");
    }
    // check how many trees inside the user input file
    int k = checkCharInFile(params.user_file, ';');
    if (k <= 1) {
        outError("Tree mixture model only supports at least 2 trees inside the tree file: " + string(params.user_file) + ". Each tree must be followed by the character ';'.");
    }
    if (n == 0) {
        n = k;
        cout << "Number of input trees: " << n << endl;
    } else if (n < k) {
        cout << "Note: Only " << n << " trees are considered, although there are more than " << n << " trees in the tree file: " << params.user_file << endl;
    } else if (n > k) {
        outError("The number of trees inside the tree file '" + string(params.user_file) + "' is less than " + convertIntToString(n));
    }
    return n;
}


IQTreeMix::IQTreeMix() : IQTree() {
    ptn_freq = NULL; // double frequencies of each pattern (can be changed)
    patn_isconst = NULL;
    patn_parsimony = NULL;
    ptn_like_cat = NULL;
    _ptn_like_cat = NULL;
    ptn_scale_cat = NULL;
    single_ptn_tree_like = NULL;
    ptn_like = NULL;
    _pattern_scaling = NULL;
    // initialize the variables
    isTreeWeightFixed = false;
    weightGrpExist = false;
    isEdgeLenRestrict = false;
    optim_type = 1;
    parsi_computed = false;
    logl_variance = 0.0;
    isLinkModel = true;
    isLinkSiteRate = true;
    anySiteRate = false;
    isNestedOpenmp = false;
    rhas_var = NULL;
    ntree = 0;
}

IQTreeMix::IQTreeMix(Params &params, Alignment *aln) : IQTree(aln) {
    size_t i;

    clear();
    weights.clear();
    weight_logs.clear();
    
    // checking on params
    if (params.compute_ml_tree_only) {
        outError("option compute_ml_tree_only cannot be set for tree-mixture model");
    }
    // change to 0.04 for tree mixture model as 0.02 and 0.03 cause numerical problems
    if (params.min_gamma_shape < MIN_GAMMA_SHAPE_TREEMIX) {
        if (params.min_gamma_shape != MIN_GAMMA_SHAPE)
            cout << "The minimum value for Gamma shape is changed to " << MIN_GAMMA_SHAPE_TREEMIX << endl;
        params.min_gamma_shape = MIN_GAMMA_SHAPE_TREEMIX;
    }
    if (params.user_file == NULL) {
        outError("To use tree-mixture model, use an option: -te <newick file with multiple trees>");
    }

    // get the number of trees for tree-mixture model
    ntree = getTreeMixNum(params);
    
    // create the trees and initialize tree-weights
    double init_weight = 1.0 / (double) ntree;
    double init_weight_log = log(init_weight);
    for (i=0; i<ntree; i++) {
        push_back(new IQTree(aln));
        weights.push_back(init_weight);
        weight_logs.push_back(init_weight_log);
    }
    
    // allocate memory for the arrays
    nptn = aln->getNPattern();
    size_t mem_size = get_safe_upper_limit(nptn);
    size_t mem_size32 = get_safe_upper_limit_float(nptn);
    size_t block_size = get_safe_upper_limit(nptn) * ntree;
    size_t block_size32 = get_safe_upper_limit_float(nptn) * ntree;

    ptn_freq = aligned_alloc<double>(mem_size);
    patn_isconst = aligned_alloc<int>(mem_size32);
    ptn_like_cat = aligned_alloc<double>(block_size);
    _ptn_like_cat = aligned_alloc<double>(block_size);
    ptn_scale_cat = aligned_alloc<double>(block_size);
    patn_parsimony = aligned_alloc<int>(block_size32);
    single_ptn_tree_like = aligned_alloc<double>(get_safe_upper_limit(ntree));
    ptn_like = aligned_alloc<double>(mem_size);
    _pattern_scaling = aligned_alloc<double>(mem_size);

    // get the pattern frequencies
    for (i=0; i<nptn; i++)
        ptn_freq[i] = aln->at(i).frequency;
    
    // get whether the pattern is constant
    for (i=0; i<nptn; i++) {
        patn_isconst[i] = aln->at(i).isConst();
    }
    
    // get number of parsimony informative sites
    ninformsite = 0;
    for (i=0; i<nptn; i++) {
        if (aln->at(i).isInformative())
            ninformsite+=aln->at(i).frequency;
    }
    
    // number of optimization steps, default: number of Trees * 2
    // optimize_steps = 2 * size();
    optimize_steps = OPTIMIZE_STEPS;
    
    // initialize the variables
    isTreeWeightFixed = false;
    weightGrpExist = false;
    isEdgeLenRestrict = false;
    parsi_computed = false;
    optim_type = 1;
    ntip = aln->getNSeq();
    if (ntree > 0)
        nbranch = 2 * ntip - 3; // number of branches in a unrooted tree
    else
        nbranch = 0;
    logl_variance = 0.0;
    isLinkModel = true;
    isLinkSiteRate = true;
    anySiteRate = false;
    isNestedOpenmp = false;
    rhas_var = NULL;
}

IQTreeMix::~IQTreeMix() {
    size_t i;
    
    // restore back the initial models and site rates
    for (i=0; i<size(); i++) {
        at(i)->getModelFactory()->model = initial_models[i];
        at(i)->setModel(initial_models[i]);
        at(i)->getModel()->setTree(at(i));
        at(i)->getModelFactory()->site_rate = initial_site_rates[i];
        at(i)->setRate(initial_site_rates[i]);
        at(i)->getRate()->setTree(at(i));
    }

    model_factory = NULL;
    model = NULL;
    site_rate = NULL;
    
    for (i=0; i<size(); i++) {
        delete (at(i));
    }
    if (ptn_like_cat != NULL) {
        aligned_free(ptn_like_cat);
    }
    if (_ptn_like_cat != NULL) {
        aligned_free(_ptn_like_cat);
    }
    if (ptn_scale_cat != NULL){
        aligned_free(ptn_scale_cat);
    }
    if (ptn_freq != NULL) {
        aligned_free(ptn_freq);
    }
    if (patn_isconst != NULL) {
        aligned_free(patn_isconst);
    }
    if (patn_parsimony != NULL) {
        aligned_free(patn_parsimony);
    }
    if (single_ptn_tree_like != NULL) {
        aligned_free(single_ptn_tree_like);
    }
    if (ptn_like != NULL) {
        aligned_free(ptn_like);
    }
    if (_pattern_scaling != NULL) {
        aligned_free(_pattern_scaling);
    }
    if (rhas_var != NULL) {
        aligned_free(rhas_var);
    }
}

void separateStr(string str, vector<string>& substrs, char separator) {
    int startpos = 0;
    int pos = 0;
    int brac_num = 0;
    while (pos < str.length()) {
        if (str[pos] == '{')
            brac_num++;
        else if (str[pos] == '}')
            brac_num--;
        else if (str[pos] == separator && brac_num<=0) {
            if (pos - startpos > 0)
                substrs.push_back(str.substr(startpos, pos-startpos));
            brac_num = 0;
            startpos = pos+1;
        }
        pos++;
    }
    if (pos - startpos > 0)
        substrs.push_back(str.substr(startpos, pos-startpos));
}

void divideModelNSiteRate(string name, string& model, string& siteRate) {
    string orig, s;
    size_t plus_pos;
    int i;
    // assume model always exists
    
    orig = name;
    model = "";
    siteRate = "";
    i = 0;
    while (name.length() > 0) {
        plus_pos = name.find("+");
        if (plus_pos != string::npos) {
            s = name.substr(0,plus_pos);
            name = name.substr(plus_pos+1);
        } else {
            s = name;
            name = "";
        }
        if (s.length() == 0) {
            outError(orig + " is not a valid model");
        }
        if (isRHS(s)) {
            if (siteRate.length() > 0)
                siteRate.append("+");
            siteRate.append(s);
        } else {
            if (model.length() > 0)
                model.append("+");
            model.append(s);
        }
        i++;
    }
    if (siteRate.length() == 0)
        siteRate = "E";
}

void rmSpace(string& s) {
    size_t i,k;
    k=0;
    for (i=0; i<s.length(); i++) {
        if (s[i]!=' ') {
            if (k < i) {
                s[k] = s[i];
            }
            k++;
        }
    }
    if (k == 0)
        s = "";
    else if (k < s.length())
        s = s.substr(0,k);
}

// to separate the submodel names and the site rate names from the full model name
void IQTreeMix::separateModel(string modelName) {
    size_t t_pos, t_pos2, i, k, d_pos, d_pos2;
    string s;
    vector<string> model_array, submodel_array, weight_array;
    
    rmSpace(modelName);
    treemix_model = modelName;
    model_names.clear();
    siterate_names.clear();
    
    // check how many trees
    t_pos = modelName.rfind("+T");
    if (t_pos == string::npos) {
        outError("This model is not a tree mixture model, because there is no '+T'");
    }
    t_pos2 = t_pos+2;
    if (t_pos2 < modelName.length()) {
        // if "+TR", then it is a edge-length-restricted model in which the edges of different trees having the same partition have the same lengths
        if (modelName[t_pos2] == 'R') {
            cout << "This is a edge-length-constrained model, in which the edges of different trees having the same partition are constrained to have the same lengths" << endl;
            isEdgeLenRestrict = true;
            t_pos2++;
        }
        d_pos = t_pos2;
        while (d_pos < modelName.length() && isdigit(modelName[d_pos]))
            d_pos++;
        if (d_pos - t_pos2 > 0) {
            ntree = atoi(modelName.substr(t_pos2, d_pos-t_pos2).c_str());
            if (ntree < size()) {
                cout << endl << "NOTE: Only the first " << ntree << " trees in the treefile are considered, because '" <<  modelName.substr(t_pos+1) << "' has been specified!" << endl;
                resize(ntree);
                weights.resize(ntree);
                weight_logs.resize(ntree);
                double init_weight = 1.0 / (double) ntree;
                double init_weight_log = log(init_weight);
                for (i=0; i<ntree; i++) {
                    weights[i] = init_weight;
                    weight_logs[i] = init_weight_log;
                }
            }
        }
        // check whether any tree weight is specified
        if (d_pos < modelName.length() && modelName[d_pos] == '{') {
            d_pos2 = modelName.find_first_of('}', d_pos+1);
            if (d_pos2 != string::npos) {
                s = modelName.substr(d_pos+1, d_pos2-d_pos-1);
                separateStr(s, weight_array, ',');
                if (weight_array.size() != ntree) {
                    outError("The number of specified tree weights does not match with the number of trees");
                }
                double sum = 0.0;
                for (i=0; i<ntree; i++) {
                    weights[i] = atof(weight_array[i].c_str());
                    sum += weights[i];
                }
                // normalize the weights
                for (i=0; i<ntree; i++) {
                    weights[i] = weights[i] / sum;
                    weight_logs[i] = log(weights[i]);
                }
                isTreeWeightFixed = true;
                cout << "The fixed tree weights:";
                for (i=0; i<ntree; i++) {
                    cout << " " << weights[i];
                }
                cout << endl;
            }
        }
        // check whether any grouping of tree weights is specified
        if (d_pos < modelName.length() && modelName[d_pos] == '[') {
            map<string,int> weight_grps;
            map<string,int>::iterator itr;
            weight_group_member.clear();
            d_pos2 = modelName.find_first_of(']', d_pos+1);
            if (d_pos2 != string::npos) {
                s = modelName.substr(d_pos+1, d_pos2-d_pos-1);
                separateStr(s, weight_array, ',');
                if (weight_array.size() != ntree) {
                    outError("The number of specified tree weights does not match with the number of trees");
                }
                for (i=0; i<ntree; i++) {
                    itr = weight_grps.find(weight_array[i]);
                    int grpID;
                    if (itr != weight_grps.end()) {
                        grpID = itr->second;
                        weightGrpExist = true; // more than one tree weights belong to the same group
                    } else {
                        grpID = weight_grps.size();
                        weight_grps.insert(pair<string,int>(weight_array[i],grpID));
                        weight_group_member.push_back(vector<int>());
                    }
                    weight_group_member[grpID].push_back(i);
                }
            }
        }
    }
    
    if (weight_group_member.size() == 0) {
        for (i=0; i<ntree; i++) {
            vector<int> new_grp;
            new_grp.push_back(i);
            weight_group_member.push_back(new_grp);
        }
    }
    
    if (ntree <= 1) {
        outError("For tree mixture model, number of trees has to be at least 2.");
    }
    
    // remove the '+Txx'
    modelName = modelName.substr(0, t_pos);
    
    // break the whole name according to '+'
    separateStr(modelName, model_array, '+');
    
    // check each model/siterate
    for (i=0; i<model_array.size(); i++) {
        s = model_array[i];
        if (s.length() == 0) {
            continue;
        } else if (s.length() > 6 && s.substr(0,5) == "TMIX{" && s.substr(s.length()-1,1) == "}") {
            // mixture model
            s = s.substr(5,s.length()-6); // remove the beginning "MIX{" and the ending "}"
            if (i==0) {
                // unlinked substitution models (while site rates may or may not be linked)
                bool siteRateAppear = false;
                string curr_model, curr_siterate;
                separateStr(s, submodel_array, ',');
                for (k=0; k<submodel_array.size(); k++) {
                    divideModelNSiteRate(submodel_array[k], curr_model, curr_siterate);
                    model_names.push_back(curr_model);
                    siterate_names.push_back(curr_siterate);
                    if (curr_siterate.length() > 0 && curr_siterate != "E") {
                        siteRateAppear = true; // some site rates appear
                    }
                }
                if (!siteRateAppear) {
                    // all siterate_names are empty, thus remove them
                    siterate_names.clear();
                }
            } else if (siterate_names.size()==0) {
                // unlinked site rates
                separateStr(s, submodel_array, ',');
                string curr_model, curr_siterate;
                for (k=0; k<submodel_array.size(); k++) {
                    if (submodel_array[k].length() == 0) {
                        outError("There is an empty submodel.");
                    }
                    divideModelNSiteRate(submodel_array[k], curr_model, curr_siterate);
                    if (curr_model.length() > 0) {
                        outError("The model: " + treemix_model + " is not correctly specified.");
                    }
                    siterate_names.push_back(submodel_array[k]);
                }
            } else {
                outError("The model: " + treemix_model + " is not correctly specified. Are you using more than one 'TMIX'?");
            }
        } else if (i==0) {
            // linked substitution model
            // assuming the first is the substitution model
            model_names.push_back(s);
        } else if (s.length() > 0 && isRHS(s)) {
            // linked RHS model
            if (siterate_names.size() > 1) {
                outError("'" + s + "' is linked, but the site rates are unlinked");
            } else if (siterate_names.size() == 1) {
                siterate_names[0].append("+" + s);
            } else {
                siterate_names.push_back(s);
            }
        } else {
            // other linked element which belongs to a part of the substutation model
            if (model_names.size() > 1) {
                outError("'" + s + "' is linked, but the model is unlinked");
            } else if (model_names.size() == 1){
                model_names[0].append("+" + s);
            } else {
                outError("'" + s + "' appears before the model defined");
            }
        }
    }
    if (model_names.size() == 0) {
        outError("It seems no model is defined.");
    }
    if (siterate_names.size() == 0) {
        anySiteRate = false;
    } else {
        anySiteRate = true;
    }
    
    // check correctness
    if (model_names.size() > 1 && model_names.size() != ntree) {
        outError("The number of submodels specified in the mixture does not match with the tree number");
    }
    if (siterate_names.size() > 1 && siterate_names.size() != ntree) {
        outError("The number of site rates specified in the mixture does not match with the tree number");
    }

    // build the branch ID
    computeBranchID();
    
    // check whether it is linked or unlinked model or rate
    if (model_names.size() == 1) {
        isLinkModel = true;
    } else {
        isLinkModel = false;
    }
    if (anySiteRate) {
        if (siterate_names.size() == 1) {
            isLinkSiteRate = true;
        } else {
            isLinkSiteRate = false;
        }
    }

    /*
    // show summary
    cout << endl;
    if (isLinkModel) {
        cout << "Linked substitution model:" << endl;
    } else {
        cout << "Unlinked substitution models:" << endl;
    }
    for (i = 0; i < model_names.size(); i++) {
        cout << "   " << model_names[i] << endl;
    }
    cout << endl;
    if (anySiteRate) {
        if (isLinkSiteRate) {
            cout << "Linked RHS model:" << endl;
        } else {
            cout << "Unlinked RHS models:" << endl;
        }
        for (i = 0; i < siterate_names.size(); i++) {
            cout << "   " << siterate_names[i] << endl;
        }
    }
    cout << endl;
    */
}

void IQTreeMix::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    size_t i;
    string curr_model;

    models.clear();
    site_rates.clear();
    site_rate_trees.clear();
    separateModel(model_name);

    // initialize the models
    for (i=0; i<ntree; i++) {
        if (isLinkModel) {
            curr_model = model_names[0];
        } else {
            curr_model = model_names[i];
        }
        if (anySiteRate) {
            params.optimize_alg_gammai = "BFGS";
            params.optimize_alg_freerate = "2-BFGS";
            if (isLinkSiteRate) {
                if (siterate_names[0] != "E")
                    curr_model += "+" + siterate_names[0];
            } else {
                params.optimize_alg_mixlen = "BFGS";
                if (siterate_names[i] != "E")
                    curr_model += "+" + siterate_names[i];
            }
        }
        // cout << "model: " << curr_model << endl;
        at(i)->initializeModel(params, curr_model, models_block);
    }
    
    // save the initial models and site rates
    for (i=0; i<ntree; i++) {
        initial_models.push_back(at(i)->getModelFactory()->model);
        initial_site_rates.push_back(at(i)->getModelFactory()->site_rate);
    }
    
    // handle the linked or unlinked substitution model(s)
    if (isLinkModel) {
        models.push_back(at(0)->getModelFactory()->model);
        for (i=1; i<ntree; i++) {
            at(i)->getModelFactory()->model = models[0];
            at(i)->setModel(models[0]);
        }
        // for linked subsitution model, set its tree to this tree
        models[0]->setTree(this);
    } else {
        for (i=0; i<ntree; i++) {
            models.push_back(at(i)->getModelFactory()->model);
        }
        // for unlinked subsitution models, set their trees to the corresponding trees
        for (i=0; i<ntree; i++) {
            models[i]->setTree(at(i));
        }
    }
    
    // handle the linked or unlinked site rate(s)
    if (anySiteRate) {
        for (i=0; i<ntree; i++) {
            site_rates.push_back(at(i)->getModelFactory()->site_rate);
        }
        if (isLinkSiteRate) {
            // for linked site rate model, set their trees to this tree
            for (i=0; i<ntree; i++) {
                site_rates[i]->setTree(this);
                site_rate_trees.push_back(at(i));
            }
        } else {
            // for unlinked site rate model, set their trees to corresponding trees
            for (i=0; i<ntree; i++) {
                site_rates[i]->setTree(at(i));
                site_rate_trees.push_back(at(i));
            }
        }
    }
}

int IQTreeMix::getNumLhCat(SiteLoglType wsl) {
    int ncat = 0;
    switch (wsl) {
    case WSL_NONE: ASSERT(0 && "is not WSL_NONE"); return 0;
    case WSL_SITE: ASSERT(0 && "is not WSL_SITE"); return 0;
    case WSL_MIXTURE_RATECAT:
        ncat = getRate()->getNDiscreteRate();
        if (getModel()->isMixture() && !getModelFactory()->fused_mix_rate)
            ncat *= getModel()->getNMixtures();
        return ncat;
    case WSL_RATECAT:
        return getRate()->getNDiscreteRate();
    case WSL_MIXTURE:
        return getModel()->getNMixtures();
    case WSL_TMIXTURE:
        return size();
    }
    return ncat;
}

// compute the log-likelihood values for every site and tree
// updated array: _ptn_like_cat
// update_which_tree: only that tree has been updated
void IQTreeMix::computeSiteTreeLogLike(int update_which_tree) {
    // cout << "enter IQTreeMix::computeSiteTreeLogLike" << endl;
    // cout << "update_which_tree = " << update_which_tree << endl;
    PhyloTree* ptree;
    int k,t;

    t = update_which_tree;
    if (t == -1) {
        IQTreeMix::computeLikelihood();
        return;
    }
    if (isLinkSiteRate && t > 0) {
        // Store the RHAS variables of tree 0 to the array rhas_var
        storeTree0RHAS();
    }
    
    // compute likelihood for each tree
    double* patternlh_tree = _ptn_like_cat + t*nptn;
    ptree = at(t)->getRate()->getTree();
    // set the tree t as the site rate's tree
    // and compute the likelihood values
    at(t)->getRate()->setTree(at(t));
    if (isLinkSiteRate && t > 0) {
        // Replace the RHAS variables of tree t by those of tree 0
        copyRHASfrTree0(t);
    }
    at(t)->initializeAllPartialLh();
    at(t)->clearAllPartialLH();
    at(t)->computeLikelihood(patternlh_tree, false);
    // set back the prevoius site rate's tree
    at(t)->getRate()->setTree(ptree);

    // reorganize the array
    k=t;
    for (size_t ptn=0; ptn<nptn; ptn++) {
        ptn_like_cat[k] = patternlh_tree[ptn];
        k+=ntree;
    }

    // synchonize the scaling factor between the updated tree and the other trees for every pattern
    #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
    for (size_t ptn=0; ptn<nptn; ptn++) {
        double* pattern_lh_tree = ptn_like_cat + ntree * ptn;
        // check whether the scaling factor of the updated tree has much difference
        double scale_diff = at(t)->_pattern_scaling[ptn] - this->_pattern_scaling[ptn];
        double abs_scale_diff = fabs(scale_diff);
        if (abs_scale_diff > TINY_SCALE_DIFF) {
            // the difference in the scaling factor is significant
            if (scale_diff > 0) {
                // the scaling factor of the updated tree is larger
                if (abs_scale_diff > ONE_LOG_SCALE_DIFF) {
                    // the difference is more than one scaling threshold, which is 2^{256}
                    // then the likelihoods of the other trees are insignificant
                    for (int j=0; j<ntree; j++) {
                        if (j == t)
                            continue;
                        pattern_lh_tree[j] = 0.0;
                    }
                } else {
                    for (int j=0; j<ntree; j++) {
                        if (j == t)
                            continue;
                        pattern_lh_tree[j] *= SCALING_THRESHOLD;
                    }
                }
                this->_pattern_scaling[ptn] = at(t)->_pattern_scaling[ptn];
            } else {
                // the scaling factor of the updated tree is smaller
                if (abs_scale_diff > ONE_LOG_SCALE_DIFF) {
                    // the difference is more than one scaling threshold, which is 2^{256}
                    // then the likelihood is insignificant
                    pattern_lh_tree[t] = 0.0;
                } else {
                    pattern_lh_tree[t] *= SCALING_THRESHOLD;
                }
            }
        }
    }

    // combining all the existing likelihood values of the trees
    // computeLikelihood_combine();
}

double IQTreeMix::computeLikelihood(double *pattern_lh, bool save_log_value) {
    // size_t i,j;
    double logLike = 0.0;
    double subLike;
    double score;
    
    if (isLinkSiteRate) {
        // Store the RHAS variables of tree 0 to the array rhas_var
        storeTree0RHAS();
    }

    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(1);
        omp_set_max_active_levels(2);
    }
    #endif

    // compute likelihood for each tree
    #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
    for (size_t t=0; t<ntree; t++) {
        if (isNestedOpenmp) {
            #ifdef _OPENMP
            omp_set_num_threads(at(t)->num_threads);
            #endif
        }
        double* pattern_lh_tree = _ptn_like_cat + nptn * t;
        // save the site rate's tree
        PhyloTree* ptree = at(t)->getRate()->getTree();
        // set the tree t as the site rate's tree
        // and compute the likelihood values
        at(t)->getRate()->setTree(at(t));
        
        if (isLinkSiteRate && t > 0) {
            // Replace the RHAS variables of tree t by those of tree 0
            copyRHASfrTree0(t);
        }
        
        at(t)->initializeAllPartialLh();
        at(t)->clearAllPartialLH();
        at(t)->computeLikelihood(pattern_lh_tree, false);
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
    }

    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(num_threads);
    }
    #endif

    // reorganize the array
    // #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
    for (size_t t=0; t<ntree; t++) {
        size_t i = t * nptn;
        size_t j = t;
        for (size_t ptn=0; ptn<nptn; ptn++) {
            // cout << setprecision(7) << t << "\t" << ptn << "\t" << _ptn_like_cat[i] << endl;
            ptn_like_cat[j] = _ptn_like_cat[i];
            ptn_scale_cat[j] = at(t)->_pattern_scaling[ptn];
            i++;
            j+=ntree;
        }
    }

    // synchonize the scaling factor among the trees for every pattern
    #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
    for (size_t ptn=0; ptn<nptn; ptn++) {
        double* pattern_lh_tree = ptn_like_cat + ntree * ptn;
        double* pattern_scale_tree = ptn_scale_cat + ntree * ptn;
        // find the max scaling factor among the trees
        double max_scale = pattern_scale_tree[0];
        int max_tree = 0;
        for (size_t t=1; t<ntree; t++) {
            if (max_scale < pattern_scale_tree[t]) {
                max_scale = pattern_scale_tree[t];
                max_tree = t;
            }
        }
        // synchonize the scaling factor among the trees
        for (size_t t=0; t<ntree; t++) {
            if (t == max_tree)
                continue;
            double scale_diff = max_scale - pattern_scale_tree[t];
            if (scale_diff > TINY_SCALE_DIFF) {
                if (scale_diff > ONE_LOG_SCALE_DIFF) {
                    // the difference is more than one scaling threshold, which is 2^{256}
                    // then the likelihood is insignificant
                    pattern_lh_tree[t] = 0.0;
                } else {
                    pattern_lh_tree[t] *= SCALING_THRESHOLD;
                }
            }
        }
        _pattern_scaling[ptn] = max_scale;
    }
    
    // compute the overall likelihood value by combining all the existing likelihood values of the trees
    return computeLikelihood_combine(pattern_lh, save_log_value);
}

double IQTreeMix::computePatternLhCat(SiteLoglType wsl) {
    // cout << "enter IQTreeMix::computePatternLhCat" << endl;
    // this function only supports wsl = WSL_MIXTURE or WSL_TMIXTURE
    
    // double* pattern_lh_tree;
    // double* pattern_lh_cat;
    //PhyloTree* ptree;
    // size_t nmix, p, m, i, j, ptn;
    size_t nmix;
    size_t tot_mem_size;
    double score;
    int index, indexlh;
    bool save_log_value = false;
    // double scaling_factor;

    int numStates = at(0)->model->num_states;
    size_t mem_size = get_safe_upper_limit(getAlnNPattern()) + max(get_safe_upper_limit(numStates),
        get_safe_upper_limit(at(0)->model_factory->unobserved_ptns.size()));
    if (!_pattern_lh_cat) {
        if (wsl == WSL_TMIXTURE)
            tot_mem_size = mem_size * ntree;
        else
            tot_mem_size = mem_size * at(0)->site_rate->getNDiscreteRate() * ((at(0)->model_factory->fused_mix_rate)? 1 : at(0)->model->getNMixtures());
        _pattern_lh_cat = aligned_alloc<double>(tot_mem_size);
    }

    if (isLinkSiteRate) {
        // Store the RHAS variables of tree 0 to the array rhas_var
        storeTree0RHAS();
    }

    if (wsl == WSL_TMIXTURE) {
        // compute likelihood for each tree
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            // omp_set_nested(1);
            omp_set_max_active_levels(2);
        }
        #endif
        #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
        for (size_t t=0; t<ntree; t++) {
            #ifdef _OPENMP
            if (isNestedOpenmp) {
                omp_set_num_threads(at(t)->num_threads);
            }
            #endif
            double* pattern_lh_tree = _ptn_like_cat + t * nptn;
            // save the site rate's tree
            PhyloTree* ptree = at(t)->getRate()->getTree();
            // set the tree t as the site rate's tree
            // and compute the likelihood values
            at(t)->getRate()->setTree(at(t));
            if (isLinkSiteRate && t > 0) {
                // Replace the RHAS variables of tree t by those of tree 0
                copyRHASfrTree0(t);
            }
            at(t)->initializeAllPartialLh();
            at(t)->clearAllPartialLH();
            at(t)->computeLikelihood(pattern_lh_tree, save_log_value);
            // set back the prevoius site rate's tree
            at(t)->getRate()->setTree(ptree);
        }
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            // omp_set_nested(0);
            omp_set_max_active_levels(1);
            omp_set_num_threads(num_threads);
        }
        #endif

        // reorganize the array
        // #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
        for (size_t t=0; t<ntree; t++) {
            int i = t * nptn;
            int j = t;
            for (size_t ptn=0; ptn<nptn; ptn++) {
                ptn_like_cat[j] = _ptn_like_cat[i];
                ptn_scale_cat[j] = at(t)->_pattern_scaling[ptn];
                i++;
                j += ntree;
            }
        }

        // synchonize the scaling factor among the trees for every pattern
        #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
        for (size_t ptn=0; ptn<nptn; ptn++) {
            // find the max scaling factor among the trees
            double* pattern_lh_tree = ptn_like_cat + ptn * ntree;
            double* pattern_scale_tree = ptn_scale_cat + ptn * ntree;
            double max_scale = at(0)->_pattern_scaling[ptn];
            int max_tree = 0;
            for (size_t t=1; t<ntree; t++) {
                if (max_scale < at(t)->_pattern_scaling[ptn]) {
                    max_scale = at(t)->_pattern_scaling[ptn];
                    max_tree = t;
                }
            }
            // synchonize the scaling factor among the trees
            for (size_t t=0; t<ntree; t++) {
                if (t == max_tree)
                    continue;
                double scale_diff = max_scale - pattern_scale_tree[t];
                if (scale_diff > TINY_SCALE_DIFF) {
                    if (scale_diff > ONE_LOG_SCALE_DIFF) {
                        // the difference is more than one scaling threshold, which is 2^{256}
                        // then the likelihood is insignificant
                        // WARNING!! There must be a minimum threshold for the tree weight
                        // this approximation will be incorrect if there are some tree weights close to zero
                        pattern_lh_tree[t] = 0.0;
                    } else {
                        pattern_lh_tree[t] *= SCALING_THRESHOLD;
                    }
                }
            }
            _pattern_scaling[ptn] = max_scale;
        }
        
        // computing the array _pattern_lh_cat
        #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
        for (size_t ptn=0; ptn<nptn; ptn++) {
            size_t idx = ptn * ntree;
            double* pattern_lh_tree = ptn_like_cat + idx;
            double* pattern_lh_cat = _pattern_lh_cat + idx;
            double scaling_factor = exp(_pattern_scaling[ptn]);
            for (size_t t=0; t<ntree; t++) {
                pattern_lh_cat[t] = pattern_lh_tree[t] * weights[t] * scaling_factor;
            }
        }
        
        // combining all the existing likelihood values of the trees
        score = computeLikelihood_combine();

    } else {
        // compute _pattern_lh_cat for each tree
        nmix = at(0)->getModel()->getNMixtures();
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            // omp_set_nested(1);
            omp_set_max_active_levels(2);
        }
        #endif
        #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
        for (size_t t = 0; t < ntree; t++) {
            #ifdef _OPENMP
            if (isNestedOpenmp) {
                omp_set_num_threads(at(t)->num_threads);
            }
            #endif
            if (isLinkSiteRate && t > 0) {
                // Replace the RHAS variables of tree t by those of tree 0
                copyRHASfrTree0(t);
            }
            if (t > 0 && nmix != at(t)->getModel()->getNMixtures()) {
                outError("The number of classes inside each tree is not the same!");
            }
            at(t)->computePatternLhCat(wsl);
        }
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            // omp_set_nested(0);
            omp_set_max_active_levels(1);
            omp_set_num_threads(num_threads);
        }
        #endif

        // compute the overall _pattern_lh_cat
        #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
        for (size_t p = 0; p < nptn * nmix; p++) {
            double ptnLike = 0.0;
            for (size_t t = 0; t < ntree; t++) {
                ptnLike += at(t)->_pattern_lh_cat[p] * weights[t];
            }
            _pattern_lh_cat[p] = ptnLike;
        }
        
        score = computeLikelihood();
    }


    return score;
}

// compute the overall likelihood value by combining all the existing likelihood values of the trees
double IQTreeMix::computeLikelihood_combine(double *pattern_lh, bool save_log_value) {
    double* pattern_lh_tree;
    // size_t i,ptn,t;
    double logLike = 0.0;
    // double subLike;
    double ptnLike;
    PhyloTree* ptree;
    
    // compute the total likelihood
    #pragma omp parallel for schedule(static) num_threads(num_threads) reduction(+:logLike) if (num_threads > 1)
    for (size_t ptn=0; ptn<nptn; ptn++) {
        size_t i = ptn * ntree;
        double subLike = 0.0;
        for (size_t t=0; t<ntree; t++) {
            subLike += ptn_like_cat[i] * weights[t];
            i++;
        }
        if (pattern_lh != NULL && !save_log_value) {
            pattern_lh[ptn] = subLike;
        }
        // cout << ptn << "\t" << log(subLike) << "\t" << ptn_freq[ptn] << endl;
        double ptnLike = log(subLike) + _pattern_scaling[ptn];
        if (pattern_lh != NULL && save_log_value) {
            pattern_lh[ptn] = ptnLike;
        }
        ptnLike = ptnLike * ptn_freq[ptn];
        logLike += ptnLike;
    }
    return logLike;
}

/*
// compute the overall likelihood value by combining all the existing likelihood values of the trees
double IQTreeMix::computeLikelihood_combine(double *pattern_lh) {
    double* pattern_lh_tree;
    size_t i,ptn,t;
    double logLike = 0.0;
    double subLike;
    double ptnLike;
    double maxLnLike;
    int max_t;
    PhyloTree* ptree;
    
    // compute the total likelihood
    i=0;
    for (ptn=0; ptn<nptn; ptn++) {
        for (t=0; t<ntree; t++) {
            single_ptn_tree_like[t] = ptn_like_cat[i] + weight_logs[t];
            i++;
        }
        // for the max
        maxLnLike = single_ptn_tree_like[0];
        max_t = 0;
        for (t=1; t<ntree; t++) {
            if (single_ptn_tree_like[t] > maxLnLike) {
                maxLnLike = single_ptn_tree_like[t];
                max_t = t;
            }
        }
        // sum of all values of x, exp(x - max_t)
        ptnLike = 1.0;
        for (t=0; t<ntree; t++) {
            if (t==max_t)
                continue;
            ptnLike += exp(single_ptn_tree_like[t] - maxLnLike);
        }
        ptnLike = log(ptnLike) + maxLnLike;
        
        // cout << ptn << "\t" << ptnlike << "\t" << ptn_freq[ptn] << endl;
        logLike += ptnLike * (double) ptn_freq[ptn];
        ptn_like[ptn] = ptnLike;
        if (pattern_lh != NULL) {
            pattern_lh[ptn] = ptnLike;
        }
    }

    return logLike;
}
*/


/**
        compute pattern likelihoods only if the accumulated scaling factor is non-zero.
        Otherwise, copy the pattern_lh attribute
        @param pattern_lh (OUT) pattern log-likelihoods,
                        assuming pattern_lh has the size of the number of patterns
        @param cur_logl current log-likelihood (for sanity check)
        @param pattern_lh_cat (OUT) if not NULL, store all pattern-likelihood per category
 */
void IQTreeMix::computePatternLikelihood(double *pattern_lh, double *cur_logl,
                                         double *pattern_lh_cat, SiteLoglType wsl) {
    // cout << "enter IQTreeMix::computePatternLikelihood" << endl;
    // double* pattern_lh_tree;
    // double* pattern_scale_tree;
    // double logLike = 0.0;
    // double subLike;
    double score;
    //PhyloTree* ptree;
    // bool rhas_changed;

    if (isLinkSiteRate) {
        // Store the RHAS variables of tree 0 to the array rhas_var
        storeTree0RHAS();
    }
    
    // compute likelihood for each tree
    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(1);
        omp_set_max_active_levels(2);
    }
    #endif
    #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
    for (size_t t=0; t<ntree; t++) {
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            omp_set_num_threads(at(t)->num_threads);
        }
        #endif
        double* pattern_lh_tree = _ptn_like_cat + t * nptn;
        // save the site rate's tree
        PhyloTree* ptree = at(t)->getRate()->getTree();
        // set the tree t as the site rate's tree
        // and compute the likelihood values
        at(t)->getRate()->setTree(at(t));
        if (isLinkSiteRate && t > 0) {
            // Replace the RHAS variables of tree t by those of tree 0
            copyRHASfrTree0(t);
        }
        at(t)->initializeAllPartialLh();
        at(t)->clearAllPartialLH();
        at(t)->computeLikelihood(pattern_lh_tree, false);
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
    }
    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(num_threads);
    }
    #endif

    // reorganize the array
    #pragma omp parallel for schedule(dynamic) num_threads(num_threads) if (num_threads > 1)
    for (size_t t=0; t<ntree; t++) {
        int i = t * nptn;
        int j = t;
        for (size_t ptn=0; ptn<nptn; ptn++) {
            ptn_like_cat[j] = _ptn_like_cat[i];
            ptn_scale_cat[j] = at(t)->_pattern_scaling[ptn];
            if (pattern_lh_cat)
                pattern_lh_cat[j] = log(ptn_like_cat[j]) + ptn_scale_cat[j];
            i++;
            j+=ntree;
        }
    }

    // synchonize the scaling factor among the trees for every pattern
    #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
    for (size_t ptn=0; ptn<nptn; ptn++) {
        // find the max scaling factor among the trees
        double* pattern_lh_tree = ptn_like_cat + ptn * ntree;
        double* pattern_scale_tree = ptn_scale_cat + ptn * ntree;
        double max_scale = pattern_scale_tree[0];
        int max_tree = 0;
        for (size_t t=1; t<ntree; t++) {
            if (max_scale < pattern_scale_tree[t]) {
                max_scale = pattern_scale_tree[t];
                max_tree = t;
            }
        }
        // synchonize the scaling factor among the trees
        for (size_t t=0; t<ntree; t++) {
            if (t == max_tree)
                continue;
            double scale_diff = max_scale - pattern_scale_tree[t];
            if (scale_diff > TINY_SCALE_DIFF) {
                if (scale_diff > ONE_LOG_SCALE_DIFF) {
                    // the difference is more than one scaling threshold, which is 2^{256}
                    // then the likelihood is insignificant
                    pattern_lh_tree[t] = 0.0;
                } else {
                    pattern_lh_tree[t] *= SCALING_THRESHOLD;
                }
            }
        }
        _pattern_scaling[ptn] = max_scale;
    }

    score = computeLikelihood_combine(pattern_lh);
    
    // compute the overall likelihood value by combining all the existing likelihood values of the trees
    if (cur_logl != NULL) {
        *cur_logl = score;
    }
}

void IQTreeMix::initializeAllPartialLh() {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->initializeAllPartialLh();
    }
}

void IQTreeMix::deleteAllPartialLh() {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->deleteAllPartialLh();
    }
}

void IQTreeMix::clearAllPartialLH(bool make_null) {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->clearAllPartialLH(make_null);
    }
}

/**
        compute pattern posterior probabilities per rate/mixture category
        @param pattern_prob_cat (OUT) all pattern-probabilities per category
        @param wsl either WSL_RATECAT, WSL_MIXTURE or WSL_MIXTURE_RATECAT
 */
void IQTreeMix::computePatternProbabilityCategory(double *ptn_prob_cat, SiteLoglType wsl) {

    if (wsl != WSL_TMIXTURE) {
        // Not a tree mixture model
        PhyloTree::computePatternProbabilityCategory(ptn_prob_cat, wsl);
        return;
    }

    bool need_computeLike = true;
    int update_which_tree = -1;
    bool need_multiplyFreq = false;
    getPostProb(ptn_prob_cat, need_computeLike, update_which_tree, need_multiplyFreq);
}

/**
        optimize all branch lengths of one tree
        @param iterations number of iterations to loop through all branches
 */
void IQTreeMix::optimizeAllBranchesOneTree(int whichtree, int my_iterations, double tolerance, int maxNRStep) {
    PhyloTree* ptree;

    // save the tree of the site rate
    ptree = at(whichtree)->getRate()->getTree();
    at(whichtree)->getRate()->setTree(at(whichtree));
    at(whichtree)->optimizeAllBranches(my_iterations, tolerance, maxNRStep);
    // restore the tree of the site rate
    at(whichtree)->getRate()->setTree(ptree);
}

/**
        optimize all branch lengths of all trees
        @param iterations number of iterations to loop through all branches
        @return the likelihood of the tree
 */
double IQTreeMix::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {

    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(1);
        omp_set_max_active_levels(2);
    }
    #endif
    #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
    for (size_t i=0; i<ntree; i++) {
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            omp_set_num_threads(at(i)->num_threads);
        }
        #endif
        optimizeAllBranchesOneTree(i, my_iterations, tolerance, maxNRStep);
    }
    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(num_threads);
    }
    #endif

    return computeLikelihood();
}

/**
        compute the updated tree weights according to the likelihood values along each site
        prerequisite: computeLikelihood() has been invoked

 */
double IQTreeMix::optimizeBranchLensByBFGS(double gradient_epsilon) {
    int ndim = branch_group.size();
    size_t i;
    double *variables; // used for BFGS numerical recipes
    double *upper_bound;
    double *lower_bound;
    bool *bound_check;
    double score;

    // optimization on which variable
    // 2 - branch lengths
    optim_type = 2;
    
    // collect the branch lengths of the tree
    getBranchLengths(branch_len);

    // compute the likelihood of each tree
    computeLikelihood();
    
    // allocate memory to the arrays
    variables = new double[ndim+1]; // used for BFGS numerical recipes
    upper_bound = new double[ndim+1];
    lower_bound = new double[ndim+1];
    bound_check = new bool[ndim+1];

    // by BFGS algorithm
    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);
    score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, gradient_epsilon);
    getVariables(variables);

    delete[] variables;
    delete[] upper_bound;
    delete[] lower_bound;
    delete[] bound_check;
    
    return score;
}

// save branch lengths of all trees
// node and dad are always NULL
void IQTreeMix::getBranchLengths(vector<DoubleVector> &len, Node *node, Node *dad) {
    size_t i;

    if (len.size() < ntree) {
        len.resize(ntree);
    }
    for (i=0; i<ntree; i++) {
        at(i)->saveBranchLengths(len[i]);
    }
}

// restore branch lengths of all trees
// node and dad are always NULL
void IQTreeMix::setBranchLengths(vector<DoubleVector> &len, Node *node, Node *dad) {
    ASSERT(len.size() == ntree);
    size_t i;

    for (i=0; i<ntree; i++) {
        at(i)->restoreBranchLengths(len[i]);
    }
}

/**
         If there are multiple tree weights belonging to the same group
         set all the tree weights of the same group to their average
 */
void IQTreeMix::checkWeightGrp() {
    size_t i,j;
    double grp_mean;
    if (weightGrpExist) {
        // set all the tree weights in the same group as their average
        for (i = 0; i < weight_group_member.size(); i++) {
            if (weight_group_member[i].size() > 0) {
                grp_mean = 0.0;
                for (j = 0; j < weight_group_member[i].size(); j++) {
                    grp_mean += weights[weight_group_member[i].at(j)];
                }
                grp_mean = grp_mean / weight_group_member[i].size();
                for (j = 0; j < weight_group_member[i].size(); j++) {
                    weights[weight_group_member[i].at(j)] = grp_mean;
                }
            }
        }
    }
}


/**
         If there are multiple branches belonging to the same group
         set all the branches of the same group to their average
 */
void IQTreeMix::checkBranchGrp() {
    size_t i,j;
    size_t treeIdx,branchIdx;
    double grp_len;
    double sum_treeweight;
    
    // collect the branch lengths of the tree
    getBranchLengths(branch_len);
    
    for (i = 0; i < branch_group.size(); i++) {
        sum_treeweight = 0.0;
        for (j = 0; j < branch_group[i].size(); j++) {
            treeIdx = branch_group[i].at(j) / nbranch;
            sum_treeweight += weights[treeIdx];
        }
        grp_len = 0.0;
        for (j = 0; j < branch_group[i].size(); j++) {
            treeIdx = branch_group[i].at(j) / nbranch;
            branchIdx = branch_group[i].at(j) % nbranch;
            grp_len += branch_len[treeIdx].at(branchIdx) * weights[treeIdx] / sum_treeweight;
        }
        for (j = 0; j < branch_group[i].size(); j++) {
            treeIdx = branch_group[i].at(j) / nbranch;
            branchIdx = branch_group[i].at(j) % nbranch;
            branch_len[treeIdx].at(branchIdx) = grp_len;
        }
    }
    
    // save the updated branch lengths of the tree
    setBranchLengths(branch_len);
}


/**
        compute the updated tree weights according to the likelihood values along each site
        prerequisite: computeLikelihood() has been invoked

 */
double IQTreeMix::optimizeTreeWeightsByEM(double* pattern_mix_lh, double logl_epsilon, int max_steps, bool& tree_weight_converge) {
    size_t ptn, c;
    double *this_lk_cat;
    double lk_ptn;
    double prev_score, score;
    int step;

    prev_score = computeLikelihood();
    
    tmp_weights.resize(ntree);
    // save the previous weights
    for (c=0; c<ntree; c++)
        tmp_weights[c] = weights[c];

    for (step = 0; step < max_steps || max_steps == -1; step++) {
        
        getPostProb(pattern_mix_lh, false);

        // reset the weights
        for (c = 0; c < ntree; c++) {
            weights[c] = 0.0;
        }

        // E-step
        for (ptn = 0; ptn < nptn; ptn++) {
            this_lk_cat = pattern_mix_lh + ptn*ntree;
            for (c = 0; c < ntree; c++) {
                weights[c] += this_lk_cat[c];
            }
        }

        // M-step
        for (c = 0; c < ntree; c++) {
            weights[c] = weights[c] / getAlnNSite();
            if (weights[c] < 1e-10) weights[c] = 1e-10;
        }
        
        // for weight-restricted condition when weightGrpExist
        checkWeightGrp();
        
        // update the weight_logs
        for (c = 0; c < ntree; c++) {
            weight_logs[c] = log(weights[c]);
        }
    
        // if all differences between the new weights and the old weights <= WEIGHT_EPSILON, then tree_weight_converge = true
        tree_weight_converge = true;
        for (c = 0; c< ntree && tree_weight_converge; c++)
            if (fabs(weights[c] - tmp_weights[c]) > WEIGHT_EPSILON)
                tree_weight_converge = false;

        score = computeLikelihood_combine();

        if (score < prev_score + logl_epsilon) {
            // converged
            break;
        }

        prev_score = score;

        // save the previous weights
        for (c=0; c<ntree; c++)
            tmp_weights[c] = weights[c];

    }
    
    /*
    // show the values of weights
    cout << "weights: ";
    for (c=0; c<ntree; c++) {
        if (c>0)
            cout << ",";
        cout << weights[c];
    }
    cout << endl;
    */
    
    return score;
}

/**
        compute the updated tree weights according to the likelihood values along each site
        prerequisite: computeLikelihood() has been invoked

 */
double IQTreeMix::optimizeTreeWeightsByBFGS(double gradient_epsilon) {
    int ndim = weight_group_member.size();
    size_t i;
    double *variables; // used for BFGS numerical recipes
    double *upper_bound;
    double *lower_bound;
    bool *bound_check;
    double score;

    // optimization on which variable
    // 1 - tree weights
    optim_type = 1;

    // special case: ndim = 1, i.e. all tree weights are forced the same
    if (ndim == 1) {
        double w = 1.0 / size();
        double w_log = log(w);
        for (i=0; i<size(); i++) {
            weights[i] = w;
            weight_logs[i] = w_log;
        }
        return computeLikelihood();
    }
    
    // compute the likelihood of each tree
    computeLikelihood();
    
    // allocate memory to the arrays
    variables = new double[ndim+1]; // used for BFGS numerical recipes
    upper_bound = new double[ndim+1];
    lower_bound = new double[ndim+1];
    bound_check = new bool[ndim+1];

    // initialize tmp_weights for optimzation
    tmp_weights.resize(ndim);
    for (i=0; i<ndim; i++) {
        tmp_weights[i] = weights[weight_group_member[i].at(0)];
    }
    // by BFGS algorithm
    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);
    score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, gradient_epsilon);
    getVariables(variables);

    delete[] variables;
    delete[] upper_bound;
    delete[] lower_bound;
    delete[] bound_check;
    
    /*
    //show the tree weights
    cout << "Tree weights: ";
    for (i=0; i<size(); i++) {
        cout << weights[i] << ";";
    }
    cout << endl;
    */
    
    return score;
}

void IQTreeMix::showTree() {
    size_t i;
    for (i=0; i<size(); i++) {
        cout << "Tree " << i+1 << ": ";
        at(i)->printTree(cout);
        cout << endl;
    }
}

void IQTreeMix::setRootNode(const char *my_root, bool multi_taxa) {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->setRootNode(my_root, multi_taxa);
    }
}

/**
    set checkpoint object
    @param checkpoint
*/
void IQTreeMix::setCheckpoint(Checkpoint *checkpoint) {
    size_t i;
    IQTree::setCheckpoint(checkpoint);
    for (i=0; i<size(); i++) {
        at(i)->setCheckpoint(checkpoint);
    }
}

void IQTreeMix::startCheckpoint() {
    checkpoint->startStruct("IQTreeMix" + convertIntToString(size()));
}

void IQTreeMix::saveCheckpoint() {
    size_t i;
    startCheckpoint();
    ASSERT(weights.size() == size());
    double* relative_weights = new double[size()];
    for (i=0; i<size(); i++) {
        relative_weights[i]=weights[i];
    }
    CKP_ARRAY_SAVE(size(), relative_weights);
    for (i=0; i<size(); i++) {
        checkpoint->startStruct("Tree" + convertIntToString(i+1));
        at(i)->saveCheckpoint();
        checkpoint->endStruct();
    }
    endCheckpoint();
    delete[] relative_weights;
}

void IQTreeMix::restoreCheckpoint() {
    size_t i;
    startCheckpoint();
    ASSERT(weights.size() == size());
    double* relative_weights = new double[size()];
    if (CKP_ARRAY_RESTORE(size(), relative_weights)) {
        for (i = 0; i < size(); i++) {
            this->weights[i] = relative_weights[i];
            this->weight_logs[i] = log(relative_weights[i]);
        }
    }
    for (i=0; i<size(); i++) {
        checkpoint->startStruct("Tree" + convertIntToString(i+1));
        at(i)->restoreCheckpoint();
        checkpoint->endStruct();
    }
    endCheckpoint();
    clearAllPartialLH();
    delete[] relative_weights;
}

void IQTreeMix::setMinBranchLen(Params& params) {
    size_t i;
    int num_prec;
    
    if (params.min_branch_length <= 0.0) {
        params.min_branch_length = MAST_MIN_BRANCH_LEN;
        /*
        if (size() > 0) {
            if (!at(0)->isSuperTree() && at(0)->getAlnNSite() >= 100000) {
                params.min_branch_length = MAST_MIN_BRANCH_LEN; // 0.1 / (at(0)->getAlnNSite());
                // params.min_branch_length = 0.1 / (at(0)->getAlnNSite());
                num_prec = max((int)ceil(-log10(Params::getInstance().min_branch_length))+1, 6);
                for (i=0; i<size(); i++)
                    at(i)->num_precision = num_prec;
                cout.precision(12);
                // cout << "NOTE: minimal branch length is reduced to " << params.min_branch_length << " for long alignment" << endl;
                cout.precision(3);
            }
        }
        // Increase the minimum branch length if PoMo is used.
        if (aln->seq_type == SEQ_POMO) {
            params.min_branch_length *= aln->virtual_pop_size * aln->virtual_pop_size;
            cout.precision(12);
            cout << "NOTE: minimal branch length is increased to " << params.min_branch_length << " because PoMo infers number of mutations and frequency shifts" << endl;
            cout.precision(3);
        }
        */
    }
    cout << setprecision(7) << "Minimum branch length is set to " << params.min_branch_length << endl;
}

/** set pointer of params variable */
void IQTreeMix::setParams(Params* params) {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->setParams(params);
    }
    this->params = params;
};

/*
 * Generate the branch IDs
 * Branches of different trees with the same partition share the same ID
 */
void IQTreeMix::computeBranchID() {
    StrVector taxname;
    NodeVector nodes1, nodes2;
    Split* split;
    string split_str0, split_str1, split_str;
    size_t i,j,l;
    size_t branch_idx;
    UINT k;
    int curr_id;
    
    map<string,int> split_id;
    map<string,int>::iterator itr;
    branch_id.resize(nbranch * size());
    DoubleVector lenvec;
    IntVector ids;
    
    // convert the trees to unrooted if they are rooted
    for (i=0; i<size(); i++) {
        if (at(i)->rooted) {
            at(i)->convertToUnrooted();
            if (verbose_mode >= VB_MED)
                cout << "Convert tree " << i+1 << " into an unrooted tree" << endl;
        }
    }
    
    for (j=0; j<size(); j++) {
        
        nodes1.clear();
        nodes2.clear();
        taxname.clear();
        lenvec.clear();
        ids.clear();

        at(j)->getTaxaName(taxname);
        at(j)->getBranches(nodes1,nodes2,ids);
        at(j)->buildNodeSplit();
        for (i=0; i<nodes1.size(); i++) {
            branch_idx = j * nbranch + ids[i];
            split = at(j)->getSplit(nodes2[i], nodes1[i]);
            split_str0 = "";
            split_str1 = "";
            for (l=0; l<split->size(); l++) {
                for (k = 0; k < UINT_BITS && (l*UINT_BITS+k < split->getNTaxa()); k++) {
                    if (split->at(l) & (1 << k)) {
                        if (split_str1.length() > 0) {
                            split_str1.append(",");
                        }
                        split_str1.append(taxname[l * UINT_BITS + k]);
                    } else {
                        if (split_str0.length() > 0) {
                            split_str0.append(",");
                        }
                        split_str0.append(taxname[l * UINT_BITS + k]);
                    }
                }
            }
            // select the smaller one according to lexical order
            if (split_str0 < split_str1)
                split_str = split_str0;
            else
                split_str = split_str1;

            itr = split_id.find(split_str);
            if (itr == split_id.end()) {
                // new split_str
                curr_id = split_id.size();
                split_id.insert(pair<string,int> (split_str, curr_id) );
                // cout << split_str << " -> " << curr_id << endl;
                branch_group.push_back(vector<int>(1,branch_idx));
            } else {
                curr_id = itr->second;
                branch_group[curr_id].push_back(branch_idx);
            }
            branch_id[branch_idx] = curr_id;
            // cout << "branch_idx:" << branch_idx << " -> " << curr_id << endl;
        }
    }
    
    /*
    // for debugging
    // list branch_id
    cout << "branch_id:";
    for (i=0; i<branch_id.size(); i++) {
        if (i>0)
            cout << ",";
        cout << branch_id[i];
    }
    cout << endl;
    
    // list branch_group
    cout << "branch_group" << endl;
    for (i=0; i<branch_group.size(); i++) {
        cout << "group " << i << ":";
        for (j=0; j<branch_group[i].size(); j++) {
            if (j>0)
                cout << ",";
            cout << branch_group[i].at(j);
        }
        cout << endl;
    }
    */
}

/**
 * Generate the initial tree (usually used for model parameter estimation)
 */
void IQTreeMix::computeInitialTree(LikelihoodKernel kernel, istream* in) {
    size_t i;
    ifstream fin;

    if (size() == 0)
        outError("No tree is inputted for the tree-mixture model");
    if (params->user_file == NULL) {
        outError("Tree file has to be inputed (using the option -te) for tree-mixture model");
    }
    
    fin.open(params->user_file);
    
    for (i=0; i<size(); i++) {
        at(i)->computeInitialTree(kernel, &fin);
    }
    
    fin.close();
    
    // show trees
    // showTree();
}

/**
 * setup all necessary parameters
 */
void IQTreeMix::initSettings(Params &params) {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->initSettings(params);
    }
    IQTree::initSettings(params);
}

uint64_t IQTreeMix::getMemoryRequired(size_t ncategory, bool full_mem) {
    uint64_t mem_size = 0;
    size_t i;
    for (i=0; i<size(); i++) {
        mem_size += at(i)->getMemoryRequired(ncategory, full_mem);
    }
    return mem_size;
}

// get memory requirement for ModelFinder
uint64_t IQTreeMix::getMemoryRequiredThreaded(size_t ncategory, bool full_mem) {
    // only get the largest k partitions (k=#threads)
    int threads = (params->num_threads != 0) ? params->num_threads : params->num_threads_max;
    threads = min(threads, countPhysicalCPUCores());
    threads = min(threads, (int)size());
    
    // sort partition by computational cost for OpenMP effciency
    uint64_t *part_mem = new uint64_t[size()];
    int i;
    for (i = 0; i < size(); i++) {
        part_mem[i] = at(i)->getMemoryRequired(ncategory, full_mem);
    }
    
    // sort partition memory in increasing order
    quicksort<uint64_t, int>(part_mem, 0, size()-1);
    
    uint64_t mem = 0;
    for (i = size()-threads; i < size(); i++)
        mem += part_mem[i];
    
    delete [] part_mem;
    
    return mem;
}

void IQTreeMix::setNumThreads(int num_threads) {

    PhyloTree::setNumThreads(num_threads);

    if (num_threads >= size()) {
        if (num_threads % size() != 0)
            cout << "Warnings: setting number of threads as the multiples of the number of trees will be more efficient!";
        // distribute threads among the trees
        int threads_per_tree = num_threads / size();
        int* nthreads = new int[size()];
        int i;
        for (i = 0; i < size(); i++)
            nthreads[i] = threads_per_tree;
        i = 0;
        num_threads -= (size() * threads_per_tree);
        while (num_threads > 0) {
            nthreads[i]++;
            num_threads--;
            i++;
        }
        cout << "Number of threads for trees:";
        for (i = 0; i < size(); i++) {
            cout << " " << nthreads[i];
            at(i)->setNumThreads(nthreads[i]);
        }
        cout << endl;
        delete[] nthreads;
        isNestedOpenmp = true;
    } else {
        for (size_t i = 0; i < size(); i++)
            at(i)->setNumThreads(num_threads);
    }
}

/**
    test the best number of threads
*/
int IQTreeMix::testNumThreads() {
    // int bestNThres = at(0)->testNumThreads();
    int bestNThres = size();
    setNumThreads(bestNThres);
    return bestNThres;
}

/**
    optimize tree k separately
 */
void IQTreeMix::optimizeTreeSeparately(int k, bool printInfo, double logl_epsilon, double gradient_epsilon) {
    size_t step;
    double prev_score, score, t_score;
    PhyloTree* ptree;
    int n = 1;
    
    // save the trees
    ptree = at(k)->getRate()->getTree();
    at(k)->getRate()->setTree(at(k));
    if (anySiteRate)
        at(k)->getModelFactory()->site_rate->setTree(at(k));
    at(k)->clearAllPartialLH();
    prev_score = score = at(k)->computeLikelihood();

    for (step = 0; step < optimize_steps; step++) {
        
        // optimize the unlinked substitution models one by one
        if (!isLinkModel) {
             models[k]->optimizeParameters(gradient_epsilon);
        }

        // optimize the tree branches
        if (params->fixed_branch_length != BRLEN_FIX)
            score = at(k)->optimizeAllBranches(n, logl_epsilon);
        else
            score = at(k)->computeLikelihood();

        // for unlinked substitution models, optimize the unlinked site-rate models one by one
        if (anySiteRate && !isLinkSiteRate && !isLinkModel) {
            site_rates[k]->optimizeParameters(gradient_epsilon);
        }

        if (score < prev_score + gradient_epsilon) {
            // converged
            break;
        }
        prev_score = score;
    }
    
    // restore the trees
    at(k)->getRate()->setTree(ptree);
    if (anySiteRate)
        at(k)->getModelFactory()->site_rate->setTree(site_rate_trees[k]);
}


/**
    optimize each tree separately
    prerequisite: the array weights should be set
 */
void IQTreeMix::optimizeTreesSeparately(bool printInfo, double logl_epsilon, double gradient_epsilon) {

    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(1);
        omp_set_max_active_levels(2);
    }
    #endif
    #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
    for (size_t i=0; i<ntree; i++) {
        // optimize tree i
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            omp_set_num_threads(at(i)->num_threads);
        }
        #endif
        optimizeTreeSeparately(i, printInfo, logl_epsilon, gradient_epsilon);
    }
    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(num_threads);
    }
    #endif
}

void showDoubleArrayContent(string name, int dim, double* arr) {
    int i;
    // show the values of array
    cout << name << ":";
    for (i=1; i<=dim; i++)
        cout << " " << arr[i];
    cout << endl;
}

/**
    Initialize the tree weights using parsimony scores
    Idea:
    1. Check the parsimony score for each tree along all the sites
    2. Select the sites with different parsimony score between the trees.
    3. For each selected site, we check which parsimony score of the tree is minimum, and assign the site to the tree.
    4. The tree weights are estimated according to the proportion of the sites assigned to each tree.
 */
void IQTreeMix::initializeTreeWeights() {
    size_t i, j, noptn;
    int* parsimony_scores;
    UINT min_par_score;
    vector<int> tree_with_min_pars;
    double weight_sum;

    if (!parsi_computed) {
        // compute parsimony scores for each tree along the patterns
        // results are stored in the array patn_parsimony
        computeParsimony();
    }

    // reset the tree weights
    for (i=0; i<ntree; i++) {
        weights[i] = 0.0;
    }
    
    // estimate the tree weights
    for (i=0; i<nptn; i++) {
        if (!aln->at(i).isInformative())
            continue;
        parsimony_scores = &(patn_parsimony[i * ntree]);
        min_par_score = parsimony_scores[0];
        if (min_par_score < 0) // not a parsimony informative site
            continue;
        tree_with_min_pars.clear();
        tree_with_min_pars.push_back(0);
        for (j=1; j<ntree; j++) {
            if (parsimony_scores[j] < min_par_score) {
                tree_with_min_pars.clear();
                tree_with_min_pars.push_back(j);
                min_par_score = parsimony_scores[j];
            } else if (parsimony_scores[j] == min_par_score) {
                tree_with_min_pars.push_back(j);
            }
        }
        // if (tree_with_min_pars.size() == 1) {
        if (tree_with_min_pars.size() < ntree) {
            for (j=0; j<tree_with_min_pars.size(); j++) {
                weights[tree_with_min_pars[j]] += ptn_freq[i] / (double) tree_with_min_pars.size();
            }
        }
    }
    
    // normalize the tree weights
    weight_sum = 0.0;
    for (i=0; i<ntree; i++) {
        weight_sum += weights[i];
    }
    if (weight_sum > 0.0) {
        // there are parsimony informative sites
        for (i=0; i<ntree; i++) {
            weights[i] = weights[i] / weight_sum;
        }
    } else {
        // there is no parsimony informative sites
        for (i=0; i<ntree; i++) {
            weights[i] = 1.0/ntree;
        }
    }

    // If there are multiple tree weights belonging to the same group
    // set all the tree weights of the same group to their average
    // when weightGrpExist
    checkWeightGrp();

    // update the weight_logs
    for (i = 0; i < ntree; i++) {
        weight_logs[i] = log(weights[i]);
    }

    // show the initial tree weights
    cout << "According to the parsimony scores along the sites, the tree weights are initialized to:" << endl;

    for (i=0; i<ntree; i++) {
        if (i>0)
            cout << ",";
        cout << weights[i];
    }
    cout << endl;
}

/**
    Initialize the tree weights using parsimony scores
    Idea:
    1. Check the parsimony score for each tree along all the parsimony informative sites
    2. Initialize equal tree weights
    3. For each sites, compute the posterior probabilities of the trees with the same minimum parsimony scores
    4. Update the tree weights according to the posterior probabilities of trees across all the sites
    5. Repeat the steps 3 and 4 until the sum of the absolute changes in the tree weights <= 0.0001
    6. Report the final tree weights
 */
void IQTreeMix::initializeTreeWeights2() {
    size_t i, j, k;
    int* parsimony_scores;
    int min_par_score;
    vector<vector<int> > tree_with_min_pars;
    vector<double> pre_weights;
    vector<int> freqs;
    double weight_sum;
    double diff_sum;

    if (!parsi_computed) {
        // compute parsimony scores for each tree along the patterns
        // results are stored in the array patn_parsimony
        computeParsimony();
    }
    
    // cout << "the parimony scores:" << endl;
    // collect the tree(s) with min parsimony scores along sites
    tree_with_min_pars.clear();
    freqs.clear();
    for (i=0; i<nptn; i++) {
        if (!aln->at(i).isInformative())
            continue;
        parsimony_scores = &(patn_parsimony[i * ntree]);
        min_par_score = parsimony_scores[0];
        if (min_par_score < 0) // not a parsimony informative site
            continue;
        // cout << parsimony_scores[0];
        vector<int> min_pars;
        min_pars.push_back(0);
        for (j=1; j<ntree; j++) {
            // cout << ", " << parsimony_scores[j];
            if (parsimony_scores[j] < min_par_score) {
                min_pars.clear();
                min_pars.push_back(j);
                min_par_score = parsimony_scores[j];
            } else if (parsimony_scores[j] == min_par_score) {
                min_pars.push_back(j);
            }
        }
        // cout << ", " << ptn_freq[i] << endl;
            tree_with_min_pars.push_back(min_pars);
            freqs.push_back(ptn_freq[i]);
    }
    // cout << "--------------------" << endl;
    // exit(1);
    
    // show the parsimony statistics
    if (verbose_mode >= VB_MED) {
        map<vector<int>,int> pars_freqs;
        map<vector<int>,int>::iterator itr;
        for (i=0; i<tree_with_min_pars.size(); i++) {
            itr = pars_freqs.find(tree_with_min_pars[i]);
            if (itr == pars_freqs.end()) {
                // not found
                pars_freqs.insert(pair<vector<int>,int>(tree_with_min_pars[i],freqs[i]));
            } else {
                itr->second += freqs[i];
            }
        }
        // print out the summary
        cout << endl;
        cout << "Statistics on the sites with maximum parsimony score:" << endl;
        cout << "-----------------------------------------------------" << endl;
        cout << "Topologies; Number of parsimony informatic sites" << endl;
        for (itr = pars_freqs.begin(); itr != pars_freqs.end(); itr++) {
            vector<int> pars = itr->first;
            for (i=0; i<pars.size(); i++) {
                if (i>0)
                    cout << ",";
                cout << pars[i] + 1;
            }
            cout << "; " << itr->second << endl;
        }
        cout << endl;
    }

    // initialize the tree weights
    double w = 1.0 / ntree;
    double w_log = log(w);
    for (i=0; i<ntree; i++) {
        weights[i] = w;
        weight_logs[i] = w_log;
        pre_weights.push_back(weights[i]);
    }
    
    if (tree_with_min_pars.size() == 0)
        return;

    // estimate the tree weights
    diff_sum = ntree; // initialize diff_sum as a large value
    while (diff_sum > 0.0001) {
        // set pre_weights equal to weights
        for (i=0; i<ntree; i++)
            pre_weights[i] = weights[i];
        // reset weights
        for (i=0; i<ntree; i++)
            weights[i] = 0.0;
        for (i=0; i<tree_with_min_pars.size(); i++) {
            weight_sum = 0.0;
            if (tree_with_min_pars[i].size() == 1) {
                k = tree_with_min_pars[i].at(0);
                weights[k] += freqs[i];
            } else {
                for (j=0; j<tree_with_min_pars[i].size(); j++) {
                    k = tree_with_min_pars[i].at(j);
                    weight_sum += pre_weights[k];
                }
                for (j=0; j<tree_with_min_pars[i].size(); j++) {
                    k = tree_with_min_pars[i].at(j);
                    weights[k] += pre_weights[k] * freqs[i] / weight_sum; // posterior prob
                }
            }
        }
        // normalize weights
        weight_sum = 0.0;
        for (i=0; i<ntree; i++)
            weight_sum += weights[i];
        for (i=0; i<ntree; i++)
            weights[i] = weights[i] / weight_sum;
        /*
        // list out the intermediate results
        cout << "weights:";
        for (i=0; i<ntree; i++)
            cout << " " << weights[i];
        cout << endl;
         */
        // compute the sum of difference
        diff_sum = 0.0;
        for (i=0; i<ntree; i++)
            diff_sum += fabs(weights[i] - pre_weights[i]);
    }

    // If there are multiple tree weights belonging to the same group
    // set all the tree weights of the same group to their average
    // for weightGrpExist
    checkWeightGrp();

    // update the weight_logs
    for (i = 0; i < ntree; i++) {
        weight_logs[i] = log(weights[i]);
    }

    // show the initial tree weights
    cout << "According to the parsimony scores along the sites, the tree weights are initialized to:" << endl;

    for (i=0; i<ntree; i++) {
        if (i>0)
            cout << ",";
        cout << weights[i];
    }
    cout << endl;
}

// reset the ptn_freq array to the original frequencies of the patterns
void IQTreeMix::resetPtnOrigFreq() {
    size_t i, ptn;
    for (i = 0; i < ntree; i++) {
        memcpy(at(i)->ptn_freq, ptn_freq, sizeof(double)*nptn);
        /*
        for (ptn = 0; ptn < nptn; ptn++) {
            at(i)->ptn_freq[ptn] = ptn_freq[ptn];
        }*/
    }
}

string IQTreeMix::optimizeModelParameters(bool printInfo, double logl_epsilon) {
    size_t i, ptn;
    int step, n, m, substep1, nsubstep1, nsubstep1_start, nsubstep1_max, nsubstep2_start, nsubstep2_max, substep2, nsubstep2, substep2_tot;
    double* pattern_mix_lh;
    double curr_epsilon;
    double prev_score, prev_score1, prev_score2, score, t_score;
    double* prev_ptn_invar;
    bool tree_weight_converge = false;
    bool firsttime_substmodel = true;
    bool firsttime_RHASmodel = true;
    bool firsttime_branchlen = true;
    double gradient_epsilon = 0.0001;
    PhyloTree *ptree;
    
    n = 1;
    nsubstep1_start = 5;
    nsubstep1_max = 10;
    nsubstep2_start = 5;
    nsubstep2_max = 10;
    substep2 = 0;
    
    // allocate memory
    pattern_mix_lh = new double[ntree * nptn];
    
    if (!isTreeWeightFixed) {
        // initialize the tree weights according to parsimony scores along the sites
        initializeTreeWeights();
    }
    
    // initialize the parameters
    optimizeTreesSeparately(printInfo, logl_epsilon, gradient_epsilon);
    prev_ptn_invar = ptn_invar;
    ptn_invar = at(0)->ptn_invar;

    if (params->fixed_branch_length != BRLEN_FIX) {
        // set all the branches of the same group to their weighted average for initialization of the branch lengths
        checkBranchGrp();
    }

    // show trees
    // cout << "Initial trees:" << endl;
    // showTree();
    // cout << setprecision(5) << "Estimate model parameters (logl epsilon = " << logl_epsilon << ", gradient epsilon = " << gradient_epsilon << ")" << endl;
    cout << setprecision(5) << "Estimate model parameters (epsilon = " << logl_epsilon << ")" << endl;
    score = computeLikelihood();
    cout << "1. Initial log-likelihood: " << score;
    if (verbose_mode >= VB_MED) {
        cout << " weights=(";
        for (i=0; i<ntree; i++) {
            if (i>0)
                cout << ",";
            cout << weights[i];
        }
        cout << ")";
    }
    cout << endl;

    bool is_ptnfrq_posterior = false;

    nsubstep1 = nsubstep1_start;
    nsubstep2 = nsubstep2_start;
    prev_score = score;

    for (step = 0; step < optimize_steps; step++) {
        
        prev_score1 = score;
        substep2_tot = 0;

        for (substep1 = 0; substep1<nsubstep1; substep1++) {
            
            prev_score2 = score;

            for (substep2=0; substep2<nsubstep2; substep2++) {

                // optimize the linked subsitution model
                if (isLinkModel) {
                    // reset the ptn_freq array to the original frequencies of the patterns
                    if (is_ptnfrq_posterior) {
                        resetPtnOrigFreq();
                        is_ptnfrq_posterior = false;
                    }
                    models[0]->optimizeParameters(gradient_epsilon);
                    if (verbose_mode >= VB_MED) {
                        score = computeLikelihood();
                        cout << "after optimizing linked subsitution model, likelihood = " << score << endl;
                    }
                }

                // optimize the unlinked subsitution models one by one
                if (!isLinkModel) {
                    if (!is_ptnfrq_posterior) {
                        // cout << "before calling computeFreqArray()" << endl;
                        computeFreqArray(pattern_mix_lh, true);
                        // cout << "after calling computeFreqArray()" << endl;
                        is_ptnfrq_posterior = true;
                    }
                    if (substep2 == 0) {
                        // call computeFreqArray() after each optimizeParameters()
                        for (i=0; i<ntree; i++) {
                            models[i]->optimizeParameters(gradient_epsilon);
                            computeFreqArray(pattern_mix_lh, true, i);
                        }
                    } else {
                        // call computeFreqArray() once after all optimizeParameters()
                        #ifdef _OPENMP
                        if (isNestedOpenmp) {
                            // omp_set_nested(1);
                            omp_set_max_active_levels(2);
                        }
                        #endif
                        #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
                        for (int k=0; k<ntree; k++) {
                            #ifdef _OPENMP
                            if (isNestedOpenmp) {
                                omp_set_num_threads(at(k)->num_threads);
                            }
                            #endif
                            models[k]->optimizeParameters(gradient_epsilon);
                        }
                        #ifdef _OPENMP
                        if (isNestedOpenmp) {
                            // omp_set_nested(0);
                            omp_set_max_active_levels(1);
                            omp_set_num_threads(num_threads);
                        }
                        #endif
                        computeFreqArray(pattern_mix_lh, true);
                    }
                    if (verbose_mode >= VB_MED) {
                        score = computeLikelihood();
                        cout << "after optimizing unlinked subsitution model, likelihood = " << score << endl;
                    }
                }

                // optimize tree branches for non-branch-length-restricted model
                if (!isEdgeLenRestrict && params->fixed_branch_length != BRLEN_FIX) {
                    if (!is_ptnfrq_posterior) {
                        computeFreqArray(pattern_mix_lh, true);
                        is_ptnfrq_posterior = true;
                    }
                    if (substep2 == 0) {
                        // call computeFreqArray() after each optimizeAllBranchesOneTree()
                        for (i=0; i<size(); i++) {
                            optimizeAllBranchesOneTree(i, 1, logl_epsilon);
                            computeFreqArray(pattern_mix_lh, true, i);
                        }
                    } else {
                        // call computeFreqArray() once after optimizeAllBranches()
                        optimizeAllBranches(1, logl_epsilon);
                        computeFreqArray(pattern_mix_lh, true);
                    }
                    if (verbose_mode >= VB_MED) {
                        score = computeLikelihood();
                        cout << "after optimizing branches, likelihood = " << score << endl;
                    }
                }

                // optimize tree branches for branch-length-restricted model
                if (isEdgeLenRestrict && params->fixed_branch_length != BRLEN_FIX) {
                    if (is_ptnfrq_posterior) {
                        resetPtnOrigFreq();
                        is_ptnfrq_posterior = false;
                    }
                    optimizeBranchLensByBFGS(gradient_epsilon);
                    if (verbose_mode >= VB_MED) {
                        score = computeLikelihood();
                        cout << "after optimizing branches by BFGS, likelihood = " << score << endl;
                    }
                }

                score = computeLikelihood();
                if (score < prev_score2 + logl_epsilon) {
                    // converged
                    break;
                }
                prev_score2 = score;
            }

            substep2_tot += substep2;

            // optimize the linked site rate model
            if (anySiteRate && isLinkSiteRate) {
                if (is_ptnfrq_posterior) {
                    resetPtnOrigFreq();
                    is_ptnfrq_posterior = false;
                }
                site_rates[0]->setTree(this);
                site_rates[0]->optimizeParameters(gradient_epsilon);
                if (verbose_mode >= VB_MED) {
                    score = computeLikelihood();
                    cout << "after optimizing linked site rate model, likelihood = " << score << endl;
                }
            }

            // optimize the unlinked site rate model
            if (anySiteRate && !isLinkSiteRate) {
                if (!is_ptnfrq_posterior) {
                    computeFreqArray(pattern_mix_lh, true);
                    is_ptnfrq_posterior = true;
                }
                // if (substep1 == 0) {
                    // call computeFreqArray() after each optimizeParameters()
                    for (i=0; i<ntree; i++) {
                        site_rates[i]->optimizeParameters(gradient_epsilon);
                        computeFreqArray(pattern_mix_lh, true, i);
                    }
                /* } else {
                    // call computeFreqArray() once after all optimizeParameters()
                    if (isNestedOpenmp)
                        omp_set_nested(1);
                    #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
                    for (int k=0; k<ntree; k++) {
                        if (isNestedOpenmp) {
                            omp_set_num_threads(at(k)->num_threads);
                        }
                        site_rates[k]->optimizeParameters(gradient_epsilon);
                    }
                    if (isNestedOpenmp) {
                        omp_set_nested(0);
                        omp_set_num_threads(num_threads);
                    }
                    computeFreqArray(pattern_mix_lh, true);
                }*/
                if (verbose_mode >= VB_MED) {
                    score = computeLikelihood();
                    cout << "after optimizing unlinked site-rate model, likelihood = " << score << endl;
                }
            }

            score = computeLikelihood();

            if (score < prev_score1 + logl_epsilon) {
                // converged
                break;
            }
            prev_score1 = score;
        }
        
        // optimize tree weights
        if (!isTreeWeightFixed) {
            // if (weightGrpExist || params->optimize_alg_treeweight == "BFGS") {
            if (params->optimize_alg_treeweight == "BFGS") {
                score = optimizeTreeWeightsByBFGS(gradient_epsilon);
                tree_weight_converge = true;
            } else {
                m = 1 + step / 100;
                score = optimizeTreeWeightsByEM(pattern_mix_lh, gradient_epsilon, m, tree_weight_converge);  // loop max n times
            }
            if (verbose_mode >= VB_MED) {
                score = computeLikelihood();
                cout << "after optimizing tree weights, likelihood = " << score << endl;
            }
        }

        cout << step+2 << ". Current log-likelihood: " << score;
        if (verbose_mode >= VB_MED) {
            cout << " weights=(";
            for (i=0; i<ntree; i++) {
                if (i>0)
                    cout << ",";
                cout << weights[i];
            }
            cout << ")";
        }
        cout << endl;

        if (score < prev_score + logl_epsilon) {
            // converged
            break;
        }
        
        if (nsubstep1 < nsubstep1_max)
            nsubstep1++;
        if (nsubstep2 < nsubstep2_max)
            nsubstep2++;

        prev_score = score;
    }

    setCurScore(score);
    stop_rule.setCurIt(step);
    
    // compute the proportion of sites for each tree with the maximum posterior probability
    // computeMaxPosteriorRatio(pattern_mix_lh, false, false);
    
    // compute the proportion of sites for each tree with the maximum high-enough likelihood values
    // computeMaxLikeRatio();

    delete[] pattern_mix_lh;
    ptn_invar = prev_ptn_invar;
    
    // show the weights
    cout << "Final estimation on weights:" << endl;
    for (i = 0; i < ntree; i++) {
        if (i > 0)
            cout << ",";
        cout << weights[i];
    }
    cout << endl;
    
    return getTreeString();
}

/**
        print tree to .treefile
        @param params program parameters, field root is taken
 */
void IQTreeMix::printResultTree(string suffix) {
    ofstream fout;
    size_t i;
    
    if (MPIHelper::getInstance().isWorker()) {
        return;
    }
    if (params->suppress_output_flags & OUT_TREEFILE)
        return;
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    if (suffix.compare("") != 0) {
        tree_file_name += "." + suffix;
    }
    fout.open(tree_file_name.c_str());
    setRootNode(params->root, true);
    for (i=0; i<size(); i++) {
        at(i)->printTree(fout);
        fout << endl;
    }
    setRootNode(params->root, false);
    fout.close();
    if (verbose_mode >= VB_MED)
        cout << "Best tree printed to " << tree_file_name << endl;
}

string IQTreeMix::getTreeString() {
    stringstream tree_stream;
    size_t i;
    
    for (i=0; i<size(); i++) {
        at(i)->printTree(tree_stream, WT_TAXON_ID + WT_BR_LEN + WT_SORT_TAXA);
        tree_stream << "\n";
    }
    return tree_stream.str();
}

// return the weighted sum of the tree lengths over all trees
double IQTreeMix::treeLength(Node *node, Node *dad) {
    double len = 0.0;
    size_t i;
    for (i=0; i<size(); i++)
        len += at(i)->treeLength() * weights[i];
    return len;
}

// return the weighted sum of the lengths of all internal branches over all trees
double IQTreeMix::treeLengthInternal( double epsilon, Node *node, Node *dad) {
    double len = 0.0;
    size_t i;
    for (i = 0; i < size(); i++)
        len += at(i)->treeLengthInternal(epsilon) * weights[i];
    return len;
}

int IQTreeMix::getNParameters() {
    int df = 0;
    int k;
    size_t i;
    
    if (verbose_mode >= VB_MED)
        cout << endl << "Number of parameters:" << endl;
    for (i=0; i<models.size(); i++) {
        k = models[i]->getNDim() + models[i]->getNDimFreq();
        if (verbose_mode >= VB_MED) {
            if (models.size() == 1)
                cout << " linked model : " << k << endl;
            else
                cout << " model " << i+1 << " : " << k << endl;
        }
        df += k;
    }
    for (i=0; i<site_rates.size(); i++) {
        if (verbose_mode >= VB_MED) {
            if (site_rates.size() == 1)
                cout << " linked site rate : " << site_rates[i]->getNDim() << endl;
            else
                cout << " siterate " << i+1 << " : " << site_rates[i]->getNDim() << endl;
        }
        df += site_rates[i]->getNDim();
    }
    // for branch parameters
    if (params->fixed_branch_length != BRLEN_FIX) {
        if (isEdgeLenRestrict) {
            if (verbose_mode >= VB_MED)
                cout << " branch groups (for branch-len-restricted) : " << branch_group.size() << endl;
            df += branch_group.size();
        } else {
            for (i=0; i<size(); i++) {
                k = at(i)->getNBranchParameters(BRLEN_OPTIMIZE);
                if (verbose_mode >= VB_MED)
                    cout << " branches of tree " << i+1 << " : " << k << endl;
                df += k;
            }
        }
    }
    // for tree weights
    if (!isTreeWeightFixed) {
        if (weightGrpExist) {
            if (verbose_mode >= VB_MED)
                cout << " tree weight groups (for weight-restricted) : " << (weight_group_member.size() - 1) << endl;
            df += (weight_group_member.size() - 1);
        } else {
            if (verbose_mode >= VB_MED)
                cout << " tree weights : " << (size() - 1) << endl;
            df += (size() - 1);
        }
    }
    if (verbose_mode >= VB_MED)
        cout << " == Total : " << df << " == " << endl << endl;
    return df;
}

void IQTreeMix::drawTree(ostream &out, int brtype, double zero_epsilon) {
    size_t i;
    for (i=0; i<size(); i++) {
        out << "Tree " << i+1 << ":" << endl;
        at(i)->drawTree(out, brtype, zero_epsilon);
    }
}

/**
        print the tree to the output file in newick format
        @param out the output file.
        @param node the starting node, NULL to start from the root
        @param dad dad of the node, used to direct the search
        @param brtype type of branch to print
        @return ID of the taxon with smallest ID
 */
void IQTreeMix::printTree(ostream & out, int brtype) {
    for (int i=0; i<size(); i++) {
        at(i)->printTree(out, brtype);
    }
}
    
/*
int IQTreeMix::printTree(ostream &out, int brtype, Node *node, Node *dad) {
    size_t i;
    int value = 0;
    for (i=0; i<size(); i++) {
        // out << "Tree " << i+1 << ":" << endl;
        value = at(i)->printTree(out, brtype, node, dad);
        if (i<size()-1)
            out << ";" << endl;
    }
    return value;
}
*/
    
/**
 *  Return best tree string from the candidate set
 *
 *  @param numTrees
 *      Number of best trees to return
 *  @return
 *      A string vector of trees
 */
vector<string> IQTreeMix::getBestTrees(int numTrees) {
    vector<string> res;
    size_t i,j;
    for (i=0; i<size() && numTrees>0; i++) {
        vector<string> curr_res = at(i)->getBestTrees(numTrees);
        if (curr_res.size() < numTrees) {
            numTrees = curr_res.size();
        }
        if (numTrees > 0) {
            if (i==0) {
                // the array res is empty
                for (j=0; j<numTrees; j++) {
                    res.push_back(curr_res[j]);
                }
            } else {
                // the array res is not empty
                for (j=0; j<numTrees; j++) {
                    res[j].append(";"+curr_res[j]);
                }
            }
        }
    }
    if (numTrees == 0)
        res.clear();
    else if (numTrees < res.size())
        res.resize(numTrees);
    return res;
}

/**
        Read the tree saved with Taxon IDs and branch lengths.
        @param tree_string tree string to read from
        @param updatePLL if true, tree is read into PLL
 */
void IQTreeMix::readTreeString(const string &tree_string) {
    vector<string> substrs;
    size_t i;
    
    separateStr(tree_string, substrs, ';');
    ASSERT(substrs.size() == size());
    for (i=0; i<size(); i++) {
        at(i)->readTreeString(substrs[i]);
    }
}

/*
// get posterior probabilities along each site for each tree
void IQTreeMix::getPostProb(double* pattern_mix_lh, bool need_computeLike, int updated_which_tree, bool need_multiplyFreq) {
    size_t i, ptn, c;
    double* this_lk_cat;
    double lk_ptn;
    double ln_freq_minus_ln_like;
    double v;

    if (need_computeLike) {
        computeSiteTreeLogLike(updated_which_tree);
        // computeLikelihood_combine();
        // computeLikelihood();
    }

    memcpy(pattern_mix_lh, ptn_like_cat, nptn*ntree*sizeof(double));

    // pattern_mix_lh = pattern_mix_lh + logarithm of tree weights - log-likelihood of the pattern + logarithm of pattern frequency
    i = 0;
    for (ptn = 0; ptn < nptn; ptn++) {
        if (need_multiplyFreq)
            v = ptn_freq_logs[ptn] - ptn_like[ptn];
        else
            v = - ptn_like[ptn];
        for (c = 0; c < ntree; c++) {
            pattern_mix_lh[i] += (weight_logs[c] + v);
            i++;
        }
    }
    c = nptn * ntree;
    for (i = 0; i < c; i++)
        pattern_mix_lh[i] = exp(pattern_mix_lh[i]);
}
*/

// get posterior probabilities along each site for each tree
void IQTreeMix::getPostProb(double* pattern_mix_lh, bool need_computeLike, int update_which_tree, bool need_multiplyFreq) {
    // size_t i, ptn, c;
    // double* this_lk_cat;
    // double lk_ptn;

    if (need_computeLike) {
        computeSiteTreeLogLike(update_which_tree);
    }
    /*
    cout << "ptn_like_cat:";
    for (size_t i = 0; i < nptn*ntree; i++) {
        cout << " " << ptn_like_cat[i];
    }
    cout << endl;
    cout << "end of ptn_like_cat" << endl;
    */
    memcpy(pattern_mix_lh, ptn_like_cat, nptn*ntree*sizeof(double));

    /*
    cout << "2" << endl;
    cout << "ntree = " << ntree << endl;
    cout << "nptn = " << nptn << endl;
    cout << "weights:";
    for (size_t c = 0; c < ntree; c++) {
        cout << " " << weights[c];
    }
    cout << endl;
    cout << "ptn_freq:";
    for (size_t ptn = 0; ptn < nptn; ptn++) {
        cout << " " << ptn_freq[ptn];
    }
    cout << endl;
    cout << "pattern_mix_lh:";
    for (size_t i = 0; i < nptn*ntree; i++) {
        cout << " " << pattern_mix_lh[i];
    }
    cout << endl;
    cout << "here 0" << endl << flush;
    */

    // multiply pattern_mix_lh by tree weights
    #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
    for (size_t ptn = 0; ptn < nptn; ptn++) {
        size_t i = ptn * ntree;
        for (size_t c = 0; c < ntree; c++) {
            pattern_mix_lh[i] *= weights[c];
            i++;
        }
    }

    // cout << "here 1" << endl << flush;

    if (need_multiplyFreq) {
        #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
        for (size_t ptn = 0; ptn < nptn; ptn++) {
            double* this_lk_cat = pattern_mix_lh + ptn * ntree;
            double lk_ptn = 0.0;
            for (size_t c = 0; c < ntree; c++) {
                lk_ptn += this_lk_cat[c];
            }
            ASSERT(lk_ptn != 0.0);
            lk_ptn = ptn_freq[ptn] / lk_ptn;
            
            // transform pattern_mix_lh into posterior probabilities of each category
            for (size_t c = 0; c < ntree; c++) {
                this_lk_cat[c] *= lk_ptn;
            }
        }
    } else {
        #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
        for (size_t ptn = 0; ptn < nptn; ptn++) {
            double* this_lk_cat = pattern_mix_lh + ptn * ntree;
            double lk_ptn = 0.0;
            for (size_t c = 0; c < ntree; c++) {
                lk_ptn += this_lk_cat[c];
            }
            ASSERT(lk_ptn != 0.0);
            lk_ptn = 1.0 / lk_ptn;
            
            // transform pattern_mix_lh into posterior probabilities of each category
            for (size_t c = 0; c < ntree; c++) {
                this_lk_cat[c] *= lk_ptn;
            }
        }
    }
    // cout << "here 2" << endl << flush;
}

// update the ptn_freq array according to the posterior probabilities along each site for each tree
void IQTreeMix::computeFreqArray(double* pattern_mix_lh, bool need_computeLike, int update_which_tree) {

    // get posterior probabilities along each site for each tree
    getPostProb(pattern_mix_lh, need_computeLike, update_which_tree);

    // #pragma omp parallel for schedule(dynamic) num_threads(num_threads) if (num_threads > 1)
    for (size_t i = 0; i < ntree; i++) {
        PhyloTree* tree = at(i);
        double *this_lk_cat = pattern_mix_lh+i;
        for (size_t ptn = 0; ptn < nptn; ptn++) {
            tree->ptn_freq[ptn] = this_lk_cat[0];
            this_lk_cat += ntree;
        }
    }
}

// For the linked RHAS model
// Store the RHAS variables of tree 0 to the array rhas_var
void IQTreeMix::storeTree0RHAS() {
    RateFree* rfmodel = NULL;
    int rf_orig_optim_params;
    
    if (at(0)->getRate()->isFreeRate()) {
        // for Free Rate model
        rfmodel = dynamic_cast<RateFree*>(at(0)->getRate());
        rf_orig_optim_params = rfmodel->optimizing_params;
        // change the optimizing params to 0 so that the function getVariables() replaces all variables
        rfmodel->optimizing_params = 0;
    }

    if (rhas_var == NULL) {
        // allocate memory to the array rhas_var
        rhas_var = aligned_alloc<double>(at(0)->getRate()->getNDim() + 1);
    }
    // Store the RHAS variables for tree 0 to the array rhas_var
    at(0)->getRate()->setVariables(rhas_var);

    if (at(0)->getRate()->isFreeRate()) {
        // for Free Rate model
        // change the optimizing params back to the original value
        rfmodel->optimizing_params = rf_orig_optim_params;
    }
}

// For the linked RHAS model
// Replace the RHAS variables of tree t by those of tree 0
// The array rhas_var should have stored the updated RHAS variables of tree 0
void IQTreeMix::copyRHASfrTree0(int t) {
    RateFree* rfmodel = NULL;
    int rf_orig_optim_params;
    
    
    if (at(t)->getRate()->isFreeRate()) {
        // for Free Rate model
        rfmodel = dynamic_cast<RateFree*>(at(t)->getRate());
        rf_orig_optim_params = rfmodel->optimizing_params;
        // change the optimizing params to 0 so that the function getVariables() replaces all variables
        rfmodel->optimizing_params = 0;
    }

    // copy the RHAS variables of tree 0 to this
    at(t)->getRate()->getVariables(rhas_var);

    if (at(t)->getRate()->isFreeRate()) {
        // for Free Rate model
        // change the optimizing params back to the original value
        rfmodel->optimizing_params = rf_orig_optim_params;
    }
}

double IQTreeMix::targetFunk(double x[]) {
    double score;
    getVariables(x);
    if (optim_type==1)
        score = -computeLikelihood_combine();
    else
        score = -computeLikelihood();
    return score;
}

// read the tree weights and write into "variables"
void IQTreeMix::setVariables(double *variables) {
    // for tree weights
    size_t i;
    size_t ndim;
    size_t treeidx, branchidx;
    
    if (optim_type==1) {
        // optimization on tree weight
        ndim = weight_group_member.size();
        for (i=0; i<ndim; i++) {
            variables[i+1] = tmp_weights[i];
        }
    } else {
        // optimization on branch length
        ndim = branch_group.size();
        for (i=0; i<ndim; i++) {
            treeidx = branch_group[i].at(0) / nbranch;
            branchidx = branch_group[i].at(0) % nbranch;
            variables[i+1] = branch_len[treeidx].at(branchidx);
        }
    }
}

// read the "variables" and write into tree
void IQTreeMix::getVariables(double *variables) {
    size_t i,j;
    size_t ndim;
    size_t treeidx,branchidx;
    double sum;
    double w, w_log;
    
    if (optim_type==1) {
        // for tree weight
        ndim = weight_group_member.size();
        sum = 0.0;
        for (i=0; i<ndim; i++) {
            tmp_weights[i] = variables[i+1];
            sum += tmp_weights[i] * weight_group_member[i].size();
        }
        for (i=0; i<ndim; i++) {
            w = tmp_weights[i] / sum;
            w_log = log(w);
            for (j=0; j<weight_group_member[i].size(); j++) {
                weights[weight_group_member[i].at(j)] = w;
                weight_logs[weight_group_member[i].at(j)] = w_log;
            }
        }
    } else {
        // for branch length
        ndim = branch_group.size();
        for (i=0; i<ndim; i++) {
            for (j=0; j<branch_group[i].size(); j++) {
                treeidx = branch_group[i].at(j) / nbranch;
                branchidx = branch_group[i].at(j) % nbranch;
                branch_len[treeidx].at(branchidx) = variables[i+1];
            }
        }
        setBranchLengths(branch_len);
    }
}

// set the bounds
void IQTreeMix::setBounds(double *lower_bound, double *upper_bound, bool* bound_check) {
    size_t i;
    size_t ndim;
    
    if (optim_type == 1) {
        // optimization on tree weight
        ndim = weight_group_member.size();
        for (i=0; i<ndim; i++) {
            lower_bound[i+1] = MIN_PROP;
            upper_bound[i+1] = MAX_PROP;
            bound_check[i+1] = false;
        }
    } else {
        // optimization on branch length
        ndim = branch_group.size();
        for (i=0; i<ndim; i++) {
            lower_bound[i+1] = params->min_branch_length;
            upper_bound[i+1] = MAX_LEN;
            bound_check[i+1] = false;
        }
    }
}

// get the dimension of the variables
int IQTreeMix::getNDim() {
    size_t s;
    if (optim_type == 1) {
        // optimization on tree weight
        if (weight_group_member.size() > 0)
            s = weight_group_member.size();
        else
            s = size();
    } else {
        // optimization on branch length
        s = branch_group.size();
    }
    return s;
}

/**
        get the name of the model
 */
string IQTreeMix::getModelName() {
    size_t i;
    string model = "";

    if (isLinkModel) {
        // for linked model
        model += model_names[0];
        
        if (anySiteRate) {
            model += "+";
            if (isLinkSiteRate) {
                // for linked site rate
                model += siterate_names[0];
            } else {
                // for unlinked site rate
                model += "TMIX{";
                for (i=0; i<size(); i++) {
                    if (i>0)
                        model += ",";
                    model += siterate_names[i];
                }
                model += "}";
            }
        }
    } else {
        // for unlinked model
        model += "TMIX{";
        if (anySiteRate && !isLinkSiteRate) {
            // show both unlinked model and unlinked site rate
            for (i=0; i<size(); i++) {
                if (i>0)
                    model += ",";
                model += model_names[i] + "+" + siterate_names[i];
            }
        } else {
            // show unlinked model only
            for (i=0; i<size(); i++) {
                if (i>0)
                    model += ",";
                model += model_names[i];
            }
        }
        model += "}";
        if (anySiteRate && isLinkSiteRate) {
            // for linked site rate
            model += "+" + siterate_names[0];
        }
    }
    if (isEdgeLenRestrict)
        model += "+TR";
    else if (model.find("TMIX") == string::npos)
        model += "+T";
    return model;
}

/**
        get the name of the model
 */
string IQTreeMix::getModelNameParams(bool show_fixed_params) {
    size_t i;
    string model = "";

    if (isLinkModel) {
        // for linked model
        model += at(0)->getModel()->getNameParams(show_fixed_params);
        
        if (anySiteRate) {
            if (isLinkSiteRate) {
                // for linked site rate
                model += at(0)->getRate()->getNameParams();
            } else {
                // for unlinked site rate
                model += "+TMIX{";
                for (i=0; i<size(); i++) {
                    if (i>0)
                        model += ",";
                    string ratenameparam = at(i)->getRate()->getNameParams();
                    if (ratenameparam.length() > 0 && ratenameparam[0] == '+') {
                        ratenameparam = ratenameparam.substr(1);
                    }
                    model += ratenameparam;
                }
                model += "}";
            }
        }
    } else {
        // for unlinked model
        model += "TMIX{";
        if (anySiteRate && !isLinkSiteRate) {
            // show both unlinked model and unlinked site rate
            for (i=0; i<size(); i++) {
                if (i>0)
                    model += ",";
                model += at(i)->getModel()->getNameParams(show_fixed_params)
                            + at(i)->getRate()->getNameParams();
            }
        } else {
            // show unlinked model only
            for (i=0; i<size(); i++) {
                if (i>0)
                    model += ",";
                model += at(i)->getModel()->getNameParams(show_fixed_params);
            }
        }
        model += "}";
        if (anySiteRate && isLinkSiteRate) {
            // for linked site rate
            model += at(0)->getRate()->getNameParams();
        }
    }
    if (isEdgeLenRestrict)
        model += "+TR";
    else if (model.find("TMIX") == string::npos)
        model += "+T";
    return model;
}
    
struct classcomp {
  bool operator() (const Pattern& lhs, const Pattern& rhs) const {
      for (size_t i=0; i<lhs.size() && i<rhs.size(); i++) {
          if (lhs[i] != rhs[i]) {
              return (lhs[i] < rhs[i]);
          }
      }
      return false;
  }
};

// compute parsimony scores for each tree along the patterns
// results are stored in the array patn_parsimony
void IQTreeMix::computeParsimony() {
    map<Pattern,int,classcomp> opattern2id;
    size_t noptn;

    noptn = aln->ordered_pattern.size();
    UINT* ptn_scores = new UINT[ntree * noptn];
    
    deleteAllPartialLh();
    
    // compute the parsimony scores along patterns for each tree
    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(1);
        omp_set_max_active_levels(2);
    }
    #endif
    #pragma omp parallel for schedule(static) num_threads(ntree) if (isNestedOpenmp)
    for (size_t t=0; t<ntree; t++) {
        #ifdef _OPENMP
        if (isNestedOpenmp) {
            omp_set_num_threads(at(t)->num_threads);
        }
        #endif
        UINT* curr_ptn_scores = ptn_scores + t * noptn;
        at(t)->initCostMatrix(CM_UNIFORM);
        at(t)->setParsimonyKernel(params->SSE);
        at(t)->initializeAllPartialPars();
        at(t)->computeTipPartialParsimony();
        at(t)->computeParsimonyOutOfTreeSankoff(curr_ptn_scores);
    }
    #ifdef _OPENMP
    if (isNestedOpenmp) {
        // omp_set_nested(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(num_threads);
    }
    #endif

    // build opattern2id
    for (size_t optn=0; optn<noptn; optn++) {
        if (aln->ordered_pattern[optn].frequency == 0)
            continue;
        opattern2id.insert(pair<Pattern,int>(aln->ordered_pattern.at(optn),optn));
        
        /*
        cout << optn+1;
        cout << "," << aln->ordered_pattern[optn].frequency;

        // show the parsimony scores
        curr_ptn_scores = ptn_scores;
        for (t=0; t<ntree; t++) {
            cout << "," << curr_ptn_scores[optn];
            curr_ptn_scores += noptn;
        }

        // show nucleotide (for ordered pattern)
        for (t=0; t<ntip; t++) {
            cout << "," << aln->convertStateBackStr(aln->ordered_pattern[optn].at(t));
        }

        cout << endl;
         */
    }

    #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
    for (size_t ptn=0; ptn<nptn; ptn++) {
        map<Pattern,int>::iterator itr = opattern2id.find(aln->at(ptn));
        int* curr_pars_scores = &patn_parsimony[ptn*ntree];
        if (itr != opattern2id.end()) {
            size_t optn = itr->second;
            // the parsimony scores
            UINT* curr_ptn_scores = ptn_scores;
            for (size_t t=0; t<ntree; t++) {
                curr_pars_scores[t] = curr_ptn_scores[optn];
                curr_ptn_scores += noptn;
            }

            /*
            cout << ptn+1;
            cout << "," << aln->at(ptn).frequency;

            // show the parsimony scores
            for (t=0; t<ntree; t++) {
                cout << "," << curr_pars_scores[t];
            }

            // show nucleotide (for ordered pattern)
            for (t=0; t<ntip; t++) {
                cout << "," << aln->convertStateBackStr(aln->at(ptn).at(t));
            }

            cout << endl;
            */
        } else {
            for (size_t t=0; t<ntree; t++) {
                curr_pars_scores[t] = -1; // parsimony score is not available
            }
        }
    }
    
    parsi_computed = true;

    // free the memory of the array
    delete[] ptn_scores;
}
