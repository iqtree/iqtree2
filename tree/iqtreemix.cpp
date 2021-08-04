//
//  iqtreemix.cpp
//  tree
//
//  Created by Thomas Wong on 14/12/20.
//

#include "iqtreemix.h"
const double MIN_PROP = 0.001;
const double MAX_PROP = 1000.0;

// Input formats for the tree-mixture model
// 1. linked models and site rates: GTR+G4+T
// 2. unlinked models and linked site rates: MIX{GTR,GTR}+G4+T
// 3. linked models and unlinked site rates: GTR+MIX{G4,E}+T
// 4. unlinked models and unlinked site rates: MIX{GTR+G4,GTR}+T
// The situation that a part of the model is linked while another part is unlinked is not allowed.
//    For example, MIX{GTR,GTR}+FO+T or GTR+MIX{FO+F0}+T is not be accepted
// Similarly, the situation that a part of the site rate is linked while another part is unlinked is also not allowed.
//    For example, GTR+MIX{I,I}+G4+T or GTR+I+MIX{G4+G4}+T is not be accepted
// The tree weights can be fixed by: T{0.2,0.3,0.5}, otherwise, optimization on tree weights will be performed.

IQTreeMix::IQTreeMix() : IQTree() {
    patn_freqs = NULL;
    patn_isconst = NULL;
    ptn_like_cat = NULL;
    _ptn_like_cat = NULL;
    weightGrpExist = false;
}

IQTreeMix::IQTreeMix(Params &params, Alignment *aln, vector<IQTree*> &trees) : IQTree(aln) {
    size_t i;

    clear();
    weights.clear();

    // store the trees and initialize tree-weights
    ntree = trees.size();
    double init_weight = 1.0 / (double) ntree;
    for (i=0; i<ntree; i++) {
        push_back(trees[i]);
        weights.push_back(init_weight);
    }
    
    // allocate memory for the arrays
    ptn_like_cat = new double[size() * aln->getNPattern()];
    _ptn_like_cat = new double[size() * aln->getNPattern()];
    patn_freqs = new int[aln->getNPattern()];
    patn_isconst = new int[aln->getNPattern()];
    
    // get the pattern frequencies
    aln->getPatternFreq(patn_freqs);
    
    // get whether the pattern is constant
    for (i=0; i<aln->getNPattern(); i++) {
        patn_isconst[i] = aln->at(i).isConst();
    }
    
    // number of optimization steps, default: number of Trees * 2
    // optimize_steps = 2 * size();
    // optimize_steps = 100;
    optimize_steps = 1000;
    
    // initialize the tree weights as non-fixed, that means it needs to be optimized
    isTreeWeightFixed = false;
    
    weightGrpExist = false;
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
    
    for (i=0; i<size(); i++) {
        delete (at(i));
    }
    if (ptn_like_cat != NULL) {
        delete[] ptn_like_cat;
    }
    if (_ptn_like_cat != NULL) {
        delete[] _ptn_like_cat;
    }
    if (patn_freqs != NULL) {
        delete[] patn_freqs;
    }
    if (patn_isconst != NULL) {
        delete[] patn_isconst;
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
        if (i==0 || (s[0]=='F')) {
            if (model.length() > 0)
                model.append("+");
            model.append(s);
        } else {
            if (siteRate.length() > 0)
                siteRate.append("+");
            siteRate.append(s);
        }
        i++;
    }
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
    size_t t_pos, i, k, d_pos, d_pos2;
    string s;
    vector<string> model_array, submodel_array, weight_array;
    
    rmSpace(modelName);
    treemix_model = modelName;
    model_names.clear();
    siterate_names.clear();
    isLinkSiteRate = true; // initialize to true
    
    // check how many trees
    //cout << "[IQTreeMix::separateModel] modelName = " << modelName << endl;
    t_pos = modelName.rfind("+T");
    if (t_pos == string::npos) {
        outError("This model is not a tree mixture model, because there is no '+T'");
    }
    if (t_pos < modelName.length()-2) {
        d_pos = t_pos+2;
        while (d_pos < modelName.length() && isdigit(modelName[d_pos]))
            d_pos++;
        if (d_pos - t_pos > 2) {
            ntree = atoi(modelName.substr(t_pos+2, d_pos-t_pos-2).c_str());
            if (ntree < size()) {
                cout << endl << "NOTE: Only the first " << ntree << " trees in the treefile are considered, because '" <<  modelName.substr(t_pos+1) << "' has been specified!" << endl;
                resize(ntree);
                weights.resize(ntree);
                double init_weight = 1.0 / (double) ntree;
                for (i=0; i<ntree; i++) {
                    weights[i] = init_weight;
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
                for (i=0; i<ntree; i++)
                    weights[i] = weights[i] / sum;
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
    
    /*
    if (ntree <= 1) {
        outError("For tree mixture model, number of trees has to be at least 2.");
    }*/
    
    // remove the '+Txx'
    modelName = modelName.substr(0, t_pos);
    
    // break the whole name according to '+'
    separateStr(modelName, model_array, '+');
    
    // check each model/siterate
    for (i=0; i<model_array.size(); i++) {
        s = model_array[i];
        if (s.length() == 0) {
            continue;
        } else if (s.length() > 5 && s.substr(0,4) == "MIX{" && s.substr(s.length()-1,1) == "}") {
            // mixture model
            s = s.substr(4,s.length()-5); // remove the beginning "MIX{" and the ending "}"
            if (i==0) {
                // unlinked models (while site rates may or may not be linked)
                bool siteRateAppear = false;
                string curr_model, curr_siterate;
                separateStr(s, submodel_array, ',');
                for (k=0; k<submodel_array.size(); k++) {
                    divideModelNSiteRate(submodel_array[k], curr_model, curr_siterate);
                    model_names.push_back(curr_model);
                    siterate_names.push_back(curr_siterate);
                    if (curr_siterate.length() > 0) {
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
                for (k=0; k<submodel_array.size(); k++) {
                    siterate_names.push_back(submodel_array[k]);
                }
            } else {
                outError("Error! The model: " + treemix_model + " is not correctly specified. Are you using too many 'MIX'?");
            }
        } else if (i==0) {
            // not a mixture model
            model_names.push_back(s);
        } else if (s.length() <= 2 && s[0] == 'F') {
            // F or FO
            if (model_names.size() > 1) {
                outError("Error! 'F' is linked, but the model is unlinked");
            } else if (model_names.size() == 1){
                model_names[0].append("+" + s);
            } else {
                outError("Error! 'F' appears before the model does");
            }
        } else {
            // assume this is the site rate model
            if (siterate_names.size() > 1) {
                outError("Error! '" + s + "' is linked, but the site rates are unlinked");
            } else if (siterate_names.size() == 1) {
                siterate_names[0].append("+" + s);
            } else {
                siterate_names.push_back(s);
            }
        }
    }
    if (model_names.size() == 0) {
        outError("Error! It seems no model is defined.");
    }
    isLinkModel = (model_names.size() == 1);
    if (siterate_names.size() == 0) {
        anySiteRate = false;
    } else {
        anySiteRate = true;
        isLinkSiteRate = (siterate_names.size() == 1);
    }
    
    // check correctness
    if (model_names.size() > 1 && model_names.size() != ntree) {
        outError("Error! The number of submodels specified in the mixture does not match with the tree number");
    }
    if (siterate_names.size() > 1 && siterate_names.size() != ntree) {
        outError("Error! The number of site rates specified in the mixture does not match with the tree number");
    }
}

void IQTreeMix::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    size_t i;
    string curr_model;

    models.clear();
    site_rates.clear();
    separateModel(model_name);

    // initialize the models
    for (i=0; i<ntree; i++) {
        if (isLinkModel) {
            curr_model = model_names[0];
        } else {
            curr_model = model_names[i];
        }
        if (anySiteRate) {
            if (isLinkSiteRate) {
                if (siterate_names[0] != "E")
                    curr_model += "+" + siterate_names[0];
                params.optimize_alg_gammai = "BFGS";
                params.optimize_alg_freerate = "2-BFGS";
                params.optimize_alg_mixlen = "BFGS";
            } else {
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
            ModelSubst *m = at(i)->getModelFactory()->model;
            // delete m (need to be addressed);
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
        if (isLinkSiteRate) {
            site_rates.push_back(at(0)->getModelFactory()->site_rate);
            for (i=1; i<ntree; i++) {
                RateHeterogeneity *r = at(i)->getModelFactory()->site_rate;
                // delete r (need to be addressed);
                at(i)->getModelFactory()->site_rate = site_rates[0];
                at(i)->setRate(site_rates[0]);
            }
            // for linked site rate model, set its tree to this tree
            site_rates[0]->setTree(this);
        } else {
            for (i=0; i<ntree; i++) {
                site_rates.push_back(at(i)->getModelFactory()->site_rate);
            }
            // for unlinked site rate model, set their trees to the corresponding trees
            for (i=0; i<ntree; i++) {
                site_rates[i]->setTree(at(i));
            }
        }
    }
    
    /*
    // initialize the tree weights using parsimony scores
    double w_sum = 0.0;
    cout << "parismony scores:";
    for (i=0; i<ntree; i++) {
        int par_score = at(i)->computeParsimony();
        cout << " " << par_score;
        weights[i] = 1.0 / (double) par_score;
        w_sum += weights[i];
    }
    cout << endl;
    cout << "tree weights are initialized as:";
    for (i=0; i<ntree; i++) {
        weights[i] = weights[i] / w_sum;
        cout << " " << weights[i];
    }
    cout << endl;
    */
}

double IQTreeMix::computeLikelihood(double *pattern_lh) {
    double* pattern_lh_tree;
    size_t i,j,ptn,t;
    size_t nptn,ntree;
    double logLike = 0.0;
    double subLike;
    double score;
    PhyloTree* ptree;
    
    nptn = aln->getNPattern();
    ntree = size();

    // compute likelihood for each tree
    pattern_lh_tree = _ptn_like_cat;
    for (t=0; t<ntree; t++) {
        // save the site rate's tree
        ptree = at(t)->getRate()->getTree();
        // set the tree t as the site rate's tree
        // and compute the likelihood values
        at(t)->getRate()->setTree(at(t));
        at(t)->initializeAllPartialLh();
        score = at(t)->computeLikelihood(pattern_lh_tree);
        at(t)->clearAllPartialLH();
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
        // cout << "[IQTreeMix::computeLikelihood] Tree " << t+1 << " : " << score << endl;
        pattern_lh_tree += nptn;
    }

    // reorganize the array
    i=0;
    for (t=0; t<ntree; t++) {
        j=t;
        for (ptn=0; ptn<nptn; ptn++) {
            ptn_like_cat[j] = exp(_ptn_like_cat[i]);
            i++;
            j+=ntree;
        }
    }
    
    // compute the total likelihood
    i=0;
    for (ptn=0; ptn<nptn; ptn++) {
        subLike = 0.0;
        for (t=0; t<ntree; t++) {
            subLike += ptn_like_cat[i] * weights[t];
            i++;
        }
        if (pattern_lh) {
            pattern_lh[ptn] = subLike;
        }
        // cout << ptn << "\t" << log(subLike) << "\t" << patn_freqs[ptn] << endl;
        logLike += log(subLike) * (double) patn_freqs[ptn];
    }
    // cout << "[IQTreeMix::computeLikelihood] log-likelihood: " << logLike << endl;

    return logLike;
}

// compute the overall likelihood value by combining all the existing likelihood values of the trees
double IQTreeMix::computeLikelihood_combine() {
    double* pattern_lh_tree;
    size_t i,j,ptn,t;
    size_t nptn,ntree;
    double logLike = 0.0;
    double subLike;
    double score;
    PhyloTree* ptree;
    
    nptn = aln->getNPattern();
    ntree = size();

    // compute the total likelihood
    i=0;
    for (ptn=0; ptn<nptn; ptn++) {
        subLike = 0.0;
        for (t=0; t<ntree; t++) {
            subLike += ptn_like_cat[i] * weights[t];
            i++;
        }
        // cout << ptn << "\t" << log(subLike) << "\t" << patn_freqs[ptn] << endl;
        logLike += log(subLike) * (double) patn_freqs[ptn];
    }
    // cout << "[IQTreeMix::computeLikelihood] log-likelihood: " << logLike << endl;

    return logLike;
}

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
    size_t i,ptn,t;
    size_t nptn,ntree;
    double subLike;

    computeLikelihood(pattern_lh);
}

void IQTreeMix::initializeAllPartialLh() {
    size_t i;
    // IQTree::initializeAllPartialLh();
    for (i=0; i<size(); i++) {
        at(i)->initializeAllPartialLh();
    }
}

void IQTreeMix::deleteAllPartialLh() {
    size_t i;
    // IQTree::deleteAllPartialLh();
    for (i=0; i<size(); i++) {
        at(i)->deleteAllPartialLh();
    }
}

void IQTreeMix::clearAllPartialLH(bool make_null) {
    size_t i;
    // IQTree::clearAllPartialLH(make_null);
    for (i=0; i<size(); i++) {
        at(i)->clearAllPartialLH(make_null);
    }
}

/**
        optimize all branch lengths of the tree
        @param iterations number of iterations to loop through all branches
        @return the likelihood of the tree
 */
double IQTreeMix::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
    size_t i;
    PhyloTree* ptree;
    
    for (i=0; i<size(); i++) {
        // cout << "[IQTreeMix::optimizeAllBranches] i = " << i << endl;
        // because the tree of the site rate is pointing to this tree
        // save the tree of the site rate
        ptree = at(i)->getRate()->getTree();
        at(i)->getRate()->setTree(at(i));
        at(i)->optimizeAllBranches(my_iterations, tolerance, maxNRStep);
        // restore the tree of the site rate
        at(i)->getRate()->setTree(ptree);
    }
    return computeLikelihood();
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
        compute the updated tree weights according to the likelihood values along each site
        prerequisite: computeLikelihood() has been invoked

 */
double IQTreeMix::optimizeTreeWeightsByEM(double* pattern_mix_lh, double gradient_epsilon, int max_steps) {
    size_t nptn, ntree, ptn, c;
    double *this_lk_cat;
    double lk_ptn;
    double prev_score, score;
    int step;

    nptn = aln->getNPattern();
    ntree = size();
    
    initializeAllPartialLh();
    prev_score = computeLikelihood();
    clearAllPartialLH();
    
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
        
        score = computeLikelihood_combine();

        if (score < prev_score + gradient_epsilon) {
            // converged
            break;
        }
        prev_score = score;

    }
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

    // special case: ndim = 1, i.e. all tree weights are forced the same
    if (ndim == 1) {
        for (i=0; i<size(); i++) {
            weights[i] = 1.0 / size();
        }
        initializeAllPartialLh();
        score = computeLikelihood();
        clearAllPartialLH();
        return score;
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
        for (i = 0; i < size(); i++)
            this->weights[i] = relative_weights[i];
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
        params.min_branch_length = 1e-6;
        if (size() > 0) {
            if (!at(0)->isSuperTree() && at(0)->getAlnNSite() >= 100000) {
                params.min_branch_length = 0.1 / (at(0)->getAlnNSite());
                num_prec = max((int)ceil(-log10(Params::getInstance().min_branch_length))+1, 6);
                for (i=0; i<size(); i++)
                    at(i)->num_precision = num_prec;
                cout.precision(12);
                cout << "NOTE: minimal branch length is reduced to " << params.min_branch_length << " for long alignment" << endl;
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
    }
}

/** set pointer of params variable */
void IQTreeMix::setParams(Params* params) {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->setParams(params);
    }
    this->params = params;
};

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
    showTree();
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

/**
    test the best number of threads
*/
int IQTreeMix::testNumThreads() {
    return at(0)->testNumThreads();
}

/**
    optimize each tree separately
 */
void IQTreeMix::OptimizeTreesSeparately(bool printInfo, double gradient_epsilon) {
    size_t step, i, ptn, nptn;
    string result = "";
    string curr_res;
    PhyloTree* ptree;
    double t_score, score, prev_score;
    
    int n = 1;

    score = computeLikelihood();
    prev_score = score;
    for (step = 0; step < optimize_steps; step++) {
        
        // optimize the linked site rate model
        if (anySiteRate && isLinkSiteRate) {
            t_score = site_rates[0]->optimizeParameters(gradient_epsilon);
            if (t_score < 0.0)
                score = t_score;
        }

        // optimize the linked subsitution model
        if (isLinkModel) {
            t_score = models[0]->optimizeParameters(gradient_epsilon);
            if (t_score < 0.0)
                score = t_score;
        }

        // optimize the unlinked site-rate models one by one
        if (anySiteRate && !isLinkSiteRate) {
            for (i=0; i<site_rates.size(); i++) {
                site_rates[i]->optimizeParameters(gradient_epsilon);
            }
        }

        // optimize the unlinked subsitution models one by one
        if (!isLinkModel) {
            for (i=0; i<models.size(); i++) {
                models[i]->optimizeParameters(gradient_epsilon);
            }
        }

        // optimize tree branches
        score = optimizeAllBranches(n, gradient_epsilon);  // loop max n times
        
        score = computeLikelihood();
            
        if (score < prev_score + gradient_epsilon) {
            // converged
            break;
        }
        prev_score = score;
    }
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
    size_t i, j, ntree, nptn;
    
    ntree = size();
    nptn = aln->ordered_pattern.size();
    UINT* ptn_scores = new UINT[ntree * nptn];
    UINT* curr_ptn_scores;
    UINT min_par_score;
    vector<int> tree_with_min_pars;
    double weight_sum;
    
    // compute the parsimony scores along patterns for each tree
    for (i=0; i<ntree; i++) {
        curr_ptn_scores = ptn_scores + i * nptn;

        at(i)->initCostMatrix(CM_UNIFORM);
        at(i)->setParsimonyKernel(params->SSE);
        at(i)->computeTipPartialParsimony();
        at(i)->computeParsimonyOutOfTreeSankoff(curr_ptn_scores);
    }
    
    // reset the tree weights
    for (i=0; i<ntree; i++) {
        weights[i] = 0.0;
    }
    
    // estimate the tree weights
    for (i=0; i<nptn; i++) {
        min_par_score = ptn_scores[i];
        tree_with_min_pars.clear();
        tree_with_min_pars.push_back(0);
        for (j=1; j<ntree; j++) {
            if (ptn_scores[i+j*nptn] < min_par_score) {
                tree_with_min_pars.clear();
                tree_with_min_pars.push_back(j);
                min_par_score = ptn_scores[i+j*nptn];
            } else if (ptn_scores[i+j*nptn] == min_par_score) {
                tree_with_min_pars.push_back(j);
            }
        }
        if (tree_with_min_pars.size() < ntree) {
            for (j=0; j<tree_with_min_pars.size(); j++) {
                weights[tree_with_min_pars[j]] += aln->ordered_pattern[i].frequency;
            }
        }
    }
    
    // normalize the tree weights
    weight_sum = 0.0;
    for (i=0; i<ntree; i++) {
        weight_sum += weights[i];
    }
    for (i=0; i<ntree; i++) {
        weights[i] = weights[i] / weight_sum;
    }
    
    // If there are multiple tree weights belonging to the same group
    // set all the tree weights of the same group to their average
    checkWeightGrp();
    
    // show the initial tree weights
    cout << "According to the parsimony scores along the sites, the tree weights are initialized to:";
    for (i=0; i<ntree; i++) {
        cout << " " << weights[i];
    }
    cout << endl;

    delete[] ptn_scores;
}

string IQTreeMix::optimizeModelParameters(bool printInfo, double logl_epsilon) {
    
    size_t i, ntree, nptn, ptn;
    int step, n, substep;//, nsubstep;
    double* pattern_mix_lh;
    double gradient_epsilon = 0.0001;
    double epsilon_start = 0.0001;
    double epsilon_step = 0.1;
    double curr_epsilon = epsilon_start;
    double prev_score, prev_score2, score, t_score;
    PhyloTree *ptree;

    n = 1;
    ntree = size();
    nptn = aln->getNPattern();
    
    // allocate memory
    pattern_mix_lh = new double[ntree * nptn];
    
    if (!isTreeWeightFixed) {
        // initialize the tree weights according to parsimony scores along the sites
        initializeTreeWeights();
    }
    
    prev_score = score = -DBL_MAX;

    for (step = 0; step < optimize_steps; step++) {
        
        if (step > 0) {
            // reset the ptn_freq array to the original frequencies of the patterns
            for (i = 0; i < ntree; i++) {
                for (ptn = 0; ptn < nptn; ptn++) {
                    at(i)->ptn_freq[ptn] = patn_freqs[ptn];
                }
            }
        }

        // optimize the linked site rate model
        if (anySiteRate && isLinkSiteRate) {
            site_rates[0]->optimizeParameters(curr_epsilon);
            if (siterate_names[0].find("R") != string::npos) {
                site_rates[0]->rescaleRates();
            }
            // cout << "after optimizing linked site rate model, likelihood = " << score << "(t_score=" << t_score << ")" << endl;
        }

        // optimize the linked subsitution model
        if (isLinkModel) {
            models[0]->optimizeParameters(curr_epsilon);
            // cout << "after optimizing linked subsitution model, likelihood = " << score << "(t_score=" << t_score << ")" << endl;
        }
        
        score = computeLikelihood();
        prev_score2 = score;
        
        for (substep = 0; substep<step+1; substep++) {
            
            // compute the ptn_freq array according to the posterior probabilities along each site for each tree
            computeFreqArray(pattern_mix_lh, false);

            if (params->fixed_branch_length != BRLEN_FIX) {
                // optimize tree branches
                score = optimizeAllBranches(1, curr_epsilon);  // loop max n times
                // cout << "after optimizing branches, likelihood = " << score << endl;
            }

            // optimize the unlinked subsitution models one by one
            if (!isLinkModel) {
                for (i=0; i<models.size(); i++) {
                    models[i]->optimizeParameters(curr_epsilon);
                }
                score = computeLikelihood();
                // cout << "after optimizing unlinked subsitution model, likelihood = " << score << endl;
            }

            // optimize the unlinked site-rate models one by one
            if (anySiteRate && !isLinkSiteRate) {
                for (i=0; i<site_rates.size(); i++) {
                    site_rates[i]->optimizeParameters(curr_epsilon);
                    if (siterate_names[i].find("R") != string::npos) {
                        site_rates[i]->rescaleRates();
                    }
                }
                score = computeLikelihood();
                // cout << "after optimizing unlinked site-rate model, likelihood = " << score << endl;
            }

            // optimize tree weights
            if (!isTreeWeightFixed) {
                if (weightGrpExist || params->optimize_alg_treeweight == "BFGS") {
                    score = optimizeTreeWeightsByBFGS(curr_epsilon);
                } else {
                    score = optimizeTreeWeightsByEM(pattern_mix_lh, curr_epsilon, 1);  // loop max n times
                }
            }

            score = computeLikelihood();
            if (score < prev_score2 + curr_epsilon) {
                // converged
                break;
            }
            prev_score2 = score;
        }
        
        cout << "step= " << step << " substep = " << substep << " curr_epsilon=" << curr_epsilon << " score=" << score << endl;

        if (score < prev_score + curr_epsilon && curr_epsilon > gradient_epsilon) {
            // update the epsilon value
            curr_epsilon = curr_epsilon * epsilon_step;
        } else if (score < prev_score + gradient_epsilon) {
            // converged
            break;
        }

        prev_score = score;
    }

    setCurScore(score);

    delete[] pattern_mix_lh;

    // show the weights
    cout << "Final estimation on weights:";
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
    for (i=0; i<size(); i++)
        at(i)->printTree(fout);
    setRootNode(params->root, false);
    fout.close();
    if (verbose_mode >= VB_MED)
        cout << "Best tree printed to " << tree_file_name << endl;
}

string IQTreeMix::getTreeString() {
    stringstream tree_stream;
    size_t i;
    
    for (i=0; i<size(); i++)
        at(i)->printTree(tree_stream, WT_TAXON_ID + WT_BR_LEN + WT_SORT_TAXA);
    return tree_stream.str();
}

// return the average of the tree lengths
double IQTreeMix::treeLength(Node *node, Node *dad) {
    double len = 0.0;
    size_t i;
    for (i=0; i<size(); i++)
        len += at(i)->treeLength();
    return len / size();
}

// return the average length of all internal branches
double IQTreeMix::treeLengthInternal( double epsilon, Node *node, Node *dad) {
    double len = 0.0;
    size_t i;
    for (i = 0; i < size(); i++)
        len += at(i)->treeLengthInternal(epsilon);
    return len / size();
}

int IQTreeMix::getNParameters() {
    int df = 0;
    size_t i;
    for (i=0; i<models.size(); i++) {
        df += (models[i]->getNDim() + models[i]->getNDimFreq());
    }
    for (i=0; i<site_rates.size(); i++) {
        df += site_rates[i]->getNDim();
    }
    for (i=0; i<size(); i++) {
        df += at(i)->getNBranchParameters(BRLEN_OPTIMIZE);
    }
    // for tree weights
    df += (size() - 1);
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
int IQTreeMix::printTree(ostream &out, int brtype, Node *node, Node *dad) {
    size_t i;
    int value = 0;
    for (i=0; i<size(); i++) {
        out << "Tree " << i+1 << ":" << endl;
        value = at(i)->printTree(out, brtype, node, dad);
    }
    return value;
}

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

// get posterior probabilities along each site for each tree
void IQTreeMix::getPostProb(double* pattern_mix_lh, bool need_computeLike) {
    size_t ntree, nptn, i, ptn, c;
    double* this_lk_cat;
    double lk_ptn;

    ntree = size();
    nptn = aln->getNPattern();

    if (need_computeLike) {
        initializeAllPartialLh();
        computeLikelihood();
        clearAllPartialLH();
    }

    memcpy(pattern_mix_lh, ptn_like_cat, nptn*ntree*sizeof(double));

    // multiply pattern_mix_lh by tree weights
    i = 0;
    for (ptn = 0; ptn < nptn; ptn++) {
        for (c = 0; c < ntree; c++) {
            pattern_mix_lh[i] *= weights[c];
            i++;
        }
    }

    for (ptn = 0; ptn < nptn; ptn++) {
        this_lk_cat = pattern_mix_lh + ptn*ntree;
        lk_ptn = 0.0;
        for (c = 0; c < ntree; c++) {
            lk_ptn += this_lk_cat[c];
        }
        ASSERT(lk_ptn != 0.0);
        lk_ptn = patn_freqs[ptn] / lk_ptn;

        // transform pattern_mix_lh into posterior probabilities of each category
        for (c = 0; c < ntree; c++) {
            this_lk_cat[c] *= lk_ptn;
        }
    }
}

// update the ptn_freq array according to the posterior probabilities along each site for each tree
void IQTreeMix::computeFreqArray(double* pattern_mix_lh, bool need_computeLike) {
    size_t i, ptn;
    PhyloTree* tree;
    size_t ntree = size();
    size_t nptn = aln->getNPattern();

    // get posterior probabilities along each site for each tree
    getPostProb(pattern_mix_lh, need_computeLike);

    for (i = 0; i < ntree; i++) {
        tree = at(i);
        // initialize likelihood
        // tree->initializeAllPartialLh();
        // copy posterior probability into ptn_freq
        // tree->computePtnFreq();
        double *this_lk_cat = pattern_mix_lh+i;
        for (ptn = 0; ptn < nptn; ptn++)
            tree->ptn_freq[ptn] = this_lk_cat[ptn*ntree];
    }
}


double IQTreeMix::targetFunk(double x[]) {
    getVariables(x);
    return -computeLikelihood_combine();
}

// read the tree weights and write into "variables"
void IQTreeMix::setVariables(double *variables) {
    // for tree weights
    size_t i;
    size_t ndim = weight_group_member.size();
    for (i=0; i<ndim; i++) {
        variables[i+1] = tmp_weights[i];
    }
}

// read the "variables" and write into tree weights
void IQTreeMix::getVariables(double *variables) {
    // for tree weights
    size_t i,j;
    size_t ndim = weight_group_member.size();
    double sum = 0.0;
    double w;
    for (i=0; i<ndim; i++) {
        tmp_weights[i] = variables[i+1];
        sum += tmp_weights[i] * weight_group_member[i].size();
    }
    for (i=0; i<ndim; i++) {
        w = tmp_weights[i] / sum;
        for (j=0; j<weight_group_member[i].size(); j++) {
            weights[weight_group_member[i].at(j)] = w;
        }
    }
}

// set the bounds
void IQTreeMix::setBounds(double *lower_bound, double *upper_bound, bool* bound_check) {
    size_t i;
    size_t ndim = weight_group_member.size();
    for (i=0; i<ndim; i++) {
        lower_bound[i+1] = MIN_PROP;
        upper_bound[i+1] = MAX_PROP;
        bound_check[i+1] = false;
    }
}

// get the dimension of the variables (for tree weights)
int IQTreeMix::getNDim() {
    size_t s;
    if (weight_group_member.size() > 0)
        s = weight_group_member.size();
    else
        s = size();
    return s;
}

// show the log-likelihoods and posterior probabilties for each tree along the sites
void IQTreeMix::showLhProb(ofstream& out) {
    double* pattern_lh_tree;
    double* curr_ptn_lh;
    double* post_prob;
    size_t t,site,idx;
    size_t nsite,nptn,ntree;
    PhyloTree* ptree;
    double sum;
    
    // pattern_index[i]
    
    nsite = aln->getNSite();
    nptn = aln->getNPattern();
    ntree = size();

    IntVector pattern_index;
    aln->getSitePatternIndex(pattern_index);

    // compute likelihood for each tree
    pattern_lh_tree = new double[nptn * ntree];
    curr_ptn_lh = pattern_lh_tree;
    for (t=0; t<ntree; t++) {
        // save the site rate's tree
        ptree = at(t)->getRate()->getTree();
        // set the tree t as the site rate's tree
        // and compute the likelihood values
        at(t)->getRate()->setTree(at(t));
        at(t)->initializeAllPartialLh();
        at(t)->computeLikelihood(curr_ptn_lh);
        at(t)->clearAllPartialLH();
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
        curr_ptn_lh += nptn;
    }
    
    // for posterior probabilities
    post_prob = new double[ntree];
    
    // print out the log-likelihoods and posterior probabilties for each tree along the sites
    out << "site,log-like";
    for (t=0; t<ntree; t++) {
        out << ",log-like tree " << t+1;
    }
    for (t=0; t<ntree; t++) {
        out << ",post-prob tree " << t+1;
    }
    out << endl;
    for (site=0; site<nsite; site++) {
        out << site+1;
        idx = pattern_index[site];
        curr_ptn_lh = pattern_lh_tree;
        sum = 0.0;
        for (t=0; t<ntree; t++) {
            post_prob[t] = exp(curr_ptn_lh[idx]) * weights[t];
            sum += post_prob[t];
            curr_ptn_lh += nptn;
        }
        out << "," << log(sum); // log-likelihood of the site
        curr_ptn_lh = pattern_lh_tree;
        for (t=0; t<ntree; t++) {
            out << "," << curr_ptn_lh[idx]; // log-likelihood of the site for tree t
            curr_ptn_lh += nptn;
        }
        for (t=0; t<ntree; t++) {
            post_prob[t] = post_prob[t] / sum;
            out << "," << post_prob[t]; // posterior probability of the site for tree t
        }
        out << endl;
    }
    
    // free the memory of the array
    delete[] pattern_lh_tree;
    delete[] post_prob;
}
