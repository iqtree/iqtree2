//
//  iqtreemix.cpp
//  tree
//
//  Created by Thomas Wong on 14/12/20.
//

#include "iqtreemix.h"
const double MIN_PROP = 0.001;
const double MAX_PROP = 1000.0;
const double MIN_LEN = 1e-4;
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
    int i;
    if (m=="I" || m=="G" || m=="R")
        return true;
    if (m.length() > 1 && (m[0]=='G' || m[0]=='R')) {
        for (i=1; i<m.length(); i++) {
            if (!isdigit(m[i]))
                return false;
        }
        return true;
    }
    return false;
}

IQTreeMix::IQTreeMix() : IQTree() {
    patn_freqs = NULL;
    ptn_freq = NULL;
    patn_isconst = NULL;
    patn_parsimony = NULL;
    ptn_like_cat = NULL;
    _ptn_like_cat = NULL;
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
    nptn = aln->getNPattern();
    size_t mem_size = get_safe_upper_limit(nptn);
    size_t block_size = get_safe_upper_limit(nptn) * ntree;
    patn_freqs = aligned_alloc<int>(mem_size);
    ptn_freq = aligned_alloc<double>(mem_size);
    patn_isconst = aligned_alloc<int>(mem_size);
    ptn_like_cat = aligned_alloc<double>(block_size);
    _ptn_like_cat = aligned_alloc<double>(block_size);
    patn_parsimony = aligned_alloc<int>(block_size);

    // get the pattern frequencies
    aln->getPatternFreq(patn_freqs);
    for (i=0; i<nptn; i++) {
        ptn_freq[i] = (double) patn_freqs[i];
    }
    
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
        aligned_free(ptn_like_cat);
    }
    if (_ptn_like_cat != NULL) {
        aligned_free(_ptn_like_cat);
    }
    if (patn_freqs != NULL) {
        aligned_free(patn_freqs);
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

    // if it is a edge-length-restricted model, build the branch ID
    if (isEdgeLenRestrict) {
        // build the branch ID
        computeBranchID();
        // show trees
        // showTree();
    }
    
    // show summary
    cout << endl;
    if (model_names.size() == 1) {
        cout << "Linked substitution model:" << endl;
        isLinkModel = true;
    } else {
        cout << "Unlinked substitution models:" << endl;
        isLinkModel = false;
    }
    for (i = 0; i < model_names.size(); i++) {
        cout << "   " << model_names[i] << endl;
    }
    cout << endl;
    if (anySiteRate) {
        if (siterate_names.size() == 1) {
            cout << "Linked RHS model:" << endl;
            isLinkSiteRate = true;
        } else {
            cout << "Unlinked RHS models:" << endl;
            isLinkSiteRate = false;
        }
        for (i = 0; i < siterate_names.size(); i++) {
            cout << "   " << siterate_names[i] << endl;
        }
    }
    cout << endl;
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
                at(i)->getModelFactory()->site_rate = site_rates[0];
                at(i)->setRate(site_rates[0]);
            }
            // for linked site rate model, set its tree to this tree
            site_rates[0]->setTree(this);
            for (i=0; i<ntree; i++) {
                site_rate_trees.push_back(this);
            }
        } else {
            for (i=0; i<ntree; i++) {
                site_rates.push_back(at(i)->getModelFactory()->site_rate);
            }
            // for unlinked site rate model, set their trees to corresponding trees
            for (i=0; i<ntree; i++) {
                site_rates[i]->setTree(at(i));
                site_rate_trees.push_back(at(i));
            }
        }
    }
}

// this function is designed for the situation that
// only one tree has been updated.
double IQTreeMix::computeLikelihood_oneTreeUpdated(int whichTree) {
    double* pattern_lh_tree;
    size_t i,j,ptn;
    double logLike = 0.0;
    double subLike;
    double score;
    size_t t = whichTree;
    PhyloTree* ptree;
    
    // compute likelihood for each tree
    pattern_lh_tree = _ptn_like_cat;
    pattern_lh_tree += t*nptn;
    
    ptree = at(t)->getRate()->getTree();
    // set the tree t as the site rate's tree
    // and compute the likelihood values
    at(t)->getRate()->setTree(at(t));
    at(t)->initializeAllPartialLh();
    at(t)->clearAllPartialLH();
    at(t)->computeLikelihood(pattern_lh_tree);
    // set back the prevoius site rate's tree
    at(t)->getRate()->setTree(ptree);
    // cout << "[IQTreeMix::computeLikelihood] Tree " << t+1 << " : " << score << endl;

    // reorganize the array
    j=t;
    for (ptn=0; ptn<nptn; ptn++) {
        ptn_like_cat[j] = exp(pattern_lh_tree[ptn]);
        j+=ntree;
    }

    // compute the overall likelihood value by combining all the existing likelihood values of the trees
    return computeLikelihood_combine();
}

double IQTreeMix::computeLikelihood(double *pattern_lh) {
    double* pattern_lh_tree;
    size_t i,j,ptn,t;
    double logLike = 0.0;
    double subLike;
    double score;
    PhyloTree* ptree;
    
    // compute likelihood for each tree
    pattern_lh_tree = _ptn_like_cat;
    for (t=0; t<ntree; t++) {
        // save the site rate's tree
        ptree = at(t)->getRate()->getTree();
        // set the tree t as the site rate's tree
        // and compute the likelihood values
        at(t)->getRate()->setTree(at(t));
        at(t)->initializeAllPartialLh();
        at(t)->clearAllPartialLH();
        at(t)->computeLikelihood(pattern_lh_tree);
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
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
    
    // compute the overall likelihood value by combining all the existing likelihood values of the trees
    return computeLikelihood_combine(pattern_lh);
}

double IQTreeMix::computePatternLhCat(SiteLoglType wsl) {
    // this function is only for the linked mixture model
    size_t nmix, t, p, m, i;
    size_t tot_mem_size;
    double ptnLike, score;
    int index, indexlh;

    int numStates = at(0)->model->num_states;
    size_t mem_size = get_safe_upper_limit(getAlnNPattern()) + max(get_safe_upper_limit(numStates),
        get_safe_upper_limit(at(0)->model_factory->unobserved_ptns.size()));
    if (!_pattern_lh_cat) {
        tot_mem_size = mem_size * at(0)->site_rate->getNDiscreteRate() * ((at(0)->model_factory->fused_mix_rate)? 1 : at(0)->model->getNMixtures());
        _pattern_lh_cat = aligned_alloc<double>(tot_mem_size);
    }
    
    // compute _pattern_lh_cat for each tree
    for (t = 0; t < ntree; t++) {
        if (t==0) {
            nmix = at(t)->getModel()->getNMixtures();
        } else if (nmix != at(t)->getModel()->getNMixtures()) {
            outError("The number of classes inside each tree is not the same!");
        }
        at(t)->computePatternLhCat(wsl);
    }
    
    // compute the overall _pattern_lh_cat
    for (p = 0; p < nptn * nmix; p++) {
        ptnLike = 0.0;
        for (t = 0; t < ntree; t++) {
            ptnLike += at(t)->_pattern_lh_cat[p] * weights[t];
        }
        _pattern_lh_cat[p] = ptnLike;
    }
    
    score = computeLikelihood();
    return score;
}

// compute the overall likelihood value by combining all the existing likelihood values of the trees
double IQTreeMix::computeLikelihood_combine(double *pattern_lh) {
    double* pattern_lh_tree;
    size_t i,ptn,t;
    double logLike = 0.0;
    double subLike;
    double ptnLike;
    PhyloTree* ptree;
    
    // compute the total likelihood
    i=0;
    for (ptn=0; ptn<nptn; ptn++) {
        subLike = 0.0;
        for (t=0; t<ntree; t++) {
            subLike += ptn_like_cat[i] * weights[t];
            i++;
        }
        // cout << ptn << "\t" << log(subLike) << "\t" << patn_freqs[ptn] << endl;
        ptnLike = log(subLike);
        logLike += ptnLike * (double) patn_freqs[ptn];
        if (pattern_lh != NULL) {
            pattern_lh[ptn] = ptnLike;
        }
    }

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
    computeLikelihood(pattern_lh);
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
    size_t i;
    
    for (i=0; i<size(); i++)
        optimizeAllBranchesOneTree(i, my_iterations, tolerance, maxNRStep);

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
double IQTreeMix::optimizeTreeWeightsByEM(double* pattern_mix_lh, double gradient_epsilon, int max_steps, bool& tree_weight_converge) {
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
        
        // for weight-restricted condition
        checkWeightGrp();
        
        // if all differences between the new weights and the old weights <= WEIGHT_EPSILON, then tree_weight_converge = true
        tree_weight_converge = true;
        for (c = 0; c< ntree && tree_weight_converge; c++)
            if (fabs(weights[c] - tmp_weights[c]) > WEIGHT_EPSILON)
                tree_weight_converge = false;

        score = computeLikelihood_combine();

        if (score < prev_score + gradient_epsilon) {
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
        for (i=0; i<size(); i++) {
            weights[i] = 1.0 / size();
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

/**
    test the best number of threads
*/
int IQTreeMix::testNumThreads() {
    return at(0)->testNumThreads();
}

/**
    optimize tree k separately
 */
void IQTreeMix::optimizeTreeSeparately(int k, bool printInfo, double gradient_epsilon) {
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
            score = at(k)->optimizeAllBranches(n, gradient_epsilon);
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
void IQTreeMix::optimizeTreesSeparately(bool printInfo, double gradient_epsilon) {
    int i;

    for (i=0; i<ntree; i++) {
        // optimize tree i
        optimizeTreeSeparately(i, printInfo, gradient_epsilon);
    }
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
                weights[tree_with_min_pars[j]] += patn_freqs[i];
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
    checkWeightGrp();
    
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
        // cout << ", " << patn_freqs[i] << endl;
            tree_with_min_pars.push_back(min_pars);
            freqs.push_back(patn_freqs[i]);
    }
    // cout << "--------------------" << endl;

    // initialize the tree weights
    for (i=0; i<ntree; i++) {
        weights[i] = 1.0 / ntree;
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
    checkWeightGrp();
    
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
        for (ptn = 0; ptn < nptn; ptn++) {
            at(i)->ptn_freq[ptn] = patn_freqs[ptn];
        }
    }
}

string IQTreeMix::optimizeModelParameters(bool printInfo, double logl_epsilon) {
    
    size_t i, ptn;
    int step, n, m, substep1, nsubstep1, nsubstep1_start, nsubstep1_max, nsubstep2_start, nsubstep2_max, substep2, nsubstep2, substep2_tot;
    double* pattern_mix_lh;
    double gradient_epsilon = 0.001;
    double gradient2_epsilon = 0.0001;
    double curr_epsilon;
    double prev_score, prev_score1, prev_score2, score, t_score;
    double* prev_ptn_invar;
    bool tree_weight_converge = false;
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
        initializeTreeWeights2();
    }
    
    // initialize the parameters
    optimizeTreesSeparately(printInfo, gradient2_epsilon);
    prev_ptn_invar = ptn_invar;
    ptn_invar = at(0)->ptn_invar;
    
    if (isEdgeLenRestrict) {
        // it is a highly-restricted model in which the edges of different trees having the same partition have the same lengths

        // If there are multiple branches belonging to the same group
        // set all the branches of the same group to their average
        checkBranchGrp();
    }
    
    // show trees
    // cout << "Initial trees:" << endl;
    // showTree();
    // score = computeLikelihood();
    // cout << "Initial likelihood = " << score << endl;

    bool is_ptnfrq_posterior = false;

    nsubstep1 = nsubstep1_start;
    nsubstep2 = nsubstep2_start;
    score = computeLikelihood();
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
                    models[0]->optimizeParameters(gradient2_epsilon);
                    // score = computeLikelihood();
                    // cout << "after optimizing linked subsitution model, likelihood = " << score << endl;
                }

                // optimize the unlinked subsitution models one by one
                if (!isLinkModel) {
                    if (!is_ptnfrq_posterior) {
                        computeFreqArray(pattern_mix_lh, true);
                        is_ptnfrq_posterior = true;
                    }
                    for (i=0; i<models.size(); i++) {
                        models[i]->optimizeParameters(gradient2_epsilon);
                        computeFreqArray(pattern_mix_lh, true, i);
                    }
                    // score = computeLikelihood();
                    // cout << "after optimizing unlinked subsitution model, likelihood = " << score << endl;
                }

                // optimize tree branches for non-branch-length-restricted model
                if (!isEdgeLenRestrict && params->fixed_branch_length != BRLEN_FIX) {
                    if (!is_ptnfrq_posterior) {
                        computeFreqArray(pattern_mix_lh, true);
                        is_ptnfrq_posterior = true;
                    }
                    for (i=0; i<size(); i++) {
                        optimizeAllBranchesOneTree(i, 1, gradient2_epsilon);
                        computeFreqArray(pattern_mix_lh, true, i);
                    }
                    // score = computeLikelihood();
                    // cout << "after optimizing branches by EM, likelihood = " << score << endl;
                }

                // optimize tree branches for branch-length-restricted model
                if (isEdgeLenRestrict && params->fixed_branch_length != BRLEN_FIX) {
                    if (is_ptnfrq_posterior) {
                        resetPtnOrigFreq();
                        is_ptnfrq_posterior = false;
                    }
                    optimizeBranchLensByBFGS(gradient2_epsilon);
                    // score = computeLikelihood();
                    // cout << "after optimizing branches by BFGS, likelihood = " << score << endl;
                }

                score = computeLikelihood();
                if (score < prev_score2 + gradient2_epsilon) {
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
                site_rates[0]->optimizeParameters(gradient2_epsilon);
                if (siterate_names[0].find("R") != string::npos) {
                    site_rates[0]->rescaleRates();
                }
                // score = computeLikelihood();
                // cout << "after optimizing linked site rate model, likelihood = " << score << endl;
            }

            // optimize the unlinked site rate model
            if (anySiteRate && !isLinkSiteRate) {
                if (!is_ptnfrq_posterior) {
                    computeFreqArray(pattern_mix_lh, true);
                    is_ptnfrq_posterior = true;
                }
                for (i=0; i<site_rates.size(); i++) {
                    site_rates[i]->optimizeParameters(gradient2_epsilon);
                    if (siterate_names[i].find("R") != string::npos) {
                        site_rates[i]->rescaleRates();
                    }
                    computeFreqArray(pattern_mix_lh, true, i);
                }
                // score = computeLikelihood();
                // cout << "after optimizing unlinked site-rate model, likelihood = " << score << endl;
            }

            score = computeLikelihood();

            if (score < prev_score1 + gradient2_epsilon) {
                // converged
                break;
            }
            prev_score1 = score;
        }
        
        // optimize tree weights
        if (!isTreeWeightFixed) {
            // if (weightGrpExist || params->optimize_alg_treeweight == "BFGS") {
            if (params->optimize_alg_treeweight == "BFGS") {
                score = optimizeTreeWeightsByBFGS(gradient2_epsilon);
                tree_weight_converge = true;
            } else {
                m = 1 + step / 100;
                score = optimizeTreeWeightsByEM(pattern_mix_lh, gradient2_epsilon, m, tree_weight_converge);  // loop max n times
                // cout << "tree_weight_converge = " << tree_weight_converge << endl;
            }
            // score = computeLikelihood();
            // cout << "after optimizing tree weights, likelihood = " << score << endl;
        }

        if ((step <= 10) || (step<=100 && step%10==0) || (step%50==0)) {
            cout << "step=" << step; // << " substep1=" << substep1 << " substep2_tot=" << substep2_tot;
            cout << " score=" << score;
            /*
            cout << " weights=(";
            for (i=0; i<ntree; i++) {
                if (i>0)
                    cout << ", ";
                cout << weights[i];
            }
            cout << ")";*/
            cout << endl;
        }

        if (score < prev_score + gradient_epsilon) {
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
    // for branch parameters
    if (params->fixed_branch_length != BRLEN_FIX) {
        if (isEdgeLenRestrict) {
            df += branch_group.size();
        } else {
            for (i=0; i<size(); i++) {
                df += at(i)->getNBranchParameters(BRLEN_OPTIMIZE);
            }
        }
    }
    // for tree weights
    if (!isTreeWeightFixed) {
        if (weightGrpExist) {
            df += (weight_group_member.size() - 1);
        } else {
            df += (size() - 1);
        }
    }
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
void IQTreeMix::getPostProb(double* pattern_mix_lh, bool need_computeLike, int update_which_tree) {
    size_t i, ptn, c;
    double* this_lk_cat;
    double lk_ptn;

    if (need_computeLike) {
        computeLikelihood();
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
void IQTreeMix::computeFreqArray(double* pattern_mix_lh, bool need_computeLike, int update_which_tree) {
    size_t i, ptn;
    PhyloTree* tree;

    // get posterior probabilities along each site for each tree
    getPostProb(pattern_mix_lh, need_computeLike, update_which_tree);

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
    double w;
    
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
            for (j=0; j<weight_group_member[i].size(); j++) {
                weights[weight_group_member[i].at(j)] = w;
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
            lower_bound[i+1] = MIN_LEN;
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

// show the log-likelihoods and posterior probabilties for each tree along the sites
void IQTreeMix::showLhProb(ofstream& out) {
    double* pattern_lh_tree;
    double* curr_ptn_lh;
    double* post_prob;
    size_t t,site,idx;
    size_t nsite;
    PhyloTree* ptree;
    double sum;
    int* curr_p_scores;
    int p_score;
    int same_p_score;
    
    if (!parsi_computed) {
        // compute parsimony scores for each tree along the patterns
        // results are stored in the array patn_parsimony
        computeParsimony();
    }

    nsite = aln->getNSite();

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
        at(t)->clearAllPartialLH();
        at(t)->computeLikelihood(curr_ptn_lh);
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
        curr_ptn_lh += nptn;
    }
    
    // for posterior probabilities
    post_prob = new double[ntree];
    
    // print out the log-likelihoods and posterior probabilties for each tree along the sites
    out << "site,log-like,isConstant,isInformative,sameParsimony";
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
        // is the site constant
        if (aln->at(idx).isConst())
            out << "," << 1;
        else
            out << "," << 0;
        // is the site informative
        if (aln->at(idx).isInformative())
            out << "," << 1;
        else
            out << "," << 0;
        
        // does the site have the same parsimony score for all trees
        curr_p_scores = &(patn_parsimony[idx*ntree]);
        p_score = curr_p_scores[0];
        same_p_score = 1;
        if (p_score >= 0) {
            // this is an informative site
            for (t = 1; t < ntree; t++) {
                if (curr_p_scores[t] != p_score) {
                    same_p_score = 0;
                    break;
                }
            }
        }
        out << "," << same_p_score;
        
        curr_ptn_lh = pattern_lh_tree;
        for (t=0; t<ntree; t++) {
            out << "," << curr_ptn_lh[idx]; // log-likelihood of the site for tree t
            curr_ptn_lh += nptn;
        }
        for (t=0; t<ntree; t++) {
            post_prob[t] = post_prob[t] / sum;
            out << "," << post_prob[t]; // posterior probability of the site for tree t
        }
        /*
        // show the characters
        for (t=0; t<ntip; t++) {
            out << "," << aln->convertStateBackStr(aln->at(idx).at(t));
        }
        */
        out << endl;
    }
    
    // free the memory of the array
    delete[] pattern_lh_tree;
    delete[] post_prob;
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
    map<Pattern,int>::iterator itr;
    size_t t,optn,noptn,ptn;

    noptn = aln->ordered_pattern.size();
    UINT* ptn_scores = new UINT[ntree * noptn];
    UINT* curr_ptn_scores;
    int* curr_pars_scores;
    
    deleteAllPartialLh();
    
    // compute the parsimony scores along patterns for each tree
    for (t=0; t<ntree; t++) {
        curr_ptn_scores = ptn_scores + t * noptn;
        at(t)->initCostMatrix(CM_UNIFORM);
        at(t)->setParsimonyKernel(params->SSE);
        at(t)->initializeAllPartialPars();
        at(t)->computeTipPartialParsimony();
        at(t)->computeParsimonyOutOfTreeSankoff(curr_ptn_scores);
    }
    
    // build opattern2id
    for (optn=0; optn<noptn; optn++) {
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

    for (ptn=0; ptn<nptn; ptn++) {
        itr = opattern2id.find(aln->at(ptn));
        curr_pars_scores = &patn_parsimony[ptn*ntree];
        if (itr != opattern2id.end()) {
            optn = itr->second;
            // the parsimony scores
            curr_ptn_scores = ptn_scores;
            for (t=0; t<ntree; t++) {
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
            for (t=0; t<ntree; t++) {
                curr_pars_scores[t] = -1; // parsimony score is not available
            }
        }
    }
    
    parsi_computed = true;

    // free the memory of the array
    delete[] ptn_scores;
}

// show the log-likelihoods and posterior probabilties for each tree along the patterns
void IQTreeMix::showPatternLhProb(ofstream& out) {
    double* pattern_lh_tree;
    double* curr_ptn_lh;
    double* post_prob;
    size_t t,ptn,optn,idx;
    PhyloTree* ptree;
    double sum;
    map<Pattern,int,classcomp> opattern2id;
    map<Pattern,int>::iterator itr;

    size_t noptn = aln->ordered_pattern.size();
    UINT* ptn_scores = new UINT[ntree * noptn];
    UINT* curr_ptn_scores;

    // build opattern2id
    for (optn=0; optn<noptn; optn++) {
        opattern2id.insert(pair<Pattern,int>(aln->ordered_pattern.at(optn),optn));
    }
    
    // compute the parsimony scores along patterns for each tree
    for (t=0; t<ntree; t++) {
        curr_ptn_scores = ptn_scores + t * noptn;
        at(t)->initCostMatrix(CM_UNIFORM);
        at(t)->setParsimonyKernel(params->SSE);
        at(t)->initializeAllPartialPars();
        at(t)->computeTipPartialParsimony();
        at(t)->computeParsimonyOutOfTreeSankoff(curr_ptn_scores);
    }

    // compute likelihood for each tree
    pattern_lh_tree = new double[nptn * ntree];
    curr_ptn_lh = pattern_lh_tree;
    for (t=0; t<ntree; t++) {
        // save the site rate's tree
        ptree = at(t)->getRate()->getTree();
        // set the tree t as the site rate's tree
        // and compute the likelihood values
        at(t)->getRate()->setTree(at(t));
        at(t)->clearAllPartialLH();
        at(t)->computeLikelihood(curr_ptn_lh);
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
        curr_ptn_lh += nptn;
    }
    
    // for posterior probabilities
    post_prob = new double[ntree];
    
    // print out the log-likelihoods and posterior probabilties for each tree along the sites
    out << "pattern-id,is-informative,freq,log-like";
    for (t=0; t<ntree; t++) {
        out << ",parsimony tree " << t+1;
    }
    for (t=0; t<ntree; t++) {
        out << ",log-like tree " << t+1;
    }
    for (t=0; t<ntree; t++) {
        out << ",post-prob tree " << t+1;
    }
    // print out the sequence names
    for (t=0; t<ntip; t++) {
        out << "," << aln->getSeqName(t);
    }
    out << endl;
    
    for (ptn=0; ptn<nptn; ptn++) {
        
        out << ptn+1;
        if (aln->at(ptn).isInformative())
            out << "," << 1;
        else
            out << "," << 0;
        out << "," << patn_freqs[ptn];
        
        curr_ptn_lh = pattern_lh_tree;
        sum = 0.0;
        for (t=0; t<ntree; t++) {
            post_prob[t] = exp(curr_ptn_lh[ptn]) * weights[t];
            sum += post_prob[t];
            curr_ptn_lh += nptn;
        }
        out << "," << log(sum); // log-likelihood of the site
        
        itr = opattern2id.find(aln->at(ptn));
        if (itr != opattern2id.end()) {
            optn = itr->second;
            // show the parsimony scores
            curr_ptn_scores = ptn_scores;
            for (t=0; t<ntree; t++) {
                out << "," << curr_ptn_scores[optn];
                curr_ptn_scores += noptn;
            }
        } else {
            for (t=0; t<ntree; t++) {
                out << ",--";
            }
        }
        curr_ptn_lh = pattern_lh_tree;
        for (t=0; t<ntree; t++) {
            out << "," << curr_ptn_lh[ptn]; // log-likelihood of the pattern for tree t
            curr_ptn_lh += nptn;
        }
        for (t=0; t<ntree; t++) {
            post_prob[t] = post_prob[t] / sum;
            out << "," << post_prob[t]; // posterior probability of the site for tree t
        }
        for (t=0; t<ntip; t++) {
            out << "," << aln->convertStateBackStr(aln->at(ptn).at(t));
        }
        out << endl;
    }

    // free the memory of the array
    delete[] pattern_lh_tree;
    delete[] post_prob;
    delete[] ptn_scores;
}

// show the log-likelihoods and posterior probabilties for each tree along the patterns
void IQTreeMix::showOrderedPatternLhProb(ofstream& out) {
    double* pattern_lh_tree;
    double* curr_ptn_lh;
    double* post_prob;
    size_t t,ptn,optn,idx;
    PhyloTree* ptree;
    double sum;
    map<Pattern,int,classcomp> pattern2id;
    map<Pattern,int>::iterator itr;

    size_t noptn = aln->ordered_pattern.size();
    UINT* ptn_scores = new UINT[ntree * noptn];
    UINT* curr_ptn_scores;
    
    // compute likelihood for each tree
    pattern_lh_tree = new double[nptn * ntree];
    curr_ptn_lh = pattern_lh_tree;
    for (t=0; t<ntree; t++) {
        // save the site rate's tree
        ptree = at(t)->getRate()->getTree();
        // set the tree t as the site rate's tree
        // and compute the likelihood values
        at(t)->getRate()->setTree(at(t));
        at(t)->clearAllPartialLH();
        at(t)->computeLikelihood(curr_ptn_lh);
        // set back the prevoius site rate's tree
        at(t)->getRate()->setTree(ptree);
        curr_ptn_lh += nptn;
    }
    
    // build pattern2id
    ASSERT(aln->size() == nptn);
    for (ptn=0; ptn<nptn; ptn++) {
        pattern2id.insert(pair<Pattern,int>(aln->at(ptn),ptn));
    }
    
    // compute the parsimony scores along patterns for each tree
    for (t=0; t<ntree; t++) {
        curr_ptn_scores = ptn_scores + t * noptn;
        at(t)->initCostMatrix(CM_UNIFORM);
        at(t)->setParsimonyKernel(params->SSE);
        at(t)->initializeAllPartialPars();
        at(t)->computeTipPartialParsimony();
        at(t)->computeParsimonyOutOfTreeSankoff(curr_ptn_scores);
    }

    // for posterior probabilities
    post_prob = new double[ntree];
    
    // print out the log-likelihoods and posterior probabilties for each tree along the sites
    out << "ordered-pattern-id,is-informative,freq optn, freq ptn,log-like";
    for (t=0; t<ntree; t++) {
        out << ",parsimony tree " << t+1;
    }
    for (t=0; t<ntree; t++) {
        out << ",log-like tree " << t+1;
    }
    for (t=0; t<ntree; t++) {
        out << ",post-prob tree " << t+1;
    }
    // print out the sequence names
    for (t=0; t<ntip; t++) {
        out << "," << aln->getSeqName(t);
    }
    // print out the sequence names
    for (t=0; t<ntip; t++) {
        out << ",O" << aln->getSeqName(t);
    }
    out << endl;
    
    for (optn=0; optn<noptn; optn++) {
        if (aln->ordered_pattern[optn].frequency == 0)
            continue;
        out << optn+1;
        itr = pattern2id.find(aln->ordered_pattern[optn]);
        if (itr == pattern2id.end()) {
            // show nucleotide (for ordered pattern)
            cout << "The " << optn << "-th ordered pattern with frequency " << aln->ordered_pattern[optn].frequency << " ";
            for (t=0; t<ntip; t++) {
                cout << aln->convertStateBackStr(aln->ordered_pattern[optn].at(t));
            }
            cout << " cannot be found!" << endl;
        }
        ASSERT(itr != pattern2id.end());
        ptn = itr->second;
        if (aln->at(ptn).isInformative())
            out << "," << 1;
        else
            out << "," << 0;
        out << "," << aln->ordered_pattern[optn].frequency;
        out << "," << patn_freqs[ptn];
        curr_ptn_lh = pattern_lh_tree;
        sum = 0.0;
        for (t=0; t<ntree; t++) {
            post_prob[t] = exp(curr_ptn_lh[ptn]) * weights[t];
            sum += post_prob[t];
            curr_ptn_lh += nptn;
        }
        out << "," << log(sum); // log-likelihood of the site

        // show the parsimony scores
        curr_ptn_scores = ptn_scores;
        for (t=0; t<ntree; t++) {
            out << "," << curr_ptn_scores[optn];
            curr_ptn_scores += noptn;
        }
        
        // show log-likelihood values
        curr_ptn_lh = pattern_lh_tree;
        for (t=0; t<ntree; t++) {
            out << "," << curr_ptn_lh[ptn];
            curr_ptn_lh += nptn;
        }
        
        // show posterior probabilties
        for (t=0; t<ntree; t++) {
            post_prob[t] = post_prob[t] / sum;
            out << "," << post_prob[t];
        }
        
        // show nucleotide
        for (t=0; t<ntip; t++) {
            out << "," << aln->convertStateBackStr(aln->at(ptn).at(t));
        }

        // show nucleotide (for ordered pattern)
        for (t=0; t<ntip; t++) {
            out << "," << aln->convertStateBackStr(aln->ordered_pattern.at(optn).at(t));
        }

        out << endl;
    }

    // free the memory of the array
    delete[] pattern_lh_tree;
    delete[] post_prob;
    delete[] ptn_scores;
}
