/*
 *  alisim.h
 *  implemetation of AliSim (Alignment Simulator)
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "alisim.h"

void runAliSim(Params &params, Checkpoint *checkpoint)
{
    // inferring input parameters if inference mode is active
    if (params.alisim_inference)
    {
        inferInputParameters(params, checkpoint);
    }
    // otherwise, generate a random tree if the random mode is active
    else if (params.tree_gen != NONE)
        {
            generateRandomTree(params);
        }
    
    // run AliSim
    runAliSimWithoutInference(params);
}

/**
*  inferring input parameters for AliSim
*/
void inferInputParameters(Params &params, Checkpoint *checkpoint)
{
    // check if a tree_file is provided or not
    if (params.user_file)
    {
        // case 0: skip inferring step if all input parameters for AliSim are provided.
        if ((params.original_params.find("-m ") != std::string::npos)
            && (params.original_params.find("-m TEST") == std::string::npos)
            && (params.original_params.find("-m MF") == std::string::npos))
        {
            outWarning("Inference mode is deactive since all input parameters for AliSim are provided.");
            return;
        }
        // case 1: only a tree is supplied
        else
        {
            // run ModelFinder to infer the model
            params.model_name = "MF";
            runPhyloAnalysis(Params::getInstance(), checkpoint);
            
            // initialize dummy variables
            int sequence_length = 1000;
            string model = "";
            
            // extract sequence_length and model's parameters from the output of ModelFinder
            char *iqtree_file_path = new char[strlen(params.aln_file) + 8];
            strcpy(iqtree_file_path, params.aln_file);
            strcat(iqtree_file_path,".iqtree");
            extractInputParameters(iqtree_file_path, sequence_length, model, true);
            puts (iqtree_file_path);
            
            // transfer inferred parameters to AliSim's inputs
            params.model_name = model;
            if (params.original_params.find("--length") == std::string::npos)
                params.alisim_sequence_length = sequence_length;
        }
    }
    else
    {
        // case 2: only a model is supplied
        if ((params.original_params.find("-m ") != std::string::npos)
            && (params.original_params.find("-m TEST") == std::string::npos)
            && (params.original_params.find("-m MF") == std::string::npos))
        {
            // run IQTree to infer the tree (provided the model)
            runPhyloAnalysis(Params::getInstance(), checkpoint);
            
            // initialize the tree_file_path
            params.user_file = new char[strlen(params.aln_file) + 10];
            strcpy(params.user_file, params.aln_file);
            strcat(params.user_file,".treefile");
            
            // extract the sequence_length for AliSim from the output of IQTree (if it's yet been specified by the user)
            if (params.original_params.find("--length") == std::string::npos)
            {
                int sequence_length = 1000;
                string model = "";
                
                char *iqtree_file_path = new char[strlen(params.aln_file) + 8];
                strcpy(iqtree_file_path, params.aln_file);
                strcat(iqtree_file_path,".iqtree");
                extractInputParameters(iqtree_file_path, sequence_length, model, false);
                puts (iqtree_file_path);
                
                params.alisim_sequence_length = sequence_length;
            }
        }
        // case 3: neither a tree nor a model is supplied
        else
        {
            // run IQTree to infer the tree and the model
            //params.model_name = "MF";
            runPhyloAnalysis(Params::getInstance(), checkpoint);
            
            // initialize the tree_file_path
            params.user_file = new char[strlen(params.aln_file) + 10];
            strcpy(params.user_file, params.aln_file);
            strcat(params.user_file,".treefile");
            
            // initialize dummy variables
            int sequence_length = 1000;
            string model = "";
            
            // extract sequence_length and model's parameters from the output of IQTree
            char *iqtree_file_path = new char[strlen(params.aln_file) + 8];
            strcpy(iqtree_file_path, params.aln_file);
            strcat(iqtree_file_path,".iqtree");
            extractInputParameters(iqtree_file_path, sequence_length, model, true);
            puts (iqtree_file_path);
            
            // transfer inferred parameters to AliSim's inputs
            params.model_name = model;
            if (params.original_params.find("--length") == std::string::npos)
            {
                params.alisim_sequence_length = sequence_length;
            }
        }
    }
    
    cout << "Start AliSim with Input parameters inferred from IQTree: " << endl;
    cout << "- Treefile: " << params.user_file << endl;
    cout << "- Model: " << params.model_name << endl;
}

/**
*  generate a random tree
*/
void generateRandomTree(Params &params)
{
    if (params.sub_size < 3 && !params.aln_file) {
        outError(ERR_FEW_TAXA);
    }

    if (!params.user_file) {
        outError("Please specify an output tree file name");
    }
    ////cout << "Random number seed: " << params.ran_seed << endl << endl;

    SplitGraph sg;

    try {

        if (params.tree_gen == YULE_HARDING || params.tree_gen == CATERPILLAR ||
            params.tree_gen == BALANCED || params.tree_gen == UNIFORM || params.tree_gen == STAR_TREE || params.tree_gen == BIRTH_DEATH) {
            if (!overwriteFile(params.user_file)) return;
            ofstream out;
            out.open(params.user_file);
            MTree itree;

            if (params.second_tree) {
                cout << "Generating random branch lengths on tree " << params.second_tree << " ..." << endl;
                itree.readTree(params.second_tree, params.is_rooted);
            } else
            switch (params.tree_gen) {
            case YULE_HARDING:
                cout << "Generating random Yule-Harding tree..." << endl;
                break;
            case UNIFORM:
                cout << "Generating random uniform tree..." << endl;
                break;
            case CATERPILLAR:
                cout << "Generating random caterpillar tree..." << endl;
                break;
            case BALANCED:
                cout << "Generating random balanced tree..." << endl;
                break;
            case STAR_TREE:
                cout << "Generating star tree with random external branch lengths..." << endl;
                break;
            case BIRTH_DEATH:
                cout << "Generating random Birth-death tree..." << endl;
                break;
            default: break;
            }
            ofstream out2;
            if (params.num_zero_len) {
                cout << "Setting " << params.num_zero_len << " internal branches to zero length..." << endl;
                string str = params.user_file;
                str += ".collapsed";
                out2.open(str.c_str());
            }
            for (int i = 0; i < params.repeated_time; i++) {
                MExtTree mtree;
                if (itree.root) {
                    mtree.copyTree(&itree);
                    mtree.generateRandomBranchLengths(params);
                } else {
                    mtree.generateRandomTree(params.tree_gen, params);
                }
                if (params.num_zero_len) {
                    mtree.setZeroInternalBranches(params.num_zero_len);
                    MExtTree collapsed_tree;
                    collapsed_tree.copyTree(&mtree);
                    collapsed_tree.collapseZeroBranches();
                    collapsed_tree.printTree(out2);
                    out2 << endl;
                }
                mtree.printTree(out);
                out << endl;
            }
            out.close();
            cout << params.repeated_time << " tree(s) printed to " << params.user_file << endl;
            if (params.num_zero_len) {
                out2.close();
                cout << params.repeated_time << " collapsed tree(s) printed to " << params.user_file << ".collapsed" << endl;
            }
        }
        // Generate random trees if optioned
        else if (params.tree_gen == CIRCULAR_SPLIT_GRAPH) {
            cout << "Generating random circular split network..." << endl;
            if (!overwriteFile(params.user_file)) return;
            sg.generateCircular(params);
        } else if (params.tree_gen == TAXA_SET) {
            sg.init(params);
            cout << "Generating random taxa set of size " << params.sub_size <<
                " overlap " << params.overlap << " with " << params.repeated_time << " times..." << endl;
            if (!overwriteFile(params.pdtaxa_file)) return;
            sg.generateTaxaSet(params.pdtaxa_file, params.sub_size, params.overlap, params.repeated_time);
        }
    } catch (bad_alloc) {
        outError(ERR_NO_MEMORY);
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, params.user_file);
    }

    // calculate the distance
    if (params.run_mode == CALC_DIST) {
        if (params.tree_gen == CIRCULAR_SPLIT_GRAPH) {
            cout << "Calculating distance matrix..." << endl;
            sg.calcDistance(params.dist_file);
            cout << "Distances printed to " << params.dist_file << endl;
        }// else {
            //mtree.calcDist(params.dist_file);
        //}
    }

}

/**
*  extract input parameters for AliSim from output (*.iqtree) file after inferring
*/
void extractInputParameters(char *iqtree_file_path, int &sequence_length, string &model, bool extract_model){
    
    ostringstream err_str;
    igzstream in;
    int line_num = 1;
    string line;
    
    bool rate_parameter_begin = false;
    bool state_freqs_begin = false;
    DoubleVector state_freqs;
    DoubleVector rate_params;
    bool rate_equal = true;
    double invariant_prop = 0;
    double alpha_shape = 0;
    int gamma_num_categories = 0;
    int freerate_num_categories = 0;
    DoubleVector freerate_props;
    DoubleVector freerate_rates;
    string model_name = "";
    
    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(iqtree_file_path);
    
    // remove the failbit
    in.exceptions(ios::badbit);

    for (; !in.eof(); line_num++) {
        safeGetline(in, line);
        line = line.substr(0, line.find_first_of("\n\r"));
        if (line == "") continue;
        
        // extract sequence_length
        if (regex_match(line, regex("(Input data: )(.*)")))
        {
            string delimiter = "sequences with ";
            line.erase(0, line.find(delimiter) + delimiter.length());
            delimiter = " ";
            sequence_length = convert_int(line.substr(0, line.find(delimiter)).c_str());
            if (!extract_model)
                break;
            else
                continue;
        }
        
        // extract state_freqs
        if (regex_match(line, regex("(Rate matrix Q:)")))
        {
            state_freqs_begin = false;
            continue;
        }
        if (state_freqs_begin)
        {
            string delimiter = " = ";
            double state_freq = convert_double(line.substr(delimiter.length() + line.find(delimiter), line.length() - (delimiter.length() + line.find(delimiter))).c_str());
            state_freqs.resize(state_freqs.size() + 1);
            state_freqs[state_freqs.size() - 1] = state_freq;
            continue;
        }
        if (regex_match(line, regex("(State frequencies:)(.*)")))
        {
            rate_parameter_begin = false;
            // normalize rate_params and remove the last param
            if (rate_params.size()>0)
            {
                for (int i = 0; i < rate_params.size() - 1; i ++)
                {
                    rate_params[i] = rate_params[i]/rate_params[rate_params.size()-1];
                    if (rate_equal && rate_params[i] != 1.0)
                        rate_equal = false;
                }
                
                rate_params.resize(rate_params.size() - 1);
            }
            state_freqs_begin = true;
            continue;
        }
        
        // extract rate_params
        if (rate_parameter_begin)
        {
            string delimiter = ": ";
            double rate_param = convert_double(line.substr(delimiter.length() + line.find(delimiter), line.length() - (delimiter.length() + line.find(delimiter))).c_str());
            rate_params.resize(rate_params.size() + 1);
            rate_params[rate_params.size() - 1] = rate_param;
            continue;
        }
        if (regex_match(line, regex("(Rate parameter R:)(.*)")))
        {
            rate_parameter_begin = true;
            continue;
        }
        
        // extract model_name
        if (regex_match(line, regex("(Model of substitution: )(.*)")))
        {
            string delimiter = "Model of substitution: ";
            model_name = line.substr(delimiter.length(), line.length() - delimiter.length());
            continue;
        }
        
        // extract invariant_prop
        if (regex_match(line, regex("(Proportion of invariable sites: )(.*)")))
        {
            string delimiter = "Proportion of invariable sites: ";
            invariant_prop = convert_double(line.substr(delimiter.length(), line.length() - delimiter.length()).c_str());
            continue;
        }
        
        // extract alpha_shape
        if (regex_match(line, regex("(Gamma shape alpha: )(.*)")))
        {
            string delimiter = "Gamma shape alpha: ";
            alpha_shape = convert_double(line.substr(delimiter.length(), line.length() - delimiter.length()).c_str());
            continue;
        }
        
        // extract gamma_num_categories
        if (regex_match(line, regex("(.*)(Gamma with )(.*)")))
        {
            string delimiter = "Gamma with ";
            gamma_num_categories = convert_int(line.substr(delimiter.length() + line.find(delimiter), line.length() - (delimiter.length() + line.find(delimiter) + 11)).c_str());
            continue;
        }
        
        // extract freerate_num_categories
        if (regex_match(line, regex("(.*)(FreeRate with )(.*)")))
        {
            string delimiter = "FreeRate with ";
            freerate_num_categories = convert_int(line.substr(delimiter.length() + line.find(delimiter), line.length() - (delimiter.length() + line.find(delimiter) + 11)).c_str());
            continue;
        }
        
        // extract freerate_props and freerate_rates
        if (regex_match(line, regex("(Site proportion and rates: )(.*)")))
        {
            string delimiter = "Site proportion and rates: ";
            line = line.substr(delimiter.length(), line.length()-delimiter.length());
            while (line.length() > 0)
            {
                string prop_rate = line.substr(line.find("("), line.find(")"));
                freerate_props.resize(freerate_props.size() + 1);
                freerate_rates.resize(freerate_rates.size() + 1);
                freerate_props[freerate_props.size() - 1] = convert_double(prop_rate.substr(1, prop_rate.find(",") - 1).c_str());
                freerate_rates[freerate_rates.size() - 1] = convert_double(prop_rate.substr(prop_rate.find(",") + 1, prop_rate.length() - prop_rate.find(",") - 2).c_str());
                line.erase(0, line.find(")")+1);
            }
            continue;
        }
        
    }
    
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();
    
    // reconstruct the model_fullname
    model = model_name.substr(0, model_name.find("+"));
    // add rate_params
    if (!rate_equal)
    {
        model = model + "{";
        for (int i = 0; i < rate_params.size(); i++)
            if (i == rate_params.size() - 1)
                model = model + convertDoubleToString(rate_params[i]);
            else
                model = model + convertDoubleToString(rate_params[i])+",";
        model = model + "}";
    }
    // add +F & its paramters (if any)
    if (model_name.find("+F") != std::string::npos)
    {
        model = model + "+F{";
        for (int i = 0; i < state_freqs.size(); i++)
            if (i == state_freqs.size() - 1)
                model = model + convertDoubleToString(state_freqs[i]);
            else
                model = model + convertDoubleToString(state_freqs[i])+",";
        model = model + "}";
    }
    // add +I & its paramter (if any)
    if (model_name.find("+I") != std::string::npos)
    {
        model = model + "+I{" + convertDoubleToString(invariant_prop) + "}";
    }
    // add +G & its paramters (if any)
    if (model_name.find("+G") != std::string::npos)
    {
        model = model + "+G" + convertIntToString(gamma_num_categories) + "{" + convertDoubleToString(alpha_shape) + "}";
    }
    // add +R & its paramters (if any)
    if (model_name.find("+R") != std::string::npos)
    {
        model = model + "+R" + convertIntToString(freerate_num_categories) + "{";
        for (int i = 0; i < freerate_num_categories; i++)
            if (i == freerate_num_categories - 1)
                model = model + convertDoubleToString(freerate_props[i]) + "," + convertDoubleToString(freerate_rates[i]);
            else
                model = model + convertDoubleToString(freerate_props[i]) + "," + convertDoubleToString(freerate_rates[i]) + ",";
        convertDoubleToString(alpha_shape);
        model = model + "}";
    }
}

/**
*  execute AliSim without inference
*/
void runAliSimWithoutInference(Params params)
{
    cout << "[Alignment Simulator] Executing" <<"\n";
    
    // case 1 (default): without rate heterogeneity
    AliSimulator *alisimulator = new AliSimulator(&params);
    
    // get variables
    string rate_name = alisimulator->tree->getRateName();
    double invariant_proportion = alisimulator->tree->getRate()->getPInvar();
    
    // case 2: with rate heterogeneity
    if (!rate_name.empty())
    {
        if((rate_name.find("+G") != std::string::npos) || (rate_name.find("+R") != std::string::npos))
        {
            // case 2.1: with rate heterogeneity (gamma/freerate model with invariant sites)
            if (invariant_proportion > 0)
            {
                alisimulator = new AliSimulatorHeterogeneityInvar(alisimulator, invariant_proportion);
            }
            // case 2.2: with rate heterogeneity (gamma/freerate model without invariant sites)
            else
            {
                alisimulator = new AliSimulatorHeterogeneity(alisimulator);
                
            }
        }
        // case 2.3: without gamma/freerate model with only invariant sites
        else if (rate_name.find("+I") != std::string::npos)
        {
            alisimulator = new AliSimulatorInvar(alisimulator, invariant_proportion);
        }
    }
    
    // show parameters
    alisimulator->showParameters();
    
    // iteratively generate multiple/a single  alignment(s) for each tree
    alisimulator->generateMultipleAlignmentsFromSingleTree();
    
    cout << "[Alignment Simulator] Done"<<"\n";
    
    // delete alisimulator
    delete alisimulator;
}
