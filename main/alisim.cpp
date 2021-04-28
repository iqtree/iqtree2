/*
 *  alisim.h
 *  implemetation of AliSim (Alignment Simulator)
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "alisim.h"

void runAliSim(Params &params, Checkpoint *checkpoint)
{
    // Init variables
    IQTree *tree;
    Alignment *aln;
    
    // generate a random tree if neccessary
    if (params.tree_gen != NONE)
    {
        // draw a random num_taxa from a user-defined list
        if (!params.alisim_num_taxa_list.empty())
        {
            int rand_index = random_int(params.alisim_num_taxa_list.size());
            params.sub_size = params.alisim_num_taxa_list.at(rand_index);
        }
        // draw a random num_taxa from uniform distribution
        else if (params.alisim_num_taxa_uniform_start > 3)
        {
            int range = params.alisim_num_taxa_uniform_end - params.alisim_num_taxa_uniform_start + 1;
            params.sub_size = params.alisim_num_taxa_uniform_start + random_int(range);
        }
            
        generateRandomTree(params);
        
        // reset flags to make sure IQTree will not re-generate new random starting_tree
        params.start_tree = STT_PLL_PARSIMONY;
        params.tree_gen = NONE;
    }
    
    // inferring input parameters if inference mode is active
    if (params.aln_file)
    {
        inferInputParameters(params, checkpoint, tree, aln);
    }
    
    // run AliSim without inference
    runAliSimWithoutInference(params, tree);
    
    // aln and tree are deleted in distructor of AliSimSimulator
}

/**
*  inferring input parameters for AliSim
*/
void inferInputParameters(Params &params, Checkpoint *checkpoint, IQTree *&tree, Alignment *&aln)
{
    // if user has not specified model_name -> set model_name = "" (to override the default model JC)
    if ((params.original_params.find("-m ") == std::string::npos)
        || (params.original_params.find("-m TEST") != std::string::npos)
        || (params.original_params.find("-m MF") != std::string::npos))
    {
        params.model_name = "";
    }
    
    // infer model's parameters [and a tree]
    runPhyloAnalysis(Params::getInstance(), checkpoint, tree, aln);
    ASSERT(tree && tree->getModel() && tree->aln);
    
    // update model_name
    params.model_name = tree->getModel()->getNameParams();
    
    // initialize the tree_file if it has not been provided by the user
    if (!params.user_file)
    {
        // initialize the tree_file if it has not been provided by the user
        params.user_file = new char[strlen(params.aln_file) + 10];
        strcpy(params.user_file, params.aln_file);
        strcat(params.user_file,".treefile");
    }
    // reload tree from tree_file
    else
    {
        // do not reload the tree from file in case with heterotachy
        if (!tree->getRate()->isHeterotachy())
        {
            bool is_rooted = false;
            tree->readTree(params.user_file, is_rooted);
        }
    }
    
    // update sequence_length
    if (params.original_params.find("--length") == std::string::npos)
        params.alisim_sequence_length = tree->aln->getNSite() * params.alisim_sites_per_state;
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
*  execute AliSim without inference
*/
void runAliSimWithoutInference(Params params, IQTree *&tree)
{
    cout << "[Alignment Simulator] Executing" <<"\n";
    
    // case 1 (default): without rate heterogeneity
    AliSimulator *alisimulator;
    if (tree && params.aln_file)
        alisimulator = new AliSimulator(&params, tree);
    else
        alisimulator = new AliSimulator(&params);
    
    // get variables
    string rate_name = alisimulator->tree->getRateName();
    double invariant_proportion = alisimulator->tree->getRate()->getPInvar();
    bool is_mixture_model = alisimulator->tree->getModel()->isMixture();
    
    // case 2: with rate heterogeneity or mixture model
    if ((!rate_name.empty()) || is_mixture_model)
    {
        // case 2.3: with only invariant sites (without gamma/freerate model/mixture models)
        if (!rate_name.compare("+I"))
        {
            alisimulator = new AliSimulatorInvar(alisimulator, invariant_proportion);
        }
        else
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
        
    }
    
    // show parameters
    alisimulator->showParameters();
    
    // iteratively generate multiple/a single  alignment(s) for each tree
    alisimulator->generateMultipleAlignmentsFromSingleTree();
    
    cout << "[Alignment Simulator] Done"<<"\n";
    
    // delete alisimulator
    delete alisimulator;
}
