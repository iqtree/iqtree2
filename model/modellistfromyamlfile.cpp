//
// modellistfromyamlfile.cpp
// 
/***************************************************************************
 *   Created by James Barbetti on 30-Apr-2021                              *
 *   james_barbetti@yahoo.com                                              *
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

#include "model/modelinfo.h"
#include "modelinfofromyamlfile.h"
#include "modelsubst.h"            //for OPEN_BRACKET and CLOSE_BRACKET
#include "modelfileloader.h"
#include "modelexpression.h"       //for InterpretedExpression
#include "yamlmodelwrapper.h"      //for YAMLModelWrapper template class

#include <utils/my_assert.h>       //for ASSERT macro
#include <utils/stringfunctions.h> //for convert_int
#include <utils/tools.h>           //for outError

namespace {
    const string dummy_rate_params;
    const string dummy_freq_params;    
};

void ModelListFromYAMLFile::loadFromFile (Params& params,
                                          const char* file_path,
                                          LoggingTarget* report_to_tree) {
    try {
        YAML::Node yaml_model_list = YAML::LoadFile(file_path);
        ModelFileLoader loader(file_path);
        if (!yaml_model_list.IsSequence()) {
            throw YAML::Exception(yaml_model_list.Mark(),
                                  "list '[...]' expected");
        }
        for (auto node : yaml_model_list) {
            if (!(node["substitutionmodel"])) {
                if (node["ratemodel"]) {
                    std::string rate_model_name = node["ratemodel"].Scalar();
                    TREE_LOG_LINE(*report_to_tree, YAMLParsingVerbosity,
                                "Parsing YAML rate model " << rate_model_name);
                    if (rate_models_found.hasName(rate_model_name)) {
                        std::stringstream complaint;
                        complaint << "Duplicate definition of rate model "
                                  << rate_model_name <<".";
                        std::string complaint_string = complaint.str();
                        throw YAML::Exception(yaml_model_list.Mark(),
                                              complaint_string.c_str());
                    }
                    ModelInfoFromYAMLFile* rate = rate_models_found.insertNew(rate_model_name);
                    rate->is_rate_model          = true;
                    loader.parseYAMLModel(params, node, rate_model_name, *rate, *this,
                                          nullptr, report_to_tree);
                }
                continue;
            }
            std::string yaml_model_name = node["substitutionmodel"].Scalar();
            TREE_LOG_LINE(*report_to_tree, YAMLParsingVerbosity,
                          "Parsing YAML substitution model " << yaml_model_name);
            ModelInfoFromYAMLFile* model = models_found.insertNew(yaml_model_name);
            loader.parseYAMLModel(params, node, yaml_model_name, *model, 
                                  *this, nullptr, report_to_tree);
        }
    }
    catch (YAML::ParserException& e) {
        outError(e.what());
    }
    catch (YAML::Exception& e) {
        outError(e.what());
    }
    catch (ModelExpression::ModelException& x) {
        outError(x.getMessage());
    }
}

bool ModelListFromYAMLFile::isSubstitutionModelNameRecognized 
        (const char* model_name) {
    size_t i = 0;
    while (model_name[i]!='\0' && model_name[i]!='{') {
        ++i;
    }
    std::string model_front = std::string(model_name, i);

    auto found      = models_found.hasName(model_front);
    return found;
}

bool ModelListFromYAMLFile::isSubstitutionModelNameRecognized 
        (const std::string& model_name) {
    return isSubstitutionModelNameRecognized(model_name.c_str());
}

bool ModelListFromYAMLFile::hasModel(const std::string& model_name) const {
    return models_found.hasName(model_name) 
        || rate_models_found.hasName(model_name);
}

StrVector ModelListFromYAMLFile::getSubstitutionModelNames() const {
    StrVector answer;
    for (const std::string& name : models_found.getNames()) {
        answer.push_back(name);
    }
    return answer;
}

std::string ModelListFromYAMLFile::getListOfSubstitutionModelNames() const {
    std::stringstream answer;
    const char* separator = "";
    for (std::string name : models_found.getNames()) {
        answer << separator << name;
        separator = ", ";
    }
    return answer.str();
}

std::string ModelListFromYAMLFile::getListOfRateModelNames() const {
    std::stringstream answer;
    const char* separator = "";
    for (std::string name :rate_models_found.getNames()) {
        answer << separator << name;
        separator = ", ";
    }
    return answer.str();
}

const ModelInfoFromYAMLFile*
    ModelListFromYAMLFile::getModel(const std::string& model_name) const {
    const ModelInfoFromYAMLFile* info = models_found.getPointer(model_name);
    if (info!=nullptr) {
        return info;
    }
    info = rate_models_found.getPointer(model_name);
    return info;
}

namespace {
    void extractModelNameAndParameters(const char* model_plus_params,
                                       std::string& name,
                                       std::string& params) {
        size_t i = 0;
        while (model_plus_params[i]!='\0' && model_plus_params[i]!='{') {
            ++i;
        }
        name   = std::string(model_plus_params,i);
        params = std::string(model_plus_params+i);
        //includes the { at front and the } at back.
    }
};

bool ModelListFromYAMLFile::isRateHeterotachyRequired
        (Params& params, 
         const std::string& model_name,
         int &num_mixlen,
         LoggingTarget* logging_target) {
    std::string name;
    std::string inherit_list;    //list of classes to inherit
    std::string parameter_list;  //list of parameters to pass to model
    extractModelNameAndParameters(model_name.c_str(), name, parameter_list);
    ModelInfoFromYAMLFile* model_info = models_found.getPointer(name);
    ModelInfoFromYAMLFile  dynamic_model;
    if (model_info==nullptr) {
        //std::cout << "XX did not find substitution model pointer"
        //          << " for " << name << std::endl;
        std::string base_model = split_string(name, "+")[0];
        if (base_model != name) {
            auto base_model_info = models_found.getPointer(base_model);
            inherit_list  = name.substr(base_model.length()+1);
            //std::cout << "XX base model is " << base_model << std::endl;
            //std::cout << "XX inherit list is " << inherit_list << std::endl;

            if (base_model_info!=nullptr) {
                dynamic_model = *base_model_info;
                model_info    = &dynamic_model;
                ModelFileLoader loader(params.yaml_model_file.c_str());
                loader.handleInheritance(params, dynamic_model, *this,
                                         inherit_list, false, logging_target);
            }
        } 
        if (model_info==nullptr) {
            return false;
        }
    }    
    dynamic_model = *model_info;
    model_info    = &dynamic_model;
    if (!parameter_list.empty()) {
        dynamic_model.acceptParameterList(params, parameter_list, logging_target);
    }
    auto rate_info = model_info->specified_rate_model_info;
    
    if (rate_info==nullptr) {
        //std::cout << "XX did not find rate info for " << name << std::endl;
        return false;
    }
    if (!rate_info->hasRateHeterotachy()) {
        //std::cout << "XX Rate model wasn't heterotachic, for " << name << std::endl;
        return false;
    }
    num_mixlen = rate_info->getNumberOfProportions()
               - rate_info->getNumberOfInvariantProportions();
    //Problem.  This includes the invariant proportion

    //std::cout << "XX Number of heterotarchic proportions"
    //          << " from " << rate_info->getName()
    //          << " was " << num_mixlen << std::endl;
    return true;
}

ModelMarkov* ModelListFromYAMLFile::getModelByName
            (const char* model_name,   PhyloTree*    tree,
             const char* model_params, StateFreqType freq_type,
             const char* freq_params,  ModelsBlock*  models_block,
             PhyloTree* report_to_tree) {
    std::string name;
    std::string parameter_list;
    extractModelNameAndParameters(model_name, name, parameter_list);
    ModelInfoFromYAMLFile* model_info = models_found.getPointer(model_name);
    ASSERT(model_info!=nullptr);
    if (freq_type == StateFreqType::FREQ_UNKNOWN) {
        freq_type = model_info->frequency_type;
    }
    if (0<strlen(model_params) || 0<strlen(freq_params)) {
        TREE_LOG_LINE(*report_to_tree, YAMLFrequencyVerbosity,
                      "Model Params: " << model_params
                      << " Freq Params: " << freq_params);
        if (parameter_list.empty()) {
            parameter_list = model_params;
        }
    }
    ModelMarkov* model = nullptr;
    try {
    return getModelByReference(*model_info, tree, freq_type, models_block, 
                               parameter_list, report_to_tree);
    } catch (ModelExpression::ModelException& x) {
        std::stringstream complaint;
        complaint << "Cannot initialize model " 
                  << model_name << ": " << x.getMessage();
        outError(complaint.str());
    }
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getModelByReference
                (ModelInfoFromYAMLFile& model_info, PhyloTree*   tree,
                 StateFreqType freq_type,           ModelsBlock* models_block,
                 const std::string &parameter_list, 
                 PhyloTree* report_to_tree) {
    if (model_info.isDivergentModel()) {
        return getDivergentModel(model_info, parameter_list,
                               freq_type, models_block, tree,  
                               report_to_tree);
    }
    if (model_info.isMixtureModel()) {
        return getMixtureModel(model_info, parameter_list,
                               freq_type, models_block, tree,  
                               report_to_tree);
    }
    switch (model_info.sequence_type) {
        case SeqType::SEQ_BINARY:
            return getBinaryModel(model_info, parameter_list,
                                 freq_type, tree, 
                                 report_to_tree);
        case SeqType::SEQ_CODON:
            return getCodonModel(model_info, parameter_list,
                                 freq_type, tree, 
                                 report_to_tree);

        case SeqType::SEQ_DNA:
            return getDNAModel(model_info, parameter_list,
                               freq_type, tree, 
                               report_to_tree);

        case SeqType::SEQ_MORPH:
            return getMorphologicalModel(model_info, parameter_list,
                                         freq_type, tree, 
                                         report_to_tree);

        case SeqType::SEQ_PROTEIN:
            return getProteinModel(model_info, parameter_list,
                                   freq_type, tree, models_block,
                                   report_to_tree);

        default:
            outError("YAML model uses unsupported sequence type");
            return nullptr;
    };
}

ModelMarkov* ModelListFromYAMLFile::getBinaryModel(ModelInfoFromYAMLFile& model_info,
                                                   const std::string& parameter_list,
                                                   StateFreqType freq_type,
                                                   PhyloTree* tree,
                                                   PhyloTree* report_to_tree) {
    insistOnAlignmentSequenceType(tree->aln, SeqType::SEQ_BINARY);
    YAMLModelBinary* model;
    model = new YAMLModelBinary(model_info, model_info.parent_model==nullptr,
                                "", dummy_rate_params, freq_type,
                                dummy_freq_params, tree, report_to_tree);
    model->acceptParameterList(*tree->params, parameter_list, report_to_tree);
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getCodonModel(ModelInfoFromYAMLFile& model_info,
                                                  const std::string& parameter_list,
                                                  StateFreqType freq_type,
                                                  PhyloTree* tree,
                                                  PhyloTree* report_to_tree) {
    insistOnAlignmentSequenceType(tree->aln, SeqType::SEQ_CODON);
    YAMLModelCodon* model;

    model = new YAMLModelCodon(model_info, model_info.parent_model==nullptr, 
                               "", dummy_rate_params, freq_type,
                               dummy_freq_params, tree, report_to_tree);
    model->acceptParameterList(*tree->params, parameter_list, report_to_tree);
    return model;
}
    
ModelMarkov* ModelListFromYAMLFile::getDNAModel(ModelInfoFromYAMLFile& model_info,
                                                const std::string& parameter_list,
                                                StateFreqType freq_type,
                                                PhyloTree* tree,
                                                PhyloTree* report_to_tree) {
    ModelMarkov* model = nullptr;
    const YAMLFileParameter* p =
        model_info.findParameter("epsilon",
                                 ModelParameterType::RATE);
    if (p!=nullptr) {
        if (p->is_subscripted) {
            outError("epsilon parameter for DNA+error model"
                     "may not be subscripted");
        }
        double epsilon = p->value;
        bool epsilon_is_fixed = model_info.isVariableFixed("epsilon");
        model_info.moveParameterToBack("epsilon",
                                       ModelParameterType::RATE);
        YAMLModelDNAError* emodel;
        emodel = new YAMLModelDNAError(model_info, model_info.parent_model==nullptr,
                                       "", dummy_rate_params, freq_type,
                                       dummy_freq_params, tree, report_to_tree);
        std::string error_model = model_info.getStringProperty("errormodel", "+E");
        emodel->setEpsilon(epsilon, epsilon_is_fixed, error_model);
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      "epsilon is " << epsilon
                      << ", fixed is " << epsilon_is_fixed
                      << ", and errormodel is " << error_model);
        emodel->acceptParameterList(*tree->params, parameter_list, report_to_tree);
        model = emodel;
    }
    else {
        YAMLModelDNA* dmodel;
        dmodel = new YAMLModelDNA(model_info, model_info.parent_model==nullptr,
                                  "", dummy_rate_params, freq_type,
                                  dummy_freq_params, tree,
                                  report_to_tree);
        dmodel->acceptParameterList(*tree->params, parameter_list, report_to_tree);
        model = dmodel;
    }
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getMorphologicalModel(ModelInfoFromYAMLFile& model_info,
                                                          const std::string& parameter_list,
                                                          StateFreqType freq_type,
                                                          PhyloTree* tree,
                                                          PhyloTree* report_to_tree) {
    YAMLModelMorphology* model;
    model = new YAMLModelMorphology(model_info, model_info.parent_model==nullptr,
                                    "", dummy_rate_params, freq_type,
                                    dummy_freq_params, tree, report_to_tree);
    model->acceptParameterList(*tree->params, parameter_list, report_to_tree);
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getProteinModel(ModelInfoFromYAMLFile& model_info,
                                                    const std::string& parameter_list,
                                                    StateFreqType freq_type,
                                                    PhyloTree* tree,
                                                    ModelsBlock* models_block,
                                                    PhyloTree* report_to_tree) {
    insistOnAlignmentSequenceType(tree->aln, SeqType::SEQ_PROTEIN);
    YAMLModelProtein* model;
    model = new YAMLModelProtein(model_info, model_info.parent_model==nullptr,
                                 "", dummy_rate_params, freq_type,
                                 dummy_freq_params, models_block,
                                 tree, report_to_tree);
    model->acceptParameterList(*tree->params, parameter_list, report_to_tree);
    return model;
}

void ModelListFromYAMLFile::insistOnAlignmentSequenceType(const Alignment* alignment, 
                                                          SeqType desired_type) {
    if (alignment->seq_type != desired_type) {
        std::stringstream complaint;
        complaint << "Cannot use " << getSeqTypeName(desired_type) << " model"
                  << " with a " << getSeqTypeName(alignment->seq_type) << " alignment."
                  << "\nPlease re-run, with"
                  << " a " << getSeqTypeName(alignment->seq_type) << " model,"
                  << " or passing -seqtype " << getSeqTypeShortName(desired_type, false)
                  << " on the command-line.";
        outError(complaint.str());
    }
}

ModelMarkov* ModelListFromYAMLFile::getMixtureModel(ModelInfoFromYAMLFile& model_info,
                                                    const std::string& parameter_list,
                                                    StateFreqType freq_type,
                                                    ModelsBlock* models_block,
                                                    PhyloTree* tree,
                                                    PhyloTree* report_to_tree) {
    YAMLModelMixture* model;
    model = new YAMLModelMixture(model_info, true, model_info.getName().c_str(),
                                 freq_type,  models_block, tree, report_to_tree);
    model->acceptParameterList(*tree->params, parameter_list, report_to_tree);
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getDivergentModel
    ( ModelInfoFromYAMLFile& model_info, 
      const std::string& parameter_list,
      StateFreqType freq_type, ModelsBlock* models_block,
      PhyloTree* tree, PhyloTree* report_to_tree) {
    YAMLModelDivergent* model;
    model = new YAMLModelDivergent(model_info, true, model_info.getName().c_str(),
                                   freq_type,  models_block, tree, report_to_tree);
    model->acceptParameterList(*tree->params, parameter_list, report_to_tree);
    return model;
}


