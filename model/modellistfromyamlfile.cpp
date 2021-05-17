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

#include "modelinfofromyamlfile.h"
#include "modelsubst.h"            //for OPEN_BRACKET and CLOSE_BRACKET
#include "modelfileloader.h"
#include "modelexpression.h"       //for InterpretedExpression
#include "yamlmodelwrapper.h"      //for YAMLModelWrapper template class

#include <utils/my_assert.h>       //for ASSERT macro
#include <utils/stringfunctions.h> //for convert_int
#include <utils/tools.h>           //for outError


void ModelListFromYAMLFile::loadFromFile (const char* file_path,
                                          PhyloTree* report_to_tree) {
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
                    TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                                "Parsing YAML rate model " << rate_model_name);
                    ModelInfoFromYAMLFile &rate = rate_models_found[rate_model_name]
                                                = ModelInfoFromYAMLFile();
                    rate.is_rate_model = true;
                    loader.parseYAMLModel(node, rate_model_name, rate, *this,
                                          nullptr, report_to_tree);
                }
                continue;
            }
            std::string yaml_model_name = node["substitutionmodel"].Scalar();
            TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                          "Parsing YAML substitution model " << yaml_model_name);
            ModelInfoFromYAMLFile &model = models_found[yaml_model_name]
                                         = ModelInfoFromYAMLFile();
            loader.parseYAMLModel(node, yaml_model_name, model, *this,
                                  nullptr, report_to_tree);
        }
    }
    catch (YAML::ParserException e) {
        outError(e.what());
    }
    catch (YAML::Exception e) {
        outError(e.what());
    }
    catch (ModelExpression::ModelException x) {
        outError(x.getMessage());
    }
}

bool ModelListFromYAMLFile::isSubstitutionModelNameRecognized (const char* model_name) {
    size_t i = 0;
    while (model_name[i]!='\0' && model_name[i]!='{') {
        ++i;
    }
    std::string model_front = std::string(model_name, i);

    auto found      = models_found.find(model_front);
    bool recognized = found != models_found.end();
    return recognized;
}

bool ModelListFromYAMLFile::hasModel(const std::string& model_name) const {
    return models_found.find(model_name) != models_found.end()
        || rate_models_found.find(model_name) != rate_models_found.end();
}

StrVector ModelListFromYAMLFile::getSubstitutionModelNames() const {
    StrVector answer;
    for (auto it = models_found.begin(); it != models_found.end(); ++it) {
        answer.push_back(it->first);
    }
    return answer;
}

std::string ModelListFromYAMLFile::getListOfSubstitutionModelNames() const {
    std::stringstream answer;
    const char* separator = "";
    for (auto it = models_found.begin(); it != models_found.end(); ++it) {
        answer << separator << it->first;
        separator = ", ";
    }
    return answer.str();
}

std::string ModelListFromYAMLFile::getListOfRateModelNames() const {
    std::stringstream answer;
    const char* separator = "";
    for (auto it = rate_models_found.begin(); 
         it != rate_models_found.end(); ++it) {
        answer << separator << it->first;
        separator = ", ";
    }
    return answer.str();
}

const ModelInfoFromYAMLFile&
    ModelListFromYAMLFile::getModel(const std::string& model_name) const {
    auto it = models_found.find(model_name);
    if (it!=models_found.end()) {
        return it->second;
    }
    it = rate_models_found.find(model_name);
    ASSERT(it != rate_models_found.end());
    return it->second;
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

ModelMarkov* ModelListFromYAMLFile::getModelByName(const char* model_name,   PhyloTree *tree,
                                                   const char* model_params, StateFreqType freq_type,
                                                   const char* freq_params,  ModelsBlock* models_block,
                                                   PhyloTree* report_to_tree) {
    std::string name;
    std::string parameter_list;
    extractModelNameAndParameters(model_name, name, parameter_list);
    ModelInfoFromYAMLFile& model_info = models_found[model_name];
    if (0<strlen(model_params) || 0<strlen(freq_params)) {
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      "Model Params: " << model_params
                      << " Freq Params: " << freq_params);
        if (parameter_list.empty()) {
            parameter_list = model_params;
        }
    }
    if (freq_type == StateFreqType::FREQ_UNKNOWN) {
        freq_type = model_info.frequency_type;
    }
    switch (model_info.sequence_type) {
        case SeqType::SEQ_BINARY:
            return getBinaryModel(model_info, parameter_list,
                                 freq_type, tree, report_to_tree);
        case SeqType::SEQ_CODON:
            return getCodonModel(model_info, parameter_list,
                                 freq_type, tree, report_to_tree);

        case SeqType::SEQ_DNA:
            return getDNAModel(model_info, parameter_list,
                               freq_type, tree, report_to_tree);

        case SeqType::SEQ_MORPH:
            return getMorphologicalModel(model_info, parameter_list,
                                         freq_type, tree, report_to_tree);

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
    model = new YAMLModelBinary("", dummy_rate_params, freq_type,
                                dummy_freq_params, tree, report_to_tree,
                                model_info);
    model->acceptParameterList(parameter_list);
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getCodonModel(ModelInfoFromYAMLFile& model_info,
                                                  const std::string& parameter_list,
                                                  StateFreqType freq_type,
                                                  PhyloTree* tree,
                                                  PhyloTree* report_to_tree) {
    insistOnAlignmentSequenceType(tree->aln, SeqType::SEQ_CODON);
    YAMLModelCodon* model;

    model = new YAMLModelCodon("", dummy_rate_params, freq_type,
                               dummy_freq_params, tree, report_to_tree,
                               model_info);
    model->acceptParameterList(parameter_list);
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
        emodel = new YAMLModelDNAError("", dummy_rate_params, freq_type,
                                      dummy_freq_params, tree,
                                      report_to_tree, model_info);
        std::string error_model = model_info.getStringProperty("errormodel", "+E");
        emodel->setEpsilon(epsilon, epsilon_is_fixed, error_model);
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      "epsilon is " << epsilon
                      << ", fixed is " << epsilon_is_fixed
                      << ", and errormodel is " << error_model);
        emodel->acceptParameterList(parameter_list);
        model = emodel;
    }
    else {
        YAMLModelDNA* dmodel;
        dmodel = new YAMLModelDNA("", dummy_rate_params, freq_type,
                                  dummy_freq_params, tree,
                                  report_to_tree, model_info);
        dmodel->acceptParameterList(parameter_list);
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
    model = new YAMLModelMorphology("", dummy_rate_params, freq_type,
                                  dummy_freq_params, tree, report_to_tree,
                                  model_info);
    model->acceptParameterList(parameter_list);
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
    model = new YAMLModelProtein(models_block, "", dummy_rate_params, freq_type,
                                  dummy_freq_params, tree, report_to_tree,
                                  model_info);
    model->acceptParameterList(parameter_list);
    return model;
}

void ModelListFromYAMLFile::insistOnAlignmentSequenceType(const Alignment* alignment, 
                                                          SeqType desired_type) const {
    if (alignment->seq_type != desired_type) {
        std::stringstream complaint;
        complaint << "Cannot use " << getSeqTypeName(desired_type) << " model"
                  << " with a " << getSeqTypeName(alignment->seq_type) << "alignment."
                  << "\nPlease re-run, with"
                  << " a " << getSeqTypeName(alignment->seq_type) << " model,"
                  << " or passing -seqtype " << getSeqTypeShortName(desired_type, false)
                  << " on the command-line.";
        outError(complaint.str());
    }
}