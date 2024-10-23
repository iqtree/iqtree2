//
// Created by tamara on 3/8/21.
//

#include <assert.h>
#include "neuralnetwork.h"
#include <random>
#include <chrono>
#include <vector>

NeuralNetwork::NeuralNetwork(Alignment *alignment) {
    this->alignment = alignment;
}

NeuralNetwork::~NeuralNetwork() {}

double NeuralNetwork::doAlphaInference() {
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "alpha_find");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

    const char *model_path = Params::getInstance().nn_path_rates.c_str();

    // printf("Using Onnxruntime C++ API\n");
    Ort::Session session(env, model_path, session_options);

    size_t num_input_nodes = session.GetInputCount();
    std::vector<const char *> input_node_names(num_input_nodes);

    Ort::AllocatorWithDefaultOptions allocator;
    //std::vector<int64_t> input_node_dims;

    //printf("Number of inputs = %zu\n", num_input_nodes);
#ifndef _OLD_NN
    std::vector<Ort::AllocatedStringPtr> inputNodeNameAllocatedStrings;
#endif
    
    // iterate over all input nodes
    for (int i = 0; i < num_input_nodes; i++) {
        // print input node names
#ifndef _OLD_NN
        auto input_name = session.GetInputNameAllocated(i, allocator);
        inputNodeNameAllocatedStrings.push_back(std::move(input_name));
        input_node_names[i] = inputNodeNameAllocatedStrings.back().get();
#else
        char *input_name = session.GetInputName(i, allocator);
        input_node_names[i] = input_name;
#endif

    /*
        // print input node types
        Ort::TypeInfo type_info = session.GetInputTypeInfo(i);
        auto tensor_info = type_info.GetTensorTypeAndShapeInfo();

        ONNXTensorElementDataType type = tensor_info.GetElementType();
        printf("Input %d : type=%d\n", i, type);

        // print input dims
        input_node_dims = tensor_info.GetShape();
        printf("Input %d : num_dims=%zu\n", i, input_node_dims.size());
        for (int j = 0; j < input_node_dims.size(); j++) {
            printf("Input %d : dim %d=%jd\n", i, j, input_node_dims[j]);
        }*/
    }

    std::vector<float> input_tensor_(10000 * 4);
    std::vector<int64_t> input_shape_{1, 10000, 4};

    std::vector<const char *> output_node_names = {"alpha", "ev_model"};

    size_t num_sites = this->alignment->getNSite();
    size_t num_taxa = this->alignment->getNSeq();

    // choose 10,000 random positions (with repetition) in (0, num_sites - 1)
    //mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    //std::uniform_int_distribution<size_t> dist(0, num_sites - 1);

    for (size_t i = 0; i < 40000; i = i + 4) {
        size_t site_idx = random_int(num_sites); //dist(rng);
        vector<size_t> freqs = this->alignment->getPattern(site_idx).freqs;
        // in case of gaps, adjust number of taxa
        // size_t num_taxa = accumulate(freqs.begin(), freqs.end(), 0);

        input_tensor_[i] = (float) freqs[0] / num_taxa;
        input_tensor_[i + 1] = (float) freqs[1] / num_taxa;
        input_tensor_[i + 2] = (float) freqs[2] / num_taxa;
        input_tensor_[i + 3] = (float) freqs[3] / num_taxa;
    }

    // create input tensor object from data values
    auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    Ort::Value input_tensor = Ort::Value::CreateTensor<float>(memory_info, input_tensor_.data(), input_tensor_.size(),
                                                              input_shape_.data(), input_shape_.size());
    assert(input_tensor.IsTensor());

    // do inference
    auto output_tensors = session.Run(Ort::RunOptions{nullptr}, input_node_names.data(), &input_tensor, 1,
                                      output_node_names.data(), 2);
    assert(output_tensors.size() == 2 && output_tensors.front().IsTensor());

    // get pointer to output tensor float values
    //float *floatarr = output_tensors.front().GetTensorMutableData<float>();
    float *alpha = output_tensors[0].GetTensorMutableData<float>();
    float *check = output_tensors[1].GetTensorMutableData<float>();

    printf("Check whether heterogeneous (0) or homogeneous (1) = %f\n", check[0]);

    // print alpha value
    // printf("Alpha value =  %f\n", alpha[0] / 1000);

    if (check[0] > 0.5)
        return -1;
    return alpha[0] / 1000;
}

string NeuralNetwork::doModelInference() {
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "model_find");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

    const char *model_path = Params::getInstance().nn_path_model.c_str();

    // printf("Using Onnxruntime C++ API\n");
    // cout << "model_path: " << model_path << endl;
    Ort::Session session(env, model_path, session_options);

    size_t num_input_nodes = session.GetInputCount();
    std::vector<const char *> input_node_names(num_input_nodes);

    Ort::AllocatorWithDefaultOptions allocator;
    //std::vector<int64_t> input_node_dims;

    // printf("Number of inputs = %zu\n", num_input_nodes);

    std::vector<Ort::AllocatedStringPtr> inputNodeNameAllocatedStrings;

    // iterate over all input nodes
    for (int i = 0; i < num_input_nodes; i++) {
#ifndef _OLD_NN
        auto input_name = session.GetInputNameAllocated(i, allocator);
        inputNodeNameAllocatedStrings.push_back(std::move(input_name));
        input_node_names[i] = inputNodeNameAllocatedStrings.back().get();
#else
        char *input_name = session.GetInputName(i, allocator);
        input_node_names[i] = input_name;
#endif
    /*
        // print input node types
        Ort::TypeInfo type_info = session.GetInputTypeInfo(i);
        auto tensor_info = type_info.GetTensorTypeAndShapeInfo();

        ONNXTensorElementDataType type = tensor_info.GetElementType();
        printf("Input %d : type=%d\n", i, type);

        // print input dims
        input_node_dims = tensor_info.GetShape();
        printf("Input %d : num_dims=%zu\n", i, input_node_dims.size());
        for (int j = 0; j < input_node_dims.size(); j++) {
            printf("Input %d : dim %d=%jd\n", i, j, input_node_dims[j]);
        }*/
    }

    std::vector<float> input_tensor_(40 * 250 * 26);
    std::vector<int64_t> input_shape_{1, 40, 250, 26};

    std::vector<const char *> output_node_names = {"dense"};

    const size_t num_taxa = this->alignment->getNSeq();

    // choose 10,000 random sequence pairs (with repetition) in (0, num_sites - 1)
    //mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    //std::uniform_int_distribution<size_t> dist_taxa(0, num_taxa - 1);

    vector<float> summary_stats(26);

    this->alignment->ungroupSitePattern();

    for (size_t i = 0; i < 260000; i = i + 26) {
        size_t seq1_idx;
        size_t seq2_idx;
        while (true) {
            seq1_idx = random_int(num_taxa); //dist_taxa(rng);
            seq2_idx = random_int(num_taxa); //dist_taxa(rng);
            if (seq1_idx != seq2_idx)
                break;
        }

        summary_stats = this->alignment->computeSummaryStats(seq1_idx, seq2_idx);
        std::copy(std::begin(summary_stats), std::end(summary_stats), std::begin(input_tensor_) + i); // replace part of the vector with this
    }

    // create input tensor object from data values
    auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    Ort::Value input_tensor = Ort::Value::CreateTensor<float>(memory_info, input_tensor_.data(), input_tensor_.size(),
                                                              input_shape_.data(), input_shape_.size());
    assert(input_tensor.IsTensor());

    // do inference
    auto output_tensors = session.Run(Ort::RunOptions{nullptr}, input_node_names.data(), &input_tensor, 1,
                                      output_node_names.data(), 1);
    assert(output_tensors.size() == 1 && output_tensors.front().IsTensor());

    // get pointer to output tensor float values
    float *floatarr = output_tensors.front().GetTensorMutableData<float>();

    // print values for JC,K2P,F81,HKY,Tn,GTR
    float max_val = 0.0;
    size_t chosen_model;

    for (size_t i = 0; i < 6; i++) {
        if (floatarr[i] > max_val) {
            max_val = floatarr[i];
            chosen_model = i;
        }
        //printf("Model value [%zu] =  %f\n", i, floatarr[i]);
    }

    // return chosen model
    switch(chosen_model) {
        case 0: return "JC";
        case 1: return "K2P";
        case 2: return "F81";
        case 3: return "HKY";
        case 4: return "TN";
        case 5: return "GTR";
        default: throw "Model not known";

    }

}
