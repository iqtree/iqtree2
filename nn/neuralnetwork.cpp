//
// Created by tamara on 3/8/21.
//

#include <assert.h>
#include "neuralnetwork.h"
#include <random>
#include <chrono>

NeuralNetwork::NeuralNetwork(Alignment *alignment) {

    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "alpha_find");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

    const char *model_path = "../nn_models/alpha_tmp_model.onnx";

    printf("Using Onnxruntime C++ API\n");
    Ort::Session session(env, model_path, session_options);
    Ort::AllocatorWithDefaultOptions allocator;

    size_t num_input_nodes = session.GetInputCount();
    std::vector<const char *> input_node_names(num_input_nodes);
    std::vector<int64_t> input_node_dims;

    printf("Number of inputs = %zu\n", num_input_nodes);

    // iterate over all input nodes
    for (int i = 0; i < num_input_nodes; i++) {
        // print input node names
        char *input_name = session.GetInputName(i, allocator);
        printf("Input %d : name=%s\n", i, input_name);
        input_node_names[i] = input_name;

        // print input node types
        Ort::TypeInfo type_info = session.GetInputTypeInfo(i);
        auto tensor_info = type_info.GetTensorTypeAndShapeInfo();

        ONNXTensorElementDataType type = tensor_info.GetElementType();
        printf("Input %d : type=%d\n", i, type);

        // print input dims
        input_node_dims = tensor_info.GetShape();
        printf("Input %d : num_dims=%zu\n", i, input_node_dims.size());
        for (int j = 0; j < input_node_dims.size(); j++) {
            if (input_node_dims[j] == -1) // TODO: check with Sebastian if there is a nicer way to deal with wildcard dim
                input_node_dims[j] = 1;
            printf("Input %d : dim %d=%jd\n", i, j, input_node_dims[j]);
        }
    }

    //std::vector<float> input_tensor_(10000 * 4);
    std::array<float, 10000 * 4> input_tensor_{};
    std::array<int64_t, 3> input_shape_{1, 10000, 4};

    //std::vector<const char*> output_node_names = {"ev_model"};
    std::vector<const char *> output_node_names = {"alpha"};

    const size_t num_taxa = alignment->getNSeq();
    const size_t num_sites = alignment->getNSite();

    // choose 10,000 random positions (with repetition) in (0, num_sites - 1)
    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<size_t> dist(0, num_sites - 1);

    // row major
    for (size_t i = 0; i < 10000; i += 4) {
        size_t site_idx = dist(rng);
        vector<size_t> freqs = alignment->getPattern(site_idx).freqs;
        input_tensor_[i] = (float) freqs[0] / num_taxa;
        input_tensor_[i + 1] = (float) freqs[1] / num_taxa;
        input_tensor_[i + 2] = (float) freqs[2] / num_taxa;
        input_tensor_[i + 3] = (float) freqs[3] / num_taxa;
    }

    // column major
//    for (size_t i = 0; i < 10000; i++) {
//        size_t site_idx = dist(rng);
//        auto mimi = alignment->getPattern(site_idx);
//        vector<size_t> freqs = alignment->getPattern(site_idx).freqs;
//        input_tensor_[i] = (float)freqs[0] / num_taxa;
//        input_tensor_[i+10000] = (float)freqs[1] / num_taxa;
//        input_tensor_[i+20000] = (float)freqs[2] / num_taxa;
//        input_tensor_[i+30000] = (float)freqs[3] / num_taxa;
//    }

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

    // print alpha value
    for (int i = 0; i < 1; i++)
        printf("Alpha value [%d] =  %f\n", i, floatarr[i]);

}

NeuralNetwork::~NeuralNetwork() {

}

// TODO: move inference from constructor to this function
int NeuralNetwork::doInference(int input) {
    return 0;
}
