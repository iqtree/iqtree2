//
//modelinfo.h
//Created by James Barbetti on 28-Jan-2021
//

#ifndef modelinfo_h
#define modelinfo_h

#include <string>
#include <utils/tools.h> //for ASCType
#include <nclextra/modelsblock.h>

extern VerboseMode YAMLModelVerbosity;
extern VerboseMode YAMLVariableVerbosity;
extern VerboseMode YAMLFrequencyVerbosity;
extern VerboseMode YAMLMatrixVerbosity;


class ModelMarkov;
class PhyloTree;
class ModelFileLoader;

class ModelInfo {
public:
    ModelInfo()                     = default;
    ModelInfo(const ModelInfo& rhs) = default;
    virtual ~ModelInfo()            = default;
    
    virtual std::string getFreeRateParameters(int& num_rate_cats,
                                              bool &fused_mix_rate) const = 0;
    virtual std::string getFrequencyMixtureParams(std::string& freq_str) const  = 0;
    virtual void getFrequencyOptions(std::string& freq_str,
                                     StateFreqType& freq_type,
                                     std::string& freq_params,
                                     bool& optimize_mixmodel_weight) const = 0;

    virtual void   getGammaParameters(int& num_rate_cats,
                                      double& gamma_shape) const = 0;
    virtual std::string getHeterotachyParameters(bool is_mixture_model,
                                                 int& num_rate_cats,
                                                 bool& fused_mix_rate) const = 0;
    virtual double getProportionOfInvariantSites() const = 0;
    
    virtual bool hasAscertainmentBiasCorrection()  const = 0;
    virtual bool hasRateHeterotachy()              const = 0;
    
    virtual bool isFreeRate()                      const = 0;
    virtual bool isFrequencyMixture()              const = 0;
    virtual bool isGammaModel()                    const = 0;
    virtual bool isInvariantModel()                const = 0;
    virtual bool isMixtureModel()                  const = 0;
    virtual bool isModelFinder()                   const = 0;
    virtual bool isModelFinderOnly()               const = 0;
    virtual bool isPolymorphismAware()             const = 0;
    virtual bool isWeissAndVonHaeselerTest()       const = 0;
    
    virtual ASCType     extractASCType(std::string& leftover_name) const = 0;
    virtual std::string extractMixtureModelList(std::string& leftover_name) const = 0;
    virtual std::string extractPolymorphicHeterozygosity(std::string& leftover_name) const = 0;
    virtual void        updateName(const std::string& name) = 0;
};

class ModelInfoFromName: public ModelInfo {
private:
    std::string model_name;
public:
    friend class ModelFileLoader;
    
    explicit ModelInfoFromName(std::string name);
    explicit ModelInfoFromName(const char* name);
    virtual ~ModelInfoFromName() = default;
    
    virtual std::string getFreeRateParameters(int& num_rate_cats,
                                              bool &fused_mix_rate) const;
    virtual std::string getFrequencyMixtureParams(std::string& freq_str) const;
    virtual void        getFrequencyOptions(std::string& freq_str,
                                            StateFreqType& freq_type,
                                            std::string& freq_params,
                                            bool& optimize_mixmodel_weight) const;
    virtual void        getGammaParameters(int& num_rate_cats,
                                           double& gamma_shape) const;
    virtual std::string getHeterotachyParameters(bool is_mixture_model,
                                                 int& num_rate_cats,
                                                 bool& fused_mix_rate) const;
    virtual double getProportionOfInvariantSites() const;

    virtual bool hasAscertainmentBiasCorrection()  const;
    virtual bool hasRateHeterotachy()              const;
    
    virtual bool isFreeRate()                      const;
    virtual bool isFrequencyMixture()              const;
    virtual bool isGammaModel()                    const;
    virtual bool isInvariantModel()                const;
    virtual bool isMixtureModel()                  const;
    virtual bool isModelFinder()                   const;
    virtual bool isModelFinderOnly()               const;
    virtual bool isPolymorphismAware()             const;
    virtual bool isWeissAndVonHaeselerTest()       const;

    ASCType     extractASCType(std::string& leftover_name) const;
    std::string extractMixtureModelList(std::string& leftover_name) const;
    std::string extractPolymorphicHeterozygosity(std::string& leftover_name) const;
    void        updateName(const std::string& name);
};


#endif /* modelinfo_h */
