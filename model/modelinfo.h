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
extern VerboseMode YAMLParsingVerbosity;
extern VerboseMode YAMLRateVerbosity;
extern VerboseMode YAMLVariableVerbosity;
extern VerboseMode YAMLFrequencyVerbosity;
extern VerboseMode YAMLMatrixVerbosity;
extern VerboseMode YAMLWarningVerbosity;


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
    
    virtual bool isDivergentModel()                const = 0;
    virtual bool isFreeRate()                      const = 0;
    virtual bool isFrequencyMixture()              const = 0;
    virtual bool isGammaModel()                    const = 0;
    virtual bool isInvariantModel()                const = 0;
    virtual bool isKategoryModel()                 const = 0;
    virtual bool isMixtureModel()                  const = 0;
    virtual bool isModelFinder()                   const = 0;
    virtual bool isModelFinderOnly()               const = 0;
    virtual bool isPolymorphismAware()             const = 0;
    virtual bool isWeissAndVonHaeselerTest()       const = 0;

    virtual int  getKategoryRateCount(int rate_count, int min_count) const = 0;
    
    virtual ASCType     extractASCType(std::string& leftover_name) const = 0;
    virtual std::string extractMixtureModelList(std::string& leftover_name) const = 0;
    virtual std::string extractPolymorphicHeterozygosity(std::string& leftover_name) const = 0;
    virtual void        updateName(const std::string& name) = 0;
};

class ModelInfoFromName: public ModelInfo {
private:
    std::string model_name;
    int getNumberOfCategories(std::string model_name,  string::size_type posH, 
                              bool is_mixture_model,   bool fused_mix_rate, 
                              int &end_pos) const;
    std::string getOtherHeterotachyParameters
                (std::string model_name, string::size_type posH, 
                 int end_pos) const;
    void extractUserDefinedFrequency(std::string fstr,  StateFreqType& freq_type,
                                     std::string& freq_params) const;

public:
    friend class ModelFileLoader;
    
    explicit ModelInfoFromName(const std::string& name);
    explicit ModelInfoFromName(const char* name);
    virtual ~ModelInfoFromName() = default;
    
    virtual std::string getFreeRateParameters(int& num_rate_cats,
                                              bool &fused_mix_rate) const override;
    virtual std::string getFrequencyMixtureParams(std::string& freq_str) const override;
    virtual void        getFrequencyOptions(std::string& freq_str,
                                            StateFreqType& freq_type,
                                            std::string& freq_params,
                                            bool& optimize_mixmodel_weight) const override;
    virtual void        getGammaParameters(int& num_rate_cats,
                                           double& gamma_shape) const override;
    virtual std::string getHeterotachyParameters(bool is_mixture_model,
                                                 int& num_rate_cats,
                                                 bool& fused_mix_rate) const override;
    virtual double getProportionOfInvariantSites() const override;

    virtual bool hasAscertainmentBiasCorrection()  const override;
    virtual bool hasRateHeterotachy()              const override;
    
    virtual bool isDivergentModel()                const override;
    virtual bool isFreeRate()                      const override;
    virtual bool isFrequencyMixture()              const override;
    virtual bool isGammaModel()                    const override;
    virtual bool isInvariantModel()                const override;
    virtual bool isKategoryModel()                 const override;
    virtual bool isMixtureModel()                  const override;
    virtual bool isModelFinder()                   const override;
    virtual bool isModelFinderOnly()               const override;
    virtual bool isPolymorphismAware()             const override;
    virtual bool isWeissAndVonHaeselerTest()       const override;

    virtual int  getKategoryRateCount(int rate_count, int min_count) const override;

    ASCType     extractASCType(std::string& leftover_name) const override;
    std::string extractMixtureModelList(std::string& leftover_name) const override;
    std::string extractPolymorphicHeterozygosity(std::string& leftover_name) const override;
    void        updateName(const std::string& name) override;
};


#endif /* modelinfo_h */
