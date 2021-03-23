//
//  modelinfo.h
//  Created by James Barbetti on 28-Jan-2021
//

#ifndef modelinfo_h
#define modelinfo_h

#include <string>
#include <utils/tools.h> //for ASCType

class ModelInfo {
public:
    ModelInfo()                     = default;
    ModelInfo(const ModelInfo& rhs) = delete;
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

class ModelInfoFromYAMLFile: public ModelInfo {
private:
    std::string model_name;
    std::string model_file_path;
public:
    explicit ModelInfoFromYAMLFile(const std::string& file_path);
    ~ModelInfoFromYAMLFile() = default;
    
    virtual std::string getFreeRateParameters(int& num_rate_cats,
                                              bool &fused_mix_rate) const {
        FUNCTION_NOT_IMPLEMENTED;
        return "";
    }
    virtual std::string getFrequencyMixtureParams(std::string& freq_str) const {
        FUNCTION_NOT_IMPLEMENTED;
        return "";
    }
    virtual void        getFrequencyOptions(std::string& freq_str,
                                            StateFreqType& freq_type,
                                            std::string& freq_params,
                                            bool& optimize_mixmodel_weight) const {
        FUNCTION_NOT_IMPLEMENTED;
    }
    virtual void        getGammaParameters(int& num_rate_cats,
                                           double& gamma_shape) const {
        FUNCTION_NOT_IMPLEMENTED;
    }
    virtual std::string getHeterotachyParameters(bool is_mixture_model,
                                                 int& num_rate_cats,
                                                 bool& fused_mix_rate) const {
        FUNCTION_NOT_IMPLEMENTED;
        return "";
    }
    
    virtual double getProportionOfInvariantSites() const { return 0.0; /*stub*/ }

    virtual bool hasAscertainmentBiasCorrection()  const { return false; /*stub*/ }
    virtual bool hasRateHeterotachy()              const { return false; /*stub*/ }
    
    virtual bool isFreeRate()                      const { return false; /*stub*/ }
    virtual bool isFrequencyMixture()              const { return false; /*stub*/ }
    virtual bool isGammaModel()                    const { return false; /*stub*/ }
    virtual bool isInvariantModel()                const { return false; /*stub*/ }
    virtual bool isMixtureModel()                  const { return false; /*stub*/ }
    virtual bool isModelFinder()                   const { return false; /*stub*/ }
    virtual bool isModelFinderOnly()               const { return false; /*stub*/ }
    virtual bool isPolymorphismAware()             const { return false; /*stub*/ }
    virtual bool isWeissAndVonHaeselerTest()       const { return false; /*stub*/ }

    ASCType     extractASCType(std::string& leftover_name) const {
        FUNCTION_NOT_IMPLEMENTED;
        return ASC_VARIANT;
    }
    std::string extractMixtureModelList(std::string& leftover_name) const {
        FUNCTION_NOT_IMPLEMENTED;
        return "";
    }
    std::string extractPolymorphicHeterozygosity(std::string& leftover_name) const {
        FUNCTION_NOT_IMPLEMENTED;
        return "";
    }
    void        updateName(const std::string& name);
};


#endif /* modelinfo_hpp */
