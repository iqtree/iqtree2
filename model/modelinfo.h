//
//modelinfo.h
//Created by James Barbetti on 28-Jan-2021
//

#ifndef modelinfo_h
#define modelinfo_h

#include <string>
#include <utils/tools.h> //for ASCType

extern VerboseMode YAMLModelVerbosity;

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

class ModelParameterRange: public std::pair<double,double> {
public:
    typedef std::pair<double,double> super;
    bool is_set;
    ModelParameterRange(): super(0,0), is_set(false) {}
};

enum ModelParameterType {
    RATE, FREQUENCY, WEIGHT, OTHER
};

class YAMLFileParameter {
public:
    //Member variables
    std::string         name;
    std::string         description;
    bool                is_subscripted;
    int                 minimum_subscript;
    int                 maximum_subscript;
    std::string         type_name;
    ModelParameterType  type;
    ModelParameterRange range;
    double              value;

    //Constructors
    YAMLFileParameter();
    YAMLFileParameter( const YAMLFileParameter& ) = default;
    YAMLFileParameter& operator= ( const YAMLFileParameter& ) = default;

    //Member functions
    std::string getSubscriptedVariableName(int subscript) const;
    bool isMatchFor(const std::string& match_name /* assumed: lower-case */,
                    ModelParameterType match_type) const;
};

class ModelVariable {
protected:
    ModelParameterRange range;
    ModelParameterType  type;
    double              value;
    bool                is_fixed;
public:
    ModelVariable();
    ModelVariable(ModelParameterType t, const ModelParameterRange& r, double v);
    ~ModelVariable() = default;
    ModelVariable(const ModelVariable& rhs) = default;
    ModelVariable& operator=(const ModelVariable& rhs) = default;
    void   setValue(double v);
    void   markAsFixed();
    double getValue() const;
    bool   isFixed() const;
};

class StringMatrix: public std::vector<StrVector> {
public:
    void makeRectangular(size_t column_count) {
        size_t row_count = size();
        for (int row_num=0; row_num<row_count; ++row_num) {
            StrVector& row = at(row_num);
            if (row.size() != column_count) {
                row.resize(row_count, "");
            }
        }
    }
};

class ModelInfoFromYAMLFile: public ModelInfo {
public:
    typedef std::vector<YAMLFileParameter>               Parameters;
    typedef std::map<std::string, ModelVariable>         Variables;
    typedef std::map<std::string, ModelInfoFromYAMLFile> MapOfModels;
    typedef  std::map<std::string, std::string>          StringMap;
private:
    std::string   model_name;       //
    std::string   model_file_path;  //
    std::string   citation;         //Citation string
    std::string   DOI;              //DOI for the publication (optional)
    std::string   data_type_name;   //
    int           num_states;       //number of states
    bool          reversible;       //don't trust this; check against rate matrix
    int           rate_matrix_rank; //
    StringMatrix  rate_matrix_expressions;    //row major (expression strings)
    std::string   rate_matrix_formula;
    int           tip_likelihood_rank;
    StringMatrix  tip_likelihood_expressions; //likewise (for tip likelihoods)
    std::string   tip_likelihood_formula;
    Parameters    parameters;      //parameters
    StateFreqType frequency_type;
    Variables     variables;
    MapOfModels*  mixed_models;
    StringMap     string_properties;
    mutable StrVector variable_names;
    
    friend class ModelListFromYAMLFile;
    friend class ModelFileLoader;
public:
    ModelInfoFromYAMLFile(); //Only ModelListFromYAMLFile uses it.
    ModelInfoFromYAMLFile(const ModelInfoFromYAMLFile& rhs);
    explicit ModelInfoFromYAMLFile(const std::string& file_path);
    ~ModelInfoFromYAMLFile();
    
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
    virtual bool isMixtureModel()                  const;
    virtual bool isModelFinder()                   const;
    virtual bool isModelFinderOnly()               const;
    virtual bool isPolymorphismAware()             const { return false; /*stub*/ }
    virtual bool isReversible()                    const;
    virtual bool isWeissAndVonHaeselerTest()       const { return false; /*stub*/ }

    ASCType extractASCType
        (std::string& leftover_name) const {
        FUNCTION_NOT_IMPLEMENTED;
        return ASC_VARIANT;
    }
    std::string extractMixtureModelList
        (std::string& leftover_name) const {
        FUNCTION_NOT_IMPLEMENTED;
        return "";
    }
    std::string extractPolymorphicHeterozygosity
        (std::string& leftover_name) const {
        FUNCTION_NOT_IMPLEMENTED;
        return "";
    }
    void updateName  (const std::string&       name);
    void addParameter(const YAMLFileParameter& p);
    
public:
    const std::string& getName ()                                       const;
    std::string getLongName    ()                                       const;
    
    bool hasDot    (const char* name)                                   const;
    void breakAtDot(const char* name, std::string& sub_model_name,
                    const char* &remainder)                             const;
    MapOfModels::const_iterator findMixedModel(const std::string& name) const;
    MapOfModels::iterator       findMixedModel(const std::string& name);

    //Parameters
    const YAMLFileParameter* findParameter(const char* name,
                                           ModelParameterType type) const;
    void  moveParameterToBack(const char* name,
                              ModelParameterType type);
    bool   isFrequencyParameter(const std::string& param_name)          const;
    std::string        getParameterList(ModelParameterType param_type)  const;
    void               appendParameterList(ModelParameterType param_type,
                                           std::stringstream& list) const;

    //Rate matrices
    int                getRateMatrixRank()                              const;
    const std::string& getRateMatrixExpression(int row, int col)        const;
    int                getNumStates() const;
    
    //Tip Likelihood matrices
    int  getTipLikelihoodMatrixRank() const;
    void computeTipLikelihoodsForState(int state, int num_states, double* likelihoods);
    
    //Variables
    bool   hasVariable         (const char* name)                       const;
    bool   hasVariable         (const std::string& name)                const;
    double getVariableValue    (const std::string& name)                const;
    double getVariableValue    (const char* name)                       const;
    bool   isVariableFixed     (const std::string& name)                const;
    void   setBounds           (int bound_count, double *lower_bound,
                                double *upper_bound, bool *bound_check) const;
    void   updateVariables     (const double* variables,
                                int param_count);
    void   logVariablesTo      (PhyloTree& report_to_tree)              const;
    ModelVariable& assign      (const std::string& var_name,
                                double value);
    ModelVariable& forceAssign (const std::string& var_name,
                                double value);
    const StrVector&            getVariableNamesByPosition() const;
    ModelVariable& assignByPosition(size_t position,
                                    double value_to_set);

    bool   assignLastFrequency (double value);

    //String Properties
    std::string getStringProperty(const char* name,
                                  const char* default_value)            const;
    
};

typedef std::map<std::string, ModelInfoFromYAMLFile> MapOfModels;

class ModelListFromYAMLFile {
protected:
    MapOfModels models_found;

public:
    friend class ModelFileLoader;
    ModelListFromYAMLFile()  = default;
    ~ModelListFromYAMLFile() = default;

    void loadFromFile          (const char* file_path,    PhyloTree* report_to_tree);
    bool isModelNameRecognized (const char* model_name);
    ModelMarkov* getModelByName(const char* model_name,   PhyloTree *tree,
                                const char* model_params, StateFreqType freq_type,
                                const char* freq_params,  PhyloTree* report_to_tree);
    
    bool hasModel(const std::string& model_name) const;
    const ModelInfoFromYAMLFile& getModel(const std::string& model_name) const;
};

#endif /* modelinfo_h */
