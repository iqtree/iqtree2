#pragma once
#ifndef modelinfofromyamlfile_h
#define modelinfofromyamlfile_h

#include "modelinfo.h"
#include <utils/stringfunctions.h> //for string_to_lower

class Alignment;
class RateHeterogeneity;

class ModelParameterRange : public std::pair<double, double> {
public:
    typedef std::pair<double, double> super;
    bool is_set;
    ModelParameterRange() : super(0, 0), is_set(false) {}
    ModelParameterRange& operator= (const ModelParameterRange&) = default;
    ~ModelParameterRange() = default;
};

enum class ModelParameterType {
    FREQUENCY, PROPORTION, INVARIANT_PROPORTION, 
    RATE, SHAPE, WEIGHT, OTHER
};

inline ModelParameterType modelParameterTypeFromString(const std::string& type_name) {
    std::string s = string_to_lower(type_name);
    if (s=="frequency")            return ModelParameterType::FREQUENCY;
    if (s=="other")                return ModelParameterType::OTHER;
    if (s=="proportion")           return ModelParameterType::PROPORTION;
    if (s=="invariant_proportion") return ModelParameterType::INVARIANT_PROPORTION;
    if (s=="rate")                 return ModelParameterType::RATE;
    if (s=="shape")                return ModelParameterType::SHAPE;
    if (s=="weight")               return ModelParameterType::WEIGHT;
    return ModelParameterType::RATE;
}

inline std::string modelParameterTypeToString(ModelParameterType type_code) {
    if (type_code == ModelParameterType::FREQUENCY)            return "frequency";
    if (type_code == ModelParameterType::INVARIANT_PROPORTION) return "proportion";
    if (type_code == ModelParameterType::OTHER)                return "other";
    if (type_code == ModelParameterType::PROPORTION)           return "invariant_proportion";
    if (type_code == ModelParameterType::RATE)                 return "rate";
    if (type_code == ModelParameterType::SHAPE)                return "shape";
    if (type_code == ModelParameterType::WEIGHT)               return "weight";
    return "";
}


class YAMLFileParameter {
public:
    //Member variables
    std::string         name;
    std::string         description;
    bool                is_subscripted;
    std::string         subscript_expression; //Rate models often have parameters
                                              //subscripted by the number of categories
    int                 minimum_subscript;
    int                 maximum_subscript;
    std::string         type_name;
    ModelParameterType  type;
    ModelParameterRange range;
    double              tolerance;            //Tolerance
    std::string         init_expression;      //Expression for initializing parameter
                                              //(Rate models may need to evaluate)
    double              value;

    //Constructors
    YAMLFileParameter();
    ~YAMLFileParameter() = default;
    YAMLFileParameter(const YAMLFileParameter&) = default;
    YAMLFileParameter& operator= (const YAMLFileParameter&) = default;

    //Member functions
    std::string getSubscriptedVariableName(int subscript) const;
    bool isMatchFor(const std::string& match_name) const; /* assumed: lower-case */
    bool isMatchFor(const std::string& match_name,        /* assumed: lower-case */
                    ModelParameterType match_type) const;

    void logParameterState(const char* verb, PhyloTree* report_to_tree) const;
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

    ModelVariable&     operator=(const ModelVariable& rhs) = default;

    void               setTypeName(const std::string& type_name);
    void               setMinimum (double minimum_value);
    void               setMaximum (double maximum_value);
    void               setValue   (double v);
    void               markAsFixed();
    bool               constrainValueToRange();

    std::string        getTypeName() const;
    ModelParameterType getType    () const;
    double             getValue   () const;
    bool               isFixed    () const;

};

class StringMatrix : public std::vector<StrVector> {
public:
    void makeRectangular(size_t column_count);

    /**
        makes a matrix (typically a rate matrix square), possibly  by
        reflecting expressions across the diagonal.
        @param reflect true if expressions to be reflected
     */
    void makeSquare(bool reflect);
};

class Checkpoint;

namespace ModelExpression { class Expression; }

class ModelInfoFromYAMLFile;

template <class R> class NamedPointerMap: protected std::vector< R* > {
protected:
    std::map<std::string, int> name_map;
    std::vector<std::string>   names;
public:
    typedef std::vector< R* > super;
    using super::at;
    using super::begin;
    using super::clear;
    using super::push_back;
    using super::size;

    //Constructors, Assignment operators, Destructors
    //===============================================
    NamedPointerMap() = default;
    NamedPointerMap(const NamedPointerMap& rhs) {
        copyFrom(rhs);
    }
    NamedPointerMap& operator=(const NamedPointerMap& rhs) {
        if (this==&rhs) {
            return *this;
        }
        super::clear();
        name_map = rhs.name_map;
        names    = rhs.names;
        for (const R* pointer: rhs) {
            push_back(pointer);
        }
        return *this;
    }
    NamedPointerMap& copyFrom(const NamedPointerMap& rhs) {
        for (size_t i=0; i< rhs.size(); ++i) {
            std::string name = rhs.names[i];
            ASSERT(!hasName(name));
            insert(name, new R(*(rhs[i])));
        }
        return *this;
    }
    virtual ~NamedPointerMap() {
        name_map.clear();   
        for (R* pointer: *this) {
            delete pointer;
        }
        super::clear();
    }

    //Adding
    //======
    R* insertNew(const std::string& name) {
        R* newPointer = new R();
        insert(name, newPointer);
        return newPointer;
    }
    void insert(const std::string& name, R* pointer) {
        int index = static_cast<int>(size());
        push_back(pointer);
        names.push_back(name);
        name_map[name] = index;
    }

    //Querying
    //========
    const std::vector<std::string>& getNames() const {
        return names;
    }
    typename super::iterator find(const std::string& name) {
        auto it = name_map.find(name);
        ASSERT(it!=name_map.end());
        return begin() + it->second;

    }
    typename super::const_iterator find(const std::string& name) const {
        auto it = name_map.find(name);
        ASSERT(it!=name_map.end());
        return begin() + it->second;
    }
    bool hasName(const std::string& name) const {
        return name_map.find(name) != name_map.end();
    }
    R* getPointer(const std::string& name) const {
        R* pointer_found = nullptr;
        auto it = name_map.find(name);
        if (it!=name_map.end()) {
            pointer_found = at(it->second);
        }
        return pointer_found;
    }

    //Iterators (I would prefer to disable the non-cost ones,
    //but Clang doesn't seem to want to play ball).)
    typedef typename super::const_iterator const_iterator;
    typedef typename super::iterator       iterator;

    iterator       begin()       {  return super::begin(); }
    iterator       end()         {  return super::end();   }

    const_iterator begin() const { return super::begin(); }
    const_iterator end()   const { return super::end(); }
};

class MapOfModels: public NamedPointerMap<ModelInfoFromYAMLFile> {
public:
    typedef NamedPointerMap<ModelInfoFromYAMLFile> super;
    typedef typename super::const_iterator const_iterator;
    MapOfModels()  = default;
    MapOfModels(const MapOfModels& rhs) = default;
    ~MapOfModels() = default;
};

#define PROPERTY_NAME_ERROR_MODEL "errormodel"
#define PROPERTY_NAME_ASC         "ascertainmentbiascorrection"

class ModelInfoFromYAMLFile : public ModelInfo {

    friend class ModelListFromYAMLFile;
    friend class ModelFileLoader;

public:
    typedef std::vector<YAMLFileParameter>                Parameters;
    typedef std::map<std::string, ModelVariable>          Variables;
    typedef std::map<std::string, std::string>            StringMap;

private:
    std::string   model_name;         //
    std::string   model_file_path;    //
    std::string   superclass_model_name; //The name of the superclass model (if any)
    bool          is_modifier_model;  //True for models that work by modifying
                                      //or extending an existing model.  False for
                                      //others.
    bool          is_rate_model;      //True if this is a *rate*, rather than a
                                      //substitution model.
    std::string   rate_distribution;  //Indicates what sort of rate distribution it is
                                      //(Even though, I'm not sure this is really necessary!)

    std::string   citation;         //Citation string
    std::string   DOI;              //DOI for the publication (optional)
    std::string   url;              //URL for the model (optional).
    std::string   data_type_name;   //
    SeqType       sequence_type;    //sequence type (perhaps from data_type_name ?)
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
    StringMap     string_properties;
    mutable StrVector variable_names;

    std::string   opt_algorithm;  //optimization algorithm

    //Only of relevance for child modules in mixtures
    ModelInfoFromYAMLFile* parent_model;
    MapOfModels*           mixed_models;   //nullptr, except for model mixtures
    MapOfModels*           linked_models;  //nullptr, except for tree mixtures
    std::string            weight_formula; //formula for deciding on weight
    double                 model_weight;   //weight

    //Set if a subsitution model has a specified rate model
    ModelInfoFromYAMLFile* specified_rate_model_info;

protected:
    void appendTo(const std::string& append_me,
                  const char* with_sep,
                  std::string& to_me);

    bool checkIntConsistent(const std::string& value_source,
                            const char* int_name,
                            int new_value, int &old_value,
                            std::stringstream& complaint);
    void copyMixedAndLinkedModels (const ModelInfoFromYAMLFile& rhs);
    void changeParameterSubscriptRange(int new_min, int new_max, 
                                       YAMLFileParameter& param);
    void setSubscriptedVariable   (const YAMLFileParameter& p, int i);
    bool removeSubscriptedVariable(const YAMLFileParameter&p, int i);

public:
    ModelInfoFromYAMLFile(); //Only ModelListFromYAMLFile uses it.
    ModelInfoFromYAMLFile(const ModelInfoFromYAMLFile& rhs);
    explicit ModelInfoFromYAMLFile(const std::string& file_path);
    ~ModelInfoFromYAMLFile();
    ModelInfoFromYAMLFile& operator=(const ModelInfoFromYAMLFile& rhs);
    void specifyRateModel(const ModelInfoFromYAMLFile& ancestor);

    virtual std::string getFreeRateParameters(int& num_rate_cats,
        bool& fused_mix_rate) const {
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

    virtual bool hasRateHeterotachy()              const;

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

    virtual ASCType extractASCType                  (std::string& leftover_name)        const;
    bool            checkAscertainmentBiasCorrection(bool warnIfInvalid, ASCType &type) const;
    virtual bool    hasAscertainmentBiasCorrection  ()                                  const;

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
    void updateName  (const std::string& name);
    void addParameter(const YAMLFileParameter& p);

/***
 * @param p
 * @param min_subscript
 * @param max_subsciprt
 * @note  Assumes: at least one of the subscripts has changed
 *        Doesnt assume: The new subscripts make sense.
 ***/
    void updateParameterSubscriptRange(YAMLFileParameter& p, 
                                       int min_subscript,
                                       int max_subscript,
                                       PhyloTree* report_to_tree);

public:
    const std::string& getName()                                        const;
    std::string        getLongName()                                    const;

    bool hasDot    (const char* name)                                   const;
    void breakAtDot(const char* name, std::string& sub_model_name,
                    const char*& remainder)                             const;
    MapOfModels::const_iterator findMixedModel(const std::string& name) const;

    //Initialization helper functions
    void setNumberOfStatesAndSequenceType(int requested_num_states,
                                          PhyloTree* report_to_tree);
    double evaluateExpression(std::string& expression, std::string context);

    //Parameters and variables
    const YAMLFileParameter* findParameter(const char* name)            const;
    const YAMLFileParameter* findParameter(const std::string& name)     const;
    const YAMLFileParameter* findParameter(const char* name,
                                           ModelParameterType type)     const;
    const YAMLFileParameter* findParameter(const std::string& name,
                                           ModelParameterType type)     const;
    bool  hasFrequencyParameters(int min_variable_count)                const;
    void  moveParameterToBack   (const char* name,
                                 ModelParameterType type);
    bool  isFrequencyParameter  (const std::string& param_name)         const;
    StateFreqType      getFrequencyType()                               const;
    void               setFrequencyType(StateFreqType new_type);
    std::string        getParameterList(ModelParameterType param_type)  const;
    bool               appendParameterList(ModelParameterType param_type,
                                           std::stringstream& list)     const;
    double             getModelWeight()                                 const;
    double             getModelWeight(); //This version re-evaluates weight_formula
    bool               isModelWeightFixed(); //Can't be const.

    //Rate matrices and associated rate models
    int                getRateMatrixRank()                              const;
    const std::string& getRateMatrixFormula()                           const;
    const std::string& getRateMatrixExpression(int row, int col)        const;
    int                getNumStates()                                   const;
	RateHeterogeneity* getSpecifiedRateModel(PhyloTree* tree)           const;
    RateHeterogeneity* getRateHeterogeneity(PhyloTree* tree)            const;

    //Tip Likelihood matrices
    int  getTipLikelihoodMatrixRank() const;
    void computeTipLikelihoodsForState(int state, int num_states, double* likelihoods);

    //Variables
    const ModelVariable* getVariableByName(const char*        name) const;
    const ModelVariable* getVariableByName(const std::string& name) const;
    bool   hasVariable     (const char* name)        const;
    bool   hasVariable     (const std::string& name) const;
    double getVariableValue(const std::string& name) const;
    double getVariableValue(const char* name)        const;
    bool   isVariableFixed (const std::string& name) const;
    void   setBounds       (int bound_count, 
                            const std::vector<ModelParameterType>& types,
                            double* lower_bound, double* upper_bound, 
                            bool* bound_check) const;
    void   updateVariables (const double* variables,
                            int first_freq_index, int last_param_index);
    bool   updateModelVariablesByType(const double* updated_values,
                                      int param_count, bool even_fixed_ones,
                                      ModelParameterType param_type,
                                      int &i);
    void   readModelVariablesByType  (double* write_them_here,
                                      int param_count, bool even_fixed_ones,
                                      ModelParameterType param_type,
                                      int &i) const;                                          


    void   logVariablesTo(PhyloTree& report_to_tree)           const;
    ModelVariable& assign(const std::string& var_name,
                          double value);
    ModelVariable& assign(const std::string& var_name,  
                          ModelExpression::Expression *x,
                          bool fix, const char* how,
                          PhyloTree* report_tree);


    ModelVariable& forceAssign(const std::string& var_name,
                               double value);
    ModelVariable& forceAssign(const std::string& var_name,
                               int value);
    ModelVariable& forceAssign(const char* var_name,
                               int value);

    const StrVector& getVariableNamesByPosition() const;
    /**
        note:   Assumes getVariableNamesByPosition() called in advance
    */
    std::string getVariableNameByPosition(int position) const;

    ModelVariable& assignByPosition(size_t position,
        double value_to_set);

    bool   assignLastFrequency(double value);
    const  ModelVariable* getInvariantProportionVariable() const;

    //String Properties
    std::string getStringProperty(const char* name,
        const char* default_value)            const;
    bool hasStringProperty(const char* name, std::string& value) const;

    //Inheriting
    void inheritModel(const ModelInfoFromYAMLFile &);
    void inheritModelProperties(const ModelInfoFromYAMLFile& mummy,
                                std::stringstream& complaint);
    void inheritModelParameters(const ModelInfoFromYAMLFile& mummy,
                                std::stringstream& complaint);
    void inheritModelVariables (const ModelInfoFromYAMLFile& mummy,
                                std::stringstream& complaint);
    void inheritModelMatrices  (const ModelInfoFromYAMLFile& mummy,
                                std::stringstream& complaint);

    //Accepting parameters supplied on command-line
    bool acceptParameterList   (std::string parameter_list,
                                PhyloTree* report_tree);

    //Mixture Model stuff
    MapOfModels& getMixedModels();
    const MapOfModels& getMixedModels() const;
    void copyVariablesFrom(const ModelInfoFromYAMLFile* original);
    void setParentModel(ModelInfoFromYAMLFile* parent);

    //Rate-Only stuff
    const std::string& getOptimizationAlgorithm() const;
    int getNumberOfRateCategories()      const;
    int getNumberOfVariableRates()       const;
    int getNumberOfVariableShapes()      const;
    int getNumberOfProportions()         const;
    int getNumberOfVariableProportions() const; 

    //Output
    void writeInfo(const char* caption, ModelParameterType type,
                   std::ostream& out) const;

    //Checkpoints
    void saveToCheckpoint     (Checkpoint* checkpoint) const;
    void restoreFromCheckpoint(Checkpoint* checkpoint);
};

class ModelListFromYAMLFile {
protected:
    MapOfModels  models_found;
    MapOfModels  rate_models_found;

public:
    friend class ModelFileLoader;
    ModelListFromYAMLFile() = default;
    ~ModelListFromYAMLFile() = default;

    void loadFromFile(const char* file_path, PhyloTree* report_to_tree);
    bool isSubstitutionModelNameRecognized(const char* model_name);

    ModelMarkov* getModelByName(const char* model_name,   PhyloTree* tree,
                                const char* model_params, StateFreqType freq_type,
                                const char* freq_params,  ModelsBlock* blocks_model,
                                PhyloTree* report_to_tree);

    static ModelMarkov* getModelByReference
                    (ModelInfoFromYAMLFile& model_info, PhyloTree *tree,
                     StateFreqType freq_type,           ModelsBlock* models_block,
                     const std::string &parameter_list, 
                     PhyloTree* report_to_tree);

    static ModelMarkov* getBinaryModel(ModelInfoFromYAMLFile& model_info,
                                       const std::string& parameter_list,
                                       StateFreqType freq_type, PhyloTree* tree, 
                                       PhyloTree* report_to_tree);

    static ModelMarkov* getCodonModel(ModelInfoFromYAMLFile& model_info,
                                      const std::string& parameter_list,
                                      StateFreqType freq_type, PhyloTree* tree,
                                      PhyloTree* report_to_tree);

    static ModelMarkov* getDNAModel(ModelInfoFromYAMLFile& model_info,
                                    const std::string& parameter_list,
                                    StateFreqType freq_type, PhyloTree* tree,
                                    PhyloTree* report_to_tree);

    static ModelMarkov* getMorphologicalModel(ModelInfoFromYAMLFile& model_info,
                                              const std::string& parameter_list,
                                              StateFreqType freq_type, PhyloTree* tree,
                                              PhyloTree* report_to_tree);

    static ModelMarkov* getProteinModel(ModelInfoFromYAMLFile& model_info,
                                        const std::string& parameter_list,
                                        StateFreqType freq_type, 
                                        PhyloTree*    tree,
                                        ModelsBlock*  models_block, 
                                        PhyloTree*    report_to_trees);

    static ModelMarkov* getMixtureModel(ModelInfoFromYAMLFile& model_info,
                                        const std::string& parameter_list,
                                        StateFreqType freq_type, 
                                        ModelsBlock*  models_block, 
                                        PhyloTree*    tree, 
                                        PhyloTree*    report_to_tree);

    static void insistOnAlignmentSequenceType(const Alignment* alignment, 
                                              SeqType desired_type);

    bool hasModel(const std::string& model_name)  const;
    StrVector   getSubstitutionModelNames()       const;
    std::string getListOfSubstitutionModelNames() const;
    std::string getListOfRateModelNames()         const;
    
    const ModelInfoFromYAMLFile* getModel(const std::string& model_name) const;
};

#endif //modelinfofromyamlfile_h
