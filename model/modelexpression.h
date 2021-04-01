//
//  modelexpression.h
//  alignment
//
//  Created by James Barbetti on 1/4/21.
//

#ifndef modelexpression_h
#define modelexpression_h

#include "modelinfo.h"  //for ModelInfoFromYAMLFile

class ModelExpression {
protected:
    const ModelInfoFromYAMLFile& model;
public:
    ModelExpression(ModelInfoFromYAMLFile& for_model);
    virtual ~ModelExpression() = default;
    virtual double evaluate() const;
};

class ModelVariable: public ModelExpression {
    std::string variable_name;
public:
    typedef ModelExpression super;
    ModelVariable(ModelInfoFromYAMLFile& for_model, const char* name);
    virtual ~ModelVariable() = default;
    virtual double evaluate() const;
};

class ModelInfix: public ModelExpression {
protected:
    ModelExpression* lhs;
    ModelExpression* rhs;
public:
    typedef ModelExpression super;
    ModelInfix(ModelInfoFromYAMLFile& for_model,
               ModelExpression* left, ModelExpression* right );
    virtual ~ModelInfix();
};

class ModelMultiplication: public ModelInfix {
public:
    typedef ModelInfix super;
    ModelMultiplication(ModelInfoFromYAMLFile& for_model,
                        ModelExpression* left, ModelExpression* right );
    virtual ~ModelMultiplication() = default;
    virtual double evaluate() const;
};

class ModelSubtraction: public ModelInfix {
    typedef ModelInfix super;
    ModelSubtraction(ModelInfoFromYAMLFile& for_model,
                     ModelExpression* left, ModelExpression* right );
    virtual ~ModelSubtraction() = default;
    virtual double evaluate() const;
};

#endif /* modelexpression_hpp */
