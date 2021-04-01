//
//  modelexpression.cpp
//  alignment
//
//  Created by James Barbetti on 1/4/21.
//

#include "modelexpression.h"

ModelExpression::ModelExpression(ModelInfoFromYAMLFile& for_model) : model(for_model) {
}

double ModelExpression::evaluate() const {
    return 0;
}

ModelVariable::ModelVariable(ModelInfoFromYAMLFile& for_model, const char* name)
    : super(for_model), variable_name(name) {
    if (!for_model.hasVariable(name)) {
        std::stringstream complaint;
        complaint << "Could not evaluate variable " << name
                  << " for model " << for_model.getLongName();
    }
}

double ModelVariable::evaluate() const {
    return model.getVariableValue(variable_name);
}

ModelInfix::ModelInfix(ModelInfoFromYAMLFile& for_model,
          ModelExpression* left, ModelExpression* right )
    : super(for_model), lhs(left), rhs(right) {
    ASSERT(lhs!=nullptr && rhs!=nullptr);
}

ModelInfix::~ModelInfix() {
    delete lhs;
    lhs = nullptr;
    delete rhs;
    rhs = nullptr;
}

ModelMultiplication::ModelMultiplication(ModelInfoFromYAMLFile& for_model,
                   ModelExpression* left, ModelExpression* right )
    : super(for_model, left, right) { }
double ModelMultiplication::evaluate() const { return lhs->evaluate() * rhs->evaluate(); }

ModelSubtraction::ModelSubtraction(ModelInfoFromYAMLFile& for_model,
                                   ModelExpression* left, ModelExpression* right )
    : super(for_model, left, right) { }
double ModelSubtraction::evaluate() const {
    return lhs->evaluate() - rhs->evaluate();
}
