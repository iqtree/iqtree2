//
// modelexpression.cpp
// Implementations of classes used for analyzing expressions
// (e.g. for values of rate matrix entries) in YAML model files.
// Created by James Barbetti on 01-Apr-2021
//

#include "modelexpression.h"
#include <utils/stringfunctions.h> //for convert_double
#include <model/rategamma.h>       //for RateGamma::cmpPointChi2

namespace ModelExpression {

    ModelException::ModelException(const char* s) : message(s) {
    }
    ModelException::ModelException(const std::string& s) : message(s) {        
    }
    ModelException::ModelException(const std::stringstream& s) : message(s.str()) {
    }
    const std::string& ModelException::getMessage() const {
        return message;
    }
    /*static*/ void ModelException::throwIfNonBlank(const std::stringstream& complaint) {
        std::string problem = complaint.str();
        if (!problem.empty()) {
            throw ModelException(problem);
        }
    }

    class BuiltIns {
    public:
        class Exp : public UnaryFunctionImplementation {
            double callFunction(ModelInfoFromYAMLFile& mf,
                                double parameter) const {
                return exp(parameter);
            }
        } exp_body;
        class Logarithm: public UnaryFunctionImplementation {
            double callFunction(ModelInfoFromYAMLFile& mf,
                                double parameter) const {
                return log(parameter);
            }
        } ln_body;
        class ChiSquared: public MultiFunctionImplementation {
        public:
            typedef MultiFunctionImplementation super;
            ChiSquared(): super(2,2) {}
            double callFunction(ModelInfoFromYAMLFile& mf,
                                ModelExpression::ListOperator* parameter_list) const {
                double p1 = parameter_list->evaluateEntry(0);
                double p2 = parameter_list->evaluateEntry(1);
                return RateGamma::cmpPointChi2(p1, p2);
            } 
        } chi2_body;

        std::map<std::string, UnaryFunctionImplementation*> unary_functions;
        std::map<std::string, MultiFunctionImplementation*> multi_functions;
        
        BuiltIns() {
            unary_functions["exp"]  = &exp_body;
            unary_functions["ln"]   = &ln_body;
            multi_functions["chi2"] = &chi2_body;
        }

        bool isUnaryFunction(const std::string &name) {
            return unary_functions.find(name) != unary_functions.end();
        }

        bool isMultiFunction(const std::string &name) {
            return multi_functions.find(name) != multi_functions.end();
        }

        bool isBuiltIn(const std::string &name) {
            return isUnaryFunction(name) || isMultiFunction(name);
        }

        Expression* asBuiltIn(ModelInfoFromYAMLFile& mf,
                              const std::string& name) {
            if (isUnaryFunction(name)) {
                return new UnaryFunction(mf, name.c_str(), unary_functions[name]);
            } else if (isMultiFunction(name) ) {
                return new MultiFunction(mf, name.c_str(), multi_functions[name]);
            }
            std::stringstream complaint;
            complaint << "Function " << name << " not recognized";
            throw ModelExpression::ModelException(complaint.str());
        }
    } built_in_functions;

    Expression::Expression(ModelInfoFromYAMLFile& for_model) : model(for_model) {
    }

    double Expression::evaluate()           const { return 0; }
    bool   Expression::isAssignment()       const { return false; }
    bool   Expression::isBoolean()          const { return false; }
    bool   Expression::isConstant()         const { return false; }
    bool   Expression::isFunction()         const { return false; }
    bool   Expression::isList()             const { return false; }
    bool   Expression::isOperator()         const { return false; }
    bool   Expression::isRange()            const { return false; }
    bool   Expression::isRightAssociative() const { return false; }
    bool   Expression::isToken(char c)      const { return false; }
    bool   Expression::isVariable()         const { return false; }
    int    Expression::getPrecedence() const { return 0; }
    int    Expression::evaluateAsInteger() const {
        //
        //Todo: better error handling.  Probably want to throw
        //      a ModelExpression::ModelException if not an integral
        //      integer (in a sensible range).
        //
        return (int)floor(evaluate());
    }

    std::string Expression::getText() const {
        std::stringstream s;
        writeTextTo(s);
        return s.str();
    }


    Token::Token(ModelInfoFromYAMLFile& for_model, char c): super(for_model) {
        token_char = c;
    }

    bool Token::isToken(char c) const {
        return (token_char == c);
    }

    void Token::writeTextTo(std::stringstream &text) const {
        text << token_char;
    }

    Constant::Constant(ModelInfoFromYAMLFile& for_model,
                       double v) : super(for_model), value(v) {
    }

    double Constant::evaluate() const {
        return value;
    }

    void   Constant::writeTextTo(std::stringstream &text) const {
        text << value;
    }

    bool Constant::isConstant() const {
        return true;
    }

    Variable::Variable(ModelInfoFromYAMLFile& for_model) 
        : super(for_model) {
    }

    Variable::Variable(ModelInfoFromYAMLFile& for_model,
                       const std::string& name)
        : super(for_model), variable_name(name) {
        if (!for_model.hasVariable(name)) {
            std::stringstream complaint;
            complaint << "Could not evaluate variable " << name
                      << " for model " << for_model.getLongName();
            throw ModelException(complaint.str());
        }
    }

    double Variable::evaluate() const {
        double v =  model.getVariableValue(variable_name);
        return v;
    }

    void Variable::writeTextTo(std::stringstream &text) const {
        text << variable_name;
    }

    bool Variable::isVariable() const {
        return true;
    }

    const std::string& Variable::getName() const {
        return variable_name;
    }

    ParameterSubscript::ParameterSubscript(ModelInfoFromYAMLFile& for_model,
                                          const YAMLFileParameter* param,
                                          Expression* expr) 
        : super(for_model)
        , parameter_to_subscript(param->name)
        , subscript_expression(expr) {
        if (!param->is_subscripted) {
            delete expr;
            std::stringstream complaint;
            complaint << "Cannot subscript parameter " << param->name
                      << " (it is not subscripted).";
            throw ModelException(complaint.str());
        }
        std::stringstream long_name;
        long_name <<  param->name + "(" ;
        expr->writeTextTo(long_name);
        long_name << ")";
        variable_name = long_name.str();
        //std::cout << "PS varname is " << variable_name << std::endl;
    }

    ParameterSubscript::~ParameterSubscript() {
        delete subscript_expression;
    }

    double ParameterSubscript::evaluate() const {
        return model.getVariableValue(getName());
    }

    std::string ParameterSubscript::getName() const {
        const YAMLFileParameter* param 
            = model.findParameter(parameter_to_subscript);
        if (param==nullptr) {
            std::stringstream complaint;
            complaint << "Parameter " << parameter_to_subscript
                      << ", not found.";
            throw ModelException(complaint.str());
        }
        double x = subscript_expression->evaluate();
        int    i = (int)floor(x);
        int   lo = param->minimum_subscript;
        int   hi = param->maximum_subscript;
        if (x<lo || hi<x) {
            std::stringstream complaint;
            complaint << "Invalid subscript " << i
                      << " for parameter " << param->name
                      << ", for which the valid subscript range is " << lo 
                      << " through " << hi << " inclusive.";
            throw ModelException(complaint.str());
        }
        std::stringstream var_name_stream;
        var_name_stream << param->name << "(" << i << ")";
        return var_name_stream.str();
    }

    Function::Function(ModelInfoFromYAMLFile& for_model,
                       const char* name)
        : super(for_model), function_name(name) {
    }

    bool   Function::isFunction() const {
        return true;
    }


    UnaryFunction::UnaryFunction(ModelInfoFromYAMLFile& for_model,
                                 const char* name,
                                 const UnaryFunctionImplementation* implementation)
        : super(for_model, name)
        , body(implementation), parameter(nullptr) {
    }

    void UnaryFunction::setParameter(Expression* param) {
        parameter = param;
    }

    double UnaryFunction::evaluate() const {
        double parameter_value = parameter->evaluate();
        return body->callFunction(model, parameter_value);
    }

    void   UnaryFunction::writeTextTo(std::stringstream &text) const {
        text << function_name << "(";
        parameter->writeTextTo(text);
        text << ")";
    }

    UnaryFunction::~UnaryFunction() {
        delete parameter;
        parameter = nullptr;
    }

    InfixOperator::InfixOperator(ModelInfoFromYAMLFile& for_model )
        : super(for_model), lhs(nullptr), rhs(nullptr) {
    }

    void InfixOperator::setOperands(Expression* left, Expression* right) {
        lhs = left;
        rhs = right;
    }

    void InfixOperator::writeInfixTextTo(const char* operator_text, 
                                         std::stringstream& text) const {
        text << "(";
        lhs->writeTextTo(text);
        text << ")" << operator_text << "(";
        rhs->writeTextTo(text);
        text << ")";
    }

    bool InfixOperator::isOperator() const {
        return true;
    }

    InfixOperator::~InfixOperator() {
        delete lhs;
        lhs = nullptr;
        delete rhs;
        rhs = nullptr;
    }

    Exponentiation::Exponentiation(ModelInfoFromYAMLFile& for_model)
        : super(for_model) { }

    double Exponentiation::evaluate() const {
        double v1 = lhs->evaluate();
        double v2 = rhs->evaluate();
        return pow(v1, v2);
    }

    void  Exponentiation::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo("^", text);
    }

    int    Exponentiation::getPrecedence() const { return 12; }

    Multiplication::Multiplication(ModelInfoFromYAMLFile& for_model)
        : super(for_model) { }

    double Multiplication::evaluate() const {
        double v1 = lhs->evaluate();
        double v2 = rhs->evaluate();
        return v1 * v2;
    }

    void  Multiplication::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo("*", text);
    }

    int    Multiplication::getPrecedence() const { return 11; }

    Division::Division(ModelInfoFromYAMLFile& for_model )
        : super(for_model) { }

    double Division::evaluate() const {
        return lhs->evaluate() / rhs->evaluate();
    }

    void  Division::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo("/", text);
    }

    int    Division::getPrecedence() const { return 11; }

    Addition::Addition(ModelInfoFromYAMLFile& for_model)
        : super(for_model) { }

    double Addition::evaluate() const {
        double v1 = lhs->evaluate();
        double v2 = rhs->evaluate();
        return v1 + v2;
    }

    void  Addition::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo("+", text);
    }

    int    Addition::getPrecedence() const { return 10; }

    Subtraction::Subtraction(ModelInfoFromYAMLFile& for_model )
        : super(for_model) { }

    double Subtraction::evaluate() const {
        return lhs->evaluate() - rhs->evaluate();
    }

    void  Subtraction::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo("-", text);
    }

    int Subtraction::getPrecedence() const { return 10; }

    Assignment::Assignment(ModelInfoFromYAMLFile& for_model )
    : super(for_model) { }

    void Assignment::setOperands(Expression* left, 
                                 Expression* right) { //takes ownership
        super::setOperands(left, right);
        if (rhs==nullptr) {
            throw ModelExpression::ModelException
                  ("Assignment has null right-hand side");
        }
        if (lhs==nullptr) {
            throw ModelExpression::ModelException
                  ("Assignment has null left-hand side");
        }
        if (!lhs->isVariable()) {
            throw ModelExpression::ModelException
                  ("Can only assign to variables");
        }
    }

    double Assignment::evaluate() const {
        if (rhs==nullptr) {
            throw ModelExpression::ModelException
                  ("Assignment has null right-hand side");
        }
        double eval = rhs->evaluate();
        if (lhs==nullptr) {
            throw ModelExpression::ModelException
                  ("Assignment has null left-hand side");
        }
        if (!lhs->isVariable()) {
            throw ModelExpression::ModelException
                  ("Can only assign to variables");
        }
        Variable* v = dynamic_cast<ModelExpression::Variable*>(lhs);
        model.assign(v->getName(), eval);
        //std::cout << "VSET: " << v->getName() << "=" << eval << std::endl;
        return eval;
    }

    void  Assignment::writeTextTo(std::stringstream &text) const {
        if (lhs==nullptr || !lhs->isVariable()) {
            text << "Invalid assignment ";
            if (lhs==nullptr) {
                text << "(null)";
            } else {
                lhs->writeTextTo(text);
            }
            text << "=";
        } else {
            Variable* v = dynamic_cast<ModelExpression::Variable*>(lhs);
            text << v->getName() << "=";
        }
        rhs->writeTextTo(text);
    }

    bool Assignment::isAssignment() const {
        return true;
    }

    bool Assignment::isRightAssociative() const {
        return true;
    }

    int Assignment::getPrecedence()     const { return 9; }

    Expression* Assignment::getTarget()         const {
        return lhs;
    }

    Variable*   Assignment::getTargetVariable() const {
        return lhs->isVariable() ? dynamic_cast<Variable*>(lhs) : nullptr;
    }

    Expression* Assignment::getExpression()    const {
        return rhs;
    }

    BooleanOperator::BooleanOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    bool BooleanOperator::isBoolean() const {
        return true;
    }

    LessThanOperator::LessThanOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    int LessThanOperator::getPrecedence() const {
        return 8;
    }

    double LessThanOperator::evaluate()      const {
        return lhs->evaluate() < rhs->evaluate() ? 1.0 : 0.0;
    }

    void  LessThanOperator::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo(" < ", text);
    }

    GreaterThanOperator::GreaterThanOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    int GreaterThanOperator::getPrecedence() const {
        return 8;
    }

    double GreaterThanOperator::evaluate()      const {
        return lhs->evaluate() > rhs->evaluate() ? 1.0 : 0.0;
    }

    void  GreaterThanOperator::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo(" > ", text);
    }

    EqualityOperator::EqualityOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    int EqualityOperator::getPrecedence() const {
        return 7;
    }

    double EqualityOperator::evaluate()      const {
        return lhs->evaluate() == rhs->evaluate() ? 1.0 : 0.0;
    }

    void  EqualityOperator::writeTextTo(std::stringstream &text) const {
        lhs->writeTextTo(text);
        text << "==";
        rhs->writeTextTo(text);
    }

    InequalityOperator::InequalityOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    int InequalityOperator::getPrecedence() const {
        return 7;
    }

    double InequalityOperator::evaluate() const {
        return lhs->evaluate() != rhs->evaluate() ? 1.0 : 0.0;
    }

    void  InequalityOperator::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo("!=", text);
    }

    ShortcutAndOperator::ShortcutAndOperator(ModelInfoFromYAMLFile& for_model)
    : super(for_model) {}

    int ShortcutAndOperator::getPrecedence() const {
        return 6;
    }

    double ShortcutAndOperator::evaluate() const {
        return ( lhs->evaluate() != 0 && rhs->evaluate() !=0 ) ? 1.0 : 0.0;
    }

    void  ShortcutAndOperator::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo(" && ", text);
    }

    ShortcutOrOperator::ShortcutOrOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    int ShortcutOrOperator::getPrecedence() const {
        return 5;
    }

    double ShortcutOrOperator::evaluate() const {
        return ( lhs->evaluate() != 0 || rhs->evaluate() !=0 ) ? 1.0 : 0.0;
    }

    void  ShortcutOrOperator::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo(" || ", text);
    }

    ListOperator::ListOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {
    }

    ListOperator::ListOperator(ModelInfoFromYAMLFile& for_model, 
                               Expression*            single_param)
        : super(for_model) {
            list_entries.push_back(single_param);
        }

    ListOperator::~ListOperator() {
        for (Expression* exp : list_entries) {
            delete exp;
        }
        list_entries.clear();
    }

    int ListOperator::getPrecedence() const {
        return 4;
    }

    double ListOperator::evaluate() const {
        double v = 0;
        for (Expression* expr: list_entries) {
            v = expr->evaluate();
        }
        return v;
    }

    void  ListOperator::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo(" : ", text);
    }

    void   ListOperator::setOperands(Expression* left, Expression* right) {
        //The precedence comparison ensures that 1 : 2 , 3 
        //is treated as ( 1 : 2), 3, rather than as 4, 3, 5,
        //and that (1, 2) : (3, 4) is not treated as 1 : 2 : (3, 4).
        //
        if (left->isList() && left->getPrecedence()==this->getPrecedence()) {
            ListOperator* old = dynamic_cast<ListOperator*>(left);
            std::swap(old->list_entries, list_entries);
            list_entries.push_back(right);
            delete left;
        } else {
            list_entries.push_back(left);
            list_entries.push_back(right);
        }
    }

    bool ListOperator::isList() const {
        return true;
    }

    int  ListOperator::getEntryCount()   const {
        return static_cast<int>(list_entries.size());
    }

    double  ListOperator::evaluateEntry(int index) const {
        if (index<0) {
            std::stringstream complaint;
            complaint << "Cannot select list element"
                      <<" with zero-based index " << index << ".";
            throw ModelException(complaint.str());
        }
        if (list_entries.size()<=index) {
            std::stringstream complaint;
            complaint << "Cannot select list element"
                      << " with zero-based index " << index
                      << " from a list of " << list_entries.size() << "entries.";
            throw ModelException(complaint.str());

        }
        return list_entries[index]->evaluate();
    }

    CommaOperator::CommaOperator(ModelInfoFromYAMLFile& for_model): super(for_model) {
    }
    
    int CommaOperator::getPrecedence() const {
        return 1;
    }

    void  CommaOperator::writeTextTo(std::stringstream &text) const {
        //Doesn't use  writeInfixTextTo() because don't want
        //parentheses around all of the entries in the list.
        lhs->writeTextTo(text);
        text << ", ";
        rhs->writeTextTo(text);
    }

    MultiFunctionImplementation::MultiFunctionImplementation(int min_pars, int max_pars)
        : min_param_count(min_pars), max_param_count(max_pars) {
    }

    MultiFunction::MultiFunction(ModelInfoFromYAMLFile& for_model,
                      const char* name,
                      const MultiFunctionImplementation* implementation)
        : super(for_model, name), body(implementation) {
    }

    MultiFunction::~MultiFunction() {
        delete parameter_list;
        parameter_list = nullptr;
    }

    void   MultiFunction::setParameter(Expression* param) /*takes ownership*/ {
        if (!param->isList())
        {
            Expression* inner_param = param;
            param = new ListOperator(this->model, inner_param);              
        }
        parameter_list = dynamic_cast<ListOperator*>(param);
        if (parameter_list==nullptr) {
            delete param;
            throw ModelExpression::ModelException
                  ("Parameter could not be cast to parameter_list");
        }
        int param_count = parameter_list->getEntryCount();
        if (param_count < body->min_param_count) {
            std::stringstream complaint;
            complaint << "Too few parameters (" << param_count << ")"
                      << " for function " << function_name
                      << " (which must have at least "
                      << body->min_param_count << ")";
            throw  ModelExpression::ModelException(complaint.str());     
        }
        if (body->max_param_count < param_count) {
            std::stringstream complaint;
            complaint << "Too many parameters (" << param_count << ")"
                      << " for function " << function_name
                      << " (which can take at most "
                      << body->max_param_count << ")";
            throw  ModelExpression::ModelException(complaint.str());
        }
    }

    double MultiFunction::evaluate() const {
        return body->callFunction(model, parameter_list);
    }

    void   MultiFunction::writeTextTo(std::stringstream &text) const {
        text << function_name << "(";
        parameter_list->writeTextTo(text);
        text << ")";
    }

    SelectOperator::SelectOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    int    SelectOperator::getPrecedence() const {
        return 3;
    }

    double SelectOperator::evaluate()      const {
        if (lhs->isBoolean()) {
            if (rhs->isList()) {
                int index = (lhs->evaluate() == 0) ? 1: 0;
                ListOperator* r = dynamic_cast<ListOperator*>(rhs);
                return r->evaluateEntry(index);
            } else {
                return lhs->evaluate() ? rhs->evaluate() : 0 ;
            }
        } else {
            double index = lhs->evaluate();
            if (index<0) {
                std::stringstream complaint;
                complaint << "Cannot select list element"
                          << " with zero-based index " << index
                          << " from a list.";
                throw ModelException(complaint.str());
            }
            if (rhs->isList()) {
                ListOperator* r = dynamic_cast<ListOperator*>(rhs);
                if (r->getEntryCount() <= index) {
                    std::stringstream complaint;
                    complaint << "Cannot select list element"
                              << " with zero-based index " << index
                              << " from a list of " << r->getEntryCount()
                              << " entries.";
                    throw ModelException(complaint.str());
                }
                int int_index = (int) floor(index);
                return r->evaluateEntry(int_index);
            }
            if (index==0) {
                return 0.0;
            }
            return rhs->evaluate();
        }
    }

    void   SelectOperator::writeTextTo(std::stringstream &text) const {
        writeInfixTextTo(" ? ", text);
    }

    RangeOperator::RangeOperator(ModelInfoFromYAMLFile& for_model)
        : super(for_model) {}

    int    RangeOperator::getPrecedence()        const { return 2; }

    bool   RangeOperator::isRange()              const { return true; }

    int    RangeOperator::getIntegerLowerBound() const {
        //Todo: better error-checking (out of range for an int...
        return (int)floor(lhs->evaluate());
    }

    int    RangeOperator::getIntegerUpperBound() const {
        //Todo: better error-checking (out of range for an int...
        return (int)floor(rhs->evaluate());
    }

    void   RangeOperator::writeTextTo(std::stringstream &text) const {
        //Doesn't use  writeInfixTextTo() because don't want
        //parentheses around the lower bound, or around the
        //upper bound.
        lhs->writeTextTo(text);
        text << "..";
        rhs->writeTextTo(text);
    }

    InterpretedExpression::InterpretedExpression(ModelInfoFromYAMLFile& for_model,
                                                 const std::string& text): super(for_model) {
        is_unset  = text.empty();
        size_t ix = 0;
        root      = parseExpression(text, ix);
    }

    class ExpressionStack: public std::vector<Expression*> {
    public:
        ExpressionStack& operator << (Expression* x) {
            push_back(x);
            return *this;
        }
        Expression* pop() {
            Expression* x = back();
            pop_back();
            return x;
        }
        bool topIsOperator() const {
            return (0<size() && back()->isOperator());
        }
        bool topIsToken(char c) const {
            return (0<size() && back()->isToken(c));
        }
        bool topIsNotToken(char c) const {
            return (0<size() && !back()->isToken(c));
        }
        bool topIsFunction() const {
            return (0<size() && back()->isFunction());
        }
    };

    namespace {
        void logExpressionAndValue(const char* prefix, Expression* expr) {
            std::cout << prefix << ": ";
            std::stringstream expr_text;
            expr->writeTextTo(expr_text);
            std::cout << expr_text.str();
            try {
                auto value =  expr->evaluate();
                std::cout << " evaluated to " << value;
            }
            catch (ModelExpression::ModelException x) {
                std::cout << " could not be evaluated (" << x.getMessage() << ")";
            }
            std::cout << std::endl;
        }
    };

    Expression* InterpretedExpression::parseExpression(const std::string& text,
                                                       size_t& ix) {
        //
        //Notes: 1. Now supports unary and multiple-parameter functions
        //          (multiple-parameter functions take ONE parameter but it
        //           is a list constructed via CommaOperator or ListOperator).
        //       2. Based on the pseudocode quoted at
        //          https://en.wikipedia.org/wiki/Shunting-yard_algorithm
        //
        Expression*     token;
        ExpressionStack output;
        ExpressionStack operator_stack;
        
        while (parseToken(text, ix, token)) {
            if (token->isConstant()) {
                output << token;
            } else if (token->isFunction()) {
                operator_stack << token;
            } else if (token->isOperator()) {
                auto prec = token->getPrecedence();
                bool isRightAssociative = token->isRightAssociative();
                while (operator_stack.topIsOperator() &&
                       ( operator_stack.back()->getPrecedence() > prec ||
                        ( operator_stack.back()->getPrecedence() == prec
                          && !isRightAssociative )
                       ) ) {
                    output << operator_stack.pop();
                }
                operator_stack << token;
            } else if (token->isVariable()) {
                output << token;
            } else if (token->isToken('(')) {
                operator_stack << token;
            } else if (token->isToken(')')) {
                while (operator_stack.topIsNotToken('(')) {
                    output << operator_stack.pop();
                }
                if (operator_stack.topIsToken('(')) {
                    delete operator_stack.pop();
                }
                if (operator_stack.topIsFunction()) {
                    output << operator_stack.pop();
                }
                delete token;
            }
        } //no more tokens
        while ( !operator_stack.empty() ) {
            output << operator_stack.pop();
        }

        return parseTokenizedExpressions(output);
    }

    Expression* InterpretedExpression::parseTokenizedExpressions(ExpressionStack& tokenized) {
        ExpressionStack operand_stack;
        for (size_t i=0; i<tokenized.size(); ++i) {
            Expression* token = tokenized[i];
            if (token->isOperator()) {
                Expression*    rhs = operand_stack.pop();
                Expression*    lhs = operand_stack.pop();
                InfixOperator* op  = dynamic_cast<InfixOperator*>(token);
                op->setOperands(lhs, rhs);
                operand_stack << op;
            } else if (token->isFunction()) {
                Expression* param = operand_stack.pop();
                Function*   fn    = dynamic_cast<Function*>(token);
                if (fn==nullptr) {
                    throw ModelException("Internal Logic Error:"
                                            " Unrecognized function type.");
                }
                fn->setParameter(param);
                operand_stack << fn;
            } else {
                operand_stack << token;
            }
        }
        if (operand_stack.size()!=1) {
            std::stringstream complaint;
            complaint << "Malformed expression appeared to be " 
                      << operand_stack.size() << " space-separated expressions:";
            int parameter_number = 1;
            for (auto entry : operand_stack ) {
                complaint << "\n" << parameter_number << " is: ";
                entry->writeTextTo(complaint);
                ++parameter_number;
            }
            throw ModelException(complaint.str());
        }
        return operand_stack[0];
    }

    bool InterpretedExpression::parseVariable(const std::string& text,
                                              std::string& var_name,
                                              size_t& ix,
                                              Expression*& expr) {
        if (text[ix]=='(') {
            //Subscripted: Parse the subscript expression
            size_t subscript_start = ix;
            size_t text_length = text.length();
            int    bracket_depth = 0;
            do {
                auto ch = text[ix];
                if      (ch==')')  --bracket_depth;
                else if (ch=='(')  ++bracket_depth;
                ++ix;
            } while (ix<text_length && 0<bracket_depth);
            //ix = index of first character after the closing bracket.
            std::string subscript_expr;
            if (bracket_depth==0) {
                subscript_expr = text.substr(subscript_start+1, ix-subscript_start-2);
                if (is_string_all_digits(subscript_expr)) {
                    //std::cout << "VS expr: " << subscript_expr << std::endl;
                    var_name = var_name + "(" + subscript_expr + ")";
                } else {
                    //Oh, boy.  Subscript expression!
                    InterpretedExpression x(model, subscript_expr);
                    var_name = string_to_lower(var_name);
                    const YAMLFileParameter* param = model.findParameter(var_name) ;
                    if (param == nullptr ) {
                        std::stringstream complaint;
                        complaint << "Subscripted parameter " << var_name
                                    << " not defined, for model " << model.getLongName();
                        throw ModelException(complaint.str());
                    }
                    //std::cout << "PS expr: " << subscript_expr << std::endl;
                    expr = new ParameterSubscript(model, param, x.detatchExpression());
                    if (false) {
                        logExpressionAndValue("PS", expr);
                    }
                    return true;
                }
            } else {
                std::stringstream complaint;
                complaint << "Subscripted reference for " << var_name
                            << " for model " << model.getLongName()
                            << " was not terminated by a closing bracket:\n" << text;
                throw ModelException(complaint.str());
            }
        }
        expr = new Variable(model, var_name);
        return true;
    }

    std::string InterpretedExpression::parseIdentifier(const std::string& text,
                                                       size_t& ix) {
        size_t var_start = ix;
        for (++ix; ix<text.size(); ++ix ) {
            char ch = text[ix];
            if (!isalpha(ch) && (ch < '0' || '9' < ch)
                && ch != '.' && ch != '_' ) {
                break;
            }
        }
        return text.substr(var_start, ix-var_start);
    }

    bool InterpretedExpression::parseNumericConstant(const std::string& text,
                                                     size_t& ix,
                                                     Expression*& expr) {
        //Number
        int    endpos    = 0;
        double v         = convert_double(text.c_str()+ix, endpos);
        expr             = new Constant(model, v);
        //".." sequences are read a bit funny.  We don't want them
        //skipped over, because then expressions like "1..5" would not
        //be read as "a range with lower bound 1 and upper bound 5".
        for (; endpos>0 && (text[ix]!='.' || text[ix+1]!='.'); --endpos) {
            ++ix;
        }
        return true;
    }

    void InterpretedExpression::skipWhiteSpace(const std::string& text,
                                               size_t& ix) {
        while (ix<text.size() && text[ix]==' ') {
            ++ix;
        }
    }

    bool InterpretedExpression::parseToken(const std::string& text,
                                           size_t& ix,
                                           Expression*& expr) {
        skipWhiteSpace(text, ix);
        if (text.size()<=ix) {
            expr = nullptr;
            return false;
        }
        auto ch = tolower(text[ix]);
        if (isalpha(ch)) {
            //Variable, or function call
            std::string var_name = parseIdentifier(text, ix);
            skipWhiteSpace(text, ix);
            if (built_in_functions.isBuiltIn(var_name)) {
                expr = built_in_functions.asBuiltIn(model, var_name);
                return true;
            }
            return parseVariable(text, var_name, ix, expr);
        }
        char nextch = ((ix+1)<text.length()) ? text[ix+1] : '\0';
        if (('0'<=ch && ch<='9') || (ch=='.' && '0'<=nextch && nextch<='9')) {
            return parseNumericConstant(text, ix, expr);
        }
        return parseOtherToken(text,ix,expr);
    }

    bool InterpretedExpression::parseOtherToken(const std::string& text,
                                                size_t&            ix,
                                                Expression*&       expr) {
        char ch = text[ix];
        char nextch = ((ix+1)<text.length()) ? text[ix+1] : '\0';
        switch (ch) {
            case '(': expr = new Token(model, ch);      break;
            case ')': expr = new Token(model, ch);      break;
            case '!': if (nextch=='=') {
                          expr = new InequalityOperator(model);
                          ++ix;
                      } else {
                          throw ModelException("unary not (!) operator not supported (yet)");
                      }
                      break;
            case '^': expr = new Exponentiation(model); break;
            case '*': expr = new Multiplication(model); break;
            case '/': expr = new Division(model);       break;
            case '+': expr = new Addition(model);       break;
            case '-': expr = new Subtraction(model);    break;
            case '<': expr = new LessThanOperator(model);    break;
            case '>': expr = new GreaterThanOperator(model); break;
            case '=': if (nextch=='=') {
                          expr = new EqualityOperator(model);
                          ++ix;
                      } else {
                          expr = new Assignment(model);
                      }
                      break;
            case '&': if (nextch=='&') {
                          expr = new ShortcutAndOperator(model);
                          ++ix;
                      } else {
                          throw ModelException("bitwise-and & operator not supported");
                      }
                      break;
            case '|': if (nextch=='|') {
                          expr = new ShortcutOrOperator(model);
                          ++ix;
                      } else {
                            throw ModelException("bitwise-or | operator not supported");
                      }
                      break;
            case ':': expr = new ListOperator(model);   break;
            case '?': expr = new SelectOperator(model); break;
            case '.': if (nextch=='.') {
                          expr = new RangeOperator(model);
                          ++ix;
                      } else {
                          throw ModelException("period that wasn't part of .. or a number,"
                                                   " was not understood");
                      }
                      break;
            case ',': expr = new CommaOperator(model); break;
            default:
                throw ModelException(std::string("unrecognized character '") +
                                     std::string(1, ch) + "'' in expression");
        }
        ++ix;
        return true;
    }

    bool InterpretedExpression::isSet() const {
        return !is_unset;
    }

    InterpretedExpression::~InterpretedExpression() {
        delete root;
        root = nullptr;
    }

    double InterpretedExpression::evaluate() const {
        ASSERT(root != nullptr);
        return root->evaluate();
    }

    void InterpretedExpression::writeTextTo(std::stringstream &text) const {
        root->writeTextTo(text);
    }

    Expression* InterpretedExpression::expression() const {
        return root;
    }

    Expression* InterpretedExpression::detatchExpression()  {
        //Note: This yields ownership
        Expression* expr = root;
        root = nullptr;
        return expr;
    }

    bool InterpretedExpression::evaluateIntegerRange(std::pair<int,int>& range) const {
        if (!root->isRange()) {
            return false;
        }
        RangeOperator* expr = dynamic_cast<RangeOperator*>(root);
        const char* for_what;
        try {
            for_what = "lower bound";
            range.first  = expr->getIntegerLowerBound();
            for_what = "upper bound";
            range.second = expr->getIntegerUpperBound();
        }
        catch (ModelException x) {
            throw ModelException("Error evaluating" + std::string(for_what) + 
                                 " " + x.getMessage());
        }
        return true;
    }
}
