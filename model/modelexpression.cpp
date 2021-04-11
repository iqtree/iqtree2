//
// modelexpression.cpp
// Implementations of classes used for analyzing expressions
// (e.g. for values of rate matrix entries) in YAML model files.
// Created by James Barbetti on 01-Apr-2021
//

#include "modelexpression.h"
#include <utils/stringfunctions.h> //for convert_double

namespace ModelExpression {

    Expression::Expression(const ModelInfoFromYAMLFile& for_model) : model(for_model) {
    }

    double Expression::evaluate()      const { return 0; }
    bool   Expression::isConstant()    const { return false; }
    bool   Expression::isFunction()    const { return false; }
    bool   Expression::isOperator()    const { return false; }
    bool   Expression::isToken(char c) const { return false; }
    bool   Expression::isVariable()    const { return false; }
    int    Expression::getPrecedence() const { return false; }

    Constant::Constant(const ModelInfoFromYAMLFile& for_model,
                       double v) : super(for_model), value(v) {
    }

    double Constant::evaluate() const {
        return value;
    }

    Variable::Variable(const ModelInfoFromYAMLFile& for_model, const std::string& name)
        : super(for_model), variable_name(name) {
        if (!for_model.hasVariable(name)) {
            std::stringstream complaint;
            complaint << "Could not evaluate variable " << name
                      << " for model " << for_model.getLongName();
        }
    }

    double Variable::evaluate() const {
        return model.getVariableValue(variable_name);
    }

    InfixOperator::InfixOperator(const ModelInfoFromYAMLFile& for_model )
        : super(for_model), lhs(nullptr), rhs(nullptr) {
    }

    void InfixOperator::setOperands(Expression* left, Expression* right) {
        lhs = left;
        rhs = right;
    }

    InfixOperator::~InfixOperator() {
        delete lhs;
        lhs = nullptr;
        delete rhs;
        rhs = nullptr;
    }

    Multiplication::Multiplication(const ModelInfoFromYAMLFile& for_model)
        : super(for_model) { }

    double Multiplication::evaluate() const {
        return lhs->evaluate() * rhs->evaluate();
    }
    int    Multiplication::getPrecedence() const { return 2; }


    Division::Division(const ModelInfoFromYAMLFile& for_model )
        : super(for_model) { }

    double Division::evaluate() const {
        return lhs->evaluate() / rhs->evaluate();
    }
    int    Division::getPrecedence() const { return 2; }

    Addition::Addition(const ModelInfoFromYAMLFile& for_model)
        : super(for_model) { }

    double Addition::evaluate() const {
        return lhs->evaluate() + rhs->evaluate();
    }
    int    Addition::getPrecedence() const { return 1; }

    Subtraction::Subtraction(const ModelInfoFromYAMLFile& for_model )
        : super(for_model) { }

    double Subtraction::evaluate() const {
        return lhs->evaluate() - rhs->evaluate();
    }

    InterpretedExpression::InterpretedExpression(const ModelInfoFromYAMLFile& for_model,
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

    Expression* InterpretedExpression::parseExpression(const std::string& text,
                                                       size_t& ix) {
        //
        //Notes: 1. As yet, only supports unary functions
        //       2. Based on the pseudocode quoted at
        //          https://en.wikipedia.org/wiki/Shunting-yard_algorithm
        //
        Expression* token;
        ExpressionStack output;
        ExpressionStack operator_stack;
        
        while (parseToken(text, ix, token)) {
            if (token->isConstant()) {
                output << token;
            } else if (token->isFunction()) {
                operator_stack << token;
            } else if (token->isOperator()) {
                while (operator_stack.topIsOperator() &&
                       operator_stack.back()->getPrecedence() >
                       token->getPrecedence()) {
                    output << operator_stack.pop();
                }
                operator_stack << token;
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
        ExpressionStack operand_stack;
        for (size_t i=0; i<output.size(); ++i) {
            token = output[i];
            if (token->isOperator()) {
                Expression* rhs = operand_stack.pop();
                Expression* lhs = operand_stack.pop();
                InfixOperator* op = dynamic_cast<InfixOperator*>(token);
                op->setOperands(lhs, rhs);
                operand_stack << op;
            } else {
                operand_stack << token;
            }
        }
        ASSERT(operand_stack.size()==1);
        return operand_stack[0];
    }

    bool InterpretedExpression::parseToken(const std::string& text,
                                                  size_t& ix,
                                                  Expression*& expr) {
        while (ix<text.size() && text[ix]==' ') {
            ++ix;
        }
        if (text.size()<=ix) {
            expr = nullptr;
            return false;
        }
        auto ch = tolower(text[ix]);
        if (isalpha(ch)) {
            //Variable
            size_t var_start = ix;
            for (++ix; ix<text.size(); ++ix ) {
                ch = text[ix];
                if (!isalpha(ch) && !isnumber(ch)) {
                    break;
                }
            }
            expr = new Variable(model, text.substr(var_start, ix-var_start));
            return true;
        }
        if (isnumber(ch)) {
            //Number
            int    endpos    = 0;
            double v         = convert_double(text.c_str()+ix, endpos);
            expr             = new Constant(model, v);
            ix              += endpos;
            return true;
        }
        switch (ch) {
            case '(': expr = new Token(model, ch);
            case ')': expr = new Token(model, ch);
            case '-': expr = new Subtraction(model);
            case '+': expr = new Addition(model);
            case '*': expr = new Multiplication(model);
            case '/': expr = new Division(model);
            default:
                throw std::string("unrecognized character '") +
                      std::string(1, ch) + "'' in expression";
        }
        return true;
    }

    bool InterpretedExpression::isSet() const {
        return !is_unset;
    }

    InterpretedExpression::~InterpretedExpression() {
        delete root;
    }

    double InterpretedExpression::evaluate() const {
        return root->evaluate();
    }
}
