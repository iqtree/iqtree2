//
// modelexpression.h
// Implementations of classes used for analyzing expressions
// (e.g. for values of rate matrix entries) in YAML model files.
//
// ModelExpression
//   Expression
//   Variable
//   InfixOperator
//     Multiplication
//     Subtraction
//
// Created by James Barbetti on 01-Apr-2021
//

#ifndef modelexpression_h
#define modelexpression_h

#include "modelinfo.h"  //for ModelInfoFromYAMLFile

namespace ModelExpression {

    class Expression {
    protected:
        const ModelInfoFromYAMLFile& model;
    public:
        Expression(const ModelInfoFromYAMLFile& for_model);
        virtual ~Expression() = default;
        virtual double evaluate()      const;
        virtual bool   isConstant()    const;
        virtual bool   isFunction()    const;
        virtual bool   isOperator()    const;
        virtual bool   isToken(char c) const;
        virtual bool   isVariable()    const;
        virtual int    getPrecedence() const;
    };

    class Token: public Expression {
    protected:
        char token_char;
    public:
        typedef Expression super;
        Token(const ModelInfoFromYAMLFile& for_model, char c);
        virtual ~Token() = default;
        virtual bool isToken(char c) const;
    };

    class Variable: public Expression {
        std::string variable_name;
    public:
        typedef Expression super;
        Variable(const ModelInfoFromYAMLFile& for_model,
                 const std::string& name);
        virtual ~Variable() = default;
        virtual double evaluate() const;
        virtual bool isVariable() const;
    };

    class Constant: public Expression {
        double value;
    public:
        typedef Expression super;
        Constant(const ModelInfoFromYAMLFile& for_model,
                 double v);
        virtual ~Constant() = default;
        virtual double evaluate() const;
        virtual bool isConstant() const;
    };

    class UnaryFunctionImplementation {
    public:
        virtual double callFunction(const  ModelInfoFromYAMLFile&,
                                    double parameter) const = 0;
        virtual ~UnaryFunctionImplementation() = default;
    };

    class UnaryFunction: public Expression {
        const UnaryFunctionImplementation* body;
        Expression*                        parameter;
    public:
        typedef Expression super;
        UnaryFunction(const ModelInfoFromYAMLFile& for_model,
                      const UnaryFunctionImplementation* implementation);
        virtual ~UnaryFunction();
        virtual void   setParameter(Expression* param); //takes ownership
        virtual double evaluate() const;
        virtual bool   isFunction() const;
    };

    class InfixOperator: public Expression {
    protected:
        Expression* lhs;
        Expression* rhs;
    public:
        typedef Expression super;
        InfixOperator(const ModelInfoFromYAMLFile& for_model);
        virtual ~InfixOperator();
        virtual bool isOperator() const;
        void setOperands(Expression* left, Expression* right); //takes ownership
    };

    class Exponentiation: public InfixOperator {
        public:
            typedef InfixOperator super;
            Exponentiation(const ModelInfoFromYAMLFile& for_model);
            virtual ~Exponentiation() = default;
            virtual double evaluate() const;
            virtual int    getPrecedence() const;
    };

    class Multiplication: public InfixOperator {
    public:
        typedef InfixOperator super;
        Multiplication(const ModelInfoFromYAMLFile& for_model);
        virtual ~Multiplication() = default;
        virtual double evaluate() const;
        virtual int    getPrecedence() const;
    };

    class Division: public InfixOperator {
    public:
        typedef InfixOperator super;
        Division(const ModelInfoFromYAMLFile& for_model);
        virtual ~Division() = default;
        virtual double evaluate() const;
        virtual int    getPrecedence() const;
    };

    class Addition: public InfixOperator {
    public:
        typedef InfixOperator super;
        Addition(const ModelInfoFromYAMLFile& for_model );
        virtual ~Addition() = default;
        virtual double evaluate() const;
        virtual int    getPrecedence() const;
    };

    class Subtraction: public InfixOperator {
    public:
        typedef InfixOperator super;
        Subtraction(const ModelInfoFromYAMLFile& for_model);
        virtual ~Subtraction() = default;
        virtual double evaluate() const;
        virtual int    getPrecedence() const;
    };

    class InterpretedExpression: public Expression {
    protected:
        bool        is_unset;
        Expression* root;
        Expression* parseExpression(const   std::string& expression_text,
                                    size_t& index);
        bool        parseToken     (const std::string& text,
                                    size_t& ix, Expression*& expr);
    public:
        typedef Expression super;
        using super::model;
        InterpretedExpression(const ModelInfoFromYAMLFile& for_model,
                              const std::string& expression_text);
        virtual ~InterpretedExpression();
        bool    isSet() const;
        virtual double evaluate() const;
    };
} //ModelExpression namespace

#endif /* modelexpression_hpp */
