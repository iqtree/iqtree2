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
//     Division
//     Addition
//     Subtraction
//     Assignment
//
// Created by James Barbetti on 01-Apr-2021
//

#ifndef modelexpression_h
#define modelexpression_h

#include "modelinfofromyamlfile.h"  //for ModelInfoFromYAMLFile

namespace ModelExpression {

    class ModelException {
        protected:
            std::string message;
        public: 
            explicit ModelException(const char* s);
            explicit ModelException(const std::string& s);
            explicit ModelException(const std::stringstream& s);
            const std::string& getMessage() const;
    };

    class Expression {
    protected:
        ModelInfoFromYAMLFile& model;
    public:
        Expression(ModelInfoFromYAMLFile& for_model);
        virtual ~Expression() = default;
        virtual double evaluate()           const;
        virtual bool   isAssignment()       const;
        virtual bool   isBoolean()          const;
        virtual bool   isConstant()         const;
        virtual bool   isFunction()         const;
        virtual bool   isList()             const;
        virtual bool   isOperator()         const;
        virtual bool   isRange()            const;
        virtual bool   isRightAssociative() const;
        virtual bool   isToken(char c)      const;
        virtual bool   isVariable()         const;
        virtual int    getPrecedence()      const;
        int            evaluateAsInteger()  const;
        virtual void   writeTextTo(std::stringstream &text) const = 0;
    };

    class Token: public Expression {
    protected:
        char token_char;
    public:
        typedef Expression super;
        Token(ModelInfoFromYAMLFile& for_model, char c);
        virtual ~Token() = default;
        virtual void writeTextTo(std::stringstream &text) const;
        virtual bool isToken(char c) const;
    };

    class Variable: public Expression {
    protected:
        std::string variable_name;
        explicit Variable(ModelInfoFromYAMLFile& for_model);
        //The single parameter constructor should only be called
        //from subclasses (notably ParameterSubscript).
    public:
        typedef Expression super;
        Variable(ModelInfoFromYAMLFile& for_model,
                 const std::string& name);
        virtual ~Variable() = default;
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual bool   isVariable() const;
        const std::string& getName() const;
    };

    class ParameterSubscript: public Variable {
        std::string parameter_to_subscript;
        Expression* subscript_expression;
    public:
        typedef Variable super;
        ParameterSubscript(ModelInfoFromYAMLFile& for_model,
                            const YAMLFileParameter* param,
                            Expression* subscript_expr);
        virtual ~ParameterSubscript();
        virtual double evaluate() const;
        std::string getName() const;
    };

    class Constant: public Expression {
        double value;
    public:
        typedef Expression super;
        Constant(ModelInfoFromYAMLFile& for_model,
                 double v);
        virtual ~Constant() = default;
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual bool   isConstant() const;
    };

    class UnaryFunctionImplementation {
    public:
        virtual double callFunction(ModelInfoFromYAMLFile&,
                                    double parameter) const = 0;
        virtual ~UnaryFunctionImplementation() = default;
    };

    class UnaryFunction: public Expression {
        std::string                        function_name;
        const UnaryFunctionImplementation* body;
        Expression*                        parameter;
    public:
        typedef Expression super;
        UnaryFunction(ModelInfoFromYAMLFile& for_model,
                      const char* name,
                      const UnaryFunctionImplementation* implementation);
        virtual ~UnaryFunction();
        virtual void   setParameter(Expression* param); //takes ownership
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual bool   isFunction() const;
    };

    class InfixOperator: public Expression {
    protected:
        Expression* lhs;
        Expression* rhs;
    public:
        typedef Expression super;
        InfixOperator(ModelInfoFromYAMLFile& for_model);
        virtual ~InfixOperator();
        virtual bool isOperator() const;
        virtual void setOperands(Expression* left, Expression* right); //takes ownership
        virtual void writeInfixTextTo(const char* operator_text, std::stringstream& text) const;
    };

    class Exponentiation: public InfixOperator {
        public:
            typedef InfixOperator super;
            Exponentiation(ModelInfoFromYAMLFile& for_model);
            virtual ~Exponentiation() = default;
            virtual double evaluate() const;
            virtual void   writeTextTo(std::stringstream &text) const;
            virtual int    getPrecedence() const;
    };

    class Multiplication: public InfixOperator {
    public:
        typedef InfixOperator super;
        Multiplication(ModelInfoFromYAMLFile& for_model);
        virtual ~Multiplication() = default;
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual int    getPrecedence() const;
    };

    class Division: public InfixOperator {
    public:
        typedef InfixOperator super;
        Division(ModelInfoFromYAMLFile& for_model);
        virtual ~Division() = default;
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual int    getPrecedence() const;
    };

    class Addition: public InfixOperator {
    public:
        typedef InfixOperator super;
        Addition(ModelInfoFromYAMLFile& for_model );
        virtual ~Addition() = default;
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual int    getPrecedence() const;
    };

    class Subtraction: public InfixOperator {
    public:
        typedef InfixOperator super;
        Subtraction(ModelInfoFromYAMLFile& for_model);
        virtual ~Subtraction() = default;
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual int    getPrecedence() const;
    };

    class Assignment: public InfixOperator {
    public:
        typedef InfixOperator super;
        Assignment(ModelInfoFromYAMLFile& for_model);
        virtual ~Assignment() = default;
        virtual double evaluate()           const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual bool   isAssignment()       const;
        virtual bool   isRightAssociative() const;

        virtual int    getPrecedence()      const;
        Expression*    getTarget()          const;
        Variable*      getTargetVariable()  const;
        Expression*    getExpression()      const;
    };

    class BooleanOperator: public InfixOperator {
    public:
        typedef InfixOperator super;
        BooleanOperator(ModelInfoFromYAMLFile& for_model);
        virtual bool isBoolean() const;
    };

    class LessThanOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        LessThanOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const;
        virtual double evaluate()      const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class GreaterThanOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        GreaterThanOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const;
        virtual double evaluate()      const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class EqualityOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        EqualityOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const;
        virtual double evaluate()      const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class InequalityOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        InequalityOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const;
        virtual double evaluate()      const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class ShortcutAndOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        ShortcutAndOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const;
        virtual double evaluate()      const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class ShortcutOrOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        ShortcutOrOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const;
        virtual double evaluate()      const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class ListOperator: public InfixOperator {
        std::vector<Expression*> list_entries;
    public:
        typedef InfixOperator super;
        ListOperator(ModelInfoFromYAMLFile& for_model);
        virtual ~ListOperator();
        virtual int    getPrecedence()   const;
        virtual double evaluate()        const;
        virtual void   writeTextTo(std::stringstream &text) const;
        virtual void   setOperands(Expression* left, Expression* right); //takes ownership
        virtual bool   isList()          const;
        virtual int    getEntryCount()   const;
        virtual double evaluateEntry(int index) const;
    };

    class CommaOperator: public ListOperator {
    public:
        typedef ListOperator super;
        CommaOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence()   const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class SelectOperator: public InfixOperator {
    public:
        typedef InfixOperator super;
        SelectOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const;
        virtual double evaluate()      const;
        virtual void   writeTextTo(std::stringstream &text) const;
    };

    class RangeOperator: public InfixOperator {
    public:
        typedef InfixOperator super;
        RangeOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence()        const;
        virtual bool   isRange()              const;
        virtual int    getIntegerLowerBound() const;
        virtual int    getIntegerUpperBound() const;
        virtual void   writeTextTo(std::stringstream &text) const;
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
        InterpretedExpression(ModelInfoFromYAMLFile& for_model,
                              const std::string& expression_text);
        virtual ~InterpretedExpression();
        bool    isSet() const;
        virtual double evaluate() const;
        virtual void   writeTextTo(std::stringstream &text) const;
        Expression* expression() const; //Does *not* yield ownership
        Expression* detatchExpression(); //*Does* yield ownership
    };
} //ModelExpression namespace

#endif /* modelexpression_hpp */
