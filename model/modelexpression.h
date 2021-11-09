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
            static void throwIfNonBlank(const std::stringstream& s);
    };

    class Expression {
    protected:
        ModelInfoFromYAMLFile& model;
    public:
        Expression(const Expression& rhs) = delete;
        explicit Expression(ModelInfoFromYAMLFile& for_model);
        virtual ~Expression() = default;
        virtual double evaluate()           const;
        virtual bool   isAssignment()       const;
        virtual bool   isBoolean()          const;
        virtual bool   isConstant()         const;
        virtual bool   isEstimate()         const;
        virtual bool   isFixed()            const = 0;
        virtual bool   isFunction()         const;
        virtual bool   isInitializedList()  const;
        virtual bool   isList()             const;
        virtual bool   isOperator()         const;
        virtual bool   isRange()            const;
        virtual bool   isRightAssociative() const;
        virtual bool   isToken(char c)      const;
        virtual bool   isVariable()         const;
        virtual int    getPrecedence()      const;
        int            evaluateAsInteger()  const;
        virtual void   writeTextTo(std::stringstream &text) const = 0;
        std::string    getText()            const;
    };

    class Token: public Expression {
    protected:
        char token_char;
    public:
        typedef Expression super;
        Token(ModelInfoFromYAMLFile& for_model, char c);
        virtual ~Token() = default;
        virtual void writeTextTo(std::stringstream &text) const override;
        virtual bool isFixed() const override;
        virtual bool isToken(char c) const override;
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
        virtual double evaluate()    const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual bool   isFixed()     const override;
        virtual bool   isVariable()  const override;
        const std::string& getName() const;
    };

    class ParameterSubscript: public Variable {
    protected:
        std::string parameter_to_subscript;
        Expression* subscript_expression;
        int checkSubscript(const YAMLFileParameter* param,
                           double x) const;
        const YAMLFileParameter* checkParameter()         const;

    public:
        typedef Variable super;
        ParameterSubscript(ModelInfoFromYAMLFile& for_model,
                            const YAMLFileParameter* param,
                            Expression* subscript_expr);
        virtual ~ParameterSubscript();
        virtual double evaluate()   const override;
        virtual bool   isFixed()    const override;
        std::string    getName()    const;
    };

    class Constant: public Expression {
        double value;
    public:
        typedef Expression super;
        Constant(ModelInfoFromYAMLFile& for_model,
                 double v);
        virtual ~Constant() = default;
        virtual double evaluate()   const override;
        virtual bool   isConstant() const override;
        virtual bool   isFixed()    const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

 
    class UnaryFunctionImplementation {
    public:
        virtual double callFunction(ModelInfoFromYAMLFile&,
                                    double parameter) const = 0;
        virtual ~UnaryFunctionImplementation() = default;
    };

    class Function: public Expression {
    protected:
        std::string    function_name;
    public:
        typedef Expression super;
        Function() = delete;
        Function(ModelInfoFromYAMLFile& for_model,
                      const char* name);
        virtual ~Function() = default;

        virtual bool   isFunction() const override;
        virtual void   setParameter(Expression* param) = 0; //takes ownership
    };

   class Estimate: public Function {
    protected:
        Expression* expression;
    public:
        typedef Function super;
        explicit Estimate(ModelInfoFromYAMLFile& for_model);
        virtual ~Estimate();
        virtual double evaluate() const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual bool   isConstant() const override;
        virtual bool   isEstimate() const override;
        virtual bool   isFixed()    const override;
        virtual void   setParameter(Expression* param) override; //takes ownership
    };

    class UnaryFunction: public Function {
        const UnaryFunctionImplementation* body;
        Expression*                        parameter;
    public:
        typedef Function super;
        UnaryFunction(ModelInfoFromYAMLFile& for_model,
                      const char* name,
                      const UnaryFunctionImplementation* implementation);
        virtual ~UnaryFunction();

        virtual void   setParameter(Expression* param) override; //takes ownership
        virtual double evaluate() const override;
        virtual bool   isFixed()  const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class InfixOperator: public Expression {
    protected:
        Expression* lhs;
        Expression* rhs;
    public:
        typedef Expression super;
        explicit InfixOperator(ModelInfoFromYAMLFile& for_model);
        virtual ~InfixOperator();
        virtual bool   isFixed()    const override;
        virtual bool   isOperator() const override;
        virtual void   setOperands(Expression* left, 
                                   Expression* right); //takes ownership
        virtual void   writeInfixTextTo(const char* operator_text, 
                                        std::stringstream& text) const;
    };

    class Exponentiation: public InfixOperator {
        public:
            typedef InfixOperator super;
            explicit Exponentiation(ModelInfoFromYAMLFile& for_model);
            virtual ~Exponentiation() = default;
            virtual double evaluate() const override;
            virtual void   writeTextTo(std::stringstream &text) const override;
            virtual int    getPrecedence() const override;
    };

    class Multiplication: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit Multiplication(ModelInfoFromYAMLFile& for_model);
        virtual ~Multiplication() = default;
        virtual double evaluate() const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual int    getPrecedence() const override;
    };

    class Division: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit Division(ModelInfoFromYAMLFile& for_model);
        virtual ~Division() = default;
        virtual double evaluate() const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual int    getPrecedence() const override;
    };

    class Addition: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit Addition(ModelInfoFromYAMLFile& for_model );
        virtual ~Addition() = default;
        virtual double evaluate() const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual int    getPrecedence() const override;
    };

    class Subtraction: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit Subtraction(ModelInfoFromYAMLFile& for_model);
        virtual ~Subtraction() = default;
        virtual double evaluate() const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual int    getPrecedence() const override;
    };

    class Assignment: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit Assignment(ModelInfoFromYAMLFile& for_model);
        virtual ~Assignment() = default;
        virtual void   setOperands(Expression* left, 
                                   Expression* right) override; //takes ownership


        virtual double evaluate()           const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual bool   isAssignment()       const override;
        virtual bool   isFixed()            const override;
        virtual bool   isRightAssociative() const override;
        virtual int    getPrecedence()      const override;
        Expression*    getTarget()          const;
        Variable*      getTargetVariable()  const;
        Expression*    getExpression()      const;
    protected:
        void    checkAssignment(Expression* left, 
                                Expression* right);
        /* called from setOperands */
        double  doAssignment(Expression* left,
                             Expression *right) const;
        /* called from evaluate */

    };

    class BooleanOperator: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit BooleanOperator(ModelInfoFromYAMLFile& for_model);
        virtual bool isBoolean() const override;
    };

    class LessThanOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        explicit LessThanOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual double evaluate()      const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class GreaterThanOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        explicit GreaterThanOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual double evaluate()      const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class EqualityOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        explicit EqualityOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual double evaluate()      const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class InequalityOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        explicit InequalityOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual double evaluate()      const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class ShortcutAndOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        explicit ShortcutAndOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual double evaluate()      const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class ShortcutOrOperator: public BooleanOperator {
    public:
        typedef BooleanOperator super;
        explicit ShortcutOrOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual double evaluate()      const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class ListOperator: public InfixOperator {
    protected:
        std::vector<Expression*> list_entries;
    public:
        typedef InfixOperator super;
        explicit ListOperator(ModelInfoFromYAMLFile& for_model);
        ListOperator(ModelInfoFromYAMLFile& for_model, 
                     Expression* only); //1-entry list (takes ownership)
        virtual ~ListOperator();
        virtual int    getPrecedence()     const override;
        virtual double evaluate()          const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        virtual void   setOperands(Expression* left, Expression* right) override; //takes ownership
        virtual void   addOperand(Expression* op); //takes ownership
        virtual bool   isFixed()           const override;
        virtual bool   isInitializedList() const override;
        virtual bool   isList()            const override;
        virtual int    getEntryCount()     const;
        virtual double evaluateEntry(int index) const;
        virtual Expression* getEntryExpression(int index) const;
        std::string    getEntryExpressionText(int ix) const;

    };

    class CommaOperator: public ListOperator {
    public:
        typedef ListOperator super;
        explicit CommaOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class MultiFunction; //So it can be friend of MultiFunctionImplementation

    class MultiFunctionImplementation {
    protected:
        int min_param_count;
        int max_param_count;
    public:
        friend class MultiFunction;

        MultiFunctionImplementation(int min_params, int max_params);
        ~MultiFunctionImplementation() = default;
        virtual double callFunction(ModelInfoFromYAMLFile& model,
                                    ListOperator* parameter_list) const = 0;
    };

    class MultiFunction: public Function {
        const MultiFunctionImplementation* body;
        ListOperator*                      parameter_list;
    public:
        typedef Function super;
        MultiFunction(ModelInfoFromYAMLFile& for_model,
                      const char* name,
                      const MultiFunctionImplementation* implementation);
        virtual ~MultiFunction();
        virtual void   setParameter(Expression* param) override; //takes ownership
        virtual double evaluate() const override;
        virtual bool   isFixed()  const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class SelectOperator: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit SelectOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence() const override;
        virtual double evaluate()      const  override;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class RangeOperator: public InfixOperator {
    public:
        typedef InfixOperator super;
        explicit RangeOperator(ModelInfoFromYAMLFile& for_model);
        virtual int    getPrecedence()     const override;
        virtual bool   isRange()           const override;
        virtual int    getIntegerMinimum() const;
        virtual int    getIntegerMaximum() const;
        virtual double getMinimum()        const;
        virtual double getMaximum()        const;
        virtual void   writeTextTo(std::stringstream &text) const override;
    };

    class ExpressionStack;

    class InterpretedExpression: public Expression {
    protected:
        bool        is_unset;
        Expression* root;
        Expression* parseExpression          (const   std::string& expression_text,
                                              size_t& index);                        
        Expression* parseTokenizedExpressions(ExpressionStack& tokenized);
        void        checkOperandStack        (ExpressionStack& operand_stack);
        bool        parseToken               (const std::string& text,
                                              size_t& ix, Expression*& expr);
        bool        parseVariable            (const std::string& text,
                                              std::string& var_name,
                                              size_t& ix, Expression*& expr);
        std::string parseIdentifier          (const std::string& text,
                                              size_t& ix);
        bool        parseNumericConstant     (const std::string& text,
                                              size_t& ix,
                                              Expression*& expr);
        bool        parseTwoCharacterToken   (const std::string& text,
                                              size_t& ix, Expression*& expr);                                    
        bool        parseOtherToken          (const std::string& text,
                                              size_t& ix, Expression*& expr);
        void        skipWhiteSpace           (const std::string& text,
                                              size_t& ix);

    public:
        typedef Expression super;
        using super::model;
        InterpretedExpression(ModelInfoFromYAMLFile& for_model,
                              const std::string& expression_text);
        virtual ~InterpretedExpression();
        bool           isSet()      const;
        bool           isFixed()    const override;
        virtual double evaluate()   const override;
        virtual void   writeTextTo(std::stringstream &text) const override;
        Expression*    expression() const; //Does *not* yield ownership
        Expression*    detatchExpression(); //*Does* yield ownership
        bool           evaluateIntegerRange(std::pair<int,int>& range) const;
        bool           evaluatesToAList    () const;     
        void           convertToListLookup (int ix);
    };
} //ModelExpression namespace

#endif /* modelexpression_hpp */
