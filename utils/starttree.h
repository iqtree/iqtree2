//
//  starttree.h
//  iqtree
//
//  Created by James Barbetti on 26/6/20.
//

#ifndef starttree_h
#define starttree_h

#include <string>
#include <map>
#include <iostream>
#include "timeutil.h"       //for getRealTime()

namespace StartTree
{
    class BuilderInterface
    {
    public:
        BuilderInterface() {}
        virtual ~BuilderInterface() {}
        virtual void constructTree
            ( const std::string &distanceMatrixFilePath
             , const std::string & newickTreeFilePath) = 0;
        virtual const std::string& getName() = 0;
        virtual const std::string& getDescription() = 0;
    };

    class Factory
    {
    private:
        std::map<std::string, BuilderInterface*> mapOfTreeBuilders;
        //Note: Owned by the Factory, and will be deleted in ~Factory.
        
    protected:
        Factory();
        ~Factory();

        size_t getBuilderCount();
        void addBuilder(const std::string& name, BuilderInterface* builder);
        BuilderInterface* getBuilder(const std::string& name);
        BuilderInterface* getBuilder(const char* name);
    public:
        static Factory& getInstance();
        void   advertiseTreeBuilder(BuilderInterface* builder);
        static BuilderInterface* getTreeBuilderByName(const std::string& name);
    };

    template <class B> class Builder: public BuilderInterface
    {
        //Note: B must have:
        //      1. a constructor that takes the name of an ".mldist"
        //         distance matrix file as a parameter;
        //      2. a constructTree() member function; and
        //      3. a writeTreeFile() membr function.
        //
    protected:
        const std::string name;
        const std::string description;
    public:
        Builder(const char* nameToUse, const char *descriptionToGive)
        : name(nameToUse), description(descriptionToGive) {
        }
        virtual const std::string& getName() {
            return name;
        }
        virtual const std::string& getDescription() {
            return description;
        }
        virtual void constructTree
            ( const std::string &distanceMatrixFilePath
             , const std::string & newickTreeFilePath) {
                B builder(distanceMatrixFilePath);
                double buildStart = getRealTime();
                builder.constructTree();
                double buildElapsed = getRealTime() - buildStart;
                std::cout.precision(6);
                std::cout << "Elapsed time for constructing initial tree"
                    << " (with algorithm " << name << "), "
                    << buildElapsed << std::endl;
                std::cout.precision(3);
                builder.writeTreeFile(newickTreeFilePath);
        }
    };
}

#define START_TREE_RECOGNIZED(name) \
    ( StartTree::Factory::getTreeBuilderByName(name) != nullptr )

#endif /* starttree_h */
