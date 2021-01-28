//
//  modelinfo.h
//  Created by James Barbetti on 28-Jan-2021
//

#ifndef modelinfo_h
#define modelinfo_h

#include <string>

class ModelInfo {
public:
    ModelInfo()                     = default;
    ModelInfo(const ModelInfo& rhs) = delete;
    virtual ~ModelInfo()            = default;
    
    virtual bool isModelFinder()             const = 0;
    virtual bool isModelFinderOnly()         const = 0;
    virtual bool isPolymorphismAware()       const = 0;
    virtual bool isWeissAndVonHaeselerTest() const = 0;
    virtual bool hasRateHeterotachy()        const = 0;
};

class ModelInfoFromName: public ModelInfo {
private:
    std::string model_name;
public:
    explicit ModelInfoFromName(std::string name);
    explicit ModelInfoFromName(const char* name);
    virtual ~ModelInfoFromName();
    
    virtual bool isModelFinder()             const;
    virtual bool isModelFinderOnly()         const;
    virtual bool isPolymorphismAware()       const;
    virtual bool isWeissAndVonHaeselerTest() const;
    virtual bool hasRateHeterotachy()        const;

    std::string extractPolymorphicHeterozygosity(std::string& leftover_name) const;
    void        updateName(const std::string& name);
};

#endif /* modelinfo_hpp */
