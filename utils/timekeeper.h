//
// timekeeper.h
// Created by James Barbetti on 27-Oct-2020.
//

#ifndef timekeeper_h
#define timekeeper_h

#include <string>

class TimeKeeper {
protected:
    mutable std::string activity;
    mutable bool        stopped;
    mutable double      elapsed_wallclock_time;
    mutable double      elapsed_cpu_time;
    
public:
    explicit TimeKeeper(const char* activity_name);
    explicit TimeKeeper(std::string name);
    ~TimeKeeper();

    const std::string& getActivity() const;
    void start()       const;
    void stop()        const;
    void report()      const;
    std::string getElapsedDescription() const;

    void setActivity(std::string& new_activity);
    void setActivity(const char*  new_activity);

    
};

class KeptTime {
public:
    const TimeKeeper& keeper;
    KeptTime() = delete;
    KeptTime(const KeptTime& rhs) = delete;
    KeptTime(const TimeKeeper& tk);
    ~KeptTime();
};

#endif /* timekeepr_hpp */
