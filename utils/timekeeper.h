//
//  timekeeper.h
//  alignment
//
//  Created by James Barbetti on 27/10/20.
//

#ifndef timekeeper_h
#define timekeeper_h

#include <string>

class TimeKeeper {
public:
    mutable std::string activity;
    mutable double elapsed_wallclock_time;
    mutable double elapsed_cpu_time;
    TimeKeeper(const char* activity_name);
    ~TimeKeeper();
    void start()  const;
    void stop()   const;
    void report() const;
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
