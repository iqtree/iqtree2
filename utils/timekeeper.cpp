//
//  timekeeper.cpp
//  alignment
//
//  Created by James Barbetti on 27/10/20.
//

#include "timekeeper.h"
#include "timeutil.h"    //for getRealTime() and getCPUTime()
#include <iostream>
#include <math.h>        //for floor()
#include <sstream>       //for std::stringstream

TimeKeeper::TimeKeeper(const char* activity_name)
    : activity(activity_name)
    , elapsed_wallclock_time(0.0), elapsed_cpu_time(0.0) {
}

TimeKeeper::TimeKeeper(std::string activity_name)
    : activity(activity_name)
    , elapsed_wallclock_time(0.0), elapsed_cpu_time(0.0) {
}

void TimeKeeper::start()  const {
    elapsed_wallclock_time -= getRealTime();
    elapsed_cpu_time       -= getCPUTime();
}

void TimeKeeper::stop()   const {
    elapsed_cpu_time       += getCPUTime();
    elapsed_wallclock_time += getRealTime();
}

void TimeKeeper::report() const {
    std::stringstream percent;
    if (0<elapsed_wallclock_time) {
        percent.precision(4);
        percent << "(" << floor(elapsed_cpu_time * 100.0 / elapsed_wallclock_time + .5)
                << "%)";
    }
    std::cout << activity << " took " << elapsed_wallclock_time
        << " seconds of wall-clock time and " << elapsed_cpu_time
        << " seconds of CPU time" << percent.str() << "." << std::endl;
}

TimeKeeper::~TimeKeeper() = default;

KeptTime::KeptTime(const TimeKeeper& tk): keeper(tk) {
    keeper.start();
}

KeptTime:: ~KeptTime() {
    keeper.stop();
}

