//
// timekeeper.cpp
// Created by James Barbetti on 27/10/20.
//

#include "timekeeper.h"
#include "timeutil.h"    //for getRealTime() and getCPUTime()
#include "progress.h"    //for appendTimeDescription
#include <iostream>
#include <math.h>        //for floor()
#include <sstream>       //for std::stringstream

TimeKeeper::TimeKeeper(const char* activity_name)
    : activity(activity_name), stopped(true)
    , elapsed_wallclock_time(0.0), elapsed_cpu_time(0.0) {
}

TimeKeeper::TimeKeeper(std::string activity_name)
    : activity(activity_name), stopped(true)
    , elapsed_wallclock_time(0.0), elapsed_cpu_time(0.0) {
}

const std::string& TimeKeeper::getActivity() const {
    return activity;
}

void TimeKeeper::setActivity(std::string& new_activity) {
    activity = new_activity;
}

void TimeKeeper::setActivity(const char*  new_activity) {
    activity = new_activity;
}

void TimeKeeper::start()  const {
    if (stopped) {
        elapsed_wallclock_time -= getRealTime();
        elapsed_cpu_time       -= getCPUTime();
        stopped = false;
    }
}

void TimeKeeper::stop()   const {
    if (!stopped) {
        elapsed_cpu_time       += getCPUTime();
        elapsed_wallclock_time += getRealTime();
        stopped = true;
    }
}

void TimeKeeper::report() const {
    std::stringstream percent;
    double elapsed  = elapsed_wallclock_time + (stopped ? 0 : getRealTime());
    double cpu_time = elapsed_cpu_time       + (stopped ? 0 : getCPUTime());
    
    std::stringstream desc;
    desc << activity << " took ";
    appendTimeDescription(elapsed, desc);
    desc << " of wall-clock time and ";
    appendTimeDescription(elapsed_cpu_time, desc);
    desc << " seconds of CPU time";
    if (0<elapsed) {
        desc.precision(4);
        desc << " (" << floor(cpu_time * 100.0 / elapsed + .5)
             << "%)";
    }
    desc << "." << std::endl;
    std::cout << desc.str();
}

std::string TimeKeeper::getElapsedDescription() const {
    double elapsed  = elapsed_wallclock_time + (stopped ? 0 : getRealTime());
    std::stringstream desc;
    appendTimeDescription(elapsed, desc);
    return desc.str();
}

TimeKeeper::~TimeKeeper() = default;

KeptTime::KeptTime(const TimeKeeper& tk): keeper(tk) {
    keeper.start();
}

KeptTime:: ~KeptTime() {
    keeper.stop();
}

