//
// progress.cpp - Displaying colorized progress bars on the terminal
//                (based on Boost's old progress_display class).
//
//  Copyright (C) 2020, James Barbetti.
//
//  LICENSE:
//* This program is free software; you can redistribute it and/or modify
//* it under the terms of the GNU General Public License as published by
//* the Free Software Foundation; either version 2 of the License, or
//* (at your option) any later version.
//*
//* This program is distributed in the hope that it will be useful,
//* but WITHOUT ANY WARRANTY; without even the implied warranty of
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//* GNU General Public License for more details.
//*
//* You should have received a copy of the GNU General Public License
//* along with this program; if not, write to the
//* Free Software Foundation, Inc.,
//* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//

#include "progress.h"
#include <sstream>  //for std::ostringstream
#include <iostream> //for std::cout
#include <math.h>   //for floor
#include "operatingsystem.h" //for isStandardOutputATerminal


namespace {
bool displayingProgress = true;
    //You can turn off progress displays via progress_display::setProgressDisplay.
bool isTerminal = false;
}

progress_display::progress_display( double workToDo, const char* doingWhat
                                   , const char* verb, const char* unitName)
    : startTime(getRealTime()), startCPUTime(getCPUTime())
    , totalWorkToDo(workToDo),  workDone(0.0)
    , taskDescription(doingWhat), isDone(false)
    , workVerb(verb),            workUnitName(unitName) {
        lastReportedWork    = 0.0;
        lastReportedTime    = startTime;
        lastReportedCPUTime = startCPUTime;
}

progress_display & progress_display::operator ++ () {
    return (*this) += 1.0;
}

progress_display& progress_display::operator =  (double workDoneNow) {
    double increment;
    {
        #if _OPENMP
        #pragma omp critical (io)
        #endif
        increment = workDoneNow - workDone;
    }
    if (0<increment) {
        (*this) += increment;
    }
    return *this;
}

progress_display & progress_display::operator += (double incrementalWork) {
    if (incrementalWork==0.0) {
        return *this;
    }
    double time = getRealTime();
    double cpu  = getCPUTime();
    {
        #if _OPENMP
        #pragma omp critical (io)
        #endif
        workDone += incrementalWork;
    }
    if ( lastReportedWork == 0.0 || time - lastReportedTime > 1.0 ) {
        reportProgress(time, cpu, false);
    }
    return *this;
}

void progress_display::reportProgress(double time, double cpu, bool newline) {
    double elapsedTime = time - startTime;
    std::ostringstream progress;
    if (!taskDescription.empty()) {
        progress << taskDescription << ":";
    }
    bool verbed = false;
    if (totalWorkToDo <= workDone) {
        if (!taskDescription.empty()) {
            progress << " done";
        }
        verbed = true;
    } else if (!workVerb.empty() && !workUnitName.empty()) {
        if (!progress.str().empty()) {
            progress << " ";
        }
        progress << workVerb << " " << workUnitName;
        verbed = true;
    }
    if (verbed) {
        if (workDone < totalWorkToDo) {
            if (!progress.str().empty()) {
                progress << " ";
            }
            progress << workDone << " (of " << totalWorkToDo << ")";
        }
    } else if (0<totalWorkToDo) {
        double percentDone = 100.0 * ( workDone / totalWorkToDo );
        progress.precision(3);
        progress << " " << percentDone << "% done";
    }
    progress.precision(6);
    progress << " in " << elapsedTime << " secs";
    if (0<elapsedTime && lastReportedCPUTime < cpu) {
        progress.precision(4);
        double percentCPU = 100.0 * ( (cpu-startCPUTime) / elapsedTime);
        progress << " using " << percentCPU << "% CPU";
    }
    double estimatedTime = 0.0; //Estimated work still to do in seconds
    if (0.0 < workDone && 0.0 < elapsedTime && workDone < totalWorkToDo ) {
        estimatedTime = ((totalWorkToDo - workDone) / workDone) * elapsedTime;
        if (estimatedTime < 10.0) {
            progress.precision(5);
        } else if (estimatedTime < 100.0) {
            progress.precision(4);
        }
        if (estimatedTime < 600.0 ) {
            progress << " (" << estimatedTime << " secs to go)";
        } else if (estimatedTime < 7200.0 ) {
            progress << " (" << (size_t)(floor(estimatedTime/60.0)) << " min to go)";
        } else {
            progress << " (" << (size_t)(floor(estimatedTime/3600.0)) << " hrs, "
                << (((size_t)(floor(estimatedTime/60.0)))%60) << " min to go)";
        }
    }
    std::string message = progress.str();
    #if _OPENMP
    #pragma omp critical (io)
    #endif
    {
        lastReportedWork = workDone;
        lastReportedTime = time;
        lastReportedCPUTime = cpu;
        if (isTerminal) {
            std::cout << "\33[2K\r";
        }
        if (displayingProgress || newline) {
            int barLen = 80;
            if (!newline && workDone < totalWorkToDo) {
                if (message.length() < barLen ) {
                    message += std::string(barLen-message.length(), ' ');
                }
                size_t charsInGreen = (size_t) floor( workDone * barLen / totalWorkToDo );
                if (isTerminal && charsInGreen < message.length()) {
                    size_t charsInGreenOrCyan //number of chars in green or blue
                        = (( message.length() < barLen) ? message.length() : barLen);
                    message = "\33[1;30;102m" + message.substr(0, charsInGreen)
                    + "\33[1;30;106m" + message.substr(charsInGreen, charsInGreenOrCyan - charsInGreen)
                    + "\33[0m" + message.substr(charsInGreenOrCyan, message.length() - charsInGreenOrCyan);
                }
            }
            std::cout << message;
            if (newline) {
                std::cout << std::endl;
            } else {
                std::cout.flush();
            }
        }
    }
}

progress_display& progress_display::done() {
    #if _OPENMP
    #pragma omp critical (io)
    #endif
    workDone = totalWorkToDo;
    isDone = true;
    reportProgress(getRealTime(), getCPUTime(), true);
    return *this;
}


progress_display::~progress_display() {
    if (!isDone) {
        workDone = totalWorkToDo;
        reportProgress(getRealTime(), getCPUTime(), true);
    }
}

progress_display& progress_display::hide() {
    if (!isTerminal) {
        return *this;
    }
    #if _OPENMP
    #pragma omp critical (io)
    #endif
    {
        std::cout << "\33[2K\r";
        std::cout.flush();
    }
    return *this;
}

progress_display& progress_display::show() {
    if (!isTerminal) {
        return *this;
    }
    reportProgress(getRealTime(), getCPUTime(), false);
    return *this;
}

void progress_display::setProgressDisplay(bool displayIt) {
    displayingProgress = displayIt;
    isTerminal = displayIt && isStandardOutputATerminal();
}

bool progress_display::getProgressDisplay() {
    return displayingProgress;
}
