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
#include "operatingsystem.h" //for isStandardOutputATerminal and CONSOLE_FILE


namespace {
    bool displayingProgress = true;
        //You can turn off progress displays via progress_display::setProgressDisplay.
    bool isTerminal = false;
}

progress_display::progress_display( double workToDo, const char* doingWhat
                                   , const char* verb, const char* unitName)
    : startTime(getRealTime()),   startCPUTime(getCPUTime())
    , totalWorkToDo(workToDo),    workDone(0.0)
    , taskDescription(doingWhat), isDone(false)
    , workVerb(verb),             workUnitName(unitName)
    , atMost(false),              hidden(0)
    , termout(CONSOLE_FILE, std::ios_base::out)
     {
        lastReportedWork    = 0.0;
        lastReportedTime    = startTime;
        lastReportedCPUTime = startCPUTime;
}

progress_display & progress_display::operator ++ () {
    return (*this) += 1.0;
}

progress_display& progress_display::incrementBy (double increment) {
    return (*this) += increment;
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
    bool justASec = floor(time) > floor(lastReportedTime);
    if ( ( lastReportedWork == 0.0 || justASec ) && !taskDescription.empty() ) {
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
            progress << workDone << " (of " << (atMost ? "at most " : "") << totalWorkToDo << ")";
        }
    } else if (0<totalWorkToDo) {
        double percentDone = 100.0 * ( workDone / totalWorkToDo );
        progress.precision(3);
        progress << " " << percentDone << "% done";
    }
    if (elapsedTime < 60.0) /* less than a minute */ {
        progress.precision(6);
        progress << " in " << elapsedTime << " secs";
    } else {
        int64_t seconds = (int)floor(elapsedTime);
        int64_t minutes = seconds / 60;
        int64_t hours   = minutes / 60;
        int64_t days    = hours / 24;
        seconds -= minutes * 60;
        minutes -= hours * 60;
        hours   -= days * 24;
        progress << " in ";
        if (0<days) {
            progress << days << "day ";
        }
        if (0<hours) {
            progress << hours << " hr ";
        }
        progress << minutes << " min " << seconds << " sec";
    }
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
        const char* leadIn = ( atMost && !verbed) ? " (at most " : " (";
        if (estimatedTime < 600.0 ) {
            progress << leadIn << estimatedTime << " secs to go)";
        } else if (estimatedTime < 7200.0 ) {
            progress << leadIn << (size_t)(floor(estimatedTime/60.0)) << " min to go)";
        } else {
            progress << leadIn  << (size_t)(floor(estimatedTime/3600.0)) << " hrs, "
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
        if (isTerminal && !termout.fail()) {
            termout << "\33[2K\r";
        }
        if (displayingProgress) {
            int barLen = 80;
            if (newline) {
                if (!termout.fail()) {
                    termout.flush();
                }
                std::cout << message << std::endl;
            } else {
                if (workDone < totalWorkToDo) {
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
                if (!termout.fail()) {
                    termout << message;
                    termout.flush();
                }
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
    if (!taskDescription.empty()) {
        reportProgress(getRealTime(), getCPUTime(), true);
    }
    return *this;
}

progress_display::~progress_display() {
    if (!isDone) {
        done();
    }
    termout.close();
}

progress_display& progress_display::hide() {
    if (!isTerminal || termout.fail()) {
        return *this;
    }
    #if _OPENMP
    #pragma omp critical (io)
    #endif
    {
        if (++hidden == 1) {
            termout << "\33[2K\r";
            termout.flush();
        }
    }
    return *this;
}

progress_display& progress_display::show() {
    if (!isTerminal) {
        return *this;
    }
    if ( --hidden == 0 && !taskDescription.empty() ) {
        reportProgress(getRealTime(), getCPUTime(), false);
    }
    return *this;
}

void progress_display::setTaskDescription(const  char* newDescription) {
    if (this->taskDescription == newDescription) {
        return;
    }
    if (isTerminal && hidden <= 0 && !taskDescription.empty() ) {
        reportProgress(getRealTime(), getCPUTime(), false);
    }
}

void progress_display::setTaskDescription(const  std::string& newDescription) {
    setTaskDescription(newDescription.c_str());
}

void progress_display::setWorkRemaining(double newEstimate) {
    if (newEstimate < 0) {
        return; //Nonsense!
    }
    double oldWorkToDo, newWorkToDo;
    #if _OPENMP
    #pragma omp critical (io)
    #endif
    {
        oldWorkToDo = totalWorkToDo;
        totalWorkToDo = workDone + newEstimate;
        newWorkToDo = totalWorkToDo;
    }
    if (isTerminal && hidden <= 0 && 1.0 < abs(oldWorkToDo - newWorkToDo)) {
        reportProgress(getRealTime(), getCPUTime(), false);
    }
}

void progress_display::setIsEstimateABound(bool isEstimateAnUpperBound) {
    if (atMost == isEstimateAnUpperBound) {
        return;
    }
    atMost = isEstimateAnUpperBound;
    if (hidden <= 0 && !taskDescription.empty() ) {
        reportProgress(getRealTime(), getCPUTime(), false);
    }
}

void progress_display::setProgressDisplay(bool displayIt) {
    displayingProgress = displayIt;
    isTerminal = displayIt && isStandardOutputATerminal();
}

bool progress_display::getProgressDisplay() {
    return displayingProgress;
}
