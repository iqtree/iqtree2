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

#ifdef  USE_PROGRESS_DISPLAY
#undef  USE_PROGRESS_DISPLAY
#endif
#define  USE_PROGRESS_DISPLAY (1)

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

progress_display::progress_display( double      workToDo, const char* doingWhat
                                  , const char* verb,     const char* unitName)
    : startTime(getRealTime()),   startCPUTime(getCPUTime())
    , totalWorkToDo(workToDo),    workDone(0.0)
    , taskDescription(doingWhat), isDone(false)
    , workVerb(verb),             workUnitName(unitName)
    , atMost(false),              hidden(0)
    , termout(CONSOLE_FILE, std::ios_base::out)
    , failed(false) {
        lastReportedWork    = 0.0;
        lastReportedTime    = startTime;
        lastReportedCPUTime = startCPUTime;
}

progress_display::progress_display(size_t workToDo, const char* doingWhat, 
                                   const char* verb, const char* unitName)
    : progress_display((double)workToDo, doingWhat, verb, unitName) {}

progress_display::progress_display(intptr_t workToDo, const char* doingWhat,
    const char* verb, const char* unitName)
    : progress_display((double)workToDo, doingWhat, verb, unitName) {}

progress_display & progress_display::operator ++ () {
    return (*this) += 1.0;
}

progress_display& progress_display::incrementBy (double increment) {
    return (*this) += increment;
}

progress_display& progress_display::incrementBy(size_t increment) {
    return (*this) += increment;
}

progress_display& progress_display::incrementBy(intptr_t increment) {
    return (*this) += increment;
}

/**
 * @brief  Set the amount of progress made, carrying out the task
 * @param  workDoneNow - what the amount of progress made should be, now
 * @return progress_display& - reference to *this
 * @note   Any attempt to retrospectively lower the amount of progress made
 *         will be ignored.
 */
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

/**
 * @brief  Increment the amount of progress made on a task
 * @param  incrementalWork  - how much additional progress has been made
 * @return progress_display& - reference to *this
 * @note   If the task has no name, the progress bar won't be updated.
 * @note   The progress bar is only updated at most once per second.
 * @note   If no progress has been made (if incrementalWork is zero),
 *         the progress bar will not be updated.
 */
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

progress_display& progress_display::operator += (size_t incrementalWork) {
    return *this += (double)incrementalWork;
}

progress_display& progress_display::operator += (intptr_t incrementalWork) {
    return *this += (double)incrementalWork;
}

/**
 * @brief Update progress bar and report progress
 * @param time    - the current time (as returned by getRealTime())
 * @param cpu     - the current accumulated count of CPU time used, by this
 *                  process, as returned by getCPUTime().
 * @param newline - true if a new line is to be appended to the report of
 *                  progress (typically, when the task has completed, or
 *                  has failed).)
 */
void progress_display::reportProgress
    (double time, double cpu, bool newline) {
    double elapsedTime = time - startTime;
    std::ostringstream progress;
    bool verbed = false;

    appendHowMuchDone(verbed, progress);
    progress << " in ";
    appendTimeDescription(elapsedTime, progress);
    appendUsage(verbed, elapsedTime, cpu, progress);
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
            if (newline) {
                if (!termout.fail()) {
                    termout.flush();
                }
                std::cout << message << std::endl;
                #if defined(CLANG_UNDER_VS)
                    OutputDebugStringA((message + "\n").c_str());
                #endif
            } else {
                if (workDone < totalWorkToDo) {
                    formatAsProgressBar(message);
                }
                if (!termout.fail()) {
                    termout << message;
                    termout.flush();
                }
            }
        }
    }
}

/**
 * @brief Record that the task, for which progress is reported by this
 *        instance, has failed outright (e.g. couldn't write file, I/O
 *        error part way through writing it).
 */
void progress_display::markAsFailed() {
    bool already_failed = failed;
    failed = true;
    if (!already_failed) {
        double time = getRealTime();
        double cpu  = getCPUTime();
        reportProgress(time, cpu, false);
    }
}

/**
 * @brief append a description of how much progress has been made to a 
 *        stringstream (that is uspplied by reference).
 * @param verbed   - reference to a boolean that will be set to true if
 *                   the word "done" or the word "failed" has been added
 * @param progress - reference to a stringstream to append
 */
void progress_display::appendHowMuchDone(bool &verbed, std::ostringstream& progress) {
    if (!taskDescription.empty()) {
        progress << taskDescription << ":";
    }
    if (totalWorkToDo <= workDone || failed) {
        if (!taskDescription.empty()) {
            if (failed) {
                progress << " failed";
            }
            else {
                progress << " done";
            }
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
}

/**
 * @brief Append a description of how much CPU time has been spent
 *        (and an estimate of the percent CPU), and if the task isn't done
 *        or failed, and estimate of how much more time it will take to
 *        finish the task.
 * @param verbed      - true if the string under construction (in progress)
 *                      has been described as done (if so, don't want to say
 *                      how much longer it is likely to take), or failed (ditto).
 *                      otherwise, false.
 * @param elapsedTime - elapsed time.
 * @param cpu         - accumulated CPU time.
 * @param progress    - reference to the stringstream to be appended.
 */
void progress_display::appendUsage(bool verbed, double elapsedTime, double cpu, 
                                   std::ostringstream& progress) {
    if (0<elapsedTime && lastReportedCPUTime < cpu) {
        progress.precision(4);
        double percentCPU = 100.0 * ( (cpu-startCPUTime) / elapsedTime);
        progress << " using " << percentCPU << "% CPU";
    }
    if (0.0 < workDone && 0.0 < elapsedTime && workDone < totalWorkToDo ) {
        //Estimated work still to do in seconds
        double estimatedTime = ((totalWorkToDo - workDone) / workDone) * elapsedTime;
        const char* leadIn = ( atMost && !verbed) ? " (at most " : " (";
        progress << leadIn;
        appendTimeDescription(estimatedTime, progress);
        progress << " to go)";
    }
}

/**
 * @brief Given a reference to a message, reformat that message so that
 *        it looks like a progress bar.
 * @param message - reference to the message to be reformatted
 */
void progress_display::formatAsProgressBar(std::string& message) {
    size_t barLen = 80;
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

/**
 * @brief  Report that the task, whose progress is being reported, is done.
 * @return progress_display& - reference to *this
 */
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

/**
 * @brief  Hide a progress display progress bar
 * @return progress_display& - reference to *this
 * @note   works by going back to the start of the line 
 *         and clearing what follows.
 * @note   progress_display tracks how many "times" it has been hidden, so
 *         if it is has been hidden n times it must later be shown n times
 *         before a progress bar will be displayed again.
 */
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

/**
 * @brief  Show a progress display progress bar
 * @return progress_display& - reference to *this
 * @note   Basically decrements a counter that tracks "how hidden" the progress
 *         bar is, and if that goes to zero, redisplays a progress bar line.
 */
progress_display& progress_display::show() {
    if (!isTerminal) {
        return *this;
    }
    if ( --hidden == 0 && !taskDescription.empty() ) {
        reportProgress(getRealTime(), getCPUTime(), false);
    }
    return *this;
}

/**
 * @brief Change the task description for a task that is in progress
 * @param newDescription = The new description (if it differs from the
 *        the previous one, and the progress_display isn't "hidden",
 *        the progress bar will be recalculated and redisplayed with 
 *        the new description).
 */
void progress_display::setTaskDescription(const  char* newDescription) {
    if (this->taskDescription == newDescription) {
        return;
    }
    if (isTerminal && hidden <= 0 && !taskDescription.empty() ) {
        reportProgress(getRealTime(), getCPUTime(), false);
    }
}

/**
 * @brief Change the task description for a task that is in progress
 * @param newDescription - the new description 
 */
void progress_display::setTaskDescription(const  std::string& newDescription) {
    setTaskDescription(newDescription.c_str());
}

/**
 * @brief Change the estimate of how much progress is as yet unmaed
 *        (or, if you prefer, the estimate of how much work remains to be done)
 * @param newEstimate - the new estimate (in units of work)
 * @note  if newEstimate is less than zero it is ignored.
 */
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

/**
 * @brief Indicates whether an estimate of what is left to do is
 *        an upper bound on how much is left to do (how much progress
 *        is to be made, or how long it is expected to take)
 * @param isEstimateAnUpperBound - true if it's a lower bound
 *        (if it is, progress bars will say " at most" ... whatever);
 *        false if it is not (if it's not, progress bars won't say that).
 */
void progress_display::setIsEstimateABound(bool isEstimateAnUpperBound) {
    if (atMost == isEstimateAnUpperBound) {
        return;
    }
    atMost = isEstimateAnUpperBound;
    if (hidden <= 0 && !taskDescription.empty() ) {
        reportProgress(getRealTime(), getCPUTime(), false);
    }
}

/**
 * @brief Set whether progress should be displayed (globally, for all tasks)
 * @param displayIt - true if it should, false if it should not
 */
void progress_display::setProgressDisplay(bool displayIt) {
    displayingProgress = displayIt;
    isTerminal = displayIt && isStandardOutputATerminal();
}

/**
 * @brief  Indicate whether progress bars are being displayed.
 * @return true  - if progress bars are being displayed
 * @return false - if they are not (e.g. program running with a -q for quiet
 *                 command-line parameter, or something like that)
 */
bool progress_display::getProgressDisplay() {
    return displayingProgress;
}
