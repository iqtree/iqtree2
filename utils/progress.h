//
//  progress.hpp
//  alignment
//
//  Created by James Barbetti on 28/7/20.
//

#ifndef progress_h
#define progress_h

#include <string>
#include "timeutil.h"

class progress_display {
public:
    double startTime;
    double startCPUTime;
    double totalWorkToDo;
    double workDone;             //Amount of work done so far
    std::string taskDescription; //A description of the task
    bool        isDone;          //True if the task has been completed
    std::string workVerb;        //A past tense verb (e.g. sorted)
    std::string workUnitName;    //A singular noun   (e.g. row)
    double lastReportedWork;     //Last reported work done
    double lastReportedTime;     //Last reported wall-clock time
    double lastReportedCPUTime;  //Last reported CPU time
    
    explicit progress_display ( double workToDo, const char* doingWhat=""
                              , const char* verb="", const char* unitName="");
    ~progress_display();
    progress_display& operator += (double incrementalWork);
    progress_display& operator ++ ();
    progress_display& operator =  (double workDoneNow);
    progress_display& hide ();
    progress_display& show ();
    progress_display& done ();
    void reportProgress(double time, double cpu, bool newline);
    static void setProgressDisplay(bool displayIt);
    static bool getProgressDisplay();
};


#endif /* progress_hpp */
