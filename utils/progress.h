//
//  progress.h
//  alignment
//
//  Created by James Barbetti on 28/7/20.
//

#ifndef progress_h
#define progress_h

#if USE_PROGRESS_DISPLAY

#include <string>       //for std::string
#include <fstream>      //for std::fstream
#include <math.h>       //for floor()
#include "timeutil.h"

class progress_display {
public:
    explicit progress_display ( double workToDo, const char* doingWhat=""
                              , const char* verb="", const char* unitName="");
    explicit progress_display ( size_t workToDo, const char* doingWhat = ""
                              , const char* verb = "", const char* unitName = "");
    explicit progress_display ( intptr_t workToDo, const char* doingWhat = ""
                              , const char* verb = "", const char* unitName = "");
    ~progress_display();
    progress_display& operator += (double   incrementalWork);
    progress_display& operator += (size_t   incrementalWork);
    progress_display& operator += (intptr_t incrementalWork);
    progress_display& operator ++ ();
    progress_display& operator =  (double   workDoneNow);
    progress_display& incrementBy (double   increment);
    progress_display& incrementBy (size_t   increment);
    progress_display& incrementBy (intptr_t increment);
    progress_display& hide ();
    progress_display& show ();
    progress_display& done ();
    void setTaskDescription  (const  char* newDescription);
    void setTaskDescription  (const  std::string& newDescription);
    void setWorkRemaining    (double newEstimate);
    void setIsEstimateABound (bool   isEstimateAnUpperBound);
    void reportProgress      (double time, double cpu, bool newline);
    void markAsFailed        ();
    static void setProgressDisplay(bool displayIt);
    static bool getProgressDisplay();
protected:
    double       startTime;
    double       startCPUTime;
    double       totalWorkToDo;
    double       workDone;            //Amount of work done so far
    std::string  taskDescription;     //A description of the task
    bool         isDone;              //True if the task has been completed
    std::string  workVerb;            //A past tense verb (e.g. sorted)
    std::string  workUnitName;        //A singular noun   (e.g. row)
    double       lastReportedWork;    //Last reported work done
    double       lastReportedTime;    //Last reported wall-clock time
    double       lastReportedCPUTime; //Last reported CPU time
    bool         atMost;              //Estimate of work to do is an upper bound
    int          hidden;              //How many more times has hide() been called than show()
    std::fstream termout;             //Terminal output stream (not standard output!)
    bool         failed;              //Becomes true if markAsFailed() has been called.

    void appendHowMuchDone  (bool &verbed, std::ostringstream& progress);
    void appendUsage        (bool verbed,  double elapsedTime, double cpu,
                             std::ostringstream& progress);
    void formatAsProgressBar(std::string&        message);
};

typedef progress_display* progress_display_ptr;
inline void progressLocal(bool wantProgress, double workToDo,
                          const char* description, const char* verb,
                          const char* noun, progress_display_ptr& progress,
                          progress_display_ptr& progress_here ) {
    if (progress==nullptr && wantProgress) {
        #if USE_PROGRESS_DISPLAY
        progress_here = new progress_display( workToDo, description, verb, noun);
        progress = progress_here;
        #endif
    }
}
inline void progressHide(progress_display_ptr p) { 
    if (p!=nullptr) {
        p->hide(); 
    }
}
inline void progressShow(progress_display_ptr p) { 
    if (p!=nullptr) {
        p->show(); 
    }
}
inline void progressDone(progress_display_ptr p) { 
    if (p!=nullptr) {
        p->done(); 
    }
}
inline void progressDelete(progress_display_ptr& p) {
    delete p;
    p = nullptr;
}
inline void progressLocalDone(progress_display_ptr& progress,
                              progress_display_ptr& progress_here) {
    if (progress_here!=nullptr) {
        progressDone(progress_here);
        progressDelete(progress_here);
        progress = nullptr;
    }
}

#else
typedef double  progress_display;
typedef double* progress_display_ptr;
inline void progressLocal(bool wantProgress, double workToDo,
                          const char* description, const char* verb,
                          const char* noun, progress_display_ptr& progress,
                          progress_display_ptr& progress_here ) {}
inline void progressHide(progress_display_ptr) {}
inline void progressShow(progress_display_ptr) {}
inline void progressDone(progress_display_ptr) {}
inline void progressDelete(progress_display_ptr p) {}
inline void progressLocalDone(progress_display_ptr& progress,
                              progress_display_ptr& progress_here) {}
#endif

template <class S>
void appendTimeDescription(double elapsed_time, S &s) {
    if (elapsed_time < 60.0) /* less than a minute */ {
        s.precision(4);
        s << elapsed_time << " secs";
    }
    else {
        int64_t seconds = (int)floor(elapsed_time);
        int64_t minutes = seconds / 60;
        int64_t hours   = minutes / 60;
        int64_t days    = hours / 24;
        seconds -= minutes * 60;
        minutes -= hours * 60;
        hours   -= days * 24;
        if (28<days) {
            s << (days/7) << " weeks ";
            days = days % 7;
        }
        if (0<days) {
            s << days << " days ";
        }
        if (0<hours) {
            s << hours << " hrs ";
        }
        s << minutes << " min "
          << seconds << " sec";
    }
}


#endif /* progress_h */
