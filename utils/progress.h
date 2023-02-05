//
//  progress.h
//
//  Created by James Barbetti on 28/7/20.
//

#ifndef progress_h
#define progress_h

#if USE_PROGRESS_DISPLAY

#include <string>       //for std::string
#include <fstream>      //for std::fstream
#include <math.h>       //for floor()
#include "timeutil.h"   //for getRealTime()

/**
 * @brief Reports progress, for a task, by displaying a progress bar on
 *        the terminal, that is periodically updated, as it is informed
 *        how much progress has been made.
 */
class progress_display {
public:
    explicit progress_display ( double workToDo, const char* doingWhat=""
                              , const char* verb="", const char* unitName="");
    explicit progress_display ( size_t workToDo, const char* doingWhat = ""
                              , const char* verb = "", const char* unitName = "");
    explicit progress_display ( intptr_t workToDo, const char* doingWhat = ""
                              , const char* verb = "", const char* unitName = "");
    ~progress_display();
    progress_display& operator += (double   incrementalWork); //report progress
    progress_display& operator += (size_t   incrementalWork); //report progress
    progress_display& operator += (intptr_t incrementalWork); //report progress
    progress_display& operator ++ ();                         //increment progress
    progress_display& operator =  (double   workDoneNow); //set total progress
    progress_display& incrementBy (double   increment);   //report progress
    progress_display& incrementBy (size_t   increment);   //report progress
    progress_display& incrementBy (intptr_t increment);   //report progress
    progress_display& hide (); //Undisplay progress bar
    progress_display& show (); //Display progress bar
    progress_display& done (); //Report that the task, whose progress was
                               //being reported, is completed.
    void setTaskDescription  (const  char* newDescription);
    void setTaskDescription  (const  std::string& newDescription);
    void setWorkRemaining    (double newEstimate);
    void setIsEstimateABound (bool   isEstimateAnUpperBound);
    void reportProgress      (double time, double cpu, bool newline);
    void markAsFailed        (); //Report that the task, whose progress is/was
                                 //being reported, has failed to complete.
                                 //(e.g. writing file, file I/O failed)
    static void setProgressDisplay(bool displayIt);
    static bool getProgressDisplay();
protected:
    double       startTime;           //The time at which the task began
    double       startCPUTime;        //The CPU time that had accumulated
                                      //when the task began.
    double       totalWorkToDo;       //Amount of work to be done for this task
                                      //(e.g. how many bytes of file to read)
    double       workDone;            //Amount of work done so far
                                      //(e.g. how many bytes read thus far)
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

    /**
    * @brief Append a string, indicating how much progress has been made,
    *        to a stringstream (used to construct the text to report progress)
    * @param verbed   - reference to a boolean, into which true will be
    *                   written, if and only if the task has either finished
    *                   (all the work, or more than that, has been done)
    *                   or failed.
    * @param progress - reference to the stringstream to append
    */
    void appendHowMuchDone  (bool &verbed, std::ostringstream& progress);
    /**
     * @brief Append a string, indicating how much time and effort 
     *        (elapsed time and accumulated CPU) have been invested, in the
     *        task.
     * @note  It is assumed that ALL CPU being used, by any thread, is
     *        furthering the execution of the task.
     * @param verbed      - true if the task has been described as
     *                      "done" or "failed" already. false otherwise.
     * @param elapsedTime - how much time has been spent.
     * @param cpu         - how much CPU time has been spent.
     * @param progress    - reference to the stringstream to append.
     */
    void appendUsage        (bool verbed,  double elapsedTime, double cpu,
                             std::ostringstream& progress);
    /**
     * @brief Given a string, reformat it (using background and foreground
     *        colours, to look like a status bar).
     * @param message - a reference the message (on entry, as text)
     *                  on output, the message will be overwritten with
     *                  a string that will make the message appear to be
     *                  superimposed over a coloured progress bar.
     */
    void formatAsProgressBar(std::string&        message);
};

typedef progress_display* progress_display_ptr;

/**
 * @brief (Possibly) set up a pointer to a progress_display instance
 * @param wantProgress - true if progress is to be displayed
 * @param workToDo     - the total amount of work that is to be done by the task
 * @param description   - a description of the task "Reading file X"
 * @param verb          - the verb to prefix the "work done" unit with,
 *                        can be either past participle (e.g. read) or a 
 *                        gerund (e.g. reading). Example, in "read line 4",
 *                        "read" is the verb.
 * @param noun          - the noun for the work done unit (e.g. in 
 *                        "read line 4 of 1400", "line" is the noun).
 * @param progress      - a reference, to the "currently in scope task".
 *                        if we're creating a subtask, this'll be overwritten.
 * @param progress_here - a reference, to a pointer (to display_progress), that
 *                        will either be set to null (if no display_progress  
 *                        instance was created for a task or subtask),
 *                        or to a pointer, to the newly allocated instance
 *                        for the task or subtask created.
 */
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

/**
 * @brief hide progress
 * @param p a pointer to a progress_display instance
 */
inline void progressHide(progress_display_ptr p) { 
    if (p!=nullptr) {
        p->hide(); 
    }
}

/**
 * @brief show progress
 * @param p a pointer to a progress_display instance
 */
inline void progressShow(progress_display_ptr p) { 
    if (p!=nullptr) {
        p->show(); 
    }
}

/**
 * @brief record that the task is done
 * @param p a pointer to a progress_display instance
 */
inline void progressDone(progress_display_ptr p) { 
    if (p!=nullptr) {
        p->done(); 
    }
}

/**
 * @brief   delete (and stop referencing) a progress_display instance.
 * @param p reference to pointer to a progress_display instance.
 *          the pointed-to instance (if there is one) will be deleted,
 *          and the pointer (that is referenced) will be set to null.
 */
inline void progressDelete(progress_display_ptr& p) {
    delete p;
    p = nullptr;
}

/**
 * @brief Called when a subtask is compelted
 * @param progress      - reference to a display instance pointer
 * @param progress_here - reference to the currently in-scope 
 *                        progress_display instance (if any) that
 *                        is to be marked as done, and then deleted.
 */
inline void progressLocalDone(progress_display_ptr& progress,
                              progress_display_ptr& progress_here) {
    if (progress_here!=nullptr) {
        progressDone(progress_here);
        progressDelete(progress_here);
        progress = nullptr;
    }
}

/**
 * @brief  Append a description of elapsed (or estimated) time to a stream
 * @tparam S - the type of the stream (usually std::string)
 * @param  elapsed_time  - the elapsed (or estimated!) time (in seconds)
 * @param  s - reference to the stream, to be appended
 * @note   In a string like "6 min, 4 sec spent (5 sec to go)"
 *         the "6 min, 4 sec" and "5 sec" substrings are descriptions of
 *         elapsed (the first) and estimated remaining (the second) times.
 */
template <class S>
void appendTimeDescription(double elapsed_time, S &s) {
    if (elapsed_time < 60.0) /* less than a minute */ {
        s.precision(4);
        s << elapsed_time << " secs";
    }
    else {
        int64_t seconds = static_cast<int>(floor(elapsed_time));
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
#else
//Dummied versions of progress_display and progress reporting functions
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
#endif /* progress_h */
