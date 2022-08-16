/*----------------------------------------------------------------------------
VARS.H

Note: #define VAR_DECLS 1 before including this file to DECLARE and INITIALIZE
global variables.  Include this file without defining VAR_DECLS to extern
these variables.
----------------------------------------------------------------------------*/
#ifndef VAR_DEFS          // Make sure this file is included only once
#define VAR_DEFS 1

#include <memory>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>

/*----------------------------------------------
Setup variable declaration macros.
----------------------------------------------*/
#ifndef VAR_DECLS
# define _DECL extern
# define _INIT(x)
#else
# define _DECL
# define _INIT(x)  = x
#endif

/*----------------------------------------------
Declare variables as follows:

_DECL [standard variable declaration] _INIT(x);

where x is the value with which to initialize
the variable.  If there is no initial value for
the variable, it may be declared as follows:

_DECL [standard variable declaration];
----------------------------------------------*/
class ParallelContext;
struct Task;

_DECL std::atomic<int> global_terrace_trees;
_DECL std::atomic<int> global_intermediate_trees;
_DECL std::atomic<int> global_dead_ends;

static std::shared_ptr<ParallelContext> pContext;
static std::queue<Task*> taskQueue;

template<typename T, typename... Args>
std::shared_ptr<T> make_shared_global(Args&&... args) {
    return std::shared_ptr<T>(new T(std::forward<Args>(args)...));
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique_global(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template< typename T >
struct delete_pointer_element
{
    void operator()( T element ) const
    {
        delete element;
    }
};

//std::once_flag flag_1, flag_2;

#endif
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
