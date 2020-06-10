
#include <sched.h>
#include <unistd.h>
#include <limits.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#include <sys/syscall.h>

#include <cerrno>
#include <string> 
#include <iostream>
#include <sstream>
#include <iomanip>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "affinity.hpp"

using namespace girih;

Affinity::Affinity() {

    // get processor number
    nb_cores_ = get_nprocs();

    // get setsize
    setsize_ = CPU_ALLOC_SIZE(nb_cores_);
     
    // allocate a mask set
    set_ = CPU_ALLOC(nb_cores_);
    if (set_ == NULL) throw std::system_error(errno, std::generic_category());

    // clean up
    CPU_ZERO_S(setsize_,set_);

}

Affinity::~Affinity() {
    // clean up mask set
    CPU_FREE(set_);
}

void Affinity::set_openmp(const int nb_threads,
                          const int thread0,
                          const int thread_stride) {
    int tid_omp,rc;
    pid_t pid,tid;

    if (thread0 < 0 || thread0 > nb_cores_) {
        throw std::invalid_argument("thread min is not valid");
    }

    if (thread0 + (nb_threads - 1) * thread_stride < 0 || 
        thread0 + (nb_threads - 1) * thread_stride > nb_cores_) {
        throw std::invalid_argument("thread max is not valid");
    }

#ifdef _OPENMP
    omp_set_max_active_levels(1);
    omp_set_num_threads(nb_threads);

    pid = getpid();
    
    CPU_ZERO_S(setsize_,set_);
    CPU_SET_S(thread0,setsize_,set_);
                
    rc = sched_setaffinity(pid,setsize_,set_);
    if (rc != 0) throw std::system_error(errno, std::generic_category());

    #pragma omp parallel default(shared) private(tid_omp,tid)
    {
        tid_omp = omp_get_thread_num();
        tid     = syscall(SYS_gettid);

        #pragma omp for ordered schedule(static,1)
        for (int i = 0; i < omp_get_num_threads(); i++) {
            #pragma omp ordered
            if (i == tid_omp) {
                CPU_ZERO_S(setsize_,set_);
                CPU_SET_S(thread0+i*thread_stride,setsize_,set_);
                
                rc = sched_setaffinity(tid,setsize_,set_);
                if (rc != 0) throw std::system_error(errno, std::generic_category());
            }
        }
    }
#endif

}
/*

///
/// \brief set nested openmp (currently disactivated since threads are not reused)
///
void Affinity::set_openmp_nested(const int nb_thread_groups,
                                 const int thread_group_size,
                                 const int thread0,
                                 const int thread_stride) {
    int rc,nb_threads;
    pid_t pid;

    if (thread0 < 0 || thread0 > nb_cores_) {
        throw std::invalid_argument("thread min is not valid");
    }

    nb_threads = nb_thread_groups * thread_group_size;

    if (thread0 + (nb_threads - 1) * thread_stride < 0 || 
        thread0 + (nb_threads - 1) * thread_stride > nb_cores_) {
        throw std::invalid_argument("thread max is not valid");
    }

#ifdef _OPENMP

    omp_set_max_active_levels(2);
    pid = getpid();
    
    CPU_ZERO_S(setsize_,set_);
    CPU_SET_S(thread0,setsize_,set_);
                
    rc = sched_setaffinity(pid,setsize_,set_);
    if (rc != 0) throw std::system_error(errno, std::generic_category());

    omp_set_num_threads(nb_thread_groups);

    #pragma omp parallel default(shared) //num_threads(nb_thread_groups)
    {
        int groupid = omp_get_thread_num();
        int tid     = syscall(SYS_gettid);

        #pragma omp for ordered schedule(static,1)
        for (int i = 0; i < omp_get_num_threads(); i++) {
            #pragma omp ordered
            if (i == groupid) {
                CPU_ZERO_S(setsize_,set_);
                CPU_SET_S(thread0+i*thread_group_size*thread_stride,setsize_,set_);

                rc = sched_setaffinity(tid,setsize_,set_);
                if (rc != 0) throw std::system_error(errno, std::generic_category());

                omp_set_num_threads(thread_group_size);

                #pragma omp parallel default(shared) //num_threads(thread_group_size)
                {   
                    int threadid = omp_get_thread_num();
                    int ttid     = syscall(SYS_gettid);

                    #pragma omp for ordered schedule(static,1)
                    for (int j = 0; j < omp_get_num_threads(); j++) {
                        #pragma omp ordered
                        if (j == threadid) {
                            CPU_ZERO_S(setsize_,set_);
                            CPU_SET_S(thread0+i*thread_group_size*thread_stride + j,setsize_,set_);
                            std::cout << ttid<<" "<<groupid << threadid << omp_get_num_threads() << i*thread_group_size*thread_stride + j << std::endl;
                            rc = sched_setaffinity(ttid,setsize_,set_);
                            if (rc != 0) throw std::system_error(errno, std::generic_category());                            
                        }
                    }
                }
            }
        }
    }
#endif
}
*/

///
/// \brief convert a cpuset to string
///
void Affinity::cpuset_get_string(std::string& str) const {
    int start, stop;
    std::ostringstream oss;

    start = -1;
    stop  = -1;

    for (int i = 0; i < nb_cores_; i++) {
        if (CPU_ISSET_S(i,setsize_,set_)) {
            stop = i;
            if (start == -1) start = i;
        } else {
            if (start != -1) {
                if (oss.tellp() != 0) oss << ",";

                oss << start;
                if (start != stop) oss << " - " << stop;

                start = -1;
                stop  = -1;
            }
        }
    }

    if (start != -1) {
        if (oss.tellp() != 0) oss << ",";

        oss << start;
        if (start != stop) oss << " - " << stop;
                
        start = -1;
        stop  = -1;
    }

    str = oss.str();

}

///
/// \brief set cpuset with pid
///
void Affinity::pid_get_affinity(const pid_t tid,
                                std::string& str,
                                int& nb_cpu_allowed) {
    int rc;

    CPU_ZERO_S(setsize_,set_);

    // get cpuset of tid
    rc = sched_getaffinity(tid,setsize_,set_);
    if (rc != 0) throw std::system_error(errno, std::generic_category());

    // set allowed cpu number
    nb_cpu_allowed = CPU_COUNT_S(setsize_,set_);

    // convert cpuset to string
    cpuset_get_string(str);
}

///
/// \brief get a string to describe the affinity configuration
///
void Affinity::get_string(std::string& str,
                          int rank) {
    int rc,nb_cpu_allowed;

    pid_t pid;
    std::string affinity;
    char hostname[HOST_NAME_MAX + 1];
    std::ostringstream oss;

    // get hostname
    rc = gethostname(hostname,HOST_NAME_MAX);
    if (rc != 0) throw std::system_error(errno, std::generic_category());

    // processus affinity
    pid = getpid();
    pid_get_affinity(pid,affinity,nb_cpu_allowed);

    oss << std::thread::hardware_concurrency() 
        << " concurrent threads are supported." << std::endl;

    oss << "rank " << std::setw(3) << rank << " "
        << "maps to " << std::setw(3) << nb_cpu_allowed << " core "
        << "[" << affinity << "] on " << hostname 
        << "(" << pid << ")" << std::endl;
    

#ifdef _OPENMP
    // thread affinity
    if (omp_get_max_active_levels() ==1 ) {
        #pragma omp parallel default(shared) 
        {
            int groupid = omp_get_thread_num();
            int tid     = syscall(SYS_gettid);

            #pragma omp for ordered schedule(static,1) 
            for (int i = 0; i < omp_get_num_threads(); i++) {
                #pragma omp ordered
                if (i == groupid) {
                    pid_get_affinity(tid,affinity,nb_cpu_allowed);

                    oss << "          thread " << std::setw(3) << groupid << " "
                        << "maps to " << std::setw(3) << nb_cpu_allowed << " core "
                        << "[" << affinity << "] on " << hostname 
                        << "(" << tid << ")" << std::endl;
                }
            }
        }    
    } else {
        #pragma omp parallel default(shared)
        {
            int groupid = omp_get_thread_num();
            int tid     = syscall(SYS_gettid);

            #pragma omp for ordered schedule(static,1)
            for (int i = 0; i < omp_get_num_threads(); i++) {
                #pragma omp ordered
                if (i == groupid) {
                    pid_get_affinity(tid,affinity,nb_cpu_allowed);

                    oss << "          thread group " << std::setw(3) << groupid << " "
                        << "maps to " << std::setw(3) << nb_cpu_allowed << " core "
                        << "[" << affinity << "] on " << hostname 
                        << "(" << tid << ")" << std::endl;
                   
                    #pragma omp parallel default(shared)
                    {
                        int threadid = omp_get_thread_num();
                        int tid     = syscall(SYS_gettid);

                        #pragma omp for ordered schedule(static,1)
                        for (int j = 0; j < omp_get_num_threads(); j++) {
                            #pragma omp ordered
                            if (j == threadid) {
                                pid_get_affinity(tid,affinity,nb_cpu_allowed);

                                oss << "          thread (" << std::setw(3) << groupid << "," << std::setw(3) << threadid << ") "
                                    << "maps to " << std::setw(3) << nb_cpu_allowed << " core "
                                    << "[" << affinity << "] on " << hostname 
                                    << "(" << tid << ")" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    std::cout << oss.str() << std::endl;
}