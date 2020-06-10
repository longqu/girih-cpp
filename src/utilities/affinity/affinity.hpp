#ifndef GIRIH_AFFINITY_
#define GIRIH_AFFINITY_

#include <sched.h>
#include <exception>

namespace girih {

class Affinity : public std::exception {
public: 
    Affinity();
    ~Affinity();

    void set_openmp(const int nb_threads,
                    const int thread0 = 0,
                    const int thread_stride = 1);

    void set_openmp_nested(const int nb_thread_groups,
                           const int thread_group_size,
                           const int thread0 = 0,
                           const int thread_stride = 1);

    void get_string(std::string& str,
                    int rank);

private:
  int nb_cores_;
  size_t setsize_;
  cpu_set_t* set_;

  void cpuset_get_string(std::string& str) const;
  void pid_get_affinity(const pid_t tid,
                        std::string& str,
                        int& nb_cpu_allowed);
};

} // namespace girih

#endif