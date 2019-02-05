//
// Created by Nikita Kruk on 25.07.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_THREAD_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_THREAD_HPP

#include "../Definitions.hpp"

#include <vector>

class Thread
{
 public:

  Thread(int argc, char **argv);
  ~Thread();

  bool IsRoot();
  int GetNumberOfMpichThreads();
  int GetRank();

 private:

  int root_rank_;
  int rank_;
  int number_of_mpich_threads_;

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_THREAD_HPP
