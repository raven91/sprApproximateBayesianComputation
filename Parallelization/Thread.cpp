//
// Created by Nikita Kruk on 25.07.18.
//

#include "Thread.hpp"

#include <cassert>

#if defined(MPI_PARALLELIZATION)
#include <mpi.h>
#endif

#if defined(MPI_PARALLELIZATION)
Thread::Thread(int argc, char **argv)
{
  root_rank_ = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_mpich_threads_);
}
#else
Thread::Thread(int argc, char **argv)
{
  root_rank_ = 0;
  rank_ = 0;
  number_of_mpich_threads_ = 1;
}
#endif

#if defined(MPI_PARALLELIZATION)
Thread::~Thread()
{

}
#else
Thread::~Thread()
{

}
#endif

#if defined(MPI_PARALLELIZATION)
bool Thread::IsRoot()
{
  return (rank_ == root_rank_);
}
#else
bool Thread::IsRoot()
{
  return true;
}
#endif

#if defined(MPI_PARALLELIZATION)
int Thread::GetNumberOfMpichThreads()
{
  return number_of_mpich_threads_;
}
#else
int Thread::GetNumberOfMpichThreads()
{
  return number_of_mpich_threads_;
}
#endif

#if defined(MPI_PARALLELIZATION)
int Thread::GetRank()
{
  return rank_;
}
#else
int Thread::GetRank()
{
  return 0;
}
#endif