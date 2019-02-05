//
// Created by Nikita Kruk on 25.07.18.
//

#include "Parallelization.hpp"

#if defined(MPI_PARALLELIZATION)
#include <mpi.h>
#endif

#if defined(MPI_PARALLELIZATION)
void LaunchParallelSession(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
}
#else
void LaunchParallelSession(int argc, char **argv)
{

}
#endif

#if defined(MPI_PARALLELIZATION)
void FinalizeParallelSession()
{
  MPI_Finalize();
}
#else
void FinalizeParallelSession()
{

}
#endif