//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_DEFINITIONS_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_DEFINITIONS_HPP

//#define MPI_PARALLELIZATION
//#define BCS_CLUSTER

typedef double Real;

const int kS = 4; // number of state variables <-> (x,y,v_x,v_y)
const int kSe = 8; // number of extended state variables <-> (x,y,v_x,v_y,area,slope,width,height)
const Real kMicroMetersPerPixel = 0.0639;
const Real kSecondsPerImage = 1e-3;
const Real kMetersPerMicroMeter = 1e-6;

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_DEFINITIONS_HPP




