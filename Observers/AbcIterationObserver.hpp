//
// Created by Nikita Kruk on 26.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONRUNOBSERVER_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONRUNOBSERVER_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/AbstractParameterSet.hpp"
#include "SyntheticObserver.hpp"

#include <string>
#include <fstream>
#include <chrono>

class AbcIterationObserver : public SyntheticObserver
{
 public:

  explicit AbcIterationObserver(AbstractParameterSet *parameter_set,
                                const std::vector<Real> &velocity_distribution,
                                int abc_iteration = 0);

  virtual void operator()(std::vector<Real> &system_state, Real t);
  std::string GetDataFileName();

 private:

  int abc_iteration_;

  virtual std::string ConstructDataFileName();
  virtual std::string ConstructParameterFileName();
  virtual std::string ConstructActivePassiveFileName();

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONRUNOBSERVER_HPP
