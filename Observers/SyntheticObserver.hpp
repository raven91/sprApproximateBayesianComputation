//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SYNTHETICOBSERVER_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SYNTHETICOBSERVER_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/AbstractParameterSet.hpp"

#include <string>
#include <fstream>
#include <chrono>
#include <vector>

class SyntheticObserver
{
 public:

  explicit SyntheticObserver(AbstractParameterSet *parameter_set,
                             const std::vector<Real> &velocity_distribution);
  explicit SyntheticObserver(AbstractParameterSet *parameter_set,
                             const std::vector<Real> &velocity_distribution,
                             int no_file_creation);
  ~SyntheticObserver();

  virtual void operator()(std::vector<Real> &system_state, Real t);
  void SaveAllParametersIntoFile(std::string &parameters_file_name);

 protected:

  std::string data_file_name_;
  std::ofstream data_file_;
  AbstractParameterSet *parameter_set_;
  std::chrono::time_point<std::chrono::system_clock> integration_step_timer_;
  int ouput_time_counter_;
  int output_time_threshold_;
  std::vector<Real> velocity_distribution_;

  virtual std::string ConstructDataFileName();
  virtual std::string ConstructParameterFileName();
  virtual std::string ConstructActivePassiveFileName();

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SYNTHETICOBSERVER_HPP
