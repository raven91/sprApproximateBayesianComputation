//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCREJECTIONENGINE_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCREJECTIONENGINE_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/SyntheticParameterSet.hpp"
#include "../Other/PeriodicBoundaryConditionsConfiguration.hpp"

#include <fstream>
#include <random>

class AbcRejectionEngine
{
 public:

  AbcRejectionEngine();
  ~AbcRejectionEngine();

  void PrepareSyntheticData();
  void RunAbcTowardsSyntheticData();
  void RunAbcTowardsExperimentalData();

 private:

  SyntheticParameterSet synthetic_parameter_set_;
  PeriodicBoundaryConditionsConfiguration pbc_config_;

  // containers for data comparison
  std::vector<Real> synthetic_data_;
  std::vector<Real> candidate_data_;

  void SampleNewParameter(std::uniform_real_distribution<Real> &unif_rnd, Real &parameter);
  void ComputeDistanceBetweenSyntheticAndCandidateDatasets(const std::string &synthetic_data_file_name,
                                                           const std::string &candidate_dataset_file_name,
                                                           Real &distance);
  void SaveAcceptedParameters(const std::vector<Real> &sampled_parameters);

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCREJECTIONENGINE_HPP
