//
// Created by Nikita Kruk on 03.12.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCSMCENGINE_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCSMCENGINE_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/SyntheticParameterSet.hpp"
#include "../Other/PeriodicBoundaryConditionsConfiguration.hpp"

#include <fstream>
#include <random>

class AbcSmcEngineSynthetic
{
 public:

  AbcSmcEngineSynthetic();
  ~AbcSmcEngineSynthetic();

  void PrepareSyntheticData();
  void RunAbcTowardsSyntheticData();

 private:

  SyntheticParameterSet synthetic_parameter_set_;
  PeriodicBoundaryConditionsConfiguration pbc_config_;

  // containers for data comparison
  std::vector<Real> synthetic_data_;
  std::vector<Real> candidate_data_;

  void ComputeDistanceBetweenSyntheticAndCandidateDatasets(const std::string &synthetic_data_file_name,
                                                           const std::string &candidate_dataset_file_name,
                                                           Real &distance);
  void SaveAcceptedParameters(int population,
                              const std::vector<Real> &sampled_parameters,
                              const std::vector<Real> &weights);
  void SampleParameterFromPreviousPopulation(int population_indicator,
                                             std::uniform_real_distribution<Real> &unif_rnd,
                                             const std::vector<Real> &weights,
                                             const std::vector<Real> &sampled_parameters,
                                             Real &new_parameter);
  void NormalizeWeights(std::vector<Real> &weights);

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCSMCENGINE_HPP
