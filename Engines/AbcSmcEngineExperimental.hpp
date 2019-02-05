//
// Created by Nikita Kruk on 08.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCSMCENGINEEXPERIMENTAL_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCSMCENGINEEXPERIMENTAL_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/ExperimentalParameterSet.hpp"
#include "../Other/FreeBoundaryConditionsConfiguration.hpp"

#include <fstream>
#include <random>
#include <map>
#include <eigen3/Eigen/Dense>

class AbcSmcEngineExperimental
{
 public:

  AbcSmcEngineExperimental();
  ~AbcSmcEngineExperimental();

  void RunAbcTowardsExperimentalData();

 private:

  ExperimentalParameterSet experimental_parameter_set_;
  FreeBoundaryConditionsConfiguration fbc_config_;

  // containers for data comparison
  std::map<int, Eigen::VectorXd> experimental_data_;
  std::vector<Real> candidate_data_;
  std::vector<int> initial_indexes_;

  void ComputeDistanceBetweenExperimentalAndCandidateDatasets(const std::string &experimental_data_file_name,
                                                              const std::string &candidate_dataset_file_name,
                                                              Real &distance);
  void RetrieveNumberOfAgentsAndInitialIndexes(const std::string &file_name);
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

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_ABCSMCENGINEEXPERIMENTAL_HPP
