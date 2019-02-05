//
// Created by Nikita Kruk on 08.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONENGINEFOREXPERIMENTALDATA_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONENGINEFOREXPERIMENTALDATA_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/ExperimentalParameterSet.hpp"
#include "../Other/FreeBoundaryConditionsConfiguration.hpp"
#include "../Parallelization/Thread.hpp"

#include <random>

class SimulationEngineForExperimentalData
{
 public:

  explicit SimulationEngineForExperimentalData(ExperimentalParameterSet &parameter_set,
                                               FreeBoundaryConditionsConfiguration &fbc_config);
  ~SimulationEngineForExperimentalData();

  // Integrate the SPR model and save the dynamics
  void SimulateCandidateDataset(int population_idx,
                                const std::string &experimental_data_file_name,
                                std::string &simulation_data_file_name,
                                const std::vector<int> &initial_indexes);

 private:

  ExperimentalParameterSet &parameter_set_;
  FreeBoundaryConditionsConfiguration &fbc_config_;
  std::vector<Real> system_state_;
  std::vector<Real> velocity_distribution_;

  void ReconstructVelocityDistribution(const std::string &experimental_data_file_name,
                                         const std::vector<int> &initial_indexes);
  void InitializeSystemStateFromExperimentalData(const std::string &experimental_data_file_name,
                                                   const std::vector<int> &initial_indexes);

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONENGINEFOREXPERIMENTALDATA_HPP
