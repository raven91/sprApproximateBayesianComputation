//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONENGINE_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONENGINE_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/SyntheticParameterSet.hpp"
#include "../Other/PeriodicBoundaryConditionsConfiguration.hpp"
#include "../Parallelization/Thread.hpp"

#include <random>

class SimulationEngineForSyntheticData
{
 public:

  explicit SimulationEngineForSyntheticData(SyntheticParameterSet &parameter_set,
                                            PeriodicBoundaryConditionsConfiguration &pbc_config);
  ~SimulationEngineForSyntheticData();

  // Integrate the SPR model and save the dynamics
  void GenerateSyntheticDataWithUnitVelocity();
  void GenerateSyntheticDataForDifferentParameters(Thread *thread = nullptr);
  void GenerateSyntheticDataWithArbitraryStationaryVelocity();
  void SimulateCandidateDataset(int population_idx,
                                const std::string &synthetic_data_file_name,
                                std::string &simulation_data_file_name);
  const std::vector<Real> &GetVelocityDistribution();

 private:

  SyntheticParameterSet &parameter_set_;
  PeriodicBoundaryConditionsConfiguration &pbc_config_;
  std::vector<Real> system_state_;
  std::vector<Real> velocity_distribution_;

  void InitializeRandomSystemState();
  void InitializeSystemStateFromSyntheticData(const std::string &synthetic_data_file_name);
  void InitializeRandomSystemStateWithVelocityDistribution();
  void DrawVelocitiesFromPriorDistribution();
  Real DrawRayleighSample(std::uniform_real_distribution<Real> &uniform_distribution,
                          std::mt19937 &mersenne_twister_generator,
                          Real b);
};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SIMULATIONENGINE_HPP
