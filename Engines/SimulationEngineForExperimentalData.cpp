//
// Created by Nikita Kruk on 08.08.18.
//

#include "SimulationEngineForExperimentalData.hpp"
#include "../DynamicalSystems/SprSystemWithFreeBoundaries.hpp"
#include "../Steppers/StochasticEulerStepperWithArbitraryStationaryVelocity.hpp"
#include "../Observers/AbcIterationObserver.hpp"

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm> // std::for_each
#include <map>

#include <boost/numeric/odeint.hpp>
#include <eigen3/Eigen/Dense>

SimulationEngineForExperimentalData::SimulationEngineForExperimentalData(ExperimentalParameterSet &parameter_set,
                                                                         FreeBoundaryConditionsConfiguration &fbc_config)
    :
    parameter_set_(parameter_set),
    fbc_config_(fbc_config),
    system_state_(kS * parameter_set_.Get_N(), 0.0)
{

}

SimulationEngineForExperimentalData::~SimulationEngineForExperimentalData()
{

}

void SimulationEngineForExperimentalData::ReconstructVelocityDistribution(const std::string &experimental_data_file_name,
                                                                          const std::vector<int> &initial_indexes)
{
  velocity_distribution_ = std::vector<Real>(parameter_set_.Get_N(), 0.0);

  std::ifstream experimental_data_file(experimental_data_file_name, std::ios::in);
  assert(experimental_data_file.is_open());
  std::string line;
  for (int i = 0; i < parameter_set_.Get_first_image(); ++i)
  {
    std::getline(experimental_data_file, line);
  }

  int time_idx = 0;
  int number_of_agents = 0;
  int agent_idx = 0;
  Eigen::VectorXd r_ext(kSe);
  std::map<int, Real> accumulated_velocities;
  std::map<int, int> accumulated_times;
  while (experimental_data_file >> time_idx >> number_of_agents)
  {
    for (int alpha = 0; alpha < number_of_agents; ++alpha)
    {
      experimental_data_file >> agent_idx >> r_ext(0) >> r_ext(1) >> r_ext(2) >> r_ext(3)
                             >> r_ext(4) >> r_ext(5) >> r_ext(6) >> r_ext(7);
      if (accumulated_times.find(agent_idx) != accumulated_times.end())
      {
        accumulated_velocities[agent_idx] += r_ext.segment(2, 2).norm();
        ++accumulated_times[agent_idx];
      } else
      {
        accumulated_velocities[agent_idx] = r_ext.segment(2, 2).norm();
        accumulated_times[agent_idx] = 1;
      }
    }
  }

  for (int i_alpha = 0; i_alpha < initial_indexes.size(); ++i_alpha)
  {
    int alpha = initial_indexes[i_alpha];
    velocity_distribution_[i_alpha] =
        accumulated_velocities[alpha] / accumulated_times[alpha]
            * kMicroMetersPerPixel * kMetersPerMicroMeter / kSecondsPerImage;
  }

  std::ofstream output("/Users/nikita/Documents/spr/sprApproximateBayesianComputation/CandidateDatasets/velocity.bin",
                       std::ios::out | std::ios::trunc | std::ios::binary);
  output.write((char *)&velocity_distribution_[0], velocity_distribution_.size() * sizeof(Real));
}

void SimulationEngineForExperimentalData::InitializeSystemStateFromExperimentalData(const std::string &experimental_data_file_name,
                                                                                    const std::vector<int> &initial_indexes)
{
  ReconstructVelocityDistribution(experimental_data_file_name, initial_indexes);

  std::ifstream experimental_data_file(experimental_data_file_name, std::ios::in);
  assert(experimental_data_file.is_open());
  std::string line;
  for (int i = 0; i < parameter_set_.Get_first_image(); ++i)
  {
    std::getline(experimental_data_file, line);
  }
  int time_idx = 0;
  int number_of_agents = 0;
  int agent_idx = 0;
  Eigen::VectorXd r_ext(kSe);
  experimental_data_file >> time_idx >> number_of_agents;
  assert(number_of_agents == parameter_set_.Get_N());
  for (int alpha = 0; alpha < number_of_agents; ++alpha)
  {
    int vec_idx = kS * alpha;
    experimental_data_file >> agent_idx >> r_ext(0) >> r_ext(1) >> r_ext(2) >> r_ext(3)
                           >> r_ext(4) >> r_ext(5) >> r_ext(6) >> r_ext(7);
    system_state_[vec_idx] = r_ext(0) * kMicroMetersPerPixel * kMetersPerMicroMeter;  // m <- px
    system_state_[vec_idx + 1] = r_ext(1) * kMicroMetersPerPixel * kMetersPerMicroMeter;  // m <- px
    if (r_ext.segment(2, 2).squaredNorm() < 0.001) // if velocity is less than 0.001 px/image
    {
//      std::cout << r_ext.segment(2, 2).squaredNorm() << "|";
      r_ext(2) = std::cos(r_ext(5));
      r_ext(3) = std::sin(r_ext(5));
    } else
    {
      r_ext.segment(2, 2).normalize();
    }
    system_state_[vec_idx + 2] =
        r_ext(2) * kMicroMetersPerPixel * kMetersPerMicroMeter / kSecondsPerImage; // m/s <- px/image
    system_state_[vec_idx + 3] =
        r_ext(3) * kMicroMetersPerPixel * kMetersPerMicroMeter / kSecondsPerImage; // m/s <- px/image
  }
  experimental_data_file.close();
}

void SimulationEngineForExperimentalData::SimulateCandidateDataset(int population_idx,
                                                                   const std::string &experimental_data_file_name,
                                                                   std::string &simulation_data_file_name,
                                                                   const std::vector<int> &initial_indexes)
{
  InitializeSystemStateFromExperimentalData(experimental_data_file_name, initial_indexes);
  SprSystemWithFreeBoundaries system(&parameter_set_, &fbc_config_, velocity_distribution_);
  AbcIterationObserver observer(&parameter_set_, velocity_distribution_, population_idx);
  StochasticEulerStepperWithArbitraryStationaryVelocity stepper(&fbc_config_, &parameter_set_, velocity_distribution_);

  Real t = parameter_set_.Get_t_0();
  observer(system_state_, t);
  while (t <= parameter_set_.Get_t_1())
  {
    std::vector<Real> system_state_prev(system_state_);
    t += parameter_set_.Get_delta_t();
    stepper.do_step(&system, system_state_, t, parameter_set_.Get_delta_t());
    observer(system_state_, t);
  }
  simulation_data_file_name = observer.GetDataFileName();
}