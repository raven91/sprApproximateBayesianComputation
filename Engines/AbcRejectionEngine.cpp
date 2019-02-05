//
// Created by Nikita Kruk on 25.11.17.
//

#include "AbcRejectionEngine.hpp"
#include "SimulationEngineForSyntheticData.hpp"

#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm> // std::for_each, std::find
#include <cmath>

AbcRejectionEngine::AbcRejectionEngine() :
    synthetic_parameter_set_(),
    pbc_config_(synthetic_parameter_set_.Get_L(), synthetic_parameter_set_.Get_L()),
    synthetic_data_(std::vector<Real>(kS * synthetic_parameter_set_.Get_N(), 0.0)),
    candidate_data_(std::vector<Real>(kS * synthetic_parameter_set_.Get_N(), 0.0))
{
  std::cout << "Approximate Bayesian Computation - Rejection Started" << std::endl;
}

AbcRejectionEngine::~AbcRejectionEngine()
{
  std::cout << "Approximate Bayesian Computation - Rejection Complete" << std::endl;
}

void AbcRejectionEngine::PrepareSyntheticData()
{
  SimulationEngineForSyntheticData simulation_engine(synthetic_parameter_set_, pbc_config_);
  simulation_engine.GenerateSyntheticDataWithArbitraryStationaryVelocity();
}

void AbcRejectionEngine::RunAbcTowardsSyntheticData()
{
  std::ostringstream synthetic_data_file_name_buffer;
  synthetic_data_file_name_buffer << synthetic_parameter_set_.GetSyntheticDataFolderName()
                                  << "/spr_simulation_N_" << synthetic_parameter_set_.Get_N()
                                  << "_phi_" << synthetic_parameter_set_.Get_phi()
                                  << "_a_" << synthetic_parameter_set_.Get_a()
                                  << "_U0_" << synthetic_parameter_set_.Get_U_0()
                                  << "_k_" << synthetic_parameter_set_.Get_kappa() << ".bin";

  int population_size = 100;
  Real tolerance = 0.0000125;
  std::uniform_real_distribution<Real> unif_rnd(0.0, 1e-18);
  int particle_idx = 0;
  std::vector<Real> sampled_parameters;
  while (particle_idx < population_size)
  {
    Real U_0 = 0.0;
    SampleNewParameter(unif_rnd, U_0);
    synthetic_parameter_set_.Set_U_0(U_0);
    synthetic_parameter_set_.Set_t_0(0.0);
    synthetic_parameter_set_.Set_t_1(10.0);
    SimulationEngineForSyntheticData simulation_engine(synthetic_parameter_set_, pbc_config_);
    std::string candidate_dataset_file_name;
    simulation_engine.SimulateCandidateDataset(particle_idx,
                                               synthetic_data_file_name_buffer.str(),
                                               candidate_dataset_file_name);
    Real distance = 0.0;
    ComputeDistanceBetweenSyntheticAndCandidateDatasets(synthetic_data_file_name_buffer.str(),
                                                        candidate_dataset_file_name,
                                                        distance);
    if (distance <= tolerance)
    {
      sampled_parameters.push_back(U_0);
    }
  }
  SaveAcceptedParameters(sampled_parameters);
}

void AbcRejectionEngine::RunAbcTowardsExperimentalData()
{

}

void AbcRejectionEngine::SampleNewParameter(std::uniform_real_distribution<Real> &unif_rnd, Real &parameter)
{
  std::mt19937 mersenne_twister_generator(std::random_device{}());
  parameter = unif_rnd(mersenne_twister_generator);
}

void AbcRejectionEngine::ComputeDistanceBetweenSyntheticAndCandidateDatasets(const std::string &synthetic_data_file_name,
                                                                             const std::string &candidate_dataset_file_name,
                                                                             Real &distance)
{
  std::ifstream synthetic_data_file(synthetic_data_file_name, std::ios::binary | std::ios::in);
  assert(synthetic_data_file.is_open());
  synthetic_data_file.seekg(synthetic_parameter_set_.Get_t_0() * (1 + kS * synthetic_parameter_set_.Get_N()) *
      sizeof(Real), std::ios::beg);
  std::ifstream candidate_dataset_file(candidate_dataset_file_name, std::ios::binary | std::ios::in);
  assert(candidate_dataset_file.is_open());

  Real t = synthetic_parameter_set_.Get_t_0();
  Real root_mean_square_error = 0.0, mean_square_error_per_time_unit = 0.0;
  Vector2D r_alpha, r_beta, dr_alpha_beta;
  std::mt19937 mersenne_twister_generator(std::random_device{}());
  std::normal_distribution<Real> norm_dist(0.0, synthetic_parameter_set_.Get_lambda() / 3.0);

  int time_count = 0;
  while (t <= synthetic_parameter_set_.Get_t_1())
  {
    Real tmp = 0.0;
    synthetic_data_file.read((char *) &tmp, sizeof(Real));
    candidate_dataset_file.read((char *) &tmp, sizeof(Real));

    synthetic_data_file.read((char *) &synthetic_data_[0], kS * synthetic_parameter_set_.Get_N() * sizeof(Real));
    // exert observational noise
    std::for_each(synthetic_data_.begin(),
                  synthetic_data_.end(),
                  [&](Real &n)
                  {
                    n += norm_dist(mersenne_twister_generator);
                  });
    candidate_dataset_file.read((char *) &candidate_data_[0], kS * synthetic_parameter_set_.Get_N() * sizeof(Real));

    mean_square_error_per_time_unit = 0.0;
    for (int alpha = 0; alpha < synthetic_parameter_set_.Get_N(); ++alpha)
    {
      r_alpha.x = synthetic_data_[kS * alpha];
      r_alpha.y = synthetic_data_[kS * alpha + 1];
      r_beta.x = candidate_data_[kS * alpha];
      r_beta.y = candidate_data_[kS * alpha + 1];

      pbc_config_.ClassAEffectiveParticleDistance(r_alpha, r_beta, dr_alpha_beta);
      mean_square_error_per_time_unit += NormSquared(dr_alpha_beta);
    }
    root_mean_square_error += mean_square_error_per_time_unit / synthetic_parameter_set_.Get_N();
    t += synthetic_parameter_set_.Get_delta_t() * synthetic_parameter_set_.Get_output_interval();
    ++time_count;
  } // t

  distance = std::sqrt(root_mean_square_error / time_count);
  std::cout << "RMSE:" << distance << std::endl;
}

void AbcRejectionEngine::SaveAcceptedParameters(const std::vector<Real> &sampled_parameters)
{
  std::ostringstream accepted_parameters_file_name_buffer;
  accepted_parameters_file_name_buffer << synthetic_parameter_set_.GetAbcFolderName()
                                       << synthetic_parameter_set_.GetPosteriorDistributionsSubfolderName()
                                       << "accepted_parameters.txt";
  std::remove(accepted_parameters_file_name_buffer.str().c_str());
  std::ofstream accepted_parameters_file(accepted_parameters_file_name_buffer.str(), std::ios::out | std::ios::app);
  for (int particle_idx = 0; particle_idx < sampled_parameters.size(); ++particle_idx)
  {
    accepted_parameters_file << sampled_parameters[particle_idx] << std::endl;
  }
  accepted_parameters_file.close();
}