//
// Created by Nikita Kruk on 03.12.17.
//

#include "AbcSmcEngineSynthetic.hpp"
#include "SimulationEngineForSyntheticData.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>    // std::for_each

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/normal.hpp>

AbcSmcEngineSynthetic::AbcSmcEngineSynthetic() :
    synthetic_parameter_set_(),
    pbc_config_(synthetic_parameter_set_.Get_L(), synthetic_parameter_set_.Get_L()),
    synthetic_data_(std::vector<Real>(kS * synthetic_parameter_set_.Get_N(), 0.0)),
    candidate_data_(std::vector<Real>(kS * synthetic_parameter_set_.Get_N(), 0.0))
{
  std::cout << "Approximate Bayesian Computation - Sequential Monte Carlo Started" << std::endl;
}

AbcSmcEngineSynthetic::~AbcSmcEngineSynthetic()
{
  std::cout << "Approximate Bayesian Computation - Sequential Monte Carlo Complete" << std::endl;
}

void AbcSmcEngineSynthetic::PrepareSyntheticData()
{
  SimulationEngineForSyntheticData simulation_engine(synthetic_parameter_set_, pbc_config_);
  simulation_engine.GenerateSyntheticDataWithArbitraryStationaryVelocity();
}

void AbcSmcEngineSynthetic::RunAbcTowardsSyntheticData()
{
  std::ostringstream synthetic_data_file_name_buffer;
  synthetic_data_file_name_buffer << synthetic_parameter_set_.GetSyntheticDataFolderName()
                                  << "/spr_simulation_N_" << synthetic_parameter_set_.Get_N()
                                  << "_phi_" << synthetic_parameter_set_.Get_phi()
                                  << "_a_" << synthetic_parameter_set_.Get_a()
                                  << "_U0_" << synthetic_parameter_set_.Get_U_0()
                                  << "_k_" << synthetic_parameter_set_.Get_kappa() << ".bin";
  synthetic_parameter_set_.Set_t_0(0.0);
  synthetic_parameter_set_.Set_t_1(1.0);

  // Initialize tolerances
  int number_of_populations = 4;
  std::vector<Real> tolerances(number_of_populations);
  tolerances[0] = 12.5 * 1e-6;
  tolerances[1] = 12 * 1e-6;
  tolerances[2] = 11.5 * 1e-6;
  tolerances[3] = 11.4 * 1e-6;

  int population_size = 25;
  std::uniform_real_distribution<Real> unif_rnd(0.0, 1e-19);
  boost::math::uniform_distribution<Real> unif_pdf(0.0, 1e-19);
  std::vector<Real> sampled_parameters(population_size);
  std::vector<Real> weights(population_size, 1.0 / population_size);
  int population_idx = 0;
  while (population_idx < number_of_populations)
  {
    std::vector<Real> new_sampled_parameters(sampled_parameters);
    std::vector<Real> new_weights(weights);

    int particle_idx = 0;
    while (particle_idx < population_size)
    {
      // S2.1 Sample new parameters/particles
      Real new_parameter = 0;
      SampleParameterFromPreviousPopulation(population_idx, unif_rnd, weights, sampled_parameters, new_parameter);
      // if \pi(\theta^{**})=0, return to S2.1
      if (boost::math::pdf(unif_pdf, new_parameter) == 0.0)
      {
        continue;
      }

      // simulate candidate dataset
      int number_of_simulations = 1;
      int number_of_successful_simulations = 0;
      for (int simulation_idx = 0; simulation_idx < number_of_simulations; ++simulation_idx)
      {
        synthetic_parameter_set_.Set_U_0(new_parameter);
        SimulationEngineForSyntheticData simulation_engine(synthetic_parameter_set_, pbc_config_);
        std::string candidate_dataset_file_name;
        simulation_engine.SimulateCandidateDataset(population_idx,
                                                   synthetic_data_file_name_buffer.str(),
                                                   candidate_dataset_file_name);

        Real distance = 0.0;
        ComputeDistanceBetweenSyntheticAndCandidateDatasets(synthetic_data_file_name_buffer.str(),
                                                            candidate_dataset_file_name,
                                                            distance);
        if (distance <= tolerances[population_idx])
        {
          ++number_of_successful_simulations;
        }
      } // simulation_idx
      if (number_of_successful_simulations == 0)
      {
        continue;
      }
      std::cout << "population:" << population_idx << ", particle:" << particle_idx << std::endl;

      // set the new particle
      new_sampled_parameters[particle_idx] = new_parameter;
      // calculate the weight for that particle
      Real weight = 0.0;
      if (population_idx == 0)
      {
        weight = number_of_successful_simulations;
      } else
      {
        for (int j = 0; j < population_size; ++j)
        {
          boost::math::normal_distribution<Real>
              perturbation_kernel(sampled_parameters[j], synthetic_parameter_set_.Get_U_0() / 3.0);
          weight += weights[j] * boost::math::pdf(perturbation_kernel, new_parameter);
        }
        weight = boost::math::pdf(unif_pdf, new_parameter) * number_of_successful_simulations / weight;
      }
      new_weights[particle_idx] = weight;

      ++particle_idx;
    } // particle_idx

    sampled_parameters = new_sampled_parameters;
    weights = new_weights;
    NormalizeWeights(weights);
    SaveAcceptedParameters(population_idx, sampled_parameters, weights);
    ++population_idx;
  } // population_idx
}

void AbcSmcEngineSynthetic::ComputeDistanceBetweenSyntheticAndCandidateDatasets(const std::string &synthetic_data_file_name,
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
  Real root_mean_square_error = 0.0;
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

    Real mean_square_error_per_time_unit = 0.0;
    for (int alpha = 0; alpha < synthetic_parameter_set_.Get_N(); ++alpha)
    {
      r_alpha.x = synthetic_data_[kS * alpha];
      r_alpha.y = synthetic_data_[kS * alpha + 1];
      r_beta.x = candidate_data_[kS * alpha];
      r_beta.y = candidate_data_[kS * alpha + 1];

      pbc_config_.ClassAEffectiveParticleDistance(r_alpha, r_beta, dr_alpha_beta);
      mean_square_error_per_time_unit += NormSquared(dr_alpha_beta);
//      mean_square_error_per_time_unit += (synthetic_parameter_set_.Get_L() / 2.0 * M_SQRT2) * (synthetic_parameter_set_.Get_L() / 2.0 * M_SQRT2);
//      mean_square_error_per_time_unit += std::pow(synthetic_parameter_set_.Get_L() / 2.0 * M_SQRT2, -1.0);
    } // alpha
    root_mean_square_error += mean_square_error_per_time_unit / synthetic_parameter_set_.Get_N();
    t += synthetic_parameter_set_.Get_delta_t() * synthetic_parameter_set_.Get_output_interval();
    ++time_count;
  } // t

  distance = std::sqrt(root_mean_square_error / time_count);
//  distance = std::pow(root_mean_square_error / time_count, -1.0);
  std::cout << "RMSE:" << distance << std::endl;
}

void AbcSmcEngineSynthetic::SaveAcceptedParameters(int population,
                                          const std::vector<Real> &sampled_parameters,
                                          const std::vector<Real> &weights)
{
  std::ostringstream accepted_parameters_file_name_buffer;
  accepted_parameters_file_name_buffer << synthetic_parameter_set_.GetAbcFolderName()
                                       << synthetic_parameter_set_.GetPosteriorDistributionsSubfolderName()
                                       << "accepted_parameters_"
                                       << population << ".txt";
  std::remove(accepted_parameters_file_name_buffer.str().c_str());
  std::ofstream accepted_parameters_file(accepted_parameters_file_name_buffer.str(), std::ios::out | std::ios::app);
  for (int particle_idx = 0; particle_idx < sampled_parameters.size(); ++particle_idx)
  {
    accepted_parameters_file << sampled_parameters[particle_idx] << '\t' << weights[particle_idx] << std::endl;
  }
  accepted_parameters_file.close();
}

void AbcSmcEngineSynthetic::SampleParameterFromPreviousPopulation(int population_indicator,
                                                         std::uniform_real_distribution<Real> &unif_rnd,
                                                         const std::vector<Real> &weights,
                                                         const std::vector<Real> &sampled_parameters,
                                                         Real &new_parameter)
{
  std::mt19937 mersenne_twister_generator(std::random_device{}());
  if (population_indicator == 0)
  {
    new_parameter = unif_rnd(mersenne_twister_generator);
  } else
  {
    std::discrete_distribution<int> previous_population_distribution(weights.begin(),
                                                                     weights.end());
    int new_parameter_index = previous_population_distribution(mersenne_twister_generator);
    new_parameter = sampled_parameters[new_parameter_index];
    std::normal_distribution<Real> perturbation_kernel(new_parameter, synthetic_parameter_set_.Get_U_0() / 3.0);
    new_parameter = perturbation_kernel(mersenne_twister_generator);
  }
}

void AbcSmcEngineSynthetic::NormalizeWeights(std::vector<Real> &weights)
{
  Real sum = std::accumulate(weights.begin(), weights.end(), 0.0);
  for (Real &w : weights)
  {
    w /= sum;
  }
}