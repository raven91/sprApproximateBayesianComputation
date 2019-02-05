//
// Created by Nikita Kruk on 08.08.18.
//

#include "AbcSmcEngineExperimental.hpp"
#include "SimulationEngineForExperimentalData.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>    // std::for_each

#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/normal.hpp>

AbcSmcEngineExperimental::AbcSmcEngineExperimental() :
    experimental_parameter_set_(),
    fbc_config_(experimental_parameter_set_.Get_L(), experimental_parameter_set_.Get_L())
{

}

AbcSmcEngineExperimental::~AbcSmcEngineExperimental()
{

}

void AbcSmcEngineExperimental::RunAbcTowardsExperimentalData()
{
  std::ostringstream experimental_data_file_name_buffer;
  experimental_data_file_name_buffer << experimental_parameter_set_.GetExperimentalDataFolderName()
                                     << experimental_parameter_set_.GetExperimentalDataFileName();
  experimental_parameter_set_.Set_first_image(11);
  experimental_parameter_set_.Set_last_image(20);
  RetrieveNumberOfAgentsAndInitialIndexes(experimental_data_file_name_buffer.str());

  // Initialize tolerances
  int number_of_populations = 5;
  std::vector<Real> tolerances(number_of_populations);
  tolerances[0] = 2.0 * 1e-06;
  tolerances[1] = 1.75 * 1e-06;
  tolerances[2] = 1.7 * 1e-06;
  tolerances[3] = 1.65 * 1e-06;
  tolerances[4] = 1.6 * 1e-06;

  int population_size = 1;
  std::uniform_real_distribution<Real> unif_rnd(0.0, 1e-18);
  boost::math::uniform_distribution<Real> unif_pdf(0.0, 1e-18);
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
//      new_parameter = 1e-19;
      // if \pi(\theta^{**})=0, return to S2.1
      if (boost::math::pdf(unif_pdf, new_parameter) == 0.0)
      {
        continue;
      }

      // simulate candidate dataset
      int number_of_simulations = 5;
      int number_of_successful_simulations = 0;
      for (int simulation_idx = 0; simulation_idx < number_of_simulations; ++simulation_idx)
      {
        experimental_parameter_set_.Set_U_0(new_parameter);
        SimulationEngineForExperimentalData simulation_engine(experimental_parameter_set_, fbc_config_);
        std::string candidate_dataset_file_name;
        simulation_engine.SimulateCandidateDataset(population_idx,
                                                   experimental_data_file_name_buffer.str(),
                                                   candidate_dataset_file_name, initial_indexes_);
        Real distance = 0.0;
        ComputeDistanceBetweenExperimentalAndCandidateDatasets(experimental_data_file_name_buffer.str(),
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
              perturbation_kernel(sampled_parameters[j], experimental_parameter_set_.Get_U_0() / 3.0);
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

void AbcSmcEngineExperimental::RetrieveNumberOfAgentsAndInitialIndexes(const std::string &file_name)
{
  std::ifstream experimental_data_file(file_name, std::ios::in);
  assert(experimental_data_file.is_open());
  std::string line;
  for (int i = 0; i < experimental_parameter_set_.Get_first_image(); ++i)
  {
    std::getline(experimental_data_file, line);
  }

  int time_idx = 0;
  int agent_idx = 0;
  int number_of_agents = 0;
  Eigen::VectorXd agent(kSe);
  experimental_data_file >> time_idx >> number_of_agents;
  candidate_data_.resize(kS * number_of_agents, 0.0);
  initial_indexes_.resize(number_of_agents, 0);
  for (int alpha = 0; alpha < number_of_agents; ++alpha)
  {
    experimental_data_file >> agent_idx >> agent(0) >> agent(1) >> agent(2) >> agent(3)
                           >> agent(4) >> agent(5) >> agent(6) >> agent(7);
    initial_indexes_[alpha] = agent_idx;
  }
  experimental_data_file.close();
  experimental_parameter_set_.Set_N(number_of_agents);
}

void AbcSmcEngineExperimental::SampleParameterFromPreviousPopulation(int population_indicator,
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
    std::normal_distribution<Real> perturbation_kernel(new_parameter, 1e-20 / 3.0);
    new_parameter = perturbation_kernel(mersenne_twister_generator);
  }
}

void AbcSmcEngineExperimental::ComputeDistanceBetweenExperimentalAndCandidateDatasets(const std::string &experimental_data_file_name,
                                                                                      const std::string &candidate_dataset_file_name,
                                                                                      Real &distance)
{
  std::ifstream experimental_data_file(experimental_data_file_name, std::ios::in);
  assert(experimental_data_file.is_open());
  std::string line;
  for (int i = 0; i < experimental_parameter_set_.Get_first_image(); ++i)
  {
    std::getline(experimental_data_file, line);
  }
  std::ifstream candidate_dataset_file(candidate_dataset_file_name, std::ios::binary | std::ios::in);
  assert(candidate_dataset_file.is_open());

  Real t = experimental_parameter_set_.Get_t_0();
  Real root_mean_square_error = 0.0;
//  std::mt19937 mersenne_twister_generator(std::random_device{}());
//  std::normal_distribution<Real> norm_dist(0.0, experimental_parameter_set_.Get_lambda() / 3.0);

  int time_count = 0;
  while (t <= experimental_parameter_set_.Get_t_1())
  {
    experimental_data_.clear();

    int time_idx = 0;
    int number_of_agents = 0;
    int agent_idx = 0;
    experimental_data_file >> time_idx >> number_of_agents;
    for (int alpha = 0; alpha < number_of_agents; ++alpha)
    {
      Eigen::VectorXd r_ext(kSe);
      experimental_data_file >> agent_idx >> r_ext(0) >> r_ext(1) >> r_ext(2) >> r_ext(3)
                             >> r_ext(4) >> r_ext(5) >> r_ext(6) >> r_ext(7);
      r_ext(0) *= kMicroMetersPerPixel * kMetersPerMicroMeter;
      r_ext(1) *= kMicroMetersPerPixel * kMetersPerMicroMeter;
      r_ext(2) *= kMicroMetersPerPixel * kMetersPerMicroMeter / kSecondsPerImage;
      r_ext(3) *= kMicroMetersPerPixel * kMetersPerMicroMeter / kSecondsPerImage;
      experimental_data_[agent_idx] = r_ext.head(kS);
    }
    Real tmp = 0.0;
    candidate_dataset_file.read((char *) &tmp, sizeof(Real));
    candidate_dataset_file.read((char *) &candidate_data_[0], kS * experimental_parameter_set_.Get_N() * sizeof(Real));

    Real mean_square_error_per_time_unit = 0.0;
    int mean_square_error_per_time_unit_times = 0;
    for (int i_alpha = 0; i_alpha < initial_indexes_.size(); ++i_alpha)
    {
      int alpha = initial_indexes_[i_alpha];
      Eigen::VectorXd r_alpha(kS);
      r_alpha(0) = candidate_data_[kS * i_alpha];
      r_alpha(1) = candidate_data_[kS * i_alpha + 1];
      r_alpha(2) = candidate_data_[kS * i_alpha + 2];
      r_alpha(3) = candidate_data_[kS * i_alpha + 3];
      // TODO: Normalized velocities?
      if (r_alpha(0) >= 0.0 && r_alpha(0) < experimental_parameter_set_.Get_L()
          && r_alpha(1) >= 0.0 && r_alpha(1) < experimental_parameter_set_.Get_L())
      {
        if (experimental_data_.find(alpha) != experimental_data_.end())
        {
          Eigen::VectorXd r_beta(kS);
          r_beta = experimental_data_[alpha];
          Eigen::VectorXd dr_alpha_beta(2);
          dr_alpha_beta = r_beta.head(2) - r_alpha.head(2);
          mean_square_error_per_time_unit += dr_alpha_beta.squaredNorm();
//          mean_square_error_per_time_unit += dr_alpha_beta.lpNorm<1>();
//          mean_square_error_per_time_unit += 1.0 / dr_alpha_beta.norm();
          ++mean_square_error_per_time_unit_times;
        }
      }
    } // i_alpha
    if (mean_square_error_per_time_unit_times != 0)
    {
      root_mean_square_error += mean_square_error_per_time_unit / mean_square_error_per_time_unit_times;
    }
    t += experimental_parameter_set_.Get_delta_t() * experimental_parameter_set_.Get_output_interval();
    ++time_count;
  } // t

  distance = std::sqrt(root_mean_square_error / time_count);
//  distance = root_mean_square_error / time_count;
//  distance = 1.0 / (root_mean_square_error / time_count);
  std::cout << "RMSE:" << distance << "m|" << distance / kMetersPerMicroMeter / kMicroMetersPerPixel << "px"
            << std::endl;
}

void AbcSmcEngineExperimental::NormalizeWeights(std::vector<Real> &weights)
{
  Real sum = std::accumulate(weights.begin(), weights.end(), 0.0);
  for (Real &w : weights)
  {
    w /= sum;
  }
}

void AbcSmcEngineExperimental::SaveAcceptedParameters(int population,
                                                      const std::vector<Real> &sampled_parameters,
                                                      const std::vector<Real> &weights)
{
  std::ostringstream accepted_parameters_file_name_buffer;
  accepted_parameters_file_name_buffer << experimental_parameter_set_.GetAbcFolderName()
                                       << experimental_parameter_set_.GetPosteriorDistributionsSubfolderName()
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