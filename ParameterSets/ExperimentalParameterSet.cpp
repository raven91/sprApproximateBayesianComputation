//
// Created by Nikita Kruk on 25.11.17.
//

#include "ExperimentalParameterSet.hpp"

#include <fstream>
#include <unordered_map>
#include <random>
#include <cassert>
#include <cmath>
#include <iostream>
#include <boost/math/constants/constants.hpp>

ExperimentalParameterSet::ExperimentalParameterSet() :
    AbstractParameterSet()
{
#if defined(__linux__) && defined(LICHTENBERG)
  std::ifstream parameters_file("/home/nk59zoce/cpp/sprApproximateBayesianComputation/ParameterSets/ExperimentalParameters.cfg", std::ios::in);
#elif defined(__linux__) && defined(BCS_CLUSTER)
  std::ifstream parameters_file("/home/nkruk/cpp/sprApproximateBayesianComputation/ParameterSets/ExperimentalParameters.cfg", std::ios::in);
#elif defined(__linux__)
#elif defined(_WIN32)
#elif defined(__APPLE__)
  std::ifstream parameters_file(
      "/Users/nikita/CLionProjects/sprApproximateBayesianComputation/ParameterSets/ExperimentalParameters.cfg",
      std::ios::in);
#endif
  assert(parameters_file.is_open());

  //read string values
  parameters_file >> abc_folder_name_;//description of the parameter is omitted
  parameters_file >> abc_folder_name_;
  parameters_file >> synthetic_data_folder_name_;//description of the parameter is omitted
  parameters_file >> synthetic_data_folder_name_;
  parameters_file >> experimental_data_folder_name_;//description of the parameter is omitted
  parameters_file >> experimental_data_folder_name_;
  parameters_file >> experimental_configuration_file_name_;// description of the parameter is omitted
  parameters_file >> experimental_configuration_file_name_;
  parameters_file >> experimental_data_file_name_;// description of the parameter is omitted
  parameters_file >> experimental_data_file_name_;
  parameters_file >> posterior_distributions_subfolder_name_;//description of the parameter is omitted
  parameters_file >> posterior_distributions_subfolder_name_;
  parameters_file >> candidate_datasets_subfolder_name_;//description of the parameter is omitted
  parameters_file >> candidate_datasets_subfolder_name_;

  //read Real values
  std::unordered_map<std::string, Real> parameters_dictionary;
  std::string key;
  Real value;
  while (parameters_file >> key >> value)
  {
    parameters_dictionary[key] = value;
  }
  parameters_file.close();

  //initialize class independent parameters
  a_ = parameters_dictionary["a"];
  lambda_ = parameters_dictionary["lambda"];
  delta_t_ = parameters_dictionary["delta_t"];
  output_interval_ = (int) parameters_dictionary["output_interval"];
  T_ = parameters_dictionary["T"];
  eta_S_ = parameters_dictionary["eta_S"];
  kappa_ = parameters_dictionary["kappa"];
  F_ = parameters_dictionary["F"];
  f_0_ = parameters_dictionary["f_0"];
  D_ = parameters_dictionary["D"];

  // read experimenal configuration
  std::ifstream experimental_configuration_file
      (experimental_data_folder_name_ + experimental_configuration_file_name_, std::ios::in);
  assert(experimental_configuration_file.is_open());
  std::unordered_map<std::string, std::string> experimental_configuration;
  std::string exp_key, exp_value;
  while (experimental_configuration_file >> exp_key >> exp_value)
  {
    experimental_configuration[exp_key] = exp_value;
  }
  experimental_configuration_file.close();
  first_image_ = (int) std::stod(experimental_configuration["first_image"]);
  last_image_ = (int) std::stod(experimental_configuration["last_image"]);
  t_0_ = first_image_ * kSecondsPerImage;
  t_1_ = last_image_ * kSecondsPerImage;
  L_ = std::stod(experimental_configuration["subimage_x_size"]) * kMicroMetersPerPixel * kMetersPerMicroMeter;
  l_ = lambda_ * a_;
  A_ = L_ * L_;

  InitializeDependentParameters();
}

ExperimentalParameterSet::~ExperimentalParameterSet()
{

}

void ExperimentalParameterSet::InitializeDependentParameters()
{
  N_ = 0;
  phi_ = 0.0;
  rho_ = 0.0;

  ComputeNumberOfSegments();
  ComputeDimentionlessGeometricFactors();
  ComputeSelfDiffusionCoefficients();
}

void ExperimentalParameterSet::RecalculateDependentParameters()
{
  const Real pi = boost::math::constants::pi<Real>();

  phi_ = (lambda_ * (l_ - lambda_) + pi * lambda_ * lambda_ / 4.0) * N_ / A_;
  rho_ = N_ * lambda_ * lambda_ / A_;
  active_passive_.resize(N_, 1.0);

  ComputeNumberOfSegments();
  ComputeDimentionlessGeometricFactors();
  ComputeSelfDiffusionCoefficients();
}

void ExperimentalParameterSet::ComputeNumberOfSegments()
{
  assert(a_ >= Real(1.0));

  if (a_ == Real(1.0))
  {
    n_ = 1;
  } else if (a_ <= Real(3.0))
  {
    n_ = 3;
  } else
  {
    n_ = (int) std::nearbyint(9.0 * a_ / 8.0);
  }
}

void ExperimentalParameterSet::ComputeDimentionlessGeometricFactors()
{
  const Real pi = boost::math::constants::pi<Real>();
  f_par_ = 2.0 * pi / (std::log(a_) - 0.207 + 0.980 / a_ - 0.133 / (a_ * a_));
  f_orth_ = 4.0 * pi / (std::log(a_) + 0.839 + 0.185 / a_ + 0.233 / (a_ * a_));
  f_R_ = pi * a_ * a_ / (3 * (std::log(a_) - 0.662 + 0.917 / a_ - 0.050 / (a_ * a_)));
}

void ExperimentalParameterSet::ComputeSelfDiffusionCoefficients()
{
  const Real pi = boost::math::constants::pi<Real>();
  const Real k_B = 1.38064852 * std::pow(10.0, -23.0);
  D_0_ = k_B * T_ / (eta_S_ * l_);
  D_parallel_ = D_0_ / (2.0 * pi) * (std::log(a_) - 0.207 + 0.980 / a_ - 0.133 / (a_ * a_));
  D_perp_ = D_0_ / (4.0 * pi) * (std::log(a_) + 0.839 + 0.185 / a_ + 0.233 / (a_ * a_));
  D_R_ = 3.0 * D_0_ / (pi * l_ * l_) * (std::log(a_) - 0.662 + 0.917 / a_ - 0.050 / (a_ * a_));
}

int ExperimentalParameterSet::Get_N()
{
  return N_;
}

Real ExperimentalParameterSet::Get_l()
{
  return l_;
}

Real ExperimentalParameterSet::Get_L()
{
  return L_;
}

Real ExperimentalParameterSet::Get_A()
{
  return A_;
}

int ExperimentalParameterSet::Get_n()
{
  return n_;
}

Real ExperimentalParameterSet::Get_U_0()
{
  return U_0_;
}

Real ExperimentalParameterSet::Get_lambda()
{
  return lambda_;
}

Real ExperimentalParameterSet::Get_a()
{
  return a_;
}

Real ExperimentalParameterSet::Get_f_0()
{
  return f_0_;
}

Real ExperimentalParameterSet::Get_f_par()
{
  return f_par_;
}

Real ExperimentalParameterSet::Get_f_orth()
{
  return f_orth_;
}

Real ExperimentalParameterSet::Get_f_R()
{
  return f_R_;
}

Real ExperimentalParameterSet::Get_F()
{
  return F_;
}

Real ExperimentalParameterSet::Get_phi()
{
  return phi_;
}

Real ExperimentalParameterSet::Get_kappa()
{
  return kappa_;
}

const std::vector<bool> &ExperimentalParameterSet::Get_active_passive()
{
  return active_passive_;
}

Real ExperimentalParameterSet::Get_rho()
{
  return rho_;
}

Real ExperimentalParameterSet::Get_burn_in_time()
{
  return burn_in_time_;
}

Real ExperimentalParameterSet::Get_t_0()
{
  return t_0_;
}

Real ExperimentalParameterSet::Get_t_1()
{
  return t_1_;
}

Real ExperimentalParameterSet::Get_delta_t()
{
  return delta_t_;
}

int ExperimentalParameterSet::Get_output_interval()
{
  return output_interval_;
}

Real ExperimentalParameterSet::Get_D()
{
  return D_;
}

Real ExperimentalParameterSet::Get_T()
{
  return T_;
}

Real ExperimentalParameterSet::Get_eta_S()
{
  return eta_S_;
}

Real ExperimentalParameterSet::Get_D_0()
{
  return D_0_;
}

Real ExperimentalParameterSet::Get_D_parallel()
{
  return D_parallel_;
}

Real ExperimentalParameterSet::Get_D_perp()
{
  return D_perp_;
}

Real ExperimentalParameterSet::Get_D_R()
{
  return D_R_;
}

int ExperimentalParameterSet::Get_first_image()
{
  return first_image_;
}

int ExperimentalParameterSet::Get_last_image()
{
  return last_image_;
}

const std::string &ExperimentalParameterSet::GetAbcFolderName()
{
  return abc_folder_name_;
}

const std::string &ExperimentalParameterSet::GetSyntheticDataFolderName()
{
  return synthetic_data_folder_name_;
}

const std::string &ExperimentalParameterSet::GetExperimentalDataFolderName()
{
  return experimental_data_folder_name_;
}

const std::string &ExperimentalParameterSet::GetExperimentalDataFileName()
{
  return experimental_data_file_name_;
}

const std::string &ExperimentalParameterSet::GetPosteriorDistributionsSubfolderName()
{
  return posterior_distributions_subfolder_name_;
}

const std::string &ExperimentalParameterSet::GetCandidateDatasetsSubfolderName()
{
  return candidate_datasets_subfolder_name_;
}

void ExperimentalParameterSet::Set_N(int N)
{
  N_ = N;
  RecalculateDependentParameters();
}

void ExperimentalParameterSet::Set_a(Real a)
{
  a_ = a;
  RecalculateDependentParameters();
}

void ExperimentalParameterSet::Set_t_0(Real t_0)
{
  t_0_ = t_0;
}

void ExperimentalParameterSet::Set_t_1(Real t_1)
{
  t_1_ = t_1;
}

void ExperimentalParameterSet::Set_U_0(Real U_0)
{
  U_0_ = U_0;
  RecalculateDependentParameters();
}

void ExperimentalParameterSet::Set_first_image(int first_image)
{
  first_image_ = first_image;
  t_0_ = first_image_ * kSecondsPerImage;
}

void ExperimentalParameterSet::Set_last_image(int last_image)
{
  last_image_ = last_image;
  t_1_ = last_image_ * kSecondsPerImage;
}