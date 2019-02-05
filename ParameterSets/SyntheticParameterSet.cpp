//
// Created by Nikita Kruk on 25.11.17.
//

#include "SyntheticParameterSet.hpp"

#include <fstream>
#include <unordered_map>
#include <random>
#include <cassert>
#include <cmath>
#include <boost/math/constants/constants.hpp>

SyntheticParameterSet::SyntheticParameterSet() :
    AbstractParameterSet()
{
#if defined(__linux__) && defined(LICHTENBERG)
  std::ifstream parameters_file("/home/nk59zoce/cpp/sprApproximateBayesianComputation/ParameterSets/SyntheticParameters.cfg", std::ios::in);
#elif defined(__linux__) && defined(BCS_CLUSTER)
  std::ifstream parameters_file("/home/nkruk/cpp/sprApproximateBayesianComputation/ParameterSets/SyntheticParameters.cfg", std::ios::in);
#elif defined(__linux__)
#elif defined(_WIN32)
#elif defined(__APPLE__)
  std::ifstream parameters_file(
      "/Users/nikita/CLionProjects/sprApproximateBayesianComputation/ParameterSets/SyntheticParameters.cfg",
      std::ios::in);
#endif
  assert(parameters_file.is_open());

  //read string values
  parameters_file >> abc_folder_name_;//description of the parameter is omitted
  parameters_file >> abc_folder_name_;
  parameters_file >> synthetic_data_folder_name_;//description of the parameter is omitted
  parameters_file >> synthetic_data_folder_name_;
  parameters_file >> experimental_data_file_name_;//description of the parameter is omitted
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
  N_ = (int) parameters_dictionary["N"];
  a_ = parameters_dictionary["a"];
  U_0_ = parameters_dictionary["U_0"];
  lambda_ = parameters_dictionary["lambda"];
  F_ = parameters_dictionary["F"];
  f_0_ = parameters_dictionary["f_0"];
  phi_ = parameters_dictionary["phi"];
  kappa_ = parameters_dictionary["kappa"];
  burn_in_time_ = parameters_dictionary["burn_in_time"];
  t_0_ = parameters_dictionary["t_0"];
  t_1_ = parameters_dictionary["t_1"];
  delta_t_ = parameters_dictionary["delta_t"];
  output_interval_ = (int) parameters_dictionary["output_interval"];
  D_ = parameters_dictionary["D"];
  T_ = parameters_dictionary["T"];
  eta_S_ = parameters_dictionary["eta_S"];

  std::mt19937 mersenne_twister_generator(std::random_device{}());
  std::bernoulli_distribution bern_dist(1.0 - kappa_);
  active_passive_.resize(N_);
//	std::for_each(active_passive_.begin(), active_passive_.end(), [&](bool &b){ b = bern_dist(mersenne_twister_generator); });
  for (int alpha = 0; alpha < N_; ++alpha)
  {
    active_passive_[alpha] = bern_dist(mersenne_twister_generator);
  }

  CalculateDependentParameters();
}

SyntheticParameterSet::~SyntheticParameterSet()
{

}

void SyntheticParameterSet::CalculateDependentParameters()
{
  const Real pi = boost::math::constants::pi<Real>();

  l_ = lambda_ * a_;
  L_ = std::sqrt((lambda_ * (l_ - lambda_) + pi * lambda_ * lambda_ / 4.0) * N_ / phi_);
  A_ = L_ * L_;
  rho_ = N_ * lambda_ * lambda_ / A_;

  ComputeNumberOfSegments();
  ComputeDimentionlessGeometricFactors();
  ComputeSelfDiffusionCoefficients();
}

void SyntheticParameterSet::ComputeNumberOfSegments()
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

void SyntheticParameterSet::ComputeDimentionlessGeometricFactors()
{
  const Real pi = boost::math::constants::pi<Real>();
  f_par_ = 2.0 * pi / (std::log(a_) - 0.207 + 0.980 / a_ - 0.133 / (a_ * a_));
  f_orth_ = 4.0 * pi / (std::log(a_) + 0.839 + 0.185 / a_ + 0.233 / (a_ * a_));
  f_R_ = pi * a_ * a_ / (3 * (std::log(a_) - 0.662 + 0.917 / a_ - 0.050 / (a_ * a_)));
}

void SyntheticParameterSet::ComputeSelfDiffusionCoefficients()
{
  const Real pi = boost::math::constants::pi<Real>();
  const Real k_B = 1.38064852 * std::pow(10.0, -23.0);
  D_0_ = k_B * T_ / (eta_S_ * l_);
  D_parallel_ = D_0_ / (2.0 * pi) * (std::log(a_) - 0.207 + 0.980 / a_ - 0.133 / (a_ * a_));
  D_perp_ = D_0_ / (4.0 * pi) * (std::log(a_) + 0.839 + 0.185 / a_ + 0.233 / (a_ * a_));
  D_R_ = 3.0 * D_0_ / (pi * l_ * l_) * (std::log(a_) - 0.662 + 0.917 / a_ - 0.050 / (a_ * a_));
}

int SyntheticParameterSet::Get_N()
{
  return N_;
}

Real SyntheticParameterSet::Get_l()
{
  return l_;
}

Real SyntheticParameterSet::Get_L()
{
  return L_;
}

Real SyntheticParameterSet::Get_A()
{
  return A_;
}

int SyntheticParameterSet::Get_n()
{
  return n_;
}

Real SyntheticParameterSet::Get_U_0()
{
  return U_0_;
}

Real SyntheticParameterSet::Get_lambda()
{
  return lambda_;
}

Real SyntheticParameterSet::Get_a()
{
  return a_;
}

Real SyntheticParameterSet::Get_f_0()
{
  return f_0_;
}

Real SyntheticParameterSet::Get_f_par()
{
  return f_par_;
}

Real SyntheticParameterSet::Get_f_orth()
{
  return f_orth_;
}

Real SyntheticParameterSet::Get_f_R()
{
  return f_R_;
}

Real SyntheticParameterSet::Get_F()
{
  return F_;
}

Real SyntheticParameterSet::Get_phi()
{
  return phi_;
}

Real SyntheticParameterSet::Get_kappa()
{
  return kappa_;
}

const std::vector<bool> &SyntheticParameterSet::Get_active_passive()
{
  return active_passive_;
}

Real SyntheticParameterSet::Get_rho()
{
  return rho_;
}

Real SyntheticParameterSet::Get_burn_in_time()
{
  return burn_in_time_;
}

Real SyntheticParameterSet::Get_t_0()
{
  return t_0_;
}

Real SyntheticParameterSet::Get_t_1()
{
  return t_1_;
}

Real SyntheticParameterSet::Get_delta_t()
{
  return delta_t_;
}

int SyntheticParameterSet::Get_output_interval()
{
  return output_interval_;
}

Real SyntheticParameterSet::Get_D()
{
  return D_;
}

Real SyntheticParameterSet::Get_T()
{
  return T_;
}

Real SyntheticParameterSet::Get_eta_S()
{
  return eta_S_;
}

Real SyntheticParameterSet::Get_D_0()
{
  return D_0_;
}

Real SyntheticParameterSet::Get_D_parallel()
{
  return D_parallel_;
}

Real SyntheticParameterSet::Get_D_perp()
{
  return D_perp_;
}

Real SyntheticParameterSet::Get_D_R()
{
  return D_R_;
}

int SyntheticParameterSet::Get_first_image()
{
  return t_0_ / kSecondsPerImage;
}

int SyntheticParameterSet::Get_last_image()
{
  return t_1_ / kSecondsPerImage;
}

const std::string &SyntheticParameterSet::GetAbcFolderName()
{
  return abc_folder_name_;
}

const std::string &SyntheticParameterSet::GetSyntheticDataFolderName()
{
  return synthetic_data_folder_name_;
}

const std::string &SyntheticParameterSet::GetExperimentalDataFileName()
{
  return experimental_data_file_name_;
}

const std::string &SyntheticParameterSet::GetPosteriorDistributionsSubfolderName()
{
  return posterior_distributions_subfolder_name_;
}

const std::string &SyntheticParameterSet::GetCandidateDatasetsSubfolderName()
{
  return candidate_datasets_subfolder_name_;
}

Real SyntheticParameterSet::Get_alpha()
{
  return alpha_;
}

void SyntheticParameterSet::Set_N(int N)
{
  N_ = N;
  CalculateDependentParameters();
}

void SyntheticParameterSet::Set_a(Real a)
{
  a_ = a;
  CalculateDependentParameters();
}

void SyntheticParameterSet::Set_phi(Real phi)
{
  phi_ = phi;
  CalculateDependentParameters();
}

void SyntheticParameterSet::Set_t_0(Real t_0)
{
  t_0_ = t_0;
}

void SyntheticParameterSet::Set_t_1(Real t_1)
{
  t_1_ = t_1;
}

void SyntheticParameterSet::Set_U_0(Real U_0)
{
  U_0_ = U_0;
  CalculateDependentParameters();
}

void SyntheticParameterSet::Set_F(Real F)
{
  F_ = F;
  CalculateDependentParameters();
}

void SyntheticParameterSet::Set_f_0(Real f_0)
{
  f_0_ = f_0;
  CalculateDependentParameters();
}

void SyntheticParameterSet::Set_D(Real D)
{
  D_ = D;
}