//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SYNTHETICPARAMETERSET_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SYNTHETICPARAMETERSET_HPP

#include "../Definitions.hpp"
#include "AbstractParameterSet.hpp"

#include <string>
#include <vector>

class SyntheticParameterSet : public AbstractParameterSet
{
 public:

  SyntheticParameterSet();
  ~SyntheticParameterSet();

  int Get_N();
  Real Get_l();
  Real Get_L();
  Real Get_A();
  int Get_n();
  Real Get_U_0();
  Real Get_lambda();
  Real Get_a();
  Real Get_f_0();
  Real Get_f_par();
  Real Get_f_orth();
  Real Get_f_R();
  Real Get_F();
  Real Get_phi();
  Real Get_kappa();
  const std::vector<bool> &Get_active_passive();
  Real Get_rho();
  Real Get_burn_in_time();
  Real Get_t_0();
  Real Get_t_1();
  Real Get_delta_t();
  int Get_output_interval();
  Real Get_D();
  Real Get_T();
  Real Get_eta_S();
  Real Get_D_0();
  Real Get_D_parallel();
  Real Get_D_perp();
  Real Get_D_R();
  int Get_first_image();
  int Get_last_image();
  const std::string &GetAbcFolderName();
  const std::string &GetSyntheticDataFolderName();
  const std::string &GetExperimentalDataFileName();
  const std::string &GetPosteriorDistributionsSubfolderName();
  const std::string &GetCandidateDatasetsSubfolderName();

  Real Get_alpha();

  void Set_N(int N);
  void Set_a(Real a);
  void Set_phi(Real phi);
  void Set_t_0(Real t_0);
  void Set_t_1(Real t_1);
  void Set_U_0(Real U_0);
  void Set_F(Real F);
  void Set_f_0(Real f_0);
  void Set_D(Real D);

 private:

  int N_;//number of rods in the system
  Real l_;//length of a rod
  Real L_;//size of a simulation box in one dimention
  Real A_;//area of a simulation box
  int n_;//number of segments per rod
  Real U_0_;//potential amplitude
  Real lambda_;//screening length
  Real a_;//aspect ratio between length and width of a rod
  Real f_0_;//Stokesian friction coefficient
  Real f_par_; //
  Real f_orth_;//geometric factors for rod-like macromolecules
  Real f_R_;   //
  Real F_;//self-motility force
  Real phi_;//effective volume fraction
  Real kappa_;//fraction of passive particles
  std::vector<bool> active_passive_;//active or passive bit
  Real rho_;//rod density
  Real burn_in_time_;
  Real t_0_;//start of integration
  Real t_1_;//end of integration
  Real delta_t_;//integration time-step
  int output_interval_;//frequency of integration output
  Real D_;// diffusion constant

  Real eta_S_;// dynamic viscosity
  Real T_;// temperature
  Real D_0_;
  Real D_parallel_;
  Real D_perp_;
  Real D_R_;

  Real alpha_;//tunable parameter for clusterization

  std::string abc_folder_name_;
  std::string synthetic_data_folder_name_;
  std::string experimental_data_file_name_;
  std::string posterior_distributions_subfolder_name_;
  std::string candidate_datasets_subfolder_name_;

  void CalculateDependentParameters();
  void ComputeNumberOfSegments();
  void ComputeDimentionlessGeometricFactors();
  void ComputeSelfDiffusionCoefficients();

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SYNTHETICPARAMETERSET_HPP
