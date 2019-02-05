//
// Created by Nikita Kruk on 08.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_COMMONPARAMETERS_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_COMMONPARAMETERS_HPP

#include "../Definitions.hpp"

#include <string>
#include <vector>

class AbstractParameterSet
{
 public:

  AbstractParameterSet();
  virtual ~AbstractParameterSet();

  virtual int Get_N() = 0;
  virtual Real Get_lambda() = 0;
  virtual Real Get_L() = 0;
  virtual int Get_n() = 0;
  virtual Real Get_delta_t() = 0;
  virtual Real Get_l() = 0;
  virtual Real Get_U_0() = 0;
  virtual const std::vector<bool> &Get_active_passive() = 0;
  virtual Real Get_T() = 0;
  virtual Real Get_D_parallel() = 0;
  virtual Real Get_D_perp() = 0;
  virtual Real Get_D_R() = 0;
  virtual int Get_output_interval() = 0;
  virtual Real Get_a() = 0;
  virtual Real Get_kappa() = 0;
  virtual Real Get_phi() = 0;
  virtual Real Get_A() = 0;
  virtual Real Get_rho() = 0;
  virtual Real Get_t_0() = 0;
  virtual Real Get_t_1() = 0;
  virtual Real Get_eta_S() = 0;
  virtual Real Get_D_0() = 0;
  virtual int Get_first_image() = 0;
  virtual int Get_last_image() = 0;

  virtual Real Get_f_0() = 0;
  virtual Real Get_f_par() = 0;
  virtual Real Get_f_orth() = 0;
  virtual Real Get_f_R() = 0;
  virtual Real Get_F() = 0;
  virtual Real Get_D() = 0;

  virtual const std::string &GetAbcFolderName() = 0;
  virtual const std::string &GetSyntheticDataFolderName() = 0;
  virtual const std::string &GetCandidateDatasetsSubfolderName() = 0;

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_COMMONPARAMETERS_HPP
