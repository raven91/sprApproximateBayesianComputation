//
// Created by Nikita Kruk on 25.11.17.
//

#include "SyntheticObserver.hpp"

#include <sstream>
#include <cassert>
#include <iostream>
#include <algorithm> // std::copy
#include <iterator>  // std::ostream_iterator

SyntheticObserver::SyntheticObserver(AbstractParameterSet *parameter_set,
                                     const std::vector<Real> &velocity_distribution) :
    parameter_set_(parameter_set),
    ouput_time_counter_(0),
    output_time_threshold_(parameter_set->Get_output_interval()), // mod 1 - save at every dt
    velocity_distribution_(velocity_distribution)
{
  integration_step_timer_ = std::chrono::system_clock::now();

  data_file_name_ = ConstructDataFileName();
  std::remove(data_file_name_.c_str());
  data_file_.open(data_file_name_, std::ios::binary | std::ios::out | std::ios::app);
  assert(data_file_.is_open());

  std::string parameters_file_name = ConstructParameterFileName();
  SaveAllParametersIntoFile(parameters_file_name);
}

SyntheticObserver::SyntheticObserver(AbstractParameterSet *parameter_set,
                                     const std::vector<Real> &velocity_distribution,
                                     int no_file_creation) :
    parameter_set_(parameter_set),
    ouput_time_counter_(0),
    output_time_threshold_(parameter_set->Get_output_interval()), // mod 1 - save at every dt
    velocity_distribution_(velocity_distribution)
{
  integration_step_timer_ = std::chrono::system_clock::now();
}

SyntheticObserver::~SyntheticObserver()
{
  if (data_file_.is_open())
  {
    data_file_.close();
  }
}

std::string SyntheticObserver::ConstructDataFileName()
{
  std::ostringstream data_file_name_buffer;
  data_file_name_buffer << parameter_set_->GetSyntheticDataFolderName()
                        << "/spr_simulation_N_" << parameter_set_->Get_N()
                        << "_phi_" << parameter_set_->Get_phi()
                        << "_a_" << parameter_set_->Get_a() << "_U0_" << parameter_set_->Get_U_0() << "_k_"
                        << parameter_set_->Get_kappa() << ".bin";
  return data_file_name_buffer.str();
}

void SyntheticObserver::operator()(std::vector<Real> &system_state, Real t)
{
  if (!(ouput_time_counter_ % output_time_threshold_))
  {
    std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
    std::cout << "t=" << t << " integrated in " << elapsed_seconds.count() << "s" << std::endl;
    integration_step_timer_ = std::chrono::system_clock::now();

    data_file_.write((char *) &t, sizeof(Real));
    std::vector<Real> system_state_with_arbitrary_velocities(system_state);
    for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
    {
      system_state_with_arbitrary_velocities[kS * alpha + 2] *= velocity_distribution_[alpha];
      system_state_with_arbitrary_velocities[kS * alpha + 3] *= velocity_distribution_[alpha];
    }
    data_file_.write((char *) &system_state_with_arbitrary_velocities[0], kS * parameter_set_->Get_N() * sizeof(Real));
  }
  ++ouput_time_counter_;
}

std::string SyntheticObserver::ConstructParameterFileName()
{
  std::ostringstream parameters_file_name_buffer;
  parameters_file_name_buffer << parameter_set_->GetSyntheticDataFolderName()
                              << "/parameters_N_" << parameter_set_->Get_N()
                              << "_phi_" << parameter_set_->Get_phi()
                              << "_a_" << parameter_set_->Get_a() << "_U0_" << parameter_set_->Get_U_0() << "_k_"
                              << parameter_set_->Get_kappa() << ".txt";
  return parameters_file_name_buffer.str();
}

std::string SyntheticObserver::ConstructActivePassiveFileName()
{
  std::ostringstream active_passive_file_name_buffer;
  active_passive_file_name_buffer << parameter_set_->GetSyntheticDataFolderName()
                                  << "/active_passive_N_" << parameter_set_->Get_N()
                                  << "_phi_" << parameter_set_->Get_phi()
                                  << "_a_" << parameter_set_->Get_a() << "_U0_" << parameter_set_->Get_U_0() << "_k_"
                                  << parameter_set_->Get_kappa() << ".txt";
  return active_passive_file_name_buffer.str();
}

void SyntheticObserver::SaveAllParametersIntoFile(std::string &parameters_file_name)
{
  std::ofstream parameters_file(parameters_file_name, std::ios::out | std::ios::trunc);

  parameters_file << "N " << parameter_set_->Get_N() << std::endl;
  parameters_file << "l " << parameter_set_->Get_l() << std::endl;
  parameters_file << "L " << parameter_set_->Get_L() << std::endl;
  parameters_file << "A " << parameter_set_->Get_A() << std::endl;
  parameters_file << "n " << parameter_set_->Get_n() << std::endl;
  parameters_file << "U_0 " << parameter_set_->Get_U_0() << std::endl;
  parameters_file << "lambda " << parameter_set_->Get_lambda() << std::endl;
  parameters_file << "a " << parameter_set_->Get_a() << std::endl;
  parameters_file << "f_0 " << parameter_set_->Get_f_0() << std::endl;
  parameters_file << "f_par " << parameter_set_->Get_f_par() << std::endl;
  parameters_file << "f_orth " << parameter_set_->Get_f_orth() << std::endl;
  parameters_file << "f_R " << parameter_set_->Get_f_R() << std::endl;
  parameters_file << "F " << parameter_set_->Get_F() << std::endl;
  parameters_file << "phi " << parameter_set_->Get_phi() << std::endl;
  parameters_file << "kappa " << parameter_set_->Get_kappa() << std::endl;
  parameters_file << "rho " << parameter_set_->Get_rho() << std::endl;
  parameters_file << "t_0 " << parameter_set_->Get_t_0() << std::endl;
  parameters_file << "t_1 " << parameter_set_->Get_t_1() << std::endl;
  parameters_file << "delta_t " << parameter_set_->Get_delta_t() << std::endl;
  parameters_file << "output_interval " << parameter_set_->Get_output_interval() << std::endl;
  parameters_file << "output_delta_t " << parameter_set_->Get_delta_t() * output_time_threshold_ << std::endl;
  parameters_file << "eta_S " << parameter_set_->Get_eta_S() << std::endl;
  parameters_file << "T " << parameter_set_->Get_T() << std::endl;
  parameters_file << "D_0 " << parameter_set_->Get_D_0() << std::endl;
  parameters_file << "D_parallel " << parameter_set_->Get_D_parallel() << std::endl;
  parameters_file << "D_perp " << parameter_set_->Get_D_perp() << std::endl;
  parameters_file << "D_R " << parameter_set_->Get_D_R() << std::endl;
  parameters_file.close();

//  std::string active_passive_file_name = ConstructActivePassiveFileName();
//  std::ofstream active_passive_file(active_passive_file_name, std::ios::out | std::ios::trunc);
//  std::copy(parameter_set_.Get_active_passive().begin(),
//            parameter_set_.Get_active_passive().end(),
//            std::ostream_iterator<float>(active_passive_file, " "));
//  active_passive_file.close();
}