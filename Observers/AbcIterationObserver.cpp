//
// Created by Nikita Kruk on 26.11.17.
//

#include "AbcIterationObserver.hpp"

#include <sstream>
#include <cassert>
#include <iostream>
#include <algorithm> // std::copy
#include <iterator>  // std::ostream_iterator

AbcIterationObserver::AbcIterationObserver(AbstractParameterSet *parameter_set,
                                           const std::vector<Real> &velocity_distribution,
                                           int abc_iteration) :
    SyntheticObserver(parameter_set, velocity_distribution, 0),
    abc_iteration_(abc_iteration)
{
  data_file_name_ = ConstructDataFileName();
  std::remove(data_file_name_.c_str());
  data_file_.open(data_file_name_, std::ios::binary | std::ios::out | std::ios::app);
  assert(data_file_.is_open());

  std::string parameters_file_name = ConstructParameterFileName();
  SaveAllParametersIntoFile(parameters_file_name);
}

void AbcIterationObserver::operator()(std::vector<Real> &system_state, Real t)
{
  if (!(ouput_time_counter_ % output_time_threshold_))
  {
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

std::string AbcIterationObserver::ConstructDataFileName()
{
  std::ostringstream data_file_name_buffer;
  data_file_name_buffer << parameter_set_->GetAbcFolderName() << parameter_set_->GetCandidateDatasetsSubfolderName()
                        << "candidate_dataset.bin";
  return data_file_name_buffer.str();
}

std::string AbcIterationObserver::ConstructParameterFileName()
{
  std::ostringstream parameters_file_name_buffer;
  parameters_file_name_buffer << parameter_set_->GetAbcFolderName() << parameter_set_->GetCandidateDatasetsSubfolderName()
                              << "candidate_dataset_parameters.txt";
  return parameters_file_name_buffer.str();
}

std::string AbcIterationObserver::ConstructActivePassiveFileName()
{
  std::ostringstream active_passive_file_name_buffer;
  active_passive_file_name_buffer << parameter_set_->GetAbcFolderName()
                                  << parameter_set_->GetCandidateDatasetsSubfolderName()
                                  << "candidate_dataset_active_passive.txt";
  return active_passive_file_name_buffer.str();
}

std::string AbcIterationObserver::GetDataFileName()
{
  return data_file_name_;
}
