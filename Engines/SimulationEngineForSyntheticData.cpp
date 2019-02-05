//
// Created by Nikita Kruk on 25.11.17.
//

#include "SimulationEngineForSyntheticData.hpp"
#include "../DynamicalSystems/SprSystemWithUnitVelocity.hpp"
#include "../DynamicalSystems/SprSystemWithArbitraryStationaryVelocity.hpp"
#include "../Observers/SyntheticObserver.hpp"
#include "../Steppers/StochasticEulerStepperWithUnitVelocity.hpp"
#include "../Steppers/StochasticEulerStepperWithArbitraryStationaryVelocity.hpp"
#include "../Observers/AbcIterationObserver.hpp"

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm> // std::for_each

#include <boost/numeric/odeint.hpp>

SimulationEngineForSyntheticData::SimulationEngineForSyntheticData(SyntheticParameterSet &parameter_set,
                                   PeriodicBoundaryConditionsConfiguration &pbc_config) :
    parameter_set_(parameter_set),
    pbc_config_(pbc_config)
{
  system_state_ = std::vector<Real>(kS * parameter_set_.Get_N(), 0.0);
}

SimulationEngineForSyntheticData::~SimulationEngineForSyntheticData()
{

}

void SimulationEngineForSyntheticData::InitializeRandomSystemState()
{
  velocity_distribution_ = std::vector<Real>(parameter_set_.Get_N(), 1.0);

  std::mt19937 mersenne_twister_generator(std::random_device{}());
  std::uniform_real_distribution<Real> uni_real_dist(0, parameter_set_.Get_L());
  std::bernoulli_distribution bern_dist(0.5);

  bool placement_failed = false;
  Vector2D r_alpha, u_alpha;
  Vector2D r_beta, u_beta;
  Vector2D dr;

  std::cout << "rod initialization: ";
  for (int alpha = 0; alpha < parameter_set_.Get_N(); ++alpha)
  {
    do
    {
      placement_failed = false;

      r_alpha.x = uni_real_dist(mersenne_twister_generator);
      r_alpha.y = uni_real_dist(mersenne_twister_generator);

      for (int beta = 0; beta < alpha; ++beta)
      {
        r_beta.x = system_state_[kS * beta];
        r_beta.y = system_state_[kS * beta + 1];

        pbc_config_.ClassBEffectiveParticleDistance_unsigned(r_alpha, r_beta, dr);
        if ((fabs(dr.y) < parameter_set_.Get_l()) && (std::fabs(dr.x) < parameter_set_.Get_lambda() / 2.0))
        {
          placement_failed = true;
          break;
        }
      }
    } while (placement_failed);

    system_state_[kS * alpha] = r_alpha.x;
    system_state_[kS * alpha + 1] = r_alpha.y;
    system_state_[kS * alpha + 2] = 0.0;
    system_state_[kS * alpha + 3] = Real(bern_dist(mersenne_twister_generator)) * 2.0 - 1.0;

    std::cout << alpha << " ";
  }
  std::cout << std::endl;
  std::cout << "Rod Initialization Complete" << std::endl;
}

void SimulationEngineForSyntheticData::InitializeSystemStateFromSyntheticData(const std::string &synthetic_data_file_name)
{
  velocity_distribution_ = std::vector<Real>(parameter_set_.Get_N(), 0.0);

//  std::mt19937 mersenne_twister_generator(std::random_device{}());
//  std::normal_distribution<Real> norm_dist(0.0, parameter_set_.Get_lambda() / 3.0);

  std::ifstream synthetic_data_file(synthetic_data_file_name, std::ios::binary | std::ios::in);
  assert(synthetic_data_file.is_open());
  synthetic_data_file.seekg(parameter_set_.Get_t_0() * (1 + kS * parameter_set_.Get_N()) * sizeof(Real), std::ios::beg);
  Real time = 0.0;
  synthetic_data_file.read((char *) &time, sizeof(Real));
  synthetic_data_file.read((char *) &system_state_[0], kS * parameter_set_.Get_N() * sizeof(Real));
  synthetic_data_file.close();

  for (int alpha = 0; alpha < parameter_set_.Get_N(); ++alpha)
  {
    velocity_distribution_[alpha] = std::sqrt(system_state_[kS * alpha + 2] * system_state_[kS * alpha + 2]
                                                  + system_state_[kS * alpha + 3] * system_state_[kS * alpha + 3]);
  }
  // exert observational noise
/*  std::for_each(system_state_.begin(),
                system_state_.end(),
                [&](Real &s)
                {
                  s += norm_dist(mersenne_twister_generator);
                });
                */
}

void SimulationEngineForSyntheticData::InitializeRandomSystemStateWithVelocityDistribution()
{
  std::mt19937 mersenne_twister_generator(std::random_device{}());
  std::uniform_real_distribution<Real> uni_real_dist(0, parameter_set_.Get_L());
  std::bernoulli_distribution bern_dist(0.5);

  bool placement_failed = false;
  Vector2D r_alpha, u_alpha;
  Vector2D r_beta, u_beta;
  Vector2D dr;

  std::cout << "rod initialization: ";
  for (int alpha = 0; alpha < parameter_set_.Get_N(); ++alpha)
  {
    do
    {
      placement_failed = false;

      r_alpha.x = uni_real_dist(mersenne_twister_generator);
      r_alpha.y = uni_real_dist(mersenne_twister_generator);

      for (int beta = 0; beta < alpha; ++beta)
      {
        r_beta.x = system_state_[kS * beta];
        r_beta.y = system_state_[kS * beta + 1];

        pbc_config_.ClassBEffectiveParticleDistance_unsigned(r_alpha, r_beta, dr);
        if ((fabs(dr.y) < parameter_set_.Get_l()) && (std::fabs(dr.x) < parameter_set_.Get_lambda() / 2.0))
        {
          placement_failed = true;
          break;
        }
      }
    } while (placement_failed);

    system_state_[kS * alpha] = r_alpha.x;
    system_state_[kS * alpha + 1] = r_alpha.y;
    system_state_[kS * alpha + 2] = 0.0;
    system_state_[kS * alpha + 3] =
        velocity_distribution_[alpha] * (Real(bern_dist(mersenne_twister_generator)) * 2.0 - 1.0);

    std::cout << alpha << " ";
  }
  std::cout << std::endl;
  std::cout << "Rod Initialization Complete" << std::endl;
}

void SimulationEngineForSyntheticData::DrawVelocitiesFromPriorDistribution()
{
  velocity_distribution_ = std::vector<Real>(parameter_set_.Get_N(), 0.0);
  std::mt19937 mersenne_twister_generator(std::random_device{}());
  std::uniform_real_distribution<Real> uniform_distribution(0.0, 1.0);
  for (int alpha = 0; alpha < parameter_set_.Get_N(); ++alpha)
  {
    Real velocity = 0.0;
    do
    {
      velocity = DrawRayleighSample(uniform_distribution, mersenne_twister_generator, 300.0 * 0.000001);
    } while (velocity <= 0.0);
    velocity_distribution_[alpha] = velocity;
  }
}

Real SimulationEngineForSyntheticData::DrawRayleighSample(std::uniform_real_distribution<Real> &uniform_distribution,
                                          std::mt19937 &mersenne_twister_generator,
                                          Real b)
{
  return std::sqrt(-2.0 * b * b * std::log(uniform_distribution(mersenne_twister_generator)));
}

void SimulationEngineForSyntheticData::GenerateSyntheticDataWithUnitVelocity()
{
  InitializeRandomSystemState();

  SprSystemWithUnitVelocity system(&parameter_set_, &pbc_config_);
  SyntheticObserver observer(&parameter_set_, velocity_distribution_);
//	typedef boost::numeric::odeint::euler<state_type, Real, state_type, Real> euler_stepper_type;

  std::cout << "Synthetic Data Integration Started" << std::endl;
  // integrate_const calls first observer then stepper!
  boost::numeric::odeint::integrate_const(StochasticEulerStepperWithUnitVelocity(pbc_config_, parameter_set_),
                                          boost::ref(system),
                                          system_state_,
                                          parameter_set_.Get_t_0(),
                                          parameter_set_.Get_t_1(),
                                          parameter_set_.Get_delta_t(),
                                          boost::ref(observer));
  std::cout << "Synthetic Data Integration Complete" << std::endl;
}

void SimulationEngineForSyntheticData::GenerateSyntheticDataForDifferentParameters(Thread *thread)
{
  for (int i_a = 1; i_a <= 15; ++i_a)
  {
    parameter_set_.Set_a(i_a);
    for (int i_phi = 1; i_phi <= 9; ++i_phi)
    {
      parameter_set_.Set_phi(i_phi * 0.1);
      pbc_config_.SetNewBoundaries(parameter_set_.Get_L(), parameter_set_.Get_L());
      std::cout << "Synthetic Data Integration Started: " << "a=" << i_a << ", phi=" << i_phi * 0.1 << std::endl;

      InitializeRandomSystemState();
      SprSystemWithUnitVelocity system(&parameter_set_, &pbc_config_);
      // run through burn-in time
      // integrate_const calls first observer then stepper!
      boost::numeric::odeint::integrate_const(StochasticEulerStepperWithUnitVelocity(pbc_config_, parameter_set_),
                                              boost::ref(system),
                                              system_state_,
                                              parameter_set_.Get_t_0(),
                                              parameter_set_.Get_t_1() * 0.1,
                                              parameter_set_.Get_delta_t());
      // integrate_const calls first observer then stepper!
      SyntheticObserver observer(&parameter_set_, velocity_distribution_);
      boost::numeric::odeint::integrate_const(StochasticEulerStepperWithUnitVelocity(pbc_config_, parameter_set_),
                                              boost::ref(system),
                                              system_state_,
                                              parameter_set_.Get_t_0(),
                                              parameter_set_.Get_t_1(),
                                              parameter_set_.Get_delta_t(),
                                              boost::ref(observer));

      std::cout << "Synthetic Data Integration Complete: " << "a=" << i_a << ", phi=" << i_phi * 0.1 << std::endl;
    } // i_phi
  } // i_a
}

void SimulationEngineForSyntheticData::GenerateSyntheticDataWithArbitraryStationaryVelocity()
{
  DrawVelocitiesFromPriorDistribution();
  InitializeRandomSystemStateWithVelocityDistribution();

  SprSystemWithArbitraryStationaryVelocity system(&parameter_set_, &pbc_config_, velocity_distribution_);
  SyntheticObserver observer(&parameter_set_, velocity_distribution_);
  StochasticEulerStepperWithArbitraryStationaryVelocity stepper(&pbc_config_, &parameter_set_, velocity_distribution_);
  std::cout << "Synthetic Data Integration Started" << std::endl;
  Real t = 0.0;
  while (t <= parameter_set_.Get_burn_in_time())
  {
    t += parameter_set_.Get_delta_t();
    stepper.do_step(&system, system_state_, t, parameter_set_.Get_delta_t());
  }

  t = parameter_set_.Get_t_0();
  observer(system_state_, t);
  while (t <= parameter_set_.Get_t_1())
  {
    t += parameter_set_.Get_delta_t();
    stepper.do_step(&system, system_state_, t, parameter_set_.Get_delta_t());
    observer(system_state_, t);
  }
//  boost::numeric::odeint::integrate_const(StochasticEulerStepperWithArbitraryStationaryVelocity(&pbc_config_,
//                                                                                                &parameter_set_,
//                                                                                                velocity_distribution_),
//                                          &system,
//                                          system_state_,
//                                          0.0,
//                                          parameter_set_.Get_burn_in_time(),
//                                          parameter_set_.Get_delta_t());
//  boost::numeric::odeint::integrate_const(StochasticEulerStepperWithArbitraryStationaryVelocity(&pbc_config_,
//                                                                                                &parameter_set_,
//                                                                                                velocity_distribution_),
//                                          &system,
//                                          system_state_,
//                                          parameter_set_.Get_t_0(),
//                                          parameter_set_.Get_t_1(),
//                                          parameter_set_.Get_delta_t(),
//                                          boost::ref(observer));
  std::cout << "Synthetic Data Integration Complete" << std::endl;
}

void SimulationEngineForSyntheticData::SimulateCandidateDataset(int population_idx,
                                                const std::string &synthetic_data_file_name,
                                                std::string &simulation_data_file_name)
{
  InitializeSystemStateFromSyntheticData(synthetic_data_file_name);
  SprSystemWithArbitraryStationaryVelocity system(&parameter_set_, &pbc_config_, velocity_distribution_);
  AbcIterationObserver observer(&parameter_set_, velocity_distribution_, population_idx);
  StochasticEulerStepperWithArbitraryStationaryVelocity stepper(&pbc_config_, &parameter_set_, velocity_distribution_);

//  std::cout << "ABC Iteration:" << population_idx << " | Simulation Started" << std::endl;
  Real t = parameter_set_.Get_t_0();
  observer(system_state_, t);
  while (t <= parameter_set_.Get_t_1())
  {
    t += parameter_set_.Get_delta_t();
    stepper.do_step(&system, system_state_, t, parameter_set_.Get_delta_t());
    observer(system_state_, t);
  }
//  boost::numeric::odeint::integrate_const(StochasticEulerStepperWithArbitraryStationaryVelocity(pbc_config_,
//                                                                                                parameter_set_,
//                                                                                                velocity_distribution_),
//                                          &system,
//                                          system_state_,
//                                          parameter_set_.Get_t_0(),
//                                          parameter_set_.Get_t_1(),
//                                          parameter_set_.Get_delta_t(),
//                                          boost::ref(observer));
//  std::cout << "ABC Iteration:" << population_idx << " | Simulation Complete" << std::endl;
  simulation_data_file_name = observer.GetDataFileName();
}

const std::vector<Real> &SimulationEngineForSyntheticData::GetVelocityDistribution()
{
  return velocity_distribution_;
}