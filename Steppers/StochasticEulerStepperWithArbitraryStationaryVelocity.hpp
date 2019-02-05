//
// Created by Nikita Kruk on 02.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_STOCHASTICEULERSTEPPERWITHARBITRARYSTATIONARYVELOCITY_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_STOCHASTICEULERSTEPPERWITHARBITRARYSTATIONARYVELOCITY_HPP

#include "../Definitions.hpp"
#include "../Other/BoundaryConditionsConfiguration.hpp"
#include "../ParameterSets/AbstractParameterSet.hpp"

#include <cmath>

#include <boost/numeric/odeint.hpp>

class StochasticEulerStepperWithArbitraryStationaryVelocity
{
 public:

  typedef std::vector<Real> state_type;
  typedef std::vector<Real> deriv_type;
  typedef Real value_type;
  typedef Real time_type;
  typedef unsigned short order_type;

  typedef boost::numeric::odeint::stepper_tag stepper_category;

  explicit StochasticEulerStepperWithArbitraryStationaryVelocity(BoundaryConditionsConfiguration *bc_config,
                                                                 AbstractParameterSet *parameter_set,
                                                                 std::vector<Real> &velocity_distribution) :
      bc_config_(bc_config),
      parameter_set_(parameter_set),
      velocity_distribution_(velocity_distribution)
  {}
  ~StochasticEulerStepperWithArbitraryStationaryVelocity()
  {}

  static order_type order()
  { return 1; }

  template<class System>
  void do_step(System *system, state_type &x, time_type t, time_type dt) const
  {
    deriv_type det(x.size(), 0.0), stoch(x.size(), 0.0);

    system->OperatorVerletNeighborList(x, det, t);
    system->OperatorNoise(x, stoch, t);

    for (size_t i = 0; i < x.size(); ++i)
    {
      x[i] += (dt * det[i] + sqrt(dt) * stoch[i]);
    }

    // keep all positions under periodic boundary conditions
    bc_config_->ApplyPeriodicBoundaryConditions(x, parameter_set_->Get_N());
    // keep all velocity vectors unitary
    for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
    {
      Normalize(x[kS * alpha + 2], x[kS * alpha + 3]);
    } // TODO: circumvent the artificial velocity normalization
  }

 private:

  BoundaryConditionsConfiguration *bc_config_;
  AbstractParameterSet *parameter_set_;
  std::vector<Real> &velocity_distribution_;

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_STOCHASTICEULERSTEPPERWITHARBITRARYSTATIONARYVELOCITY_HPP
