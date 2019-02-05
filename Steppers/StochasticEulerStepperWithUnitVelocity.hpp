//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_STOCHASTICEULERSTEPPER_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_STOCHASTICEULERSTEPPER_HPP

#include "../Definitions.hpp"
#include "../Other/PeriodicBoundaryConditionsConfiguration.hpp"
#include "../ParameterSets/SyntheticParameterSet.hpp"

#include <cmath>

#include <boost/numeric/odeint.hpp>

class StochasticEulerStepperWithUnitVelocity
{
 public:

  typedef std::vector<Real> state_type;
  typedef std::vector<Real> deriv_type;
  typedef Real value_type;
  typedef Real time_type;
  typedef unsigned short order_type;

  typedef boost::numeric::odeint::stepper_tag stepper_category;

  explicit StochasticEulerStepperWithUnitVelocity(PeriodicBoundaryConditionsConfiguration &pbc_config,
                                  SyntheticParameterSet &parameter_set) :
      pbc_config_(pbc_config),
      parameter_set_(parameter_set)
  {}
  ~StochasticEulerStepperWithUnitVelocity()
  {}

  static order_type order()
  { return 1; }

  template<class System>
  void do_step(System system, state_type &x, time_type t, time_type dt) const
  {
    deriv_type det(x.size(), 0.0), stoch(x.size(), 0.0);

    system.get().OperatorVerletNeighborList(x, det, t);
    system.get().OperatorNoise(x, stoch, t);

    for (size_t i = 0; i < x.size(); ++i)
    {
      x[i] += (dt * det[i] + sqrt(dt) * stoch[i]);
    }

    // keep all positions under periodic boundary conditions
    pbc_config_.ApplyPeriodicBoundaryConditions(x, parameter_set_.Get_N());
    // keep all velocities unitary
    for (int alpha = 0; alpha < parameter_set_.Get_N(); ++alpha)
    {
      Normalize(x[kS * alpha + 2], x[kS * alpha + 3]);
    } // TODO: circumvent the artificial velocity normalization
  }

 private:

  PeriodicBoundaryConditionsConfiguration &pbc_config_;
  SyntheticParameterSet &parameter_set_;

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_STOCHASTICEULERSTEPPER_HPP
