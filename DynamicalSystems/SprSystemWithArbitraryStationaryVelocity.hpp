//
// Created by Nikita Kruk on 02.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEMWITHARBITRARYSTATIONARYVELOCITY_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEMWITHARBITRARYSTATIONARYVELOCITY_HPP

#include "../Definitions.hpp"
#include "SprSystemWithUnitVelocity.hpp"
#include "../ParameterSets/AbstractParameterSet.hpp"
#include "../Other/BoundaryConditionsConfiguration.hpp"
#include "../Other/Vector2D.hpp"
#include "TwoParticleInteractionForce.hpp"

class SprSystemWithArbitraryStationaryVelocity : public SprSystemWithUnitVelocity
{
 public:

  explicit SprSystemWithArbitraryStationaryVelocity(AbstractParameterSet *parameter_set,
                                                    BoundaryConditionsConfiguration *bc_config,
                                                    std::vector<Real> &velocity_distribution);
  virtual ~SprSystemWithArbitraryStationaryVelocity();

  virtual void OperatorNoise(std::vector<Real> &system_state, std::vector<Real> &derivative, Real t);

 protected:

  std::vector<Real> &velocity_distribution_;

  virtual void CalculateDerivatives(const std::vector<Real> &system_state, std::vector<Real> &derivative);

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEMWITHARBITRARYSTATIONARYVELOCITY_HPP
