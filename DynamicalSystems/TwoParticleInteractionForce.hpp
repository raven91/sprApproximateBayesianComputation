//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_TWOPARTICLEINTERACTIONFORCE_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_TWOPARTICLEINTERACTIONFORCE_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/AbstractParameterSet.hpp"
#include "../Other/Vector2D.hpp"

class TwoParticleInteractionForce
{
 public:

  explicit TwoParticleInteractionForce(AbstractParameterSet *parameter_set);
  ~TwoParticleInteractionForce();

  // \nabla_{\vec{r}^\alpha} U^{\alpha\beta}_{ij} = -f^{\alpha\beta}_{ij}
  Vector2D Yukawa(const Vector2D &vec_r, Real r);
  Vector2D LennardJones(const Vector2D &vec_r, Real r);
  Vector2D QuasiMorse(const Vector2D &vec_r, Real r);

 private:

  AbstractParameterSet *parameter_set_;
};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_TWOPARTICLEINTERACTIONFORCE_HPP
