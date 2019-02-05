//
// Created by Nikita Kruk on 25.11.17.
//

#include "TwoParticleInteractionForce.hpp"

#include <cmath>

TwoParticleInteractionForce::TwoParticleInteractionForce(AbstractParameterSet *parameter_set) :
    parameter_set_(parameter_set)
{

}

TwoParticleInteractionForce::~TwoParticleInteractionForce()
{

}

Vector2D TwoParticleInteractionForce::Yukawa(const Vector2D &vec_r, Real r)
{
  Real lambda = parameter_set_->Get_lambda();
  return (lambda + r) / (lambda * r * r * r) * std::exp(-r / lambda) * vec_r;
}

Vector2D TwoParticleInteractionForce::LennardJones(const Vector2D &vec_r, Real r)
{
  Real lambda = parameter_set_->Get_lambda();
  Real tmp = std::pow(lambda / r, 6.0);
  return 12.0 * (tmp * tmp - tmp) / (r * r) * vec_r;
}

Vector2D TwoParticleInteractionForce::QuasiMorse(const Vector2D &vec_r, Real r)
{
  return Vector2D();
}
