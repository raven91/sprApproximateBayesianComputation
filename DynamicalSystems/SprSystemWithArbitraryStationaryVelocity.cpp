//
// Created by Nikita Kruk on 02.08.18.
//

#include "SprSystemWithArbitraryStationaryVelocity.hpp"

#include <iostream>

SprSystemWithArbitraryStationaryVelocity::SprSystemWithArbitraryStationaryVelocity(AbstractParameterSet *parameter_set,
                                                                                   BoundaryConditionsConfiguration *bc_config,
                                                                                   std::vector<Real> &velocity_distribution)
    : SprSystemWithUnitVelocity(parameter_set, bc_config),
      velocity_distribution_(velocity_distribution)
{

}

SprSystemWithArbitraryStationaryVelocity::~SprSystemWithArbitraryStationaryVelocity()
{

}

void SprSystemWithArbitraryStationaryVelocity::CalculateDerivatives(const std::vector<Real> &system_state,
                                                                    std::vector<Real> &derivative)
{
  int n = parameter_set_->Get_n();
  const Real k_B = 1.38064852 * std::pow(10.0, -23.0);

  Vector2D u_alpha_normalized;
  Matrix2D diffusion_tensor, identity(1.0, 0.0, 0.0, 1.0);
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    u_alpha_normalized.x = system_state[kS * alpha + 2];
    u_alpha_normalized.y = system_state[kS * alpha + 3];
    diffusion_tensor = OuterProduct(u_alpha_normalized, u_alpha_normalized);
    diffusion_tensor = parameter_set_->Get_D_perp() * (identity - diffusion_tensor)
        + parameter_set_->Get_D_parallel() * diffusion_tensor;

    position_interaction_force_[alpha] *= (parameter_set_->Get_U_0() / (n * n));
    position_interaction_force_[alpha] = velocity_distribution_[alpha] * u_alpha_normalized
        - 1.0 / (k_B * parameter_set_->Get_T()) * diffusion_tensor * position_interaction_force_[alpha];
    velocity_interaction_force_[alpha] *=
        (-1.0 / (k_B * parameter_set_->Get_T()) * parameter_set_->Get_D_R() * parameter_set_->Get_U_0() / (n * n));

    derivative[kS * alpha] = position_interaction_force_[alpha].x;
    derivative[kS * alpha + 1] = position_interaction_force_[alpha].y;
    derivative[kS * alpha + 2] = velocity_interaction_force_[alpha].x;
    derivative[kS * alpha + 3] = velocity_interaction_force_[alpha].y;
  }
}

void SprSystemWithArbitraryStationaryVelocity::OperatorNoise(std::vector<Real> &system_state,
                                                             std::vector<Real> &derivative,
                                                             Real t)
{
  std::normal_distribution<Real> standard_normal_distribution(0.0, 1.0);
  Matrix2D rotate_90_cw;
  rotate_90_cw.xx = 0.0;
  rotate_90_cw.xy = 1.0;
  rotate_90_cw.yx = -1.0;
  rotate_90_cw.yy = 0.0;
  Vector2D u_alpha, u_perp_alpha;
  Real noise_1 = 0.0, noise_2 = 0.0, noise_3 = 0.0;

  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    if (parameter_set_->Get_active_passive().at(alpha))
    {
      u_alpha.x = system_state[kS * alpha + 2];
      u_alpha.y = system_state[kS * alpha + 3];
      u_perp_alpha = rotate_90_cw * u_alpha;

      noise_1 = standard_normal_distribution(mersenne_twister_generator_);
      noise_2 = standard_normal_distribution(mersenne_twister_generator_);
      noise_3 = standard_normal_distribution(mersenne_twister_generator_);

      derivative[kS * alpha] = std::sqrt(2.0 * parameter_set_->Get_D_parallel()) * noise_1 * u_alpha.x
          + std::sqrt(2.0 * parameter_set_->Get_D_perp()) * noise_2 * u_perp_alpha.x;
      derivative[kS * alpha + 1] = std::sqrt(2.0 * parameter_set_->Get_D_parallel()) * noise_1 * u_alpha.y
          + std::sqrt(2.0 * parameter_set_->Get_D_perp()) * noise_2 * u_perp_alpha.y;
      derivative[kS * alpha + 2] = std::sqrt(2.0 * parameter_set_->Get_D_R()) * noise_3 * u_perp_alpha.x;
      derivative[kS * alpha + 3] = std::sqrt(2.0 * parameter_set_->Get_D_R()) * noise_3 * u_perp_alpha.y;
    }
  }
}