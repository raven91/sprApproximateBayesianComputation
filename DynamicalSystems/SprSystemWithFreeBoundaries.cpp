//
// Created by Nikita Kruk on 08.08.18.
//

#include "SprSystemWithFreeBoundaries.hpp"

SprSystemWithFreeBoundaries::SprSystemWithFreeBoundaries(AbstractParameterSet *parameter_set,
                                                         BoundaryConditionsConfiguration *bc_config,
                                                         std::vector<Real> &velocity_distribution)
    : SprSystemWithArbitraryStationaryVelocity(parameter_set, bc_config, velocity_distribution)
{
  neighboring_cells_ =
      {
          {-1, -1}, {0, -1}, {1, -1},
          {-1, 0}, {0, 0}, {1, 0},
          {-1, 1}, {0, 1}, {1, 1}
      };
}

SprSystemWithFreeBoundaries::~SprSystemWithFreeBoundaries()
{

}

void SprSystemWithFreeBoundaries::CalculateMaxDisplacement(const std::vector<Real> &system_state)
{
  Real L = parameter_set_->Get_L();
  Real max = 0.0;
  Real dist = 0.0;
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    if (system_state[kS * alpha] >= 0.0 && system_state[kS * alpha] < L
        && system_state[kS * alpha + 1] >= 0.0 && system_state[kS * alpha + 1] < L)
    {
      dist = Norm(position_interaction_force_[alpha]) * parameter_set_->Get_delta_t();
      if (max < dist)
      {
        max = dist;
      }
    }
  }
  accumulated_displacement_ += max;
  if (accumulated_displacement_ > 0.5 * (r_max_ - r_min_))
  {
    accumulated_displacement_ = 0.0;
    should_update_lists_ = true;
  }
}

void SprSystemWithFreeBoundaries::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y)
{
  // No adjustment needed for free boundary conditions
}

bool SprSystemWithFreeBoundaries::ShouldConsiderFreeBoundaries()
{
  return true;
}

bool SprSystemWithFreeBoundaries::ShouldConsiderHalfOfNeighbors()
{
  return false;
}