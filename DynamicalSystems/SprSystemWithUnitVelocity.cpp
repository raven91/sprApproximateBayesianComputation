//
// Created by Nikita Kruk on 25.11.17.
//

#include "SprSystemWithUnitVelocity.hpp"

#include <iostream>

SprSystemWithUnitVelocity::SprSystemWithUnitVelocity(AbstractParameterSet *parameter_set,
                                                     BoundaryConditionsConfiguration *bc_config) :
    parameter_set_(parameter_set),
    bc_config_(bc_config),
    two_particle_interaction_force_(parameter_set),
    position_interaction_force_(std::vector<Vector2D>(parameter_set->Get_N())),
    velocity_interaction_force_(std::vector<Vector2D>(parameter_set->Get_N())),
    mersenne_twister_generator_(std::random_device{}())
{
  //linked list and cell division
  x_size_ = parameter_set_->Get_L();
  y_size_ = parameter_set_->Get_L();

  num_subcells_x_ = int(parameter_set_->Get_L() / (2.0 * parameter_set_->Get_lambda()));
  num_subcells_y_ = num_subcells_x_;

  linked_list_ = std::vector<std::vector<SegmentIndex>>(num_subcells_x_ * num_subcells_y_, std::vector<SegmentIndex>());

  if (ShouldConsiderHalfOfNeighbors())
  {
    //half of neighbors
    neighboring_cells_ =
        {
            {0, 0}, {1, 0},
            {-1, 1}, {0, 1}, {1, 1}
        };
  } else
  {
    //all neighbors
    neighboring_cells_ =
        {
            {-1, -1}, {0, -1}, {1, -1},
            {-1, 0}, {0, 0}, {1, 0},
            {-1, 1}, {0, 1}, {1, 1}
        };
  }

  //Verlet neighbor list
  r_min_ = parameter_set_->Get_lambda() * 1.5;
  r_max_ = r_min_ + parameter_set_->Get_lambda();
  should_update_lists_ = true;
  accumulated_displacement_ = 0.0;

  verlet_list_ = std::vector<std::vector<SegmentIndex>>(parameter_set_->Get_N() * parameter_set_->Get_n(),
                                                        std::vector<SegmentIndex>());
}

SprSystemWithUnitVelocity::~SprSystemWithUnitVelocity()
{

}

void SprSystemWithUnitVelocity::OperatorLinkedList(std::vector<Real> &system_state,
                                                   std::vector<Real> &derivative,
                                                   Real t)
{
  ConstructNewLinkedList(system_state);
  CalculateInteractionForcesUsingLinkedList(system_state);
  CalculateDerivatives(system_state, derivative);

//	Real max = 0.0;
//	Real dist = 0.0;
//	for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
//	{
//		dist = Norm(position_interaction_force_[alpha]) * parameter_set_->Get_delta_t();
//		if (max < dist)
//		{
//			max = dist;
//		}
//	}
//	std::cout << "dist: " << max << std::endl;
}

void SprSystemWithUnitVelocity::OperatorVerletNeighborList(std::vector<Real> &system_state,
                                                           std::vector<Real> &derivative,
                                                           Real t)
{
  //if it is time to update the linked list and the Verlet list
  if (should_update_lists_)
  {
    ConstructNewLinkedList(system_state);
    ConstructNewVerletList(system_state);
    //reset the counter of the passed time steps
    should_update_lists_ = false;
  }
  CalculateInteractionForcesUsingVerletList(system_state);
  CalculateDerivatives(system_state, derivative);

  CalculateMaxDisplacement(system_state);
}

void SprSystemWithUnitVelocity::CalculateMaxDisplacement(const std::vector<Real> &system_state)
{
  Real max = 0.0;
  Real dist = 0.0;
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    dist = Norm(position_interaction_force_[alpha]) * parameter_set_->Get_delta_t();
    if (max < dist)
    {
      max = dist;
    }
  }
  accumulated_displacement_ += max;
  if (accumulated_displacement_ > 0.5 * (r_max_ - r_min_))
  {
    accumulated_displacement_ = 0.0;
    should_update_lists_ = true;
  }
}

void SprSystemWithUnitVelocity::ConstructNewLinkedList(const std::vector<Real> &system_state)
{
  int n = parameter_set_->Get_n();
  Real L = parameter_set_->Get_L();
  Real l = parameter_set_->Get_l();
  Real lambda = parameter_set_->Get_lambda();

#pragma unroll
  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  }

  Vector2D r_alpha, u_alpha, r_i_alpha;
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  Real l_i = 0.0;
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    r_alpha.x = system_state[kS * alpha];
    r_alpha.y = system_state[kS * alpha + 1];
    u_alpha.x = system_state[kS * alpha + 2];
    u_alpha.y = system_state[kS * alpha + 3];

    for (int i = 0; i < n; ++i)
    {
      if (1 == n)
      {
        l_i = 0.0;
      } else
      {
        l_i = (l - lambda) * (-0.5 + Real(i) / Real(n - 1));
      }
      r_i_alpha = r_alpha + l_i * u_alpha;
      if (r_i_alpha.x < 0.0 || r_i_alpha.x >= L || r_i_alpha.y < 0.0 || r_i_alpha.y >= L)
      {
        bc_config_->ApplyPeriodicBoundaryConditions(r_i_alpha, r_i_alpha);
        if (ShouldConsiderFreeBoundaries())
        {
          continue; // Ignore bacteria outside the ROI with free boundaries
        }
      }

      i_cell_x = int(r_i_alpha.x / x_size_ * num_subcells_x_);
      i_cell_y = int(r_i_alpha.y / y_size_ * num_subcells_y_);
      i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

      linked_list_[i_cell].push_back(SegmentIndex(alpha, i));
    }
  }
}

void SprSystemWithUnitVelocity::ConstructNewVerletList(const std::vector<Real> &system_state)
{
  int n = parameter_set_->Get_n();
  Real L = parameter_set_->Get_L();
  Real l = parameter_set_->Get_l();
  Real lambda = parameter_set_->Get_lambda();

#pragma unroll
  for (int i = 0; i < verlet_list_.size(); ++i)
  {
    verlet_list_[i].clear();
  }

  int beta = 0, j = 0;
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  Vector2D r_alpha, u_alpha, r_i_alpha;
  Vector2D r_beta, u_beta, r_j_beta;
  Real l_i = 0.0, l_j = 0.0;
  Vector2D dr_ij_ab;
  Real r_ij_ab = 0.0, r_ij_ab_squared = 0.0;
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    r_alpha.x = system_state[kS * alpha];
    r_alpha.y = system_state[kS * alpha + 1];
    u_alpha.x = system_state[kS * alpha + 2];
    u_alpha.y = system_state[kS * alpha + 3];

    for (int i = 0; i < n; ++i)
    {
      if (1 == n)
      {
        l_i = 0.0;
      } else
      {
        l_i = (l - lambda) * (-0.5 + Real(i) / Real(n - 1));
      }
      r_i_alpha = r_alpha + l_i * u_alpha;
      if (r_i_alpha.x < 0.0 || r_i_alpha.x >= L || r_i_alpha.y < 0.0 || r_i_alpha.y >= L)
      {
        bc_config_->ApplyPeriodicBoundaryConditions(r_i_alpha, r_i_alpha);
        if (ShouldConsiderFreeBoundaries())
        {
          continue; // Ignore bacteria outside the ROI with free boundaries
        }
      }

      i_cell_x = int(r_i_alpha.x / x_size_ * num_subcells_x_);
      i_cell_y = int(r_i_alpha.y / y_size_ * num_subcells_y_);
      i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

      for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
      {
        j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
        j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
        if (j_cell_x < 0 || j_cell_x >= num_subcells_x_ || j_cell_y < 0 || j_cell_y >= num_subcells_y_)
        {
          if (ShouldConsiderFreeBoundaries())
          {
            continue; // Ignore bacteria outside the ROI with free boundaries
          }
        }
        AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
        j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

        for (std::vector<SegmentIndex>::iterator neighbr_iter = linked_list_[j_cell].begin();
             neighbr_iter != linked_list_[j_cell].end(); ++neighbr_iter)
        {
          beta = neighbr_iter->alpha;
          j = neighbr_iter->i;
          if ((i_cell != j_cell && beta != alpha) || (i_cell == j_cell &&
              ((beta < alpha && ShouldConsiderHalfOfNeighbors())
                  || (beta != alpha && !ShouldConsiderHalfOfNeighbors()))))
          {
            r_beta.x = system_state[kS * beta];
            r_beta.y = system_state[kS * beta + 1];
            u_beta.x = system_state[kS * beta + 2];
            u_beta.y = system_state[kS * beta + 3];

            if (1 == n)
            {
              l_j = 0.0;
            } else
            {
              l_j = (l - lambda) * (-0.5 + Real(j) / Real(n - 1));
            }
            r_j_beta = r_beta + l_j * u_beta;
            if (r_j_beta.x < 0.0 || r_j_beta.x >= L || r_j_beta.y < 0.0 || r_j_beta.y >= L)
            {
              bc_config_->ApplyPeriodicBoundaryConditions(r_j_beta, r_j_beta);
              if (ShouldConsiderFreeBoundaries())
              {
                continue; // Ignore bacteria outside the ROI with free boundaries
              }
            } // TODO: is this condition necessary if class B distance is computed?

            bc_config_->ClassBEffectiveParticleDistance_signed(r_i_alpha, r_j_beta, dr_ij_ab);
            r_ij_ab_squared = NormSquared(dr_ij_ab);

            if (r_ij_ab_squared <= r_max_ * r_max_)
            {
              verlet_list_[alpha * n + i].push_back(SegmentIndex(beta, j));
            }
          }
        } // neighbr_iter
      } // neighboring_cell
    } // i
  } // alpha
}

void SprSystemWithUnitVelocity::CalculateInteractionForcesUsingLinkedList(const std::vector<Real> &system_state)
{
  int n = parameter_set_->Get_n();
  Real L = parameter_set_->Get_L();
  Real l = parameter_set_->Get_l();
  Real lambda = parameter_set_->Get_lambda();

  int beta = 0, j = 0;
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  Vector2D r_alpha, u_alpha, r_i_alpha;
  Vector2D r_beta, u_beta, r_j_beta;
  Real l_i = 0.0, l_j = 0.0;
  Vector2D dr_ij_ab, df_ij_ab;
  Real r_ij_ab = 0.0;

  //loop using the linked list
  std::fill(position_interaction_force_.begin(), position_interaction_force_.end(), Vector2D());
  std::fill(velocity_interaction_force_.begin(), velocity_interaction_force_.end(), Vector2D());
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    r_alpha.x = system_state[kS * alpha];
    r_alpha.y = system_state[kS * alpha + 1];
    u_alpha.x = system_state[kS * alpha + 2];
    u_alpha.y = system_state[kS * alpha + 3];

    for (int i = 0; i < n; ++i)
    {
      if (1 == n)
      {
        l_i = 0.0;
      } else
      {
        l_i = (l - lambda) * (-0.5 + Real(i) / Real(n - 1));
      }
      r_i_alpha = r_alpha + l_i * u_alpha;
      if (r_i_alpha.x < 0.0 || r_i_alpha.x >= L || r_i_alpha.y < 0.0 || r_i_alpha.y >= L)
      {
        bc_config_->ApplyPeriodicBoundaryConditions(r_i_alpha, r_i_alpha);
        if (ShouldConsiderFreeBoundaries())
        {
          continue; // Ignore bacteria outside the ROI with free boundaries
        }
      }

      i_cell_x = int(r_i_alpha.x / x_size_ * num_subcells_x_);
      i_cell_y = int(r_i_alpha.y / y_size_ * num_subcells_y_);
      i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

      for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
      {
        j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
        j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
        if (j_cell_x < 0 || j_cell_x >= num_subcells_x_ || j_cell_y < 0 || j_cell_y >= num_subcells_y_)
        {
          if (ShouldConsiderFreeBoundaries())
          {
            continue; // Ignore bacteria outside the ROI with free boundaries
          }
        }
        AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
        j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

        for (std::vector<SegmentIndex>::iterator neighbr_iter = linked_list_[j_cell].begin();
             neighbr_iter != linked_list_[j_cell].end(); ++neighbr_iter)
        {
          beta = neighbr_iter->alpha;
          j = neighbr_iter->i;
          if ((i_cell != j_cell && beta != alpha) || (i_cell == j_cell &&
              ((beta < alpha && ShouldConsiderHalfOfNeighbors())
                  || (beta != alpha && !ShouldConsiderHalfOfNeighbors()))))
          {
            r_beta.x = system_state[kS * beta];
            r_beta.y = system_state[kS * beta + 1];
            u_beta.x = system_state[kS * beta + 2];
            u_beta.y = system_state[kS * beta + 3];

            l_j = (l - lambda) * (-0.5 + Real(j) / Real(n - 1));
            r_j_beta = r_beta + l_j * u_beta;
            if (r_j_beta.x < 0.0 || r_j_beta.x >= L || r_j_beta.y < 0.0 || r_j_beta.y >= L)
            {
              bc_config_->ApplyPeriodicBoundaryConditions(r_j_beta, r_j_beta);
              if (ShouldConsiderFreeBoundaries())
              {
                continue; // Ignore bacteria outside the ROI with free boundaries
              }
            }

            bc_config_->ClassBEffectiveParticleDistance_signed(r_i_alpha, r_j_beta, dr_ij_ab);
            r_ij_ab = Norm(dr_ij_ab);

            if (r_ij_ab <= lambda)
            {
              df_ij_ab = two_particle_interaction_force_.Yukawa(dr_ij_ab, r_ij_ab);

              position_interaction_force_[alpha] += df_ij_ab;
              velocity_interaction_force_[alpha] += l_i * df_ij_ab;

              if (ShouldConsiderHalfOfNeighbors())
              {
                position_interaction_force_[beta] -= df_ij_ab;
                velocity_interaction_force_[beta] -= l_j * df_ij_ab;
              }
            }
          }
        } // neighbr_iter
      } // neighboring_cell
    } // i
  } // alpha
}

void SprSystemWithUnitVelocity::CalculateInteractionForcesUsingVerletList(const std::vector<Real> &system_state)
{
  int n = parameter_set_->Get_n();
  Real L = parameter_set_->Get_L();
  Real l = parameter_set_->Get_l();
  Real lambda = parameter_set_->Get_lambda();

  int beta = 0, j = 0;
  Vector2D r_alpha, u_alpha, r_i_alpha;
  Vector2D r_beta, u_beta, r_j_beta;
  Real l_i = 0.0, l_j = 0.0;
  Vector2D dr_ij_ab, df_ij_ab;
  Real r_ij_ab = 0.0;

  //loop using the verlet list
  std::fill(position_interaction_force_.begin(), position_interaction_force_.end(), Vector2D());
  std::fill(velocity_interaction_force_.begin(), velocity_interaction_force_.end(), Vector2D());
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    r_alpha.x = system_state[kS * alpha];
    r_alpha.y = system_state[kS * alpha + 1];
    u_alpha.x = system_state[kS * alpha + 2];
    u_alpha.y = system_state[kS * alpha + 3];

    for (int i = 0; i < n; ++i)
    {
      if (1 == n)
      {
        l_i = 0.0;
      } else
      {
        l_i = (l - lambda) * (-0.5 + Real(i) / Real(n - 1));
      }
      r_i_alpha = r_alpha + l_i * u_alpha;
      if (r_i_alpha.x < 0.0 || r_i_alpha.x >= L || r_i_alpha.y < 0.0 || r_i_alpha.y >= L)
      {
        bc_config_->ApplyPeriodicBoundaryConditions(r_i_alpha, r_i_alpha);
        if (ShouldConsiderFreeBoundaries())
        {
          continue; // Ignore bacteria outside the ROI with free boundaries
        }
      }

      for (std::vector<SegmentIndex>::iterator neighbr_iter = verlet_list_[alpha * n + i].begin();
           neighbr_iter != verlet_list_[alpha * n + i].end(); ++neighbr_iter)
      {
        beta = neighbr_iter->alpha;
        j = neighbr_iter->i;

        r_beta.x = system_state[kS * beta];
        r_beta.y = system_state[kS * beta + 1];
        u_beta.x = system_state[kS * beta + 2];
        u_beta.y = system_state[kS * beta + 3];

        if (1 == n)
        {
          l_j = 0.0;
        } else
        {
          l_j = (l - lambda) * (-0.5 + Real(j) / Real(n - 1));
        }
        r_j_beta = r_beta + l_j * u_beta;
        if (r_j_beta.x < 0.0 || r_j_beta.x >= L || r_j_beta.y < 0.0 || r_j_beta.y >= L)
        {
          bc_config_->ApplyPeriodicBoundaryConditions(r_j_beta, r_j_beta);
          if (ShouldConsiderFreeBoundaries())
          {
            continue; // Ignore bacteria outside the ROI with free boundaries
          }
        }

        bc_config_->ClassBEffectiveParticleDistance_signed(r_i_alpha, r_j_beta, dr_ij_ab);
        r_ij_ab = Norm(dr_ij_ab);

        if (r_ij_ab <= lambda)
        {
          df_ij_ab = two_particle_interaction_force_.Yukawa(dr_ij_ab, r_ij_ab);

          position_interaction_force_[alpha] += df_ij_ab;
          velocity_interaction_force_[alpha] += l_i * df_ij_ab;

          if (ShouldConsiderHalfOfNeighbors())
          {
            position_interaction_force_[beta] -= df_ij_ab;
            velocity_interaction_force_[beta] -= l_j * df_ij_ab;
          }
        }
      } // neighbr_iter
    } // i
  } // alpha
}

void SprSystemWithUnitVelocity::CalculateDerivatives(const std::vector<Real> &system_state,
                                                     std::vector<Real> &derivative)
{
  int n = parameter_set_->Get_n();

  Vector2D u_alpha;
  Matrix2D inverse_tensor;
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    u_alpha.x = system_state[kS * alpha + 2];
    u_alpha.y = system_state[kS * alpha + 3];

    inverse_tensor.xx =
        parameter_set_->Get_f_par() * u_alpha.y * u_alpha.y + parameter_set_->Get_f_orth() * u_alpha.x * u_alpha.x;
    inverse_tensor.xy = (parameter_set_->Get_f_orth() - parameter_set_->Get_f_par()) * u_alpha.x * u_alpha.y;
    inverse_tensor.yx = inverse_tensor.xy;
    inverse_tensor.yy =
        parameter_set_->Get_f_par() * u_alpha.x * u_alpha.x + parameter_set_->Get_f_orth() * u_alpha.y * u_alpha.y;
    inverse_tensor /= (parameter_set_->Get_f_0() * parameter_set_->Get_f_par() * parameter_set_->Get_f_orth());

    position_interaction_force_[alpha] *= (parameter_set_->Get_U_0() / (n * n));
    position_interaction_force_[alpha] = parameter_set_->Get_active_passive().at(alpha)
        * (parameter_set_->Get_F() / parameter_set_->Get_f_0() / parameter_set_->Get_f_par()) * u_alpha
        - inverse_tensor * position_interaction_force_[alpha];
    velocity_interaction_force_[alpha] *=
        (-1.0 / parameter_set_->Get_f_0() / parameter_set_->Get_f_R() * parameter_set_->Get_U_0() / (n * n));

    derivative[kS * alpha] = position_interaction_force_[alpha].x;
    derivative[kS * alpha + 1] = position_interaction_force_[alpha].y;
    derivative[kS * alpha + 2] = velocity_interaction_force_[alpha].x;
    derivative[kS * alpha + 3] = velocity_interaction_force_[alpha].y;
  }
}

void SprSystemWithUnitVelocity::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y)
{
  if (cell_x < 0)
  {
    cell_x = num_subcells_x_ - 1;
  } else if (cell_x >= num_subcells_x_)
  {
    cell_x = 0;
  }

  if (cell_y < 0)
  {
    cell_y = num_subcells_y_ - 1;
  } else if (cell_y >= num_subcells_y_)
  {
    cell_y = 0;
  }
}

//Calculate the stochastic part of SDE
void SprSystemWithUnitVelocity::OperatorNoise(std::vector<Real> &system_state, std::vector<Real> &derivative, Real t)
{
//  std::normal_distribution<Real> norm_dist(0.0, 2.0 / (parameter_set_->Get_f_0() * parameter_set_->Get_f_R()));
  std::normal_distribution<Real> norm_dist(0.0, 1.0);
  for (int alpha = 0; alpha < parameter_set_->Get_N(); ++alpha)
  {
    if (parameter_set_->Get_active_passive().at(alpha))
    {
      derivative[kS * alpha + 2] = parameter_set_->Get_D() * norm_dist(mersenne_twister_generator_);
      derivative[kS * alpha + 3] = parameter_set_->Get_D() * norm_dist(mersenne_twister_generator_);
    }
  }
}

bool SprSystemWithUnitVelocity::ShouldConsiderFreeBoundaries()
{
  return false;
}

bool SprSystemWithUnitVelocity::ShouldConsiderHalfOfNeighbors()
{
  return false;
}