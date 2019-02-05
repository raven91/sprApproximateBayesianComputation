//
// Created by Nikita Kruk on 25.11.17.
//

#include "PeriodicBoundaryConditionsConfiguration.hpp"

#include <cmath>
#include <iostream>

PeriodicBoundaryConditionsConfiguration::PeriodicBoundaryConditionsConfiguration(Real x_size, Real y_size) :
    BoundaryConditionsConfiguration(x_size, y_size)
{

}

PeriodicBoundaryConditionsConfiguration::~PeriodicBoundaryConditionsConfiguration()
{

}

void PeriodicBoundaryConditionsConfiguration::ClassAEffectiveParticleDistance(Real x_i,
                                                                              Real y_i,
                                                                              Real x_j,
                                                                              Real y_j,
                                                                              Real &dx,
                                                                              Real &dy)
{
  dx = x_j - x_i;
  dx -= static_cast<int>(dx * 2.0 * x_rsize_) * x_size_;

  dy = y_j - y_i;
  dy -= static_cast<int>(dy * 2.0 * y_rsize_) * y_size_;
}

void PeriodicBoundaryConditionsConfiguration::ClassAEffectiveParticleDistance(const Vector2D &r_i,
                                                                              const Vector2D &r_j,
                                                                              Vector2D &dr_ij)
{
  dr_ij.x = r_j.x - r_i.x;
  dr_ij.x -= static_cast<int>(dr_ij.x * 2.0 / x_size_) * x_size_;

  dr_ij.y = r_j.y - r_i.y;
  dr_ij.y -= static_cast<int>(dr_ij.y * 2.0 / y_size_) * y_size_;
}

void PeriodicBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_signed(Real x_i,
                                                                                     Real y_i,
                                                                                     Real x_j,
                                                                                     Real y_j,
                                                                                     Real &dx,
                                                                                     Real &dy)
{
  dx = x_j - x_i;
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  } else if (dx < -0.5 * x_size_)
  {
    dx += x_size_;
  }

  dy = y_j - y_i;
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  } else if (dy < -0.5 * y_size_)
  {
    dy += y_size_;
  }
}

void PeriodicBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_unsigned(Real x_i,
                                                                                       Real y_i,
                                                                                       Real x_j,
                                                                                       Real y_j,
                                                                                       Real &dx,
                                                                                       Real &dy)
{
  dx = std::fabs(x_j - x_i);
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  }

  dy = std::fabs(y_j - y_i);
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  }
}

void PeriodicBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_signed(const Vector2D &r_i,
                                                                                     const Vector2D &r_j,
                                                                                     Vector2D &dr_ij)
{
  dr_ij.x = r_j.x - r_i.x;
  if (dr_ij.x > 0.5 * x_size_)
  {
    dr_ij.x -= x_size_;
  } else if (dr_ij.x < -0.5 * x_size_)
  {
    dr_ij.x += x_size_;
  }

  dr_ij.y = r_j.y - r_i.y;
  if (dr_ij.y > 0.5 * y_size_)
  {
    dr_ij.y -= y_size_;
  } else if (dr_ij.y < -0.5 * y_size_)
  {
    dr_ij.y += y_size_;
  }
}

void PeriodicBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_unsigned(const Vector2D &r_i,
                                                                                       const Vector2D &r_j,
                                                                                       Vector2D &dr_ij)
{
  dr_ij.x = std::fabs(r_j.x - r_i.x);
  if (dr_ij.x > 0.5 * x_size_)
  {
    dr_ij.x -= x_size_;
  }

  dr_ij.y = std::fabs(r_j.y - r_i.y);
  if (dr_ij.y > 0.5 * y_size_)
  {
    dr_ij.y -= y_size_;
  }
}

//The minimum image convention for the calculation of effective particle distances
//if the sign of the distance is relevant
void PeriodicBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_signed(Real x_i,
                                                                                     Real y_i,
                                                                                     Real x_j,
                                                                                     Real y_j,
                                                                                     Real &dx,
                                                                                     Real &dy)
{
  dx = x_j - x_i;
  dx -= x_size_ * std::nearbyint(dx * x_rsize_);

  dy = y_j - y_i;
  dy -= y_size_ * std::nearbyint(dy * y_rsize_);
}

//The minimum image convention for the calculation of effective particle distances
//if the sign of the distances is not relevant
void PeriodicBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_unsigned(Real x_i,
                                                                                       Real y_i,
                                                                                       Real x_j,
                                                                                       Real y_j,
                                                                                       Real &dx,
                                                                                       Real &dy)
{
  dx = std::fabs(x_j - x_i);
  dx -= static_cast<int>(dx * x_rsize_ + 0.5) * x_size_;

  dy = std::fabs(y_j - y_i);
  dy -= static_cast<int>(dy * y_rsize_ + 0.5) * y_size_;
}

void PeriodicBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_signed(const Vector2D &r_i,
                                                                                     const Vector2D &r_j,
                                                                                     Vector2D &dr_ij)
{
  dr_ij.x = r_j.x - r_i.x;
  dr_ij.x -= x_size_ * std::nearbyint(dr_ij.x * x_rsize_);

  dr_ij.y = r_j.y - r_i.y;
  dr_ij.y -= y_size_ * std::nearbyint(dr_ij.y * y_rsize_);
}

void PeriodicBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_unsigned(const Vector2D &r_i,
                                                                                       const Vector2D &r_j,
                                                                                       Vector2D &dr_ij)
{
  dr_ij.x = std::fabs(r_j.x - r_i.x);
  dr_ij.x -= static_cast<int>(dr_ij.x * x_rsize_ + 0.5) * x_size_;

  dr_ij.y = std::fabs(r_j.y - r_i.y);
  dr_ij.y -= static_cast<int>(dr_ij.y * y_rsize_ + 0.5) * y_size_;
}

//Restrict particle coordinates to the simulation box
void PeriodicBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc)
{
  x_pbc = x - std::floor(x * x_rsize_) * x_size_;
  y_pbc = y - std::floor(y * y_rsize_) * y_size_;
}

//Restrict particle coordinates to the simulation box
void PeriodicBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state, int N)
{
#pragma unroll
  for (int alpha = 0; alpha < N; ++alpha)
  {
    system_state[kS * alpha] -= std::floor(system_state[kS * alpha] * x_rsize_) * x_size_;
    system_state[kS * alpha + 1] -= std::floor(system_state[kS * alpha + 1] * y_rsize_) * y_size_;
  }
}

void PeriodicBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(const Vector2D &r, Vector2D &r_pbc)
{
  r_pbc.x = r.x - std::floor(r.x * x_rsize_) * x_size_;
  r_pbc.y = r.y - std::floor(r.y * y_rsize_) * y_size_;
}