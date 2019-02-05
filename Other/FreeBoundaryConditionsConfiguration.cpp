//
// Created by Nikita Kruk on 08.08.18.
//

#include "FreeBoundaryConditionsConfiguration.hpp"

#include <cmath>
#include <iostream>

FreeBoundaryConditionsConfiguration::FreeBoundaryConditionsConfiguration(Real x_size, Real y_size) :
    BoundaryConditionsConfiguration(x_size, y_size)
{

}

FreeBoundaryConditionsConfiguration::~FreeBoundaryConditionsConfiguration()
{

}

void FreeBoundaryConditionsConfiguration::BoundedParticleDistance(Real x_i,
                                                                  Real y_i,
                                                                  Real x_j,
                                                                  Real y_j,
                                                                  Real &dx,
                                                                  Real &dy)
{
  dx = x_j - x_i;
  dy = y_j - y_i;
}

void FreeBoundaryConditionsConfiguration::BoundedParticleDistance(const Vector2D &r_i,
                                                                  const Vector2D &r_j,
                                                                  Vector2D &dr_ij)
{
  dr_ij.x = r_j.x - r_i.x;
  dr_ij.y = r_j.y - r_i.y;
}

void FreeBoundaryConditionsConfiguration::ClassAEffectiveParticleDistance(Real x_i,
                                                                          Real y_i,
                                                                          Real x_j,
                                                                          Real y_j,
                                                                          Real &dx,
                                                                          Real &dy)
{
  BoundedParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
}

void FreeBoundaryConditionsConfiguration::ClassAEffectiveParticleDistance(const Vector2D &r_i,
                                                                          const Vector2D &r_j,
                                                                          Vector2D &dr_ij)
{
  BoundedParticleDistance(r_i, r_j, dr_ij);
}

void FreeBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_signed(Real x_i,
                                                                                 Real y_i,
                                                                                 Real x_j,
                                                                                 Real y_j,
                                                                                 Real &dx,
                                                                                 Real &dy)
{
  BoundedParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
}

void FreeBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_unsigned(Real x_i,
                                                                                   Real y_i,
                                                                                   Real x_j,
                                                                                   Real y_j,
                                                                                   Real &dx,
                                                                                   Real &dy)
{
  BoundedParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
}

void FreeBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_signed(const Vector2D &r_i,
                                                                                 const Vector2D &r_j,
                                                                                 Vector2D &dr_ij)
{
  BoundedParticleDistance(r_i, r_j, dr_ij);
}

void FreeBoundaryConditionsConfiguration::ClassBEffectiveParticleDistance_unsigned(const Vector2D &r_i,
                                                                                   const Vector2D &r_j,
                                                                                   Vector2D &dr_ij)
{
  BoundedParticleDistance(r_i, r_j, dr_ij);
}

void FreeBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_signed(Real x_i,
                                                                                 Real y_i,
                                                                                 Real x_j,
                                                                                 Real y_j,
                                                                                 Real &dx,
                                                                                 Real &dy)
{
  BoundedParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
}

void FreeBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_unsigned(Real x_i,
                                                                                   Real y_i,
                                                                                   Real x_j,
                                                                                   Real y_j,
                                                                                   Real &dx,
                                                                                   Real &dy)
{
  BoundedParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
}

void FreeBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_signed(const Vector2D &r_i,
                                                                                 const Vector2D &r_j,
                                                                                 Vector2D &dr_ij)
{
  BoundedParticleDistance(r_i, r_j, dr_ij);
}

void FreeBoundaryConditionsConfiguration::ClassCEffectiveParticleDistance_unsigned(const Vector2D &r_i,
                                                                                   const Vector2D &r_j,
                                                                                   Vector2D &dr_ij)
{
  BoundedParticleDistance(r_i, r_j, dr_ij);
}

// Does not do anything in the free boundary conditions setup
void FreeBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc)
{
  x_pbc = x;
  y_pbc = y;
}

// Does not do anything in the free boundary conditions setup
void FreeBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state, int N)
{

}

// Does not do anything in the free boundary conditions setup
void FreeBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(const Vector2D &r, Vector2D &r_pbc)
{
  r_pbc = r;
}