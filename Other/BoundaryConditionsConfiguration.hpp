//
// Created by Nikita Kruk on 08.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_BOUNDARYCONDITIONSCONFIGURATION_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_BOUNDARYCONDITIONSCONFIGURATION_HPP

#include "../Definitions.hpp"
#include "Vector2D.hpp"

#include <vector>

class BoundaryConditionsConfiguration
{
 public:

  BoundaryConditionsConfiguration(Real x_size, Real y_size);
  virtual ~BoundaryConditionsConfiguration();

  void SetNewBoundaries(Real x_size, Real y_size);

  virtual void ClassAEffectiveParticleDistance(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) = 0;
  virtual void ClassAEffectiveParticleDistance(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij) = 0;

  virtual void ClassBEffectiveParticleDistance_signed(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) = 0;
  virtual void ClassBEffectiveParticleDistance_unsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) = 0;
  virtual void ClassBEffectiveParticleDistance_signed(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij) = 0;
  virtual void ClassBEffectiveParticleDistance_unsigned(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij) = 0;

  virtual void ClassCEffectiveParticleDistance_signed(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) = 0;
  virtual void ClassCEffectiveParticleDistance_unsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy) = 0;
  virtual void ClassCEffectiveParticleDistance_signed(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij) = 0;
  virtual void ClassCEffectiveParticleDistance_unsigned(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij) = 0;

  virtual void ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc) = 0;
  virtual void ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state, int N) = 0;
  virtual void ApplyPeriodicBoundaryConditions(const Vector2D &r, Vector2D &r_pbc) = 0;

 protected:

  Real x_size_;
  Real y_size_;
  Real x_rsize_;
  Real y_rsize_;

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_BOUNDARYCONDITIONSCONFIGURATION_HPP
