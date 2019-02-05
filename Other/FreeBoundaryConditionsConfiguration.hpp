//
// Created by Nikita Kruk on 08.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_FREEBOUNDARYCONDITIONSCONFIGURATION_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_FREEBOUNDARYCONDITIONSCONFIGURATION_HPP

#include "../Definitions.hpp"
#include "BoundaryConditionsConfiguration.hpp"
#include "Vector2D.hpp"

#include <vector>

class FreeBoundaryConditionsConfiguration : public BoundaryConditionsConfiguration
{
 public:

  FreeBoundaryConditionsConfiguration(Real x_size, Real y_size);
  ~FreeBoundaryConditionsConfiguration();

  void BoundedParticleDistance(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  void BoundedParticleDistance(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij);

  virtual void ClassAEffectiveParticleDistance(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  virtual void ClassAEffectiveParticleDistance(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij);

  virtual void ClassBEffectiveParticleDistance_signed(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  virtual void ClassBEffectiveParticleDistance_unsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  virtual void ClassBEffectiveParticleDistance_signed(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij);
  virtual void ClassBEffectiveParticleDistance_unsigned(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij);

  virtual void ClassCEffectiveParticleDistance_signed(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  virtual void ClassCEffectiveParticleDistance_unsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  virtual void ClassCEffectiveParticleDistance_signed(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij);
  virtual void ClassCEffectiveParticleDistance_unsigned(const Vector2D &r_i, const Vector2D &r_j, Vector2D &dr_ij);

  virtual void ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc);
  virtual void ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state, int N);
  virtual void ApplyPeriodicBoundaryConditions(const Vector2D &r, Vector2D &r_pbc);

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_FREEBOUNDARYCONDITIONSCONFIGURATION_HPP
