//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEM_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEM_HPP

#include "../Definitions.hpp"
#include "../ParameterSets/AbstractParameterSet.hpp"
#include "../Other/BoundaryConditionsConfiguration.hpp"
#include "../Other/Vector2D.hpp"
#include "TwoParticleInteractionForce.hpp"

#include <vector>
#include <set>
#include <random>

class SprSystemWithUnitVelocity
{
 public:

  explicit SprSystemWithUnitVelocity(AbstractParameterSet *parameter_set, BoundaryConditionsConfiguration *bc_config);
  virtual ~SprSystemWithUnitVelocity();

  void OperatorLinkedList(std::vector<Real> &system_state, std::vector<Real> &derivative, Real t);
  void OperatorVerletNeighborList(std::vector<Real> &system_state, std::vector<Real> &derivative, Real t);
  virtual void OperatorNoise(std::vector<Real> &system_state, std::vector<Real> &derivative, Real t);

 protected:

  AbstractParameterSet *parameter_set_;
  BoundaryConditionsConfiguration *bc_config_;
  TwoParticleInteractionForce two_particle_interaction_force_;

  //linked list and cell division
  struct SegmentIndex
  {
    int alpha;//index of a rod which a segment belongs to
    int i;//index of a segment inside a rod

    SegmentIndex() : alpha(0), i(0)
    {}
    SegmentIndex(int _alpha, int _i) : alpha(_alpha), i(_i)
    {}
  };

  Real x_size_;
  Real y_size_;
  int num_subcells_x_;
  int num_subcells_y_;
  std::vector<std::vector<SegmentIndex>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;//list of index increments for all of the 9 neighboring subsells
  std::vector<Vector2D> position_interaction_force_;
  std::vector<Vector2D> velocity_interaction_force_;

  //Verlet neighbor list
  Real r_max_;//skin region
  Real r_min_;//cut-off distance
  std::vector<std::vector<SegmentIndex>> verlet_list_;
  Real accumulated_displacement_;
  bool should_update_lists_;
  std::mt19937 mersenne_twister_generator_;

  virtual void CalculateMaxDisplacement(const std::vector<Real> &system_state);
  void ConstructNewLinkedList(const std::vector<Real> &system_state);
  void ConstructNewVerletList(const std::vector<Real> &system_state);
  void CalculateInteractionForcesUsingLinkedList(const std::vector<Real> &system_state);
  void CalculateInteractionForcesUsingVerletList(const std::vector<Real> &system_state);
  virtual void CalculateDerivatives(const std::vector<Real> &system_state, std::vector<Real> &derivative);
  virtual void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y);
  virtual bool ShouldConsiderFreeBoundaries();
  virtual bool ShouldConsiderHalfOfNeighbors();

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEM_HPP
