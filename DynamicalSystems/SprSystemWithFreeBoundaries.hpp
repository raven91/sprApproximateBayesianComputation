//
// Created by Nikita Kruk on 08.08.18.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEMWITHFREEBOUNDARIES_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEMWITHFREEBOUNDARIES_HPP

#include "../Definitions.hpp"
#include "SprSystemWithArbitraryStationaryVelocity.hpp"
#include "../ParameterSets/AbstractParameterSet.hpp"
#include "../Other/BoundaryConditionsConfiguration.hpp"

class SprSystemWithFreeBoundaries : public SprSystemWithArbitraryStationaryVelocity
{
 public:

  SprSystemWithFreeBoundaries(AbstractParameterSet *parameter_set,
                              BoundaryConditionsConfiguration *bc_config,
                              std::vector<Real> &velocity_distribution);
  virtual ~SprSystemWithFreeBoundaries();

 protected:

  virtual void CalculateMaxDisplacement(const std::vector<Real> &system_state);
  virtual void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y);
  virtual bool ShouldConsiderFreeBoundaries();
  virtual bool ShouldConsiderHalfOfNeighbors();

};

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_SPRSYSTEMWITHFREEBOUNDARIES_HPP
