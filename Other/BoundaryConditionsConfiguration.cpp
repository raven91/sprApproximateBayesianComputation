//
// Created by Nikita Kruk on 08.08.18.
//

#include "BoundaryConditionsConfiguration.hpp"

BoundaryConditionsConfiguration::BoundaryConditionsConfiguration(Real x_size, Real y_size) :
    x_size_(x_size),
    y_size_(y_size),
    x_rsize_(1.0 / x_size),
    y_rsize_(1.0 / y_size)
{

}

BoundaryConditionsConfiguration::~BoundaryConditionsConfiguration()
{

}

void BoundaryConditionsConfiguration::SetNewBoundaries(Real x_size, Real y_size)
{
  x_size_ = x_size;
  y_size_ = y_size;

  x_rsize_ = 1.0 / x_size;
  y_rsize_ = 1.0 / y_size;
}