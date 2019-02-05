//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_MATRIX2D_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_MATRIX2D_HPP

#include "../Definitions.hpp"
#include "Vector2D.hpp"

class Matrix2D
{
 public:

  Matrix2D();
  Matrix2D(Real _xx, Real _xy, Real _yx, Real _yy);
  Matrix2D(const Matrix2D &mat);
  ~Matrix2D();

  Matrix2D &operator+=(const Matrix2D &v);
  Matrix2D &operator-=(const Matrix2D &v);
  Matrix2D &operator*=(Real a);
  Matrix2D &operator/=(Real a);

  Real xx;
  Real xy;
  Real yx;
  Real yy;

};

Matrix2D operator*(Real a, const Matrix2D &v);
Matrix2D operator*(const Matrix2D &v, Real a);
Matrix2D operator+(const Matrix2D &v, const Matrix2D &u);
Matrix2D operator-(const Matrix2D &v, const Matrix2D &u);

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_MATRIX2D_HPP
