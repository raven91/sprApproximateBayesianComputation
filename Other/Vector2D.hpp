//
// Created by Nikita Kruk on 25.11.17.
//

#ifndef SPRAPPROXIMATEBAYESIANCOMPUTATION_VECTOR2D_HPP
#define SPRAPPROXIMATEBAYESIANCOMPUTATION_VECTOR2D_HPP

#include "../Definitions.hpp"
#include "Matrix2D.hpp"

struct Matrix2D;

class Vector2D
{
 public:

  Vector2D();
  Vector2D(Real _x, Real _y);
  Vector2D(const Vector2D &vec);
  ~Vector2D();

  Vector2D &operator=(const Vector2D &v);
  Vector2D &operator+=(const Vector2D &v);
  Vector2D &operator-=(const Vector2D &v);
  Vector2D &operator*=(Real a);

  Real x;
  Real y;

};

Vector2D operator*(Real a, const Vector2D &v);
Vector2D operator*(const Vector2D &v, Real a);
Vector2D operator+(const Vector2D &v, const Vector2D &u);
Vector2D operator-(const Vector2D &v, const Vector2D &u);
Real Norm(const Vector2D &v);
Real NormSquared(const Vector2D &v);
void Normalize(Real &x, Real &y);
void Normalize(Vector2D &v);
Real InnerProduct(const Vector2D &v, const Vector2D &u);
Real CrossProduct3(const Vector2D &v, const Vector2D &u);
Matrix2D OuterProduct(const Vector2D &v, const Vector2D &u);

Vector2D operator*(const Matrix2D &mat, const Vector2D &vec);

#endif //SPRAPPROXIMATEBAYESIANCOMPUTATION_VECTOR2D_HPP
