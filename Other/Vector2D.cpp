//
// Created by Nikita Kruk on 25.11.17.
//

#include "Vector2D.hpp"

#include <cmath>

Vector2D::Vector2D() :
    x(0.0),
    y(0.0)
{

}

Vector2D::Vector2D(Real _x, Real _y) :
    x(_x),
    y(_y)
{

}

Vector2D::Vector2D(const Vector2D &v) :
    x(v.x),
    y(v.y)
{

}

Vector2D::~Vector2D()
{

}

Vector2D &Vector2D::operator=(const Vector2D &v)
{
  if (&v == this)
  {
    return *this;
  }
  this->x = v.x;
  this->y = v.y;

  return *this;
}

Vector2D &Vector2D::operator+=(const Vector2D &v)
{
  this->x += v.x;
  this->y += v.y;

  return *this;
}

Vector2D &Vector2D::operator-=(const Vector2D &v)
{
  this->x -= v.x;
  this->y -= v.y;

  return *this;
}

Vector2D &Vector2D::operator*=(Real a)
{
  this->x *= a;
  this->y *= a;

  return *this;
}

Vector2D operator*(Real a, const Vector2D &v)
{
  Vector2D u;
  u.x = a * v.x;
  u.y = a * v.y;

  return u;
}

Vector2D operator*(const Vector2D &v, Real a)
{
  Vector2D u;
  u.x = v.x * a;
  u.y = v.y * a;

  return u;
}

Vector2D operator+(const Vector2D &v, const Vector2D &u)
{
  Vector2D w;
  w.x = v.x + u.x;
  w.y = v.y + u.y;

  return w;
}

Vector2D operator-(const Vector2D &v, const Vector2D &u)
{
  Vector2D w;
  w.x = v.x - u.x;
  w.y = v.y - u.y;

  return w;
}

Real Norm(const Vector2D &v)
{
  return sqrt(v.x * v.x + v.y * v.y);
}

Real NormSquared(const Vector2D &v)
{
  return v.x * v.x + v.y * v.y;
}

void Normalize(Real &x, Real &y)
{
  Real r = std::sqrt(x * x + y * y);
  if (r != 0.0)
  {
    x /= r;
    y /= r;
  }
}

void Normalize(Vector2D &v)
{
  Real r = Norm(v);
  if (r != 0.0)
  {
    v.x /= r;
    v.y /= r;
  }
}

Real InnerProduct(const Vector2D &v, const Vector2D &u)
{
  return v.x * u.x + v.y * u.y;
}

Real CrossProduct3(const Vector2D &v, const Vector2D &u)
{
  return v.x * u.y - v.y * u.x;
}

Matrix2D OuterProduct(const Vector2D &v, const Vector2D &u)
{
  Matrix2D outer_product;
  outer_product.xx = v.x * u.x;
  outer_product.xy = v.x * u.y;
  outer_product.yx = v.y * u.x;
  outer_product.yy = v.y * u.y;
  return outer_product;
}

Vector2D operator*(const Matrix2D &mat, const Vector2D &vec)
{
  Vector2D res;
  res.x = mat.xx * vec.x + mat.xy * vec.y;
  res.y = mat.yx * vec.x + mat.yy * vec.y;

  return res;
}