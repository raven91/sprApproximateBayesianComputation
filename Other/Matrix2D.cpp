//
// Created by Nikita Kruk on 25.11.17.
//

#include "Matrix2D.hpp"

Matrix2D::Matrix2D() :
    xx(0.0),
    xy(0.0),
    yx(0.0),
    yy(0.0)
{

}

Matrix2D::Matrix2D(Real _xx, Real _xy, Real _yx, Real _yy) :
    xx(_xx),
    xy(_xy),
    yx(_yx),
    yy(_yy)
{

}

Matrix2D::Matrix2D(const Matrix2D &mat) :
    xx(mat.xx),
    xy(mat.xy),
    yx(mat.yx),
    yy(mat.yy)
{

}

Matrix2D::~Matrix2D()
{

}

Matrix2D &Matrix2D::operator+=(const Matrix2D &v)
{
  this->xx += v.xx;
  this->xy += v.xy;
  this->yx += v.yx;
  this->yy += v.yy;

  return *this;
}

Matrix2D &Matrix2D::operator-=(const Matrix2D &v)
{
  this->xx -= v.xx;
  this->xy -= v.xy;
  this->yx -= v.yx;
  this->yy -= v.yy;

  return *this;
}

Matrix2D &Matrix2D::operator*=(Real a)
{
  this->xx *= a;
  this->xy *= a;
  this->yx *= a;
  this->yy *= a;

  return *this;
}

Matrix2D &Matrix2D::operator/=(Real a)
{
  this->xx /= a;
  this->xy /= a;
  this->yx /= a;
  this->yy /= a;

  return *this;
}

Matrix2D operator*(Real a, const Matrix2D &v)
{
  Matrix2D w;
  w.xx = a * v.xx;
  w.xy = a * v.xy;
  w.yx = a * v.yx;
  w.yy = a * v.yy;

  return w;
}

Matrix2D operator*(const Matrix2D &v, Real a)
{
  Matrix2D w;
  w.xx = v.xx * a;
  w.xy = v.xy * a;
  w.yx = w.yx * a;
  w.yy = w.yy * a;

  return w;
}

Matrix2D operator+(const Matrix2D &v, const Matrix2D &u)
{
  Matrix2D w;
  w.xx = v.xx + u.xx;
  w.xy = v.xy + u.xy;
  w.yx = v.yx + u.yx;
  w.yy = v.yy + u.yy;

  return w;
}

Matrix2D operator-(const Matrix2D &v, const Matrix2D &u)
{
  Matrix2D w;
  w.xx = v.xx - u.xx;
  w.xy = v.xy - u.xy;
  w.yx = v.yx - u.yx;
  w.yy = v.yy - u.yy;

  return w;
}