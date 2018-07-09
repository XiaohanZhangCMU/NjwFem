#ifndef __SVECTOR_H__
#define __SVECTOR_H__

#include<iostream>
#include<cmath>

// Small vector (svector) class (everything inlined)

template<const int Ndim> class svector;
template<const int Ndim> class smatrix;
template<const int Ndim> class sbmatrix;
template<const int Ndim> 
           svector<Ndim> operator*(const double a, const svector<Ndim>& v);

template<const int Ndim>
class svector 
{
public:
  svector();
  svector(const svector&);
  explicit svector(double *);
  svector(double);        // only for Ndim = 1, otherwise double -> svectors
  svector(double,double);
  svector(double,double,double);

  svector& zero();
  double   dot(const svector&) const;
  double norm() const;

  svector  operator+(const svector&) const;
  svector& operator+=(const svector&);
  svector& operator-=(const svector&);
  svector  operator-() const;
  svector  operator-(const svector&) const;
  double   operator[](int) const;
  double&  operator[](int);

private:
  double x[Ndim];

  friend class smatrix<Ndim>;
  friend class sbmatrix<Ndim>;
  friend svector<Ndim> operator*<Ndim>(const double a, const svector<Ndim>& v);
};

template<const int Ndim> 
std::ostream& operator<<(std::ostream& out, svector<Ndim> x);

svector<3> cross(svector<3> &a, svector<3> &b);   // 3d cross product

#endif                          // !defined(__SVECTOR_H__)


