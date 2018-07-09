#ifndef __SMATRIX_H__
#define __SMATRIX_H__
#include<cassert>
#include<iostream>

// Small matrix (smatrix) class (everything inlined)

template<const int Ndim> class svector;
template<const int Ndim> class smatrix;

template<const int Ndim>
smatrix<Ndim> operator*(const double c, const smatrix<Ndim>& a);

template<const int Ndim>
class smatrix 
{
public:
  smatrix();
  smatrix(const smatrix&);
  smatrix(const double *);  // a[i][j] = arg[i*Ndim + j]
  smatrix(const double a[Ndim][Ndim]);
  smatrix(const svector<Ndim>&, const svector<Ndim>&);  // outer product

  smatrix(const double, const double, 
	  const double, const double);  // 2x2 matrix, rows
  smatrix(const double, const double, const double, 
	  const double, const double, const double, 
	  const double, const double, const double);  // 3x3 matrix, rows

  smatrix& zero();                // reset to zero
  smatrix& identity();            // reset to identity
  smatrix transpose() const;      // return the transpose
  smatrix inverse(double&) const; // return the inverse and determinant
  double  determinant() const;    // return the determinant

  smatrix operator+(const smatrix&) const;
  smatrix& operator+=(const smatrix&);
  smatrix operator-() const;
  smatrix operator-(const smatrix&) const;
  svector<Ndim> operator*(const svector<Ndim>&) const;
  smatrix operator*(const smatrix&) const;

  double trace() const; 
  double frobNorm() const;
  double frobDot(const smatrix&) const;

  const svector<Ndim> operator[](int i) const;
  svector<Ndim>& operator[](int i); 

    //friend smatrix<Ndim> operator* <> (const double , const smatrix<Ndim> &);


  svector<Ndim> m[Ndim];   // m[i] = i-th row of the matrix
  private:

};

#endif                          // !defined(__SMATRIX_H__)
