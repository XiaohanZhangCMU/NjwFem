#ifndef __SMATRIX_C__
#define __SMATRIX_C__
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
  smatrix cofactor() const;       // return the cofactor

  smatrix operator+(const smatrix&) const;
  smatrix& operator+=(const smatrix&);
  smatrix& operator-=(const smatrix&);
  //xiaohan
  smatrix operator*(double scalar) const;
  smatrix& operator*=(double scalar);
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


template<const int Ndim>
inline double smatrix<Ndim>::frobNorm() const
{
  return sqrt(frobDot(*this));  // compute the froebenius norm of this matrix
}

template<const int Ndim>
inline double smatrix<Ndim>::trace() const
{
  double tr = 0.0;

  for(int i = 0; i < Ndim; i++) tr += m[i][i];

  return(tr);
}

template<const int Ndim>
inline double smatrix<Ndim>::frobDot(const smatrix<Ndim>& a) const
{
  double fDot = 0.0;

// compute the froebenius inner produc of this matrix

  for (int i = 0; i < Ndim; i++) fDot += m[i].dot(a[i]);

  return fDot;
}

template<const int Ndim>
inline smatrix<Ndim>::smatrix()
{
  for(int i = 0; i < Ndim; i++) m[i].zero();
}

template<const int Ndim>
inline smatrix<Ndim>::smatrix(const smatrix<Ndim>& a)
{
  for(int i = 0; i < Ndim; i++) m[i] = a.m[i];
}

template<const int Ndim>
inline smatrix<Ndim>::smatrix(const double *a)
{
  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim; j++) m[i][j] = a[i*Ndim + j];
}

template<const int Ndim>
inline smatrix<Ndim>::smatrix(const double a[Ndim][Ndim])
{
  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim; j++) m[i][j] = a[i][j];
}

template<const int Ndim>
inline smatrix<Ndim>::smatrix(const svector<Ndim>& v, const svector<Ndim>& w)
{
  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim; j++) m[i][j] = v.x[i] * w.x[j];
}


template<const int Ndim>
inline smatrix<Ndim>::smatrix(const double a00, const double a01, 
			      const double a10, const double a11)
{
  assert(Ndim == 2);

  m[0][0] = a00; m[0][1] = a01;
  m[1][0] = a10; m[1][1] = a11;
}

template<const int Ndim>
inline 
smatrix<Ndim>::smatrix(const double a00, const double a01, const double a02, 
		       const double a10, const double a11, const double a12, 
		       const double a20, const double a21, const double a22)
{
  assert(Ndim == 3);

  m[0][0] = a00; m[0][1] = a01; m[0][2] = a02;
  m[1][0] = a10; m[1][1] = a11; m[1][2] = a12;
  m[2][0] = a20; m[2][1] = a21; m[2][2] = a22;
}

template<const int Ndim>
inline smatrix<Ndim>& smatrix<Ndim>::zero()
{
  for(int i = 0; i < Ndim; i++) m[i].zero();

  return(*this);
}

template<const int Ndim>
inline smatrix<Ndim>& smatrix<Ndim>::identity()
{
  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim; j++) 
    {
      if(i == j) m[i][j] = 1.0;
      else       m[i][j] = 0.0;
    }

  return(*this);
}

template<const int Ndim>
inline smatrix<Ndim> smatrix<Ndim>::transpose() const
{
  smatrix<Ndim> a;

  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim; j++) a.m[i][j] = m[j][i];
  return(a);
}

template<const int Ndim>
inline smatrix<Ndim> smatrix<Ndim>::operator+(const smatrix<Ndim>& a) const
{
  smatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = m[i] + a.m[i];

  return(b);
}

template<const int Ndim>
inline smatrix<Ndim>& smatrix<Ndim>::operator+=(const smatrix<Ndim>& a)
{
  for(int i = 0; i < Ndim; i++) m[i] += a.m[i];

  return(*this);
}
//xiaohan


template<const int Ndim>
inline smatrix<Ndim> smatrix<Ndim>::operator*(double scal) const
{
  smatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = m[i]*scal;

  return(b);
}


template<const int Ndim>
inline smatrix<Ndim>& smatrix<Ndim>::operator-=(const smatrix<Ndim>& a)
{
  for(int i = 0; i < Ndim; i++) m[i] -= a.m[i];

  return(*this);
}

template<const int Ndim>
inline smatrix<Ndim>& smatrix<Ndim>::operator*=(double scalar)
{
  for(int i = 0;i<Ndim;i++) m[i] *=scalar;

  return (*this);
}
template<const int Ndim>
inline smatrix<Ndim> smatrix<Ndim>::operator-() const
{
  smatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = -m[i];

  return(b);
}

template<const int Ndim>
inline smatrix<Ndim> smatrix<Ndim>::operator-(const smatrix<Ndim>& a) const
{
  smatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = m[i] - a.m[i];

  return(b);
}

template<const int Ndim>
smatrix<Ndim> operator*(double c, const smatrix<Ndim>& a)
{
  smatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = c * a.m[i];

  return(b);
}

template<const int Ndim>
svector<Ndim> smatrix<Ndim>::operator*(const svector<Ndim>& v) const
{
  svector<Ndim> w;

  for(int i = 0; i < Ndim; i++) w.x[i] += m[i].dot(v);

  return(w);
}

template<const int Ndim>
smatrix<Ndim> smatrix<Ndim>::operator*(const smatrix<Ndim>& a) const
{
  int i,j,k;
  smatrix<Ndim> b;

  for(i = 0; i < Ndim; i++)
  for(j = 0; j < Ndim; j++) 
    {
      b.m[i][j] = 0.0;
      for(k = 0; k < Ndim; k++) b.m[i][j] += m[i][k] * a.m[k][j];
    }
  return(b);
}

template<const int Ndim>
inline svector<Ndim>& smatrix<Ndim>::operator[](int i)
{
  return(m[i]);
}

template<const int Ndim>
inline const svector<Ndim> smatrix<Ndim>::operator[](int i) const
{
  // Note: Returning "const svector<Ndim>" guarantees that the
  //  syntax "a[i][j] = c" will not just update a copy of m[i]

  return(m[i]);
}

template<const int Ndim>
smatrix<Ndim> smatrix<Ndim>::inverse(double &det) const
{
  smatrix<Ndim> a;

if(Ndim == 1)
  {
  if(m[0][0] == 0.0) 
    goto error;
  else 
    {
      a.m[0][0] = 1 / m[0][0];
      det = m[0][0];
    }
  }
else if(Ndim == 2)
 {
 det = m[0][0] * m[1][1] - m[0][1] * m[1][0];

 if( det == 0.0 ) goto error;

 a.m[0][0] = m[1][1] / det;
 a.m[1][1] = m[0][0] / det;

 a.m[0][1] = -m[0][1] / det;
 a.m[1][0] = -m[1][0] / det;
 }
else if(Ndim == 3)
  {
  det =  m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1]
       - m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1]
       + m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1];
 
  if( det == 0.0 ) goto error;
 
  a.m[0][0] = ( m[1][1]*m[2][2] - m[1][2]*m[2][1]) / det;
  a.m[0][1] = (-m[0][1]*m[2][2] + m[0][2]*m[2][1]) / det;
  a.m[0][2] = ( m[0][1]*m[1][2] - m[0][2]*m[1][1]) / det;
  a.m[1][0] = (-m[1][0]*m[2][2] + m[1][2]*m[2][0]) / det;
  a.m[1][1] = ( m[0][0]*m[2][2] - m[0][2]*m[2][0]) / det;
  a.m[1][2] = (-m[0][0]*m[1][2] + m[0][2]*m[1][0]) / det;
  a.m[2][0] = ( m[1][0]*m[2][1] - m[1][1]*m[2][0]) / det;
  a.m[2][1] = (-m[0][0]*m[2][1] + m[0][1]*m[2][0]) / det;
  a.m[2][2] = ( m[0][0]*m[1][1] - m[0][1]*m[1][0]) / det;
  }
else
  {
    std::cerr 
      << "smatrix<>::invert() not yet implemented in dimension Ndim \n";

  throw("smatrix<>::invert() not yet implemented in dimension Ndim \n");
  }

return(a);

error:
 {
   std::cerr << "Error in invert(): zero determinant \n";
   throw("Error in invert(): zero determinant \n");
 }
}

template<const int Ndim>
smatrix<Ndim> smatrix<Ndim>::cofactor() const
{
  smatrix<Ndim> a;

if(Ndim == 1)
  {
    a.m[0][0] = 1 / m[0][0];
  }
else if(Ndim == 2)
 {
 a.m[0][0] = m[1][1];
 a.m[1][1] = m[0][0];

 a.m[1][0] = -m[0][1];
 a.m[0][1] = -m[1][0];
 }
else if(Ndim == 3)
  {
  a.m[0][0] = ( m[1][1]*m[2][2] - m[1][2]*m[2][1]);
  a.m[1][0] = (-m[0][1]*m[2][2] + m[0][2]*m[2][1]);
  a.m[2][0] = ( m[0][1]*m[1][2] - m[0][2]*m[1][1]);
  a.m[0][1] = (-m[1][0]*m[2][2] + m[1][2]*m[2][0]);
  a.m[1][1] = ( m[0][0]*m[2][2] - m[0][2]*m[2][0]);
  a.m[2][1] = (-m[0][0]*m[1][2] + m[0][2]*m[1][0]);
  a.m[0][2] = ( m[1][0]*m[2][1] - m[1][1]*m[2][0]);
  a.m[1][2] = (-m[0][0]*m[2][1] + m[0][1]*m[2][0]);
  a.m[2][2] = ( m[0][0]*m[1][1] - m[0][1]*m[1][0]);
  }
else
  {
    std::cerr 
      << "smatrix<>::invert() not yet implemented in dimension Ndim \n";

  throw("smatrix<>::invert() not yet implemented in dimension Ndim \n");
  }

return(a);
}

template<const int Ndim>
double smatrix<Ndim>::determinant() const
{

  if(Ndim == 1)      return(m[0][0]);
  else if(Ndim == 2) return(m[0][0] * m[1][1] - m[0][1] * m[1][0]);
  else if(Ndim == 3)
    return( m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] -
	    m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1] +
	    m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1]);
  else
    {
      std::cerr 
	<< "smatrix<>::determiant() not yet implemented in dimension Ndim \n";

    throw("smatrix<>::determinant() not yet implemented in dimension Ndim \n");
    }
}

template<const int Ndim> 
std::ostream& operator<<(std::ostream& out, smatrix<Ndim> m)
{
  out << "[";

  for(int i = 0; i < Ndim; i++)
    {
      out << m[i];
      if(i < Ndim-1) out << ";";
    }

  out << "]";

  return(out);
}

#endif                          // !defined(__SMATRIX_C__)
