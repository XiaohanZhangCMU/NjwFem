#ifndef __SBMATRIX_C__
#define __SBMATRIX_C__

#include<iostream>
#include<cassert>
#include<cmath>

// Small "boundary" matrix (sbmatrix) class (everything inlined)
// Construct an (Ndim x Ndim-1) matrix used for parameterizing 
// surfaces.

template<const int Ndim> class svector;
template<const int Ndim> class smatrix;
template<const int Ndim> class sbmatrix;
template<const int Ndim>
           sbmatrix<Ndim> operator*(double c, const sbmatrix<Ndim>& a);

template<const int Ndim>
           sbmatrix<Ndim> operator*(const smatrix<Ndim>& a, 
				    const sbmatrix<Ndim>& b);

template<const int Ndim>
class sbmatrix 
{
public:
  sbmatrix();
  sbmatrix(const sbmatrix&);
  sbmatrix(const double a[Ndim][Ndim-1]);
  sbmatrix(const svector<Ndim>&, const svector<Ndim-1>&);  // outer product

  sbmatrix& zero();                // reset to zero

  sbmatrix operator+(const sbmatrix&) const;
  sbmatrix& operator+=(const sbmatrix&);
  sbmatrix operator-() const;
  sbmatrix operator-(const sbmatrix&) const;
  svector<Ndim> operator*(const svector<Ndim-1>&) const;

  double frobNorm() const;

  const svector<Ndim-1> operator[](int i) const;
  svector<Ndim-1>& operator[](int i); 

  svector<Ndim> normal(double &det) const;

  // private:
  svector<Ndim-1> m[Ndim];   // m[i] = i-th row of the matrix

  // friend sbmatrix<Ndim> operator*<Ndim>(double c, const sbmatrix<Ndim> &a);
  // friend sbmatrix<Ndim> operator*<Ndim>(const smatrix<Ndim> &a, 
  //		  			   const sbmatrix<Ndim> &b);
};


template<const int Ndim>
inline double sbmatrix<Ndim>::frobNorm() const
{
  double norm = 0.0;

  for (int i = 0; i < Ndim; i++) norm += m[i].dot(m[i]);

  return sqrt(norm);  // compute the froebenius norm of this matrix
}

template<const int Ndim>
inline sbmatrix<Ndim>::sbmatrix()
{
  for(int i = 0; i < Ndim; i++) m[i].zero();
}

template<const int Ndim>
inline sbmatrix<Ndim>::sbmatrix(const sbmatrix<Ndim>& a)
{
  for(int i = 0; i < Ndim; i++) m[i] = a.m[i];
 }

template<const int Ndim>
inline sbmatrix<Ndim>::sbmatrix(const double a[Ndim][Ndim-1])
{
  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim-1; j++) m[i][j] = a[i][j];
}

template<const int Ndim>
inline sbmatrix<Ndim>::sbmatrix(const svector<Ndim>& v, 
				const svector<Ndim-1>& w)
{
  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim-1; j++) m[i][j] = v[i] * w[j];
}

template<const int Ndim>
inline sbmatrix<Ndim>& sbmatrix<Ndim>::zero()
{
  for(int i = 0; i < Ndim; i++) m[i].zero();

  return(*this);
}

template<const int Ndim>
inline sbmatrix<Ndim> sbmatrix<Ndim>::operator+(const sbmatrix<Ndim>& a) const
{
  sbmatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = m[i] + a.m[i];

  return(b);
}

template<const int Ndim>
inline sbmatrix<Ndim>& sbmatrix<Ndim>::operator+=(const sbmatrix<Ndim>& a)
{
  for(int i = 0; i < Ndim; i++) m[i] += a.m[i];

  return(*this);
}

template<const int Ndim>
inline sbmatrix<Ndim> sbmatrix<Ndim>::operator-() const
{
  sbmatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = -m[i];

  return(b);
}

template<const int Ndim>
inline sbmatrix<Ndim> sbmatrix<Ndim>::operator-(const sbmatrix<Ndim>& a) const
{
  sbmatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = m[i] - a.m[i];

  return(b);
}

template<const int Ndim>
inline sbmatrix<Ndim> operator*(double c, const sbmatrix<Ndim>& a)
{
  sbmatrix<Ndim> b;

  for(int i = 0; i < Ndim; i++) b.m[i] = c * a.m[i];

  return(b);
}

template<const int Ndim>
inline sbmatrix<Ndim> operator*(const smatrix<Ndim>& a,
				const sbmatrix<Ndim>& b)
{
  sbmatrix<Ndim> c;

  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim-1; j++)
    {
      c[i][j] = 0.0;

      for(int k = 0; k < Ndim; k++) c[i][j] += a[i][k]*b[k][j];
    }

  return(c);
}

template<const int Ndim>
inline svector<Ndim> sbmatrix<Ndim>::operator*(const svector<Ndim-1>& v) const
{
  svector<Ndim> w;

  for(int i = 0; i < Ndim; i++) 
    {
      w.x[i] = 0.0;

      for(int j = 0; j < Ndim-1; j++) w.x[i] += m[i][j] * v.x[j];
    }

  return(w);
}


template<const int Ndim>
inline svector<Ndim-1>& sbmatrix<Ndim>::operator[](int i)
{
  return(m[i]);
}

template<const int Ndim>
inline const svector<Ndim-1> sbmatrix<Ndim>::operator[](int i) const
{
  // Note: Returning "const svector<Ndim>" guarantees that the
  //  syntax "a[i][j] = c" will not just update a copy of m[i]

  return(m[i]);
}

template<const int Ndim>
svector<Ndim> sbmatrix<Ndim>::normal(double &det) const
{
  if(Ndim == 2)
    {
      svector<Ndim> nn(m[1][0], -m[0][0]);
      
      det = nn.norm();

      assert(det != 0.0);

      return( (1.0/det) * nn );
    }
  else if(Ndim == 3)
    {
      svector<Ndim> nn(m[1][0]*m[2][1]-m[2][0]*m[1][1],
		       m[2][0]*m[0][1]-m[0][0]*m[2][1],
		       m[0][0]*m[1][1]-m[1][0]*m[0][1]);  // cross product

      det = nn.norm();

      assert(det != 0.0);

      return( (1/det)*nn );
    }
  else
    {
      throw("sbmatrix<Ndim>::normal() only implemented for Ndim = 2 and 3");
    }
}
    
#endif                          // !defined(__SBMATRIX_C__)
