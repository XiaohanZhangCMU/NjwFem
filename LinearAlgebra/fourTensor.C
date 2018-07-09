#ifndef __FOURTENSOR_C__
#define __FOURTENSOR_C__

#include<iostream>

#include"../LinearAlgebra/svector.C"
#include"../LinearAlgebra/smatrix.C"

template<const int Ndim> class fourTensor;
template<const int Ndim>
fourTensor<Ndim> operator*(const double a, const fourTensor<Ndim>& v);

template<const int Ndim>
class fourTensor
{
public:
  fourTensor();
  fourTensor(const fourTensor<Ndim>&);
  fourTensor(const smatrix<Ndim>&, const smatrix<Ndim>&); //tensor product
  fourTensor(double mu, double lambda);
  fourTensor(double c11, double c12, double c44);

  fourTensor<Ndim>& zero();                  //reset to zero
  fourTensor<Ndim> identity();        //return identity

  fourTensor<Ndim> operator+(const fourTensor<Ndim>&) const;
  fourTensor<Ndim>& operator+=(const fourTensor<Ndim>&);
  fourTensor<Ndim> operator-() const;
  fourTensor<Ndim> operator-(const fourTensor<Ndim>&) const;
  fourTensor<Ndim>& operator=(const fourTensor<Ndim>& f);
  // double& operator()(int i, int j, int k, int l);
  double  operator()(int i, int j, int k, int l) const;
  double & operator()(int i, int j, int k, int l);
  smatrix<Ndim> operator*(const smatrix<Ndim>&) const;  // fourTensor * matrix
  smatrix<Ndim> contract(const smatrix<Ndim>& s, int idx1, int idx2) const;
  void setValue(int i, int j, int k, int l, double val);
  
  fourTensor<Ndim> compose(const fourTensor<Ndim>&);  // composition

  //  std::ostream& operator<<(std::ostream& out, fourTensor<Ndim> cc);

private:
  double cc[Ndim][Ndim][Ndim][Ndim];

  // I don't know how to make this work when the above operator*() function
  // exists
  //  friend fourTensor<Ndim> operator*<Ndim>(const double a, 
  //					  const fourTensor<Ndim>& v);
};

template<const int  Ndim>
inline fourTensor<Ndim>::fourTensor()   // default constructor
{
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0;k<Ndim;k++){
	for(int l = 0;l<Ndim;l++){
	  cc[i][j][k][l] = 0;
	}}}}
}

template<const int Ndim>                 // copy constructor
inline fourTensor<Ndim>::fourTensor(const fourTensor<Ndim>& f)
{
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0;k<Ndim;k++){
	for(int l = 0;l<Ndim;l++){
	  cc[i][j][k][l] = f(i,j,k,l);
	}}}}
}


template<const int Ndim>  // tensor product
inline fourTensor<Ndim>::fourTensor(const smatrix<Ndim>& s1, 
				    const smatrix<Ndim>& s2)
{
  for(int i = 0; i < Ndim; i++)
    for(int j = 0; j < Ndim; j++) 
      for(int k = 0; k < Ndim; k++)
	for(int l = 0; l < Ndim; l++) 
	  cc[i][j][k][l] = s1[i][j]*s2[k][l];
      
}

template<const int Ndim>
inline fourTensor<Ndim> ::fourTensor(double mu, double lambda)
{
  smatrix<Ndim> delta;
  delta.identity();
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++) {
      for(int k = 0;k<Ndim;k++){
	for(int l = 0;l<Ndim;l++){
	  cc[i][j][k][l] = lambda*delta[i][j]*delta[k][l]+mu*(delta[i][k]*delta[j][l]+ delta[i][l]*delta[j][k]);
	}}}}
}


template<const int Ndim>
inline fourTensor<Ndim> ::fourTensor(double c11, double c12, double c44)
{
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++) {
      for(int k = 0;k<Ndim;k++){
	for(int l = 0;l<Ndim;l++){
	  cc[i][j][k][l] = 0;

	  if(i == 0 && j == 0 && k == 0 && l == 0 )
	    cc[i][j][k][l] = c11;
	  if(i == 1 && j == 1 && k == 1 && l == 1 )
	    cc[i][j][k][l] = c11;
	  if(i == 2 && j == 2 && k == 2 && l == 2 )
	    cc[i][j][k][l] = c11;

	  if(i == 0 && j == 0 && k == 1 && l == 1 )
	    cc[i][j][k][l] = c12;
	  if(i == 1 && j == 1 && k == 0 && l == 0 )
	    cc[i][j][k][l] = c12;
	  if(i == 0 && j == 0 && k == 2 && l == 2 )
	    cc[i][j][k][l] = c12;
	  if(i == 2 && j == 2 && k == 0 && l == 0 )
	    cc[i][j][k][l] = c12;
	  if(i == 1 && j == 1 && k == 2 && l == 2 )
	    cc[i][j][k][l] = c12;
	  if(i == 2 && j == 2 && k == 1 && l == 1 )
	    cc[i][j][k][l] = c12;

	  if(i == 1 && j == 2 && k == 1 && l == 2 )
	    cc[i][j][k][l] = c44;
	  if(i == 1 && j == 2 && k == 2 && l == 1 )
	    cc[i][j][k][l] = c44;
	  if(i == 2 && j == 1 && k == 1 && l == 2 )
	    cc[i][j][k][l] = c44;
	  if(i == 2 && j == 1 && k == 2 && l == 1 )
	    cc[i][j][k][l] = c44;
	 

	  if(i == 0 && j == 2 && k == 0 && l == 2 )
	    cc[i][j][k][l] = c44;
	  if(i == 0 && j == 2 && k == 2 && l == 0 )
	    cc[i][j][k][l] = c44;
	  if(i == 2 && j == 0 && k == 0 && l == 2 )
	    cc[i][j][k][l] = c44;
	  if(i == 2 && j == 0 && k == 2 && l == 0 )
	    cc[i][j][k][l] = c44;

	  if(i == 0 && j == 1 && k == 0 && l == 1 )
	    cc[i][j][k][l] = c44;
	  if(i == 0 && j == 1 && k == 1 && l == 0 )
	    cc[i][j][k][l] = c44;
	  if(i == 1 && j == 0 && k == 0 && l == 1 )
	    cc[i][j][k][l] = c44;
	  if(i == 1 && j == 0 && k == 1 && l == 0 )
	    cc[i][j][k][l] = c44;

	  
	}}}}
}

template<const int Ndim> 
inline fourTensor<Ndim>& fourTensor<Ndim>::zero()
{
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0;k<Ndim;k++){
	for(int l = 0;l<Ndim;l++){
	  cc[i][j][k][l] = 0;
	}}}}
  return(*this);
}

template<const int Ndim>
inline fourTensor<Ndim> fourTensor<Ndim>::identity()
{
  fourTensor<Ndim> id;

  id.zero();

  for(int i = 0; i < Ndim; i++)
    for(int j = 0; j < Ndim; j++)
      for(int k = 0; k < Ndim; k++)
	for(int l = 0; l < Ndim; l++)
	  if((i == k) && (j == l)) id.cc[i][j][k][l] = 1;

  return id;
}


// template<const int Ndim>
// inline fourTensor<Ndim>& fourTensor<Ndim>::identity()
// {
//   smatrix<Ndim> delta;
//   delta.identity();
//   for(int i = 0; i < Ndim; i++){
//     for(int j = 0; j < Ndim; j++){ 
//       for(int k = 0; k < Ndim; k++){
// 	for(int l = 0; l < Ndim; l++){
// 	  {
// 	    cc[i][j][k][l] = delta[i][j]*delta[k][l];
// 	  }}}}}
//   return (*this);
// }

template<const int Ndim>
inline fourTensor<Ndim> fourTensor<Ndim>::operator+(const fourTensor<Ndim>& f) const
{
  fourTensor<Ndim> cpf;
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++){
	  double scal =cc[i][j][k][l]+f(i,j,k,l);
	  cpf.setValue(i,j,k,l,scal);
	}}}}
  return cpf;
}

template<const int Ndim>
inline fourTensor<Ndim>& fourTensor<Ndim>::operator+=(const fourTensor<Ndim>& f)
{
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++){
	  double scal =cc[i][j][k][l]+f(i,j,k,l);
	  cc[i][j][k][l] = scal;
	}}}}
  return(*this);
}

template<const int Ndim> 
inline fourTensor<Ndim> fourTensor<Ndim>::operator-(const fourTensor<Ndim>& f) const
{
  fourTensor<Ndim> cmf;
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++){
	  double scal =cc[i][j][k][l]-f(i,j,k,l);
	  cmf.setValue(i,j,k,l,scal);
	}}}}
  return cmf;
}

template<const int Ndim>
inline fourTensor<Ndim> fourTensor<Ndim>::operator-() const
{
  fourTensor<Ndim> mc;
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++){
	  double scal =-1.0*cc[i][j][k][l];
	  mc.setValue(i,j,k,l,scal);
	}}}}
  return mc;
}

template<const int Ndim>
inline fourTensor<Ndim>& fourTensor<Ndim>::operator=(const fourTensor<Ndim>& f)
{
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++){
	  cc[i][j][k][l] = f(i,j,k,l);
	}}}}
  return (*this);
}

template<const int Ndim>
inline double fourTensor<Ndim>::operator()(int i, int j, int k, int l) const
{
  assert((0 <= i) && (i < Ndim));
  assert((0 <= j) && (j < Ndim));
  assert((0 <= k) && (k < Ndim));
  assert((0 <= l) && (l < Ndim));

  return( cc[i][j][k][l] );
}



template<const int Ndim>
inline double& fourTensor<Ndim>::operator()(int i, int j, int k, int l)
{
  assert((0 <= i) && (i < Ndim));
  assert((0 <= j) && (j < Ndim));
  assert((0 <= k) && (k < Ndim));
  assert((0 <= l) && (l < Ndim));

  return( cc[i][j][k][l] );
}


template<const int Ndim>
inline smatrix<Ndim> fourTensor<Ndim>::operator*(const smatrix<Ndim>& s) const
{
  double arg[Ndim][Ndim];
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      double scal = 0;
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++){
	  scal += cc[i][j][k][l]*s[k][l] ;
	}}
      //std::cout<<"scal = "<<scal<<std::endl;
      arg[i][j] = scal;
    }}
  smatrix<Ndim> cts(arg);
  return(cts);
}

template<const int Ndim>
inline smatrix<Ndim> fourTensor<Ndim>::contract(const smatrix<Ndim>& s, int idx1, int idx2) const
{
  if(idx1 == 2 && idx2 == 3)
    {
      double arg[Ndim][Ndim];
      for(int i = 0; i < Ndim; i++){
	for(int j = 0; j < Ndim; j++){ 
	  double scal = 0;
	  for(int k = 0; k < Ndim; k++){
	    for(int l = 0; l < Ndim; l++){
	      scal += cc[i][j][k][l]*s[k][l] ;
	    }}
	  // std::cout<<"scal = "<<scal<<std::endl;
	  arg[i][j] = scal;
	}}
      smatrix<Ndim> cts(arg);
      return(cts);
    }
  else if(idx1 == 1 && idx2 == 3)
    {   
      double arg[Ndim][Ndim];
      for(int i = 0; i < Ndim; i++){
	for(int k = 0; k < Ndim; k++){
	  double scal = 0;
	  for(int j = 0; j < Ndim; j++){ 
	    for(int l = 0; l < Ndim; l++){
	      scal += cc[i][j][k][l]*s[j][l] ;
	    }}
	  // std::cout<<"scal = "<<scal<<std::endl;
	  arg[i][k] = scal;
	}}
      smatrix<Ndim> cts(arg);
      return(cts);
    }
  else
    throw("not valid contract indexes, to be implemented.");
}


template<const int Ndim>
inline void fourTensor<Ndim>::setValue(int i, int j, int k, int l, double val)
{
  if(i>=Ndim || j>=Ndim || k>=Ndim || l>=Ndim) 
    throw("fourTensor->setvalue(),index out of range!");
  if(i<0||j<0||k<0||l<0)
    throw("fourTensor->setvalue(),index out of range!");
  cc[i][j][k][l] = val;
}  


  template<const int Ndim>
  inline fourTensor<Ndim>
  fourTensor<Ndim>::compose(const fourTensor<Ndim>& dd)
  {
    fourTensor<Ndim> ee;

    ee.zero();

    for(int i = 0; i < Ndim; i++)
      for(int j = 0; j < Ndim; j++)
	for(int k = 0; k < Ndim; k++)
	  for(int l = 0; l < Ndim; l++) 
	    for(int m = 0;m<Ndim;m++)
	      for(int n = 0;n<Ndim;n++)
		ee.cc[i][j][k][l] += cc[i][j][m][n] * dd.cc[m][n][k][l];

    return ee;
  }


// utility functions
/*
  template<const int Ndim> 
  inline fourTensor<Ndim> operator*(const double a, const fourTensor<Ndim>& v)
  {
  fourTensor<Ndim> w;

  for(int i = 0; i < Ndim; i++)
  for(int j = 0; j < Ndim; j++)
  for(int k = 0; k < Ndim; k++)
  for(int l = 0; l < Ndim; l++) w(i,j,k,l) = a * v(i,j,k,l);

  return(w);
  }
*/

template<const int Ndim>
std::ostream& operator<<(std::ostream& out, fourTensor<Ndim> cc)
{
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){
      out<<std::endl;
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++)
	  out<<"F["<<i<<"]"<<"["<<j<<"]"<<"["<<k<<"]"<<"["<<l<<"] = "<<setw(6)<<cc(i,j,k,l)<<";";
      }
    }
  }
  return(out);
}




template<const int Ndim>
inline fourTensor<Ndim> operator*(const double a, const fourTensor<Ndim>& f)
{
  fourTensor<Ndim> cpf;
  for(int i = 0; i < Ndim; i++){
    for(int j = 0; j < Ndim; j++){ 
      for(int k = 0; k < Ndim; k++){
	for(int l = 0; l < Ndim; l++){
	  cpf.setValue(i,j,k,l, a*f(i,j,k,l));
	}}}}
  return cpf;
}

#endif                      // !defined (__FOURTENSOR_C__)
