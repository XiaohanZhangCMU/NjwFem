#ifndef __GELEMENT_C__
#define __GELEMENT_C__

#include"masterCell.C"
#include"../LinearAlgebra/svector.C"

template<const int Ndim, const int Nchi>
class GeometricElement
{
public:
  GeometricElement() : type(0), master(NULL), basis(NULL) {return;}

  int type;

  int nodes[Nchi];
  
  int tag;

  masterCell<Ndim> *master;

  svector<Ndim> map(svector<Ndim> xi);
  svector<Ndim> map(svector<Ndim> xi, smatrix<Ndim> &dxds);

  svector<Ndim> xx[Nchi];
  void (*basis)(svector<Ndim>, double phi[Nchi], svector<Ndim> dphi[Nchi]);
};

template<const int Ndim, const int Nchi>
class BoundaryGelement
{
public:
  BoundaryGelement() : type(0), master(NULL), basis(NULL) {return;}

  int type;

  int nodes[Nchi];

  int tag;

  void print(std::ostream &);

  masterCell<Ndim-1> *master;

  svector<Ndim> map(svector<Ndim-1> xi);
  svector<Ndim> map(svector<Ndim-1> xi, sbmatrix<Ndim> &dxds);

  svector<Ndim> xx[Nchi];
  void (*basis)(svector<Ndim-1>, double *phi, svector<Ndim-1> *dphi);
};


template<const int Ndim, const int Nchi>
svector<Ndim> GeometricElement<Ndim,Nchi>::map(svector<Ndim> xi)
{
  double psi[Nchi];
  svector<Ndim> pt;

  basis(xi, psi, NULL);

  pt.zero();

  for(int i = 0; i < Nchi; i++) pt += psi[i] * xx[i];

  return(pt);
}

template<const int Ndim, const int Nchi>
svector<Ndim> GeometricElement<Ndim,Nchi>::map(svector<Ndim> xi, smatrix<Ndim> &dxds)
{
  double psi[Nchi];
  svector<Ndim> pt, dpsi[Nchi];

  basis(xi, psi, dpsi);

  pt.zero();
  dxds.zero();

  for(int i = 0; i < Nchi; i++) 
    {
      pt   += psi[i] * xx[i];
      dxds += smatrix<Ndim>(xx[i], dpsi[i]);   // outer product
    }

  return(pt);
}

template<const int Ndim, const int Nchi>
svector<Ndim> BoundaryGelement<Ndim,Nchi>::map(svector<Ndim-1> xi)
{
  double psi[Nchi];
  svector<Ndim> pt;

  basis(xi, psi, NULL);

  pt.zero();

  for(int i = 0; i < Nchi; i++) pt += psi[i] * xx[i];

  return(pt);
}

template<const int Ndim, const int Nchi>
svector<Ndim> BoundaryGelement<Ndim,Nchi>::map(svector<Ndim-1> xi, 
					       sbmatrix<Ndim> &dxds)
{
  double psi[Nchi];
  svector<Ndim> pt;
  svector<Ndim-1> dpsi[Nchi];

  basis(xi, psi, dpsi);

  pt.zero();
  dxds.zero();

  for(int i = 0; i < Nchi; i++) 
    {
      pt += psi[i] * xx[i];
      dxds += sbmatrix<Ndim>(xx[i], dpsi[i]);
    }

  return(pt);
}
#endif                          // !defined(__GELEMENT_C__)
