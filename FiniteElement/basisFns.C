// Base class for a basis function

#ifndef __BASISFNS_C__
#define __BASISFNS_C__

#include<cassert>

#include "../MeshGeneration/masterCell.C"
#include "../LinearAlgebra/svector.C"
#include "../LinearAlgebra/smatrix.C"
#include "../LinearAlgebra/sbmatrix.C"

template<const int Ndim>
struct basisFunction 
{
  const int ndof;           // total number of degrees of freedom
  int subDof[Ndim+1];       // number of degrees of freedom on sub-cells

  basisFunction<Ndim-1> *boundary; // boundary basis function  

  masterCell<Ndim> *master;

  const int npts;           // number of intrpolation points

  svector<Ndim> *pts;       // local coordinates of interpolation points

  void (*permuteDof)(const int *nodes, int *dof);
};

// ******** empty classes used when no boundary values exist ***********
// ******** empty classes used when no boundary values exist ***********
// ******** empty classes used when no boundary values exist ***********

void identityPermutation(const int *nodes, int *dof) { return; }

basisFunction<0> EmptyVertexClass =   
  {0,
   {0},
   NULL,
   NULL,
   0,
   NULL,
   identityPermutation
  }, *EmptyVertex = &EmptyVertexClass;

basisFunction<1> EmptyIntervalClass =   
  {0,
   {0,0},
   EmptyVertex,
   &Interval,
   0,
   NULL,
   identityPermutation
  }, *EmptyInterval = &EmptyIntervalClass;

basisFunction<2> EmptyTriangleClass =   // empty triangle
  {0,
   {0,0,0},
   EmptyInterval,
   &Triangle,
   0,
   NULL,
   identityPermutation
  }, *EmptyTriangle = &EmptyTriangleClass;

basisFunction<2> EmptySquareClass =  // empty square 
  {0,
   {0,0,0},
   EmptyInterval,
   &Square,
   0,
   NULL,
   identityPermutation
  }, *EmptySquare = &EmptySquareClass;

// ************************* classical families of elements **************
// ************************* classical families of elements **************
// ************************* classical families of elements **************

#include"lagrangeBasisFns.h"
#include"discontinuousFns.C"   // currently uses lagrangeBasisFns <----
// #include"hdivBasisFns.C"
// #include"awBasisFns.C"
// #include"hermiteBasisFns.h"

#endif                          // !defined(__BASISFNS_C__)
