#include<map>
#include<utility>

// Discontinuous basis functions are similar to Lagrange basis functions
// but don't have (well defined) boundary basis functions; all of their
// dof's are interior to the cell. *boundary should point to an empty 
// boundary (no dof's) to avoid screwing up the boundary FE itterator.

// For DG methods it may be useful to include a function to
// compute the trace on the boundary:
//
// void (*boundaryBasis)(svector<Ndim-1>, double *phi, svector<Ndim> *dphi);
	
template<const int Ndim>
struct DiscontinuousBasisFunction : public basisFunction<Ndim>
{
  DiscontinuousBasisFunction(basisFunction<Ndim> &bf, 
	       void (*bb)(svector<Ndim>, double *phi, svector<Ndim> *dphi))
    : basisFunction<Ndim>(bf), basis(bb)
  {return;}

  void (*basis)(svector<Ndim>, double *phi, svector<Ndim> *dphi);
};

// **************** discontinous finite elements ***************
// **************** discontinous finite elements ***************
// **************** discontinous finite elements ***************

void tri1(svector<2> xi, double psi[1], svector<2> dpsi[1]);

void discontinuousTriangle3Permutation(const int *nodes, int *dof);
void discontinuousTet4Permutation (const int *nodes, int *dof);
void discontinuousTet10Permutation(const int *nodes, int *dof);

void tet1(svector<3> xi, double phi[1], svector<3> dphi[1]);

// ******************* Interval ****************************
// ******************* Interval ****************************
// ******************* Interval ****************************

svector<1> Line1Pts[] = {svector<1>(0.0)};

basisFunction<1> Line1BF =   // one-d p/w constant basis fn
  {1,
   {0,1},
   EmptyVertex,
   &Interval,
   1,
   Line1Pts,
   identityPermutation
  };

DiscontinuousBasisFunction<1> Line1class(Line1BF, line1),
  *Line1 = &Line1class;


// ******************* Triangles ****************************
// ******************* Triangles ****************************
// ******************* Triangles ****************************

svector<2> DiscontinuousTriangle1Pts[] = {svector<2>(1.0/3.0, 1.0/3.0)};

basisFunction<2> DiscontinuousTriangle1BF = // piecewise constant triangle
  {1,
   {0,0,1},
   EmptyInterval,
   &Triangle,
   1,
   DiscontinuousTriangle1Pts,
   identityPermutation
  };

DiscontinuousBasisFunction<2> 
  DiscontinuousTriangle1Class(DiscontinuousTriangle1BF, tri1),
  *DiscontinuousTriangle1 = &DiscontinuousTriangle1Class;

// use the same interpolation as Lagrange triangle
// so that we can use the same basis functions 

basisFunction<2> DiscontinuousTriangle3BF =  // piecewise linear triangle
  {3,
   {0,0,3},
   EmptyInterval,
   &Triangle,
   3,
   Triangle6Pts, 
   discontinuousTriangle3Permutation
  };

DiscontinuousBasisFunction<2> 
  DiscontinuousTriangle3Class(DiscontinuousTriangle3BF, tri3),
 *DiscontinuousTriangle3 = &DiscontinuousTriangle3Class;

// ******************* Squares ****************************
// ******************* Squares ****************************
// ******************* Squares ****************************

svector<2> DiscontinuousSquare1Pts[] = {svector<2>(0.0, 0.0)};

basisFunction<2> DiscontinuousSquare1BF = // piecewise constant square
  {1,
   {0,0,1},
   EmptyInterval,
   &Square,
   1,
   DiscontinuousSquare1Pts,
   identityPermutation
  };

// tri1() is the constant function phi(x) = 1 for x in 2d

DiscontinuousBasisFunction<2> 
  DiscontinuousSquare1Class(DiscontinuousSquare1BF, tri1),
 *DiscontinuousSquare1 = &DiscontinuousSquare1Class;

// ******************* Tetrahedra ****************************
// ******************* Tetrahedra ****************************
// ******************* Tetrahedra ****************************

svector<3> DiscontinuousTet1Pts[] ={svector<3>(1.0/4.0, 1.0/4.0, 1.0/4.0)};

basisFunction<3> DiscontinuousTet1BF = // piecewise constant tet
  {1,
   {0,0,0,1},
   EmptyTriangle,
   &Tetrahedron,
   1,
   DiscontinuousTet1Pts,
   identityPermutation
  };

DiscontinuousBasisFunction<3> 
  DiscontinuousTet1Class(DiscontinuousTet1BF, tet1),
 *DiscontinuousTet1 = &DiscontinuousTet1Class;

basisFunction<3> DiscontinuousTet4BF = // piecewise linear tet
  {4,
   {0,0,0,4},
   EmptyTriangle,
   &Tetrahedron,
   4,
   Tet4BubblePts,
   discontinuousTet4Permutation 
  };

DiscontinuousBasisFunction<3> 
  DiscontinuousTet4Class(DiscontinuousTet4BF, tet4),
 *DiscontinuousTet4 = &DiscontinuousTet4Class;

basisFunction<3> DiscontinuousTet10BF = // piecewise quadratic tet
  {10,
   {0,0,0,10},
   EmptyTriangle,
   &Tetrahedron,
   10,
   Tet10Pts,
   discontinuousTet10Permutation 
  };

DiscontinuousBasisFunction<3> 
  DiscontinuousTet10Class(DiscontinuousTet10BF, tet10),
 *DiscontinuousTet10 = &DiscontinuousTet10Class;

// ******************* Cubes ****************************
// ******************* Cubes ****************************
// ******************* Cubes ****************************

svector<3> DiscontinuousCube1Pts[] ={svector<3>(0.0, 0.0, 0.0)};

basisFunction<3> DiscontinuousCube1BF = // piecewise constant tet
  {1,
   {0,0,0,1},
   EmptySquare,
   &Cube,
   1,
   DiscontinuousCube1Pts,
   identityPermutation
  };

// tet1() is the constant function phi(x) = 1 for x in 3d

DiscontinuousBasisFunction<3> 
  DiscontinuousCube1Class(DiscontinuousCube1BF, tet1),
 *DiscontinuousCube1 = &DiscontinuousCube1Class;

// ************** Basis Functions and Permutaionts *******************
// ************** Basis Functions and Permutaionts *******************
// ************** Basis Functions and Permutaionts *******************

void tri1(svector<2> xi, double phi[1], svector<2> dphi[1])
{
  phi[0] = 1.0;

 if(dphi == NULL) return;

 dphi[0] = svector<2>(0,0);
 
 return;
}

void discontinuousTriangle3Permutation(const int nodes[3], int dof[3])
{
  int nn[3];

  for(int i = 0; i < 3; i++) nn[i] = nodes[i];

  sortDof(3, nn, dof);

  return;
}

void tet1(svector<3> xi, double phi[1], svector<3> dphi[1])
{
  phi[0] = 1.0;
  
  if(dphi == NULL) return;
  
  dphi[0] = svector<3>(0,0,0);
  
  return;
}

void discontinuousTet4Permutation(const int nodes[4], int dof[4])
{
  int nn[4];

  for(int i = 0; i < 4; i++) nn[i] = nodes[i];

  sortDof(4, nn, dof);

  return;
}

void discontinuousTet10Permutation(const int nodes[10], int dof[10])
{
  int nn[4];

  for(int i = 0; i < 4; i++) nn[i] = nodes[i];

  sortDof(4, nn, dof);        // sorts vertex dof

  std::map< std::pair<int,int>, int > eDof;      // store edge dof

  for(int i = 0; i < 6; i++) 
    {
    std::pair<int, int> 
	e1(nn[Tetrahedron.subs(1,i,0)], nn[Tetrahedron.subs(1,i,1)]),
	e2(nn[Tetrahedron.subs(1,i,1)], nn[Tetrahedron.subs(1,i,0)]);


      eDof[e1] = eDof[e2] = dof[4+i];   // dof at mid-edge
    }

  for(int i = 0; i < 6; i++) 
    {
      std::pair<int, int> ee(nodes[Tetrahedron.subs(1,i,0)], 
			nodes[Tetrahedron.subs(1,i,1)]);

      assert(eDof.find(ee) != eDof.end());

      dof[4+i] = eDof[ee];
    }

  return;
}
