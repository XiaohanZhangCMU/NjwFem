// Classical Lagrange basis functions are scalar valued and all
// derivatives are defined.  The boundary is a Lagrange basis
// function of one dimension less.

#ifndef __LAGRANGEBASISFNS_H__
#define __LAGRANGEBASISFNS_H__

template<const int Ndim>
struct LagrangeBasisFunction : public basisFunction<Ndim>
{
  LagrangeBasisFunction(basisFunction<Ndim> &bf, 
	       LagrangeBasisFunction<Ndim-1> *bdy,
	       void (*bb)(svector<Ndim>, double *phi, svector<Ndim> *dphi))
    : basisFunction<Ndim>(bf), boundary(bdy), basis(bb)
  {return;}

  LagrangeBasisFunction<Ndim-1> *boundary;
  void (*basis)(svector<Ndim>, double *phi, svector<Ndim> *dphi);
};

template<const int Ndim>
struct LagrangeBasisFunctionBar : public basisFunction<Ndim>
{
  LagrangeBasisFunctionBar(basisFunction<Ndim> &bf, 
   LagrangeBasisFunction<Ndim-1> *bdy,
   void (*bb)(svector<Ndim>, double *phi, svector<Ndim> *dphi, double *phiBar))
    : basisFunction<Ndim>(bf), boundary(bdy), basis(bb)
  {return;}

  LagrangeBasisFunction<Ndim-1> *boundary;
  void (*basis)(svector<Ndim>, double *phi, 
		svector<Ndim> *dphi, double *phiBar);
};

// **************** classical Lagrange finite elements ***************
// **************** classical Lagrange finite elements ***************
// **************** classical Lagrange finite elements ***************

void vertex1(svector<0> eta, double psi[1], svector<0> dpsi[1])
{
  psi[0] = 1;

  return;
}

void line1(svector<1> eta, double psi[1], svector<1> dpsi[1]);
void line2(svector<1> eta, double psi[2], svector<1> dpsi[2]);
void line3(svector<1> eta, double psi[3], svector<1> dpsi[3]);
void line4(svector<1>  xi, double psi[4], svector<1> dpsi[4]);
void line5(svector<1>  xi, double psi[5], svector<1> dpsi[5]);
void line6(svector<1>  xi, double psi[6], svector<1> dpsi[6]);

// projections of lineX() basis functions onto polynomials
// of one degree less. The last argument should be NULL

void line2bar(svector<1> peta, double psi[2], 
	      svector<1> dpsi[2], double psiBar[2]);
void line3bar(svector<1> peta, double psi[3], 
	      svector<1> dpsi[3], double psiBar[3]);

void tri3(svector<2> xi, double psi[3], svector<2> dpsi[3]);
void tri6(svector<2> xi, double psi[6], svector<2> dpsi[6]);
void tri10(svector<2> pt, double psi[10], svector<2> dpsi[10]);

void tet4( const svector<3> xi, double phi[4] , svector<3> dphi[4]);
void tet10(const svector<3> xi, double phi[10], svector<3> dphi[10]);
void tet20(const svector<3> pt, double psi[20], svector<3> dpsi[20]);  

void quad4(svector<2> xi, double psi[4], svector<2> dsi[4]);
void quad8(svector<2> xi, double phi[8], svector<2> dphi[8]);
void quad9(svector<2> xi, double phi[9], svector<2> dphi[9]);
void quad12(svector<2> xi, double phi[12], svector<2> dphi[12]);
void quad16(svector<2> xi, double phi[16], svector<2> dphi[16]);
void quad17(svector<2> xi, double phi[17], svector<2> dphi[17]);
void quad24(svector<2> xi, double phi[24], svector<2> dphi[24]);

void cube8(svector<3>, double psi[8], svector<3> dpsi[8]);
void cube20(svector<3>, double psi[20], svector<3> dpsi[20]);
void cube27(svector<3>, double psi[27], svector<3> dpsi[27]);
void cube32(svector<3>, double psi[32], svector<3> dpsi[32]);
void cube50(svector<3>, double psi[50], svector<3> dpsi[50]);

void tri3Bubble(svector<2> xi, double psi[4], svector<2> dpsi[4]);
void tet4Bubble(svector<3> xi, double phi[5] , svector<3> dphi[5]);
void quad4Bubble(svector<2> xi, double psi[5], svector<2> dsi[5]);
void cube8Bubble(svector<3>, double psi[9], svector<3> dpsi[9]);

void tri3BubbleInterp(svector<2> xi, double  sn[4], svector<2>  dsn[4]);
void tet4BubbleInterp(svector<3> xi, double phi[5], svector<3> dphi[5]);

void sortDof(int n, int *nn, int *dd);               // in utility.C
void line4Permutation(const int *nodes, int *dof);
void line5Permutation(const int *nodes, int *dof);
void line6Permutation(const int *nodes, int *dof);
void tri10Permutation(const int *nodes, int *dof);
void tet20Permutation(const int *nodes, int *dof);
void quad12Permutation(const int *nodes, int *dof);
void quad16Permutation(const int *nodes, int *dof);
void quad17Permutation(const int *nodes, int *dof);
void quad24Permutation(const int *nodes, int *dof);
void cube32Permutation(const int *nodes, int *dof);
void cube50Permutation(const int *nodes, int *dof);

// ********** Vertex terminates the boundary reccursion *********
// ********** Vertex terminates the boundary reccursion *********
// ********** Vertex terminates the boundary reccursion *********

svector<0> Vertex1Pts[] = {svector<0>()};

basisFunction<0> Vertex1BF =    // ends the boundary reccursion
  {1,
   {1},
   NULL,
   &Point,
   1,
   Vertex1Pts,
   identityPermutation,
  };

LagrangeBasisFunction<0> Vertex1class(Vertex1BF, NULL, vertex1),
  *Vertex1 = &Vertex1class;

// ***************** 1d Functions on Interval **********************
// ***************** 1d Functions on Interval **********************
// ***************** 1d Functions on Interval **********************

svector<1> Line3Pts[] = {svector<1>(-1.0), svector<1>(1.0), svector<1>(0.0)};
 
basisFunction<1> Line2BF =  // one-d linear element basis fns
  {2,
   {1,0},
   Vertex1,
   &Interval,
   2,
   Line3Pts,
   identityPermutation
  };

LagrangeBasisFunction<1> Line2class(Line2BF, Vertex1, line2),
  *Line2 = &Line2class;

basisFunction<1> Line2barBF =   // projections of Line2 onto P_0
  {2,
   {1,0},
   Vertex1,
   &Interval,
   2,
   Line3Pts,
   identityPermutation,
  };

LagrangeBasisFunctionBar<1> Line2barClass(Line2barBF, Vertex1, line2bar),
  *Line2bar = &Line2barClass;

basisFunction<1> Line3BF =    // one-d quadratic basis fns
  {3,
   {1,1},
   Vertex1,
   &Interval,
   3,
   Line3Pts,
   identityPermutation
  };

LagrangeBasisFunction<1> Line3class(Line3BF, Vertex1, line3),
  *Line3 = &Line3class;

basisFunction<1> Line3barBF =   // Project one-d quadratics onto linears
  {3,
   {1,1},
   Vertex1,
   &Interval,
   3,
   Line3Pts,
   identityPermutation
  };

LagrangeBasisFunctionBar<1> Line3barClass(Line3barBF, Vertex1, line3bar),
  *Line3bar = &Line3barClass;

// svector<1> Line4Pts[] = {svector<1>(-1.0), svector<1>(1.0), 
// 			 svector<1>(-1.0/3.0), svector<1>(1.0/3.0)};

svector<1> Line4Pts[] = {svector<1>(-1.0), svector<1>(1.0), 
	       svector<1>(-1.0/sqrt(5.0)), svector<1>(1.0/sqrt(5.0))};

basisFunction<1> Line4BF =    // one-d cubic basis fns
  {4,
   {1,2},
   Vertex1,
   &Interval,
   4,
   Line4Pts,
   line4Permutation
  };

LagrangeBasisFunction<1> Line4class(Line4BF, Vertex1, line4),
  *Line4 = &Line4class;

svector<1> Line5Pts[] = {svector<1>(-1.0), svector<1>(1.0), 
   svector<1>(-sqrt(3.0/7.0)), svector<1>(0.0), svector<1>(sqrt(3.0/7.0))};

basisFunction<1> Line5BF =    // one-d cubic basis fns
  {5,
   {1,3},
   Vertex1,
   &Interval,
   5,
   Line5Pts,
   line5Permutation
  };

LagrangeBasisFunction<1> Line5class(Line5BF, Vertex1, line5),
  *Line5 = &Line5class;

svector<1> Line6Pts[] = 
  {svector<1>(-1.0), svector<1>(1.0), 
   svector<1>(-sqrt((3+2*sqrt(6.0/5))/7.0)),
   svector<1>(-sqrt((3-2*sqrt(6.0/5))/7.0)), 
   svector<1>( sqrt((3-2*sqrt(6.0/5))/7.0)),
   svector<1>( sqrt((3+2*sqrt(6.0/5))/7.0))};

basisFunction<1> Line6BF =    // one-d cubic basis fns
  {6,
   {1,4},
   Vertex1,
   &Interval,
   6,
   Line6Pts,
   line6Permutation
  };

LagrangeBasisFunction<1> Line6class(Line6BF, Vertex1, line6),
  *Line6 = &Line6class;

// ***************** 2d Functions on Triangle **********************
// ***************** 2d Functions on Triangle **********************
// ***************** 2d Functions on Triangle **********************

svector<2> Triangle6Pts[] = 
  {svector<2>(0.0, 0.0), svector<2>(1.0, 0.0), svector<2>(0.0, 1.0), 
   svector<2>(0.5, 0.0), svector<2>(0.5, 0.5), svector<2>(0.0, 0.5)};

basisFunction<2> Triangle3BF =    // two-d linear basis fns
  {3,
   {1,0,0},
   Line2,
   &Triangle,
   3,
   Triangle6Pts,
   identityPermutation
  };

LagrangeBasisFunction<2> Triangle3class(Triangle3BF, Line2, tri3),
 *Triangle3 = &Triangle3class;

svector<2> Triangle3BubblePts[] = 
  {svector<2>(0.0, 0.0), svector<2>(1.0, 0.0), svector<2>(0.0, 1.0), 
   svector<2>(1.0/3.0, 1.0/3.0)};

basisFunction<2> Triangle3BubbleBF =    // two-d linear basis fns
  {4,
   {1,0,1},
   Line2,
   &Triangle,
   3,                        // not clear how to initialize bubble
   Triangle3BubblePts,
   identityPermutation,
  };

LagrangeBasisFunction<2> 
  Triangle3BubbleClass(Triangle3BubbleBF, Line2, tri3Bubble),
  *Triangle3Bubble = &Triangle3BubbleClass;

// Alternative, interpolate at all 4 points and use tri3BubbleInterp()

basisFunction<2> Triangle3BubbleInterpBF =
  {4,
   {1,0,1},
   Line2,
   &Triangle,
   4,                        // not clear how to initialize bubble
   Triangle3BubblePts,
   identityPermutation,
  };

LagrangeBasisFunction<2> 
  Triangle3BubbleInterpClass(Triangle3BubbleInterpBF, Line2, tri3BubbleInterp),
  *Triangle3BubbleInterp = &Triangle3BubbleInterpClass;

basisFunction<2> Triangle6BF =    // two-d quadratic basis fns
  {6,
   {1,1,0},
   Line3,
   &Triangle,
   6,
   Triangle6Pts,
   identityPermutation
  };

LagrangeBasisFunction<2> Triangle6class(Triangle6BF, Line3, tri6),
  *Triangle6 = &Triangle6class;

svector<2> Triangle10Pts[10] = 
  {svector<2>(0.0, 0.0), svector<2>(1.0, 0.0), svector<2>(0.0, 1.0), 
   svector<2>(0.5*(1-1.0/sqrt(5.0)), 0.0), 
   svector<2>(0.5*(1+1.0/sqrt(5.0)), 0.0),
   svector<2>(0.5*(1+1.0/sqrt(5.0)), 0.5*(1-1.0/sqrt(5.0))), 
   svector<2>(0.5*(1-1.0/sqrt(5.0)), 0.5*(1+1.0/sqrt(5.0))),
   svector<2>(0.0, 0.5*(1+1.0/sqrt(5.0))), 
   svector<2>(0.0, 0.5*(1-1.0/sqrt(5.0))),
   svector<2>(1.0/3.0, 1.0/3.0)};

basisFunction<2> Triangle10BF =    // two-d quadratic basis fns
  {10,
   {1,2,1},
   Line4,
   &Triangle,
   10,
   Triangle10Pts,
   tri10Permutation
  };

LagrangeBasisFunction<2> Triangle10class(Triangle10BF, Line4, tri10),
  *Triangle10 = &Triangle10class;

// ***************** 2d Functions on Square **********************
// ***************** 2d Functions on Square **********************
// ***************** 2d Functions on Square **********************

svector<2> Square9Pts[] = {svector<2>(-1.0,-1.0), svector<2>( 1.0,-1.0), 
			 svector<2>( 1.0, 1.0), svector<2>(-1.0, 1.0),
			 svector<2>( 0.0,-1.0), svector<2>( 1.0, 0.0),
			 svector<2>( 0.0, 1.0), svector<2>(-1.0, 0.0),
			 svector<2>( 0.0, 0.0)};

basisFunction<2> Square4BF =    // two-d bilinear basis fns
  {4,
   {1,0,0},
   Line2,
   &Square,
   4,
   Square9Pts,
   identityPermutation
  };

LagrangeBasisFunction<2> Square4class(Square4BF, Line2, quad4),
  *Square4 = &Square4class;

svector<2> Square4BubblePts[] = 
  {svector<2>(-1.0,-1.0), svector<2>( 1.0,-1.0), 
   svector<2>( 1.0, 1.0), svector<2>(-1.0, 1.0),
   svector<2>( 0.0, 0.0)};

basisFunction<2> Square4BubbleBF =    // two-d bilinear basis fns
  {5,
   {1,0,1},
   Line2,
   &Square,
   4,                           // not clear how to initialize bubble
   Square4BubblePts,
   identityPermutation
  };

LagrangeBasisFunction<2> 
  Square4BubbleClass(Square4BubbleBF, Line2, quad4Bubble),
  *Square4Bubble = &Square4BubbleClass;

basisFunction<2> Square8BF =    // serindipity two-d biquadratic basis fns
  {8,
   {1,1,0},
   Line3,
   &Square,
   8,
   Square9Pts,
   identityPermutation
  };

LagrangeBasisFunction<2> Square8class(Square8BF, Line3, quad8),
  *Square8 = &Square8class;

basisFunction<2> Square9BF =    // two-d biquadratic basis fns
  {9,
   {1,1,1},
   Line3,
   &Square,
   9,
   Square9Pts,
   identityPermutation
  };

LagrangeBasisFunction<2> Square9class(Square9BF, Line3, quad9),
 *Square9 = &Square9class;

svector<2> Square16Pts[] = {svector<2>(-1.0,-1.0), svector<2>( 1.0,-1.0), 
			    svector<2>( 1.0, 1.0), svector<2>(-1.0, 1.0),
			   
	svector<2>(-1.0/sqrt(5.0),-1.0), svector<2>( 1.0/sqrt(5.0), -1.0),
	svector<2>( 1.0,-1.0/sqrt(5.0)), svector<2>( 1.0,  1.0/sqrt(5.0)),
	svector<2>( 1.0/sqrt(5.0), 1.0), svector<2>(-1.0/sqrt(5.0),  1.0),
	svector<2>(-1.0, 1.0/sqrt(5.0)), svector<2>(-1.0, -1.0/sqrt(5.0)),

	svector<2>(-1.0/sqrt(5.0), -1.0/sqrt(5.0)), 
        svector<2>( 1.0/sqrt(5.0), -1.0/sqrt(5.0)),
	svector<2>( 1.0/sqrt(5.0),  1.0/sqrt(5.0)), 
        svector<2>(-1.0/sqrt(5.0),  1.0/sqrt(5.0))};

basisFunction<2> Square12BF =    // two-d serendipity cubic fns
  {12,
   {1,2,0},
   Line4,
   &Square,
   12,
   Square16Pts,                  // first 12 points are the same
   quad12Permutation
  };

LagrangeBasisFunction<2> Square12class(Square12BF, Line4, quad12),
 *Square12 = &Square12class;

basisFunction<2> Square16BF =    // two-d bicubic basis fns
  {16,
   {1,2,4},
   Line4,
   &Square,
   16,
   Square16Pts,
   quad16Permutation
  };

LagrangeBasisFunction<2> Square16class(Square16BF, Line4, quad16),
 *Square16 = &Square16class;

svector<2> Square17Pts[] = {svector<2>(-1.0,-1.0), svector<2>( 1.0,-1.0), 
			    svector<2>( 1.0, 1.0), svector<2>(-1.0, 1.0),
			   
  svector<2>(-sqrt(3.0/7.0),-1), svector<2>( 0,-1), svector<2>(sqrt(3.0/7.0),  -1),
  svector<2>( 1,-sqrt(3.0/7.0)), svector<2>( 1, 0), svector<2>( 1,  sqrt(3.0/7.0)),
  svector<2>( sqrt(3.0/7.0), 1), svector<2>( 0, 1), svector<2>(-sqrt(3.0/7.0),  1),
  svector<2>(-1, sqrt(3.0/7.0)), svector<2>(-1, 0), svector<2>(-1, -sqrt(3.0/7.0)),

  svector<2>(0.0, 0.0)};

basisFunction<2> Square17BF =    // two-d serendipity quartic basis fns
  {17,
   {1,3,1},
   Line5,
   &Square,
   17,
   Square17Pts,
   quad17Permutation
  };

LagrangeBasisFunction<2> Square17class(Square17BF, Line5, quad17),
  *Square17 = &Square17class;

svector<2> Square24Pts[] = 
  {svector<2>(-1.0,-1.0), svector<2>( 1.0,-1.0), 
   svector<2>( 1.0, 1.0), svector<2>(-1.0, 1.0),

   svector<2>(-sqrt((3+2*sqrt(6.0/5))/7.0), -1),
   svector<2>(-sqrt((3-2*sqrt(6.0/5))/7.0), -1),
   svector<2>( sqrt((3-2*sqrt(6.0/5))/7.0), -1),
   svector<2>( sqrt((3+2*sqrt(6.0/5))/7.0), -1),

   svector<2>(1, -sqrt((3+2*sqrt(6.0/5))/7.0)),
   svector<2>(1, -sqrt((3-2*sqrt(6.0/5))/7.0)),
   svector<2>(1,  sqrt((3-2*sqrt(6.0/5))/7.0)),
   svector<2>(1,  sqrt((3+2*sqrt(6.0/5))/7.0)),

   svector<2>( sqrt((3+2*sqrt(6.0/5))/7.0), 1),
   svector<2>( sqrt((3-2*sqrt(6.0/5))/7.0), 1),
   svector<2>(-sqrt((3-2*sqrt(6.0/5))/7.0), 1),
   svector<2>(-sqrt((3+2*sqrt(6.0/5))/7.0), 1),

   svector<2>(-1,  sqrt((3+2*sqrt(6.0/5))/7.0)),
   svector<2>(-1,  sqrt((3-2*sqrt(6.0/5))/7.0)),
   svector<2>(-1, -sqrt((3-2*sqrt(6.0/5))/7.0)),
   svector<2>(-1, -sqrt((3+2*sqrt(6.0/5))/7.0)),

   svector<2>(-sqrt(3.0)/4.0, -sqrt(3.0)/4.0),
   svector<2>( sqrt(3.0)/4.0, -sqrt(3.0)/4.0),
   svector<2>( sqrt(3.0)/4.0,  sqrt(3.0)/4.0),
   svector<2>(-sqrt(3.0)/4.0,  sqrt(3.0)/4.0)};

basisFunction<2> Square24BF =    // two-d serendipity quintic basis fns
  {24,
   {1,4,4},
   Line6,
   &Square,
   24,
   Square24Pts,
   quad24Permutation
  };

LagrangeBasisFunction<2> Square24class(Square24BF, Line6, quad24),
 *Square24 = &Square24class;

// ***************** 3d Functions on Cube **********************
// ***************** 3d Functions on Cube **********************
// ***************** 3d Functions on Cube **********************

svector<3> Cube27Pts[] = 
  {svector<3>(-1,-1,-1), svector<3>(1,-1,-1),  // vertices
   svector<3>(-1, 1,-1), svector<3>(1, 1,-1), 
   svector<3>(-1,-1, 1), svector<3>(1,-1, 1), 
   svector<3>(-1, 1, 1), svector<3>(1, 1, 1),
   svector<3>(0, -1,-1), svector<3>(0, 1, -1),  // mid edges
   svector<3>(0, -1, 1), svector<3>(0, 1, 1), 
   svector<3>(-1, 0,-1), svector<3>(1, 0,-1), 
   svector<3>(-1, 0, 1), svector<3>(1, 0, 1), 
   svector<3>(-1,-1, 0), svector<3>(1,-1, 0), 
   svector<3>(-1, 1, 0), svector<3>(1, 1, 0),
   svector<3>(0, 0, -1), svector<3>(0, 0, 1),  // mid faces 
   svector<3>(0, -1, 0), svector<3>(0, 1, 0), 
   svector<3>(-1, 0, 0), svector<3>(1, 0, 0),
   svector<3>( 0, 0, 0)};                      // center

basisFunction<3> Cube8BF =    // two-d biquadratic basis fns
  {8,
   {1,0,0,0},
   Square4,
   &Cube,
   8,
   Cube27Pts,           // Cube8, Cube20, Cube 27 interpolation points nest
   identityPermutation
  };

LagrangeBasisFunction<3> Cube8class(Cube8BF, Square4, cube8),
  *Cube8 = &Cube8class;

svector<3> Cube8BubblePts[] = 
  {svector<3>(-1,-1,-1), svector<3>(1,-1,-1),  // vertices
   svector<3>(-1, 1,-1), svector<3>(1, 1,-1), 
   svector<3>(-1,-1, 1), svector<3>(1,-1, 1), 
   svector<3>(-1, 1, 1), svector<3>(1, 1, 1),
   svector<3>( 0, 0, 0)};                      // center

basisFunction<3> Cube8BubbleBF =    // two-d biquadratic basis fns
  {9,
   {1,0,0,1},
   Square4,
   &Cube,
   8,                         // not clear how to initialize bubble
   Cube8BubblePts,
   identityPermutation,
  };

LagrangeBasisFunction<3> 
  Cube8BubbleClass(Cube8BubbleBF, Square4, cube8Bubble),
 *Cube8Bubble = &Cube8BubbleClass;

basisFunction<3> Cube20BF =    // two-d biquadratic basis fns
  {20,
   {1,1,0,0},
   Square8,
   &Cube,
   20,
   Cube27Pts,             // points same as first 20 tensor product
   identityPermutation
  };

LagrangeBasisFunction<3> Cube20class(Cube20BF, Square8, cube20),
  *Cube20 = &Cube20class;

basisFunction<3> Cube27BF =    // two-d biquadratic basis fns
  {27,
   {1,1,1,1},
   Square9,
   &Cube,
   27,
   Cube27Pts,
   identityPermutation
  };

LagrangeBasisFunction<3> Cube27class(Cube27BF, Square9, cube27),
  *Cube27 = &Cube27class;

svector<3> Cube32Pts[] = 
  {svector<3>(-1,-1,-1), svector<3>(1,-1,-1),  // vertices
   svector<3>(-1, 1,-1), svector<3>(1, 1,-1), 
   svector<3>(-1,-1, 1), svector<3>(1,-1, 1), 
   svector<3>(-1, 1, 1), svector<3>(1, 1, 1),

   svector<3>(-1.0/sqrt(5.0), -1, -1), svector<3>( 1.0/sqrt(5.0), -1, -1), 
   svector<3>(-1.0/sqrt(5.0),  1, -1), svector<3>( 1.0/sqrt(5.0),  1, -1), 
   svector<3>(-1.0/sqrt(5.0), -1,  1), svector<3>( 1.0/sqrt(5.0), -1,  1), 
   svector<3>(-1.0/sqrt(5.0),  1,  1), svector<3>( 1.0/sqrt(5.0),  1,  1), 

   svector<3>(-1, -1.0/sqrt(5.0), -1), svector<3>(-1,  1.0/sqrt(5.0), -1), 
   svector<3>( 1, -1.0/sqrt(5.0), -1), svector<3>( 1,  1.0/sqrt(5.0), -1), 
   svector<3>(-1, -1.0/sqrt(5.0),  1), svector<3>(-1,  1.0/sqrt(5.0),  1), 
   svector<3>( 1, -1.0/sqrt(5.0),  1), svector<3>( 1,  1.0/sqrt(5.0),  1), 

   svector<3>(-1, -1, -1.0/sqrt(5.0)), svector<3>(-1, -1,  1.0/sqrt(5.0)), 
   svector<3>( 1, -1, -1.0/sqrt(5.0)), svector<3>( 1, -1,  1.0/sqrt(5.0)),
   svector<3>(-1,  1, -1.0/sqrt(5.0)), svector<3>(-1,  1,  1.0/sqrt(5.0)), 
   svector<3>( 1,  1, -1.0/sqrt(5.0)), svector<3>( 1,  1,  1.0/sqrt(5.0))};

basisFunction<3> Cube32BF =    // serendipity cubic cube
  {32,
   {1,2,0,0},
   Square12,
   &Cube,
   32,
   Cube32Pts,
   cube32Permutation
  };

LagrangeBasisFunction<3> Cube32class(Cube32BF, Square12, cube32),
  *Cube32 = &Cube32class;

svector<3> Cube50Pts[] = 
  {svector<3>(-1,-1,-1), svector<3>(1,-1,-1),  // vertices
   svector<3>(-1, 1,-1), svector<3>(1, 1,-1), 
   svector<3>(-1,-1, 1), svector<3>(1,-1, 1), 
   svector<3>(-1, 1, 1), svector<3>(1, 1, 1),

   svector<3>(-sqrt(3.0/7), -1, -1),    // edges
   svector<3>(0, -1, -1),
   svector<3>(sqrt(3.0/7), -1, -1),
   svector<3>(sqrt(3.0/7), 1, -1),
   svector<3>(0, 1, -1),
   svector<3>(-sqrt(3.0/7), 1, -1),
   svector<3>(sqrt(3.0/7), -1, 1),
   svector<3>(0, -1, 1),
   svector<3>(-sqrt(3.0/7), -1, 1),
   svector<3>(sqrt(3.0/7), 1, 1),
   svector<3>(0, 1, 1),
   svector<3>(-sqrt(3.0/7), 1, 1),
   svector<3>(-1, sqrt(3.0/7), -1),
   svector<3>(-1, 0, -1),
   svector<3>(-1, -sqrt(3.0/7), -1),
   svector<3>(1, -sqrt(3.0/7), -1),
   svector<3>(1, 0, -1),
   svector<3>(1, sqrt(3.0/7), -1),
   svector<3>(-1, sqrt(3.0/7), 1),
   svector<3>(-1, 0, 1),
   svector<3>(-1, -sqrt(3.0/7), 1),
   svector<3>(1, -sqrt(3.0/7), 1),
   svector<3>(1, 0, 1),
   svector<3>(1, sqrt(3.0/7), 1),
   svector<3>(-1, -1, -sqrt(3.0/7)),
   svector<3>(-1, -1, 0),
   svector<3>(-1, -1, sqrt(3.0/7)),
   svector<3>(1, -1, -sqrt(3.0/7)),
   svector<3>(1, -1, 0),
   svector<3>(1, -1, sqrt(3.0/7)),
   svector<3>(-1, 1, -sqrt(3.0/7)),
   svector<3>(-1, 1, 0),
   svector<3>(-1, 1, sqrt(3.0/7)),
   svector<3>(1, 1, -sqrt(3.0/7)),
   svector<3>(1, 1, 0),
   svector<3>(1, 1, sqrt(3.0/7)),

   svector<3>(0, 0, -1), svector<3>(0, 0, 1),  // mid faces 
   svector<3>(0, -1, 0), svector<3>(0, 1, 0), 
   svector<3>(-1, 0, 0), svector<3>(1, 0, 0)};

basisFunction<3> Cube50BF =    // serendipity cubic cube
  {50,
   {1,3,1,0},
   Square17,
   &Cube,
   50,
   Cube50Pts,
   cube50Permutation
  };

LagrangeBasisFunction<3> Cube50class(Cube50BF, Square17, cube50),
  *Cube50 = &Cube50class;

// ***************** 3d Functions on Tetrahedra **********************
// ***************** 3d Functions on Tetrahedra **********************
// ***************** 3d Functions on Tetrahedra **********************

svector<3> Tet10Pts[] = 
  {svector<3>(0, 0, 0), svector<3>(1, 0, 0),  // vertices
   svector<3>(0, 1, 0), svector<3>(0, 0, 1), 
   svector<3>(0.5, 0.0, 0.0), svector<3>(0.5, 0.5, 0.0),  // mid edges
   svector<3>(0.0, 0.5, 0.0), svector<3>(0.0, 0.0, 0.5), 
   svector<3>(0.5, 0.0, 0.5), svector<3>(0.0, 0.5, 0.5)};

basisFunction<3> Tet4BF =   
  {4,
   {1,0,0,0},
   Triangle3,
   &Tetrahedron,
   4,
   Tet10Pts,
   identityPermutation
  };

LagrangeBasisFunction<3> Tet4class(Tet4BF, Triangle3, tet4),
 *Tet4 = &Tet4class;

svector<3> Tet4BubblePts[] = 
  {svector<3>(0, 0, 0), svector<3>(1, 0, 0),  // vertices
   svector<3>(0, 1, 0), svector<3>(0, 0, 1), 
   svector<3>(1.0/4.0, 1.0/4.0, 1.0/4.0)};

basisFunction<3> Tet4BubbleBF =   
  {5,
   {1,0,0,1},
   Triangle3,
   &Tetrahedron,
   4,
   Tet4BubblePts,
   identityPermutation
  };

LagrangeBasisFunction<3> 
  Tet4BubbleClass(Tet4BubbleBF, Triangle3, tet4Bubble),
  *Tet4Bubble = &Tet4BubbleClass;

// Lagrange type bubble

basisFunction<3> Tet4BubbleInterpBF =   
  {5,
   {1,0,0,1},
   Triangle3,
   &Tetrahedron,
   5,
   Tet4BubblePts,
   identityPermutation
  };

LagrangeBasisFunction<3> 
  Tet4BubbleInterpClass(Tet4BubbleInterpBF, Triangle3, tet4BubbleInterp),
  *Tet4BubbleInterp = &Tet4BubbleInterpClass;

basisFunction<3> Tet10BF =  
  {10,
   {1,1,0,0},
   Triangle6,
   &Tetrahedron,
   10,
   Tet10Pts, 
   identityPermutation
  };

LagrangeBasisFunction<3> Tet10class(Tet10BF, Triangle6, tet10),
 *Tet10 = &Tet10class;

svector<3> Tet20Pts[] = 
  {svector<3>(0, 0, 0), svector<3>(1, 0, 0),  // vertices
   svector<3>(0, 1, 0), svector<3>(0, 0, 1), 

   svector<3>(0.5*(1-1.0/sqrt(5.0)), 0.0, 0.0),   // bottom edges
   svector<3>(0.5*(1+1.0/sqrt(5.0)), 0.0, 0.0),
   svector<3>(0.5*(1+1.0/sqrt(5.0)), 0.5*(1-1.0/sqrt(5.0)), 0.0), 
   svector<3>(0.5*(1-1.0/sqrt(5.0)), 0.5*(1+1.0/sqrt(5.0)), 0.0),
   svector<3>(0.0, 0.5*(1+1.0/sqrt(5.0)), 0.0), 
   svector<3>(0.0, 0.5*(1-1.0/sqrt(5.0)), 0.0),

   svector<3>(0.0, 0.0, 0.5*(1+1.0/sqrt(5.0))),   // edges from apex 
   svector<3>(0.0, 0.0, 0.5*(1-1.0/sqrt(5.0))),
   svector<3>(0.5*(1-1.0/sqrt(5.0)), 0.0, 0.5*(1+1.0/sqrt(5.0))),
   svector<3>(0.5*(1+1.0/sqrt(5.0)), 0.0, 0.5*(1-1.0/sqrt(5.0))), 
   svector<3>(0.0, 0.5*(1-1.0/sqrt(5.0)), 0.5*(1+1.0/sqrt(5.0))),
   svector<3>(0.0, 0.5*(1+1.0/sqrt(5.0)), 0.5*(1-1.0/sqrt(5.0))), 

   svector<3>(1.0/3.0, 1.0/3.0, 1.0/3.0),     // mid-faces
   svector<3>(0.0, 1.0/3.0, 1.0/3.0),
   svector<3>(1.0/3.0, 0.0, 1.0/3.0),
   svector<3>(1.0/3.0, 1.0/3.0, 0.0)};

basisFunction<3> Tet20BF =  
  {20,
   {1,2,1,0},
   Triangle10,
   &Tetrahedron,
   20,
   Tet20Pts, 
   tet20Permutation
  };

LagrangeBasisFunction<3> Tet20class(Tet20BF, Triangle10, tet20),
 *Tet20 = &Tet20class;

#endif                          // !defined(__LAGRANGEBASISFNS_H__)
