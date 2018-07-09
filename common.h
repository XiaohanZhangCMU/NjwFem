#include"LinearAlgebra/svector.C"
#include"LinearAlgebra/smatrix.C"
#include"LinearAlgebra/sbmatrix.C"
#include"LinearAlgebra/fourTensor.C"
// Define dvector/dmatrix to be the small (Ndim) vector/matrix class

typedef fourTensor<Ndim> dfourTensor;
typedef svector<Ndim> dvector;
typedef smatrix<Ndim> dmatrix;
typedef sbmatrix<Ndim> bmatrix;
typedef svector<Ndim-1> bvector;    // boundary Gauss points

// Get the big matrix/vector class

class FeMdata
{
 public:
  int nvar;      // size of the matrix/vectors
  int nzero;     // estimated number of zeros per row
  int band;      // bandwidth (for banded solvers)
  bool pivot;    // pivoting required for direct solvers
  bool pcflag;   // seperate preconditioner for itterative solvers
  double tol;    // tolerence for itterative solvers
  feMesh *mm;
  FeMdata() 
    : nvar(0), nzero(0), band(0), pivot(false), 
    pcflag(false), tol(1.0e-10) , mm(NULL)
    {return;}

  FeMdata(int dof, feMesh *fem)
    : nvar(dof),mm(fem)
    {
      nzRowBand(*mm,nzero,band);
    }
};

#ifdef BandMatrix
  #include"LinearAlgebra/bandMatrix.C"
  typedef bandMV feMV;
#endif

#ifdef SpoolesMatrix
  #include"LinearAlgebra/spoolesMatrix.C"
  typedef spoolesMV feMV;
#endif

#ifdef PetscMatrix
// petsc 2.3
// #define PETSC23
// petsc 3.0

#define PETSC40
  #include"LinearAlgebra/petscMatrix.C"
  typedef petscMV feMV;
#endif

#ifdef MPIPetscMatrix
  #define PETSC40
  #include "LinearAlgebra/petscMatrixMPI.C"
  typedef petscMV feMV;
#endif

#ifdef LapackMatrix
  #include"LinearAlgebra/lapackMatrix.C"
  typedef lapackMV feMV;
#endif

#ifdef SuperluMatrix
  #include"LinearAlgebra/superluMatrix.C"
  typedef superluMV feMV;
#endif

#ifdef PardisoMatrix
  #include"LinearAlgebra/pardisoMatrix.C"
  typedef pardisoMV feMV;
#endif


#ifdef MumpsMatrix
  #include"LinearAlgebra/mumpsMatrix.C"
  typedef mumpsMV feMV;
#endif

#ifdef PastixMatrix
  #include"LinearAlgebra/pastixMatrix.C"
  typedef pastixMV feMV;
#endif

// Set up Gauss quadrature points

namespace QuadPtWt
{
  // Guass points on [-1,1]

  double GP3[] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
  double GW3[] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

  //  double GP3[] = {-0.774596669241483,0.0,0.774596669241483};
  //  double GW3[] = {0.555555555555556,0.88888888888889,0.555555555555556};

  double 
    gp4pt1 = sqrt(3.0 - 2.0*sqrt(6.0/5.0))/7.0, 
    gp4pt2 = sqrt(3.0 + 2.0*sqrt(6.0/5.0))/7.0,

    gp4wt1 = 0.5 + sqrt(30.0)/36.0,
    gp4wt2 = 0.5 - sqrt(30.0)/36.0;

  double GP4[] = {-gp4pt2, -gp4pt1, gp4pt1, gp4pt2};
  double GW4[] = { gp4wt2,  gp4wt1, gp4wt1, gp4wt2};

  double 
    gp5pt1 = sqrt(5.0 - 2.0*sqrt(10.0/7.0))/3.0, 
    gp5pt2 = sqrt(5.0 + 2.0*sqrt(10.0/7.0))/3.0,

    gp5wt1 = 161.0/450.0 + 13.0*sqrt(70.0)/900.0,
    gp5wt2 = 161.0/450.0 - 13.0*sqrt(70.0)/900.0;

  double GP5[] = {-gp5pt2, -gp5pt1,     0.0,     gp5pt1, gp5pt2};
  double GW5[] = { gp5wt2,  gp5wt1, 128.0/225.0, gp5wt1, gp5wt2};

  //  double GW5[] = {0.236926885056189, 0.478628670499366,0.568888888888889,
  //		      0.236926885056189, 0.478628670499366};

  //  double GP5[] = {-0.906179845938664,-0.538469310105683,0.0,
  //	  	       0.906179845938664, 0.538469310105683};
  //  double GW5[] = {0.236926885056189, 0.478628670499366,0.568888888888889,
  //		      0.236926885056189, 0.478628670499366};
  
  double GP7[] = {-0.949107912342759, -0.741531185599394, -0.405845151377397,
		  0.000000000000000,
		  0.405845151377397,  0.741531185599394, 0.949107912342759};
  
  double GW7[] = {0.129484966168870, 0.279705391489277, 0.381830050505119,
		  0.417959183673469,
		  0.381830050505119, 0.279705391489277, 0.129484966168870};
  
  // Triangle: see http://arxiv.org/abs/math/0501496v2 for up to degree 25

  // Seven point quadrature rule for the trinagle (exact on P5)
  // Note the Gauss points are  B =  2/7 (+ or -) sqrt(15)/21
  // and  the Gauss weights are W = 31/480 (+ -)  sqrt(15)/2400
  // The center point and weight below are exact
  
  double 
    tri7pt1 = 2.0/7.0 + sqrt(15.0)/21.0, tri7pt1a = 1.0-2.0*tri7pt1, 
    tri7pt2 = 2.0/7.0 - sqrt(15.0)/21.0, tri7pt2a = 1.0-2.0*tri7pt2, 
    tri7wt1 = 31.0/480.0 + sqrt(15.0)/2400.0,
    tri7wt2 = 31.0/480.0 - sqrt(15.0)/2400.0;

  //  double BB1=0.4701420641051151,AA1=1.0-2.0*BB1,W1=0.1323941527885062/2.0;
  //  double BB2=0.1012865073234563,AA2=1.0-2.0*BB2,W2=0.1259391805448272/2.0;
  
  svector<2> GpTri7[] = {svector<2>(1.0/3.0, 1.0/3.0),
			 svector<2>(tri7pt1a, tri7pt1), 
			 svector<2>(tri7pt1,  tri7pt1a), 
			 svector<2>(tri7pt1,  tri7pt1),
			 svector<2>(tri7pt2a, tri7pt2), 
			 svector<2>(tri7pt2,  tri7pt2a), 
			 svector<2>(tri7pt2,  tri7pt2)};

  double GwTri7[] = {0.1125, tri7wt1,tri7wt1,tri7wt1, tri7wt2,tri7wt2,tri7wt2};

  // Twelve point quadrature rule for the triangle (exact on P6)

  double GwTri12[] = {
	0.05839314, 0.05839314, 0.05839314,
	0.02542245, 0.02542245, 0.02542245,
	0.04142554, 0.04142554, 0.04142554,
	0.04142554, 0.04142554, 0.04142554};
      
  svector<2> GpTri12[] = {
	svector<2>(0.5014265, 0.2492867),
	svector<2>(0.2492867, 0.5014265),
	svector<2>(0.2492867, 0.2492867),
	svector<2>(0.8738220, 0.06308901),
	svector<2>(0.06308901,0.8738220),
	svector<2>(0.06308901,0.06308901),
	svector<2>(0.6365025, 0.05314505),
	svector<2>(0.6365025, 0.3103525),
	svector<2>(0.05314505,0.6365025),
	svector<2>(0.05314505,0.3103525),
	svector<2>(0.3103525, 0.6365025),
	svector<2>(0.3103525, 0.05314505)};
      
  // 15 point quadrature rule for the triangle (exact on P7)
      
  double GwTri15[] = { 
	0.25*0.0102558174092, 0.25*0.0102558174092, 0.25*0.0102558174092,
	0.25*0.1116047046647, 0.25*0.1116047046647, 0.25*0.1116047046647,
	0.25*0.1116047046647, 0.25*0.1116047046647, 0.25*0.1116047046647,
	0.25*0.1679775595335, 0.25*0.1679775595335, 0.25*0.1679775595335,
	0.25*0.2652238803946, 0.25*0.2652238803946, 0.25*0.2652238803946};
      
  svector<2> GpTri15[] = {
	svector<2>(1.0000000000000, 0.0000000000000),
	svector<2>(0.0000000000000, 0.0000000000000),
	svector<2>(0.0000000000000, 1.0000000000000),
	svector<2>(0.7839656651012, 0.0421382841642),
	svector<2>(0.1738960507345, 0.7839656651012),
	svector<2>(0.1738960507345, 0.0421382841642),
	svector<2>(0.0421382841642, 0.1738960507345),
	svector<2>(0.7839656651012, 0.1738960507345),
	svector<2>(0.0421382841642, 0.7839656651012),
	svector<2>(0.4743880861752, 0.4743880861752),
	svector<2>(0.4743880861752, 0.0512238276497),
	svector<2>(0.0512238276497, 0.4743880861752),
	svector<2>(0.2385615300181, 0.5228769399639),
	svector<2>(0.5228769399639, 0.2385615300181),
	svector<2>(0.2385615300181, 0.2385615300181)};

  double     GwCompositeTri21[21], GwCompositeTri36[36], GwCompositeTri45[45];
  svector<2> GpCompositeTri21[21], GpCompositeTri36[36], GpCompositeTri45[45];


  // Tetrahedron: 14 point quadrature rule for the tetrahedra (exact on P5)

  static double a1 = 0.31088591926330060980,  b1 = 1.0 - 3.0*a1; 
  static double a2 = 0.092735250310891226402, b2 = 1.0 - 3.0*a2;
  static double a3 = 0.045503704125649649492, b3 = 0.5 - a3;
  static double w1 = 0.018781320953002641800;
  static double w2 = 0.012248840519393658257;
  static double w3 = 0.0070910034628469110730;
  
  svector<3> GpTet14[] = 
	{svector<3>(a1,a1,a1), svector<3>(a1,a1,b1), 
	 svector<3>(a1,b1,a1), svector<3>(b1,a1,a1),
	 svector<3>(a2,a2,a2), svector<3>(a2,a2,b2), 
	 svector<3>(a2,b2,a2), svector<3>(b2,a2,a2),
	 svector<3>(a3,a3,b3), svector<3>(a3,b3,a3), svector<3>(b3,a3,a3),
	 svector<3>(b3,b3,a3), svector<3>(b3,a3,b3), svector<3>(a3,b3,b3)};
      
  double GwTet14[] = {w1,w1,w1,w1, w2,w2,w2,w2, w3,w3,w3, w3,w3,w3};

  double     GwCompositeTet56[56];
  svector<3> GpCompositeTet56[56];
}

class QuadratureRuleClass
{
 public:

  void interval(const int &n, double **gp, double **gw)
  {
    if(n == 3)      {*gp = QuadPtWt::GP3; *gw = QuadPtWt::GW3;}
    else if(n == 4) {*gp = QuadPtWt::GP4; *gw = QuadPtWt::GW4;}
    else if(n == 5) {*gp = QuadPtWt::GP5; *gw = QuadPtWt::GW5;}
    else if(n == 7) {*gp = QuadPtWt::GP7; *gw = QuadPtWt::GW7;}
    else
      assert(false);

    return;
  }

  void triangle(const int &n, svector<2> **gp, double **gw)
    {
      if(n == 7)       {*gp = QuadPtWt::GpTri7;  *gw = QuadPtWt::GwTri7;}
      else if(n == 12) {*gp = QuadPtWt::GpTri12; *gw = QuadPtWt::GwTri12;}
      else if(n == 15) {*gp = QuadPtWt::GpTri15; *gw = QuadPtWt::GwTri15;}
      else
	assert(false);
      
      return;
    }

  void compositeTriangle(const int &n, svector<2> **gp, double **gw)
    {
      // Quadrature rule piecewise smooth functions on
      // the baracentric subdivision of a trinagle

      int nt;
      svector<2> *gpt; double *gwt;

      if(n == 21)   
	{
	  nt = 7;
	  *gp = QuadPtWt::GpCompositeTri21; gpt = QuadPtWt::GpTri7;
	  *gw = QuadPtWt::GwCompositeTri21; gwt = QuadPtWt::GwTri7;
	}
      else if(n == 36) 
	{
	  nt = 12;
	  *gp = QuadPtWt::GpCompositeTri36; gpt = QuadPtWt::GpTri12;
	  *gw = QuadPtWt::GwCompositeTri36; gwt = QuadPtWt::GwTri12;
	}
      else if(n == 45) 
	{
	  nt = 15;
	  *gp = QuadPtWt::GpCompositeTri45; gpt = QuadPtWt::GpTri15;
	  *gw = QuadPtWt::GwCompositeTri45; gwt = QuadPtWt::GwTri15;
	}
      else
	assert(false);

      svector<2> cc(1.0/3.0,1.0/3.0), e0(0,0), e1(1,0), e2(0,1);
      svector<2> kk[3][3] = {{cc,e0,e1}, {cc,e1,e2}, {cc,e2,e0}};

      for(int k = 0; k < 3; k++)
	{
	  // ff = [f0,f1] 

          svector<2> f0 = kk[k][1]-kk[k][0], f1 = kk[k][2]-kk[k][0];
	  smatrix<2> ff(f0[0],f1[0], f0[1], f1[1]);

	  for(int i = 0; i < nt; i++)
	    {
	      (*gp)[k*nt+i] = cc + ff*gpt[i];  (*gw)[k*nt+i] = gwt[i]/3.0;
	    }
	}

      return;
    }

  void tetrahedra(const int &n, svector<3> **gp, double **gw)
    {
      if(n == 14) {*gp = QuadPtWt::GpTet14;  *gw = QuadPtWt::GwTet14;}
      else
	assert(false);

      return;
    }

  void compositeTetrahedra(const int &n, svector<3> **gp, double **gw)
    {
      // Quadrature rule piecewise smooth functions on
      // the baracentric subdivision of a trinagle

      int nt;
      svector<3> *gpt; double *gwt;

      if(n == 56)   
	{
	  nt = 14;
	  *gp = QuadPtWt::GpCompositeTet56; gpt = QuadPtWt::GpTet14;
	  *gw = QuadPtWt::GwCompositeTet56; gwt = QuadPtWt::GwTet14;
	}
      else
	assert(false);

      svector<3> cc(0.25, 0.25, 0.25), 
	e0(0,0,0), e1(1,0,0), 
	e2(0,1,0), e3(0,0,1);

      svector<3> kk[4][4] =
	{{cc,e1,e2,e3}, {cc,e0,e3,e2}, {cc,e0,e1,e3}, {cc,e0,e2,e1}};

      for(int k = 0; k < 4; k++)
	{
	  // ff = [f0,f1,f2] 

          svector<3> 
	    f0 = kk[k][1]-kk[k][0], 
	    f1 = kk[k][2]-kk[k][0],
	    f2 = kk[k][3]-kk[k][0];

	  smatrix<3> 
	    ff(f0[0],f1[0],f2[0], f0[1],f1[1],f2[1], f0[2],f1[2],f2[2]);

	  for(int i = 0; i < nt; i++)
	    {
	      (*gp)[k*nt+i] = cc + ff*gpt[i];  (*gw)[k*nt+i] = gwt[i]/4.0;
	    }
	}

      return;
    }

}  QuadratureRules;


