// Hermite Basis Functions
// The boundary function must return the functions for both u and du/dn
// If normal derivatives are used as dof's signs need to be determined

class felement;
class bfelement;

template<const int Ndim>
struct HermiteBasisFunction : public basisFunction<Ndim>
{
  HermiteBasisFunction(basisFunction<Ndim> &bf, 
    void (*bdy)(svector<Ndim-1>, double *phi, double *phin),
    void (*bb) (svector<Ndim>, 
		double *phi, svector<Ndim> *dphi, smatrix<Ndim> *d2phi),
    void (*bdytt)(int *nodes, sbmatrix<Ndim> &dxds,
		double *Phi, double *Phin, double *phi, double *phin),
    void (*tt)(int *nodes, smatrix<Ndim> &dxds,
	       double *Phi, svector<Ndim> *Dphi, smatrix<Ndim> *D2phi,
	       double *phi, svector<Ndim> *dphi, smatrix<Ndim> *d2phi))
    : basisFunction<Ndim>(bf), 
      bdyBasis(bdy), basis(bb), bdyTransform(bdytt), transform(tt)
  {return;}

  void (*bdyBasis)(svector<Ndim-1>, double *phi, double *phin);
  void (*basis)(svector<Ndim>, 
		double *phi, svector<Ndim> *dphi, smatrix<Ndim> *d2phi);
  void (*bdyTransform)(int *nodes, sbmatrix<Ndim> &dxds,
		       double *Phi, double *Phin, double *phi, double *phin);
  void (*transform)(int *nodes, smatrix<Ndim> &dxds,
   		    double *Phi, svector<Ndim> *Dphi, smatrix<Ndim> *D2phi,
		    double *phi, svector<Ndim> *dphi, smatrix<Ndim> *d2phi);
};

// CloughTocher Triangle
// CloughTocher Triangle
// CloughTocher Triangle

void cloughTocherTriangleBdy(svector<1> xi, double psi[4], double psin[3]);

void cloughTocherTriangle(svector<2> xi, double psi[12], 
	       svector<2> *dpsi=NULL, smatrix<2> *d2psi=NULL);

void cloughTocherBdyTransform(int *nodes, sbmatrix<2> &dxds,
		      double Phi[4], double Phin[3], 
		      double phi[7], double phin[7]);

void cloughTocherTransform(int *nodes, smatrix<2> &dxds,
		   double Phi[12], svector<2> Dphi[12], smatrix<2> D2phi[12],
		   double phi[12], svector<2> dphi[12], smatrix<2> d2phi[12]);

basisFunction<1> CloughTocherBdyBF =
  {7,
   {3,1},                     // vertex gets all 3 dof's from the triangle
   Vertex1,                   // boundary not used
   &Interval,
   3,
   Line3Pts,                  
   identityPermutation
  };

basisFunction<2> CloughTocherTriangleBF =  
  {12,
   {3,1,0},                 // u, grad(u)
   &CloughTocherBdyBF,
   &Triangle,
   6,
   Triangle6Pts,            
   identityPermutation
  };

HermiteBasisFunction<2>  CloughTocherTriangleClass(CloughTocherTriangleBF,
      &cloughTocherTriangleBdy, &cloughTocherTriangle, 
      &cloughTocherBdyTransform, &cloughTocherTransform),
  *CloughTocherTriangle = &CloughTocherTriangleClass;	    

// Bell Triangle
// Bell Triangle
// Bell Triangle

void bellTriangleBdy(svector<1> xi, double psi[6], double psin[4]);

void bellTriangle(svector<2> xi, double psi[18], 
	       svector<2> *dpsi=NULL, smatrix<2> *d2psi=NULL);

void bellBdyTransform(int *nodes, sbmatrix<2> &dxds,
		      double Phi[6],  double Phin[4], 
		      double phi[12], double phin[12]);

void bellTransform(int *nodes, smatrix<2> &dxds,
		   double Phi[18], svector<2> Dphi[18], smatrix<2> D2phi[18],
		   double phi[18], svector<2> dphi[18], smatrix<2> d2phi[18]);

basisFunction<1> BellBdyBF =
  {12,
   {6,0},                     // vertex gets all 6 dof's from the triangle
   Vertex1,                   // boundary not used
   &Interval,
   2,
   Line3Pts,                  // first two points are -1 and +1
   identityPermutation
  };

basisFunction<2> BellTriangleBF =  
  {18,
   {6,0,0},                 // u, grad(u) and D2(u) (symmetric)
   &BellBdyBF,
   &Triangle,
   3,
   Triangle6Pts,            // first 3 points are the vertices
   identityPermutation
  };

HermiteBasisFunction<2>  BellTriangleClass(BellTriangleBF,
      &bellTriangleBdy, &bellTriangle, &bellBdyTransform, &bellTransform),
  *BellTriangle = &BellTriangleClass;	    

// Argyris Triangle
// Argyris Triangle
// Argyris Triangle

void argyrisTriangleBdy(svector<1> xi, double psi[6], double psin[5]);

void argyrisTriangle(svector<2> xi, double psi[21], 
	       svector<2> *dpsi=NULL, smatrix<2> *d2psi=NULL);

void argyrisBdyTransform(int *nodes, sbmatrix<2> &dxds,
		      double Phi[6],  double Phin[5], 
		      double phi[13], double phin[13]);

void argyrisTransform(int *nodes, smatrix<2> &dxds,
		   double Phi[21], svector<2> Dphi[21], smatrix<2> D2phi[21],
		   double phi[21], svector<2> dphi[21], smatrix<2> d2phi[21]);

basisFunction<1> ArgyrisBdyBF =
  {13,
   {6,1},                     // vertex gets all 6 dof's from the triangle
   Vertex1,                   // boundary not used
   &Interval,
   3,
   Line3Pts,                  
   identityPermutation
  };

basisFunction<2> ArgyrisTriangleBF =  
  {21,
   {6,1,0},                 // u, grad(u) and D2(u) (symmetric)
   &ArgyrisBdyBF,
   &Triangle,
   6,
   Triangle6Pts,            
   identityPermutation
  };

HermiteBasisFunction<2>  ArgyrisTriangleClass(ArgyrisTriangleBF,
      &argyrisTriangleBdy, &argyrisTriangle, 
      &argyrisBdyTransform, &argyrisTransform),
  *ArgyrisTriangle = &ArgyrisTriangleClass;	    

// BFS Square
// BFS Square
// BFS Square

void bfsSquareBdy(svector<1> xi, double psi[4], double psin[4]) ;

void bfsSquare(svector<2> xi, double phi[16], 
	       svector<2> *dphi=NULL, smatrix<2> *d2phi=NULL);

void bfsSquareBdyTransform(int *nodes, sbmatrix<2> &dxds,
		     double Phi[4],  double Phin[4], 
		     double phi[8], double phin[8]);

void bfsSquareTransform(int *nodes, smatrix<2> &dxds,
		  double Phi[16], svector<2> Dphi[16], smatrix<2> D2phi[16],
		  double phi[16], svector<2> dphi[16], smatrix<2> d2phi[16]);

basisFunction<1> BFSsquareBdyBF =
  {8,
   {4,0},                     // vertex gets all 4 dof's from the square
   Vertex1,                   // boundary not used
   &Interval,
   2,
   Line3Pts,                  // first two points are -1 and +1
   identityPermutation
  };

basisFunction<2> BFSsquareBF =  
  {16,
   {4,0,0},                 // u, grad(u) and u_xy
   &BFSsquareBdyBF,
   &Square,
   4,
   Square9Pts,            // first 4 points are the vertices
   identityPermutation
  };

HermiteBasisFunction<2>  BFSsquareClass(BFSsquareBF,
      &bfsSquareBdy, &bfsSquare, &bfsSquareBdyTransform, &bfsSquareTransform),
  *BFSsquare = &BFSsquareClass;	    

// BFS Cube
// BFS Cube
// BFS Cube

void bfsCubeBdy(svector<2> xi, double psi[16], double psin[16]);

void bfsCube(svector<3> xi, double psi[64], 
	     svector<3> dpsi[64], smatrix<3> d2psi[64]);

void bfsCubeBdyTransform(int *nodes, sbmatrix<3> &dxds,
			 double Phi[16], double Phin[16], 
			 double phi[32], double phin[32]);

void bfsCubeTransform(int *nodes, smatrix<3> &dxds,
	   double Phi[64], svector<3> Dphi[64], smatrix<3> D2phi[64],
	   double phi[64], svector<3> dphi[64], smatrix<3> d2phi[64]);

basisFunction<2> BFScubeBdyBF =
  {32,
   {8,0},                     // vertex gets all 8 dof's from the cube
   &BFSsquareBdyBF,           // boundary not used
   &Square,
   4,
   Square9Pts,                // first 4 points are the vertices
   identityPermutation
  };

basisFunction<3> BFScubeBF =  
  {64,
   {8,0,0},                 // u, grad(u) and u_xy, u_xz, u_yz and u_xyz
   &BFScubeBdyBF,
   &Cube,
   8,
   Cube27Pts,               // first 8 points are the vertices
   identityPermutation
  };

HermiteBasisFunction<3>  BFScubeClass(BFScubeBF,
      &bfsCubeBdy, &bfsCube, &bfsCubeBdyTransform, &bfsCubeTransform),
  *BFScube = &BFScubeClass;	    

// Specific Basis and Transformation Functions
// Specific Basis and Transformation Functions
// Specific Basis and Transformation Functions

// CloughTocher Triangle
// CloughTocher Triangle
// CloughTocher Triangle

void cloughTocherTriangle(svector<2> xi, double psi[21], 
			  svector<2> *dpsi, smatrix<2> *d2psi)
{
  double x = xi[0], y = xi[1];

  double omxmy = 1-x-y,  omxmysq=omxmy*omxmy, sqrt2=sqrt(2.0);
  double  x2=x*x, y2=y*y, x3=x2*x, y3=y2*y, eps=1.0e-12;

  // The following basis functions were generated by Maple

  if((y > 0.5*(1-x)-eps) && (y > 2*(0.5-x)-eps)) // triangle opposite vertex 0
    {
      psi[0] = (7*x-1+7*y)*omxmysq;
      psi[1] = (1.0/6)*(13*x-1-2*y)*omxmysq;
      psi[2] = -(1.0/6)*(2*x-13*y+1)*omxmysq;
      psi[3] = -(9.0/2)*x-(21.0/2)*x2*y+9*x2-(9.0/2)*y-(5.0/2)*y3-(21.0/2)*x*y2+6*y2-(9.0/2)*x3+1+15*x*y;
      psi[4] = (7.0/4)*x2*y+(3.0/4)*y+(3.0/4)*x+(5.0/12)*y3-y2+(7.0/4)*x*y2-1.0/6+(17.0/12)*x3-2*x2-(5.0/2)*x*y;
      psi[5] = 2.0/3-(23.0/4)*x2*y-(13.0/4)*x-(25.0/4)*x*y2+(7.0/2)*y2-(11.0/4)*y-(17.0/12)*y3-(23.0/12)*x3+(19.0/2)*x*y+(9.0/2)*x2;
      psi[6] = 9*y2+1-(9.0/2)*y-(9.0/2)*x-(21.0/2)*x2*y-(9.0/2)*y3-(5.0/2)*x3-(21.0/2)*x*y2+6*x2+15*x*y;
      psi[7] = (9.0/2)*y2-(23.0/4)*x*y2+(7.0/2)*x2+(19.0/2)*x*y-(25.0/4)*x2*y-(11.0/4)*x-(13.0/4)*y+2.0/3-(23.0/12)*y3-(17.0/12)*x3;
      psi[8] = (7.0/4)*x*y2-x2-2*y2+(17.0/12)*y3+(3.0/4)*y-1.0/6+(7.0/4)*x2*y-(5.0/2)*x*y+(3.0/4)*x+(5.0/12)*x3;
      psi[9] = -(2.0/3)*(4*x-1+y)*omxmysq;
      psi[10] = -(1.0/3)*sqrt2*(5*x2+16*x*y-7*x+5*y2-7*y+2)*omxmy;
      psi[11] = -(2.0/3)*(x+4*y-1)*omxmysq;

      if(dpsi == NULL) return;

      dpsi[0][0] = -3*(7*x-3+7*y)*omxmy;
      dpsi[0][1] = -3*(7*x-3+7*y)*omxmy;
      dpsi[1][0] = -(1.0/2)*omxmy*(13*x-5+3*y);
      dpsi[1][1] = -omxmy*(4*x-y);
      dpsi[2][0] = omxmy*(x-4*y);
      dpsi[2][1] = -(1.0/2)*(3*x-5+13*y)*omxmy;
      dpsi[3][0] = -9.0/2-21*x*y+18*x-(21.0/2)*y2-(27.0/2)*x2+15*y;
      dpsi[3][1] = -(21.0/2)*x2-9.0/2-(15.0/2)*y2-21*x*y+12*y+15*x;
      dpsi[4][0] = (7.0/2)*x*y+3.0/4+(7.0/4)*y2+(17.0/4)*x2-4*x-(5.0/2)*y;
      dpsi[4][1] = (7.0/4)*x2+3.0/4+(5.0/4)*y2-2*y+(7.0/2)*x*y-(5.0/2)*x;
      dpsi[5][0] = -(23.0/2)*x*y-13.0/4-(25.0/4)*y2-(23.0/4)*x2+(19.0/2)*y+9*x;
      dpsi[5][1] = -(23.0/4)*x2-(25.0/2)*x*y+7*y-11.0/4-(17.0/4)*y2+(19.0/2)*x;
      dpsi[6][0] = -9.0/2-21*x*y-(15.0/2)*x2-(21.0/2)*y2+12*x+15*y;
      dpsi[6][1] = 18*y-9.0/2-(21.0/2)*x2-(27.0/2)*y2-21*x*y+15*x;
      dpsi[7][0] = -(23.0/4)*y2+7*x+(19.0/2)*y-(25.0/2)*x*y-11.0/4-(17.0/4)*x2;
      dpsi[7][1] = 9*y-(23.0/2)*x*y+(19.0/2)*x-(25.0/4)*x2-13.0/4-(23.0/4)*y2;
      dpsi[8][0] = (7.0/4)*y2-2*x+(7.0/2)*x*y-(5.0/2)*y+3.0/4+(5.0/4)*x2;
      dpsi[8][1] = (7.0/2)*x*y-4*y+(17.0/4)*y2+3.0/4+(7.0/4)*x2-(5.0/2)*x;
      dpsi[9][0] = 4*omxmy*(2*x-1+y);
      dpsi[9][1] = 2*omxmy*(3*x+y-1);
      dpsi[10][0] = sqrt2*(7*y2-8*x+14*x*y-10*y+3+5*x2);
      dpsi[10][1] = sqrt2*(7*x2+3+5*y2-8*y+14*x*y-10*x);
      dpsi[11][0] = 2*(x+3*y-1)*omxmy;
      dpsi[11][1] = 4*(x-1+2*y)*omxmy;

      if(d2psi == NULL) return;

      d2psi[0][0][0] = 42*x-30+42*y;
      d2psi[0][0][1] = 42*x-30+42*y;
      d2psi[0][1][0] = d2psi[0][0][1];
      d2psi[0][1][1] = 42*x-30+42*y;
      d2psi[1][0][0] = 13*x-9+8*y;
      d2psi[1][0][1] = 8*x-4+3*y;
      d2psi[1][1][0] = d2psi[1][0][1];
      d2psi[1][1][1] = 3*x+1-2*y;
      d2psi[2][0][0] = -2*x+1+3*y;
      d2psi[2][0][1] = 3*x-4+8*y;
      d2psi[2][1][0] = d2psi[2][0][1];
      d2psi[2][1][1] = 8*x-9+13*y;
      d2psi[3][0][0] = -21*y+18-27*x;
      d2psi[3][0][1] = -21*x-21*y+15;
      d2psi[3][1][0] = d2psi[3][0][1];
      d2psi[3][1][1] = -15*y-21*x+12;
      d2psi[4][0][0] = (7.0/2)*y+(17.0/2)*x-4;
      d2psi[4][0][1] = (7.0/2)*x+(7.0/2)*y-5.0/2;
      d2psi[4][1][0] = d2psi[4][0][1];
      d2psi[4][1][1] = (5.0/2)*y-2+(7.0/2)*x;
      d2psi[5][0][0] = -(23.0/2)*y-(23.0/2)*x+9;
      d2psi[5][0][1] = -(23.0/2)*x-(25.0/2)*y+19.0/2;
      d2psi[5][1][0] = d2psi[5][0][1];
      d2psi[5][1][1] = -(25.0/2)*x+7-(17.0/2)*y;
      d2psi[6][0][0] = -21*y-15*x+12;
      d2psi[6][0][1] = -21*x-21*y+15;
      d2psi[6][1][0] = d2psi[6][0][1];
      d2psi[6][1][1] = 18-27*y-21*x;
      d2psi[7][0][0] = 7-(25.0/2)*y-(17.0/2)*x;
      d2psi[7][0][1] = -(23.0/2)*y+19.0/2-(25.0/2)*x;
      d2psi[7][1][0] = d2psi[7][0][1];
      d2psi[7][1][1] = -(23.0/2)*y-(23.0/2)*x+9;
      d2psi[8][0][0] = -2+(7.0/2)*y+(5.0/2)*x;
      d2psi[8][0][1] = (7.0/2)*x+(7.0/2)*y-5.0/2;
      d2psi[8][1][0] = d2psi[8][0][1];
      d2psi[8][1][1] = (7.0/2)*x-4+(17.0/2)*y;
      d2psi[9][0][0] = -16*x+12-12*y;
      d2psi[9][0][1] = -12*x+8-8*y;
      d2psi[9][1][0] = d2psi[9][0][1];
      d2psi[9][1][1] = -8*x+4-4*y;
      d2psi[10][0][0] = 2*sqrt2*(-4+7*y+5*x);
      d2psi[10][0][1] = 2*sqrt2*(7*y+7*x-5);
      d2psi[10][1][0] = d2psi[10][0][1];
      d2psi[10][1][1] = 2*sqrt2*(5*y-4+7*x);
      d2psi[11][0][0] = -4*x+4-8*y;
      d2psi[11][0][1] = -8*x+8-12*y;
      d2psi[11][1][0] = d2psi[11][0][1];
      d2psi[11][1][1] = -12*x+12-16*y;

      return;
    }
  else if((y < 2*(0.5-x) + eps) && (y > x-eps)) // triangle opposite veretx 1
    {
      psi[0] = -3*y2+2*y3-3*x2-3*x2*y+1+3*x3;
      psi[1] = -(1.0/6)*x*(-6+7*x2+18*y+3*x-12*x*y-12*y2);
      psi[2] = (7.0/3)*x3+y3-(3.0/2)*x2-2*y2-(1.0/2)*x2*y+y;
      psi[3] = -(1.0/2)*x2*(-6-3*y+5*x);
      psi[4] = (1.0/12)*x2*(13*x-12-3*y);
      psi[5] = (1.0/12)*x2*(5*x-6+27*y);
      psi[6] = -(1.0/2)*x3+(3.0/2)*x2*y-2*y3+3*y2;
      psi[7] = -(1.0/12)*x*(13*x2-24*y2-6*x-21*x*y+12*y);
      psi[8] = (1.0/12)*x3-y2+y3-(1.0/4)*x2*y;
      psi[9] = (2.0/3)*x2*(4*x+3*y-3);
      psi[10] = (1.0/3)*sqrt2*(x-3*y)*x2;
      psi[11] = -(2.0/3)*x*(5*x2+6*y-6*x*y-3*x-6*y2);

      if(dpsi == NULL) return;

      dpsi[0][0] = 3*x*(-2-2*y+3*x);
      dpsi[0][1] = -6*y+6*y2-3*x2;
      dpsi[1][0] = 1-(7.0/2)*x2-3*y-x+4*x*y+2*y2;
      dpsi[1][1] = x*(-3+2*x+4*y);
      dpsi[2][0] = x*(7*x-3-y);
      dpsi[2][1] = 3*y2-4*y-(1.0/2)*x2+1;
      dpsi[3][0] = -(3.0/2)*x*(-4-2*y+5*x);
      dpsi[3][1] = (3.0/2)*x2;
      dpsi[4][0] = (1.0/4)*x*(13*x-8-2*y);
      dpsi[4][1] = -(1.0/4)*x2;
      dpsi[5][0] = (1.0/4)*x*(5*x-4+18*y);
      dpsi[5][1] = (9.0/4)*x2;
      dpsi[6][0] = -(3.0/2)*x*(x-2*y);
      dpsi[6][1] = (3.0/2)*x2-6*y2+6*y;
      dpsi[7][0] = -(13.0/4)*x2+2*y2+x+(7.0/2)*x*y-y;
      dpsi[7][1] = (1.0/4)*x*(16*y+7*x-4);
      dpsi[8][0] = (1.0/4)*x*(x-2*y);
      dpsi[8][1] = -2*y+3*y2-(1.0/4)*x2;
      dpsi[9][0] = 4*x*(2*x-1+y);
      dpsi[9][1] = 2*x2;
      dpsi[10][0] = sqrt2*x*(x-2*y);
      dpsi[10][1] = -sqrt2*x2;
      dpsi[11][0] = -10*x2-4*y+8*x*y+4*x+4*y2;
      dpsi[11][1] = 4*x*(x-1+2*y);

      if(d2psi == NULL) return;

      d2psi[0][0][0] = -6-6*y+18*x;
      d2psi[0][0][1] = -6*x;
      d2psi[0][1][0] = d2psi[0][0][1];
      d2psi[0][1][1] = -6+12*y;
      d2psi[1][0][0] = -7*x-1+4*y;
      d2psi[1][0][1] = -3+4*x+4*y;
      d2psi[1][1][0] = d2psi[1][0][1];
      d2psi[1][1][1] = 4*x;
      d2psi[2][0][0] = 14*x-3-y;
      d2psi[2][0][1] = -x;
      d2psi[2][1][0] = d2psi[2][0][1];
      d2psi[2][1][1] = 6*y-4;
      d2psi[3][0][0] = 6+3*y-15*x;
      d2psi[3][0][1] = 3*x;
      d2psi[3][1][0] = d2psi[3][0][1];
      d2psi[3][1][1] = 0;
      d2psi[4][0][0] = (13.0/2)*x-2-(1.0/2)*y;
      d2psi[4][0][1] = -(1.0/2)*x;
      d2psi[4][1][0] = d2psi[4][0][1];
      d2psi[4][1][1] = 0;
      d2psi[5][0][0] = (5.0/2)*x-1+(9.0/2)*y;
      d2psi[5][0][1] = (9.0/2)*x;
      d2psi[5][1][0] = d2psi[5][0][1];
      d2psi[5][1][1] = 0;
      d2psi[6][0][0] = -3*x+3*y;
      d2psi[6][0][1] = 3*x;
      d2psi[6][1][0] = d2psi[6][0][1];
      d2psi[6][1][1] = -12*y+6;
      d2psi[7][0][0] = -(13.0/2)*x+1+(7.0/2)*y;
      d2psi[7][0][1] = 4*y+(7.0/2)*x-1;
      d2psi[7][1][0] = d2psi[7][0][1];
      d2psi[7][1][1] = 4*x;
      d2psi[8][0][0] = (1.0/2)*x-(1.0/2)*y;
      d2psi[8][0][1] = -(1.0/2)*x;
      d2psi[8][1][0] = d2psi[8][0][1];
      d2psi[8][1][1] = -2+6*y;
      d2psi[9][0][0] = 16*x+4*y-4;
      d2psi[9][0][1] = 4*x;
      d2psi[9][1][0] = d2psi[9][0][1];
      d2psi[9][1][1] = 0;
      d2psi[10][0][0] = 2*sqrt2*(x-y);
      d2psi[10][0][1] = -2*sqrt2*x;
      d2psi[10][1][0] = d2psi[10][0][1];
      d2psi[10][1][1] = 0;
      d2psi[11][0][0] = -20*x+8*y+4;
      d2psi[11][0][1] = -4+8*x+8*y;
      d2psi[11][1][0] = d2psi[11][0][1];
      d2psi[11][1][1] = 8*x;

      return;
    }
  else  // triangle opposite vertex 2
    {
      assert(y < x+eps);
      assert(y < 0.5*(1-x)+eps);
       
      psi[0] = -3*x2+3*y3-3*x*y2-3*y2+1+2*x3;
      psi[1] = x3+x-(3.0/2)*y2-2*x2+(7.0/3)*y3-(1.0/2)*x*y2;
      psi[2] = (1.0/6)*y*(6+12*x2-3*y-7*y2-18*x+12*x*y);
      psi[3] = 3*x2-(1.0/2)*y3-2*x3+(3.0/2)*x*y2;
      psi[4] = -(1.0/4)*x*y2-x2+(1.0/12)*y3+x3;
      psi[5] = (1.0/12)*y*(21*x*y+6*y-13*y2+24*x2-12*x);
      psi[6] = (1.0/2)*y2*(-5*y+3*x+6);
      psi[7] = (1.0/12)*y2*(5*y+27*x-6);
      psi[8] = -(1.0/12)*y2*(-13*y+3*x+12);
      psi[9] = (2.0/3)*y*(-5*y2-6*x+6*x*y+3*y+6*x2);
      psi[10] = -(1.0/3)*sqrt2*(-y+3*x)*y2;
      psi[11] = (2.0/3)*y2*(4*y+3*x-3);

      if(dpsi == NULL) return;

      dpsi[0][0] = -6*x-3*y2+6*x2;
      dpsi[0][1] = -3*y*(-3*y+2*x+2);
      dpsi[1][0] = 3*x2+1-4*x-(1.0/2)*y2;
      dpsi[1][1] = -y*(-7*y+3+x);
      dpsi[2][0] = y*(4*x-3+2*y);
      dpsi[2][1] = 1+2*x2-y-(7.0/2)*y2-3*x+4*x*y;
      dpsi[3][0] = 6*x-6*x2+(3.0/2)*y2;
      dpsi[3][1] = (3.0/2)*y*(-y+2*x);
      dpsi[4][0] = -(1.0/4)*y2-2*x+3*x2;
      dpsi[4][1] = -(1.0/4)*y*(-y+2*x);
      dpsi[5][0] = (1.0/4)*y*(7*y+16*x-4);
      dpsi[5][1] = (7.0/2)*x*y+y-(13.0/4)*y2+2*x2-x;
      dpsi[6][0] = (3.0/2)*y2;
      dpsi[6][1] = (3.0/2)*y*(-5*y+2*x+4);
      dpsi[7][0] = (9.0/4)*y2;
      dpsi[7][1] = (1.0/4)*y*(5*y+18*x-4);
      dpsi[8][0] = -(1.0/4)*y2;
      dpsi[8][1] = -(1.0/4)*y*(-13*y+2*x+8);
      dpsi[9][0] = 4*y*(2*x-1+y);
      dpsi[9][1] = -10*y2-4*x+8*x*y+4*y+4*x2;
      dpsi[10][0] = -sqrt2*y2;
      dpsi[10][1] = -sqrt2*(-y+2*x)*y;
      dpsi[11][0] = 2*y2;
      dpsi[11][1] = 4*y*(x-1+2*y);

      if(d2psi == NULL) return;

      d2psi[0][0][0] = -6+12*x;
      d2psi[0][0][1] = -6*y;
      d2psi[0][1][0] = d2psi[0][0][1];
      d2psi[0][1][1] = 18*y-6*x-6;
      d2psi[1][0][0] = 6*x-4;
      d2psi[1][0][1] = -y;
      d2psi[1][1][0] = d2psi[1][0][1];
      d2psi[1][1][1] = -3+14*y-x;
      d2psi[2][0][0] = 4*y;
      d2psi[2][0][1] = -3+4*x+4*y;
      d2psi[2][1][0] = d2psi[2][0][1];
      d2psi[2][1][1] = -1-7*y+4*x;
      d2psi[3][0][0] = 6-12*x;
      d2psi[3][0][1] = 3*y;
      d2psi[3][1][0] = d2psi[3][0][1];
      d2psi[3][1][1] = -3*y+3*x;
      d2psi[4][0][0] = -2+6*x;
      d2psi[4][0][1] = -(1.0/2)*y;
      d2psi[4][1][0] = d2psi[4][0][1];
      d2psi[4][1][1] = -(1.0/2)*x+(1.0/2)*y;
      d2psi[5][0][0] = 4*y;
      d2psi[5][0][1] = (7.0/2)*y+4*x-1;
      d2psi[5][1][0] = d2psi[5][0][1];
      d2psi[5][1][1] = (7.0/2)*x+1-(13.0/2)*y;
      d2psi[6][0][0] = 0;
      d2psi[6][0][1] = 3*y;
      d2psi[6][1][0] = d2psi[6][0][1];
      d2psi[6][1][1] = -15*y+3*x+6;
      d2psi[7][0][0] = 0;
      d2psi[7][0][1] = (9.0/2)*y;
      d2psi[7][1][0] = d2psi[7][0][1];
      d2psi[7][1][1] = (5.0/2)*y+(9.0/2)*x-1;
      d2psi[8][0][0] = 0;
      d2psi[8][0][1] = -(1.0/2)*y;
      d2psi[8][1][0] = d2psi[8][0][1];
      d2psi[8][1][1] = (13.0/2)*y-(1.0/2)*x-2;
      d2psi[9][0][0] = 8*y;
      d2psi[9][0][1] = -4+8*x+8*y;
      d2psi[9][1][0] = d2psi[9][0][1];
      d2psi[9][1][1] = -20*y+8*x+4;
      d2psi[10][0][0] = 0;
      d2psi[10][0][1] = -2*sqrt2*y;
      d2psi[10][1][0] = d2psi[10][0][1];
      d2psi[10][1][1] = -2*sqrt2*(x-y);
      d2psi[11][0][0] = 0;
      d2psi[11][0][1] = 4*y;
      d2psi[11][1][0] = d2psi[11][0][1];
      d2psi[11][1][1] = 16*y+4*x-4;
     
      return;
    }

  return;
}

void cloughTocherTriangleBdy(svector<1> xi, double psi[4], double psin[3]) 
{

  double x = xi[0];
  double omx=1-x, omxsq=omx*omx;
  double opx=1+x, opxsq=opx*opx;
  
  psi[0] =  .25*(x+2)*omxsq;
  psi[1] =  .25*opx*omxsq;
  psi[2] = -.25*(x-2)*opxsq;
  psi[3] = -.25*omx*opxsq;

  psin[0] = -x * omx / 2; 
  psin[1] =  x * opx / 2; 
  psin[2] =  opx * omx; 

  return; 
}

void cloughTocherBdyTransform(int *nodes, sbmatrix<2> &dxds,
		      double Phi[4], double Phin[3], 
		      double phi[7], double phin[7])
{
  int i, j, n0 = 0, ii = 0, i0 = 0;
  double det;

  svector<2> normal = dxds.normal(det);

  for(i = 0; i < 2; i++)     // for each vertex
    {
    phi[ii] = Phi[i0];
    phin[ii] = 0.0;

    ii ++;
    i0 ++;
    
    for(j = 0; j < 2; j++)   // gradient dof
      {
      phi[ii] = dxds[j][0] * Phi[i0];
      
      phin[ii] = normal[j] * Phin[n0];
      
      ii ++;
      }
    
    i0 ++;
    n0 ++;
    }

  assert(ii == 6);
  assert(i0 == 4);
  assert(n0 == 2);

  phi[ii]  = 0.0;         // function values don't depend upon du/dn
  phin[ii] = Phin[n0];    // normal derivative at the middle of the edge

  if(nodes[0] > nodes[1]) phin[ii] *= -1;

  return;
}

void cloughTocherTransform(int *nodes, smatrix<2> &dxds,
		   double Phi[12], svector<2> Dphi[12], smatrix<2> D2phi[12],
		   double phi[12], svector<2> dphi[12], smatrix<2> d2phi[12])
{
  const int ndim=2;
  int i,j, alpha, ii, i0=0;

  double det, modFtn, psi[2];
  double a0, a1;
  svector<2> nn[3], dpsi[2];
  smatrix<2> d2psi[2];

  svector<2> dphiMid[9][2] = 
    {{svector<2>(-1.50,  0),     svector<2>( 0,    -1.5)},
     {svector<2>(-0.25,  0),     svector<2>( 0,     0)},
     {svector<2>( 0,     0),     svector<2>( 0,    -0.25)},
     {svector<2>( 0.75, -0.75),  svector<2>( 1.5,   0)},
     {svector<2>(-0.125, 0.125), svector<2>(-0.25,  0)},
     {svector<2>( 0.125,-0.125), svector<2>( 0,     0)},
     {svector<2>( 0,     1.5),   svector<2>(-0.75,  0.75)},
     {svector<2>( 0,     0),     svector<2>(-0.125, 0.125)},
     {svector<2>( 0,    -0.25),  svector<2>( 0.125,-0.125)}};
 
  smatrix<2> dsdx  = dxds.inverse(det);
  smatrix<2> dsdxT = dsdx.transpose();
  
  nn[0] = dsdxT * svector<2>(0.0,-1.0);
  nn[1] = dsdxT * svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0));
  nn[2] = dsdxT * svector<2>(-1.0,0.0);
  
  for(i = 0; i < 3; i++)
    {
       ii = 9 + i;

       modFtn = 1.0/nn[i].norm();

       if(nodes[i] > nodes[(i+1)%3]) modFtn *= -1;

       phi[ii]   = modFtn * Phi[ii];
       dphi[ii]  = modFtn * dsdxT * Dphi[ii];
       d2phi[ii] = modFtn * dsdxT * D2phi[ii] * dsdx;

       nn[i] = (modFtn*modFtn) * nn[i];
     }
 
  // currently Hermite elements must be affine equivalent to the master

   ii = 0;

   for(i = 0; i < 3; i++)
     {
     j = (i+2)%3;    // edge before node i

     a0 = nn[i].dot(dsdxT * dphiMid[ii][0]);
     a1 = nn[j].dot(dsdxT * dphiMid[ii][1]);

     phi[ii]  =          Phi[i0]   - a0 * Phi[9+i]   - a1 * Phi[9+j];
     dphi[ii] =  dsdxT* (Dphi[i0]  - a0 * Dphi[9+i]  - a1 * Dphi[9+j]);
     d2phi[ii] = dsdxT* (D2phi[i0] - a0 * D2phi[9+i] - a1 * D2phi[9+j])*dsdx;

     ii ++;
     i0 ++;

     for(alpha = 0; alpha < ndim; alpha++)    // gradient dof
       {
	 a0 = nn[i].dot(dsdxT * dphiMid[i0+alpha][0]);
	 a1 = nn[j].dot(dsdxT * dphiMid[i0+alpha][1]);

	 psi[alpha]   = Phi[i0+alpha]   - a0 * Phi[9+i]   - a1 * Phi[9+j];
	 dpsi[alpha]  = Dphi[i0+alpha]  - a0 * Dphi[9+i]  - a1 * Dphi[9+j];
	 d2psi[alpha] = D2phi[i0+alpha] - a0 * D2phi[9+i] - a1 * D2phi[9+j];
       }

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       for(alpha = 0; alpha < ndim; alpha++)
	 {
	   phi[ii]   += dxds[j][alpha] * psi[alpha];
	   dphi[ii]  += dxds[j][alpha] * dpsi[alpha];
	   d2phi[ii] += dxds[j][alpha] * d2psi[alpha];
	 }

       dphi[ii]  = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;

       ii ++;
       }

     i0 += ndim;
     }

   assert(ii == 9);
   assert(i0 == 9);

  return;
}

void cloughTocherBdyInterpolate(int *nodes, svector<2> &nn,
			double uu[3], svector<2> du[3], smatrix<2> d2u[3],
			double ue[7])
{
  const int ndim=2;
  int i,j, ii=0;
 
   for(i = 0; i < 2; i++)
     {
     ue[ii] = uu[i];

     ii ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       ue[ii] = du[i][j];

       ii ++;
       }
     }

   assert(ii == 6);

   ue[ii] = nn.dot(du[2]);           // du/dn at the mid point

   if(nodes[0] > nodes[1]) ue[ii] *= -1;

  return;
}

void cloughTocherInterpolate(int *nodes, smatrix<2> &dxds,
		     double uu[6], svector<2> du[6], smatrix<2> d2u[6],
		     double ue[12])
{
  const int ndim=2;
  int i,j, ii=0;
 
   for(i = 0; i < 3; i++)
     {
     ue[ii] = uu[i];

     ii ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       ue[ii] = du[i][j];

       ii ++;
       }
     }

   assert(ii == 9);

   // compute the normal derivative at the edge centers 
   // assuming the map from the parent is affine

   double det;

   smatrix<2> dsdx  = dxds.inverse(det);
   smatrix<2> dsdxT = dsdx.transpose();

   svector<2> nn[3];

   nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1.0/ nn[0].norm()) * nn[0];
   nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1.0/ nn[1].norm()) * nn[1];
   nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1.0/ nn[2].norm()) * nn[2];
  
   for(i = 0; i < 3; i++) 
     {
       ue[ii+i] = du[3+i].dot(nn[i]);

       if(nodes[i] > nodes[(i+1)%3]) ue[9 + i] *= -1;
     }

  return;
}

// Bell Triangle
// Bell Triangle
// Bell Triangle

void bellTriangle(svector<2> xi, double psi[18], 
		  svector<2> *dpsi, smatrix<2> *d2psi)
{
  double x = xi[0], y = xi[1];

  double omxmy = 1-x-y,  omxmysq=omxmy*omxmy;
  double  x2=x*x, y2=y*y, x3=x2*x, y3=y2*y;

  // The following basis functions were generated automatically by Maple

  psi[0] = -(6*x3-3*x2-12*x2*y-12*x*y2-2*x-6*x*y-3*y2+6*y3-1-2*y)*omxmysq;
  psi[1] = -x*(1+3*x)*(x-1-2*y)*omxmysq;
  psi[2] = y*(1+3*y)*(2*x+1-y)*omxmysq;
  psi[3] = -.5*x2*(x-1-2*y)*omxmysq;
  psi[4] = x*y*omxmysq;
  psi[5] = .5*y2*(2*x+1-y)*omxmysq;
  psi[6] = x2*(10*x+15*y2-15*x2-15*y3-15*x*y2+6*x3);
  psi[7] = -.5*x2*(8*x+15*y2-14*x2-15*y3-15*x*y2+6*x3);
  psi[8] = .5*x2*y*(-3*y2+3*x*y-3*y+6-4*x);
  psi[9] = .25*x2*(2*x+5*y2-4*x2-5*y3-5*x*y2+2*x3);
  psi[10] = -.5*x2*y*(-y2-y+x*y+2-2*x);
  psi[11] = .25*x2*y2*(-y+1+x);
  psi[12] = -y2*(-10*y+15*y2-15*x2-6*y3+15*x3+15*x2*y);
  psi[13] = -.5*x*y2*(-6+4*y+3*x-3*x*y+3*x2);
  psi[14] = .5*y2*(-8*y+14*y2-15*x2-6*y3+15*x3+15*x2*y);
  psi[15] = -.25*x2*y2*(-y-1+x);
  psi[16] = .5*x*y2*(-2+2*y+x-x*y+x2);
  psi[17] = -.25*y2*(-2*y+4*y2-5*x2-2*y3+5*x3+5*x2*y);

  // Derivatives

  if(dpsi == NULL) return;

  dpsi[0][0] = 30*x*omxmy*(x2-x*y-x-2*y2);
  dpsi[0][1] = -30*y*omxmy*(2*x2+x*y+y-y2);
  dpsi[1][0] = omxmy*(15*x3-17*x2-15*x2*y-12*x*y2+x+2*x*y-2*y2+1+y);
  dpsi[1][1] = -6*x*y*(1+3*x)*omxmy;
  dpsi[2][0] = -6*x*y*(1+3*y)*omxmy;
  dpsi[2][1] = -omxmy*(12*x2*y+2*x2-2*x*y+15*x*y2-x-y+17*y2-15*y3-1);
  dpsi[3][0] = .5*x*omxmy*(5*x2-7*x-5*x*y-4*y2+2+2*y);
  dpsi[3][1] = -3*x2*y*omxmy;
  dpsi[4][0] = -y*omxmy*(3*x-1+y);
  dpsi[4][1] = -x*(x-1+3*y)*omxmy;
  dpsi[5][0] = -3*x*y2*omxmy;
  dpsi[5][1] = -.5*y*omxmy*(4*x2-2*x+5*x*y-2+7*y-5*y2);
  dpsi[6][0] = 15*x*(2*x+2*y2-4*x2-2*y3-3*x*y2+2*x3);
  dpsi[6][1] = -15*x2*y*(-2+2*x+3*y);
  dpsi[7][0] = -.5*x*(24*x+30*y2-56*x2-30*y3-45*x*y2+30*x3);
  dpsi[7][1] = 7.5*x2*y*(-2+2*x+3*y);
  dpsi[8][0] = 1.5*x*y*(-2*y2-2*y+3*x*y-4*x+4);
  dpsi[8][1] = .5*x2*(-9*y2+6*x*y-6*y+6-4*x);
  dpsi[9][0] = .25*x*(6*x+10*y2-16*x2-10*y3-15*x*y2+10*x3);
  dpsi[9][1] = -1.25*x2*y*(-2+2*x+3*y);
  dpsi[10][0] = -.5*x*y*(-2*y2-2*y+3*x*y-6*x+4);
  dpsi[10][1] = -.5*x2*(-3*y2-2*y+2*x*y+2-2*x);
  dpsi[11][0] = .25*x*y2*(-2*y+2+3*x);
  dpsi[11][1] = .25*x2*y*(-3*y+2*x+2);
  dpsi[12][0] = -15*y2*x*(-2+3*x+2*y);
  dpsi[12][1] = -15*y*(-2*y+4*y2-2*x2-2*y3+2*x3+3*x2*y);
  dpsi[13][0] = -.5*y2*(-6+4*y+6*x-6*x*y+9*x2);
  dpsi[13][1] = -1.5*x*y*(-4+4*y+2*x-3*x*y+2*x2);
  dpsi[14][0] = 7.5*y2*x*(-2+3*x+2*y);
  dpsi[14][1] = .5*y*(-24*y+56*y2-30*x2-30*y3+30*x3+45*x2*y);
  dpsi[15][0] = -.25*x*y2*(-2*y-2+3*x);
  dpsi[15][1] = -.25*x2*y*(-3*y-2+2*x);
  dpsi[16][0] = .5*y2*(-2+2*y+2*x-2*x*y+3*x2);
  dpsi[16][1] = .5*x*y*(-4+6*y+2*x-3*x*y+2*x2);
  dpsi[17][0] = -1.25*y2*x*(-2+3*x+2*y);
  dpsi[17][1] = -.25*y*(-6*y+16*y2-10*x2-10*y3+10*x3+15*x2*y);

  // Second Derivatives

  if(d2psi == NULL) return;

  d2psi[0][0][0] = -60*x+180*x*y2+60*y3-120*x3-60*y2+180*x2;
  d2psi[0][1][0] = 60*x*y*(3*x+3*y-2);
  d2psi[0][0][1] = d2psi[0][1][0];
  d2psi[0][1][1] = -60*y+180*x2*y-120*y3+60*x3+180*y2-60*x2;
  d2psi[1][0][0] = -36*x+54*x*y2+12*y3-60*x3-12*y2+96*x2;
  d2psi[1][1][0] = 6*y*(-1+9*x2+6*x*y-4*x+y);
  d2psi[1][0][1] = d2psi[1][1][0];
  d2psi[1][1][1] = 6*x*(1+3*x)*(x-1+2*y);
  d2psi[2][0][0] = 6*y*(1+3*y)*(y-1+2*x);
  d2psi[2][1][0] = 6*x*(6*x*y+9*y2-4*y+x-1);
  d2psi[2][0][1] = d2psi[2][1][0];
  d2psi[2][1][1] = -36*y+54*x2*y-60*y3+12*x3+96*y2-12*x2;
  d2psi[3][0][0] = 1-9*x+9*x*y2+2*y3-10*x3-3*y2+18*x2;
  d2psi[3][1][0] = 3*x*y*(-2+3*x+2*y);
  d2psi[3][0][1] = d2psi[3][1][0];
  d2psi[3][1][1] = 3*x2*(x-1+2*y);
  d2psi[4][0][0] = 2*y*(-2+3*x+2*y);
  d2psi[4][1][0] = 1-4*x-4*y+3*x2+8*x*y+3*y2;
  d2psi[4][0][1] = d2psi[4][1][0];
  d2psi[4][1][1] = 2*x*(-2+2*x+3*y);
  d2psi[5][0][0] = 3*y2*(y-1+2*x);
  d2psi[5][1][0] = 3*x*y*(-2+2*x+3*y);
  d2psi[5][0][1] = d2psi[5][1][0];
  d2psi[5][1][1] = 1-9*y+9*x2*y-10*y3+2*x3+18*y2-3*x2;
  d2psi[6][0][0] = 60*x+30*y2-180*x2-30*y3-90*x*y2+120*x3;
  d2psi[6][1][0] = -30*x*y*(3*x+3*y-2);
  d2psi[6][0][1] = d2psi[6][1][0];
  d2psi[6][1][1] = -30*x2*(x-1+3*y);
  d2psi[7][0][0] = -24*x-15*y2+84*x2+15*y3+45*x*y2-60*x3;
  d2psi[7][1][0] = 15*x*y*(3*x+3*y-2);
  d2psi[7][0][1] = d2psi[7][1][0];
  d2psi[7][1][1] = 15*x2*(x-1+3*y);
  d2psi[8][0][0] = 3*y*(-y2+3*x*y-y-4*x+2);
  d2psi[8][1][0] = 3*x*(-3*y2+3*x*y-2*y+2-2*x);
  d2psi[8][0][1] = d2psi[8][1][0];
  d2psi[8][1][1] = 3*x2*(-3*y+x-1);
  d2psi[9][0][0] = 3*x+2.5*y2-12*x2-2.5*y3-7.5*x*y2+10*x3;
  d2psi[9][1][0] = -2.5*x*y*(3*x+3*y-2);
  d2psi[9][0][1] = d2psi[9][1][0];
  d2psi[9][1][1] = -2.5*x2*(x-1+3*y);
  d2psi[10][0][0] = -y*(-y2-y+3*x*y+2-6*x);
  d2psi[10][1][0] = -x*(-3*y2-2*y+3*x*y+2-3*x);
  d2psi[10][0][1] = d2psi[10][1][0];
  d2psi[10][1][1] = -x2*(-3*y+x-1);
  d2psi[11][0][0] = .5*y2*(-y+1+3*x);
  d2psi[11][1][0] = .5*x*y*(-3*y+2+3*x);
  d2psi[11][0][1] = d2psi[11][1][0];
  d2psi[11][1][1] = .5*x2*(-3*y+1+x);
  d2psi[12][0][0] = -30*y2*(3*x-1+y);
  d2psi[12][1][0] = -30*x*y*(3*x+3*y-2);
  d2psi[12][0][1] = d2psi[12][1][0];
  d2psi[12][1][1] = 60*y-180*y2+30*x2+120*y3-30*x3-90*x2*y;
  d2psi[13][0][0] = -3*y2*(-y+1+3*x);
  d2psi[13][1][0] = -3*y*(-2+2*y+2*x-3*x*y+3*x2);
  d2psi[13][0][1] = d2psi[13][1][0];
  d2psi[13][1][1] = -3*x*(-2+4*y+x-3*x*y+x2);
  d2psi[14][0][0] = 15*y2*(3*x-1+y);
  d2psi[14][1][0] = 15*x*y*(3*x+3*y-2);
  d2psi[14][0][1] = d2psi[14][1][0];
  d2psi[14][1][1] = -24*y+84*y2-15*x2-60*y3+15*x3+45*x2*y;
  d2psi[15][0][0] = -.5*y2*(-y-1+3*x);
  d2psi[15][1][0] = -.5*x*y*(-3*y-2+3*x);
  d2psi[15][0][1] = d2psi[15][1][0];
  d2psi[15][1][1] = -.5*x2*(-3*y+x-1);
  d2psi[16][0][0] = y2*(-y+1+3*x);
  d2psi[16][1][0] = y*(-2+3*y+2*x-3*x*y+3*x2);
  d2psi[16][0][1] = d2psi[16][1][0];
  d2psi[16][1][1] = x*(-2+6*y+x-3*x*y+x2);
  d2psi[17][0][0] = -2.5*y2*(3*x-1+y);
  d2psi[17][1][0] = -2.5*x*y*(3*x+3*y-2);
  d2psi[17][0][1] = d2psi[17][1][0];
  d2psi[17][1][1] = 3*y-12*y2+2.5*x2+10*y3-2.5*x3-7.5*x2*y;

  return;
}

void bellTriangleBdy(svector<1> xi, double psi[6], double psin[4]) 
{

  double x = xi[0], x2=x*x;
  double omx=1-x, omxsq=omx*omx;
  double opx=1+x, opxsq=opx*opx;
  
  psi[0] = .0625*(3*x2+9*x+8)*omxsq*omx;
  psi[1] = .0625*(3*x+5)*opx*omxsq*omx;
  psi[2] = .0625*opxsq*omxsq*omx;
  psi[3] = .0625*(3*x2-9*x+8)*opxsq*opx;
  psi[4] = .0625*omx*(3*x-5)*opxsq*opx;
  psi[5] = .0625*omxsq*opxsq*opx;
  
  psin[0] =  .25*(x+2)*omxsq;
  psin[1] =  .25*opx*omxsq;
  psin[2] = -.25*(x-2)*opxsq;
  psin[3] = -.25*omx*opxsq;
 
  return; 
}

void bellBdyTransform(int *nodes, sbmatrix<2> &dxds,
		      double Phi[6],  double Phin[4], 
		      double phi[12], double phin[12])
{
  int i, j, k, n0 = 0, ii = 0, i0 = 0;
  double det;

  svector<2> normal = dxds.normal(det);

  for(i = 0; i < 2; i++)     // for each vertex
    {
    phi[ii] = Phi[i0];
    phin[ii] = 0.0;

    ii ++;
    i0 ++;
    
    for(j = 0; j < 2; j++)   // gradient dof
      {
      phi[ii] = dxds[j][0] * Phi[i0];
      
      phin[ii] = normal[j] * Phin[n0];
      
      ii ++;
      }
    
    i0 ++;
    n0 ++;
    
    for(j = 0; j < 2; j++)   // Hessian dof
    for(k = j; k < 2; k++)   
      {
      phi[ii] = dxds[j][0] * dxds[k][0] * Phi[i0];

      phin[ii] = 0.5 *
	(dxds[j][0] * normal[k] + dxds[k][0] * normal[j]) * Phin[n0];

      ii++;
      }
    
    i0 ++;
    n0 ++;
    }

  assert(ii == 12);
  assert(i0 == 6);
  assert(n0 == 4);
  
  return;
}

void bellTransform(int *nodes, smatrix<2> &dxds,
		   double Phi[18], svector<2> Dphi[18], smatrix<2> D2phi[18],
		   double phi[18], svector<2> dphi[18], smatrix<2> d2phi[18])
{
  const int ndim=2;
  int i,j, k, alpha,beta, idxab, ii=0, i0=0;

  double temp, det;
  smatrix<2> dsdx  = dxds.inverse(det);
  smatrix<2> dsdxT = dsdx.transpose();
 
  // currently Hermite elements must be affine equivalent to the master

   for(i = 0; i < 3; i++)
     {
     phi[ii] = Phi[i0];
     dphi[ii] = dsdxT * Dphi[i0];
     d2phi[ii] = dsdxT * D2phi[i0] * dsdx;

     ii ++;
     i0 ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       for(alpha = 0; alpha < ndim; alpha++)
	 {
	   phi[ii]   += dxds[j][alpha] * Phi[i0+alpha];
	   dphi[ii]  += dxds[j][alpha] * Dphi[i0+alpha];
	   d2phi[ii] += dxds[j][alpha] * D2phi[i0+alpha];
	 }

       dphi[ii] = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;

       ii ++;
       }

     i0 += ndim;

     for(j = 0; j < ndim; j++)   // Hessian dof
     for(k = j; k < ndim; k++)   
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       idxab = i0;   // indexing for symmetric matrx 0 <= a <= b < ndim

       for(alpha = 0;    alpha < ndim; alpha++)
       for(beta  = alpha; beta < ndim; beta++)
	 {
	 temp = 0.5 *
	   (dxds[j][alpha] * dxds[k][beta] + dxds[k][alpha] * dxds[j][beta]);

	 phi[ii]   += temp * Phi[idxab];
	 dphi[ii]  += temp * Dphi[idxab];
	 d2phi[ii] += temp * D2phi[idxab];

	 idxab ++;
	 }

       dphi[ii] = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;
       ii++;
       }

     i0 += (ndim*(ndim+1))/2;
     }

   assert(ii == 18);
   assert(i0 == 18);

  return;
}

void bellBdyInterpolate(int *nodes,
			double uu[2], svector<2> du[2], smatrix<2> d2u[2],
			double ue[12])
{
  const int ndim=2;
  int i,j,k, ii=0;
 
   for(i = 0; i < 2; i++)
     {
     ue[ii] = uu[i];

     ii ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       ue[ii] = du[i][j];

       ii ++;
       }

     for(j = 0; j < ndim; j++)   // Hessian dof
     for(k = j; k < ndim; k++)   
       {
       ue[ii] = d2u[i][j][k];

       if(j != k) ue[ii] *= 2.0;

       ii++;
       }
     }

   assert(ii == 12);

  return;
}

void bellInterpolate(int *nodes,
		     double uu[3], svector<2> du[3], smatrix<2> d2u[3],
		     double ue[18])
{
  const int ndim=2;
  int i,j,k, ii=0;
 
   for(i = 0; i < 3; i++)
     {
     ue[ii] = uu[i];

     ii ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       ue[ii] = du[i][j];

       ii ++;
       }

     for(j = 0; j < ndim; j++)   // Hessian dof
     for(k = j; k < ndim; k++)   
       {
       ue[ii] = d2u[i][j][k];

       if(j != k) ue[ii] *= 2.0;

       ii++;
       }
     }

   assert(ii == 18);

  return;
}

// Argyris Triangle
// Argyris Triangle
// Argyris Triangle

void argyrisTriangle(svector<2> xi, double psi[21], 
		     svector<2> *dpsi, smatrix<2> *d2psi)
{
  double x = xi[0], y = xi[1];

  double omxmy = 1-x-y,  omxmysq=omxmy*omxmy;
  double  x2=x*x, y2=y*y, x3=x2*x, y3=y2*y, sqrt2=sqrt(2.0);

  // The following basis functions were generated automatically by Maple

  psi[0] = -(6*x3-3*x2-12*x2*y-12*x*y2-2*x-6*x*y-3*y2+6*y3-1-2*y)*omxmysq;
  psi[1] = -x*(3*x2-2*x-6*x*y+8*y2-1-2*y)*omxmysq;
  psi[2] = -y*(8*x2-2*x-6*x*y-2*y-1+3*y2)*omxmysq;
  psi[3] = -.5*x2*(x-1-2*y)*omxmysq;
  psi[4] = -x*y*(2*x-1+2*y)*omxmysq;
  psi[5] = .5*y2*(2*x+1-y)*omxmysq;
  psi[6] = x2*(10*x+15*y2-15*x2-15*y3-15*x*y2+6*x3);
  psi[7] = -.5*x2*(8*x+7*y2-14*x2-7*y3-7*x*y2+6*x3);
  psi[8] = -.5*x2*y*(10-37*y-28*x+27*y2+37*x*y+16*x2);
  psi[9] = .25*x2*(2*x+y2-4*x2-y3-x*y2+2*x3);
  psi[10] = .5*x2*y*(2-7*y-6*x+5*y2+7*x*y+4*x2);
  psi[11] = -.25*x2*y2*(5*y-5+3*x);
  psi[12] = -y2*(-10*y+15*y2-15*x2-6*y3+15*x2*y+15*x3);
  psi[13] = -.5*x*y2*(10-28*y-37*x+16*y2+37*x*y+27*x2);
  psi[14] = .5*y2*(-8*y+14*y2-7*x2-6*y3+7*x2*y+7*x3);
  psi[15] = -.25*x2*y2*(3*y-5+5*x);
  psi[16] = .5*x*y2*(2-6*y-7*x+4*y2+7*x*y+5*x2);
  psi[17] = -.25*y2*(-2*y+4*y2-x2-2*y3+x2*y+x3);
  psi[18] = -16*x2*y*omxmysq;
  psi[19] = -8*sqrt2*omxmy*x2*y2;
  psi[20] = -16*x*y2*omxmysq;

  // Derivatives

  if(dpsi == NULL) return;

  dpsi[0][0] = 30*x*omxmy*(x2-x-x*y-2*y2);
  dpsi[0][1] = -30*y*omxmy*(2*x2+x*y+y-y2);
  dpsi[1][0] = omxmy*(15*x3-17*x2-15*x2*y+x+12*x*y2+2*x*y-10*y2+8*y3+1+y);
  dpsi[1][1] = -2*x*y*omxmy*(x+11-16*y);
  dpsi[2][0] = 2*x*y*omxmy*(16*x-11-y);
  dpsi[2][1] = omxmy*(8*x3-10*x2+12*x2*y+x-15*x*y2+2*x*y+1+y-17*y2+15*y3);
  dpsi[3][0] = .5*x*omxmy*(5*x2-7*x-5*x*y+2-4*y2+2*y);
  dpsi[3][1] = -3*x2*y*omxmy;
  dpsi[4][0] = y*omxmy*(8*x2-7*x+10*x*y-3*y+1+2*y2);
  dpsi[4][1] = x*omxmy*(2*x2-3*x+10*x*y-7*y+1+8*y2);
  dpsi[5][0] = -3*x*y2*omxmy;
  dpsi[5][1] = -.5*y*omxmy*(4*x2-2*x+5*x*y-2+7*y-5*y2);
  dpsi[6][0] = 15*x*(2*x+2*y2-4*x2-2*y3-3*x*y2+2*x3);
  dpsi[6][1] = -15*x2*y*(-2+3*y+2*x);
  dpsi[7][0] = -.5*x*(24*x+14*y2-56*x2-14*y3-21*x*y2+30*x3);
  dpsi[7][1] = 3.5*x2*y*(-2+3*y+2*x);
  dpsi[8][0] = -.5*x*y*(20-74*y-84*x+54*y2+111*x*y+64*x2);
  dpsi[8][1] = -.5*x2*(10-74*y-28*x+81*y2+74*x*y+16*x2);
  dpsi[9][0] = .25*x*(6*x+2*y2-16*x2-2*y3-3*x*y2+10*x3);
  dpsi[9][1] = -.25*x2*y*(-2+3*y+2*x);
  dpsi[10][0] = .5*x*y*(4-14*y-18*x+10*y2+21*x*y+16*x2);
  dpsi[10][1] = .5*x2*(2-14*y-6*x+15*y2+14*x*y+4*x2);
  dpsi[11][0] = -.25*x*y2*(10*y-10+9*x);
  dpsi[11][1] = -.25*x2*y*(15*y-10+6*x);
  dpsi[12][0] = -15*y2*x*(-2+2*y+3*x);
  dpsi[12][1] = -15*y*(-2*y+4*y2-2*x2-2*y3+3*x2*y+2*x3);
  dpsi[13][0] = -.5*y2*(10-28*y-74*x+16*y2+74*x*y+81*x2);
  dpsi[13][1] = -.5*x*y*(20-84*y-74*x+64*y2+111*x*y+54*x2);
  dpsi[14][0] = 3.5*y2*x*(-2+2*y+3*x);
  dpsi[14][1] = .5*y*(-24*y+56*y2-14*x2-30*y3+21*x2*y+14*x3);
  dpsi[15][0] = -.25*x*y2*(6*y-10+15*x);
  dpsi[15][1] = -.25*x2*y*(9*y-10+10*x);
  dpsi[16][0] = .5*y2*(2-6*y-14*x+4*y2+14*x*y+15*x2);
  dpsi[16][1] = .5*x*y*(4-18*y-14*x+16*y2+21*x*y+10*x2);
  dpsi[17][0] = -.25*y2*x*(-2+2*y+3*x);
  dpsi[17][1] = -.25*y*(-6*y+16*y2-2*x2-10*y3+3*x2*y+2*x3);
  dpsi[18][0] = 32*x*y*omxmy*(2*x-1+y);
  dpsi[18][1] = 16*x2*(x-1+3*y)*omxmy;
  dpsi[19][0] = 8*sqrt2*(-2+2*y+3*x)*x*y2;
  dpsi[19][1] = 8*sqrt2*(-2+3*y+2*x)*y*x2;
  dpsi[20][0] = 16*y2*omxmy*(3*x-1+y);
  dpsi[20][1] = 32*x*y*(x-1+2*y)*omxmy;

  // Second Derivatives

  if(d2psi == NULL) return;

  d2psi[0][0][0] = -60*x-60*y2+180*x2+60*y3-120*x3+180*x*y2;
  d2psi[0][0][1] = 60*x*y*(3*x-2+3*y);
  d2psi[0][1][0] = d2psi[0][0][1];
  d2psi[0][1][1] = -60*y+180*y2-60*x2-120*y3+60*x3+180*x2*y;
  d2psi[1][0][0] = -36*x+20*y2+96*x2-20*y3-60*x3+6*x*y2;
  d2psi[1][0][1] = 2*y*(-11+27*y-16*y2+20*x-30*x*y+3*x2);
  d2psi[1][1][0] = d2psi[1][0][1];
  d2psi[1][1][1] = 2*x*(-11+10*x+54*y+x2-30*x*y-48*y2);
  d2psi[2][0][0] = -2*y*(11-54*x-10*y+48*x2+30*x*y-y2);
  d2psi[2][0][1] = -2*x*(11-27*x+16*x2-20*y-3*y2+30*x*y);
  d2psi[2][1][0] = d2psi[2][0][1];
  d2psi[2][1][1] = -36*y+96*y2+20*x2-60*y3-20*x3+6*x2*y;
  d2psi[3][0][0] = 1-9*x-3*y2+18*x2+2*y3-10*x3+9*x*y2;
  d2psi[3][0][1] = 3*x*y*(-2+2*y+3*x);
  d2psi[3][1][0] = d2psi[3][0][1];
  d2psi[3][1][1] = 3*x2*(x-1+2*y);
  d2psi[4][0][0] = -2*y*(4-10*y-15*x+6*y2+18*x*y+12*x2);
  d2psi[4][0][1] = 1-8*x-8*y+15*y2+15*x2-8*y3-8*x3+40*x*y-36*x*y2-36*x2*y;
  d2psi[4][1][0] = d2psi[4][0][1];
  d2psi[4][1][1] = -2*x*(4-10*x-15*y+6*x2+18*x*y+12*y2);
  d2psi[5][0][0] = 3*y2*(2*x-1+y);
  d2psi[5][0][1] = 3*x*y*(-2+3*y+2*x);
  d2psi[5][1][0] = d2psi[5][0][1];
  d2psi[5][1][1] = 1-9*y+18*y2-3*x2-10*y3+2*x3+9*x2*y;
  d2psi[6][0][0] = 60*x+30*y2-180*x2-30*y3-90*x*y2+120*x3;
  d2psi[6][0][1] = -30*x*y*(3*x-2+3*y);
  d2psi[6][1][0] = d2psi[6][0][1];
  d2psi[6][1][1] = -30*x2*(x-1+3*y);
  d2psi[7][0][0] = -24*x-7*y2+84*x2+7*y3+21*x*y2-60*x3;
  d2psi[7][0][1] = 7*x*y*(3*x-2+3*y);
  d2psi[7][1][0] = d2psi[7][0][1];
  d2psi[7][1][1] = 7*x2*(x-1+3*y);
  d2psi[8][0][0] = -y*(10-37*y-84*x+27*y2+111*x*y+96*x2);
  d2psi[8][0][1] = -x*(10-74*y-42*x+81*y2+111*x*y+32*x2);
  d2psi[8][1][0] = d2psi[8][0][1];
  d2psi[8][1][1] = -x2*(-37+81*y+37*x);
  d2psi[9][0][0] = 3*x+.5*y2-12*x2-.5*y3-1.5*x*y2+10*x3;
  d2psi[9][0][1] = -.5*x*y*(3*x-2+3*y);
  d2psi[9][1][0] = d2psi[9][0][1];
  d2psi[9][1][1] = -.5*x2*(x-1+3*y);
  d2psi[10][0][0] = y*(2-7*y-18*x+5*y2+21*x*y+24*x2);
  d2psi[10][0][1] = x*(2-14*y-9*x+15*y2+21*x*y+8*x2);
  d2psi[10][1][0] = d2psi[10][0][1];
  d2psi[10][1][1] = x2*(-7+15*y+7*x);
  d2psi[11][0][0] = -.5*y2*(5*y-5+9*x);
  d2psi[11][0][1] = -.5*x*y*(15*y-10+9*x);
  d2psi[11][1][0] = d2psi[11][0][1];
  d2psi[11][1][1] = -.5*x2*(15*y-5+3*x);
  d2psi[12][0][0] = -30*y2*(3*x-1+y);
  d2psi[12][0][1] = -30*x*y*(3*x-2+3*y);
  d2psi[12][1][0] = d2psi[12][0][1];
  d2psi[12][1][1] = 60*y-180*y2+30*x2+120*y3-90*x2*y-30*x3;
  d2psi[13][0][0] = -y2*(37*y-37+81*x);
  d2psi[13][0][1] = -y*(10-42*y-74*x+32*y2+111*x*y+81*x2);
  d2psi[13][1][0] = d2psi[13][0][1];
  d2psi[13][1][1] = -x*(10-84*y-37*x+96*y2+111*x*y+27*x2);
  d2psi[14][0][0] = 7*y2*(3*x-1+y);
  d2psi[14][0][1] = 7*x*y*(3*x-2+3*y);
  d2psi[14][1][0] = d2psi[14][0][1];
  d2psi[14][1][1] = -24*y+84*y2-7*x2-60*y3+21*x2*y+7*x3;
  d2psi[15][0][0] = -.5*y2*(3*y-5+15*x);
  d2psi[15][0][1] = -.5*x*y*(9*y-10+15*x);
  d2psi[15][1][0] = d2psi[15][0][1];
  d2psi[15][1][1] = -.5*x2*(9*y-5+5*x);
  d2psi[16][0][0] = y2*(7*y-7+15*x);
  d2psi[16][0][1] = y*(2-9*y-14*x+8*y2+21*x*y+15*x2);
  d2psi[16][1][0] = d2psi[16][0][1];
  d2psi[16][1][1] = x*(2-18*y-7*x+24*y2+21*x*y+5*x2);
  d2psi[17][0][0] = -.5*y2*(3*x-1+y);
  d2psi[17][0][1] = -.5*x*y*(3*x-2+3*y);
  d2psi[17][1][0] = d2psi[17][0][1];
  d2psi[17][1][1] = 3*y-12*y2+.5*x2+10*y3-1.5*x2*y-.5*x3;
  d2psi[18][0][0] = -32*y*(1-6*x-2*y+6*x2+6*x*y+y2);
  d2psi[18][0][1] = -32*x*(1-3*x-4*y+2*x2+6*x*y+3*y2);
  d2psi[18][1][0] = d2psi[18][0][1];
  d2psi[18][1][1] = -32*x2*(-2+3*y+2*x);
  d2psi[19][0][0] = 16*sqrt2*(3*x-1+y)*y2;
  d2psi[19][0][1] = 16*sqrt2*(3*x-2+3*y)*x*y;
  d2psi[19][1][0] = d2psi[19][0][1];
  d2psi[19][1][1] = 16*sqrt2*x2*(x-1+3*y);
  d2psi[20][0][0] = -32*y2*(-2+2*y+3*x);
  d2psi[20][0][1] = -32*y*(1-4*x-3*y+3*x2+6*x*y+2*y2);
  d2psi[20][1][0] = d2psi[20][0][1];
  d2psi[20][1][1] = -32*x*(1-2*x-6*y+x2+6*x*y+6*y2);

  return;
}

void argyrisTriangleBdy(svector<1> xi, double psi[6], double psin[5]) 
{

  double x = xi[0], x2=x*x;
  double omx=1-x, omxsq=omx*omx;
  double opx=1+x, opxsq=opx*opx;
  
  psi[0] = 0.0625 *(3*x2+9*x+8)*omxsq*omx;
  psi[1] = 0.0625 *(3*x+5)*opx*omxsq*omx;
  psi[2] = 0.0625 *opxsq*omxsq*omx;
  psi[3] = 0.0625 *(3*x2-9*x+8)*opxsq*opx;
  psi[4] = 0.0625 *omx*(3*x-5)*opxsq*opx;
  psi[5] = 0.0625 *omxsq*opxsq*opx;
  
  psin[0] = -0.25*x*(2*x+3)*omxsq;
  psin[1] = -0.25*x*opx*omxsq;
  psin[2] = -0.25*x*(-3+2*x)*opxsq;
  psin[3] = -0.25*x*omx*opxsq;
  psin[4] = omxsq*opxsq;
  
  return; 
}

void argyrisBdyTransform(int *nodes, sbmatrix<2> &dxds,
		      double Phi[6],  double Phin[5], 
		      double phi[13], double phin[13])
{
  int i, j, k, n0 = 0, ii = 0, i0 = 0;
  double det;

  svector<2> normal = dxds.normal(det);

  for(i = 0; i < 2; i++)     // for each vertex
    {
    phi[ii] = Phi[i0];
    phin[ii] = 0.0;

    ii ++;
    i0 ++;
    
    for(j = 0; j < 2; j++)   // gradient dof
      {
      phi[ii] = dxds[j][0] * Phi[i0];
      
      phin[ii] = normal[j] * Phin[n0];
      
      ii ++;
      }
    
    i0 ++;
    n0 ++;
    
    for(j = 0; j < 2; j++)   // Hessian dof
    for(k = j; k < 2; k++)   
      {
      phi[ii] = dxds[j][0] * dxds[k][0] * Phi[i0];

      phin[ii] = 0.5 *
	(dxds[j][0] * normal[k] + dxds[k][0] * normal[j]) * Phin[n0];

      ii++;
      }
    
    i0 ++;
    n0 ++;
    }

  assert(ii == 12);
  assert(i0 == 6);
  assert(n0 == 4);

  phi[ii]  = 0.0;         // function values don't depend upon du/dn
  phin[ii] = Phin[n0];    // normal derivative at the middle of the edge

  if(nodes[0] > nodes[1]) phin[ii] *= -1;

  return;
}

void argyrisTransform(int *nodes, smatrix<2> &dxds,
		   double Phi[21], svector<2> Dphi[21], smatrix<2> D2phi[21],
		   double phi[21], svector<2> dphi[21], smatrix<2> d2phi[21])
{
  const int ndim=2;
  int i,j, k, alpha,beta, idxab, ii, i0=0;

  double temp, det, modFtn, psi[3];
  double a0, a1;
  svector<2> nn[3], dpsi[3];
  smatrix<2> d2psi[3];

  svector<2> dphiMid[18][2] =
    {{svector<2>(-15.0/8, 0), svector<2>(0, -15.0/8)},
     {svector<2>(-7.0/16, 0), svector<2>(0, 0)},
     {svector<2>(0, 0), svector<2>(0, -7.0/16)},
     {svector<2>(-1.0/32, 0), svector<2>(0, 0)},
     {svector<2>(0, 0), svector<2>(0, 0)},
     {svector<2>(0, 0), svector<2>(0, -1.0/32)},
     
     {svector<2>(15.0/16, -15.0/16), svector<2>(15.0/8, 0)},
     {svector<2>(-7.0/32, 7.0/32), svector<2>(-7.0/16, 0)},
     {svector<2>(7.0/32, -7.0/32), svector<2>(0, 0)},
     {svector<2>(1.0/64, -1.0/64), svector<2>(1.0/32, 0)},
     {svector<2>(-1.0/32, 1.0/32), svector<2>(0, 0)},
     {svector<2>(1.0/64, -1.0/64), svector<2>(0, 0)},
     
     {svector<2>(0, 15.0/8), svector<2>(-15.0/16, 15.0/16)},
     {svector<2>(0, 0), svector<2>(-7.0/32, 7.0/32)},
     {svector<2>(0, -7.0/16), svector<2>(7.0/32, -7.0/32)},
     {svector<2>(0, 0), svector<2>(-1.0/64, 1.0/64)},
     {svector<2>(0, 0), svector<2>(1.0/32, -1.0/32)},
     {svector<2>(0, 1.0/32), svector<2>(-1.0/64, 1.0/64)}};
  
  smatrix<2> dsdx  = dxds.inverse(det);
  smatrix<2> dsdxT = dsdx.transpose();
  
  nn[0] = dsdxT * svector<2>(0.0,-1.0);
  nn[1] = dsdxT * svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0));
  nn[2] = dsdxT * svector<2>(-1.0,0.0);
  
  for(i = 0; i < 3; i++)
    {
       ii = 18 + i;

       modFtn = 1.0/nn[i].norm();

       if(nodes[i] > nodes[(i+1)%3]) modFtn *= -1;

       phi[ii]   = modFtn * Phi[ii];
       dphi[ii]  = modFtn * dsdxT * Dphi[ii];
       d2phi[ii] = modFtn * dsdxT * D2phi[ii] * dsdx;

       nn[i] = (modFtn*modFtn) * nn[i];
     }
 
  // currently Hermite elements must be affine equivalent to the master

   ii = 0;

   for(i = 0; i < 3; i++)
     {
     j = (i+2)%3;    // edge before node i

     a0 = nn[i].dot(dsdxT * dphiMid[ii][0]);
     a1 = nn[j].dot(dsdxT * dphiMid[ii][1]);

     phi[ii]  =          Phi[i0]   - a0 * Phi[18+i]   - a1 * Phi[18+j];
     dphi[ii] = dsdxT * (Dphi[i0]  - a0 * Dphi[18+i]  - a1 * Dphi[18+j]);
     d2phi[ii] = dsdxT* (D2phi[i0] - a0 * D2phi[18+i] - a1 * D2phi[18+j])*dsdx;

     ii ++;
     i0 ++;

     for(alpha = 0; alpha < ndim; alpha++)    // gradient dof
       {
	 a0 = nn[i].dot(dsdxT * dphiMid[i0+alpha][0]);
	 a1 = nn[j].dot(dsdxT * dphiMid[i0+alpha][1]);

	 psi[alpha]   = Phi[i0+alpha]   - a0 * Phi[18+i]   - a1 * Phi[18+j];
	 dpsi[alpha]  = Dphi[i0+alpha]  - a0 * Dphi[18+i]  - a1 * Dphi[18+j];
	 d2psi[alpha] = D2phi[i0+alpha] - a0 * D2phi[18+i] - a1 * D2phi[18+j];
       }

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       for(alpha = 0; alpha < ndim; alpha++)
	 {
	   phi[ii]   += dxds[j][alpha] * psi[alpha];
	   dphi[ii]  += dxds[j][alpha] * dpsi[alpha];
	   d2phi[ii] += dxds[j][alpha] * d2psi[alpha];
	 }

       dphi[ii] = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;

       ii ++;
       }

     i0 += ndim;

     j = (i+2)%3;    // edge before node i

     for(idxab = 0; idxab < 3; idxab++)    // hessian dof
       {
	 a0 = nn[i].dot(dsdxT * dphiMid[i0+idxab][0]);
	 a1 = nn[j].dot(dsdxT * dphiMid[i0+idxab][1]);

	 psi[idxab]   = Phi[i0+idxab]   - a0 * Phi[18+i]   - a1 * Phi[18+j];
	 dpsi[idxab]  = Dphi[i0+idxab]  - a0 * Dphi[18+i]  - a1 * Dphi[18+j];
	 d2psi[idxab] = D2phi[i0+idxab] - a0 * D2phi[18+i] - a1 * D2phi[18+j];
       }

     for(j = 0; j < ndim; j++)   // Hessian dof
     for(k = j; k < ndim; k++)   
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       idxab = 0;   // indexing for symmetric matrx 0 <= a <= b < ndim

       for(alpha = 0;    alpha < ndim; alpha++)
       for(beta  = alpha; beta < ndim; beta++)
	 {
	 temp = 0.5 *
	   (dxds[j][alpha] * dxds[k][beta] + dxds[k][alpha] * dxds[j][beta]);

	 phi[ii]   += temp * psi[idxab];
	 dphi[ii]  += temp * dpsi[idxab];
	 d2phi[ii] += temp * d2psi[idxab];

	 idxab ++;
	 }

       dphi[ii] = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;
       ii++;
       }

     i0 += (ndim*(ndim+1))/2;
     }

   assert(ii == 18);
   assert(i0 == 18);

  return;
}

void argyrisBdyInterpolate(int *nodes, svector<2> &nn,
			double uu[3], svector<2> du[3], smatrix<2> d2u[3],
			double ue[13])
{
  const int ndim=2;
  int i,j,k, ii=0;
 
   for(i = 0; i < 2; i++)
     {
     ue[ii] = uu[i];

     ii ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       ue[ii] = du[i][j];

       ii ++;
       }

     for(j = 0; j < ndim; j++)   // Hessian dof
     for(k = j; k < ndim; k++)   
       {
       ue[ii] = d2u[i][j][k];

       if(j != k) ue[ii] *= 2.0;

       ii++;
       }
     }

   assert(ii == 12);

   ue[ii] = nn.dot(du[2]);           // du/dn at the mid point

   if(nodes[0] > nodes[1]) ue[ii] *= -1;

  return;
}

void argyrisInterpolate(int *nodes, smatrix<2> &dxds,
		     double uu[6], svector<2> du[6], smatrix<2> d2u[6],
		     double ue[21])
{
  const int ndim=2;
  int i,j,k, ii=0;
 
   for(i = 0; i < 3; i++)
     {
     ue[ii] = uu[i];

     ii ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       ue[ii] = du[i][j];

       ii ++;
       }

     for(j = 0; j < ndim; j++)   // Hessian dof
     for(k = j; k < ndim; k++)   
       {
       ue[ii] = d2u[i][j][k];

       if(j != k) ue[ii] *= 2.0;

       ii++;
       }
     }

   assert(ii == 18);

   // compute the normal derivative at the edge centers 
   // assuming the map from the parent is affine

   double det;

   smatrix<2> dsdx  = dxds.inverse(det);
   smatrix<2> dsdxT = dsdx.transpose();

   svector<2> nn[3];

   nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1.0/ nn[0].norm()) * nn[0];
   nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1.0/ nn[1].norm()) * nn[1];
   nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1.0/ nn[2].norm()) * nn[2];
  
   for(i = 0; i < 3; i++) 
     {
       ue[18+i] = du[3+i].dot(nn[i]);

       if(nodes[i] > nodes[(i+1)%3]) ue[18 + i] *= -1;
     }

  return;
}

// BFS Element Functions 
// BFS Element Functions 
// BFS Element Functions 

void bfsSquareBdy(svector<1> xi, double psi[4], double psin[4]) 
{

  double x = xi[0];
  double omx=1-x, omxsq=omx*omx;
  double opx=1+x, opxsq=opx*opx;
  
  psi[0] = psin[0] =  .25*(x+2)*omxsq;
  psi[1] = psin[1] =  .25*opx*omxsq;
  psi[2] = psin[2] = -.25*(x-2)*opxsq;
  psi[3] = psin[3] = -.25*omx*opxsq;

  return; 
}

void bfsSquare(svector<2> xi, double phi[16], 
	       svector<2> *dphi, smatrix<2> *d2phi)
{
  double  phi00e, phi01e, phi10e, phi11e;
  double  phi00a, phi01a, phi10a, phi11a;
  double d00e,d01e,d10e,d11e,d200e,d201e,d210e,d211e;
  double d00a,d01a,d10a,d11a,d200a,d201a,d210a,d211a;

  double eta = xi[0], ata = xi[1];

  double theta = 1.5*eta, thata=1.5*ata;
  double ope = 1+eta, opa = 1+ata, ome = 1-eta, oma = 1-ata;

  phi00e= 0.25*ome*ome*(2+eta); d00e=-0.75*ome*ope;       d200e= theta;
  phi01e= 0.25*ome*ome* ope   ; d01e=-0.25*ome*(3*eta+1); d201e= theta-0.5;
  phi10e= 0.25*ope*ope*(2-eta); d10e= 0.75*ome*ope;       d210e=-theta;
  phi11e=-0.25*ope*ope* ome   ; d11e= 0.25*ope*(3*eta-1); d211e= theta+0.5;

  phi00a= 0.25*oma*oma*(2+ata); d00a=-0.75*oma*opa;       d200a= thata;
  phi01a= 0.25*oma*oma* opa   ; d01a=-0.25*oma*(3*ata+1); d201a= thata-0.5;
  phi10a= 0.25*opa*opa*(2-ata); d10a= 0.75*oma*opa;       d210a=-thata;
  phi11a=-0.25*opa*opa* oma   ; d11a= 0.25*opa*(3*ata-1); d211a= thata+0.5;

  phi[0] = phi00e*phi00a; 
  phi[1] = phi01e*phi00a; 
  phi[2] = phi00e*phi01a; 
  phi[3] = phi01e*phi01a; 
  phi[4] = phi10e*phi00a; 
  phi[5] = phi11e*phi00a; 
  phi[6] = phi10e*phi01a; 
  phi[7] = phi11e*phi01a; 
  phi[8] = phi10e*phi10a; 
  phi[9] = phi11e*phi10a; 
  phi[10]= phi10e*phi11a; 
  phi[11]= phi11e*phi11a; 
  phi[12]= phi00e*phi10a; 
  phi[13]= phi01e*phi10a; 
  phi[14]= phi00e*phi11a; 
  phi[15]= phi01e*phi11a; 

  if(dphi == NULL) return;

  dphi[0][0] = d00e*phi00a; 
  dphi[0][1] = phi00e*d00a;
  dphi[1][0] = d01e*phi00a; 
  dphi[1][1] = phi01e*d00a;
  dphi[2][0] = d00e*phi01a; 
  dphi[2][1] = phi00e*d01a;
  dphi[3][0] = d01e*phi01a; 
  dphi[3][1] = phi01e*d01a;
  dphi[4][0] = d10e*phi00a; 
  dphi[4][1] = phi10e*d00a;
  dphi[5][0] = d11e*phi00a; 
  dphi[5][1] = phi11e*d00a;
  dphi[6][0] = d10e*phi01a; 
  dphi[6][1] = phi10e*d01a;
  dphi[7][0] = d11e*phi01a; 
  dphi[7][1] = phi11e*d01a;
  dphi[8][0] = d10e*phi10a; 
  dphi[8][1] = phi10e*d10a;
  dphi[9][0] = d11e*phi10a; 
  dphi[9] [1]= phi11e*d10a;
  dphi[10][0]= d10e*phi11a; 
  dphi[10][1]= phi10e*d11a;
  dphi[11][0]= d11e*phi11a; 
  dphi[11][1]= phi11e*d11a;
  dphi[12][0]= d00e*phi10a; 
  dphi[12][1]= phi00e*d10a;
  dphi[13][0]= d01e*phi10a; 
  dphi[13][1]= phi01e*d10a;
  dphi[14][0]= d00e*phi11a; 
  dphi[14][1]= phi00e*d11a;
  dphi[15][0]= d01e*phi11a; 
  dphi[15][1]= phi01e*d11a;

  if(d2phi == NULL) return;

  d2phi[0][0][0] = d200e*phi00a; 
  d2phi[0][0][1] = d00e*d00a; 
  d2phi[0][1][1] = phi00e*d200a;
  d2phi[1][0][0] = d201e*phi00a; 
  d2phi[1][0][1] = d01e*d00a; 
  d2phi[1][1][1] = phi01e*d200a;
  d2phi[2][0][0] = d200e*phi01a; 
  d2phi[2][0][1] = d00e*d01a; 
  d2phi[2][1][1] = phi00e*d201a;
  d2phi[3][0][0] = d201e*phi01a; 
  d2phi[3][0][1] = d01e*d01a; 
  d2phi[3][1][1] = phi01e*d201a;

  d2phi[4][0][0] = d210e*phi00a; 
  d2phi[4][0][1] = d10e*d00a; 
  d2phi[4][1][1] = phi10e*d200a;
  d2phi[5][0][0] = d211e*phi00a; 
  d2phi[5][0][1] = d11e*d00a; 
  d2phi[5][1][1] = phi11e*d200a;
  d2phi[6][0][0] = d210e*phi01a; 
  d2phi[6][0][1] = d10e*d01a; 
  d2phi[6][1][1] = phi10e*d201a;
  d2phi[7][0][0] = d211e*phi01a; 
  d2phi[7][0][1] = d11e*d01a; 
  d2phi[7][1][1] = phi11e*d201a;

  d2phi[8][0][0]  =  d210e*phi10a; 
  d2phi[8][0][1] = d10e*d10a; 
  d2phi[8][1][1] = phi10e*d210a;
  d2phi[9][0][0] = d211e*phi10a; 
  d2phi[9][0][1] = d11e*d10a; 
  d2phi[9][1][1] = phi11e*d210a;
  d2phi[10][0][0]= d210e*phi11a; 
  d2phi[10][0][1]= d10e*d11a; 
  d2phi[10][1][1]= phi10e*d211a;
  d2phi[11][0][0]= d211e*phi11a; 
  d2phi[11][0][1]= d11e*d11a; 
  d2phi[11][1][1]= phi11e*d211a;

  d2phi[12][0][0] = d200e*phi10a; 
  d2phi[12][0][1] = d00e*d10a; 
  d2phi[12][1][1] = phi00e*d210a;
  d2phi[13][0][0] = d201e*phi10a; 
  d2phi[13][0][1] = d01e*d10a; 
  d2phi[13][1][1] = phi01e*d210a;
  d2phi[14][0][0] = d200e*phi11a; 
  d2phi[14][0][1] = d00e*d11a; 
  d2phi[14][1][1] = phi00e*d211a;
  d2phi[15][0][0] = d201e*phi11a; 
  d2phi[15][0][1] = d01e*d11a; 
  d2phi[15][1][1] = phi01e*d211a;

  for(int i = 0; i < 16; i++) d2phi[i][1][0] = d2phi[i][0][1];

  return;
}

void bfsSquareBdyTransform(int *nodes, sbmatrix<2> &dxds,
		     double Phi[4],  double Phin[4], 
		     double phi[8], double phin[8])
{
  int i, j, n0 = 0, ii = 0, i0 = 0;
  double det;

  svector<2> normal = dxds.normal(det);

  for(i = 0; i < 2; i++)     // for each vertex
    {
    phi[ii] = Phi[i0];
    phin[ii] = 0.0;

    ii ++;
    i0 ++;
    
    for(j = 0; j < 2; j++)   // gradient dof
      {
      phi[ii] = dxds[j][0] * Phi[i0];
      
      phin[ii] = normal[j] * Phin[n0];
      
      ii ++;
      }
    
    i0 ++;
    n0 ++;

    // Hessian dof (mixed derivative)
    
    phin[ii] = 0.0;
    phin[ii] = 0.5 *
	(dxds[0][0] * normal[1] + dxds[1][0] * normal[0]) * Phin[n0];

    ii++;
    
    n0 ++;
    }

  assert(ii == 8);
  assert(i0 == 4);
  assert(n0 == 4);
  
  return;
}

void bfsSquareTransform(int *nodes, smatrix<2> &dxds,
		   double Phi[16], svector<2> Dphi[16], smatrix<2> D2phi[16],
		   double phi[16], svector<2> dphi[16], smatrix<2> d2phi[16])
{
  const int ndim=2;
  int i,j, alpha, ii=0, i0=0;

  double temp, det;
  smatrix<2> dsdx  = dxds.inverse(det);
  smatrix<2> dsdxT = dsdx.transpose();
 
  // currently Hermite elements must be affine equivalent to the master

   for(i = 0; i < 4; i++)
     {
     phi[ii] = Phi[i0];
     dphi[ii] = dsdxT * Dphi[i0];
     d2phi[ii] = dsdxT * D2phi[i0] * dsdx;

     ii ++;
     i0 ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       for(alpha = 0; alpha < ndim; alpha++)
	 {
	   phi[ii]   += dxds[j][alpha] * Phi[i0+alpha];
	   dphi[ii]  += dxds[j][alpha] * Dphi[i0+alpha];
	   d2phi[ii] += dxds[j][alpha] * D2phi[i0+alpha];
	 }

       dphi[ii] = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;

       ii ++;
       }

     i0 += ndim;

     // Hessian dof (mixed derivative; need dxds[][] like a permutation)

     temp = 0.5 *
       (dxds[0][0] * dxds[1][1] + dxds[1][0] * dxds[0][1]);
     
     phi[ii]   = temp * Phi[i0];
     dphi[ii]  = temp * Dphi[i0];
     d2phi[ii] = temp * D2phi[i0];
     
     dphi[ii] = dsdxT * dphi[ii];
     d2phi[ii] = dsdxT * d2phi[ii] * dsdx;
     ii++;

     i0 += (ndim*(ndim-1))/2;
     }

   assert(ii == 16);
   assert(i0 == 16);

  return;
}

void bfsSquareBdyInterpolate(int *nodes,
			double uu[2], svector<2> du[2], smatrix<2> d2u[2],
			double ue[8])
{
  const int ndim=2;
  int i,j, ii=0;
  
  for(i = 0; i < 2; i++)
    {
    ue[ii++] = uu[i];

    // gradient dof
    
    for(j = 0; j < ndim; j++) ue[ii++] = du[i][j];
    
    // Hessian dof (mixed derivative)
    
    ue[ii++] = 2 * d2u[i][0][1];
    }
  
  assert(ii == 8);
  
  return;
}

void bfsSquareInterpolate(int *nodes,
		     double uu[4], svector<2> du[4], smatrix<2> d2u[4],
		     double ue[16])
{
  const int ndim=2;
  int i,j, ii=0;
  
  for(i = 0; i < 4; i++)   // each vertex
    {
    ue[ii++] = uu[i];

    // gradient dof
    
    for(j = 0; j < ndim; j++) ue[ii++] = du[i][j];
    
    // Hessian dof (mixed derivative)
    
    ue[ii++] = 2 * d2u[i][0][1];
    }
  
  assert(ii == 16);
  
  return;
}

// BFS cube -- tensor product of 1d cubic Hermite
// BFS cube -- tensor product of 1d cubic Hermite
// BFS cube -- tensor product of 1d cubic Hermite

void hermiteLine4(double x, double psi[4], double dpsi[4], double d2psi[4])
{
  double omx=1-x, omxsq=omx*omx;
  double opx=1+x, opxsq=opx*opx;

  psi[0] =  .25*(x+2)*omxsq;
  psi[1] =  .25*opx*omxsq;
  psi[2] = -.25*(x-2)*opxsq;
  psi[3] = -.25*omx*opxsq;
  
  if(dpsi == NULL) return;

  dpsi[0] = -.75*omx*opx;
  dpsi[1] = -.25*(3*x+1)*omx;
  dpsi[2] =  .75*omx*opx;
  dpsi[3] =  .25*opx*(3*x-1);
  
  if(d2psi == NULL) return;

  d2psi[0] =  1.5*x;
  d2psi[1] =  1.5*x-0.5;
  d2psi[2] = -1.5*x;
  d2psi[3] =  1.5*x+0.5;

  return;
}

void bfsCubeBdy(svector<2> xi, double psi[16], double psin[16]) 
{
  int i,j, i0, ii = 0;

  double phix[4], phiy[4];

  hermiteLine4(xi[0], phix, NULL, NULL);
  hermiteLine4(xi[1], phiy, NULL, NULL);

  for(j  = 0; j  < 3; j+=2)
  for(i0 = 0; i0 < 3; i0+=2)
    {
    if(j == 0) i = i0;
    else       i = 2 - i0;
 
    psin[ii] = psi[ii] = phix[i]   * phiy[j];   ii++;
    psin[ii] = psi[ii] = phix[i+1] * phiy[j];   ii++;
    psin[ii] = psi[ii] = phix[i]   * phiy[j+1]; ii++;
    psin[ii] = psi[ii] = phix[i+1] * phiy[j+1]; ii++;
    }

  assert(ii == 16);

  return;
}

void bfsCube(svector<3> xi, double psi[64], 
	     svector<3> dpsi[64], smatrix<3> d2psi[64])
{
  int i,j,k, ii = 0;

  double phix[4], dphix[4], d2phix[4];
  double phiy[4], dphiy[4], d2phiy[4];
  double phiz[4], dphiz[4], d2phiz[4];

  hermiteLine4(xi[0], phix, dphix, d2phix);
  hermiteLine4(xi[1], phiy, dphiy, d2phiy);
  hermiteLine4(xi[2], phiz, dphiz, d2phiz);

  for(k = 0; k < 3; k+=2)
  for(j = 0; j < 3; j+=2)
  for(i = 0; i < 3; i+=2)
    {
    psi[ii] = phix[i] * phiy[j] * phiz[k];        // function dof

    dpsi[ii][0] = dphix[i] *  phiy[j] *  phiz[k];       
    dpsi[ii][1] =  phix[i] * dphiy[j] *  phiz[k];       
    dpsi[ii][2] =  phix[i] *  phiy[j] * dphiz[k];       

    d2psi[ii][0][0] = d2phix[i] *  phiy[j] *  phiz[k];       
    d2psi[ii][0][1] =  dphix[i] * dphiy[j] *  phiz[k];       
    d2psi[ii][0][2] =  dphix[i] *  phiy[j] * dphiz[k];       

    d2psi[ii][1][1] =  phix[i] * d2phiy[j] *  phiz[k];       
    d2psi[ii][1][2] =  phix[i] *  dphiy[j] * dphiz[k];       

    d2psi[ii][2][2] =  phix[i] *  phiy[j] * d2phiz[k];       

    ii++;

    psi[ii] = phix[i+1] * phiy[j] * phiz[k];      // derivative dof

    dpsi[ii][0] = dphix[i+1] *  phiy[j] *  phiz[k];       
    dpsi[ii][1] =  phix[i+1] * dphiy[j] *  phiz[k];       
    dpsi[ii][2] =  phix[i+1] *  phiy[j] * dphiz[k];       

    d2psi[ii][0][0] = d2phix[i+1] *  phiy[j] *  phiz[k];       
    d2psi[ii][0][1] =  dphix[i+1] * dphiy[j] *  phiz[k];       
    d2psi[ii][0][2] =  dphix[i+1] *  phiy[j] * dphiz[k];       

    d2psi[ii][1][1] =  phix[i+1] * d2phiy[j] *  phiz[k];       
    d2psi[ii][1][2] =  phix[i+1] *  dphiy[j] * dphiz[k];       

    d2psi[ii][2][2] =  phix[i+1] *  phiy[j] * d2phiz[k];       

    ii++;

    psi[ii] = phix[i] * phiy[j+1] * phiz[k];

    dpsi[ii][0] = dphix[i] *  phiy[j+1] *  phiz[k];       
    dpsi[ii][1] =  phix[i] * dphiy[j+1] *  phiz[k];       
    dpsi[ii][2] =  phix[i] *  phiy[j+1] * dphiz[k];       

    d2psi[ii][0][0] = d2phix[i] *  phiy[j+1] *  phiz[k];       
    d2psi[ii][0][1] =  dphix[i] * dphiy[j+1] *  phiz[k];       
    d2psi[ii][0][2] =  dphix[i] *  phiy[j+1] * dphiz[k];       

    d2psi[ii][1][1] =  phix[i] * d2phiy[j+1] *  phiz[k];       
    d2psi[ii][1][2] =  phix[i] *  dphiy[j+1] * dphiz[k];       

    d2psi[ii][2][2] =  phix[i] *  phiy[j+1] * d2phiz[k];       

    ii++;

    psi[ii] = phix[i] * phiy[j] * phiz[k+1];

    dpsi[ii][0] = dphix[i] *  phiy[j] *  phiz[k+1];       
    dpsi[ii][1] =  phix[i] * dphiy[j] *  phiz[k+1];       
    dpsi[ii][2] =  phix[i] *  phiy[j] * dphiz[k+1];       

    d2psi[ii][0][0] = d2phix[i] *  phiy[j] *  phiz[k+1];       
    d2psi[ii][0][1] =  dphix[i] * dphiy[j] *  phiz[k+1];       
    d2psi[ii][0][2] =  dphix[i] *  phiy[j] * dphiz[k+1];       

    d2psi[ii][1][1] =  phix[i] * d2phiy[j] *  phiz[k+1];       
    d2psi[ii][1][2] =  phix[i] *  dphiy[j] * dphiz[k+1];       

    d2psi[ii][2][2] =  phix[i] *  phiy[j] * d2phiz[k+1];       

    ii++;

    psi[ii] = phix[i+1] * phiy[j+1] * phiz[k];      // mixed derivative dof

    dpsi[ii][0] = dphix[i+1] *  phiy[j+1] *  phiz[k];       
    dpsi[ii][1] =  phix[i+1] * dphiy[j+1] *  phiz[k];       
    dpsi[ii][2] =  phix[i+1] *  phiy[j+1] * dphiz[k];       

    d2psi[ii][0][0] = d2phix[i+1] *  phiy[j+1] *  phiz[k];       
    d2psi[ii][0][1] =  dphix[i+1] * dphiy[j+1] *  phiz[k];       
    d2psi[ii][0][2] =  dphix[i+1] *  phiy[j+1] * dphiz[k];       

    d2psi[ii][1][1] =  phix[i+1] * d2phiy[j+1] *  phiz[k];       
    d2psi[ii][1][2] =  phix[i+1] *  dphiy[j+1] * dphiz[k];       

    d2psi[ii][2][2] =  phix[i+1] *  phiy[j+1] * d2phiz[k];       

    ii++;

    psi[ii] = phix[i+1] * phiy[j] * phiz[k+1];

    dpsi[ii][0] = dphix[i+1] *  phiy[j] *  phiz[k+1];       
    dpsi[ii][1] =  phix[i+1] * dphiy[j] *  phiz[k+1];       
    dpsi[ii][2] =  phix[i+1] *  phiy[j] * dphiz[k+1];       

    d2psi[ii][0][0] = d2phix[i+1] *  phiy[j] *  phiz[k+1];       
    d2psi[ii][0][1] =  dphix[i+1] * dphiy[j] *  phiz[k+1];       
    d2psi[ii][0][2] =  dphix[i+1] *  phiy[j] * dphiz[k+1];       

    d2psi[ii][1][1] =  phix[i+1] * d2phiy[j] *  phiz[k+1];       
    d2psi[ii][1][2] =  phix[i+1] *  dphiy[j] * dphiz[k+1];       

    d2psi[ii][2][2] =  phix[i+1] *  phiy[j] * d2phiz[k+1];       

    ii++;

    psi[ii] = phix[i] * phiy[j+1] * phiz[k+1];

    dpsi[ii][0] = dphix[i] *  phiy[j+1] *  phiz[k+1];       
    dpsi[ii][1] =  phix[i] * dphiy[j+1] *  phiz[k+1];       
    dpsi[ii][2] =  phix[i] *  phiy[j+1] * dphiz[k+1];       

    d2psi[ii][0][0] = d2phix[i] *  phiy[j+1] *  phiz[k+1];       
    d2psi[ii][0][1] =  dphix[i] * dphiy[j+1] *  phiz[k+1];       
    d2psi[ii][0][2] =  dphix[i] *  phiy[j+1] * dphiz[k+1];       

    d2psi[ii][1][1] =  phix[i] * d2phiy[j+1] *  phiz[k+1];       
    d2psi[ii][1][2] =  phix[i] *  dphiy[j+1] * dphiz[k+1];       

    d2psi[ii][2][2] =  phix[i] *  phiy[j+1] * d2phiz[k+1];       

    ii++;

    psi[ii] = phix[i+1] * phiy[j+1] * phiz[k+1];

    dpsi[ii][0] = dphix[i+1] *  phiy[j+1] *  phiz[k+1];       
    dpsi[ii][1] =  phix[i+1] * dphiy[j+1] *  phiz[k+1];       
    dpsi[ii][2] =  phix[i+1] *  phiy[j+1] * dphiz[k+1];       

    d2psi[ii][0][0] = d2phix[i+1] *  phiy[j+1] *  phiz[k+1];       
    d2psi[ii][0][1] =  dphix[i+1] * dphiy[j+1] *  phiz[k+1];       
    d2psi[ii][0][2] =  dphix[i+1] *  phiy[j+1] * dphiz[k+1];       

    d2psi[ii][1][1] =  phix[i+1] * d2phiy[j+1] *  phiz[k+1];       
    d2psi[ii][1][2] =  phix[i+1] *  dphiy[j+1] * dphiz[k+1];       

    d2psi[ii][2][2] =  phix[i+1] *  phiy[j+1] * d2phiz[k+1];       

    ii++;
    }

  assert(ii == 64);

  for(i = 0; i < ii; i++)
    {
      for(j = 0;   j < 3; j++) 
      for(k = j+1; k < 3; k++) d2psi[i][k][j] = d2psi[i][j][k];
    }

  return;
}

void bfsCubeBdyTransform(int *nodes, sbmatrix<3> &dxds,
		      double Phi[16], double Phin[16], 
		      double phi[32], double phin[32])
{
  const int ndim=3;
  int i, j, k, n0 = 0, alpha, ii = 0, i0 = 0;
  double det;

  svector<ndim> normal = dxds.normal(det);

  for(i = 0; i < 4; i++)     // for each vertex
    {
    phi[ii] = Phi[i0];
    phin[ii] = 0.0;

    ii ++;
    i0 ++;                      // function value on face dof used
    
    for(j = 0; j < ndim; j++)   // gradient dof
      {
      phin[ii] = normal[j] * Phin[n0];

      phi[ii]  = 0.0;

      for(alpha = 0; alpha < ndim-1; alpha++)  
	phi[ii] += dxds[j][alpha] * Phi[i0+alpha];
      
      ii ++;
      }
    
    i0 += ndim-1;                 // (ndim-1) face gradients dof used
    n0 ++;                        // function face value     dof used
    
    for(j = 0;   j < ndim; j++)   // Hessian dof
    for(k = j+1; k < ndim; k++)   
      {
      phi[ii]  = 0.0;
      phin[ii] = 0.0;

      for(alpha = 0; alpha < ndim-1; alpha++)  
	{
	phin[ii] += 0.5 * ( dxds[j][alpha] * normal[k] + 
			    dxds[k][alpha] * normal[j] ) * Phin[n0+alpha];

	for(int beta = alpha+1; beta < ndim-1; beta++)  // (alpha,beta) = (0,1)
	  {
	  phi[ii]  += 0.5 * ( dxds[j][alpha] * dxds[k][beta] + 
			      dxds[k][alpha] * dxds[j][beta] ) * Phi[i0];
	  }
	}

      ii++;
      }
    
    i0 ++;          // mixed face derivative  dof used
    n0 += ndim-1;   // (ndim-1) face gradient dof used

    phi[ii] = 0.0;
    phin[ii] = (1.0/6.0) * Phin[n0] *
      ( dxds[0][0]*dxds[1][1]*normal[2] + dxds[0][0]*dxds[2][1]*normal[1] +
	dxds[1][0]*dxds[2][1]*normal[0] + dxds[1][0]*dxds[0][1]*normal[2] +
	dxds[2][0]*dxds[0][1]*normal[1] + dxds[2][0]*dxds[1][1]*normal[0] );

    ii ++;
    n0 ++;          // mixed face derivative dof used
    }

  assert(ii == 32);
  assert(i0 == 16);
  assert(n0 == 16);
  
  return;
}

void bfsCubeTransform(int *nodes, smatrix<3> &dxds,
		   double Phi[64], svector<3> Dphi[64], smatrix<3> D2phi[64],
		   double phi[64], svector<3> dphi[64], smatrix<3> d2phi[64])
{
  const int ndim=3;
  int i,j, k, alpha,beta, idxab, ii=0, i0=0;

  double temp, det;
  smatrix<3> dsdx  = dxds.inverse(det);
  smatrix<3> dsdxT = dsdx.transpose();
 
  // BFS elements require the Jacobian to be like a permutation matrix

  for(i = 0; i < 8; i++)    // each vertex
     {
     phi[ii] = Phi[i0];
     dphi[ii] = dsdxT * Dphi[i0];
     d2phi[ii] = dsdxT * D2phi[i0] * dsdx;

     ii ++;
     i0 ++;

     for(j = 0; j < ndim; j++)   // gradient dof
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       for(alpha = 0; alpha < ndim; alpha++)
	 {
	   phi[ii]   += dxds[j][alpha] * Phi[i0+alpha];
	   dphi[ii]  += dxds[j][alpha] * Dphi[i0+alpha];
	   d2phi[ii] += dxds[j][alpha] * D2phi[i0+alpha];
	 }

       dphi[ii] = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;

       ii ++;
       }

     i0 += ndim;

     for(j = 0;   j < ndim; j++)   // mixed derivative dof
     for(k = j+1; k < ndim; k++)   
       {
       phi[ii] = 0.0;
       dphi[ii].zero();
       d2phi[ii].zero();

       idxab = i0;   // indexing for symmetric matrx 0 <= a < b < ndim

       for(alpha = 0;      alpha < ndim; alpha++)
       for(beta  = alpha+1; beta < ndim; beta++)
	 {
	 temp = 0.5 *
	   (dxds[j][alpha] * dxds[k][beta] + dxds[k][alpha] * dxds[j][beta]);

	 phi[ii]   += temp * Phi[idxab];
	 dphi[ii]  += temp * Dphi[idxab];
	 d2phi[ii] += temp * D2phi[idxab];

	 idxab ++;
	 }

       dphi[ii] = dsdxT * dphi[ii];
       d2phi[ii] = dsdxT * d2phi[ii] * dsdx;
       ii++;
       }

     i0 += (ndim*(ndim-1))/2;

     temp = (1.0/6.0) *
      ( dxds[0][0]*dxds[1][1]*dxds[2][2] + dxds[0][0]*dxds[1][2]*dxds[2][1] + 
	dxds[0][1]*dxds[1][2]*dxds[2][0] + dxds[0][1]*dxds[1][0]*dxds[2][2] + 
	dxds[0][2]*dxds[1][0]*dxds[2][1] + dxds[0][2]*dxds[1][1]*dxds[2][0] ); 

     phi[ii]   = temp * Phi[i0];      // mixed tripple derivative
     dphi[ii]  = temp * Dphi[i0];
     d2phi[ii] = temp * D2phi[i0];

     dphi[ii] = dsdxT * dphi[ii];
     d2phi[ii] = dsdxT * d2phi[ii] * dsdx;

     ii ++;
     i0 ++;
     }

   assert(ii == 64);
   assert(i0 == 64);

  return;
}

void bfsCubeBdyInterpolate(int *nodes,
      double uu[4], svector<3> du[4], smatrix<3> d2u[4], smatrix<3> d3u[4][3],
      double ue[32])
{
  const int ndim=3;
  int i,j,k, ii=0;
  
  for(i = 0; i < 4; i++)   // each vertex
    {
    ue[ii++] = uu[i];

    // gradient dof
    
    for(j = 0; j < ndim; j++) ue[ii++] = du[i][j];
    
    // Hessian dofs (mixed derivative)
    
    for(j = 0; j < ndim; j++)
    for(k=j+1; k < ndim; k++) ue[ii++] = 2 * d2u[i][j][k];

    ue[ii++] = 6 * d3u[i][0][1][2];       // third mixed derivative
    }
  
  assert(ii == 32);
  
  return;
}

void bfsCubeInterpolate(int *nodes, double uu[8], 
	    	svector<3> du[8], smatrix<3> d2u[8], smatrix<3> d3u[8][3],
	       	double ue[64])
{
  const int ndim=3;
  int i,j,k, ii=0;
  
  for(i = 0; i < 8; i++)   // each vertex
    {
    ue[ii++] = uu[i];

    // gradient dof
    
    for(j = 0; j < ndim; j++) ue[ii++] = du[i][j];
    
    // Hessian dof (mixed second derivatives only)
    
    for(j = 0; j < ndim; j++)
    for(k=j+1; k < ndim; k++) ue[ii++] = 2 * d2u[i][j][k];

    ue[ii++] = 6 * d3u[i][0][1][2];       // third mixed derivative
    }
  
  assert(ii == 64);
  
  return;
}
