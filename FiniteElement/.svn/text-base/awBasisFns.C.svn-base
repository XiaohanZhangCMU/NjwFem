// S(div) baisis function
// functions are matrix valued, the only derivative is the (vector) div(S)
// their boundary values are vector valued basis functions for Sn;
// however, they need to be reconstructed from their dof's.

// Boundary basis functions are scalar, and transform() untangles 
// vertex and edge the dof's

template<const int Ndim>
struct SdivBdyBasisFunction : public basisFunction<Ndim>
{
  SdivBdyBasisFunction(basisFunction<Ndim> &bf, 
	       void (*bb)(svector<Ndim>, double *phi),
	       void (*ss)(const int *nodes, int &sign))
    : basisFunction<Ndim>(bf), basis(bb), dofSign(ss)  
  {return;}

  void (*basis)(svector<Ndim>,  double *phi);
  void (*dofSign)(const int *nodes, int &sign);
};

template<const int Ndim>
struct SdivBasisFunction : public basisFunction<Ndim>
{
  SdivBasisFunction(basisFunction<Ndim> &bf, 
	   SdivBdyBasisFunction<Ndim-1> *bdy,
	   void (*bb)(svector<Ndim>, smatrix<Ndim> *phi, svector<Ndim> *dphi),
	   void (*ss)(const int *nodes, int *signs),
	   void (*bdytt)(int *nodes, svector<Ndim> &nn,
			 double *Phi, svector<Ndim> *phi),
           void (*bdyii)(int *nodes, svector<Ndim> &nn,
		         smatrix<Ndim> *ss, double *dof),
	   void (*tt)(int *nodes, smatrix<Ndim> &dxds,
		      smatrix<Ndim> *Phi, svector<Ndim> *Dphi, 
		      smatrix<Ndim> *phi, svector<Ndim> *dphi),
           void (*interp)(int *nodes,  smatrix<Ndim> &dxds,
		          smatrix<Ndim> *ss, double *dof))
    : basisFunction<Ndim>(bf), boundary(bdy), basis(bb), dofSigns(ss),
      bdyTransform(bdytt), bdyInterpolate(bdyii),
      transform(tt), interpolate(interp)
  {return;}

  SdivBdyBasisFunction<Ndim-1> *boundary;
  void (*basis)(svector<Ndim>, smatrix<Ndim> *phi, svector<Ndim> *dphi);
  void (*dofSigns)(const int *nodes, int *signs);
  void (*bdyTransform)(int *nodes, svector<Ndim> &nn,
		       double *Phi, svector<Ndim> *phi);
  void (*bdyInterpolate)(int *nodes, svector<Ndim> &nn,
		         smatrix<Ndim> *ss, double *dof);
  void (*transform)(int *nodes, smatrix<Ndim> &dxds,
		    smatrix<Ndim> *Phi, svector<Ndim> *Dphi, 
		    smatrix<Ndim> *phi, svector<Ndim> *dphi);
  void (*interpolate)(int *nodes,  smatrix<Ndim> &dxds,
		    smatrix<Ndim> *ss, double *dof);
};

void awTriangle1Bdy(svector<1>, double phi[4]);
void awTriangleBdySign(const int *nodes, int &sign);
void awTriangle1BdyPermutation(const int nodes[2], int dof[10]);

void awTriangle1(svector<2> xx, smatrix<2> phi[24], svector<2> dphi[24]);
void awTriangle1Signs(const int *nodes, int signs[10]);
void awTriangle1Permutation(const int *nodes, int *dof);
void awTriangle1BdyTransform(int nodes[2], svector<2> &nn,
			     double Phi[4], svector<2> phi[10]);
void awTriangle1BdyInterpolate(int *nodes, svector<2> &nn,
		         smatrix<2> *ss, double *dof);
void awTriangle1Transform(int *nodes, smatrix<2> &dxds,
			  smatrix<2> *Phi, svector<2> *Dphi, 
			  smatrix<2> *phi, svector<2> *dphi);
void awTriangle1Interpolate(int nodes[3],  smatrix<2> &dxds,
			    smatrix<2> ss[10], double dof[24]);

// Construct the AW triangle with div(S) degree one
// Construct the AW triangle with div(S) degree one
// Construct the AW triangle with div(S) degree one

svector<1> AWtriangle1BdyPts[] = 
  {svector<1>(-1.0), 
   svector<1>( 1.0),
   svector<1>(-1.0/sqrt(3.0)), 
   svector<1>( 1.0/sqrt(3.0))};

basisFunction<1> AWtriangle1BdyBF =    // AW_1 boundary
  {10,
   {3,4},
   EmptyVertex,
   &Interval,
   4,
   AWtriangle1BdyPts,
   awTriangle1BdyPermutation
  };

SdivBdyBasisFunction<1> AWtriangle1BdyClass(
   AWtriangle1BdyBF, &awTriangle1Bdy, &awTriangleBdySign),
  *AWtriangle1Bdy = &AWtriangle1BdyClass;

svector<2> AWtriangle1Pts[] =  // 3 vertices, 6 edge pts, centroid 
  {svector<2>(0.0, 0.0), 
   svector<2>(1.0, 0.0), 
   svector<2>(0.0, 1.0), 
   svector<2>(0.5*(1.0 - 1.0/sqrt(3.0)), 0), 
   svector<2>(0.5*(1.0 + 1.0/sqrt(3.0)), 0),
   svector<2>(0.5*(1.0 + 1.0/sqrt(3.0)), 0.5*(1.0 - 1.0/sqrt(3.0))),
   svector<2>(0.5*(1.0 - 1.0/sqrt(3.0)), 0.5*(1.0 + 1.0/sqrt(3.0))),
   svector<2>(0, 0.5*(1.0 + 1.0/sqrt(3.0))),
   svector<2>(0, 0.5*(1.0 - 1.0/sqrt(3.0))), 
   svector<2>(1/3.0, 1/3.0)};
  
basisFunction<2> AWtriangle1BF =    // AW_1 element
  {24,
   {3,4,3},
   &AWtriangle1BdyBF,
   &Triangle,
   10,
   AWtriangle1Pts,
   awTriangle1Permutation
  };

SdivBasisFunction<2> AWtriangle1class(AWtriangle1BF, 
  AWtriangle1Bdy, 
  &awTriangle1, 
  &awTriangle1Signs,
  &awTriangle1BdyTransform, 
  &awTriangle1BdyInterpolate,
  &awTriangle1Transform,
  &awTriangle1Interpolate),
  *AWtriangle1 = &AWtriangle1class;

// ********************* Basis Functions ************************
// ********************* Basis Functions ************************
// ********************* Basis Functions ************************

void awTriangle1Bdy(svector<1> pxi, double psi[4])
{
  double xi = pxi[0], xisq = xi*xi;

  double gp = 1.0/sqrt(3.0); // (-gp, gp) = Gauss points on [-1,1]

  psi[0] = (3.0/4.0) * (xisq - 1.0/3.0) * (1.0 - xi);
  psi[1] = (3.0/4.0) * (xisq - 1.0/3.0) * (1.0 + xi);
  psi[2] = (3.0/4.0) * (1.0 - xisq) * (1.0 - xi/gp);
  psi[3] = (3.0/4.0) * (1.0 - xisq) * (1.0 + xi/gp);

  return;
}

void awTriangleBdySign(const int nodes[2], int &sign)
{
  // Convention: global edge orientation is low to high

  if(nodes[0] > nodes[1]) sign = -1; else sign = 1;
}

void awTriangle1BdyPermutation(const int nodes[2], int dof[10])
{
 if(nodes[0] > nodes[1]) 
   {
     int temp;
     temp = dof[6]; dof[6] = dof[8]; dof[8] = temp;
     temp = dof[7]; dof[7] = dof[9]; dof[9] = temp;
   }
}

void awTriangle1BdyTransform(int nodes[2], svector<2> &nn,
		    double Phi[4], svector<2> phi[10])
{
  int sign;

  awTriangleBdySign(nodes, sign);

  // vertices

  phi[0] = svector<2>(Phi[0]*nn[0], 0.0);
  phi[1] = 0.5 * svector<2>(Phi[0]*nn[1], Phi[0]*nn[0]);
  phi[2] = svector<2>(0.0, Phi[0]*nn[1]);

  phi[3] = svector<2>(Phi[1]*nn[0], 0.0);
  phi[4] = 0.5 * svector<2>(Phi[1]*nn[1], Phi[1]*nn[0]);
  phi[5] = svector<2>(0.0, Phi[1]*nn[1]);

  // mid points

  phi[6] = svector<2>(sign*Phi[2], 0.0);
  phi[7] = svector<2>(0.0, sign*Phi[2]);

  phi[8] = svector<2>(sign*Phi[3], 0.0);
  phi[9] = svector<2>(0.0, sign*Phi[3]);

  return;
}

void awTriangle1BdyInterpolate(int nodes[3],  svector<2> &normal,
			    smatrix<2> ss[4], double dof[20])
{
  // ss[] is the values of S  at the pts[]
  // compute the corresponding dof[] for the interpolant

  int sign;

  awTriangleBdySign(nodes, sign);

  // vertices

  dof[0] = ss[0][0][0];
  dof[1] = ss[0][0][1] * 2;
  dof[2] = ss[0][1][1];    

  dof[3] = ss[1][0][0];
  dof[4] = ss[1][0][1] * 2;
  dof[5] = ss[1][1][1];    

  // mid points

  svector<2> sn = ss[2] * (sign * normal);

  dof[6] = sn[0];
  dof[7] = sn[1];

  sn = ss[3] * (sign * normal);

  dof[8] = sn[0];
  dof[9] = sn[1];
  
  return;
}

void awTriangle1(svector<2> pxi, smatrix<2> phi[24], svector<2> *dphi = NULL)
{
  double xi = pxi[0], eta = pxi[1], omega = 1.0/sqrt(3.0);
  double xi2 = xi*xi, xi3 = xi2*xi, eta2 = eta*eta, eta3 = eta2*eta;
  double sqrt2 = sqrt(2.0);

  phi[0][0][0] = 1+.6*xi-7*eta-2.8*xi2+.8*xi*eta+12*eta2+1.2*xi3+5.4*xi2*eta-6*eta3;
  phi[0][0][1] = 3*xi*eta-3.6*xi2*eta-5.4*xi*eta2;
  phi[0][1][0] = 3*xi*eta-3.6*xi2*eta-5.4*xi*eta2;
  phi[0][1][1] = -.6*eta+1.2*xi*eta-1.2*eta2+3.6*xi*eta2+1.8*eta3;
  phi[1][0][0] = .5*(9*eta-7+9*xi)*xi*(-1+2*eta+xi);
  phi[1][0][1] = .5-3.5*xi-3.5*eta+6*xi2+15.5*xi*eta+6*eta2-3*xi3-13.5*xi2*eta-13.5*xi*eta2-3*eta3;
  phi[1][1][0] = .5-3.5*xi-3.5*eta+6*xi2+15.5*xi*eta+6*eta2-3*xi3-13.5*xi2*eta-13.5*xi*eta2-3*eta3;
  phi[1][1][1] = .5*(7+9*eta2-23*xi+27*xi*eta-16*eta+18*xi2)*eta;
  phi[2][0][0] = -.6*xi-1.2*xi2+1.2*xi*eta+1.8*xi3+3.6*xi2*eta;
  phi[2][0][1] = 3*xi*eta-5.4*xi2*eta-3.6*xi*eta2;
  phi[2][1][0] = 3*xi*eta-5.4*xi2*eta-3.6*xi*eta2;
  phi[2][1][1] = 1+1.2*eta3-7*xi+.6*eta-2.8*eta2+12*xi2+.8*xi*eta+5.4*xi*eta2-6*xi3;
  phi[3][0][0] = 1.4*xi+.8*xi2-2.8*xi*eta-1.2*xi3-5.4*xi2*eta;
  phi[3][0][1] = -3*xi*eta+3.6*xi2*eta+5.4*xi*eta2;
  phi[3][1][0] = -3*xi*eta+3.6*xi2*eta+5.4*xi*eta2;
  phi[3][1][1] = .6*eta-1.2*xi*eta+1.2*eta2-3.6*xi*eta2-1.8*eta3;
  phi[4][0][0] = 1.4*xi+1.3*xi2-2.8*xi*eta-2.7*xi3-5.4*xi2*eta;
  phi[4][0][1] = .5*xi-3*xi2-4*xi*eta+3*xi3+8.1*xi2*eta+5.4*xi*eta2;
  phi[4][1][0] = .5*xi-3*xi2-4*xi*eta+3*xi3+8.1*xi2*eta+5.4*xi*eta2;
  phi[4][1][1] = -.4*eta+5.3*xi*eta+2.2*eta2-9*xi2*eta-8.1*xi*eta2-1.8*eta3;
  phi[5][0][0] = 0;
  phi[5][0][1] = 0;
  phi[5][1][0] = 0;
  phi[5][1][1] = xi+eta-6*xi2-xi*eta-eta2+6*xi3;
  phi[6][0][0] = xi+eta-xi2-xi*eta-6*eta2+6*eta3;
  phi[6][0][1] = 0;
  phi[6][1][0] = 0;
  phi[6][1][1] = 0;
  phi[7][0][0] = -.1*(45*eta+18*xi-4)*xi*(-1+2*eta+xi);
  phi[7][0][1] = .1*eta*(30*eta2+81*xi*eta+54*xi2+5-30*eta-40*xi);
  phi[7][1][0] = .1*eta*(30*eta2+81*xi*eta+54*xi2+5-30*eta-40*xi);
  phi[7][1][1] = -.1*(54*xi*eta+28*xi-13*eta+27*eta2-14)*eta;
  phi[8][0][0] = .6*xi+1.2*xi2-1.2*xi*eta-1.8*xi3-3.6*xi2*eta;
  phi[8][0][1] = -3*xi*eta+5.4*xi2*eta+3.6*xi*eta2;
  phi[8][1][0] = -3*xi*eta+5.4*xi2*eta+3.6*xi*eta2;
  phi[8][1][1] = 1.4*eta-2.8*xi*eta+.8*eta2-5.4*xi*eta2-1.2*eta3;
  phi[9][0][0] = 3*xi*(9*xi*eta*omega-6*xi*omega+6*xi2*omega+eta*omega+1-2*eta-xi);
  phi[9][0][1] = -3*xi*(6*xi2*omega-11*eta*omega+3*omega+18*xi*eta*omega+9*eta2*omega-9*xi*omega-2*eta+1-xi);
  phi[9][1][0] = -3*xi*(6*xi2*omega-11*eta*omega+3*omega+18*xi*eta*omega+9*eta2*omega-9*xi*omega-2*eta+1-xi);
  phi[9][1][1] = 3*eta*(-18*xi*omega-7*eta*omega+3*eta2*omega+18*xi2*omega+4*omega+18*xi*eta*omega+1-2*xi-eta);
  phi[10][0][0] = 1.8*omega*(2*eta-2*xi-1+6*xi*eta+3*xi2)*xi;
  phi[10][0][1] = -1.8*(-5+9*xi+6*eta)*xi*eta*omega;
  phi[10][1][0] = -1.8*(-5+9*xi+6*eta)*xi*eta*omega;
  phi[10][1][1] = 3*eta-9*xi*omega+27*xi2*omega+4.8*eta*omega-8.4*eta2*omega+3.6*omega*eta3-18*xi3*omega-6.6*xi*eta*omega+16.2*omega*xi*eta2-3*xi-3*eta2+3*xi2;
  phi[11][0][0] = -3*xi*(9*xi*eta*omega-6*xi*omega+6*xi2*omega+eta*omega-1+2*eta+xi);
  phi[11][0][1] = 3*xi*(18*xi*eta*omega-11*eta*omega+6*xi2*omega+3*omega+9*eta2*omega-9*xi*omega+2*eta-1+xi);
  phi[11][1][0] = 3*xi*(18*xi*eta*omega-11*eta*omega+6*xi2*omega+3*omega+9*eta2*omega-9*xi*omega+2*eta-1+xi);
  phi[11][1][1] = -3*eta*(4*omega-7*eta*omega+18*xi*eta*omega-18*xi*omega+18*xi2*omega+3*eta2*omega+2*xi+eta-1);
  phi[12][0][0] = -1.8*omega*(2*eta-2*xi-1+6*xi*eta+3*xi2)*xi;
  phi[12][0][1] = 1.8*(-5+9*xi+6*eta)*xi*eta*omega;
  phi[12][1][0] = 1.8*(-5+9*xi+6*eta)*xi*eta*omega;
  phi[12][1][1] = 3*eta+6.6*xi*eta*omega+18*xi3*omega+9*xi*omega-27*xi2*omega-4.8*eta*omega+8.4*eta2*omega-3.6*omega*eta3-16.2*omega*xi*eta2-3*xi-3*eta2+3*xi2;
  phi[13][0][0] = .6*sqrt2*xi*(-2*omega+6*xi2*omega-4*xi*omega-eta*omega+27*xi*eta*omega-5+5*xi+10*eta);
  phi[13][0][1] = -1.8*sqrt2*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[13][1][0] = -1.8*sqrt2*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[13][1][1] = 1.8*sqrt2*eta*(-1+2*xi-2*eta+6*xi*eta+3*eta2)*omega;
  phi[14][0][0] = -1.8*sqrt2*xi*(2*eta-2*xi-1+6*xi*eta+3*xi2)*omega;
  phi[14][0][1] = 1.8*(-5+9*xi+6*eta)*sqrt2*xi*eta*omega;
  phi[14][1][0] = 1.8*(-5+9*xi+6*eta)*sqrt2*xi*eta*omega;
  phi[14][1][1] = -.6*sqrt2*eta*(-2*omega+27*xi*eta*omega-4*eta*omega+6*eta2*omega-xi*omega-10*xi-5*eta+5);
  phi[15][0][0] = -.6*sqrt2*xi*(-eta*omega-2*omega+6*xi2*omega-4*xi*omega+27*xi*eta*omega+5-5*xi-10*eta);
  phi[15][0][1] = 1.8*sqrt2*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[15][1][0] = 1.8*sqrt2*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[15][1][1] = -1.8*sqrt2*eta*(-1+2*xi-2*eta+6*xi*eta+3*eta2)*omega;
  phi[16][0][0] = 1.8*sqrt2*xi*(2*eta-2*xi-1+6*xi*eta+3*xi2)*omega;
  phi[16][0][1] = -1.8*(-5+9*xi+6*eta)*sqrt2*xi*eta*omega;
  phi[16][1][0] = -1.8*(-5+9*xi+6*eta)*sqrt2*xi*eta*omega;
  phi[16][1][1] = .6*sqrt2*eta*(-2*omega+27*xi*eta*omega-4*eta*omega+6*eta2*omega-xi*omega+10*xi+5*eta-5);
  phi[17][0][0] = -4.8*xi*omega+9*eta*omega+8.4*xi2*omega+6.6*xi*eta*omega-27*eta2*omega-3.6*xi3*omega-16.2*xi2*eta*omega+18*omega*eta3+3*xi-3*eta-3*xi2+3*eta2;
  phi[17][0][1] = 1.8*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[17][1][0] = 1.8*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[17][1][1] = -1.8*eta*(-1+2*xi-2*eta+6*xi*eta+3*eta2)*omega;
  phi[18][0][0] = -3*xi*(18*eta2*omega+4*omega+3*xi2*omega-7*xi*omega-18*eta*omega+18*xi*eta*omega-1+2*eta+xi);
  phi[18][0][1] = 3*eta*(3*omega-9*eta*omega+9*xi2*omega-11*xi*omega+6*eta2*omega+18*xi*eta*omega+2*xi+eta-1);
  phi[18][1][0] = 3*eta*(3*omega-9*eta*omega+9*xi2*omega-11*xi*omega+6*eta2*omega+18*xi*eta*omega+2*xi+eta-1);
  phi[18][1][1] = -3*eta*(-6*eta*omega+9*xi*eta*omega+xi*omega+6*eta2*omega+2*xi+eta-1);
  phi[19][0][0] = 4.8*xi*omega-9*eta*omega-8.4*xi2*omega-6.6*xi*eta*omega+27*eta2*omega+3.6*xi3*omega+16.2*xi2*eta*omega-18*omega*eta3+3*xi-3*eta-3*xi2+3*eta2;
  phi[19][0][1] = -1.8*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[19][1][0] = -1.8*xi*eta*(-5+6*xi+9*eta)*omega;
  phi[19][1][1] = 1.8*eta*(-1+2*xi-2*eta+6*xi*eta+3*eta2)*omega;
  phi[20][0][0] = 3*xi*(4*omega-7*xi*omega+18*xi*eta*omega-18*eta*omega+3*xi2*omega+18*eta2*omega+1-2*eta-xi);
  phi[20][0][1] = -3*eta*(3*omega-11*xi*omega+18*xi*eta*omega-9*eta*omega+9*xi2*omega+6*eta2*omega+1-2*xi-eta);
  phi[20][1][0] = -3*eta*(3*omega-11*xi*omega+18*xi*eta*omega-9*eta*omega+9*xi2*omega+6*eta2*omega+1-2*xi-eta);
  phi[20][1][1] = 3*eta*(xi*omega+6*eta2*omega-6*eta*omega+9*xi*eta*omega+1-2*xi-eta);
  phi[21][0][0] = 9*xi-9*xi2-9*xi*eta;
  phi[21][0][1] = 0;
  phi[21][1][0] = 0;
  phi[21][1][1] = 0;
  phi[22][0][0] = 4.5*xi-4.5*xi2-9*xi*eta;
  phi[22][0][1] = 4.5*xi*eta;
  phi[22][1][0] = 4.5*xi*eta;
  phi[22][1][1] = 4.5*eta-9*xi*eta-4.5*eta2;
  phi[23][0][0] = 0;
  phi[23][0][1] = 0;
  phi[23][1][0] = 0;
  phi[23][1][1] = 9*eta-9*xi*eta-9*eta2;
  
  if(dphi == NULL) return;

  dphi[0] = svector<2>(.6-2.6*xi+.8*eta, .6*eta-.6+1.2*xi);
  dphi[1] = svector<2>(-.5*xi+.5*eta, .5*xi-.5*eta);
  dphi[2] = svector<2>(-.6+.6*xi+1.2*eta, -2.6*eta+.6+.8*xi);
  dphi[3] = svector<2>(1.4-1.4*xi-2.8*eta, -.6*eta+.6-1.2*xi);
  dphi[4] = svector<2>(1.4-1.4*xi-2.8*eta, .1-.7*xi+.4*eta);
  dphi[5] = svector<2>(0, 1-xi-2*eta);
  dphi[6] = svector<2>(1-2*xi-eta, 0);
  dphi[7] = svector<2>(.4*xi-.7*eta+.1, -1.4*eta+1.4-2.8*xi);
  dphi[8] = svector<2>(.6-.6*xi-1.2*eta, -1.4*eta+1.4-2.8*xi);
  dphi[9] = svector<2>(-3*xi*omega+3*eta*omega+3-6*eta, -3*omega*(3*eta-1));
  dphi[10] = svector<2>(1.8*omega*(-1+2*eta+xi), -7.8*eta*omega+3+4.8*omega-6.6*xi*omega-6*eta);
  dphi[11] = svector<2>(3*xi*omega-3*eta*omega+3-6*eta, 3*omega*(3*eta-1));
  dphi[12] = svector<2>(-1.8*omega*(-1+2*eta+xi), 7.8*eta*omega+3+6.6*xi*omega-4.8*omega-6*eta);
  dphi[13] = svector<2>(-.6*sqrt2*(2*omega-7*xi*omega+eta*omega+5-10*xi-10*eta), 1.8*sqrt2*omega*(-1+2*xi+eta));
  dphi[14] = svector<2>(-1.8*sqrt2*omega*(-1+2*eta+xi), -.6*sqrt2*(7*eta*omega-2*omega-xi*omega-10*xi-10*eta+5));
  dphi[15] = svector<2>(.6*sqrt2*(eta*omega+2*omega-7*xi*omega-5+10*xi+10*eta), -1.8*sqrt2*omega*(-1+2*xi+eta));
  dphi[16] = svector<2>(1.8*sqrt2*omega*(-1+2*eta+xi), .6*sqrt2*(7*eta*omega-2*omega-xi*omega+10*xi+10*eta-5));
  dphi[17] = svector<2>(-4.8*omega+7.8*xi*omega+6.6*eta*omega+3-6*xi, -1.8*omega*(-1+2*xi+eta));
  dphi[18] = svector<2>(3*omega*(-1+3*xi), 3*eta*omega-3*xi*omega-6*xi+3);
  dphi[19] = svector<2>(4.8*omega-7.8*xi*omega-6.6*eta*omega+3-6*xi, 1.8*omega*(-1+2*xi+eta));
  dphi[20] = svector<2>(-3*omega*(-1+3*xi), -3*eta*omega+3*xi*omega+3-6*xi);
  dphi[21] = svector<2>(9-18*xi-9*eta, 0);
  dphi[22] = svector<2>(4.5-4.5*xi-9*eta, -4.5*eta+4.5-9*xi);
  dphi[23] = svector<2>(0, 9-9*xi-18*eta);
  
  return;
}

void awTriangle1Signs(const int nodes[3], int signs[24])
{
  for(int i=0; i < 24; i++) signs[i] = 1;   // only edge dof may change sign 

  // Convention: global edge orientation is low to high

  if(nodes[0] > nodes[1]) signs[9 ]=signs[10]=signs[11]=signs[12] = -1;
  if(nodes[1] > nodes[2]) signs[13]=signs[14]=signs[15]=signs[16] = -1;
  if(nodes[2] > nodes[0]) signs[17]=signs[18]=signs[19]=signs[20] = -1;

  return;
}

void awTriangle1Permutation(const int nodes[3], int dof[24])
{
  int i, temp, d0;

  for(i = 0; i < 3; i++)
    {
      if(nodes[i] > nodes[(i+1)%3])
	{
	d0 = 9 + 4*i; temp = dof[d0]; dof[d0] = dof[d0+2]; dof[d0+2] = temp;
	d0 ++;        temp = dof[d0]; dof[d0] = dof[d0+2]; dof[d0+2] = temp;
	}
    }

  return;
}

void awTriangle1Transform(int nodes[3], smatrix<2> &dxds,
		 smatrix<2> Phi[24], svector<2> Dphi[24], 
		 smatrix<2> phi[24], svector<2> dphi[24])
{
  int i,j,m, alpha,beta, mij, mab, mi, ma, signs[24];
  double det, temp;

  smatrix<2> dxdsT = dxds.transpose();
  smatrix<2> dsdx  = dxds.inverse(det);
  smatrix<2> dsdxT = dsdx.transpose();

  if(det < 0) std::cout << "awTriangle1Transform(): Warning det < 0" << "\n";

  for(m = 0; m < 3; m++)  // each vertex
    {
    for(i = 0; i <  2; i++)
    for(j = 0; j <= i; j++)
      {
      mij = 3*m + (i*(i+1))/2 + j;

      phi[mij].zero();
      dphi[mij].zero();

      for(alpha = 0; alpha < 2; alpha++)
      for(beta  = 0; beta <= alpha; beta++)
	{
	mab = 3*m + (alpha*(alpha+1))/2 + beta;

	temp = dsdx[alpha][i]*dsdx[beta][j] + dsdx[beta][i]*dsdx[alpha][j];

	if(alpha == beta) temp *= 0.5;
	
	phi[mij]  += temp * dxds*Phi[mab]*dxdsT;
	dphi[mij] += temp * dxds*Dphi[mab];
	}
      }
    }

  double fmtn[3];  // |F^{-T} nhat| for nomrals to parent edges

  awTriangle1Signs(nodes, signs);  // only edge dof's change with normal
 
  fmtn[0] = (dsdxT * svector<2>(0.0,-1.0)).norm();
  fmtn[1] = (dsdxT * svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0))).norm();
  fmtn[2] = (dsdxT * svector<2>(-1.0,0.0)).norm();
  
  for(m = 0; m < 6; m++)  // 2 per edge
    {
    for(i = 0; i < 2; i++)
      {
      mi = 9 + 2*m + i;

      phi[mi].zero();
      dphi[mi].zero();

      for(alpha = 0; alpha < 2; alpha++)
	{
	ma = 9 + 2*m + alpha;

	temp = signs[mi] * dsdx[alpha][i] * fmtn[m/2];

	phi[mi]  += temp * dxds*Phi[ma]*dxdsT;
	dphi[mi] += temp * dxds*Dphi[ma];
	}
      }
    }
    
  for(i = 0; i <  2; i++)     // centroid
  for(j = 0; j <= i; j++)
    {
    mij = 21 + (i*(i+1))/2 + j;

    phi[mij].zero();
    dphi[mij].zero();

    for(alpha = 0; alpha < 2; alpha++)
    for(beta  = 0; beta <= alpha; beta++)
      {
      mab = 21 + (alpha*(alpha+1))/2 + beta;

      temp = dsdx[alpha][i]*dsdx[beta][j] + dsdx[beta][i]*dsdx[alpha][j];
      
      if(alpha == beta) temp *= 0.5;
      
      phi[mij]  += temp * dxds*Phi[mab]*dxdsT;
      dphi[mij] += temp * dxds*Dphi[mab];
      }
    }

  return;
}

void awTriangle1Interpolate(int nodes[3],  smatrix<2> &dxds,
			    smatrix<2> ss[10], double dof[24])
{
  // ss[] is the values of S  at the pts[]
  // compute the corresponding dof[] for the interpolant

  int i,j, m, mm, mij, mi, signs[24];
  double det;

  smatrix<2> dsdx  = dxds.inverse(det);
  smatrix<2> dsdxT = dsdx.transpose();

  awTriangle1Signs(nodes, signs);  // only edge dof's change with normal

  for(m = 0; m < 3; m++)  // each vertex
    {
    for(i = 0; i <  2; i++)
    for(j = 0; j <= i; j++)
      {
      mij = 3*m + (i*(i+1))/2 + j;

      dof[mij] = ss[m][i][j];

      if(i != j) dof[mij] = 2 * dof[mij];  // off diagonals
      }
    }
  
  // compute normals for edge dof
    
  svector<2> nn[3];

  nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1.0/ nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1.0/ nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1.0/ nn[2].norm()) * nn[2];
  
  for(m = 0; m < 6; m++)  // 2 per edge
    {
    mm = m/2;
      
    svector<2> sn = ss[3+m] * nn[mm];
    
    for(i = 0; i < 2; i++)
      {
      mi = 9 + 2*m + i;

      dof[mi] = signs[mi] * sn[i];
      }
    }

  // centroid

  for(i = 0; i <  2; i++)
  for(j = 0; j <= i; j++)
    {
    mij = 21 + (i*(i+1))/2 + j;

    dof[mij] = ss[9][i][j];

    if(i != j) dof[mij] = 2 * dof[mij];  // off diagonals
    }
  
  return;
}
