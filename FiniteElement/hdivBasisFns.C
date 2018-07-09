#include<vector>

// H(div) baisis function
// functions vector valued, only derivative is the (scalar) divergence
// their boundarys have scalar valued basis functions (normal component)

template<const int Ndim, const int Nchi>
struct HdivBdyBasisFunction : public basisFunction<Ndim>
{
  HdivBdyBasisFunction(basisFunction<Ndim> &bf, 
     void (*bb)(svector<Ndim>, double *phi),
     void (*ss)(const int *nodes, int &sign),
     void (*ii)(BoundaryGelement<Ndim+1, Nchi> &ge, 
		vector<double> &bPts, double *ue))
    : basisFunction<Ndim>(bf), basis(bb), dofSign(ss), interpolate(ii)  
  {return;}

  void (*basis)(svector<Ndim>, double *phi);
  void (*dofSign)(const int *nodes, int &sign);
  void (*interpolate)(BoundaryGelement<Ndim+1, Nchi> &ge, 
		      vector<double> &bPts, double *ue);
};

template<const int Ndim, const int Nchi, const int Nchib>
struct HdivBasisFunction : public basisFunction<Ndim>
{
  HdivBasisFunction(basisFunction<Ndim> &bf, 
       HdivBdyBasisFunction<Ndim-1, Nchib> *bdy,
       void (*bb)(svector<Ndim>, svector<Ndim> *phi, double *dphi),
       void (*ss)(const int *nodes, int *signs),
       void (*tt)(GeometricElement<Ndim, Nchi> &ge, smatrix<Ndim> dxds,
		  svector<Ndim> *Phi, double *Dphi,
		  svector<Ndim> *phi, double *dphi),
       void (*ii)(GeometricElement<Ndim, Nchi> &ge,
		  vector<svector<Ndim> > &pPts, double *ue))
    : basisFunction<Ndim>(bf), boundary(bdy), basis(bb), dofSigns(ss),
      transform(tt), interpolate(ii)
  {return;}

  HdivBdyBasisFunction<Ndim-1, Nchib> *boundary;
  void (*basis)(svector<Ndim>, svector<Ndim> *phi, double *dphi);
  void (*dofSigns)(const int *nodes, int *signs);
  void (*transform)(GeometricElement<Ndim, Nchi> &ge, smatrix<Ndim> dxds,
		    svector<Ndim> *Phi, double *Dphi,
		    svector<Ndim> *phi, double *dphi);
  void (*interpolate)(GeometricElement<Ndim, Nchi> &ge, 
	     vector<svector<Ndim> > &pPts, double *ue);
};

void rtTriangle0Bdy(svector<1>, double *phi);
void rtTriangleBdySign(const int *nodes, int &sign);
void rtTriangle0BdyInterpolate(BoundaryGelement<2,2> &ge,
			       vector<double> &bPts, double ue[1]);

void rtTriangle0(svector<2> xx, svector<2> *phi, double *dphi);
void rtTriangle0Signs(const int *nodes, int *signs);
void rtTriangle0Transform(GeometricElement<2,3> &ge, smatrix<2> dxds,
		    svector<2> Phi[3], double Dphi[3],
		    svector<2> phi[3], double dphi[3]);
void rtTriangle0Interpolate(GeometricElement<2,3> &ge,
			    vector<svector<2> > &pPts, double ue[3]);

void rtTriangle1Bdy(svector<1>, double phi[2]);
void rtTriangle1BdyPermutation(const int *nodes, int *dof);
void rtTriangle1BdyInterpolate(BoundaryGelement<2,2> &ge,
			       vector<double> &bPts, double ue[3]);

void rtTriangle1(svector<2> xx, svector<2> phi[8], double dphi[8]);
void rtTriangle1Signs(const int *nodes, int signs[8]);
void rtTriangle1Permutation(const int *nodes, int *dof);
void rtTriangle1Transform(GeometricElement<2,3> &ge, smatrix<2> dxds,
		    svector<2> Phi[8], double Dphi[8],
		    svector<2> phi[8], double dphi[8]);
void rtTriangle1Interpolate(GeometricElement<2,3> &ge,
			    vector<svector<2> > &pPts, double ue[8]);

void bdmTriangle1(svector<2> xx, svector<2> phi[6], double dphi[6]);
void bdmTriangle1Signs(const int *nodes, int signs[6]);
void bdmTriangle1Transform(GeometricElement<2,3> &ge, smatrix<2> dxds,
			  svector<2> Phi[6], double Dphi[6],
			   svector<2> phi[6], double dphi[6]);
void bdmTriangle1Interpolate(GeometricElement<2,3> &ge,
			     vector<svector<2> > &pPts, double ue[6]);

void bdfmSquare0(svector<2> xi, svector<2> phi[4], double *dphi);
void bdfmSquare0Signs(const int nodes[4], int signs[4]);
void bdfmSquare0Transform(GeometricElement<2,4> &ge, smatrix<2> dxds,
		    svector<2> Phi[4], double Dphi[4],
		    svector<2> phi[4], double dphi[4]);
void bdfmSquare0Interpolate(GeometricElement<2,4> &ge,
			    vector<svector<2> > &pPts, double ue[4]);

void rtTetBdySign(const int *nodes, int &sign);

void rtTet0Bdy(svector<2>, double *phi);
void rtTet0BdyInterpolate(BoundaryGelement<3,3> &ge,
			  vector<double> &bPts, double ue[1]);

void rtTet0(svector<3> xx, svector<3> *phi, double *dphi);
void rtTet0Signs(const int *nodes, int *signs);
void rtTet0Transform(GeometricElement<3,4> &ge, smatrix<3> dxds,
		     svector<3> Phi[4], double Dphi[4],
		     svector<3> phi[4], double dphi[4]);
void rtTet0Interpolate(GeometricElement<3,4> &ge,
		       vector<svector<3> > &pPts, double ue[4]);

void rtTet1Bdy(svector<2>, double phi[3]);
void rtTet1BdyInterpolate(BoundaryGelement<3,3> &ge,
			  vector<double> &bPts, double ue[3]);

void rtTet1(svector<3> xx, svector<3> phi[15], double dphi[15]);
void rtTet1Signs(const int nodes[4], int signs[15]);
void rtTet1Permutation(const int nodes[4], int dof[15]);
void rtTet1Transform(GeometricElement<3,4> &ge, smatrix<3> dxds,
		     svector<3> Phi[15], double Dphi[15],
		     svector<3> phi[15], double dphi[15]);
void rtTet1Interpolate(GeometricElement<3,4> &ge,
		       vector<svector<3> > &pPts, double ue[15]);

void bdmTet1Bdy(svector<2>, double phi[3]);
void bdmTet1BdyInterpolate(BoundaryGelement<3,3> &ge,
			  vector<double> &bPts, double ue[3]);

void bdmTet1(svector<3> xx, svector<3> phi[12], double dphi[12]);
void bdmTet1Signs(const int nodes[4], int signs[12]);
void bdmTet1Permutation(const int nodes[4], int dof[12]);
void bdmTet1Transform(GeometricElement<3,4> &ge, smatrix<3> dxds,
		      svector<3> Phi[12], double Dphi[12],
		      svector<3> phi[12], double dphi[12]);
void bdmTet1Interpolate(GeometricElement<3,4> &ge,
			vector<svector<3> > &pPts, double ue[12]);

void bdfmCube0Bdy(svector<2> xi, double *phi);
void bdfmCubeBdySign(const int *nodes, int &sign);
void bdfmCube0BdyInterpolate(BoundaryGelement<3,4> &ge,
			     vector<double> &bPts, double ue[1]);

void bdfmCube0(svector<3> xi, svector<3> phi[6], double *dphi);
void bdfmCube0Signs(const int nodes[8], int signs[6]);
void bdfmCube0Transform(GeometricElement<3,8> &ge, smatrix<3> dxds,
			  svector<3> Phi[6], double Dphi[6],
			  svector<3> phi[6], double dphi[6]);
void bdfmCube0Interpolate(GeometricElement<3,8> &ge,
			  vector<svector<3> > &pPts, double ue[4]);

// Construct the RT degree zero triangle
// Construct the RT degree zero triangle
// Construct the RT degree zero triangle

svector<1> RTtriangle0BdyPts[] = {svector<1>(0.0)};

basisFunction<1> RTtriangle0BdyBF =    // RT_0 boundary
  {1,
   {0,1},
   EmptyVertex,
   &Interval,
   1,
   RTtriangle0BdyPts,
   identityPermutation
  };

HdivBdyBasisFunction<1,2> RTtriangle0BdyClass(
   RTtriangle0BdyBF, &rtTriangle0Bdy, &rtTriangleBdySign,
   &rtTriangle0BdyInterpolate),
  *RTtriangle0Bdy = &RTtriangle0BdyClass;

svector<2> RTtriangle0Pts[] = 
  {svector<2>(0.5, 0.0), svector<2>(0.5, 0.5), svector<2>(0.0, 0.5)};

basisFunction<2> RTtriangle0BF =    // RT_0 element
  {3,
   {0,1,0},
   &RTtriangle0BdyBF,
   &Triangle,
   3,
   RTtriangle0Pts,
   identityPermutation
  };

HdivBasisFunction<2,3,2> RTtriangle0class(RTtriangle0BF, 
   RTtriangle0Bdy, &rtTriangle0, &rtTriangle0Signs,
   &rtTriangle0Transform, &rtTriangle0Interpolate),
  *RTtriangle0 = &RTtriangle0class;

// Construct the RT degree one triangle
// Construct the RT degree one triangle
// Construct the RT degree one triangle

svector<1> RTtriangle1BdyPts[] = 
  {svector<1>(-1.0/sqrt(3.0)), svector<1>(1.0/sqrt(3.0))};

basisFunction<1> RTtriangle1BdyBF =    // RT_1 boundary
  {2,
   {0,2},
   EmptyVertex,
   &Interval,
   2,
   RTtriangle1BdyPts,
   rtTriangle1BdyPermutation
  };

HdivBdyBasisFunction<1,2> RTtriangle1BdyClass(
   RTtriangle1BdyBF, &rtTriangle1Bdy, &rtTriangleBdySign,
   &rtTriangle1BdyInterpolate),
  *RTtriangle1Bdy = &RTtriangle1BdyClass;

svector<2> RTtriangle1Pts[] = 
  {svector<2>(0.5*(1.0 - 1.0/sqrt(3.0)), 0), 
   svector<2>(0.5*(1.0 + 1.0/sqrt(3.0)), 0),
   svector<2>(0.5*(1.0 + 1.0/sqrt(3.0)), 0.5*(1.0 - 1.0/sqrt(3.0))),
   svector<2>(0.5*(1.0 - 1.0/sqrt(3.0)), 0.5*(1.0 + 1.0/sqrt(3.0))),
   svector<2>(0, 0.5*(1.0 + 1.0/sqrt(3.0))),
   svector<2>(0, 0.5*(1.0 - 1.0/sqrt(3.0))), 
   svector<2>(1/3.0, 1/3.0),  svector<2>(1/3.0, 1/3.0)}; // last is redundent

basisFunction<2> RTtriangle1BF =    // RT_1 element
  {8,
   {0,2,2},
   &RTtriangle1BdyBF,
   &Triangle,
   7,
   RTtriangle1Pts,
   rtTriangle1Permutation
  };

HdivBasisFunction<2,3,2> RTtriangle1class(RTtriangle1BF, 
   RTtriangle1Bdy, &rtTriangle1, &rtTriangle1Signs,
   &rtTriangle1Transform, &rtTriangle1Interpolate),
  *RTtriangle1 = &RTtriangle1class;

// Construct the BDM degree one triangle 
// Construct the BDM degree one triangle
// Construct the BDM degree one triangle

basisFunction<2> BDMtriangle1BF =    // BDM_1 element
  {6,
   {0,2,0},
   &RTtriangle1BdyBF,
   &Triangle,
   6,
   RTtriangle1Pts,
   rtTriangle1Permutation
  };

HdivBasisFunction<2,3,2> BDMtriangle1class(BDMtriangle1BF, 
   RTtriangle1Bdy, &bdmTriangle1, &bdmTriangle1Signs,
   &bdmTriangle1Transform, &bdmTriangle1Interpolate),
  *BDMtriangle1 = &BDMtriangle1class;


// ***************** BDFM type elements on square ***********************
// ***************** BDFM type elements on square ***********************
// ***************** BDFM type elements on square ***********************

// basisFunction<1> BDFMsquare0Bdy = RTtriangleBdy0

svector<2> BDFMsquare0Pts[] = 
  {svector<2>(0.0, -1.0), svector<2>( 1.0, 0.0), 
   svector<2>(0.0,  1.0), svector<2>(-1.0, 0.0)};

basisFunction<2> BDFMsquare0BF =   
  {4,
   {0,1,0},
   &RTtriangle0BdyBF,  
   &Square,
   4,
   BDFMsquare0Pts,
   identityPermutation
  };

HdivBasisFunction<2,4,2> BDFMsquare0class(BDFMsquare0BF, 
   RTtriangle0Bdy, &bdfmSquare0, &bdfmSquare0Signs,
   &bdfmSquare0Transform, &bdfmSquare0Interpolate),				  *BDFMsquare0 = &BDFMsquare0class;

// Construct the RT-0 degree zero tetrahedra
// Construct the RT-0 degree zero tetrahedra
// Construct the RT-0 degree zero tetrahedra

svector<2> RTtet0BdyPts[] = {svector<2>(1.0/3.0, 1.0/3.0)};

basisFunction<2> RTtet0BdyBF =    // RT_0 boundary
  {1,
   {0,0,1},
   EmptyInterval,
   &Triangle,
   1,
   RTtet0BdyPts,
   identityPermutation
  };

HdivBdyBasisFunction<2,3> RTtet0BdyClass(
   RTtet0BdyBF, &rtTet0Bdy, &rtTetBdySign, rtTet0BdyInterpolate),
  *RTtet0Bdy = &RTtet0BdyClass;

svector<3> RTtet0Pts[] = 
  {svector<3>(1.0/3.0, 1.0/3.0, 1.0/3.0), 
   svector<3>(0.0,     1.0/3.0, 1.0/3.0),
   svector<3>(1.0/3.0, 0.0,     1.0/3.0),
   svector<3>(1.0/3.0, 1.0/3.0, 0.0    )};

basisFunction<3> RTtet0BF =    // RT_0 element
  {4,
   {0,0,1,0},
   &RTtet0BdyBF,
   &Tetrahedron,
   4,
   RTtet0Pts,
   identityPermutation
  };

HdivBasisFunction<3,4,3> RTtet0class(
   RTtet0BF, RTtet0Bdy, &rtTet0, &rtTet0Signs,
   &rtTet0Transform, &rtTet0Interpolate),
  *RTtet0 = &RTtet0class;

// Construct the RT-1 degree one tetrahedra
// Construct the RT-1 degree one tetrahedra
// Construct the RT-1 degree one tetrahedra

svector<2> RTtet1BdyPts[4] = 
  {
    svector<2>(1.0/5.0, 1.0/5.0), 
    svector<2>(3.0/5.0, 1.0/5.0),
    svector<2>(1.0/5.0, 3.0/5.0),

    svector<2>(1.0/3.0, 1.0/3.0) 
  };

basisFunction<2> RTtet1BdyBF =    // RT_1 boundary
  {3,
   {0,0,3},
   EmptyInterval,
   &Triangle,
   4,
   RTtet1BdyPts,
   discontinuousTriangle3Permutation
  };

HdivBdyBasisFunction<2,3> RTtet1BdyClass(
   RTtet1BdyBF, &rtTet1Bdy, &rtTetBdySign, rtTet1BdyInterpolate),
  *RTtet1Bdy = &RTtet1BdyClass;

svector<3> RTtet1Pts[17] = 
  {
    svector<3>(3.0/5, 1.0/5, 1.0/5),
    svector<3>(1.0/5, 3.0/5, 1.0/5),
    svector<3>(1.0/5, 1.0/5, 3.0/5),
    svector<3>(1.0/3, 1.0/3, 1.0/3),

    svector<3>(0, 1.0/5, 1.0/5),
    svector<3>(0, 1.0/5, 3.0/5),
    svector<3>(0, 3.0/5, 1.0/5),
    svector<3>(0, 1.0/3, 1.0/3),

    svector<3>(1.0/5, 0, 1.0/5),
    svector<3>(3.0/5, 0, 1.0/5),
    svector<3>(1.0/5, 0, 3.0/5),
    svector<3>(1.0/3, 0, 1.0/3),

    svector<3>(1.0/5, 1.0/5, 0),
    svector<3>(1.0/5, 3.0/5, 0),
    svector<3>(3.0/5, 1.0/5, 0),
    svector<3>(1.0/3, 1.0/3, 0),

    svector<3>(1.0/4, 1.0/4, 1.0/4)
  };

basisFunction<3> RTtet1BF =    // RT_1 element
  {15,
   {0,0,3,3},
   &RTtet1BdyBF,
   &Tetrahedron,
   17,
   RTtet1Pts,
   rtTet1Permutation
  };

HdivBasisFunction<3,4,3> RTtet1class(
   RTtet1BF, RTtet1Bdy, &rtTet1, &rtTet1Signs,
   &rtTet1Transform, &rtTet1Interpolate),
  *RTtet1 = &RTtet1class;

// Construct the BDM-1 degree one tetrahedra
// Construct the BDM-1 degree one tetrahedra
// Construct the BDM-1 degree one tetrahedra

svector<2> BDMtet1BdyPts[3] = 
  {
    svector<2>(1.0/6.0, 1.0/6.0), 
    svector<2>(2.0/3.0, 1.0/6.0),
    svector<2>(1.0/6.0, 2.0/3.0),
  };

basisFunction<2> BDMtet1BdyBF =    // BDM_1 boundary
  {3,
   {0,0,3},
   EmptyInterval,
   &Triangle,
   3,
   BDMtet1BdyPts,
   discontinuousTriangle3Permutation
  };

HdivBdyBasisFunction<2,3> BDMtet1BdyClass(
   BDMtet1BdyBF, &bdmTet1Bdy, &rtTetBdySign, bdmTet1BdyInterpolate),
  *BDMtet1Bdy = &BDMtet1BdyClass;

svector<3> BDMtet1Pts[12] = 
  {
    svector<3>(2.0/3, 1.0/6, 1.0/6),
    svector<3>(1.0/6, 2.0/3, 1.0/6),
    svector<3>(1.0/6, 1.0/6, 2.0/3),
    svector<3>(0, 1.0/6, 1.0/6),
    svector<3>(0, 1.0/6, 2.0/3),
    svector<3>(0, 2.0/3, 1.0/6),
    svector<3>(1.0/6, 0, 1.0/6),
    svector<3>(2.0/3, 0, 1.0/6),
    svector<3>(1.0/6, 0, 2.0/3),
    svector<3>(1.0/6, 1.0/6, 0),
    svector<3>(1.0/6, 2.0/3, 0),
    svector<3>(2.0/3, 1.0/6, 0)
  };

basisFunction<3> BDMtet1BF =    // BDM_1 element
  {12,
   {0,0,3,0},
   &BDMtet1BdyBF,
   &Tetrahedron,
   12,
   BDMtet1Pts,
   bdmTet1Permutation
  };

HdivBasisFunction<3,4,3> BDMtet1class(
   BDMtet1BF, BDMtet1Bdy, &bdmTet1, &bdmTet1Signs,
   &bdmTet1Transform, &bdmTet1Interpolate),
  *BDMtet1 = &BDMtet1class;

// **************** BDFM-0 Cube  ************************
// **************** BDFM-0 Cube  ************************
// **************** BDFM-0 Cube  ************************

svector<2> BDFMcube0BdyPts[] = {svector<2>(0.0,0.0)};

basisFunction<2> BDFMcube0BdyBF =    // BDFM_0 boundary
  {1,
   {0,0,1},
   EmptyInterval,
   &Square,
   1,
   BDFMcube0BdyPts,
   identityPermutation
  };

HdivBdyBasisFunction<2,4> BDFMcube0BdyClass(
   BDFMcube0BdyBF, &bdfmCube0Bdy, &bdfmCubeBdySign,
   &bdfmCube0BdyInterpolate),
  *BDFMcube0Bdy = &BDFMcube0BdyClass;

svector<3> BDFMcube0Pts[] = 
  {svector<3>( 0, 0,-1), svector<3>(0, 0, 1),
   svector<3>( 0,-1, 0), svector<3>(0, 1, 0),
   svector<3>(-1, 0, 0), svector<3>(1, 0, 0)};

basisFunction<3> BDFMcube0BF =    // BDFM_0 element
  {6,
   {0,0,1,0},
   &BDFMcube0BdyBF,
   &Cube,
   6,
   BDFMcube0Pts,
   identityPermutation
  };

HdivBasisFunction<3,8,4> BDFMcube0class(
   BDFMcube0BF, BDFMcube0Bdy, &bdfmCube0, &bdfmCube0Signs,
   &bdfmCube0Transform, &bdfmCube0Interpolate),
  *BDFMcube0 = &BDFMcube0class;


// ********************* Basis Functions ************************
// ********************* Basis Functions ************************
// ********************* Basis Functions ************************

void rtTriangle0Bdy(svector<1> xi, double *phi)
{
  phi[0] = 1.0;
}

void rtTriangleBdySign(const int *nodes, int &sign)
{
  // Convention: global edge orientation is low to high

  if(nodes[0] > nodes[1]) sign = -1; else sign = 1;
}

void rtTriangle0BdyInterpolate(BoundaryGelement<2,2> &ge,
			       vector<double> &bPts, double ue[1])
{
  int ss;
  rtTriangleBdySign(ge.nodes, ss);

  ue[0] = ss * bPts[0];
}

void rtTriangle0(svector<2> xi, svector<2> *phi, double *dphi = NULL)
{
  double eta = xi[0], ata = xi[1];

 phi[0] = svector<2>(eta, ata-1);
 phi[1] = svector<2>(eta * sqrt(2.0), ata * sqrt(2.0));
 phi[2] = svector<2>(eta-1, ata);

 if(dphi == NULL) return;

 dphi[0] = 2.0;
 dphi[1] = 2.0 * sqrt(2);
 dphi[2] = 2.0;
 
 return;
}

void rtTriangle0Signs(const int nodes[3], int signs[3])
{
  // Convention: global edge orientation is low to high

  // master edges  {{0,1},{1,2},{2,0}},  (better is edge opposite 0, 1, 2)

  if(nodes[0] > nodes[1]) signs[0] = -1; else signs[0] = 1;
  if(nodes[1] > nodes[2]) signs[1] = -1; else signs[1] = 1;
  if(nodes[2] > nodes[0]) signs[2] = -1; else signs[2] = 1;

  return;
}


void rtTriangle0Transform(GeometricElement<2,3> &ge, smatrix<2> dxds,
		    svector<2> Phi[3], double Dphi[3],
		    svector<2> phi[3], double dphi[3])
{
  int i, signs[3];
  
  rtTriangle0Signs(ge.nodes, signs);  

  // Jacobian is constant

  double det;
  
  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();

  svector<2> nn[3];
  
  nn[0] = svector<2>(0.0,-1.0); 
  nn[1] = svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0)); 
  nn[2] = svector<2>(-1.0,0.0); 
  
  for(i = 0; i < 3; i++)
    {
      double modFn = signs[i] * (dsdxT * nn[i]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  return;
}

void rtTriangle0Interpolate(GeometricElement<2,3> &ge,
			      vector<svector<2> > &pPts, double ue[3])
{
  int i, signs[3];
  
  rtTriangle0Signs(ge.nodes, signs);  

  double det;
  smatrix<2> dxds;

  ge.map(svector<2>(1.0/3.0,1.0/3.0), dxds);   // Jacobian is constant
  
  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();

  svector<2> nn[3];
  
  nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1/nn[2].norm()) * nn[2];
  
  for(i = 0; i < 3; i++) ue[i] = signs[i] * pPts[i].dot(nn[i]);

  return;
}

// *********************** RT-1 Triangle ************************
// *********************** RT-1 Triangle ************************
// *********************** RT-1 Triangle ************************

void rtTriangle1Bdy(svector<1> xi, double phi[2])
{
  double gp = 1.0/sqrt(3.0); // (-gp, gp) = Gauss points on [-1,1]

  phi[0] = (gp-xi[0])/(2*gp);
  phi[1] = (gp+xi[0])/(2*gp);
}

void rtTriangle1BdyPermutation(const int *nodes, int *dof)
{
 if(nodes[0] > nodes[1]) 
   {
     int temp = dof[1]; dof[1] = dof[0]; dof[0] = temp;
   }
}

void rtTriangle1BdyInterpolate(BoundaryGelement<2,2> &ge,
			       vector<double> &bPts, double ue[2])
{
  int ss;
  rtTriangleBdySign(ge.nodes, ss);   

  for(int i = 0; i < 2; i++) ue[i] = ss * bPts[i];

  return;
}

void rtTriangle1(svector<2> pxi, svector<2> phi[8], double *dphi)
{
  double x = pxi[0], y = pxi[1], sqrt3 = sqrt(3.0), sqrt2 = sqrt(2.0);

  phi[0] = svector<2>(-(1.0/6)*sqrt3*(-sqrt3-5+8*x+4*y+4*y*sqrt3)*x, -(1.0/6)*sqrt3*(sqrt3+3-6*x-5*y*sqrt3-7*y+8*x*y+4*y*y+4*y*y*sqrt3));
  phi[1] = svector<2>(-(1.0/6)*sqrt3*(-sqrt3+5-8*x-4*y+4*y*sqrt3)*x, -(1.0/6)*sqrt3*(sqrt3-3+6*x-5*y*sqrt3+7*y-8*x*y-4*y*y+4*y*y*sqrt3));
  phi[2] = svector<2>(-(1.0/6)*sqrt2*(3+sqrt3)*(-4*x+4-sqrt3-8*y+4*y*sqrt3)*x, -(1.0/6)*sqrt2*(3+sqrt3)*(-8*y+4*y*sqrt3+5-2*sqrt3-4*x)*y);
  phi[3] = svector<2>(-(1.0/6)*sqrt2*(sqrt3-3)*(4*x-sqrt3-4+8*y+4*y*sqrt3)*x, -(1.0/6)*sqrt2*(sqrt3-3)*(8*y+4*y*sqrt3-2*sqrt3-5+4*x)*y);
  phi[4] = svector<2>(-(1.0/6)*(sqrt3-3)*(sqrt3-sqrt3*x+4*x-3*y-4*x*x+4*x*y-3*y*sqrt3+4*x*y*sqrt3), -(1.0/6)*(sqrt3-3)*(-2*sqrt3-1+4*y+4*y*sqrt3-4*x)*y);
  phi[5] = svector<2>(-(1.0/6)*(3+sqrt3)*(sqrt3-sqrt3*x-4*x+3*y+4*x*x-4*x*y-3*y*sqrt3+4*x*y*sqrt3), -(1.0/6)*(3+sqrt3)*(-2*sqrt3+1-4*y+4*y*sqrt3+4*x)*y);
  phi[6] = svector<2>(-4*x*(-2+2*x+y), -4*y*(y-1+2*x));
  phi[7] = svector<2>(-4*x*(-1+x+2*y), -4*y*(2*y-2+x));
  
  if(dphi == NULL) return;
  
  dphi[0] = -sqrt3*(4*x-sqrt3-2+2*y+2*y*sqrt3);
  dphi[1] = -sqrt3*(-4*x-sqrt3+2-2*y+2*y*sqrt3);
  dphi[2] = -(1.0/2)*sqrt2*(3+sqrt3)*(-sqrt3+4*y*sqrt3-4*x-8*y+3);
  dphi[3] = -(1.0/2)*sqrt2*(sqrt3-3)*(-sqrt3+4*y*sqrt3+4*x+8*y-3);
  dphi[4] = -(1.0/2)*(sqrt3-3)*(-sqrt3+4*y*sqrt3+1-4*x+4*y);
  dphi[5] = -(1.0/2)*(3+sqrt3)*(-sqrt3+4*y*sqrt3-1+4*x-4*y);
  dphi[6] = 12-24*x-12*y;
  dphi[7] = -24*y+12-12*x;

  return;
}

void rtTriangle1Signs(const int *nodes, int signs[8])
{
  // Convention: global edge orientation is low to high

  if(nodes[0] > nodes[1]) signs[0]=signs[1] = -1; else signs[0]=signs[1] = 1;
  if(nodes[1] > nodes[2]) signs[2]=signs[3] = -1; else signs[2]=signs[3] = 1;
  if(nodes[2] > nodes[0]) signs[4]=signs[5] = -1; else signs[4]=signs[5] = 1;

  signs[6]=signs[7] = 1;

  return;
}

void rtTriangle1Permutation(const int *nodes, int *dof)
{
  int temp;

  if(nodes[0] > nodes[1]) {temp = dof[0]; dof[0] = dof[1]; dof[1] = temp;}
  if(nodes[1] > nodes[2]) {temp = dof[2]; dof[2] = dof[3]; dof[3] = temp;}
  if(nodes[2] > nodes[0]) {temp = dof[4]; dof[4] = dof[5]; dof[5] = temp;}

  return;
}

void rtTriangle1Interpolate(GeometricElement<2,3> &ge,
			      vector<svector<2> > &pPts, double ue[8])
{
  int i, signs[8];
  double det;
  
  rtTriangle1Signs(ge.nodes, signs);  

  smatrix<2> dxds;

  ge.map(svector<2>(1.0/3.0,1.0/3.0), dxds);   // Jacobian is constant

  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();

  svector<2> nn[3], temp;
  
  nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1/nn[2].norm()) * nn[2];

  temp.zero();
  
  for(i = 0; i < 6; i++) 
    {
      ue[i] = signs[i] * pPts[i].dot(nn[i/2]);

      temp = temp + pPts[i];
    }

  temp = (1.0/12.0) * temp + 0.5 * pPts[6];

  for(i = 0; i < 2; i++) ue[i+6] = temp[i];

  return;
}

void rtTriangle1Transform(GeometricElement<2,3> &ge, smatrix<2> dxds,
			  svector<2> Phi[8], double Dphi[8],
			  svector<2> phi[8], double dphi[8])
{
  int i, signs[8];
  
  rtTriangle1Signs(ge.nodes, signs);  

  // Jacobian is constant

  double det;
  
  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();

  svector<2> nn[3];
  
  nn[0] = svector<2>(0.0,-1.0); 
  nn[1] = svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0)); 
  nn[2] = svector<2>(-1.0,0.0); 
  
  for(i = 0; i < 6; i++)
    {
      double modFn = signs[i] * (dsdxT * nn[i/2]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  svector<2> temp;

  for(i = 0; i < 2; i++)
    {
      temp.zero();
      dphi[6+i] = 0.0;

      for(int alpha = 0; alpha < 2; alpha++) 
	{
	  temp = temp + dsdxT[i][alpha] * Phi[6+alpha];

	  dphi[6+i] += dsdxT[i][alpha] * Dphi[6+alpha];
	}

      phi[6+i] = dxds * temp;
    }

  return;
}

// *********************** BDM-1 Triangle ************************
// *********************** BDM-1 Triangle ************************
// *********************** BDM-1 Triangle ************************

void bdmTriangle1(svector<2> pxi, svector<2> phi[6], double *dphi = NULL)
{
  double x = pxi[0], y = pxi[1], sqrt3 = sqrt(3.0), sqrt2 = sqrt(2.0);

  double c0 = 0.5*(1-sqrt3);
  double c1 = 0.5*(1+sqrt3);

  phi[0] = svector<2>(c0*x, -c1+sqrt3*x+c1*y);
  phi[1] = svector<2>(c1*x, -c0-sqrt3*x+c0*y);
  phi[2] = svector<2>(c1*x*sqrt2, c0*y*sqrt2);
  phi[3] = svector<2>(c0*x*sqrt2, c1*y*sqrt2);
  phi[4] = svector<2>(-c0+c0*x-sqrt3*y, c1*y);
  phi[5] = svector<2>(-c1+c1*x+sqrt3*y, c0*y);
  
  if(dphi == NULL) return;
  
  dphi[0] = 1;
  dphi[1] = 1;
  dphi[2] = sqrt2;
  dphi[3] = sqrt2;
  dphi[4] = 1;
  dphi[5] = 1;
  
  return;
}

void bdmTriangle1Signs(const int *nodes, int signs[6])
{
  // Convention: global edge orientation is low to high

  if(nodes[0] > nodes[1]) signs[0]=signs[1] = -1; else signs[0]=signs[1] = 1;
  if(nodes[1] > nodes[2]) signs[2]=signs[3] = -1; else signs[2]=signs[3] = 1;
  if(nodes[2] > nodes[0]) signs[4]=signs[5] = -1; else signs[4]=signs[5] = 1;

  return;
}

void bdmTriangle1Transform(GeometricElement<2,3> &ge, smatrix<2> dxds,
			  svector<2> Phi[6], double Dphi[6],
			  svector<2> phi[6], double dphi[6])
{
  int i, signs[6];
  
  bdmTriangle1Signs(ge.nodes, signs);  

  // Jacobian is constant

  double det;
  
  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();

  svector<2> nn[3];
  
  nn[0] = svector<2>(0.0,-1.0); 
  nn[1] = svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0)); 
  nn[2] = svector<2>(-1.0,0.0); 
  
  for(i = 0; i < 6; i++)
    {
      double modFn = signs[i] * (dsdxT * nn[i/2]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  return;
}

void bdmTriangle1Interpolate(GeometricElement<2,3> &ge,
			      vector<svector<2> > &pPts, double ue[6])
{
  int i, signs[6];
  double det;
  
  bdmTriangle1Signs(ge.nodes, signs);  

  smatrix<2> dxds;

  ge.map(svector<2>(1.0/3.0,1.0/3.0), dxds);   // Jacobian is constant
  
  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();

  svector<2> nn[3];
  
  nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1/nn[2].norm()) * nn[2];
  
  for(i = 0; i < 6; i++)  ue[i] = signs[i] * pPts[i].dot(nn[i/2]);

  return;
}

// ************************** RT0 Tetrahedron ************************
// ************************** RT0 Tetrahedron ************************
// ************************** RT0 Tetrahedron ************************

void rtTet0Bdy(svector<2> xi, double *phi)
{
  phi[0] = 1.0;
}

void rtTetBdySign(const int *nodes, int &sign)
{

  if( ((nodes[0] < nodes[1]) && (nodes[1] < nodes[2])) ||
      ((nodes[2] < nodes[0]) && (nodes[0] < nodes[1])) ||
      ((nodes[1] < nodes[2]) && (nodes[2] < nodes[0])) ) 
    sign = 1;
  else 
    sign = -1;

  return;
}

void rtTet0BdyInterpolate(BoundaryGelement<3,3> &ge,
			  vector<double> &bPts, double ue[1])
{
  int ss;
  rtTetBdySign(ge.nodes, ss);

  ue[0] = ss * bPts[0];
}

void rtTet0(svector<3> xi, svector<3> *phi, double *dphi = NULL)
{
  double eta = xi[0], ata = xi[1], chi = xi[2], sqrt3 = sqrt(3.0);

  phi[0] = sqrt3 * svector<3>(eta,  ata  , chi);
  phi[1] = svector<3>(eta-1,ata  , chi);
  phi[2] = svector<3>(eta  ,ata-1, chi);
  phi[3] = svector<3>(eta  ,ata  , chi-1);

 if(dphi == NULL) return;

 dphi[0] = 3.0 * sqrt3;
 dphi[1] = 3.0;
 dphi[2] = 3.0;
 dphi[3] = 3.0;
 
 return;
}

void rtTet0Signs(const int nodes[4], int signs[4])
{
  // master faces are {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}},

  if( ((nodes[1] < nodes[2]) && (nodes[2] < nodes[3])) ||
      ((nodes[2] < nodes[3]) && (nodes[3] < nodes[1])) ||
      ((nodes[3] < nodes[1]) && (nodes[1] < nodes[2])) )
    signs[0] = 1;
  else
    signs[0] = -1;

  if( ((nodes[0] < nodes[3]) && (nodes[3] < nodes[2])) ||
      ((nodes[3] < nodes[2]) && (nodes[2] < nodes[0])) ||
      ((nodes[2] < nodes[0]) && (nodes[0] < nodes[3])) )
    signs[1] = 1;
  else
    signs[1] = -1;

  if( ((nodes[0] < nodes[1]) && (nodes[1] < nodes[3])) ||
      ((nodes[1] < nodes[3]) && (nodes[3] < nodes[0])) ||
      ((nodes[3] < nodes[0]) && (nodes[0] < nodes[1])) )
    signs[2] = 1;
  else
    signs[2] = -1;

  if( ((nodes[0] < nodes[2]) && (nodes[2] < nodes[1])) ||
      ((nodes[2] < nodes[1]) && (nodes[1] < nodes[0])) ||
      ((nodes[1] < nodes[0]) && (nodes[0] < nodes[2])) )
    signs[3] = 1;
  else
    signs[3] = -1;

  return;
}

void rtTet0Transform(GeometricElement<3,4> &ge, smatrix<3> dxds,
		    svector<3> Phi[4], double Dphi[4],
		    svector<3> phi[4], double dphi[4])
{
  int i, signs[4];
  
  rtTet0Signs(ge.nodes, signs);  

  double det;
  
  smatrix<3> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<3> nn[4];
  
  nn[0] = (1.0/sqrt(3.0)) * svector<3>( 1.0, 1.0, 1.0);
  nn[1] = svector<3>(-1.0, 0.0, 0.0);
  nn[2] = svector<3>( 0.0,-1.0, 0.0);
  nn[3] = svector<3>( 0.0, 0.0,-1.0);
  
  for(i = 0; i < 4; i++)
    {
      double modFn = signs[i] * (dsdxT * nn[i]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  return;
}

void rtTet0Interpolate(GeometricElement<3,4> &ge,
		       vector<svector<3> > &pPts, double ue[4])
{
  int i, signs[4];
  
  rtTet0Signs(ge.nodes, signs);  

  double det;
  smatrix<3> dxds;

  ge.map(svector<3>(1.0/4.0, 1.0/4.0, 1.0/4.0), dxds); // Jacobian is constant
  
  smatrix<3> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<3> nn[4];
  
  nn[0] = dsdxT * svector<3>( 1.0, 1.0, 1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<3>(-1.0, 0.0, 0.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<3>( 0.0,-1.0, 0.0); nn[2] = (1/nn[2].norm()) * nn[2];
  nn[3] = dsdxT * svector<3>( 0.0, 0.0,-1.0); nn[3] = (1/nn[3].norm()) * nn[3];
  
  for(i = 0; i < 4; i++) ue[i] = signs[i] * pPts[i].dot(nn[i]);

  return;
}

// ************************** RT-1 Tetrahedron ************************
// ************************** RT-1 Tetrahedron ************************
// ************************** RT-1 Tetrahedron ************************

void rtTet1Bdy(svector<2> xi, double phi[3])
{
  phi[0] = 2 - 2.5*xi[0] - 2.5*xi[1];
  phi[1] = 2.5*xi[0] - 0.5;
  phi[2] = 2.5*xi[1] - 0.5;

  return;
}

void rtTet1BdyInterpolate(BoundaryGelement<3,3> &ge,
			  vector<double> &bPts, double ue[3])
{
  int j, ss;

  assert(bPts.size() == 4);

  rtTetBdySign(ge.nodes, ss);

  double w0 = (3.0/16.0) * (bPts[0] + bPts[1] + bPts[2]) 
            - (9.0/16.0) *  bPts[3];

  for(j = 0; j < 3; j++) ue[j] = ss * (bPts[j] + w0); 

  return;
}

void rtTet1(svector<3> xx, svector<3> psi[15], double dpsi[15])
{
  double x = xx[0], y = xx[1], z = xx[2], sqrt3 = sqrt(3.0);
  double xsq = x*x, ysq = y*y, zsq = z*z;

  psi[0] = svector<3>((1.0/8)*sqrt3*(-14+30*x+5*z+5*y)*x, (1.0/8)*sqrt3*(-9+5*z+5*y+30*x)*y, (1.0/8)*sqrt3*(-9+5*z+5*y+30*x)*z);
  psi[1] = svector<3>((1.0/8)*sqrt3*(-9+5*z+5*x+30*y)*x, (1.0/8)*sqrt3*(-14+30*y+5*z+5*x)*y, (1.0/8)*sqrt3*(-9+5*z+5*x+30*y)*z);
  psi[2] = svector<3>((1.0/8)*sqrt3*(-9+5*x+30*z+5*y)*x, (1.0/8)*sqrt3*(-9+5*x+30*z+5*y)*y, (1.0/8)*sqrt3*(-14+30*z+5*y+5*x)*z);
  psi[3] = svector<3>(-2+(23.0/4)*x+(5.0/2)*y+(5.0/2)*z-(15.0/4)*xsq-(25.0/8)*x*y-(25.0/8)*x*z, -(1.0/8)*y*(25*y-21+30*x+25*z), -(1.0/8)*z*(25*y-21+30*x+25*z));
  psi[4] = svector<3>(1.0/2+(1.0/8)*x-(5.0/2)*z-(5.0/8)*xsq+(25.0/8)*x*z, -(1.0/8)*y*(4+5*x-25*z), -(1.0/8)*z*(-25*z+9+5*x));
  psi[5] = svector<3>(1.0/2+(1.0/8)*x-(5.0/2)*y-(5.0/8)*xsq+(25.0/8)*x*y, -(1.0/8)*y*(-25*y+9+5*x), -(1.0/8)*z*(4+5*x-25*y));
  psi[6] = svector<3>(-(1.0/8)*x*(-21+25*x+30*y+25*z), -2+(5.0/2)*x+(23.0/4)*y+(5.0/2)*z-(25.0/8)*x*y-(15.0/4)*ysq-(25.0/8)*y*z, -(1.0/8)*z*(-21+25*x+30*y+25*z));
  psi[7] = svector<3>((1.0/8)*x*(-9+25*x-5*y), 1.0/2-(5.0/2)*x+(1.0/8)*y+(25.0/8)*x*y-(5.0/8)*ysq, (1.0/8)*z*(-4+25*x-5*y));
  psi[8] = svector<3>(-(1.0/8)*x*(4+5*y-25*z), 1.0/2+(1.0/8)*y-(5.0/2)*z-(5.0/8)*ysq+(25.0/8)*y*z, -(1.0/8)*z*(-25*z+9+5*y));
  psi[9] = svector<3>(-(1.0/8)*x*(-21+25*x+25*y+30*z), -(1.0/8)*y*(-21+25*x+25*y+30*z), -2+(5.0/2)*x+(5.0/2)*y+(23.0/4)*z-(25.0/8)*x*z-(25.0/8)*y*z-(15.0/4)*zsq);
  psi[10] = svector<3>((1.0/8)*x*(-4+25*y-5*z), (1.0/8)*y*(-9+25*y-5*z), 1.0/2-(5.0/2)*y+(1.0/8)*z+(25.0/8)*y*z-(5.0/8)*zsq);
  psi[11] = svector<3>((1.0/8)*x*(-9+25*x-5*z), (1.0/8)*y*(-4+25*x-5*z), 1.0/2-(5.0/2)*x+(1.0/8)*z+(25.0/8)*x*z-(5.0/8)*zsq);
  psi[12] = svector<3>(-5*x*(-2+2*x+y+z), -5*y*(y-1+2*x+z), -5*z*(y-1+2*x+z));
  psi[13] = svector<3>(-5*x*(-1+x+2*y+z), -5*y*(2*y-2+x+z), -5*z*(-1+x+2*y+z));
  psi[14] = svector<3>(-5*x*(-1+x+y+2*z), -5*y*(-1+x+y+2*z), -5*z*(2*z-2+x+y));
  
  if(dpsi == NULL) return;
  
  dpsi[0] = (15.0/4)*sqrt3*x+(1.0/8)*sqrt3*(-14+30*x+5*z+5*y)+(5.0/8)*sqrt3*y+(1.0/4)*sqrt3*(-9+5*z+5*y+30*x)+(5.0/8)*sqrt3*z;
  dpsi[1] = (5.0/8)*sqrt3*x+(1.0/4)*sqrt3*(-9+5*z+5*x+30*y)+(15.0/4)*sqrt3*y+(1.0/8)*sqrt3*(-14+30*y+5*z+5*x)+(5.0/8)*sqrt3*z;
  dpsi[2] = (5.0/8)*sqrt3*x+(1.0/4)*sqrt3*(-9+5*x+30*z+5*y)+(5.0/8)*sqrt3*y+(15.0/4)*sqrt3*z+(1.0/8)*sqrt3*(-14+30*z+5*y+5*x);
  dpsi[3] = 11-15*x-(25.0/2)*y-(25.0/2)*z;
  dpsi[4] = -3.0/2-(5.0/2)*x+(25.0/2)*z;
  dpsi[5] = -3.0/2-(5.0/2)*x+(25.0/2)*y;
  dpsi[6] = 11-(25.0/2)*x-15*y-(25.0/2)*z;
  dpsi[7] = -3.0/2+(25.0/2)*x-(5.0/2)*y;
  dpsi[8] = -3.0/2-(5.0/2)*y+(25.0/2)*z;
  dpsi[9] = 11-(25.0/2)*x-(25.0/2)*y-15*z;
  dpsi[10] = -3.0/2+(25.0/2)*y-(5.0/2)*z;
  dpsi[11] = -3.0/2+(25.0/2)*x-(5.0/2)*z;
  dpsi[12] = 20-40*x-20*y-20*z;
  dpsi[13] = -40*y+20-20*x-20*z;
  dpsi[14] = -20*y+20-20*x-40*z;
  
  return;
}

void rtTet1Signs(const int nodes[4], int signs[15])
{
  // master faces are {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}},

  if( ((nodes[1] < nodes[2]) && (nodes[2] < nodes[3])) ||
      ((nodes[2] < nodes[3]) && (nodes[3] < nodes[1])) ||
      ((nodes[3] < nodes[1]) && (nodes[1] < nodes[2])) )
    signs[0] = signs[1] = signs[2] = 1;
  else
    signs[0] = signs[1] = signs[2] = -1;

  if( ((nodes[0] < nodes[3]) && (nodes[3] < nodes[2])) ||
      ((nodes[3] < nodes[2]) && (nodes[2] < nodes[0])) ||
      ((nodes[2] < nodes[0]) && (nodes[0] < nodes[3])) )
    signs[3] = signs[4] = signs[5] = 1;
  else
    signs[3] = signs[4] = signs[5] = -1;

  if( ((nodes[0] < nodes[1]) && (nodes[1] < nodes[3])) ||
      ((nodes[1] < nodes[3]) && (nodes[3] < nodes[0])) ||
      ((nodes[3] < nodes[0]) && (nodes[0] < nodes[1])) )
    signs[6] = signs[7] = signs[8] = 1;
  else
    signs[6] = signs[7] = signs[8] = -1;

  if( ((nodes[0] < nodes[2]) && (nodes[2] < nodes[1])) ||
      ((nodes[2] < nodes[1]) && (nodes[1] < nodes[0])) ||
      ((nodes[1] < nodes[0]) && (nodes[0] < nodes[2])) )
    signs[9] = signs[10] = signs[11] = 1;
  else
    signs[9] = signs[10] = signs[11] = -1;

  signs[12] = signs[13]= signs[14] = 1;

  return;
}

void rtTet1Permutation(const int nodes[4], int dof[15])
{
  int nn[3];     // sort face dof's

  for(int i = 0; i < 4; i++) 
    {
      for(int j = 0; j < 3; j++) nn[j] = nodes[Tetrahedron.subs(2,i,j)];

      sortDof(3, nn, dof+3*i);
    }

  return;
}

void rtTet1Interpolate(GeometricElement<3,4> &ge,
		       vector<svector<3> > &pPts, double ue[15])
{
  int i,j, signs[15];

  assert(pPts.size() == 17);
  
  rtTet1Signs(ge.nodes, signs);  

  double det;
  smatrix<3> dxds;

  ge.map(svector<3>(1.0/4.0, 1.0/4.0, 1.0/4.0), dxds); // Jacobian is constant
  
  smatrix<3> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<3> cc, nn[4];
  
  nn[0] = dsdxT * svector<3>( 1.0, 1.0, 1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<3>(-1.0, 0.0, 0.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<3>( 0.0,-1.0, 0.0); nn[2] = (1/nn[2].norm()) * nn[2];
  nn[3] = dsdxT * svector<3>( 0.0, 0.0,-1.0); nn[3] = (1/nn[3].norm()) * nn[3];
  
  cc.zero();

  for(i = 0; i < 4; i++)   // each face
    {
      double pn[4];

      for(j = 0; j < 4; j++) pn[j] = pPts[4*i+j].dot(nn[i]);

      double w0 = (3.0/16.0) * (pn[0] + pn[1] + pn[2]) - (9.0/16.0) * pn[3];

      for(j = 0; j < 3; j++) 
	{
	  ue[3*i+j] = signs[3*i+j] * (pn[j] + w0); 

	  cc = cc + pPts[4*i + j];
	}
    }

  cc = (5.0/76.0) * cc + (4.0/19.0) * pPts[16];

  for(i = 0; i < 3; i++) ue[i+12] = cc[i];

  return;
}

void rtTet1Transform(GeometricElement<3,4> &ge, smatrix<3> dxds,
		    svector<3> Phi[15], double Dphi[15],
		    svector<3> phi[15], double dphi[15])
{
  int i, signs[15];
  
  rtTet1Signs(ge.nodes, signs);  

  double det;
  
  smatrix<3> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<3> nn[4];
  
  nn[0] = svector<3>( 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0));
  nn[1] = svector<3>(-1.0, 0.0, 0.0);
  nn[2] = svector<3>( 0.0,-1.0, 0.0);
  nn[3] = svector<3>( 0.0, 0.0,-1.0);
  
  for(i = 0; i < 12; i++)
    {
      double modFn = signs[i] * (dsdxT * nn[i/3]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  svector<3> temp;

  for(i = 0; i < 3; i++)
    {
      temp.zero();
      dphi[12+i] = 0.0;

      for(int alpha = 0; alpha < 3; alpha++) 
	{
	  temp = temp + dsdxT[i][alpha] * Phi[12+alpha];

	  dphi[12+i] += dsdxT[i][alpha] * Dphi[12+alpha];
	}

      phi[12+i] = dxds * temp;
    }

  return;
}

// ************************** BDM-1 Tetrahedron ************************
// ************************** BDM-1 Tetrahedron ************************
// ************************** BDM-1 Tetrahedron ************************

void bdmTet1Bdy(svector<2> xi, double phi[3])
{
  phi[0] = 5.0/3.0 - 2*xi[0] - 2*xi[1];
  phi[1] = 2*xi[0] - 1.0/3.0;
  phi[2] = 2*xi[1] - 1.0/3.0;

  return;
}

void bdmTet1BdyInterpolate(BoundaryGelement<3,3> &ge,
			   vector<double> &bPts, double ue[3])
{
  int j, ss;

  assert(bPts.size() == 3);

  rtTetBdySign(ge.nodes, ss);

  for(j = 0; j < 3; j++) ue[j] = ss * bPts[j]; 

  return;
}

void bdmTet1(svector<3> xx, svector<3> psi[12], double dpsi[12])
{
  double x = xx[0], y = xx[1], z = xx[2], sqrt3 = sqrt(3.0);

  psi[0] = svector<3>((5.0/3)*sqrt3*x, -(1.0/3)*sqrt3*y, -(1.0/3)*sqrt3*z);
  psi[1] = svector<3>(-(1.0/3)*sqrt3*x, (5.0/3)*sqrt3*y, -(1.0/3)*sqrt3*z);
  psi[2] = svector<3>(-(1.0/3)*sqrt3*x, -(1.0/3)*sqrt3*y, (5.0/3)*sqrt3*z);
  psi[3] = svector<3>(-5.0/3+(5.0/3)*x+2*y+2*z, -(1.0/3)*y, -(1.0/3)*z);
  psi[4] = svector<3>(1.0/3-(1.0/3)*x-2*z, -(1.0/3)*y, (5.0/3)*z);
  psi[5] = svector<3>(1.0/3-(1.0/3)*x-2*y, (5.0/3)*y, -(1.0/3)*z);
  psi[6] = svector<3>(-(1.0/3)*x, -5.0/3+2*x+(5.0/3)*y+2*z, -(1.0/3)*z);
  psi[7] = svector<3>((5.0/3)*x, 1.0/3-2*x-(1.0/3)*y, -(1.0/3)*z);
  psi[8] = svector<3>(-(1.0/3)*x, 1.0/3-(1.0/3)*y-2*z, (5.0/3)*z);
  psi[9] = svector<3>(-(1.0/3)*x, -(1.0/3)*y, -5.0/3+2*x+2*y+(5.0/3)*z);
  psi[10] = svector<3>(-(1.0/3)*x, (5.0/3)*y, 1.0/3-2*y-(1.0/3)*z);
  psi[11] = svector<3>((5.0/3)*x, -(1.0/3)*y, 1.0/3-2*x-(1.0/3)*z);
  
  if(dpsi == NULL) return;
  
  dpsi[0] = sqrt3;
  dpsi[1] = sqrt3;
  dpsi[2] = sqrt3;
  dpsi[3] = 1;
  dpsi[4] = 1;
  dpsi[5] = 1;
  dpsi[6] = 1;
  dpsi[7] = 1;
  dpsi[8] = 1;
  dpsi[9] = 1;
  dpsi[10] = 1;
  dpsi[11] = 1;
  
  return;
}

void bdmTet1Signs(const int nodes[4], int signs[12])
{
  // master faces are {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}},

  if( ((nodes[1] < nodes[2]) && (nodes[2] < nodes[3])) ||
      ((nodes[2] < nodes[3]) && (nodes[3] < nodes[1])) ||
      ((nodes[3] < nodes[1]) && (nodes[1] < nodes[2])) )
    signs[0] = signs[1] = signs[2] = 1;
  else
    signs[0] = signs[1] = signs[2] = -1;

  if( ((nodes[0] < nodes[3]) && (nodes[3] < nodes[2])) ||
      ((nodes[3] < nodes[2]) && (nodes[2] < nodes[0])) ||
      ((nodes[2] < nodes[0]) && (nodes[0] < nodes[3])) )
    signs[3] = signs[4] = signs[5] = 1;
  else
    signs[3] = signs[4] = signs[5] = -1;

  if( ((nodes[0] < nodes[1]) && (nodes[1] < nodes[3])) ||
      ((nodes[1] < nodes[3]) && (nodes[3] < nodes[0])) ||
      ((nodes[3] < nodes[0]) && (nodes[0] < nodes[1])) )
    signs[6] = signs[7] = signs[8] = 1;
  else
    signs[6] = signs[7] = signs[8] = -1;

  if( ((nodes[0] < nodes[2]) && (nodes[2] < nodes[1])) ||
      ((nodes[2] < nodes[1]) && (nodes[1] < nodes[0])) ||
      ((nodes[1] < nodes[0]) && (nodes[0] < nodes[2])) )
    signs[9] = signs[10] = signs[11] = 1;
  else
    signs[9] = signs[10] = signs[11] = -1;

  return;
}

void bdmTet1Permutation(const int nodes[4], int dof[12])
{
  int nn[3];     // sort face dof's

  for(int i = 0; i < 4; i++) 
    {
      for(int j = 0; j < 3; j++) nn[j] = nodes[Tetrahedron.subs(2,i,j)];

      sortDof(3, nn, dof+3*i);
    }

  return;
}

void bdmTet1Interpolate(GeometricElement<3,4> &ge,
			vector<svector<3> > &pPts, double ue[12])
{
  int i, signs[12];

  assert(pPts.size() == 12);
  
  bdmTet1Signs(ge.nodes, signs);  

  double det;
  smatrix<3> dxds;

  ge.map(svector<3>(1.0/4.0, 1.0/4.0, 1.0/4.0), dxds); // Jacobian is constant
  
  smatrix<3> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<3> nn[4];
  
  nn[0] = dsdxT * svector<3>( 1.0, 1.0, 1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<3>(-1.0, 0.0, 0.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<3>( 0.0,-1.0, 0.0); nn[2] = (1/nn[2].norm()) * nn[2];
  nn[3] = dsdxT * svector<3>( 0.0, 0.0,-1.0); nn[3] = (1/nn[3].norm()) * nn[3];
  
  for(i = 0; i < 12; i++) ue[i] = signs[i] * pPts[i].dot(nn[i/3]);

  return;
}

void bdmTet1Transform(GeometricElement<3,4> &ge, smatrix<3> dxds,
		      svector<3> Phi[12], double Dphi[12],
		      svector<3> phi[12], double dphi[12])
{
  int i, signs[12];
  
  bdmTet1Signs(ge.nodes, signs);  

  double det;
  
  smatrix<3> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<3> nn[4];
  
  nn[0] = svector<3>( 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0));
  nn[1] = svector<3>(-1.0, 0.0, 0.0);
  nn[2] = svector<3>( 0.0,-1.0, 0.0);
  nn[3] = svector<3>( 0.0, 0.0,-1.0);
  
  for(i = 0; i < 12; i++)
    {
      double modFn = signs[i] * (dsdxT * nn[i/3]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  return;
}


// ************************* BDFM 0 square **********************
// ************************* BDFM 0 square **********************
// ************************* BDFM 0 square **********************

void bdfmSquare0(svector<2> xi, svector<2> phi[4], double *dphi = NULL)
{
  double eta = xi[0], ata = xi[1];

  phi[0] = svector<2>(0, (ata-1)/2);
  phi[1] = svector<2>((eta+1)/2, 0);
  phi[2] = svector<2>(0, (ata+1)/2);
  phi[3] = svector<2>((eta-1)/2, 0);
 
  if(dphi == NULL) return;

  dphi[0] = 0.5;
  dphi[1] = 0.5;
  dphi[2] = 0.5;
  dphi[3] = 0.5;
  
  return;
}

void bdfmSquare0Signs(const int nodes[4], int signs[4])
{
  // Convention: global edge orientation is low to high

  // master edges {{0,1},{1,2},{2,3},{3,0}}, 

  if(nodes[0] > nodes[1]) signs[0] = -1; else signs[0] = 1;
  if(nodes[1] > nodes[2]) signs[1] = -1; else signs[1] = 1;
  if(nodes[2] > nodes[3]) signs[2] = -1; else signs[2] = 1;
  if(nodes[3] > nodes[0]) signs[3] = -1; else signs[3] = 1;

  return;
}

void bdfmSquare0Interpolate(GeometricElement<2,4> &ge,
			    vector<svector<2> > &pPts, double ue[4])
{
  int i, signs[4];
  
  bdfmSquare0Signs(ge.nodes, signs);  

  double det;
  smatrix<2> dxds, dsdxT;

  svector<2> nhat[4], nn;
  
  nhat[0] = svector<2>( 0.0,-1.0);
  nhat[1] = svector<2>( 1.0, 0.0);
  nhat[2] = svector<2>( 0.0, 1.0);
  nhat[3] = svector<2>(-1.0, 0.0);
  
  for(i = 0; i < 4; i++) 
    {
      ge.map(BDFMsquare0Pts[i], dxds);
  
      dsdxT = (dxds.inverse(det)).transpose();

      nn = dsdxT * nhat[i]; nn = (1.0/nn.norm()) * nn;

      ue[i] = signs[i] * pPts[i].dot(nn);
    }

  return;
}

void bdfmSquare0Transform(GeometricElement<2,4> &ge, smatrix<2> dxds,
			   svector<2> Phi[4], double Dphi[4],
			   svector<2> phi[4], double dphi[4])
{
  int i, signs[4];
  
  bdfmSquare0Signs(ge.nodes, signs);  

  svector<2> nhat[4], nn;
  
  nhat[0] = svector<2>( 0.0,-1.0);
  nhat[1] = svector<2>( 1.0, 0.0);
  nhat[2] = svector<2>( 0.0, 1.0);
  nhat[3] = svector<2>(-1.0, 0.0);
 
  double det = fabs(dxds.determinant());

  smatrix<2> dxdsBdy, dxdsBdyCof;

  for(i = 0; i < 4; i++)
    {
      ge.map(BDFMsquare0Pts[i], dxdsBdy);
  
      dxdsBdyCof = dxdsBdy.cofactor();

      double modFn = (signs[i]/det) * (dxdsBdyCof * nhat[i]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  return;
}

// ************************* BDFM cubes **********************
// ************************* BDFM cubes **********************
// ************************* BDFM cubes **********************

void bdfmCube0Bdy(svector<2> xi, double *phi)
{
  phi[0] = 1.0;
}

void bdfmCubeBdySign(const int *nodes, int &sign)
{
  int i0 = 0, ii[3];

  // find the lowest index and use the corresponding triangle

  for(int i = 1; i < 4; i++) if(nodes[i] < nodes[i0]) i0 = i;

  ii[0] = nodes[(i0-1+4)%4]; ii[1] = nodes[i0]; ii[2] = nodes[(i0+1)%4];

  rtTetBdySign(ii, sign);

  return; 
}

void bdfmCube0BdyInterpolate(BoundaryGelement<3,4> &ge,
			    vector<double> &bPts, double ue[1])
{
  int ss;
  bdfmCubeBdySign(ge.nodes, ss);

  ue[0] = ss * bPts[0];
}

void bdfmCube0Signs(const int nodes[8], int signs[6])
{
  int face[4];

  for(int i = 0; i < 6; i++)
    {
      for(int j = 0; j < 4; j++) face[j] = nodes[Cube.subs(2,i,j)];

      bdfmCubeBdySign(face, signs[i]);
    }

  return;
}

void bdfmCube0(svector<3> xi, svector<3> phi[6], double *dphi = NULL)
{
  double eta = xi[0], ata = xi[1], chi = xi[2];

  //  (-1,-1,-1), (1,-1,-1)   (-1, 1,-1), (1, 1,-1)   vertices
  //  (-1,-1, 1), (1,-1, 1)   (-1, 1, 1), (1, 1, 1)

  // faces  {{0,2,3,1}, {4,5,7,6}, {0,1,5,4}, {2,6,7,3}, {0,4,6,2}, {1,3,7,5}}

  phi[0] = svector<3>(0, 0, (chi-1)/2);
  phi[1] = svector<3>(0, 0, (chi+1)/2);
  phi[2] = svector<3>(0, (ata-1)/2, 0);
  phi[3] = svector<3>(0, (ata+1)/2, 0);
  phi[4] = svector<3>((eta-1)/2, 0, 0);
  phi[5] = svector<3>((eta+1)/2, 0, 0);
 
  if(dphi == NULL) return;

  dphi[0] = 0.5;
  dphi[1] = 0.5;
  dphi[2] = 0.5;
  dphi[3] = 0.5;
  dphi[4] = 0.5;
  dphi[5] = 0.5;
  
  return;
}

void bdfmCube0Interpolate(GeometricElement<3,8> &ge,
			    vector<svector<3> > &pPts, double ue[4])
{
  int i, signs[6];
  
  bdfmCube0Signs(ge.nodes, signs);  

  double det;
  smatrix<3> dxds, dsdxT;

  svector<3> nhat[6], nn;
  
  nhat[0] = svector<3>( 0.0, 0.0,-1.0);
  nhat[1] = svector<3>( 0.0, 0.0, 1.0);
  nhat[2] = svector<3>( 0.0,-1.0, 0.0);
  nhat[3] = svector<3>( 0.0, 1.0, 0.0);
  nhat[4] = svector<3>(-1.0, 0.0, 0.0);
  nhat[5] = svector<3>( 1.0, 0.0, 0.0);
  
  for(i = 0; i < 6; i++) 
    {
      ge.map(BDFMcube0Pts[i], dxds);
  
      dsdxT = (dxds.inverse(det)).transpose();

      nn = dsdxT * nhat[i]; nn = (1.0/nn.norm()) * nn;

      ue[i] = signs[i] * pPts[i].dot(nn);
    }

  return;
}

void bdfmCube0Transform(GeometricElement<3,8> &ge, smatrix<3> dxds,
			   svector<3> Phi[6], double Dphi[6],
			   svector<3> phi[6], double dphi[6])
{
  int i, signs[6];
  
  bdfmCube0Signs(ge.nodes, signs);  

  svector<3> nhat[6], nn;
  
  nhat[0] = svector<3>( 0.0, 0.0,-1.0);
  nhat[1] = svector<3>( 0.0, 0.0, 1.0);
  nhat[2] = svector<3>( 0.0,-1.0, 0.0);
  nhat[3] = svector<3>( 0.0, 1.0, 0.0);
  nhat[4] = svector<3>(-1.0, 0.0, 0.0);
  nhat[5] = svector<3>( 1.0, 0.0, 0.0);
  
  double det = fabs(dxds.determinant());

  smatrix<3> dxdsBdy, dxdsBdyCof;

  for(i = 0; i < 6; i++)
    {
      ge.map(BDFMcube0Pts[i], dxdsBdy);
  
      dxdsBdyCof = dxdsBdy.cofactor();

      double modFn = (signs[i]/det) * (dxdsBdyCof * nhat[i]).norm();

      phi[i]  = modFn * (dxds * Phi[i]);
      dphi[i] = modFn * Dphi[i];
    }

  return;
}
