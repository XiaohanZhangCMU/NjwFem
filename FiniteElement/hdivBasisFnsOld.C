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
     void (*ii)(GeometricElement<Ndim, Nchi> &ge, 
		vector<double> &bPts, double *ue))
    : basisFunction<Ndim>(bf), basis(bb), dofSign(ss), interpolate(ii)  
  {return;}

  void (*basis)(svector<Ndim>, double *phi);
  void (*dofSign)(const int *nodes, int &sign);
  void (*interpolate)(GeometricElement<Ndim, Nchi> &ge, 
		      vector<double> &bPts, double *ue);
};

template<const int Ndim, const int Nchi, const int Nchib>
struct HdivBasisFunction : public basisFunction<Ndim>
{
  HdivBasisFunction(basisFunction<Ndim> &bf, 
       HdivBdyBasisFunction<Ndim-1, Nchib> *bdy,
       void (*bb)(svector<Ndim>, svector<Ndim> *phi, double *dphi),
       void (*ss)(const int *nodes, int *signs),
       void (*ii)(GeometricElement<Ndim, Nchi> &ge,
		  vector<svector<Ndim> > &pPts, double *ue))
    : basisFunction<Ndim>(bf), boundary(bdy), basis(bb), dofSigns(ss),
      interpolate(ii)
  {return;}

  HdivBdyBasisFunction<Ndim-1, Nchib> *boundary;
  void (*basis)(svector<Ndim>, svector<Ndim> *phi, double *dphi);
  void (*dofSigns)(const int *nodes, int *signs);
  void (*ii)(GeometricElement<Ndim, Nchi> &ge, 
	     vector<svector<Ndim> > &pPts, double *ue);
};

void rtTriangle0Bdy(svector<1>, double *phi);
void rtTriangleBdySign(const int *nodes, int &sign);

void rtTriangle0(svector<2> xx, svector<2> *phi, double *dphi);
void rtTriangle0Signs(const int *nodes, int *signs);

void rtTriangle1Bdy(svector<1>, double phi[2]);
void rtTriangle1BdyPermutation(const int *nodes, int *dof);

void rtTriangle1(svector<2> xx, svector<2> phi[8], double dphi[8]);
void rtTriangle1Signs(const int *nodes, int signs[8]);
void rtTriangle1Permutation(const int *nodes, int *dof);

void bdmTriangle1(svector<2> xx, svector<2> phi[6], double dphi[6]);
void bdmTriangle1Signs(const int *nodes, int signs[6]);

void bdfmSquare0(svector<2> xi, svector<2> phi[4], double *dphi);
void bdfmSquare0Signs(const int nodes[4], int signs[4]);

void rtTetBdySign(const int *nodes, int &sign);

void rtTet0Bdy(svector<2>, double *phi);
void rtTet0(svector<3> xx, svector<3> *phi, double *dphi);
void rtTet0Signs(const int *nodes, int *signs);

void rtTet1Bdy(svector<2>, double phi[3]);
void rtTet1(svector<3> xx, svector<3> phi[15], double dphi[15]);
void rtTet1Signs(const int nodes[4], int signs[15]);
void rtTet1Permutation(const int nodes[4], int dof[15]);

void bdfmCube0Bdy(svector<2> xi, double *phi);
void bdfmCubeBdySign(const int *nodes, int &sign);
void bdfmCube0(svector<3> xi, svector<3> phi[6], double *dphi);
void bdfmCube0Signs(const int nodes[8], int signs[6]);

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

HdivBdyBasisFunction<1> RTtriangle0BdyClass(
   RTtriangle0BdyBF, &rtTriangle0Bdy, &rtTriangleBdySign),
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

HdivBasisFunction<2> RTtriangle0class(RTtriangle0BF, 
   RTtriangle0Bdy, &rtTriangle0, &rtTriangle0Signs),
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

HdivBdyBasisFunction<1> RTtriangle1BdyClass(
   RTtriangle1BdyBF, &rtTriangle1Bdy, &rtTriangleBdySign),
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
   rtTriangle1BdyPermutation
  };

HdivBasisFunction<2> RTtriangle1class(RTtriangle1BF, 
   RTtriangle1Bdy, &rtTriangle1, &rtTriangle1Signs),
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

HdivBasisFunction<2> BDMtriangle1class(BDMtriangle1BF, 
   RTtriangle1Bdy, &bdmTriangle1, &bdmTriangle1Signs),
  *BDMtriangle1 = &BDMtriangle1class;

// ***************** BDFM type elements on square ***********************
// ***************** BDFM type elements on square ***********************
// ***************** BDFM type elements on square ***********************

// basisFunction<1> BDFMsquare0Bdy ?=? RTtriangleBdy0 (check basis fn)

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

HdivBasisFunction<2> BDFMsquare0class(BDFMsquare0BF, 
   RTtriangle0Bdy, &bdfmSquare0, &bdfmSquare0Signs),
  *BDFMsquare0 = &BDFMsquare0class;

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

HdivBdyBasisFunction<2> RTtet0BdyClass(RTtet0BdyBF, &rtTet0Bdy, &rtTetBdySign),
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

HdivBasisFunction<3> RTtet0class(RTtet0BF, RTtet0Bdy, &rtTet0, &rtTet0Signs),
  *RTtet0 = &RTtet0class;

// Construct the RT-1 degree one tetrahedra
// Construct the RT-1 degree one tetrahedra
// Construct the RT-1 degree one tetrahedra

svector<2> RTtet1BdyPts[] = 
  {
    svector<2>(1.0/6.0, 1.0/6.0), 
    svector<2>(2.0/3.0, 1.0/6.0),
    svector<2>(1.0/6.0, 2.0/3.0) 
  };

basisFunction<2> RTtet1BdyBF =    // RT_0 boundary
  {3,
   {0,0,3},
   EmptyInterval,
   &Triangle,
   3,
   RTtet1BdyPts,
   discontinuousTriangle3Permutation
  };

HdivBdyBasisFunction<2> RTtet1BdyClass(RTtet1BdyBF, &rtTet1Bdy, &rtTetBdySign),
  *RTtet1Bdy = &RTtet1BdyClass;

svector<3> RTtet1Pts[] = 
  {
    svector<3>(1.0/6, 1.0/6, 2.0/3), 
    svector<3>(2.0/3, 1.0/6, 1.0/6),
    svector<3>(1.0/6, 2.0/3, 1.0/6),
    svector<3>(0, 2.0/3, 1.0/6),
    svector<3>(0, 1.0/6, 1.0/6),
    svector<3>(0, 1.0/6, 2.0/3),
    svector<3>(1.0/6, 0, 2.0/3),
    svector<3>(1.0/6, 0, 1.0/6),
    svector<3>(2.0/3, 0, 1.0/6),
    svector<3>(2.0/3, 1.0/6, 0),
    svector<3>(1.0/6, 1.0/6, 0),
    svector<3>(1.0/6, 2.0/3, 0),
    svector<3>(1.0/4, 1.0/4, 1.0/4)
  };

basisFunction<3> RTtet1BF =    // RT_1 element
  {15,
   {0,0,3,3},
   &RTtet1BdyBF,
   &Tetrahedron,
   13,
   RTtet1Pts,
   rtTet1Permutation
  };

HdivBasisFunction<3> RTtet1class(RTtet1BF, RTtet1Bdy, &rtTet1, &rtTet1Signs),
  *RTtet1 = &RTtet1class;

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

HdivBdyBasisFunction<2> 
   BDFMcube0BdyClass(BDFMcube0BdyBF, &bdfmCube0Bdy, &bdfmCubeBdySign),
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

HdivBasisFunction<3> 
  BDFMcube0class(BDFMcube0BF, BDFMcube0Bdy, &bdfmCube0, &bdfmCube0Signs),
  *BDFMcube0 = &BDFMcube0class;

// ********************* Basis Functions ************************
// ********************* Basis Functions ************************
// ********************* Basis Functions ************************

void rtTriangle0Bdy(svector<1> xi, double *phi)
{
  phi[0] = 0.5;       // reciprical of length of master edge [-1,1]
}

void rtTriangleBdySign(const int *nodes, int &sign)
{
  // Convention: global edge orientation is low to high

  if(nodes[0] > nodes[1]) sign = -1; else sign = 1;
}

void rtTriangle0(svector<2> xi, svector<2> *phi, double *dphi = NULL)
{
  double eta = xi[0], ata = xi[1];

 phi[0] = svector<2>(eta, ata-1);
 phi[1] = svector<2>(eta, ata);
 phi[2] = svector<2>(eta-1, ata);

 if(dphi == NULL) return;

 dphi[0] = 2.0;
 dphi[1] = 2.0;
 dphi[2] = 2.0;
 
 return;
}

void rtTriangle0Signs(const int *nodes, int *signs)
{
  // Convention: global edge orientation is low to high

  // master edges  {{0,1},{1,2},{2,0}},  (better is edge opposite 0, 1, 2)

  if(nodes[0] > nodes[1]) signs[0] = -1; else signs[0] = 1;
  if(nodes[1] > nodes[2]) signs[1] = -1; else signs[1] = 1;
  if(nodes[2] > nodes[0]) signs[2] = -1; else signs[2] = 1;

  return;
}

void rtTriangle0InterpolateS(GeometricElement<2,3> &ge,
			      vector<svector<2> > &pPts, double ue[3])
{
  int i, signs[3];
  
  rtTriangle0Signs(ge.nodes, signs);  

  smatrix<2> dxds;

  ge.map(svector<2>(1.0/3.0,1.0/3.0), dxds);   // Jacobian is constant
  
  smatrix<2> cof = dxds.cofactor();
  
  svector<2> nn[3];
  
  nn[0] = cof * svector<2>(0.0,-1.0);  // area * normal
  nn[1] = cof * svector<2>(1.0, 1.0);  
  nn[2] = cof * svector<2>(-1.0,0.0);
  
  for(i = 0; i < 3; i++) ue[i] = signs[i] * pPts[i].dot(nn[i]);

  return;
}

