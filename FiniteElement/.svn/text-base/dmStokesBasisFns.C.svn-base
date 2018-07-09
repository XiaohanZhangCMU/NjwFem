// Basis function for dual mixed formulation of Stokes problem
//
// 1) Construct a matrix valued H(div) baisis function for the stress
// 
// 2) Construct a matrix valued element for the skew part of the gradient
//    and provide Basis functions for the symmetric part (both discontinuous)

#include<vector>

void sortPermutation(int n, int *nn, int *pp);

template<const int Ndim, const int Nchi, const int Nchib>
struct DMstokesBasisFunctionS : public basisFunction<Ndim>
{
  DMstokesBasisFunctionS(basisFunction<Ndim> &bf, 
     void (*bb) (svector<Ndim>, smatrix<Ndim> *phi, svector<Ndim> *dphi),
     void (*tt)(int *nodes, smatrix<Ndim> &dxds,
		smatrix<Ndim> *Phi, svector<Ndim> *Dphi,
		smatrix<Ndim> *phi, svector<Ndim> *dphi),
     void (*ii)(GeometricElement<Ndim,Nchi> &ge,
		vector<smatrix<Ndim> > &sPts, double *ue),
     void (*bdybb)(svector<Ndim-1>, svector<Ndim> *phi),
     void (*bdyss)(const int *nodes, int &sign),
     void (*bdyii)(BoundaryGelement<Ndim,Nchib> &gbe,
		    vector<svector<Ndim> > &sPts, double *ue))
    : basisFunction<Ndim>(bf), basis(bb), transform(tt), interpolate(ii),
      bdyBasis(bdybb), bdySign(bdyss), bdyInterpolate(bdyii)
  {return;}

  void (*basis)(svector<Ndim>, smatrix<Ndim> *phi, svector<Ndim> *dphi);
  void (*transform)(int *nodes, smatrix<Ndim> &dxds,
		    smatrix<Ndim> *Phi, svector<Ndim> *Dphi,
		    smatrix<Ndim> *phi, svector<Ndim> *dphi);
  void (*interpolate)(GeometricElement<Ndim,Nchi> &ge,
		     vector<smatrix<Ndim> > &sPts, double *ue);
  void (*bdyBasis)(svector<Ndim-1>, svector<Ndim> *phi);
  void (*bdySign)(const int *nodes, int &sign);
  void (*bdyInterpolate)(BoundaryGelement<Ndim,Nchib> &gbe,
			 vector<svector<Ndim> > &sPts, double *ue);
};

template<const int Ndim>
struct DMstokesBasisFunctionG : public basisFunction<Ndim>
{
  DMstokesBasisFunctionG(basisFunction<Ndim> &bf, 
     void (*bb)(svector<Ndim>, smatrix<Ndim> *phi),
     void (*ii)(vector<smatrix<Ndim> > &gPts, double *gskw, double *gsym),
     void (*bbsym)(svector<Ndim>, smatrix<Ndim> *phi),
     int (*nsym)())
    : basisFunction<Ndim>(bf), basis(bb), interpolate(ii),
      basisGsym(bbsym), nGsym(nsym)
  {return;}

  void (*basis    )(svector<Ndim>, smatrix<Ndim> *phi);
  void (*interpolate)(vector<smatrix<Ndim> > &gPts, 
		      double *gskw, double *gsym);
  void (*basisGsym)(svector<Ndim>, smatrix<Ndim> *phi);
  int  (*nGsym)();
};

// AFW degree one triangle
// AFW degree one triangle
// AFW degree one triangle

void afwTriangle1BasisS(svector<2>, smatrix<2> phi[12], svector<2> dphi[12]);
void afwTriangle1PermutationS(const int nodes[6], int dof[12]);
void afwTriangle1TransformS(int *nodes, smatrix<2> &dxds,
			   smatrix<2> Phi[12], svector<2> Dphi[12],
			   smatrix<2> phi[12], svector<2> dphi[12]);
void afwTriangle1classInterpolateS(GeometricElement<2,3> &ge,
				   vector<smatrix<2> > &sPts, double ue[12]);

void afwTriangle1BdyBasisS(svector<1>, svector<2> phi[4]);
void afwTriangle1BdyPermutationS(const int *nodes, int *dof);
void afwTriangle1classBdyInterpolateS(BoundaryGelement<2,2> &gbe,
				      vector<svector<2> > &sPts, double ue[4]);

void afwTriangle1BasisGskw(svector<2>, smatrix<2> phi[1]);
void afwTriangle1classInterpolateG(vector<smatrix<2> > &gPts,
				   double gskw[1], double gsym[6]);
void afwTriangle1BasisGsym(svector<2>, smatrix<2> phi[6]);
int  afwTriangle1nGsym() { return(6); }

basisFunction<1> AFWtriangle1BdyBFS =    // AFW_1 Stress boundary
  {4,
   {0,4},
   EmptyVertex,
   &Interval,
   2,
   RTtriangle1BdyPts,
   afwTriangle1BdyPermutationS
  };

basisFunction<2> AFWtriangle1BFS =        // AFW_1 Stress (2 copies of BDM_1)
  {12,
   {0,4,0},
   &AFWtriangle1BdyBFS,
   &Triangle,
   6,
   RTtriangle1Pts,                    
   afwTriangle1PermutationS
  };

DMstokesBasisFunctionS<2,3,2> AFWtriangle1classS(AFWtriangle1BFS, 
   &afwTriangle1BasisS, &afwTriangle1TransformS, 
   &afwTriangle1classInterpolateS,
   &afwTriangle1BdyBasisS, &rtTriangleBdySign, 
   &afwTriangle1classBdyInterpolateS),
  *AFWtriangle1S = &AFWtriangle1classS;

basisFunction<2> AFWtriangle1BFG =        // AFW_1 gradient 
  {1,                                     // skew part p/w constant
   {0,0,1},
   EmptyInterval,
   &Triangle,
   4,
   Triangle3BubblePts,                    // last point is centroid 
   identityPermutation  
  };

DMstokesBasisFunctionG<2> AFWtriangle1classG(AFWtriangle1BFG, 
   &afwTriangle1BasisGskw, afwTriangle1classInterpolateG, 
   &afwTriangle1BasisGsym, &afwTriangle1nGsym),
  *AFWtriangle1G = &AFWtriangle1classG;


// ***************** Basis Functions ***************************
// ***************** Basis Functions ***************************
// ***************** Basis Functions ***************************

// AFW degree one triangle
// AFW degree one triangle
// AFW degree one triangle

void afwTriangle1BasisS(svector<2> xi, smatrix<2> phi[12], svector<2> *dphi)
{
  double dbdm[6];
  svector<2> bdm[6], zero(0,0);

  bdmTriangle1(xi, bdm, dbdm);

  for(int i = 0; i < 6; i++)   // rows of afw are bdm's
    {
      phi[2*i  ][0] = bdm[i]; phi[2*i  ][1] = zero;
      phi[2*i+1][0] = zero;   phi[2*i+1][1] = bdm[i];

      // phi[2*i  ] = smatrix<2>(bdm[i][0], bdm[i][1], 0,0);
      // phi[2*i+1] = smatrix<2>(0,0, bdm[i][0], bdm[i][1]);
    }

  if(dphi == NULL) return;

  for(int i = 0; i < 6; i++)   // rows of afw are bdm's
    {
      dphi[2*i  ] = svector<2>(dbdm[i], 0.0);
      dphi[2*i+1] = svector<2>(0.0, dbdm[i]);
    }

  return;
}

void afwTriangle1PermutationS(const int nodes[6], int dof[12])
{
  int temp;

  if(nodes[0] > nodes[1]) 
    { 
      temp = dof[0]; dof[0] = dof[2]; dof[2] = temp;
      temp = dof[1]; dof[1] = dof[3]; dof[3] = temp;
    }
  if(nodes[1] > nodes[2]) 
    {
      temp = dof[4]; dof[4] = dof[6]; dof[6] = temp;
      temp = dof[5]; dof[5] = dof[7]; dof[7] = temp;
    }
  if(nodes[2] > nodes[0]) 
    {
      temp = dof[8]; dof[8] = dof[10]; dof[10] = temp;
      temp = dof[9]; dof[9] = dof[11]; dof[11] = temp;
    }

  return;
}

void afwTriangle1TransformS(int *nodes, smatrix<2> &dxds,
			    smatrix<2> Phi[12], svector<2> Dphi[12],
			    smatrix<2> phi[12], svector<2> dphi[12])
{
  int i, signs[6];

  bdmTriangle1Signs(nodes, signs);  
 
  double det;
  smatrix<2> dxdsT = dxds.transpose();
  smatrix<2> dsdxT = dxdsT.inverse(det);

  svector<2> nn[3];
  
  nn[0] = svector<2>(0.0,-1.0); 
  nn[1] = svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0)); 
  nn[2] = svector<2>(-1.0,0.0); 
  
  if(det < 0) std::cout << "afwTriangle1Transform(): Warning det < 0" << "\n";

  for(i = 0; i < 6; i++)    // Piola transform
     {
     double modFn = signs[i] * (dsdxT * nn[i/2]).norm();

     phi[2*i  ]  = modFn *  Phi[2*i  ] * dxdsT;  // functions in H(div)
     phi[2*i+1]  = modFn *  Phi[2*i+1] * dxdsT;  // functions in H(div)

     dphi[2*i  ] = modFn * Dphi[2*i  ];          // divergences
     dphi[2*i+1] = modFn * Dphi[2*i+1];          // divergences
     }

  return;
}

void afwTriangle1classInterpolateS(GeometricElement<2,3> &ge,
				   vector<smatrix<2> > &sPts, double ue[12])
{
  int i, signs[6];
  bdmTriangle1Signs(ge.nodes, signs);  

  smatrix<2> dxds;
  ge.map(svector<2>(1.0/3.0,1.0/3.0), dxds);   // Jacobian is constant

  double det;  
  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<2> nn[3];
    
  nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1/nn[2].norm()) * nn[2];
  
  for(i = 0; i < 6; i++) 
     {
       svector<2> sn = signs[i] * ( sPts[i] * nn[i/2] );

       ue[2*i  ] = sn[0];
       ue[2*i+1] = sn[1];
     }

  return;
}

void afwTriangle1BdyBasisS(svector<1> xi, svector<2> phi[4])
{
  double psi[2];

  rtTriangle1Bdy(xi, psi);

  phi[0] = svector<2>(psi[0], 0);
  phi[1] = svector<2>(0, psi[0]);
  phi[2] = svector<2>(psi[1], 0);
  phi[3] = svector<2>(0, psi[1]);

  return;
}

void afwTriangle1BdyPermutationS(const int *nodes, int *dof)
{
  if(nodes[0] > nodes[1]) 
    {
      int temp;

      temp = dof[2]; dof[2] = dof[0]; dof[0] = temp;
      temp = dof[3]; dof[3] = dof[1]; dof[1] = temp;
   }

  return;
}

void afwTriangle1classBdyInterpolateS(BoundaryGelement<2,2> &gbe,
				       vector<svector<2> > &sPts, double ue[4])
{
  int sign;
  rtTriangleBdySign(gbe.nodes, sign);

  for(int i = 0; i < 2; i++)
    {
      ue[2*i  ] = sign * sPts[i][0];
      ue[2*i+1] = sign * sPts[i][1];
    }

  return;
}

void afwTriangle1BasisGskw(svector<2> xi, smatrix<2> phi[1])
{
  // 2d constant skew matrix

  phi[0] = smatrix<2>(0,1, -1,0);

  return;
}

void afwTriangle1BasisGsym(svector<2> xi, smatrix<2> phi[6])
{
  double psi[3];
  
  // symmetric trace free matrices p/w linear matrix

  tri3(xi, psi, NULL);

  for(int i = 0; i < 3; i++)   // symmetrix part of the matrix
    {
      phi[  i] = smatrix<2>(psi[i],0, 0,-psi[i]);   // diagonal
      phi[3+i] = smatrix<2>(0,psi[i], psi[i], 0);   // off diagonal
    }

  return;
}

void afwTriangle1classInterpolateG(vector<smatrix<2> > &gPts,
				   double gskw[1], double gsym[6])
{
  for(int i = 0; i < 3; i++)
    {
      gsym[i  ] = 0.5 * (gPts[i][0][0] - gPts[i][1][1]);
      gsym[i+3] = 0.5 * (gPts[i][0][1] + gPts[i][1][0]);
    }

  gskw[0] = 0.5*(gPts[3][0][1] - gPts[3][1][0]);    // value at centroid

  return;
}

// AFW degree one tet
// AFW degree one tet
// AFW degree one tet

void afwTet1BasisS(svector<3>, smatrix<3> phi[36], svector<3> dphi[36]);
void afwTet1PermutationS(const int nodes[4], int dof[36]);
void afwTet1TransformS(int *nodes, smatrix<3> &dxds,
		       smatrix<3> Phi[36], svector<3> Dphi[36],
		       smatrix<3> phi[36], svector<3> dphi[36]);
void afwTet1classInterpolateS(GeometricElement<3,4> &ge,
			      vector<smatrix<3> > &sPts, double ue[36]);

void afwTet1BdyBasisS(svector<2>, svector<3> phi[9]);
void afwTet1BdyPermutationS(const int *nodes, int *dof);
void afwTet1classBdyInterpolateS(BoundaryGelement<3,3> &gbe,
				 vector<svector<3> > &sPts, double ue[9]);

void afwTet1BasisGskw(svector<3>, smatrix<3> phi[1]);
void afwTet1classInterpolateG(vector<smatrix<3> > &gPts,
				   double gskw[1], double gsym[20]);
void afwTet1BasisGsym(svector<3>, smatrix<3> phi[20]);
int  afwTet1nGsym() { return(20); }

basisFunction<2> AFWtet1BdyBFS =    // AFW_1 Stress boundary
  {9,
   {0,0,9},
   EmptyInterval,
   &Triangle,
   3,
   BDMtet1BdyPts,
   afwTet1BdyPermutationS
  };

basisFunction<3> AFWtet1BFS =        // AFW_1 Stress (3 copies of BDM_1)
  {36,
   {0,0,9,0},
   &AFWtet1BdyBFS,
   &Tetrahedron,
   12,
   BDMtet1Pts,                    
   afwTet1PermutationS
  };

DMstokesBasisFunctionS<3,4,3> AFWtet1classS(AFWtet1BFS, 
   &afwTet1BasisS, &afwTet1TransformS, 
   &afwTet1classInterpolateS,
   &afwTet1BdyBasisS, &rtTetBdySign, 
   &afwTet1classBdyInterpolateS),
  *AFWtet1S = &AFWtet1classS;

basisFunction<3> AFWtet1BFG =          // AFW_1 gradient 
  {3,                                  // skew part p/w constant
   {0,0,0,3},
   EmptyTriangle,
   &Tetrahedron,
   5,
   Tet4BubblePts,                    // last point is centroid 
   identityPermutation  
  };

DMstokesBasisFunctionG<3> AFWtet1classG(AFWtet1BFG, 
   &afwTet1BasisGskw, afwTet1classInterpolateG, 
   &afwTet1BasisGsym, &afwTet1nGsym),
  *AFWtet1G = &AFWtet1classG;


// ***************** Basis Functions ***************************
// ***************** Basis Functions ***************************
// ***************** Basis Functions ***************************

// AFW degree one tet
// AFW degree one tet
// AFW degree one tet

void afwTet1BasisS(svector<3> xi, smatrix<3> phi[36], svector<3> *dphi)
{
  double dbdm[12];
  svector<3> bdm[12], zero(0,0,0);

  bdmTet1(xi, bdm, dbdm);

  for(int i = 0; i < 12; i++)   // rows of afw are bdm's
    {
    phi[3*i  ][0] = bdm[i]; phi[3*i  ][1] = zero;   phi[3*i  ][2] = zero;
    phi[3*i+1][0] = zero;   phi[3*i+1][1] = bdm[i]; phi[3*i+1][2] = zero;
    phi[3*i+2][0] = zero;   phi[3*i+2][1] = zero;   phi[3*i+2][2] = bdm[i];
    }

  if(dphi == NULL) return;

  for(int i = 0; i < 12; i++)   // rows of afw are bdm's
    {
      dphi[3*i  ] = svector<3>(dbdm[i], 0.0, 0.0);
      dphi[3*i+1] = svector<3>(0.0, dbdm[i], 0.0);
      dphi[3*i+2] = svector<3>(0.0, 0.0, dbdm[i]);
    }

  return;
}

void afwTet1PermutationS(const int nodes[4], int dof[36])
{
  int nn[3];     // sort dof face-by-face

  for(int i = 0; i < 4; i++) 
    {
      for(int j = 0; j < 3; j++) nn[j] = nodes[Tetrahedron.subs(2,i,j)];

      afwTet1BdyPermutationS(nn, dof+9*i);
    }

  return;
}

void afwTet1TransformS(int *nodes, smatrix<3> &dxds,
			    smatrix<3> Phi[36], svector<3> Dphi[36],
			    smatrix<3> phi[36], svector<3> dphi[36])
{
  int i, signs[12];

  bdmTet1Signs(nodes, signs);  
 
  double det;
  smatrix<3> dxdsT = dxds.transpose();
  smatrix<3> dsdxT = dxdsT.inverse(det);

  svector<3> nn[4];
  
  nn[0] = svector<3>( 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0));
  nn[1] = svector<3>(-1.0, 0.0, 0.0);
  nn[2] = svector<3>( 0.0,-1.0, 0.0);
  nn[3] = svector<3>( 0.0, 0.0,-1.0);
  
  if(det < 0) std::cout << "afwTet1Transform(): Warning det < 0" << "\n";

  for(i = 0; i < 12; i++)    // Piola transform
     {
     double modFn = signs[i] * (dsdxT * nn[i/3]).norm();

     phi[3*i  ]  = modFn *  Phi[3*i  ] * dxdsT;  // functions in H(div) 
     phi[3*i+1]  = modFn *  Phi[3*i+1] * dxdsT;  // functions in H(div)
     phi[3*i+2]  = modFn *  Phi[3*i+2] * dxdsT;  // functions in H(div)

     dphi[3*i  ] = modFn * Dphi[3*i  ];          // divergences
     dphi[3*i+1] = modFn * Dphi[3*i+1];          // divergences
     dphi[3*i+2] = modFn * Dphi[3*i+2];          // divergences
     }

  return;
}

void afwTet1classInterpolateS(GeometricElement<3,4> &ge,
			      vector<smatrix<3> > &sPts, double ue[36])
{
  int i, signs[12];
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
  
  for(i = 0; i < 12; i++) 
     {
       svector<3> sn = signs[i] * ( sPts[i] * nn[i/3] );

       ue[3*i  ] = sn[0];
       ue[3*i+1] = sn[1];
       ue[3*i+2] = sn[2];
     }
  
  return;
}

void afwTet1BdyBasisS(svector<2> xi, svector<3> phi[9])
{
  double psi[3];

  bdmTet1Bdy(xi, psi);

  phi[0] = svector<3>(psi[0], 0, 0);
  phi[1] = svector<3>(0, psi[0], 0);
  phi[2] = svector<3>(0, 0, psi[0]);
  phi[3] = svector<3>(psi[1], 0, 0);
  phi[4] = svector<3>(0, psi[1], 0);
  phi[5] = svector<3>(0, 0, psi[1]);
  phi[6] = svector<3>(psi[2], 0, 0);
  phi[7] = svector<3>(0, psi[2], 0);
  phi[8] = svector<3>(0, 0, psi[2]);

  return;
}

void afwTet1BdyPermutationS(const int nodes[3], int dof[9])
{
  int i, alpha, nn[3], pp[3], dof0[9];

  for(i = 0; i < 3; i++)
    {
      nn[i] = nodes[i];

      for(alpha = 0; alpha < 3; alpha++) dof0[3*i+alpha] = dof[3*i+alpha];
    }

  sortPermutation(3, nn, pp);

  for(i = 0; i < 3; i++)
    {
      for(alpha = 0; alpha < 3; alpha++) dof[3*pp[i]+alpha] = dof[3*i+alpha];
    }

  return;
}

void afwTet1classBdyInterpolateS(BoundaryGelement<3,3> &gbe,
				 vector<svector<3> > &sPts, double ue[9])
{
  int sign;
  rtTetBdySign(gbe.nodes, sign);

  for(int i = 0; i < 3; i++)
    {
      ue[3*i  ] = sign * sPts[i][0];
      ue[3*i+1] = sign * sPts[i][1];
      ue[3*i+2] = sign * sPts[i][2];
    }

  return;
}

void afwTet1BasisGskw(svector<3> xi, smatrix<3> phi[3])
{
  // 3d constant skew matrix; if a = (a0,a1,a2) then W(a) n = a x n

  phi[0] = smatrix<3>(0, 0,0,  0,0,-1,  0,1,0);
  phi[1] = smatrix<3>(0, 0,1,  0,0, 0, -1,0,0);
  phi[2] = smatrix<3>(0,-1,0,  1,0, 0,  0,0,0);

  return;
}

void afwTet1BasisGsym(svector<3> xi, smatrix<3> phi[20])
{
  double psi[4];

  svector<3> d1(-1.0/sqrt(3.0), (0.5/sqrt(3.0)+0.5), (0.5/sqrt(3.0)-0.5));
  svector<3> d2(-1.0/sqrt(3.0), (0.5/sqrt(3.0)-0.5), (0.5/sqrt(3.0)+0.5));
  
  // symmetric trace free matrices p/w linear matrix

  // G[i][j] = i(i+1)/2 + j - 1  for i=1..(d-1), j=0..i

  // Basis for the diagonal is:
  // (-1/sqrt3, (1/sqrt3+1)/2, (1/sqrt3-1)/2
  // (-1/sqrt3, (1/sqrt3-1)/2, (1/sqrt3+1)/2

  tet4(xi, psi, NULL);

  for(int i = 0; i < 4; i++)   // symmetrix part of the matrix
    {
      svector<3> d1i = psi[i] * d1, d2i = psi[i] * d2;

      phi[5*i  ] = smatrix<3>(0,psi[i],0, psi[i],0,0, 0,0,0     );   // (1,0)
      phi[5*i+1] = smatrix<3>(d1i[0],0,0, 0,d1i[1],0, 0,0,d1i[2]);   // (1,1)
      phi[5*i+2] = smatrix<3>(0,0,psi[i], 0,0,0,      psi[i],0,0);   // (2,0)
      phi[5*i+3] = smatrix<3>(0,0,0,      0,0,psi[i], 0,psi[i],0);   // (2,1)
      phi[5*i+4] = smatrix<3>(d2i[0],0,0, 0,d2i[1],0, 0,0,d2i[2]);   // (2,2)
    }

  return;
}

void afwTet1classInterpolateG(vector<smatrix<3> > &gPts,
			      double gskw[3], double gsym[20])
{
  svector<3> d1(-1.0/sqrt(3.0), (0.5/sqrt(3.0)+0.5), (0.5/sqrt(3.0)-0.5));
  svector<3> d2(-1.0/sqrt(3.0), (0.5/sqrt(3.0)-0.5), (0.5/sqrt(3.0)+0.5));

  for(int i = 0; i < 4; i++)
    {
      svector<3> dd(gPts[i][0][0], gPts[i][1][1], gPts[i][2][2]);

      gsym[5*i  ] = 0.5 * (gPts[i][0][1] + gPts[i][1][0]);   // (1,0)
      gsym[5*i+1] = dd.dot(d1);                              // (1,1)
      gsym[5*i+2] = 0.5 * (gPts[i][0][2] + gPts[i][2][0]);   // (2,0)
      gsym[5*i+3] = 0.5 * (gPts[i][1][2] + gPts[i][2][1]);   // (2,1)
      gsym[5*i+4] = dd.dot(d2);                              // (2,2)
    }

  gskw[0] = 0.5*(gPts[4][2][1] - gPts[4][1][2]);    // value at centroid
  gskw[1] = 0.5*(gPts[4][0][2] - gPts[4][2][0]);    // value at centroid
  gskw[2] = 0.5*(gPts[4][1][0] - gPts[4][0][1]);    // value at centroid

  return;
}

// SV degree one triangle
// SV degree one triangle
// SV degree one triangle

void svTriangle1BasisS(svector<2>, smatrix<2> phi[12], svector<2> dphi[12]);
void svTriangle1TransformS(int *nodes, smatrix<2> &dxds,
			   smatrix<2> Phi[12], svector<2> Dphi[12],
			   smatrix<2> phi[12], svector<2> dphi[12]);
void svTriangle1classInterpolateS(GeometricElement<2,3> &ge,
				   vector<smatrix<2> > &sPts, double ue[12]);

void svTriangle1BasisGskw(svector<2>, smatrix<2> phi[1]);
void svTriangle1classInterpolateG(vector<smatrix<2> > &gPts,
				   double gskw[1], double gsym[6]);

basisFunction<2> SVtriangle1BFS =        // SV_1 Stress (2 copies of RT_1)
  {16,
   {0,4,4},
   &AFWtriangle1BdyBFS,                  // Boundary identical to AFW
   &Triangle,
   7,
   RTtriangle1Pts,                    
   afwTriangle1PermutationS
  };

DMstokesBasisFunctionS<2,3,2> SVtriangle1classS(SVtriangle1BFS, 
   &svTriangle1BasisS, &svTriangle1TransformS, 
   &svTriangle1classInterpolateS,
   &afwTriangle1BdyBasisS, &rtTriangleBdySign, 
   &afwTriangle1classBdyInterpolateS),
  *SVtriangle1S = &SVtriangle1classS;

basisFunction<2> SVtriangle1BFG =        // SV_1 gradient 
  {3,                                    // skew part p/w linear
   {0,0,3},
   EmptyInterval,
   &Triangle,
   3,
   Triangle6Pts,                          // first 3 are vertices
   discontinuousTriangle3Permutation  
  };

DMstokesBasisFunctionG<2> SVtriangle1classG(SVtriangle1BFG, 
   &svTriangle1BasisGskw, svTriangle1classInterpolateG, 
   &afwTriangle1BasisGsym, &afwTriangle1nGsym),
  *SVtriangle1G = &SVtriangle1classG;


// ***************** Basis Functions ***************************
// ***************** Basis Functions ***************************
// ***************** Basis Functions ***************************

// SV degree one triangle
// SV degree one triangle
// SV degree one triangle

void svTriangle1BasisS(svector<2> xi, smatrix<2> phi[16], svector<2> *dphi)
{
  double drt1[8];
  svector<2> rt1[8], zero(0,0);

  rtTriangle1(xi, rt1, drt1);

  for(int i = 0; i < 8; i++)   // rows of sv are rt1's
    {
      phi[2*i  ][0] = rt1[i]; phi[2*i  ][1] = zero;
      phi[2*i+1][0] = zero;   phi[2*i+1][1] = rt1[i];

      // phi[2*i  ] = smatrix<2>(rt1[i][0], rt1[i][1], 0,0);
      // phi[2*i+1] = smatrix<2>(0,0, rt1[i][0], rt1[i][1]);
    }

  if(dphi == NULL) return;

  for(int i = 0; i < 8; i++)   // rows of sv are rt1's
    {
      dphi[2*i  ] = svector<2>(drt1[i], 0.0);
      dphi[2*i+1] = svector<2>(0.0, drt1[i]);
    }

  return;
}

void svTriangle1TransformS(int *nodes, smatrix<2> &dxds,
			    smatrix<2> Phi[16], svector<2> Dphi[16],
			    smatrix<2> phi[16], svector<2> dphi[16])
{
  int i, signs[8];

  rtTriangle1Signs(nodes, signs);  
 
  double det;
  smatrix<2> dxdsT = dxds.transpose();
  smatrix<2> dsdxT = dxdsT.inverse(det);

  svector<2> nn[3];
  
  nn[0] = svector<2>(0.0,-1.0); 
  nn[1] = svector<2>(1.0/sqrt(2.0), 1.0/sqrt(2.0)); 
  nn[2] = svector<2>(-1.0,0.0); 
  
  if(det < 0) std::cout << "svTriangle1Transform(): Warning det < 0" << "\n";

  for(i = 0; i < 6; i++)    // Piola transform
     {
     double modFn = signs[i] * (dsdxT * nn[i/2]).norm();

     phi[2*i  ]  = modFn *  Phi[2*i  ] * dxdsT;  // functions in H(div)
     phi[2*i+1]  = modFn *  Phi[2*i+1] * dxdsT;  // functions in H(div)

     dphi[2*i  ] = modFn * Dphi[2*i  ];          // divergences
     dphi[2*i+1] = modFn * Dphi[2*i+1];          // divergences
     }

  smatrix<2> temp0, temp1;

  for(i = 0; i < 2; i++)              // dof 6,7 of RT1
    {
      temp0.zero();
      temp1.zero();

      dphi[12+2*i  ].zero();
      dphi[12+2*i+1].zero();

      for(int alpha = 0; alpha < 2; alpha++) 
	{
	  temp0 = temp0 + dsdxT[i][alpha] * Phi[12+2*alpha];
	  temp1 = temp1 + dsdxT[i][alpha] * Phi[12+2*alpha+1];

	  dphi[12+2*i  ] += dsdxT[i][alpha] * Dphi[12+2*alpha];
	  dphi[12+2*i+1] += dsdxT[i][alpha] * Dphi[12+2*alpha+1];
	}

      phi[12+2*i  ] = temp0 * dxdsT;
      phi[12+2*i+1] = temp1 * dxdsT;
    }

  return;
}

void svTriangle1classInterpolateS(GeometricElement<2,3> &ge,
				   vector<smatrix<2> > &sPts, double ue[16])
{
  int i, signs[8];
  rtTriangle1Signs(ge.nodes, signs);  

  smatrix<2> dxds;
  ge.map(svector<2>(1.0/3.0,1.0/3.0), dxds);   // Jacobian is constant

  double det;  
  smatrix<2> dsdxT = (dxds.inverse(det)).transpose();
  
  svector<2> nn[3];
    
  smatrix<2> temp;

  nn[0] = dsdxT * svector<2>(0.0,-1.0); nn[0] = (1/nn[0].norm()) * nn[0];
  nn[1] = dsdxT * svector<2>(1.0, 1.0); nn[1] = (1/nn[1].norm()) * nn[1];
  nn[2] = dsdxT * svector<2>(-1.0,0.0); nn[2] = (1/nn[2].norm()) * nn[2];
  
  temp.zero();

  for(i = 0; i < 6; i++) 
     {
       svector<2> sn = signs[i] * ( sPts[i] * nn[i/2] );

       ue[2*i  ] = sn[0];
       ue[2*i+1] = sn[1];

       temp = temp + sPts[i];
     }

  temp = (1.0/12.0) * temp + 0.5 * sPts[6];

  ue[12] = temp[0][0];
  ue[13] = temp[1][0];
  ue[14] = temp[0][1];
  ue[15] = temp[1][1];

  return;
}

void svTriangle1BasisGskw(svector<2> xi, smatrix<2> phi[3])
{
  double psi[3];

  // 2d p/w linear skew matrix

  tri3(xi, psi, NULL);

  phi[0] = smatrix<2>(0,psi[0], -psi[0],0);
  phi[1] = smatrix<2>(0,psi[1], -psi[1],0);
  phi[2] = smatrix<2>(0,psi[2], -psi[2],0);

  return;
}

void svTriangle1classInterpolateG(vector<smatrix<2> > &gPts,
				   double gskw[3], double gsym[6])
{
  for(int i = 0; i < 3; i++)
    {
      gskw[i]   = 0.5 * (gPts[i][0][1] - gPts[i][1][0]);

      gsym[i  ] = 0.5 * (gPts[i][0][0] - gPts[i][1][1]);
      gsym[i+3] = 0.5 * (gPts[i][0][1] + gPts[i][1][0]);
    }

  return;
}
