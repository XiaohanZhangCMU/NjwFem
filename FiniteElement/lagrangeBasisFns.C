// Classical Lagrange basis functions are scalar valued and all
// derivatives are defined.  The boundary is a Lagrange basis
// function of one dimension less.

#include "../LinearAlgebra/svector.C"
#include "../LinearAlgebra/smatrix.C"
#include "../LinearAlgebra/sbmatrix.C"

void sortDof(int n, int *nn, int *dd);  // in utility.C

void line1(svector<1> peta, double psi[1], svector<1> dpsi[1])
{
  psi[0] = 1.0;  
  
  if(dpsi == NULL) return;

  dpsi[0][0] = 0.0; 

  return;
}

void line2(svector<1> peta, double psi[2], svector<1> dpsi[2])
{
  double eta = peta[0];

  psi[0] = (1-eta)/2.0;  
  psi[1] = (1+eta)/2.0; 
  
  if(dpsi == NULL) return;

  dpsi[0][0] = -0.5; 
  dpsi[1][0] =  0.5;

  return;
}

void line2bar(svector<1> peta, double psi[2], 
	      svector<1> dpsi[2], double psiBar[2])
{
  double eta = peta[0];

  psi[0] = (1-eta)/2.0;  
  psi[1] = (1+eta)/2.0; 

  if(dpsi == NULL) return;

  dpsi[0][0] = -0.5; 
  dpsi[1][0] =  0.5;

  psiBar[0] = psiBar[1] = 0.5;  // projections onto degree zero polys

  return;
}

void line3(svector<1> peta, double psi[3], svector<1> dpsi[3])
{
 double eta = peta[0];
 double ope = 1+eta, ome = 1-eta;

 psi[0] = -eta * ome / 2; 
 psi[1] =  eta * ope / 2; 
 psi[2] =  ope * ome    ; 

  if(dpsi == NULL) return;

  dpsi[0][0] = -0.5 + eta;
  dpsi[1][0] =  0.5 + eta;
  dpsi[2][0] = -2.0 * eta;

return;
}

void line3bar(svector<1> peta, double psi[3], 
	      svector<1> dpsi[3], double psiBar[3])
{
 double eta = peta[0];
 double ope = 1+eta, ome = 1-eta;

 psi[0] = -eta * ome / 2; 
 psi[1] =  eta * ope / 2; 
 psi[2] =  ope * ome    ; 

  if(dpsi == NULL) return;

  dpsi[0][0] = -0.5 + eta;
  dpsi[1][0] =  0.5 + eta;
  dpsi[2][0] = -2.0 * eta;

  psiBar[0] = 1.0/6.0 - 0.5*eta;     // projections onto degree one polys
  psiBar[1] = 1.0/6.0 + 0.5*eta;
  psiBar[2] = 2.0/3.0;

  return;
}

void line4(svector<1> pxi, double psi[4], svector<1> dpsi[4])
{
  double xi = pxi[0], xisq = xi*xi;
  double sqrt5 = sqrt(5.0), w1 = 5.0/8.0, w0 = sqrt5*w1;

// Use the Gauss Labotto Points -1, 1, -1/sqrt(5), 1/sqrt(5)

  psi[0] = w1 * (xisq - 1.0/5.0) * (1.0 - xi);
  psi[1] = w1 * (xisq - 1.0/5.0) * (1.0 + xi);
  psi[2] = w0 * (1.0 - xisq) * (1.0/sqrt5 - xi);
  psi[3] = w0 * (1.0 - xisq) * (1.0/sqrt5 + xi);
  
  if(dpsi == NULL) return;
  
  dpsi[0][0] = w1*(-3*xisq + 2*xi + 1.0/5.0);
  dpsi[1][0] = w1*( 3*xisq + 2*xi - 1.0/5.0);
  dpsi[2][0] = w0*( 3*xisq - (2/sqrt5)*xi - 1.0);
  dpsi[3][0] = w0*(-3*xisq - (2/sqrt5)*xi + 1.0);

// Old code used equispaced points, -1, 1, -1/3, 1/3

//   psi[0] = (9.0/16.0) * (xisq - 1.0/9.0) * (1.0 - xi);
//   psi[1] = (9.0/16.0) * (xisq - 1.0/9.0) * (1.0 + xi);
//   psi[2] = (27.0/16.0) * (1.0 - xisq) * (1.0/3.0 - xi);
//   psi[3] = (27.0/16.0) * (1.0 - xisq) * (1.0/3.0 + xi);

//  if(dpsi == NULL) return;
 
//   dpsi[0] =-(27.0/16.0)*xisq + (9.0/8.0)*xi + (1.0/16.0);
//   dpsi[1] = (27.0/16.0)*xisq + (9.0/8.0)*xi - (1.0/16.0);
//   dpsi[2] = (81.0/16.0)*xisq - (9.0/8.0)*xi - (27.0/16.0);
//   dpsi[3] =-(81.0/16.0)*xisq - (9.0/8.0)*xi + (27.0/16.0);

return;
}

void line4Permutation(const int *nodes, int *dof)
{
  if(nodes[0] > nodes[1])
    {
      int temp = dof[2]; dof[2] = dof[3]; dof[3] = temp;
    }

  return;
}

void line5(svector<1> pxi, double psi[5], svector<1> dpsi[5])
{
  double eta = pxi[0], gp =sqrt(21.0);
  double eta2 = eta*eta, eta3 = eta2*eta;
  double em1 = eta-1, ep1 = eta+1, esq7m3 = 7*eta2-3;

// Use the Gauss Labotto Points -1, 1, -sqrt(3/7), 0, sqrt(3/7)

  psi[0] = 0.125*eta*em1*esq7m3;
  psi[1] = 0.125*eta*ep1*esq7m3;
  psi[2] = (-7.0/24.0)*eta*em1*ep1*(7*eta-gp);
  psi[3] = ( 1.0/ 3.0)*em1*ep1*esq7m3;
  psi[4] = (-7.0/24.0)*eta*em1*ep1*(7*eta+gp);
  
  
  if(dpsi == NULL) return;
  
  dpsi[0] = svector<1>(3.5*eta3-0.75*eta-(21.0/8.0)*eta2+0.375);
  dpsi[1] = svector<1>(3.5*eta3+(21.0/8.0)*eta2-0.75*eta-0.375);
  dpsi[2] = svector<1>((-7.0/24.0)*(4*eta2-gp*eta+1)*(7*eta+gp));
  dpsi[3] = svector<1>(( 4.0/ 3.0)*eta*(-5+7*eta2));
  dpsi[4] = svector<1>((-7.0/24.0)*(4*eta2+gp*eta+1)*(7*eta-gp));
    
  return;
}

void line5Permutation(const int *nodes, int *dof)
{
  if(nodes[0] > nodes[1])
    {
      int temp = dof[2]; dof[2] = dof[4]; dof[4] = temp;
    }

  return;
}

void line6(svector<1> pxi, double psi[6], svector<1> dpsi[6])
{
  double s = pxi[0], sqrt30 = sqrt(30.0);
  double s2 = s*s, s3 = s2*s, s4 = s2*s2;
  double sm1 = s-1, sp1 = s+1;
  double pt1 = sqrt(525.0 + 70.0*sqrt30), pt2 = sqrt(525.0 - 70.0*sqrt30);

// Use the Gauss Labotto Points -1, 1, -pt1/35, -pt2/35, pt2/35, pt1/35
  
psi[0] = -(1.0/16)*sm1*(35*s4-30*s2+3);
psi[1] = (1.0/16)*sp1*(35*s4-30*s2+3);
psi[2] = -(1.0/117600)*pt1*(3*sqrt30-5)*(35*s2-15+2*sqrt30)*(-35*s+pt1)*sm1*sp1;
psi[3] = -(1.0/117600)*pt2*(5+3*sqrt30)*(-35*s2+15+2*sqrt30)*(-35*s+pt2)*sm1*sp1;
psi[4] = -(1.0/117600)*pt2*(5+3*sqrt30)*(-35*s2+15+2*sqrt30)*sm1*sp1*(35*s+pt2);
psi[5] = -(1.0/117600)*pt1*(3*sqrt30-5)*(35*s2-15+2*sqrt30)*(35*s+pt1)*sm1*sp1;

  if(dpsi == NULL) return;
  
dpsi[0] = svector<1>(-(175.0/16)*s4+(45.0/8)*s2-3.0/16+(35.0/4)*s3-(15.0/4)*s);
dpsi[1] = svector<1>((175.0/16)*s4-(45.0/8)*s2+3.0/16+(35.0/4)*s3-(15.0/4)*s);
dpsi[2] = svector<1>(-(1.0/117600)*pt1*(3*sqrt30-5)*(-525+70*sqrt30+140*s3*pt1+5250*s2-210*sqrt30*s2-6125*s4-100*pt1*s+4*pt1*sqrt30*s));
dpsi[3] = svector<1>(-(1.0/117600)*pt2*(5+3*sqrt30)*(70*sqrt30+525-140*s3*pt2-5250*s2-210*sqrt30*s2+6125*s4+100*s*pt2+4*s*sqrt30*pt2));
dpsi[4] = svector<1>(-(1.0/117600)*pt2*(5+3*sqrt30)*(-70*sqrt30-525-140*s3*pt2+5250*s2+210*sqrt30*s2-6125*s4+100*s*pt2+4*s*sqrt30*pt2));
dpsi[5] = svector<1> (-(1.0/117600)*pt1*(3*sqrt30-5)*(525-70*sqrt30+140*s3*pt1-5250*s2+210*s2*sqrt30+6125*s4-100*s*pt1+4*s*sqrt30*pt1));
  
  return;
}

void line6Permutation(const int *nodes, int *dof)
{
  if(nodes[0] > nodes[1])
    {
      int temp;
      temp = dof[2]; dof[2] = dof[5]; dof[5] = temp;
      temp = dof[3]; dof[3] = dof[4]; dof[4] = temp;
    }

  return;
}

void tri3(svector<2> xi, double sn[], svector<2> dsn[])
{
  double eta = xi[0], ata = xi[1], omema = 1.0-eta-ata;

  sn[0] = omema;  
  sn[1] = eta;    
  sn[2] = ata;    

 if(dsn == NULL) return;

 dsn[0] = svector<2>(-1,-1);
 dsn[1] = svector<2>(1,0);
 dsn[2] = svector<2>(0,1);
 
 return;
}

void tri3Bubble(svector<2> xi, double sn[4], svector<2> dsn[4])
{
  double eta = xi[0], ata = xi[1], omema = 1.0-eta-ata, bb = eta*ata;

  sn[0] = omema;  
  sn[1] = eta;    
  sn[2] = ata;    

  sn[3] = bb*omema;

 if(dsn == NULL) return;

 dsn[0] = svector<2>(-1,-1);
 dsn[1] = svector<2>(1,0);
 dsn[2] = svector<2>(0,1);
 
 dsn[3][0] = ata*omema - bb;
 dsn[3][1] = eta*omema - bb;

 return;
}

// This buble has the Lagrange interpolation property

void tri3BubbleInterp(svector<2> xi, double sn[4], svector<2> dsn[4])
{
  double eta = xi[0], ata = xi[1], omema = 1.0-eta-ata, bb = 9*omema*eta*ata;

  sn[0] = omema - bb;  
  sn[1] = eta   - bb;    
  sn[2] = ata   - bb;    

  sn[3] = 3.0 * bb;

 if(dsn == NULL) return;

 svector<2> db(9*ata*(omema - eta), 9*eta*(omema - ata));

 dsn[0] = svector<2>(-1,-1) - db;
 dsn[1] = svector<2>(1,0)   - db;
 dsn[2] = svector<2>(0,1)   - db;
 
 dsn[3] = 3.0 * db;

 return;
}

void tri6(svector<2> xi, double sn[], svector<2> dsn[])
{
  double eta = xi[0], ata = xi[1], omema = 1.0-eta-ata, hmema = 0.5-eta-ata;

  sn[0] = 2*omema*hmema;
  sn[1] = 2* eta *(eta-0.5);
  sn[2] = 2* ata *(ata-0.5);

  sn[3] = 4*eta*omema;
  sn[4] = 4*eta*ata;
  sn[5] = 4*ata*omema;

 if(dsn == NULL) return;
 
  dsn[0] = -2*(hmema+omema) * svector<2>(1.0,1.0);
  dsn[1] = svector<2>(4*eta-1, 0.0);
  dsn[2] = svector<2>(0.0, 4*ata-1);

  dsn[3] = 4 * svector<2>(1-2*eta-ata, -eta);
  dsn[4] = 4 * svector<2>(ata,eta);
  dsn[5] = 4 * svector<2>(-ata, 1-eta-2*ata);

  return;
}

void tri10Permutation(const int *nodes, int *dof)
{
  int temp;

  // permute the edge dof

  if(nodes[0] > nodes[1]) {temp = dof[3]; dof[3] = dof[4]; dof[4] = temp;}
  if(nodes[1] > nodes[2]) {temp = dof[5]; dof[5] = dof[6]; dof[6] = temp;}
  if(nodes[2] > nodes[0]) {temp = dof[7]; dof[7] = dof[8]; dof[8] = temp;}

  return;
}

 void tri10(svector<2> pt, double psi[10], svector<2> dpsi[10])
{
  // cubic serendipity square

  double xi = pt[0], eta = pt[1], eta2 = eta*eta, xi2 = xi*xi;
  double sqrt5 = sqrt(5.0);

psi[0] = -(xi-1+eta)*(5*xi2-5*xi+11*xi*eta-5*eta+5*eta2+1);
psi[1] = xi*(1-5*xi+eta+5*xi2-xi*eta-eta2);
psi[2] = -eta*(-1-xi+5*eta+xi2+xi*eta-5*eta2);
psi[3] = 0.5*sqrt5*(10*xi-sqrt5-5+3*eta*sqrt5+5*eta)*(xi-1+eta)*xi;
psi[4] = -0.5*sqrt5*(10*xi+sqrt5-5-3*eta*sqrt5+5*eta)*(xi-1+eta)*xi;
psi[5] = 1.25*(sqrt5+3)*(2*xi-3-3*eta*sqrt5+7*eta+sqrt5)*xi*eta;
psi[6] = -1.25*(sqrt5-3)*(2*xi-3+7*eta+3*eta*sqrt5-sqrt5)*xi*eta;
psi[7] = -1.25*(sqrt5-3)*(2*xi-3*eta*sqrt5-5*eta+sqrt5+1)*(xi-1+eta)*eta;
psi[8] = 1.25*(sqrt5+3)*(2*xi+3*eta*sqrt5-5*eta-sqrt5+1)*(xi-1+eta)*eta;
psi[9] = -27*xi*eta*(xi-1+eta);

 if(dpsi == NULL) return;

dpsi[0] = svector<2>(-15*xi2+20*xi-32*xi*eta+21*eta-16*eta2-6,
		  -16*xi2+21*xi-32*xi*eta+20*eta-15*eta2-6);
dpsi[1] = svector<2>(1-10*xi+eta+15*xi2-2*xi*eta-eta2,-xi*(-1+xi+2*eta));
dpsi[2] = svector<2>(-eta*(-1+2*xi+eta),1+xi-10*eta-xi2-2*xi*eta+15*eta2);
dpsi[3] = svector<2>(0.5*sqrt5*(30*xi2-30*xi+30*xi*eta-2*xi*sqrt5+6*xi*eta*sqrt5+sqrt5-4*eta*sqrt5+5-10*eta+3*eta2*sqrt5+5*eta2),
		  1.25*(sqrt5+1)*(6*xi-3-sqrt5+2*eta+2*eta*sqrt5)*xi);
dpsi[4] = svector<2>(-0.5*sqrt5*(30*xi2-30*xi+30*xi*eta+2*xi*sqrt5-6*xi*eta*sqrt5-sqrt5+4*eta*sqrt5+5-10*eta-3*eta2*sqrt5+5*eta2),
		  -1.25*(sqrt5-1)*(6*xi-3+sqrt5+2*eta-2*eta*sqrt5)*xi);
dpsi[5] = svector<2>(1.25*(sqrt5+3)*(7*eta-3*eta*sqrt5-3+4*xi+sqrt5)*eta,
		  1.25*(sqrt5+3)*(2*xi+sqrt5-3-6*eta*sqrt5+14*eta)*xi);
dpsi[6] = svector<2>(-1.25*(sqrt5-3)*(3*eta*sqrt5+7*eta-3+4*xi-sqrt5)*eta,
		  -1.25*(sqrt5-3)*(2*xi-3-sqrt5+14*eta+6*eta*sqrt5)*xi);
dpsi[7] = svector<2>(-1.25*(sqrt5-3)*(-3*eta-3*eta*sqrt5-1+sqrt5+4*xi)*eta,
		  -1.25*(sqrt5-3)*(-6*xi*eta-6*xi*eta*sqrt5+12*eta+8*eta*sqrt5-15*eta2-9*eta2*sqrt5+2*xi2+xi*sqrt5-xi-1-sqrt5));
dpsi[8] = svector<2>(1.25*(sqrt5+3)*(-sqrt5-1+4*xi-3*eta+3*eta*sqrt5)*eta,
		  1.25*(sqrt5+3)*(-6*xi*eta+6*xi*eta*sqrt5+12*eta-8*eta*sqrt5-15*eta2+9*eta2*sqrt5+2*xi2-xi*sqrt5-xi-1+sqrt5));
dpsi[9] = -27 * svector<2>(eta*(-1+2*xi+eta), xi*(-1+xi+2*eta));
  return;
}

void quad4(svector<2> xi, double psi[4], svector<2> dpsi[4])
{
  double eta = xi[0], ata = xi[1];
  double ope = (1+eta)/2, opa = (1+ata)/2;
  double ome = (1-eta)/2, oma = (1-ata)/2;

  psi[0] = ome * oma;
  psi[1] = ope * oma;
  psi[2] = ope * opa;
  psi[3] = ome * opa;

  if(dpsi == NULL) return;
 
  dpsi[0] = svector<2>(oma /-2, ome /-2);
  dpsi[1] = svector<2>(oma / 2, ope /-2);
  dpsi[2] = svector<2>(opa / 2, ope / 2);
  dpsi[3] = svector<2>(opa /-2, ome / 2);
  
  return;
}

void quad4Bubble(svector<2> xi, double psi[5], svector<2> dpsi[5])
{
  double eta = xi[0], ata = xi[1];
  double ope = (1+eta)/2, opa = (1+ata)/2;
  double ome = (1-eta)/2, oma = (1-ata)/2;

  double bb = (1-eta*eta)*(1-ata*ata), opepa = 1 + eta + ata;

  psi[0] = ome * oma;
  psi[1] = ope * oma;
  psi[2] = ope * opa;
  psi[3] = ome * opa;

  psi[4] = bb * opepa;  // bubble

  if(dpsi == NULL) return;
 
  dpsi[0] = svector<2>(oma /-2, ome /-2);
  dpsi[1] = svector<2>(oma / 2, ope /-2);
  dpsi[2] = svector<2>(opa / 2, ope / 2);
  dpsi[3] = svector<2>(opa /-2, ome / 2);

  dpsi[4][0] = bb - 2*eta*(1-ata*ata) * opepa;  // bubble
  dpsi[4][1] = bb - 2*ata*(1-eta*eta) * opepa;  // bubble
  
  return;
}

void quad8(svector<2> xi, double phi[8], svector<2> dphi[8])
{
  double eta = xi[0], ata = xi[1];
  double ope = 1+eta, opa = 1+ata, ome = 1-eta, oma = 1-ata;
  

  phi[0] = -0.25*ome*oma*(eta+ata+1);
  phi[1] =  0.25*ope*oma*(eta-ata-1);
  phi[2] =  0.25*ope*opa*(eta+ata-1);
  phi[3] = -0.25*ome*opa*(eta-ata+1);
  phi[4] =  0.5*ope*ome*oma;
  phi[5] =  0.5*oma*opa*ope;
  phi[6] =  0.5*ome*ope*opa;
  phi[7] =  0.5*oma*opa*ome;

  if(dphi == NULL) return;
 
  dphi[0] =  0.25* svector<2>(oma*(2*eta+ata), ome*(2*ata+eta));
  dphi[1] =  0.25* svector<2>(oma*(2*eta-ata), ope*(2*ata-eta));
  dphi[2] =  0.25* svector<2>(opa*(2*eta+ata), ope*(2*ata+eta));
  dphi[3] =  0.25* svector<2>(opa*(2*eta-ata), ome*(2*ata-eta));
  dphi[4] = -svector<2>(eta*oma, 0.5*ope*ome);
  dphi[5] =  svector<2>(0.5*oma*opa, -ata*ope);
  dphi[6] =  svector<2>(-eta*opa, 0.5*ome*ope);
  dphi[7] = -svector<2>(0.5*oma*opa,  ata*ome);
  
  return;
}   

void quad9(svector<2> xi, double phi[9], svector<2> dphi[9])
{
  double phi0e,phi1e,phi2e,dphi0e,dphi1e,dphi2e;
  double phi0a,phi1a,phi2a,dphi0a,dphi1a,dphi2a;
  double eta = xi[0], ata = xi[1];
  double ope = 1+eta, opa = 1+ata, ome = 1-eta, oma = 1-ata;
  
  phi0e =-eta * ome / 2; dphi0e =-0.5 + eta;
  phi1e = ope * ome;     dphi1e =-2.0 * eta;
  phi2e = eta * ope / 2; dphi2e = 0.5 + eta;
  
  phi0a =-ata * oma / 2; dphi0a =-0.5 + ata;
  phi1a = opa * oma;     dphi1a =-2.0 * ata;
  phi2a = ata * opa / 2; dphi2a = 0.5 + ata;
  
  phi[0] = phi0e*phi0a;
  phi[1] = phi2e*phi0a;
  phi[2] = phi2e*phi2a;
  phi[3] = phi0e*phi2a;
  
  phi[4] = phi1e*phi0a;
  phi[5] = phi2e*phi1a;
  phi[6] = phi1e*phi2a;
  phi[7] = phi0e*phi1a;
  
  phi[8] = phi1e*phi1a;  dphi[8] = svector<2>(dphi1e*phi1a, phi1e*dphi1a);
  
  if(dphi == NULL) return;
  
  dphi[0] = svector<2>(dphi0e*phi0a, phi0e*dphi0a);
  dphi[1] = svector<2>(dphi2e*phi0a, phi2e*dphi0a);
  dphi[2] = svector<2>(dphi2e*phi2a, phi2e*dphi2a);
  dphi[3] = svector<2>(dphi0e*phi2a, phi0e*dphi2a);
  
  dphi[4] = svector<2>(dphi1e*phi0a, phi1e*dphi0a);
  dphi[5] = svector<2>(dphi2e*phi1a, phi2e*dphi1a);
  dphi[6] = svector<2>(dphi1e*phi2a, phi1e*dphi2a);
  dphi[7] = svector<2>(dphi0e*phi1a, phi0e*dphi1a);
  
  dphi[8] = svector<2>(dphi1e*phi1a, phi1e*dphi1a);
  
  return;
}

void quad12(svector<2> xi, double psi[12], svector<2> dpsi[12])
{
  // cubic serendipity square

  double eta = xi[0], ata = xi[1], eta2 = eta*eta, ata2 = ata*ata;
  double one16 = 1.0/16.0, sqrt5 = sqrt(5.0), rt5o16 = sqrt5*one16;

  psi[0] =  one16*(-1+eta)*(-1+ata)*(5*ata2+5*eta2-6);
  psi[1] = -one16*( 1+eta)*(-1+ata)*(5*ata2+5*eta2-6);
  psi[2] =  one16*( 1+eta)*( 1+ata)*(5*ata2+5*eta2-6);
  psi[3] = -one16*(-1+eta)*( 1+ata)*(5*ata2+5*eta2-6);

  psi[4] =  rt5o16*(-1+ata)*(-1+eta)*(1+eta)*(-5*eta+sqrt5);
  psi[5] =  rt5o16*(-1+ata)*(-1+eta)*(1+eta)*(5*eta+sqrt5);
  psi[6] = -rt5o16*(-1+ata)*(-5*ata+sqrt5)*(1+ata)*(1+eta);
  psi[7] = -rt5o16*(-1+ata)*( 5*ata+sqrt5)*(1+ata)*(1+eta);
  psi[8] = -rt5o16*( 1+ata)*(-1+eta)*(1+eta)*(5*eta+sqrt5);
  psi[9] = -rt5o16*( 1+ata)*(-5*eta+sqrt5)*(-1+eta)*(1+eta);
  psi[10] = rt5o16*(5*ata+sqrt5)*(-1+ata)*(1+ata)*(-1+eta);
  psi[11] = rt5o16*(-1+ata)*(1+ata)*(-5*ata+sqrt5)*(-1+eta);

  if(dpsi == NULL) return;
  
  dpsi[0] = svector<2>(one16*(-1+ata)*(5*ata2+15*eta2-10*eta-6), 
		       one16*(-1+eta)*(5*eta2-6+15*ata2-10*ata));
  dpsi[1] = svector<2>(-one16*(-1+ata)*(5*ata2+15*eta2+10*eta-6), 
		       -one16*(1+eta)*(5*eta2-6+15*ata2-10*ata));
  dpsi[2] = svector<2>(one16*(1+ata)*(5*ata2+15*eta2+10*eta-6), 
		       one16*(1+eta)*(5*eta2+15*ata2-6+10*ata));
  dpsi[3] = svector<2>(-one16*(1+ata)*(5*ata2+15*eta2-10*eta-6), 
		       -one16*(-1+eta)*(5*eta2+15*ata2-6+10*ata));
  dpsi[4] = svector<2>(rt5o16*(-1+ata)*(-3*eta+sqrt5)*(5*eta+sqrt5), 
		       rt5o16*(-1+eta)*(1+eta)*(-5*eta+sqrt5));
  dpsi[5] = svector<2>(-rt5o16*(-1+ata)*(3*eta+sqrt5)*(-5*eta+sqrt5), 
		       rt5o16*(-1+eta)*(1+eta)*(5*eta+sqrt5));
  dpsi[6] = svector<2>(- rt5o16*(-1+ata)*(1+ata)*(-5*ata+sqrt5), 
		       -rt5o16*(-3*ata+sqrt5)*(5*ata+sqrt5)*(1+eta));
  dpsi[7] = svector<2>(-rt5o16*(-1+ata)*(5*ata+sqrt5)*(1+ata), 
		       rt5o16*(3*ata+sqrt5)*(-5*ata+sqrt5)*(1+eta));
  dpsi[8] = svector<2>(rt5o16*(1+ata)*(3*eta+sqrt5)*(-5*eta+sqrt5), 
		       -rt5o16*(-1+eta)*(1+eta)*(5*eta+sqrt5));
  dpsi[9] = svector<2>(-rt5o16*(1+ata)*(-3*eta+sqrt5)*(5*eta+sqrt5), 
		       -rt5o16*(-1+eta)*(1+eta)*(-5*eta+sqrt5));
  dpsi[10] = svector<2>(rt5o16*(-1+ata)*(5*ata+sqrt5)*(1+ata), 
			-rt5o16*(3*ata+sqrt5)*(-5*ata+sqrt5)*(-1+eta));
  dpsi[11] = svector<2>(rt5o16*(-1+ata)*(1+ata)*(-5*ata+sqrt5), 
			rt5o16*(-3*ata+sqrt5)*(5*ata+sqrt5)*(-1+eta));

  return;
}

void quad12Permutation(const int *nodes, int *dof)
{
  int temp;

  // permute the edge dof

  if(nodes[0] > nodes[1]) {temp = dof[ 4]; dof[ 4] = dof[ 5]; dof[ 5] = temp;}
  if(nodes[1] > nodes[2]) {temp = dof[ 6]; dof[ 6] = dof[ 7]; dof[ 7] = temp;}
  if(nodes[2] > nodes[3]) {temp = dof[ 8]; dof[ 8] = dof[ 9]; dof[ 9] = temp;}
  if(nodes[3] > nodes[0]) {temp = dof[10]; dof[10] = dof[11]; dof[11] = temp;}

  return;
}
	   
void quad16(svector<2> xi, double phi[16], svector<2> dphi[16])
{
  int i,j;
  int idx[4][4] = {{0,3,11,10},{1,2,6,7},{4,9,12,15},{5,8,13,14}};

  svector<1> eta(xi[0]), ata(xi[1]);
  double phix[4], dphix[4], phiy[4], dphiy[4];

  svector<1> vdphix[4], vdphiy[4];

  line4(eta, phix, vdphix);
  line4(ata, phiy, vdphiy);

  for(i = 0; i < 4; i++)
    {
      dphix[i] = vdphix[i][0];
      dphiy[i] = vdphiy[i][0];
    }

  for(i = 0; i < 4; i++)
  for(j = 0; j < 4; j++)
    {
      phi[idx[i][j]] = phix[i] * phiy[j];

      if(dphi != NULL)
	dphi[idx[i][j]] = svector<2>(dphix[i]*phiy[j], phix[i]*dphiy[j]);
    }

  return;
}

void quad16Permutation(const int *nodes, int *dof)
{
  int temp;

  // permute the edge dof

  if(nodes[0] > nodes[1]) {temp = dof[ 4]; dof[ 4] = dof[ 5]; dof[ 5] = temp;}
  if(nodes[1] > nodes[2]) {temp = dof[ 6]; dof[ 6] = dof[ 7]; dof[ 7] = temp;}
  if(nodes[2] > nodes[3]) {temp = dof[ 8]; dof[ 8] = dof[ 9]; dof[ 9] = temp;}
  if(nodes[3] > nodes[0]) {temp = dof[10]; dof[10] = dof[11]; dof[11] = temp;}

  // permute the cell dof

  int nn[4];
  for(int i = 0; i < 4; i++) nn[i] = nodes[i];
  sortDof(4, nn, dof+12);

  return;
}

void quad17(svector<2> xi, double psi[17], svector<2> dpsi[17])
{
  // quartic serendipity square

  double eta = xi[0], ata = xi[1], sqrt21 = sqrt(21.0);
  double eta2 = eta *eta, ata2 = ata *ata;
  double eta3 = eta2*eta, ata3 = ata2*ata;
  double ome = 1-eta, ope = 1+eta, oma = 1-ata, opa = 1+ata;

  psi[0] = -(1.0/16)*oma*ome*(7*eta3-4*ata*eta-7*eta+7*ata3-7*ata);
  psi[1] =  (1.0/16)*oma*ope*(7*eta3-4*ata*eta-7*eta-7*ata3+7*ata);
  psi[2] =  (1.0/16)*opa*ope*(7*eta3+4*ata*eta-7*eta+7*ata3-7*ata);
  psi[3] = -(1.0/16)*opa*ome*(7*eta3+4*ata*eta-7*eta-7*ata3+7*ata);
  psi[4] = -(7.0/48)*(-7*eta+sqrt21)*ope*ome*eta*oma;
  psi[5] = -(1.0/6)*oma*ome*ope*(7*eta2+3*ata);
  psi[6] =  (7.0/48)*(7*eta+sqrt21)*ope*ome*eta*oma;
  psi[7] = -(7.0/48)*ope*oma*opa*ata*(-7*ata+sqrt21);
  psi[8] =  (1.0/6)*oma*opa*ope*(3*eta-7*ata2);
  psi[9] =  (7.0/48)*ope*oma*opa*ata*(7*ata+sqrt21);
  psi[10] = (7.0/48)*ope*(7*eta+sqrt21)*ome*eta*opa;
  psi[11] = -(1.0/6)*opa*ome*ope*(7*eta2-3*ata);
  psi[12] = -(7.0/48)*(-7*eta+sqrt21)*ope*ome*eta*opa;
  psi[13] =  (7.0/48)*ome*oma*opa*(7*ata+sqrt21)*ata;
  psi[14] = -(1.0/6)*oma*opa*ome*(3*eta+7*ata2);
  psi[15] = -(7.0/48)*ome*(-7*ata+sqrt21)*oma*opa*ata;
  psi[16] =  oma*opa*ome*ope;
  
  if(dpsi == NULL) return;
 
  dpsi[0] = svector<2>(-(1.0/16)*(-1+ata)*(7*ata3-3*ata-8*ata*eta-21*eta2+7-14*eta+28*eta3),
		       -(1.0/16)*(-1+eta)*(7*eta3-8*ata*eta-3*eta+28*ata3-14*ata+7-21*ata2));
  dpsi[1] = svector<2>(-(1.0/16)*(-1+ata)*(-7*ata3+3*ata-8*ata*eta+21*eta2-7-14*eta+28*eta3), 
		       -(1.0/16)*(1+eta)*(7*eta3-8*ata*eta-3*eta+14*ata+21*ata2-7-28*ata3));
  dpsi[2] = svector<2>((1.0/16)*(1+ata)*(7*ata3-3*ata+8*ata*eta+28*eta3+21*eta2-14*eta-7), 
		       (1.0/16)*(1+eta)*(7*eta3+8*ata*eta-3*eta-14*ata+21*ata2-7+28*ata3));
  dpsi[3] = svector<2>((1.0/16)*(1+ata)*(-7*ata3+3*ata+8*ata*eta+28*eta3-21*eta2-14*eta+7), 
		       (1.0/16)*(-1+eta)*(7*eta3+8*ata*eta-3*eta+14*ata-21*ata2+7-28*ata3));
  dpsi[4] = svector<2>(-(7.0/48)*(-4*eta2+4.582575695*eta-1)*(7*eta+4.582575695)*(-1+ata), 
		       -(7.0/48)*eta*(1+eta)*(-7*eta+4.582575695)*(-1+eta));
  dpsi[5] = svector<2>(-(1.0/3)*eta*(-1+ata)*(3*ata-7+14*eta2), 
		       -(1.0/6)*(-1+eta)*(1+eta)*(7*eta2+6*ata-3));
  dpsi[6] = svector<2>(-(7.0/48)*(4*eta2+4.582575695*eta+1)*(-7*eta+4.582575695)*(-1+ata), 
	  	        (7.0/48)*eta*(1+eta)*(7*eta+4.582575695)*(-1+eta));
  dpsi[7] = svector<2>((7.0/48)*(-7*ata+4.582575695)*ata*(-1+ata)*(1+ata), 
		       (7.0/48)*(1+eta)*(-4*ata2+4.582575695*ata-1)*(7*ata+4.582575695));
  dpsi[8] = svector<2>(-(1.0/6)*(-1+ata)*(1+ata)*(-7*ata2+3+6*eta), 
		       -(1.0/3)*ata*(1+eta)*(3*eta+7-14*ata2));
  dpsi[9] = svector<2>(-(7.0/48)*ata*(7*ata+4.582575695)*(-1+ata)*(1+ata), 
		       (7.0/48)*(1+eta)*(4*ata2+4.582575695*ata+1)*(-7*ata+4.582575695));
  dpsi[10] = svector<2>((7.0/48)*(4*eta2+4.582575695*eta+1)*(-7*eta+4.582575695)*(1+ata), 
			-(7.0/48)*eta*(1+eta)*(7*eta+4.582575695)*(-1+eta));
  dpsi[11] = svector<2>((1.0/3)*eta*(1+ata)*(-3*ata+14*eta2-7), 
			(1.0/6)*(-1+eta)*(1+eta)*(7*eta2-6*ata-3));
  dpsi[12] = svector<2>((7.0/48)*(-4*eta2+4.582575695*eta-1)*(7*eta+4.582575695)*(1+ata), 
			(7.0/48)*eta*(1+eta)*(-7*eta+4.582575695)*(-1+eta));
  dpsi[13] = svector<2>((7.0/48)*ata*(7*ata+4.582575695)*(-1+ata)*(1+ata), 
			-(7.0/48)*(-1+eta)*(4*ata2+4.582575695*ata+1)*(-7*ata+4.582575695));
  dpsi[14] = svector<2>(-(1.0/6)*(-1+ata)*(1+ata)*(7*ata2-3+6*eta), 
			-(1.0/3)*ata*(-1+eta)*(3*eta-7+14*ata2));
  dpsi[15] = svector<2>(-(7.0/48)*(-7*ata+4.582575695)*ata*(-1+ata)*(1+ata), 
			-(7.0/48)*(-1+eta)*(-4*ata2+4.582575695*ata-1)*(7*ata+4.582575695));
  dpsi[16] = 2*svector<2>(eta*(-1+ata)*(1+ata), ata*(-1+eta)*(1+eta));
  
  return;
}

void quad17Permutation(const int *nodes, int *dof)
{
  int temp;

  // permute the edge dof

  if(nodes[0] > nodes[1]) {temp = dof[ 4]; dof[ 4] = dof[ 6]; dof[ 6] = temp;}
  if(nodes[1] > nodes[2]) {temp = dof[ 7]; dof[ 7] = dof[ 9]; dof[ 9] = temp;}
  if(nodes[2] > nodes[3]) {temp = dof[10]; dof[10] = dof[12]; dof[12] = temp;}
  if(nodes[3] > nodes[0]) {temp = dof[13]; dof[13] = dof[15]; dof[15] = temp;}

  return;
}

void quad24(svector<2> xi, double psi[24], svector<2> dpsi[24])
{
  // quintic serendipity square

  double s = xi[0], t = xi[1], sqrt3 = sqrt(3.0), sqrt30 = sqrt(30.0);
  double s2 = s *s,  t2 = t *t;
  double s3 = s2*s,  t3 = t2*t;
  double s4 = s2*s2, t4 = t2*t2;
  double sm1 = s-1, sp1 = s+1, tm1 = t-1, tp1 = t+1;
  double pt1 = sqrt(525.0 + 70.0*sqrt30), pt2 = sqrt(525.0 - 70.0*sqrt30);

psi[0] = (1.0/5408)*(2424-7832*s2-7832*t2+5915*s4+5915*t4+2762*s2*t2)*tm1*sm1;
psi[1] = -(1.0/5408)*(2424-7832*s2-7832*t2+5915*s4+5915*t4+2762*s2*t2)*tm1*sp1;
psi[2] = (1.0/5408)*(2424-7832*s2-7832*t2+5915*s4+5915*t4+2762*s2*t2)*tp1*sp1;
psi[3] = -(1.0/5408)*(2424-7832*s2-7832*t2+5915*s4+5915*t4+2762*s2*t2)*tp1*sm1;
psi[4] = -(1.0/218312640)*pt1*(113*sqrt30-711)*(-78*sqrt30-396+416*sqrt30*s2+1755*s2+357*t2)*tm1*sm1*sp1*(-35*s+pt1);
psi[5] = -(1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30-1755*s2+416*sqrt30*s2-357*t2)*tm1*(-35*s+pt2)*sm1*sp1;
psi[6] = -(1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30-1755*s2+416*sqrt30*s2-357*t2)*tm1*sm1*sp1*(35*s+pt2);
psi[7] = -(1.0/218312640)*pt1*(113*sqrt30-711)*(-78*sqrt30-396+416*sqrt30*s2+1755*s2+357*t2)*tm1*sm1*sp1*(35*s+pt1);
psi[8] = -(1.0/3057600)*pt1*(3*sqrt30-5)*(-6*sqrt30-60-135*s2+32*sqrt30*s2+455*t2)*tp1*tm1*(-35*t+pt1)*sp1;
psi[9] = -(1.0/3057600)*pt2*(3*sqrt30+5)*(-455*t2+60-6*sqrt30+32*sqrt30*s2+135*s2)*tp1*tm1*(-35*t+pt2)*sp1;
psi[10] = -(1.0/3057600)*pt2*(3*sqrt30+5)*(-455*t2+60-6*sqrt30+32*sqrt30*s2+135*s2)*tp1*tm1*(35*t+pt2)*sp1;
psi[11] = -(1.0/3057600)*pt1*(3*sqrt30-5)*(-6*sqrt30-60-135*s2+32*sqrt30*s2+455*t2)*tp1*(35*t+pt1)*tm1*sp1;
psi[12] = (1.0/218312640)*pt1*(113*sqrt30-711)*(-78*sqrt30-396+416*sqrt30*s2+1755*s2+357*t2)*tp1*sm1*sp1*(35*s+pt1);
psi[13] = (1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30-1755*s2+416*sqrt30*s2-357*t2)*tp1*sm1*sp1*(35*s+pt2);
psi[14] = (1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30-1755*s2+416*sqrt30*s2-357*t2)*tp1*(-35*s+pt2)*sm1*sp1;
psi[15] = (1.0/218312640)*pt1*(113*sqrt30-711)*(-78*sqrt30-396+416*sqrt30*s2+1755*s2+357*t2)*tp1*sm1*sp1*(-35*s+pt1);
psi[16] = (1.0/3057600)*pt1*(3*sqrt30-5)*(-6*sqrt30-60-135*s2+32*sqrt30*s2+455*t2)*tp1*tm1*(35*t+pt1)*sm1;
psi[17] = (1.0/3057600)*pt2*(3*sqrt30+5)*(-455*t2+60-6*sqrt30+32*sqrt30*s2+135*s2)*(35*t+pt2)*tp1*tm1*sm1;
psi[18] = (1.0/3057600)*pt2*(3*sqrt30+5)*(-455*t2+60-6*sqrt30+32*sqrt30*s2+135*s2)*tp1*tm1*(-35*t+pt2)*sm1;
psi[19] = (1.0/3057600)*pt1*(3*sqrt30-5)*(-6*sqrt30-60-135*s2+32*sqrt30*s2+455*t2)*tp1*(-35*t+pt1)*tm1*sm1;
psi[20] = (64.0/507)*tm1*tp1*(-4*t+sqrt3)*sm1*sp1*(-4*s+sqrt3);
psi[21] = (64.0/507)*tp1*(-4*t+sqrt3)*tm1*(4*s+sqrt3)*sm1*sp1;
psi[22] = (64.0/507)*tm1*tp1*(4*t+sqrt3)*sm1*sp1*(4*s+sqrt3);
psi[23] = (64.0/507)*tp1*(4*t+sqrt3)*tm1*sm1*sp1*(-4*s+sqrt3);

  if(dpsi == NULL) return;

dpsi[0] = svector<2>((1.0/5408)*tm1*(5915*t4-5524*s*t2+8286*s2*t2-7832*t2-23496*s2+15664*s+29575*s4-23660*s3+2424), (1.0/5408)*sm1*(5915*s4+8286*s2*t2-5524*s2*t-7832*s2-23496*t2-23660*t3+29575*t4+2424+15664*t));
dpsi[1] = svector<2>(-(1.0/5408)*tm1*(5915*t4+5524*s*t2+8286*s2*t2-7832*t2-23496*s2-15664*s+29575*s4+23660*s3+2424), -(1.0/5408)*sp1*(5915*s4+8286*s2*t2-5524*s2*t-7832*s2-23496*t2-23660*t3+29575*t4+2424+15664*t));
dpsi[2] = svector<2>((1.0/5408)*tp1*(5915*t4+5524*s*t2+8286*s2*t2-7832*t2-23496*s2-15664*s+29575*s4+23660*s3+2424), (1.0/5408)*sp1*(5915*s4+8286*s2*t2+5524*s2*t-7832*s2+2424-15664*t-23496*t2+23660*t3+29575*t4));
dpsi[3] = svector<2>(-(1.0/5408)*tp1*(5915*t4-5524*s*t2+8286*s2*t2-7832*t2-23496*s2+15664*s+29575*s4-23660*s3+2424), -(1.0/5408)*sm1*(5915*s4+8286*s2*t2+5524*s2*t-7832*s2+2424-15664*t-23496*t2+23660*t3+29575*t4));
dpsi[4] = svector<2>(-(1.0/218312640)*pt1*(113*sqrt30-711)*(225855*s2+12495*t2-307125*s4-2730*sqrt30-37485*s2*t2-4302*s*pt1+7020*s3*pt1-13860+1664*sqrt30*s3*pt1+714*pt1*t2*s-988*sqrt30*s*pt1-72800*sqrt30*s4+51870*sqrt30*s2)*tm1, -(1.0/218312640)*pt1*(113*sqrt30-711)*(-396-78*sqrt30+416*sqrt30*s2+1755*s2+1071*t2-714*t)*sm1*sp1*(-35*s+pt1));
dpsi[5] = svector<2>(-(1.0/218312640)*pt2*(711+113*sqrt30)*(-714*t2*s*pt2-225855*s2-12495*t2+307125*s4-2730*sqrt30+37485*s2*t2-988*sqrt30*s*pt2+13860-72800*sqrt30*s4+4302*s*pt2+51870*sqrt30*s2+1664*sqrt30*s3*pt2-7020*s3*pt2)*tm1, -(1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30+416*sqrt30*s2-1755*s2-1071*t2+714*t)*sm1*sp1*(-35*s+pt2));
dpsi[6] = svector<2>(-(1.0/218312640)*pt2*(711+113*sqrt30)*(-714*t2*s*pt2+225855*s2+12495*t2-307125*s4+2730*sqrt30-37485*s2*t2-988*sqrt30*s*pt2-13860+72800*sqrt30*s4+4302*s*pt2-51870*sqrt30*s2+1664*sqrt30*s3*pt2-7020*s3*pt2)*tm1, -(1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30+416*sqrt30*s2-1755*s2-1071*t2+714*t)*sm1*(35*s+pt2)*sp1);
dpsi[7] = svector<2>(-(1.0/218312640)*pt1*(113*sqrt30-711)*(-225855*s2-12495*t2+307125*s4+2730*sqrt30+37485*s2*t2-4302*s*pt1+7020*s3*pt1+13860+1664*sqrt30*s3*pt1+714*pt1*t2*s-988*sqrt30*s*pt1+72800*sqrt30*s4-51870*sqrt30*s2)*tm1, -(1.0/218312640)*pt1*(113*sqrt30-711)*(-396-78*sqrt30+416*sqrt30*s2+1755*s2+1071*t2-714*t)*sm1*sp1*(35*s+pt1));
dpsi[8] = svector<2>(-(1.0/3057600)*pt1*(3*sqrt30-5)*(-60-6*sqrt30-405*s2+96*sqrt30*s2+455*t2-270*s+64*sqrt30*s)*tp1*tm1*(-35*t+pt1), -(1.0/3057600)*pt1*(3*sqrt30-5)*(-2100-210*sqrt30+1120*sqrt30*s2-4725*s2+630*sqrt30*t2+54075*t2+1820*pt1*t3-79625*t4-1030*pt1*t-12*pt1*sqrt30*t-270*s2*pt1*t+64*sqrt30*s2*pt1*t+14175*s2*t2-3360*sqrt30*s2*t2)*sp1);
dpsi[9] = svector<2>(-(1.0/3057600)*pt2*(3*sqrt30+5)*(60-6*sqrt30+405*s2+96*sqrt30*s2-455*t2+64*sqrt30*s+270*s)*tp1*tm1*(-35*t+pt2), -(1.0/3057600)*pt2*(3*sqrt30+5)*(2100-210*sqrt30+4725*s2+1120*sqrt30*s2-54075*t2+630*sqrt30*t2-1820*pt2*t3+79625*t4+1030*pt2*t-12*pt2*sqrt30*t+64*sqrt30*s2*pt2*t+270*s2*pt2*t-14175*s2*t2-3360*sqrt30*s2*t2)*sp1);
dpsi[10] = svector<2>(-(1.0/3057600)*pt2*(3*sqrt30+5)*(60-6*sqrt30+405*s2+96*sqrt30*s2-455*t2+64*sqrt30*s+270*s)*tp1*(35*t+pt2)*tm1, -(1.0/3057600)*pt2*(3*sqrt30+5)*(-2100+210*sqrt30-4725*s2-1120*sqrt30*s2+54075*t2-630*sqrt30*t2-1820*pt2*t3-79625*t4+1030*pt2*t-12*pt2*sqrt30*t+64*sqrt30*s2*pt2*t+270*s2*pt2*t+14175*s2*t2+3360*sqrt30*s2*t2)*sp1);
dpsi[11] = svector<2>(-(1.0/3057600)*pt1*(3*sqrt30-5)*(-60-6*sqrt30-405*s2+96*sqrt30*s2+455*t2-270*s+64*sqrt30*s)*tp1*(35*t+pt1)*tm1, -(1.0/3057600)*pt1*(3*sqrt30-5)*(2100+210*sqrt30-1120*sqrt30*s2+4725*s2-630*sqrt30*t2-54075*t2+1820*pt1*t3+79625*t4-12*pt1*sqrt30*t-1030*pt1*t+64*sqrt30*s2*pt1*t-270*s2*pt1*t-14175*s2*t2+3360*sqrt30*s2*t2)*sp1);
dpsi[12] = svector<2>((1.0/218312640)*pt1*(113*sqrt30-711)*(-225855*s2-12495*t2+307125*s4+2730*sqrt30+37485*s2*t2-4302*s*pt1+7020*s3*pt1+13860+1664*sqrt30*s3*pt1+714*pt1*t2*s-988*sqrt30*s*pt1+72800*sqrt30*s4-51870*sqrt30*s2)*tp1, (1.0/218312640)*pt1*(113*sqrt30-711)*(-396-78*sqrt30+416*sqrt30*s2+1755*s2+1071*t2+714*t)*sm1*sp1*(35*s+pt1));
dpsi[13] = svector<2>((1.0/218312640)*pt2*(711+113*sqrt30)*(-714*t2*s*pt2+225855*s2+12495*t2-307125*s4+2730*sqrt30-37485*s2*t2-988*sqrt30*s*pt2-13860+72800*sqrt30*s4+4302*s*pt2-51870*sqrt30*s2+1664*sqrt30*s3*pt2-7020*s3*pt2)*tp1, (1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30+416*sqrt30*s2-1755*s2-1071*t2-714*t)*sm1*(35*s+pt2)*sp1);
dpsi[14] = svector<2>((1.0/218312640)*pt2*(711+113*sqrt30)*(-714*t2*s*pt2-225855*s2-12495*t2+307125*s4-2730*sqrt30+37485*s2*t2-988*sqrt30*s*pt2+13860-72800*sqrt30*s4+4302*s*pt2+51870*sqrt30*s2+1664*sqrt30*s3*pt2-7020*s3*pt2)*tp1, (1.0/218312640)*pt2*(711+113*sqrt30)*(396-78*sqrt30+416*sqrt30*s2-1755*s2-1071*t2-714*t)*sm1*sp1*(-35*s+pt2));
dpsi[15] = svector<2>((1.0/218312640)*pt1*(113*sqrt30-711)*(225855*s2+12495*t2-307125*s4-2730*sqrt30-37485*s2*t2-4302*s*pt1+7020*s3*pt1-13860+1664*sqrt30*s3*pt1+714*pt1*t2*s-988*sqrt30*s*pt1-72800*sqrt30*s4+51870*sqrt30*s2)*tp1, (1.0/218312640)*pt1*(113*sqrt30-711)*(-396-78*sqrt30+416*sqrt30*s2+1755*s2+1071*t2+714*t)*(-35*s+pt1)*sm1*sp1);
dpsi[16] = svector<2>((1.0/3057600)*pt1*(3*sqrt30-5)*(-60-6*sqrt30-405*s2+96*sqrt30*s2+455*t2+270*s-64*sqrt30*s)*tp1*(35*t+pt1)*tm1, (1.0/3057600)*pt1*(3*sqrt30-5)*(2100+210*sqrt30-1120*sqrt30*s2+4725*s2-630*sqrt30*t2-54075*t2+1820*pt1*t3+79625*t4-12*pt1*sqrt30*t-1030*pt1*t+64*sqrt30*s2*pt1*t-270*s2*pt1*t-14175*s2*t2+3360*sqrt30*s2*t2)*sm1);
dpsi[17] = svector<2>((1.0/3057600)*pt2*(3*sqrt30+5)*(60-6*sqrt30+405*s2+96*sqrt30*s2-455*t2-64*sqrt30*s-270*s)*tp1*tm1*(35*t+pt2), (1.0/3057600)*pt2*(3*sqrt30+5)*(-2100+210*sqrt30-4725*s2-1120*sqrt30*s2+54075*t2-630*sqrt30*t2-1820*pt2*t3-79625*t4+1030*pt2*t-12*pt2*sqrt30*t+64*sqrt30*s2*pt2*t+270*s2*pt2*t+14175*s2*t2+3360*sqrt30*s2*t2)*sm1);
dpsi[18] = svector<2>((1.0/3057600)*pt2*(3*sqrt30+5)*(60-6*sqrt30+405*s2+96*sqrt30*s2-455*t2-64*sqrt30*s-270*s)*tp1*tm1*(-35*t+pt2), (1.0/3057600)*pt2*(3*sqrt30+5)*(2100-210*sqrt30+4725*s2+1120*sqrt30*s2-54075*t2+630*sqrt30*t2-1820*pt2*t3+79625*t4+1030*pt2*t-12*pt2*sqrt30*t+64*sqrt30*s2*pt2*t+270*s2*pt2*t-14175*s2*t2-3360*sqrt30*s2*t2)*sm1);
dpsi[19] = svector<2>((1.0/3057600)*pt1*(3*sqrt30-5)*(-60-6*sqrt30-405*s2+96*sqrt30*s2+455*t2+270*s-64*sqrt30*s)*(-35*t+pt1)*tp1*tm1, (1.0/3057600)*pt1*(3*sqrt30-5)*(-2100-210*sqrt30+1120*sqrt30*s2-4725*s2+630*sqrt30*t2+54075*t2+1820*pt1*t3-79625*t4-1030*pt1*t-12*pt1*sqrt30*t-270*s2*pt1*t+64*sqrt30*s2*pt1*t+14175*s2*t2-3360*sqrt30*s2*t2)*sm1);
dpsi[20] = svector<2>((128.0/507)*tm1*tp1*(-4*t+sqrt3)*(s*sqrt3+2-6*s2), (128.0/507)*(2-6*t2+t*sqrt3)*sm1*sp1*(-4*s+sqrt3));
dpsi[21] = svector<2>((128.0/507)*tm1*tp1*(-4*t+sqrt3)*(s*sqrt3-2+6*s2), (128.0/507)*(2-6*t2+t*sqrt3)*sm1*sp1*(4*s+sqrt3));
dpsi[22] = svector<2>((128.0/507)*tm1*tp1*(4*t+sqrt3)*(s*sqrt3-2+6*s2), (128.0/507)*(t*sqrt3+6*t2-2)*sm1*sp1*(4*s+sqrt3));
dpsi[23] = svector<2>((128.0/507)*tm1*tp1*(4*t+sqrt3)*(s*sqrt3+2-6*s2), (128.0/507)*(t*sqrt3+6*t2-2)*sm1*sp1*(-4*s+sqrt3));
 
  return;
}

void quad24Permutation(const int *nodes, int *dof)
{
  int temp;

  // permute the edge and face dof

  if(nodes[0] > nodes[1]) 
    { 
      temp = dof[4]; dof[4] = dof[7]; dof[7] = temp;
      temp = dof[5]; dof[5] = dof[6]; dof[6] = temp;

      temp = dof[20]; dof[20] = dof[21]; dof[21] = temp;
      temp = dof[22]; dof[22] = dof[23]; dof[23] = temp;
    }
  if(nodes[1] > nodes[2])
    { 
      temp = dof[8]; dof[8] = dof[11]; dof[11] = temp;
      temp = dof[9]; dof[9] = dof[10]; dof[10] = temp;
    }
  if(nodes[2] > nodes[3])
    { 
      temp = dof[12]; dof[12] = dof[15]; dof[15] = temp;
      temp = dof[13]; dof[13] = dof[14]; dof[14] = temp;
    }
  if(nodes[3] > nodes[0])
    { 
      temp = dof[16]; dof[16] = dof[19]; dof[19] = temp;
      temp = dof[17]; dof[17] = dof[18]; dof[18] = temp;

      temp = dof[20]; dof[20] = dof[23]; dof[23] = temp;
      temp = dof[21]; dof[21] = dof[22]; dof[22] = temp;
    }

  return;
}

void tet4( const svector<3> xi, double phi[4] , svector<3> dphi[4])
{
  double xi3 = 1 - xi[0] - xi[1] - xi[2];
  
  phi[0] = xi3;  
  phi[1] = xi[0];
  phi[2] = xi[1];
  phi[3] = xi[2];
  
  if(dphi == NULL) return;
  
  dphi[0] = svector<3>(-1,-1,-1);
  dphi[1] = svector<3>( 1, 0, 0);
  dphi[2] = svector<3>( 0, 1, 0);
  dphi[3] = svector<3>( 0, 0, 1);
  
return;
}

void tet4Bubble(svector<3> xi, double phi[4] , svector<3> dphi[4])
{
  double xi3 = 1 - xi[0] - xi[1] - xi[2], bb = xi[0]*xi[1]*xi[2];
  
  phi[0] = xi3;  
  phi[1] = xi[0];
  phi[2] = xi[1];
  phi[3] = xi[2];
  
  phi[4] = bb * xi3;

  if(dphi == NULL) return;
  
  dphi[0] = svector<3>(-1,-1,-1);
  dphi[1] = svector<3>( 1, 0, 0);
  dphi[2] = svector<3>( 0, 1, 0);
  dphi[3] = svector<3>( 0, 0, 1);
  
  dphi[4][0] = xi[1]*xi[2]*xi3 - bb;
  dphi[4][1] = xi[2]*xi[0]*xi3 - bb;
  dphi[4][2] = xi[0]*xi[1]*xi3 - bb;

return;
}

// This buble has the Lagrange interpolation property

void tet4BubbleInterp(svector<3> xi, double phi[4] , svector<3> dphi[4])
{
  double xi3 = 1 - xi[0] - xi[1] - xi[2], bb = 64 * xi[0]*xi[1]*xi[2]*xi3;
  
  phi[0] = xi3   - bb;  
  phi[1] = xi[0] - bb;
  phi[2] = xi[1] - bb;
  phi[3] = xi[2] - bb;
  
  phi[4] = 4 * bb;

  if(dphi == NULL) return;
  
  svector<3> db(64*xi[1]*xi[2]*(xi3-xi[0]),
		64*xi[2]*xi[0]*(xi3-xi[1]),
		64*xi[0]*xi[1]*(xi3-xi[2]));

  dphi[0] = svector<3>(-1,-1,-1) - db;
  dphi[1] = svector<3>( 1, 0, 0) - db;
  dphi[2] = svector<3>( 0, 1, 0) - db;
  dphi[3] = svector<3>( 0, 0, 1) - db;
  
  dphi[4] = 4.0 * db;

return;
}

void tet10(const svector<3> xi, double phi[10], svector<3> dphi[10])  
{
double eta[4];

eta[0] = 1 - xi[0] - xi[1] - xi[2];   // baracentric coordiates

 if(dphi != NULL)
   {
     dphi[0] = svector<3>(-1,-1,-1);
     dphi[1] = svector<3>( 1, 0, 0);
     dphi[2] = svector<3>( 0, 1, 0);
     dphi[3] = svector<3>( 0, 0, 1);
   }

for(int i = 0; i < 3; i++)   // mid sides (0,1) (1,2) (2,0)
  {
  eta[i+1] = xi[i];   // baracentric coordiates

  phi[4+i] = 4 * eta[i] * eta[(i+1)%3];

  if(dphi != NULL)
    dphi[4+i] = (4 * eta[(i+1)%3]) * dphi[i] + (4 * eta[i]) * dphi[(i+1)%3];
  }
 
for(int i = 0; i < 3; i++)   // mid sides (0,3) (1,3) (2,3)
  {
  phi[7+i] = 4 * eta[i] * eta[3];

  if(dphi != NULL)
    dphi[7+i] = (4 * eta[3]) * dphi[i] + (4 * eta[i]) * dphi[3];
  }

for(int i = 0; i < 4; i++)   // vetriex basis functions
  {
  phi[i] = 2 * eta[i] * (eta[i] - 0.5);

  if(dphi != NULL)
    dphi[i] = (4 * eta[i] - 1) * dphi[i];
  }

return;
}

void tet20Permutation(const int *nodes, int *dof)
{
  int temp;

  // permute the edge dof

  if(nodes[0] > nodes[1]) {temp = dof[ 4]; dof[ 4] = dof[ 5]; dof[ 5] = temp;}
  if(nodes[1] > nodes[2]) {temp = dof[ 6]; dof[ 6] = dof[ 7]; dof[ 7] = temp;}
  if(nodes[2] > nodes[0]) {temp = dof[ 8]; dof[ 8] = dof[ 9]; dof[ 9] = temp;}

  if(nodes[3] > nodes[0]) {temp = dof[10]; dof[10] = dof[11]; dof[11] = temp;}
  if(nodes[3] > nodes[1]) {temp = dof[12]; dof[12] = dof[13]; dof[13] = temp;}
  if(nodes[3] > nodes[2]) {temp = dof[14]; dof[14] = dof[15]; dof[15] = temp;}

  return;
}

void tet20(const svector<3> pt, double psi[20], svector<3> dpsi[20])  
{

  double xi = pt[0], eta = pt[1], chi=pt[2];
  double  eta2 = eta*eta, xi2 = xi*xi, chi2 = chi*chi;
  double sqrt5 = sqrt(5.0);

psi[0] = -(xi-1+eta+chi)*(5*xi2+11*chi*xi-5*xi+11*xi*eta-5*eta+1+5*chi2-5*chi+11*eta*chi+5*eta2);
psi[1] = (5*xi2-chi*xi-5*xi-xi*eta+1-eta*chi+eta+chi-eta2-chi2)*xi;
psi[2] = -(-5*eta2+5*eta+eta*chi+xi*eta-1-xi-chi+xi2+chi*xi+chi2)*eta;
psi[3] = -(-5*chi2+5*chi+chi*xi+eta*chi-1-eta-xi+eta2+xi*eta+xi2)*chi;
psi[4] = 0.5*sqrt5*xi*(xi-1+eta+chi)*(10*xi-sqrt5+3*chi*sqrt5-5+5*chi+5*eta+3*eta*sqrt5);
psi[5] = -0.5*sqrt5*xi*(xi-1+eta+chi)*(10*xi-5-3*chi*sqrt5+sqrt5+5*chi+5*eta-3*eta*sqrt5);
psi[6] = 1.25*(sqrt5+3)*(2*xi-3+7*eta-3*eta*sqrt5+sqrt5)*xi*eta;
psi[7] = -1.25*(sqrt5-3)*(2*xi-3+3*eta*sqrt5+7*eta-sqrt5)*xi*eta;
psi[8] = -1.25*(sqrt5-3)*eta*(xi-1+eta+chi)*(2*xi-3*eta*sqrt5+2*chi-5*eta+sqrt5+1);
psi[9] = 1.25*(sqrt5+3)*eta*(xi-1+eta+chi)*(2*xi+2*chi-5*eta-sqrt5+3*eta*sqrt5+1);
psi[10] = -1.25*(sqrt5-3)*chi*(2*xi+1+2*eta+sqrt5-5*chi-3*chi*sqrt5)*(xi-1+eta+chi);
psi[11] = 1.25*(sqrt5+3)*chi*(xi-1+eta+chi)*(2*xi+3*chi*sqrt5-5*chi+2*eta-sqrt5+1);
psi[12] = -1.25*(sqrt5-3)*(2*xi-3+3*chi*sqrt5+7*chi-sqrt5)*xi*chi;
psi[13] = 1.25*(sqrt5+3)*(2*xi-3+7*chi-3*chi*sqrt5+sqrt5)*xi*chi;
psi[14] = -1.25*(sqrt5-3)*(2*eta-3+3*chi*sqrt5+7*chi-sqrt5)*eta*chi;
psi[15] = 1.25*(sqrt5+3)*(2*eta-3+7*chi-3*chi*sqrt5+sqrt5)*eta*chi;
psi[16] = 27*xi*eta*chi;
psi[17] = -27*(xi-1+eta+chi)*eta*chi;
psi[18] = -27*(xi-1+eta+chi)*xi*chi;
psi[19] = -27*(xi-1+eta+chi)*xi*eta;

 if(dpsi == NULL) return;

dpsi[0] = svector<3>(-15*xi2-32*chi*xi+20*xi-32*xi*eta+21*eta-6-16*chi2+21*chi-33*eta*chi-16*eta2,-16*xi2-33*chi*xi+21*xi-32*xi*eta+20*eta-6-16*chi2+21*chi-32*eta*chi-15*eta2,-16*xi2-32*chi*xi+21*xi-33*xi*eta+21*eta-6-15*chi2+20*chi-32*eta*chi-16*eta2);
dpsi[1] = svector<3>(15*xi2-2*chi*xi-10*xi-2*xi*eta+1-eta*chi+eta+chi-eta2-chi2,-(xi+chi-1+2*eta)*xi,-(xi+eta-1+2*chi)*xi);
dpsi[2] = svector<3>(-(eta-1+2*xi+chi)*eta,15*eta2-10*eta-2*eta*chi-2*xi*eta+1+xi+chi-xi2-chi*xi-chi2,-(xi+eta-1+2*chi)*eta);
dpsi[3] = svector<3>(-(eta-1+2*xi+chi)*chi,-(xi+chi-1+2*eta)*chi,15*chi2-10*chi-2*chi*xi-2*eta*chi+1+eta+xi-eta2-xi*eta-xi2);
dpsi[4] = svector<3>(0.5*sqrt5*(6*xi*eta*sqrt5+6*eta*chi*sqrt5+6*xi*chi*sqrt5-2*xi*sqrt5-4*eta*sqrt5-4*chi*sqrt5+3*eta2*sqrt5+3*chi2*sqrt5+sqrt5+5*eta2+5*chi2+30*xi2+5-10*chi-10*eta-30*xi+30*chi*xi+10*eta*chi+30*xi*eta),1.25*(sqrt5+1)*xi*(-3+2*chi-sqrt5+2*chi*sqrt5+2*eta*sqrt5+2*eta+6*xi),1.25*(sqrt5+1)*xi*(-3+2*chi-sqrt5+2*chi*sqrt5+2*eta*sqrt5+2*eta+6*xi));
dpsi[5] = svector<3>(-0.5*sqrt5*(-6*xi*eta*sqrt5-6*eta*chi*sqrt5-6*xi*chi*sqrt5+2*xi*sqrt5+4*eta*sqrt5+4*chi*sqrt5-3*eta2*sqrt5-3*chi2*sqrt5-sqrt5+5*eta2+5*chi2+30*xi2+5-10*chi-10*eta-30*xi+30*chi*xi+10*eta*chi+30*xi*eta),-1.25*(sqrt5-1)*xi*(6*xi-3-2*chi*sqrt5+sqrt5+2*chi+2*eta-2*eta*sqrt5),-1.25*(sqrt5-1)*xi*(6*xi-3-2*chi*sqrt5+sqrt5+2*chi+2*eta-2*eta*sqrt5));
dpsi[6] = svector<3>(1.25*(sqrt5+3)*(7*eta-3*eta*sqrt5-3+4*xi+sqrt5)*eta,1.25*(sqrt5+3)*(2*xi+sqrt5-3+14*eta-6*eta*sqrt5)*xi,0);
dpsi[7] = svector<3>(-1.25*(sqrt5-3)*(7*eta+3*eta*sqrt5-sqrt5+4*xi-3)*eta,-1.25*(sqrt5-3)*(2*xi-sqrt5-3+6*eta*sqrt5+14*eta)*xi,0);
dpsi[8] = svector<3>(-1.25*(sqrt5-3)*eta*(4*xi+4*chi+sqrt5-1-3*eta-3*eta*sqrt5),-1.25*(sqrt5-3)*(-6*xi*eta-6*eta*chi+4*chi*xi-xi+12*eta-chi+2*xi2-15*eta2+2*chi2-9*eta2*sqrt5-sqrt5+chi*sqrt5+8*eta*sqrt5+xi*sqrt5-1-6*eta*chi*sqrt5-6*xi*eta*sqrt5),-1.25*(sqrt5-3)*eta*(4*xi+4*chi+sqrt5-1-3*eta-3*eta*sqrt5));
dpsi[9] = svector<3>(1.25*(sqrt5+3)*eta*(4*xi+4*chi-3*eta+3*eta*sqrt5-sqrt5-1),1.25*(sqrt5+3)*(-6*xi*eta-6*eta*chi+4*chi*xi-xi+12*eta-chi+2*xi2-15*eta2+2*chi2+9*eta2*sqrt5+sqrt5-chi*sqrt5-8*eta*sqrt5-xi*sqrt5-1+6*eta*chi*sqrt5+6*xi*eta*sqrt5),1.25*(sqrt5+3)*eta*(4*xi+4*chi-3*eta+3*eta*sqrt5-sqrt5-1));
dpsi[10] = svector<3>(-1.25*(sqrt5-3)*chi*(4*xi-1+sqrt5+4*eta-3*chi-3*chi*sqrt5),-1.25*(sqrt5-3)*chi*(4*xi-1+sqrt5+4*eta-3*chi-3*chi*sqrt5),-1.25*(sqrt5-3)*(4*xi*eta-6*eta*chi-6*chi*xi-xi-eta+12*chi+2*xi2+2*eta2-15*chi2-9*chi2*sqrt5-sqrt5+8*chi*sqrt5+eta*sqrt5+xi*sqrt5-1-6*xi*chi*sqrt5-6*eta*chi*sqrt5));
dpsi[11] = svector<3>(1.25*(sqrt5+3)*chi*(4*xi-3*chi+3*chi*sqrt5+4*eta-sqrt5-1),1.25*(sqrt5+3)*chi*(4*xi-3*chi+3*chi*sqrt5+4*eta-sqrt5-1),1.25*(sqrt5+3)*(4*xi*eta-6*eta*chi-6*chi*xi-xi-eta+12*chi+2*xi2+2*eta2-15*chi2+9*chi2*sqrt5+sqrt5-8*chi*sqrt5-eta*sqrt5-xi*sqrt5-1+6*xi*chi*sqrt5+6*eta*chi*sqrt5));
dpsi[12] = svector<3>(-1.25*(sqrt5-3)*(7*chi+3*chi*sqrt5-sqrt5+4*xi-3)*chi,0,-1.25*(sqrt5-3)*(2*xi-sqrt5-3+6*chi*sqrt5+14*chi)*xi);
dpsi[13] = svector<3>(1.25*(sqrt5+3)*(7*chi-3*chi*sqrt5-3+4*xi+sqrt5)*chi,0,1.25*(sqrt5+3)*(2*xi+sqrt5-3+14*chi-6*chi*sqrt5)*xi);
dpsi[14] = svector<3>(0,-1.25*(sqrt5-3)*(7*chi+3*chi*sqrt5-sqrt5+4*eta-3)*chi,-1.25*(sqrt5-3)*(2*eta-sqrt5-3+6*chi*sqrt5+14*chi)*eta);
dpsi[15] = svector<3>(0,1.25*(sqrt5+3)*(-3*chi*sqrt5+7*chi-3+4*eta+sqrt5)*chi,1.25*(sqrt5+3)*(2*eta+sqrt5-3+14*chi-6*chi*sqrt5)*eta);
dpsi[16] = svector<3>(27*eta*chi,27*chi*xi,27*xi*eta);
dpsi[17] = svector<3>(-27*eta*chi,-27*(xi+chi-1+2*eta)*chi,-27*(xi+eta-1+2*chi)*eta);
dpsi[18] = svector<3>(-27*(eta-1+2*xi+chi)*chi,-27*chi*xi,-27*(xi+eta-1+2*chi)*xi);
dpsi[19] = svector<3>(-27*(eta-1+2*xi+chi)*eta,-27*(xi+chi-1+2*eta)*xi,-27*xi*eta);

  return;
}


void cube8(svector<3> xi, double psi[8], svector<3> dpsi[8])
{
double eta = xi[0], ata = xi[1], chi = xi[2];

double ope = (1+eta)/2, opa = (1+ata)/2, opc = (1+chi)/2;
double ome = (1-eta)/2, oma = (1-ata)/2, omc = (1-chi)/2;

 psi[0] = ome * oma * omc;
 psi[1] = ope * oma * omc;
 psi[2] = ome * opa * omc;
 psi[3] = ope * opa * omc;
 psi[4] = ome * oma * opc;
 psi[5] = ope * oma * opc;
 psi[6] = ome * opa * opc;
 psi[7] = ope * opa * opc;

 if(dpsi == NULL) return;

 dpsi[0] = -0.5 * svector<3>( oma * omc,  ome * omc,  ome * oma);
 dpsi[1] =  0.5 * svector<3>( oma * omc, -ope * omc, -ope * oma);
 dpsi[2] =  0.5 * svector<3>(-opa * omc,  ome * omc, -ome * opa);
 dpsi[3] =  0.5 * svector<3>( opa * omc,  ope * omc, -ope * opa);
 dpsi[4] =  0.5 * svector<3>(-oma * opc, -ome * opc,  ome * oma);
 dpsi[5] =  0.5 * svector<3>( oma * opc, -ope * opc,  ope * oma);
 dpsi[6] =  0.5 * svector<3>(-opa * opc,  ome * opc,  ome * opa);
 dpsi[7] =  0.5 * svector<3>( opa * opc,  ope * opc,  ope * opa);

 return;
}

void cube8Bubble(svector<3> xi, double psi[9], svector<3> dpsi[9])
{
  double eta = xi[0], ata = xi[1], chi = xi[2];

  double ope = (1+eta)/2, opa = (1+ata)/2, opc = (1+chi)/2;
  double ome = (1-eta)/2, oma = (1-ata)/2, omc = (1-chi)/2;

  double bb = (1-eta*eta)*(1-ata*ata)*(1-chi*chi);
  double bw = 1 + (eta + ata + chi) + (eta*ata + ata*chi + chi*eta);

 psi[0] = ome * oma * omc;
 psi[1] = ope * oma * omc;
 psi[2] = ome * opa * omc;
 psi[3] = ope * opa * omc;
 psi[4] = ome * oma * opc;
 psi[5] = ope * oma * opc;
 psi[6] = ome * opa * opc;
 psi[7] = ope * opa * opc;

 psi[8] = bb*bw;

 if(dpsi == NULL) return;

 dpsi[0] = -0.5 * svector<3>( oma * omc,  ome * omc,  ome * oma);
 dpsi[1] =  0.5 * svector<3>( oma * omc, -ope * omc, -ope * oma);
 dpsi[2] =  0.5 * svector<3>(-opa * omc,  ome * omc, -ome * opa);
 dpsi[3] =  0.5 * svector<3>( opa * omc,  ope * omc, -ope * opa);
 dpsi[4] =  0.5 * svector<3>(-oma * opc, -ome * opc,  ome * oma);
 dpsi[5] =  0.5 * svector<3>( oma * opc, -ope * opc,  ope * oma);
 dpsi[6] =  0.5 * svector<3>(-opa * opc,  ome * opc,  ome * opa);
 dpsi[7] =  0.5 * svector<3>( opa * opc,  ope * opc,  ope * opa);

 dpsi[8][0] = bb*(1+ata+chi) - 2*eta*(1-ata*ata)*(1-chi*chi)*bw;
 dpsi[8][1] = bb*(1+chi+eta) - 2*ata*(1-chi*chi)*(1-eta*eta)*bw;
 dpsi[8][2] = bb*(1+eta+ata) - 2*chi*(1-eta*eta)*(1-ata*ata)*bw;

 return;
}

void cube20(svector<3> xi, double psi[20], svector<3> dpsi[20])
{
  double eta = xi[0], ata = xi[1], chi = xi[2];
  double ope = (1+eta)/2, opa = (1+ata)/2, opc = (1+chi)/2;
  double ome = (1-eta)/2, oma = (1-ata)/2, omc = (1-chi)/2;
  double ome2= 1-eta*eta, oma2= 1-ata*ata, omc2= 1-chi*chi;
  
  psi[0] = ome * oma * omc * ((-eta - ata - chi) - 2);
  psi[1] = ope * oma * omc * (( eta - ata - chi) - 2);
  psi[2] = ome * opa * omc * ((-eta + ata - chi) - 2);
  psi[3] = ope * opa * omc * (( eta + ata - chi) - 2);
  psi[4] = ome * oma * opc * ((-eta - ata + chi) - 2);
  psi[5] = ope * oma * opc * (( eta - ata + chi) - 2);
  psi[6] = ome * opa * opc * ((-eta + ata + chi) - 2);
  psi[7] = ope * opa * opc * (( eta + ata + chi) - 2);
  
  psi[8] = ome2 * oma * omc;
  psi[9] = ome2 * opa * omc;
  psi[10]= ome2 * oma * opc;
  psi[11]= ome2 * opa * opc;
  
  psi[12]= ome * oma2 * omc;
  psi[13]= ope * oma2 * omc;
  psi[14]= ome * oma2 * opc;
  psi[15]= ope * oma2 * opc;
  
  psi[16]= ome * oma * omc2;
  psi[17]= ope * oma * omc2;
  psi[18]= ome * opa * omc2;
  psi[19]= ope * opa * omc2;

  if(dpsi == NULL) return;
  
  dpsi[0] = svector<3>(oma * omc * ( chi + ata + 2*eta + 1)/ 2,
		       ome * omc * ( chi + 2*ata + eta + 1)/ 2, 
		       ome * oma * ( 2*chi + ata + eta + 1)/ 2);

  dpsi[1] = svector<3>(oma * omc * ( chi + ata - 2*eta + 1)/-2,
		       ope * omc * ( chi + 2*ata - eta + 1)/ 2, 
		       ope * oma * ( 2*chi + ata - eta + 1)/ 2);
  
  dpsi[2] = svector<3>(opa * omc * ( chi - ata + 2*eta + 1)/ 2, 
		       ome * omc * ( chi - 2*ata + eta + 1)/-2, 
		       ome * opa * ( 2*chi - ata + eta + 1)/ 2);
  
  dpsi[3] = svector<3>(opa * omc * ( chi - ata - 2*eta + 1)/-2, 
		       ope * omc * ( chi - 2*ata - eta + 1)/-2, 
		       ope * opa * ( 2*chi - ata - eta + 1)/ 2);
  
  dpsi[4] = svector<3>(oma * opc * (-chi + ata + 2*eta + 1)/ 2, 
		       ome * opc * (-chi + 2*ata + eta + 1)/ 2, 
		       ome * oma * (-2*chi + ata + eta + 1)/-2);
  
  dpsi[5] = svector<3>(oma * opc * (-chi + ata - 2*eta + 1)/-2, 
		       ope * opc * (-chi + 2*ata - eta + 1)/ 2, 
		       ope * oma * (-2*chi + ata - eta + 1)/-2);
  
  dpsi[6] = svector<3>(opa * opc * (-chi - ata + 2*eta + 1)/ 2, 
		       ome * opc * (-chi - 2*ata + eta + 1)/-2, 
		       ome * opa * (-2*chi - ata + eta + 1)/-2);
  
  dpsi[7] = svector<3>(opa * opc * (-chi - ata - 2*eta + 1)/-2, 
		       ope * opc * (-chi - 2*ata - eta + 1)/-2, 
		       ope * opa * (-2*chi - ata - eta + 1)/-2);
  
  dpsi[8] = svector<3>(-2*eta*oma*omc, ome2*omc/-2, ome2*oma/-2);
  dpsi[9] = svector<3>(-2*eta*opa*omc, ome2*omc/ 2, ome2*opa/-2);
  dpsi[10] = svector<3>(-2*eta*oma*opc, ome2*opc/-2, ome2*oma/ 2);
  dpsi[11] = svector<3>(-2*eta*opa*opc, ome2*opc/ 2, ome2*opa/ 2);
  
  dpsi[12] = svector<3>(oma2*omc/-2, ome*-2*ata*omc, ome*oma2/-2);
  dpsi[13] = svector<3>(oma2*omc/ 2, ope*-2*ata*omc, ope*oma2/-2);
  dpsi[14] = svector<3>(oma2*opc/-2, ome*-2*ata*opc, ome*oma2/ 2);
  dpsi[15] = svector<3>(oma2*opc/ 2, ope*-2*ata*opc, ope*oma2/ 2);
  
  dpsi[16] = svector<3>(oma*omc2/-2, ome*omc2/-2, ome*oma*-2*chi);
  dpsi[17] = svector<3>(oma*omc2/ 2, ope*omc2/-2, ope*oma*-2*chi);
  dpsi[18] = svector<3>(opa*omc2/-2, ome*omc2/ 2, ome*opa*-2*chi);
  dpsi[19] = svector<3>(opa*omc2/ 2, ope*omc2/ 2, ope*opa*-2*chi);
  
  return;
}

void cube27(svector<3> xi, double phi[27], svector<3> dphi[27])
{
int i,j,k, ii = 0;
double onedeta[3], donedeta[3];
double onedata[3], donedata[3];
double onedchi[3], donedchi[3];

svector<1> vdonedeta[3], vdonedata[3], vdonedchi[3];

 svector<1> eta(xi[0]), ata(xi[1]), chi(xi[2]);

line3(eta,onedeta,vdonedeta);
line3(ata,onedata,vdonedata);
line3(chi,onedchi,vdonedchi);

for(int i = 0; i < 3; i++)
  {
    donedeta[i] = vdonedeta[i][0];
    donedata[i] = vdonedata[i][0];
    donedchi[i] = vdonedchi[i][0];
  } 

/* Nodes */

for(k = 0; k < 2; k ++)
for(j = 0; j < 2; j ++)
for(i = 0; i < 2; i ++)
  {
  phi[ii] = onedeta[i] * onedata[j] * onedchi[k];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[i] *  onedata[j] *  onedchi[k],
			  onedeta[i] * donedata[j] *  onedchi[k],
			  onedeta[i] *  onedata[j] * donedchi[k]);

  ii++;
  }

/* Mid Edges */

for(k = 0; k < 2; k ++)
for(j = 0; j < 2; j ++)
  {
  phi[ii] = onedeta[2] * onedata[j] * onedchi[k];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[2] *  onedata[j] *  onedchi[k],
			  onedeta[2] * donedata[j] *  onedchi[k],
			  onedeta[2] *  onedata[j] * donedchi[k]);
  
  ii++;
  }

for(k = 0; k < 2; k ++)
for(i = 0; i < 2; i ++)
  {
  phi[ii] = onedeta[i] * onedata[2] * onedchi[k];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[i] *  onedata[2] *  onedchi[k],
			  onedeta[i] * donedata[2] *  onedchi[k],
			  onedeta[i] *  onedata[2] * donedchi[k]);
  
  ii++;
  }

for(j = 0; j < 2; j ++)
for(i = 0; i < 2; i ++)
  {
  phi[ii] = onedeta[i] * onedata[j] * onedchi[2];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[i] *  onedata[j] *  onedchi[2],
			  onedeta[i] * donedata[j] *  onedchi[2],
			  onedeta[i] *  onedata[j] * donedchi[2]);
  
  ii++;
  }

/* Mid Faces */

for(k = 0; k < 2; k ++)
  {
  phi[ii] = onedeta[2] * onedata[2] * onedchi[k];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[2] *  onedata[2] *  onedchi[k],
			  onedeta[2] * donedata[2] *  onedchi[k],
			  onedeta[2] *  onedata[2] * donedchi[k]);
  
  ii++;
  }

for(j = 0; j < 2; j ++)
  {
  phi[ii] = onedeta[2] * onedata[j] * onedchi[2];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[2] *  onedata[j] *  onedchi[2],
			  onedeta[2] * donedata[j] *  onedchi[2],
			  onedeta[2] *  onedata[j] * donedchi[2]);

  ii++;
  }

for(i = 0; i < 2; i ++)
  {
  phi[ii] = onedeta[i] * onedata[2] * onedchi[2];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[i] *  onedata[2] *  onedchi[2],
			  onedeta[i] * donedata[2] *  onedchi[2],
			  onedeta[i] *  onedata[2] * donedchi[2]);
  
  ii++;
  }

/* Mid Cell */

  phi[ii] = onedeta[2] * onedata[2] * onedchi[2];

  if(dphi != NULL)
    dphi[ii] = svector<3>(donedeta[2] *  onedata[2] *  onedchi[2],
			  onedeta[2] * donedata[2] *  onedchi[2],
			  onedeta[2] *  onedata[2] * donedchi[2]);
  return;
}

void cube32(svector<3> xi, double phi[32], svector<3> dphi[32])
{
  // Serindipity cubic cube

  double x = xi[0], y = xi[1], z = xi[2];
  double x2 = x*x,  y2 = y*y,  z2 = z*z, sqrt5 = sqrt(5.0);

  // vertex basis functions

  phi[0] = -(1.0/32)*(-1+z)*(-1+y)*(x-1)*(5*x2+5*y2+5*z2-11);
  phi[1] = (1.0/32)*(-1+z)*(-1+y)*(x+1)*(5*x2+5*y2+5*z2-11);
  phi[2] = (1.0/32)*(-1+z)*(y+1)*(x-1)*(5*x2+5*y2+5*z2-11);
  phi[3] = -(1.0/32)*(-1+z)*(y+1)*(x+1)*(5*x2+5*y2+5*z2-11);
  phi[4] = (1.0/32)*(1+z)*(-1+y)*(x-1)*(5*x2+5*y2+5*z2-11);
  phi[5] = -(1.0/32)*(1+z)*(-1+y)*(x+1)*(5*x2+5*y2+5*z2-11);
  phi[6] = -(1.0/32)*(1+z)*(y+1)*(x-1)*(5*x2+5*y2+5*z2-11);
  phi[7] = (1.0/32)*(1+z)*(y+1)*(x+1)*(5*x2+5*y2+5*z2-11);

  // edge basis functions

  phi[8] = -(1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(-1+y)*(-1+z);
  phi[9] = -(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(-1+y)*(-1+z);
  phi[10] = (1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(y+1)*(-1+z);
  phi[11] = (1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(y+1)*(-1+z);
  phi[12] = (1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(-1+y)*(1+z);
  phi[13] = (1.0/32)*sqrt5*(-5*x+sqrt5)*(x-1)*(x+1)*(-1+y)*(1+z);
  phi[14] = -(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(y+1)*(1+z);
  phi[15] = -(1.0/32)*sqrt5*(x-1)*(x+1)*(-5*x+sqrt5)*(y+1)*(1+z);
  phi[16] = -(1.0/32)*sqrt5*(x-1)*(-1+y)*(5*y+sqrt5)*(y+1)*(-1+z);
  phi[17] = -(1.0/32)*sqrt5*(x-1)*(-1+y)*(-5*y+sqrt5)*(y+1)*(-1+z);
  phi[18] = (1.0/32)*sqrt5*(x+1)*(-1+y)*(y+1)*(-5*y+sqrt5)*(-1+z);
  phi[19] = (1.0/32)*sqrt5*(x+1)*(-1+y)*(5*y+sqrt5)*(y+1)*(-1+z);
  phi[20] = (1.0/32)*sqrt5*(x-1)*(-1+y)*(5*y+sqrt5)*(y+1)*(1+z);
  phi[21] = (1.0/32)*sqrt5*(x-1)*(-1+y)*(y+1)*(-5*y+sqrt5)*(1+z);
  phi[22] = -(1.0/32)*sqrt5*(x+1)*(-1+y)*(y+1)*(-5*y+sqrt5)*(1+z);
  phi[23] = -(1.0/32)*sqrt5*(x+1)*(-1+y)*(y+1)*(5*y+sqrt5)*(1+z);
  phi[24] = -(1.0/32)*sqrt5*(x-1)*(-1+y)*(-1+z)*(1+z)*(-5*z+sqrt5);
  phi[25] = -(1.0/32)*sqrt5*(x-1)*(-1+y)*(-1+z)*(5*z+sqrt5)*(1+z);
  phi[26] = (1.0/32)*sqrt5*(x+1)*(-1+y)*(-1+z)*(1+z)*(-5*z+sqrt5);
  phi[27] = (1.0/32)*sqrt5*(x+1)*(-1+y)*(-1+z)*(5*z+sqrt5)*(1+z);
  phi[28] = (1.0/32)*sqrt5*(x-1)*(y+1)*(-5*z+sqrt5)*(-1+z)*(1+z);
  phi[29] = (1.0/32)*sqrt5*(x-1)*(y+1)*(5*z+sqrt5)*(-1+z)*(1+z);
  phi[30] = -(1.0/32)*sqrt5*(x+1)*(y+1)*(-5*z+sqrt5)*(-1+z)*(1+z);
  phi[31] = -(1.0/32)*sqrt5*(x+1)*(y+1)*(-1+z)*(5*z+sqrt5)*(1+z);

  if(dphi == NULL) return;

  dphi[0] = svector<3>(-(1.0/32)*(-1+z)*(-1+y)*(5*y2+5*z2+15*x2-11-10*x),
		       -(1.0/32)*(-1+z)*(x-1)*(5*x2+5*z2-10*y+15*y2-11),
		       -(1.0/32)*(-1+y)*(x-1)*(5*x2+5*y2+15*z2-10*z-11));
  dphi[1] = svector<3>((1.0/32)*(-1+z)*(-1+y)*(5*y2+5*z2+15*x2-11+10*x), 
		       (1.0/32)*(-1+z)*(x+1)*(5*x2+5*z2-10*y+15*y2-11), 
		       (1.0/32)*(-1+y)*(x+1)*(5*x2+5*y2+15*z2-10*z-11));
  dphi[2] = svector<3>((1.0/32)*(-1+z)*(y+1)*(5*y2+5*z2+15*x2-11-10*x),
		       (1.0/32)*(-1+z)*(x-1)*(5*x2+5*z2+10*y+15*y2-11),
		       (1.0/32)*(y+1)*(x-1)*(5*x2+5*y2+15*z2-10*z-11));
  dphi[3] = svector<3>(-(1.0/32)*(-1+z)*(y+1)*(5*y2+5*z2+15*x2-11+10*x),
		       -(1.0/32)*(-1+z)*(x+1)*(5*x2+5*z2+10*y+15*y2-11),
		       -(1.0/32)*(y+1)*(x+1)*(5*x2+5*y2+15*z2-10*z-11));
  dphi[4] = svector<3>((1.0/32)*(1+z)*(-1+y)*(5*y2+5*z2+15*x2-11-10*x),
		       (1.0/32)*(1+z)*(x-1)*(5*x2+5*z2-10*y+15*y2-11),
		       (1.0/32)*(-1+y)*(x-1)*(5*x2+5*y2+15*z2+10*z-11));
  dphi[5] = svector<3>(-(1.0/32)*(1+z)*(-1+y)*(5*y2+5*z2+15*x2-11+10*x),
		       -(1.0/32)*(1+z)*(x+1)*(5*x2+5*z2-10*y+15*y2-11),
		       -(1.0/32)*(-1+y)*(x+1)*(5*x2+5*y2+15*z2+10*z-11));
  dphi[6] = svector<3>(-(1.0/32)*(1+z)*(y+1)*(5*y2+5*z2+15*x2-11-10*x),
		       -(1.0/32)*(1+z)*(x-1)*(5*x2+5*z2+10*y+15*y2-11),
		       -(1.0/32)*(y+1)*(x-1)*(5*x2+5*y2+15*z2+10*z-11));
  dphi[7] = svector<3>((1.0/32)*(1+z)*(y+1)*(5*y2+5*z2+15*x2-11+10*x),
		       (1.0/32)*(1+z)*(x+1)*(5*x2+5*z2+10*y+15*y2-11),
		       (1.0/32)*(y+1)*(x+1)*(5*x2+5*y2+15*z2+10*z-11));
  dphi[8] = svector<3>(-(1.0/32)*sqrt5*(-3*x+sqrt5)*(5*x+sqrt5)*(-1+y)*(-1+z),
		       -(1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(-1+z),
		       -(1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(-1+y));
  dphi[9] = svector<3>((1.0/32)*sqrt5*(3*x+sqrt5)*(-5*x+sqrt5)*(-1+y)*(-1+z),
		       -(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(-1+z),
		       -(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(-1+y));
  dphi[10] = svector<3>(-(1.0/32)*sqrt5*(3*x+sqrt5)*(-5*x+sqrt5)*(y+1)*(-1+z),
			(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(-1+z),
			(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(y+1));
  dphi[11] = svector<3>((1.0/32)*sqrt5*(-3*x+sqrt5)*(5*x+sqrt5)*(y+1)*(-1+z),
			(1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(-1+z),
			(1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(y+1));
  dphi[12] = svector<3>(-(1.0/32)*sqrt5*(3*x+sqrt5)*(-5*x+sqrt5)*(-1+y)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(-1+y));
  dphi[13] = svector<3>((1.0/32)*sqrt5*(-3*x+sqrt5)*(5*x+sqrt5)*(-1+y)*(1+z),
			(1.0/32)*sqrt5*(-5*x+sqrt5)*(x-1)*(x+1)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(-1+y));
  dphi[14] = svector<3>((1.0/32)*sqrt5*(3*x+sqrt5)*(-5*x+sqrt5)*(y+1)*(1+z),
			-(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(1+z),
			-(1.0/32)*sqrt5*(x-1)*(x+1)*(5*x+sqrt5)*(y+1));
  dphi[15] = svector<3>(-(1.0/32)*sqrt5*(-3*x+sqrt5)*(5*x+sqrt5)*(y+1)*(1+z),
			-(1.0/32)*sqrt5*(-5*x+sqrt5)*(x-1)*(x+1)*(1+z),
			-(1.0/32)*sqrt5*(x-1)*(-5*x+sqrt5)*(x+1)*(y+1));
  dphi[16] = svector<3>(-(1.0/32)*sqrt5*(-1+y)*(5*y+sqrt5)*(y+1)*(-1+z),
			(1.0/32)*sqrt5*(x-1)*(-5*y+sqrt5)*(3*y+sqrt5)*(-1+z),
			-(1.0/32)*sqrt5*(x-1)*(-1+y)*(5*y+sqrt5)*(y+1));
  dphi[17] = svector<3>(-(1.0/32)*sqrt5*(-1+y)*(-5*y+sqrt5)*(y+1)*(-1+z),
			-(1.0/32)*sqrt5*(x-1)*(5*y+sqrt5)*(-3*y+sqrt5)*(-1+z),
			-(1.0/32)*sqrt5*(x-1)*(-1+y)*(-5*y+sqrt5)*(y+1));
  dphi[18] = svector<3>((1.0/32)*sqrt5*(-1+y)*(-5*y+sqrt5)*(y+1)*(-1+z),
			(1.0/32)*sqrt5*(x+1)*(5*y+sqrt5)*(-3*y+sqrt5)*(-1+z),
			(1.0/32)*sqrt5*(x+1)*(-1+y)*(y+1)*(-5*y+sqrt5));
  dphi[19] = svector<3>((1.0/32)*sqrt5*(-1+y)*(5*y+sqrt5)*(y+1)*(-1+z),
			-(1.0/32)*sqrt5*(x+1)*(-5*y+sqrt5)*(3*y+sqrt5)*(-1+z),
			(1.0/32)*sqrt5*(x+1)*(-1+y)*(5*y+sqrt5)*(y+1));
  dphi[20] = svector<3>((1.0/32)*sqrt5*(-1+y)*(5*y+sqrt5)*(y+1)*(1+z),
			-(1.0/32)*sqrt5*(x-1)*(-5*y+sqrt5)*(3*y+sqrt5)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(-1+y)*(5*y+sqrt5)*(y+1));
  dphi[21] = svector<3>((1.0/32)*sqrt5*(-1+y)*(y+1)*(-5*y+sqrt5)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(5*y+sqrt5)*(-3*y+sqrt5)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(-1+y)*(-5*y+sqrt5)*(y+1));
  dphi[22] = svector<3>(-(1.0/32)*sqrt5*(-1+y)*(y+1)*(-5*y+sqrt5)*(1+z),
			-(1.0/32)*sqrt5*(x+1)*(5*y+sqrt5)*(-3*y+sqrt5)*(1+z),
			-(1.0/32)*sqrt5*(x+1)*(-1+y)*(y+1)*(-5*y+sqrt5));
  dphi[23] = svector<3>(-(1.0/32)*sqrt5*(-1+y)*(5*y+sqrt5)*(y+1)*(1+z),
			(1.0/32)*sqrt5*(x+1)*(-5*y+sqrt5)*(3*y+sqrt5)*(1+z),
			-(1.0/32)*sqrt5*(x+1)*(-1+y)*(5*y+sqrt5)*(y+1));
  dphi[24] = svector<3>(-(1.0/32)*sqrt5*(-1+y)*(-1+z)*(1+z)*(-5*z+sqrt5),
			-(1.0/32)*sqrt5*(x-1)*(-1+z)*(1+z)*(-5*z+sqrt5),
			-(1.0/32)*sqrt5*(x-1)*(-1+y)*(-3*z+sqrt5)*(5*z+sqrt5));
  dphi[25] = svector<3>(-(1.0/32)*sqrt5*(-1+y)*(-1+z)*(5*z+sqrt5)*(1+z),
			-(1.0/32)*sqrt5*(x-1)*(-1+z)*(5*z+sqrt5)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(-1+y)*(3*z+sqrt5)*(-5*z+sqrt5));
  dphi[26] = svector<3>((1.0/32)*sqrt5*(-1+y)*(-1+z)*(1+z)*(-5*z+sqrt5),
			(1.0/32)*sqrt5*(x+1)*(-1+z)*(1+z)*(-5*z+sqrt5),
			(1.0/32)*sqrt5*(x+1)*(-1+y)*(-3*z+sqrt5)*(5*z+sqrt5));
  dphi[27] = svector<3>((1.0/32)*sqrt5*(-1+y)*(-1+z)*(5*z+sqrt5)*(1+z),
			(1.0/32)*sqrt5*(x+1)*(-1+z)*(5*z+sqrt5)*(1+z),
			-(1.0/32)*sqrt5*(x+1)*(-1+y)*(3*z+sqrt5)*(-5*z+sqrt5));
  dphi[28] = svector<3>((1.0/32)*sqrt5*(y+1)*(-5*z+sqrt5)*(-1+z)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(-1+z)*(1+z)*(-5*z+sqrt5),
			(1.0/32)*sqrt5*(x-1)*(y+1)*(-3*z+sqrt5)*(5*z+sqrt5));
  dphi[29] = svector<3>((1.0/32)*sqrt5*(y+1)*(5*z+sqrt5)*(-1+z)*(1+z),
			(1.0/32)*sqrt5*(x-1)*(-1+z)*(5*z+sqrt5)*(1+z),
			-(1.0/32)*sqrt5*(x-1)*(y+1)*(3*z+sqrt5)*(-5*z+sqrt5));
  dphi[30] = svector<3>(-(1.0/32)*sqrt5*(y+1)*(-5*z+sqrt5)*(-1+z)*(1+z),
			-(1.0/32)*sqrt5*(x+1)*(-1+z)*(1+z)*(-5*z+sqrt5),
			-(1.0/32)*sqrt5*(x+1)*(y+1)*(-3*z+sqrt5)*(5*z+sqrt5));
  dphi[31] = svector<3>(-(1.0/32)*sqrt5*(y+1)*(5*z+sqrt5)*(-1+z)*(1+z),
			-(1.0/32)*sqrt5*(x+1)*(-1+z)*(5*z+sqrt5)*(1+z),
			(1.0/32)*sqrt5*(x+1)*(y+1)*(3*z+sqrt5)*(-5*z+sqrt5));
  
  return;
}

void cube32Permutation(const int *nodes, int *dof)
{
  int temp;

  int ee[12][2] = {{0,1}, {3,2}, {5,4}, {7,6}, {2,0}, {1,3}, 
		   {6,4}, {5,7}, {0,4}, {1,5}, {2,6}, {3,7}};

  // permute the edge dof

  for(int i = 0; i < 12; i++)
    {
      if(nodes[ee[i][0]] > nodes[ee[i][1]])  // swap the edge dof
	{
	  temp = dof[8 + 2*i]; 
	  dof[8 + 2*i] = dof[9 + 2*i];
	  dof[9 + 2*i] = temp;
	}
    }

  return;
}

void cube50(svector<3> xi, double phi[50], svector<3> dphi[50])
{
  // Serindipity cubic cube

  double x = xi[0],  y = xi[1], z = xi[2], sqrt21 = sqrt(21.0);
  double x2 = x*x,   y2 = y*y,   z2 = z*z; 
  double x3 = x*x2,  y3 = y*y2,  z3 = z*z2;

  // vertex basis functions

  phi[0] = (1.0/32)*(z-1)*(y-1)*(-1+x)*(7*x3-4*x*z-4*x*y-11*x+7*y3-4*y*z-11*y+7*z3-11*z-4);
  phi[1] = (1.0/32)*(z-1)*(y-1)*(1+x)*(7*x3-4*x*y-4*x*z-11*x-7*y3+4*y*z+11*y-7*z3+11*z+4);
  phi[2] = -(1.0/32)*(z-1)*(1+y)*(-1+x)*(7*x3+4*x*y-4*x*z-11*x-7*y3+4*y*z+11*y+7*z3-11*z-4);
  phi[3] = -(1.0/32)*(z-1)*(1+y)*(1+x)*(7*x3+4*x*y-11*x-4*x*z+4+7*y3-4*y*z-7*z3-11*y+11*z);
  phi[4] = -(1.0/32)*(1+z)*(y-1)*(-1+x)*(7*x3-4*x*y+4*x*z-11*x+7*y3+4*y*z-11*y-7*z3+11*z-4);
  phi[5] = -(1.0/32)*(1+z)*(y-1)*(1+x)*(7*x3+4*x*z-4*x*y-11*x-7*y3+11*y-4*y*z+7*z3+4-11*z);
  phi[6] = (1.0/32)*(1+z)*(1+y)*(-1+x)*(7*x3+4*x*y+4*x*z-11*x-7*y3-4*y*z+11*y-7*z3+11*z-4);
  phi[7] = (1.0/32)*(1+z)*(1+y)*(1+x)*(7*x3+4*x*y+4*x*z-11*x+7*y3+4*y*z-11*y+7*z3-11*z+4);

  // edge basis functions

  phi[8] = (7.0/96)*(-7*x+sqrt21)*(1+x)*x*(-1+x)*(y-1)*(z-1);
  phi[9] = (1.0/12)*(z-1)*(y-1)*(-1+x)*(1+x)*(7*x2+3*y+3*z+3);
  phi[10] = -(7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(y-1)*(z-1);
  phi[11] = (7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(1+y)*(z-1);
  phi[12] = -(1.0/12)*(z-1)*(1+y)*(-1+x)*(1+x)*(7*x2-3*y+3*z+3);
  phi[13] = -(7.0/96)*(1+x)*x*(-7*x+sqrt21)*(-1+x)*(1+y)*(z-1);
  phi[14] = (7.0/96)*(1+x)*x*(-1+x)*(7*x+sqrt21)*(y-1)*(1+z);
  phi[15] = -(1.0/12)*(1+z)*(y-1)*(-1+x)*(1+x)*(7*x2+3*y-3*z+3);
  phi[16] = -(7.0/96)*(-7*x+sqrt21)*(1+x)*x*(-1+x)*(y-1)*(1+z);
  phi[17] = -(7.0/96)*(1+x)*(7*x+sqrt21)*x*(-1+x)*(1+y)*(1+z);
  phi[18] = (1.0/12)*(1+z)*(1+y)*(-1+x)*(1+x)*(7*x2-3*y-3*z+3);
  phi[19] = (7.0/96)*(1+x)*(-7*x+sqrt21)*x*(-1+x)*(1+y)*(1+z);
  phi[20] = -(7.0/96)*(-1+x)*y*(y-1)*(7*y+sqrt21)*(1+y)*(z-1);
  phi[21] = (1.0/12)*(y-1)*(1+y)*(z-1)*(-1+x)*(3*x+7*y2+3*z+3);
  phi[22] = (7.0/96)*(-1+x)*y*(-7*y+sqrt21)*(y-1)*(1+y)*(z-1);
  phi[23] = -(7.0/96)*(1+x)*y*(y-1)*(-7*y+sqrt21)*(1+y)*(z-1);
  phi[24] = (1.0/12)*(y-1)*(1+y)*(z-1)*(1+x)*(3*x-7*y2-3*z-3);
  phi[25] = (7.0/96)*(1+x)*y*(7*y+sqrt21)*(y-1)*(1+y)*(z-1);
  phi[26] = (7.0/96)*(-1+x)*y*(7*y+sqrt21)*(y-1)*(1+y)*(1+z);
  phi[27] = -(1.0/12)*(y-1)*(1+y)*(1+z)*(-1+x)*(3*x+7*y2-3*z+3);
  phi[28] = -(7.0/96)*(-1+x)*y*(-7*y+sqrt21)*(y-1)*(1+y)*(1+z);
  phi[29] = (7.0/96)*(1+x)*y*(-7*y+sqrt21)*(y-1)*(1+y)*(1+z);
  phi[30] = -(1.0/12)*(y-1)*(1+y)*(1+z)*(1+x)*(3*x-7*y2+3*z-3);
  phi[31] = -(7.0/96)*(1+x)*y*(7*y+sqrt21)*(y-1)*(1+y)*(1+z);
  phi[32] = (7.0/96)*(-1+x)*(y-1)*(1+z)*(z-1)*(-7*z+sqrt21)*z;
  phi[33] = (1.0/12)*(z-1)*(1+z)*(y-1)*(-1+x)*(3*x+3*y+7*z2+3);
  phi[34] = -(7.0/96)*(-1+x)*(y-1)*(7*z+sqrt21)*(1+z)*(z-1)*z;
  phi[35] = -(7.0/96)*(1+x)*(y-1)*(1+z)*(z-1)*(-7*z+sqrt21)*z;
  phi[36] = (1.0/12)*(z-1)*(1+z)*(y-1)*(1+x)*(3*x-3*y-7*z2-3);
  phi[37] = (7.0/96)*(1+x)*(y-1)*(1+z)*(z-1)*(7*z+sqrt21)*z;
  phi[38] = -(7.0/96)*(-1+x)*(1+y)*(1+z)*(z-1)*(-7*z+sqrt21)*z;
  phi[39] = -(1.0/12)*(z-1)*(1+z)*(1+y)*(-1+x)*(3*x-3*y+7*z2+3);
  phi[40] = (7.0/96)*(-1+x)*(1+y)*(1+z)*(z-1)*(7*z+sqrt21)*z;
  phi[41] = (7.0/96)*(1+x)*(1+y)*(1+z)*(z-1)*(-7*z+sqrt21)*z;
  phi[42] = -(1.0/12)*(z-1)*(1+z)*(1+y)*(1+x)*(3*x+3*y-3-7*z2);
  phi[43] = -(7.0/96)*(z-1)*(1+z)*(1+y)*(1+x)*z*(7*z+sqrt21);

  // face basis functions
 
  phi[44] = -(1.0/2)*(y-1)*(1+y)*(-1+x)*(1+x)*(z-1);
  phi[45] = (1.0/2)*(y-1)*(1+y)*(-1+x)*(1+x)*(1+z);
  phi[46] = -(1.0/2)*(z-1)*(1+z)*(-1+x)*(1+x)*(y-1);
  phi[47] = (1.0/2)*(z-1)*(1+z)*(-1+x)*(1+x)*(1+y);
  phi[48] = -(1.0/2)*(z-1)*(1+z)*(y-1)*(1+y)*(-1+x);
  phi[49] = (1.0/2)*(z-1)*(1+z)*(y-1)*(1+y)*(1+x);

  if(dphi == NULL) return;

  dphi[0] = svector<3>((1.0/32)*(z-1)*(y-1)*(7*y3-4*y*z-8*x*y-7*y+7*z3-8*x*z-7*z+28*x3+7-21*x2-22*x),
		       (1.0/32)*(z-1)*(-1+x)*(7*x3-7*x-4*x*z-8*x*y-8*y*z+7*z3+7-22*y-7*z+28*y3-21*y2),
		       (1.0/32)*(y-1)*(-1+x)*(7*x3-4*x*y-7*x-8*x*z+28*z3-8*y*z+7*y3-7*y-22*z-21*z2+7));
  dphi[1] = svector<3>((1.0/32)*(z-1)*(y-1)*(-7*y3+4*y*z-8*x*y+7*y-7*z3-8*x*z+7*z+28*x3-7+21*x2-22*x),
		       (1.0/32)*(z-1)*(1+x)*(7*x3-7*x-4*x*z-8*x*y-7*z3-7+8*y*z+22*y+7*z-28*y3+21*y2),
		       (1.0/32)*(y-1)*(1+x)*(7*x3-4*x*y-7*x-8*x*z-7*y3+7*y+21*z2+8*y*z-7-28*z3+22*z));
  dphi[2] = svector<3>(-(1.0/32)*(z-1)*(1+y)*(-7*y3+4*y*z+8*x*y+7*y+7*z3-8*x*z-7*z+7+28*x3-22*x-21*x2),
		       -(1.0/32)*(z-1)*(-1+x)*(7*x3-7*x-4*x*z+8*x*y+8*y*z+7*z3+7+22*y-7*z-28*y3-21*y2),
		       -(1.0/32)*(1+y)*(-1+x)*(7*x3-7*x+4*x*y-8*x*z-7*y3-22*z+28*z3+7*y-21*z2+7+8*y*z));
  dphi[3] = svector<3>(-(1.0/32)*(z-1)*(1+y)*(7*y3-4*y*z+8*x*y-7*y-7*z3-8*x*z+7*z+28*x3-7+21*x2-22*x),
		       -(1.0/32)*(z-1)*(1+x)*(7*x3-7*x-4*x*z+8*x*y-7*z3-7-8*y*z-22*y+7*z+28*y3+21*y2),
		       -(1.0/32)*(1+y)*(1+x)*(7*x3-7*x+4*x*y-8*x*z+7*y3-7*y+22*z-8*y*z-28*z3-7+21*z2));
  dphi[4] = svector<3>(-(1.0/32)*(1+z)*(y-1)*(7*y3+4*y*z-8*x*y-7*y-7*z3+8*x*z+7*z+7+28*x3-22*x-21*x2),
		       -(1.0/32)*(1+z)*(-1+x)*(7*x3+4*x*z-8*x*y-7*x+7*z+8*y*z+7-21*y2+28*y3-22*y-7*z3),
		       -(1.0/32)*(y-1)*(-1+x)*(7*x3-4*x*y-7*x+8*x*z-28*z3+8*y*z+7*y3-7*y+22*z-21*z2+7));
  dphi[5] = svector<3>(-(1.0/32)*(1+z)*(y-1)*(-7*y3-4*y*z-8*x*y+7*y+7*z3+8*x*z-7*z+28*x3-22*x+21*x2-7),
		       -(1.0/32)*(1+z)*(1+x)*(7*x3+4*x*z-8*x*y-7*x+7*z3+22*y-8*y*z+21*y2-7*z-7-28*y3),
		       -(1.0/32)*(y-1)*(1+x)*(7*x3-4*x*y-7*x+8*x*z-7*y3+7*y+21*z2-8*y*z-7+28*z3-22*z));
  dphi[6] = svector<3>((1.0/32)*(1+z)*(1+y)*(-7*y3-4*y*z+8*x*y+7*y-7*z3+8*x*z+7*z+7+28*x3-22*x-21*x2),
		       (1.0/32)*(1+z)*(-1+x)*(7*x3+4*x*z+8*x*y-7*x+7*z-8*y*z+7-21*y2-28*y3+22*y-7*z3),
		       (1.0/32)*(1+y)*(-1+x)*(7*x3-7*x+4*x*y+8*x*z-7*y3+22*z-28*z3+7*y-21*z2+7-8*y*z));
  dphi[7] = svector<3>((1.0/32)*(1+z)*(1+y)*(7*y3+4*y*z+8*x*y-7*y+7*z3+8*x*z-7*z+28*x3-22*x+21*x2-7),
		       (1.0/32)*(1+z)*(1+x)*(7*x3+4*x*z+8*x*y-7*x+7*z3-22*y+8*y*z+21*y2-7*z-7+28*y3),
		       (1.0/32)*(1+y)*(1+x)*(7*x3-7*x+4*x*y+8*x*z+7*y3-7*y-22*z+8*y*z+28*z3-7+21*z2));
  dphi[8] = svector<3>((7.0/96)*(-4*x2+sqrt21*x-1)*(7*x+sqrt21)*(y-1)*(z-1),
		       (7.0/96)*(1+x)*x*(-7*x+sqrt21)*(-1+x)*(z-1),
		       (7.0/96)*(-7*x+sqrt21)*(1+x)*x*(-1+x)*(y-1));
  dphi[9] = svector<3>((1.0/6)*x*(z-1)*(y-1)*(3*y+3*z-4+14*x2),
		       (1.0/12)*(-1+x)*(1+x)*(z-1)*(7*x2+3*z+6*y),
		       (1.0/12)*(-1+x)*(1+x)*(y-1)*(7*x2+3*y+6*z));
  dphi[10] = svector<3>((7.0/96)*(4*x2+sqrt21*x+1)*(-7*x+sqrt21)*(y-1)*(z-1),
			-(7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(z-1),
			-(7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(y-1));
  dphi[11] = svector<3>(-(7.0/96)*(4*x2+sqrt21*x+1)*(-7*x+sqrt21)*(1+y)*(z-1),
			(7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(z-1),
			(7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(1+y));
  dphi[12] = svector<3>(-(1.0/6)*x*(z-1)*(1+y)*(-3*y+3*z-4+14*x2),
			-(1.0/12)*(-1+x)*(1+x)*(z-1)*(7*x2+3*z-6*y),
			-(1.0/12)*(-1+x)*(1+x)*(1+y)*(7*x2-3*y+6*z));
  dphi[13] = svector<3>(-(7.0/96)*(-4*x2+sqrt21*x-1)*(7*x+sqrt21)*(1+y)*(z-1),
			-(7.0/96)*(1+x)*x*(-7*x+sqrt21)*(-1+x)*(z-1),
			-(7.0/96)*(1+x)*x*(-7*x+sqrt21)*(-1+x)*(1+y));
  dphi[14] = svector<3>(-(7.0/96)*(4*x2+sqrt21*x+1)*(-7*x+sqrt21)*(y-1)*(1+z),
			(7.0/96)*(1+x)*x*(-1+x)*(7*x+sqrt21)*(1+z),
			(7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(y-1));
  dphi[15] = svector<3>(-(1.0/6)*x*(1+z)*(y-1)*(3*y-3*z-4+14*x2),
			-(1.0/12)*(-1+x)*(1+x)*(1+z)*(7*x2-3*z+6*y),
			-(1.0/12)*(-1+x)*(1+x)*(y-1)*(7*x2+3*y-6*z));
  dphi[16] = svector<3>(-(7.0/96)*(-4*x2+sqrt21*x-1)*(7*x+sqrt21)*(y-1)*(1+z),
			-(7.0/96)*(-7*x+sqrt21)*(1+x)*x*(-1+x)*(1+z),
			-(7.0/96)*(-7*x+sqrt21)*(1+x)*x*(-1+x)*(y-1));
  dphi[17] = svector<3>((7.0/96)*(4*x2+sqrt21*x+1)*(-7*x+sqrt21)*(1+y)*(1+z),
			-(7.0/96)*(1+x)*x*(-1+x)*(7*x+sqrt21)*(1+z),
			-(7.0/96)*(7*x+sqrt21)*(1+x)*x*(-1+x)*(1+y));
  dphi[18] = svector<3>((1.0/6)*x*(1+z)*(1+y)*(-3*y-3*z-4+14*x2),
			(1.0/12)*(-1+x)*(1+x)*(1+z)*(7*x2-3*z-6*y),
			(1.0/12)*(-1+x)*(1+x)*(1+y)*(7*x2-3*y-6*z));
  dphi[19] = svector<3>((7.0/96)*(-4*x2+sqrt21*x-1)*(7*x+sqrt21)*(1+y)*(1+z),
			(7.0/96)*(-7*x+sqrt21)*(1+x)*x*(-1+x)*(1+z),
			(7.0/96)*(1+x)*x*(-7*x+sqrt21)*(-1+x)*(1+y));
  dphi[20] = svector<3>(-(7.0/96)*y*(y-1)*(7*y+sqrt21)*(1+y)*(z-1),
			(7.0/96)*(-1+x)*(4*y2+sqrt21*y+1)*(-7*y+sqrt21)*(z-1),
			-(7.0/96)*(-1+x)*y*(y-1)*(7*y+sqrt21)*(1+y));
  dphi[21] = svector<3>((1.0/12)*(y-1)*(1+y)*(z-1)*(7*y2+3*z+6*x),
			(1.0/6)*y*(z-1)*(-1+x)*(3*x-4+14*y2+3*z),
			(1.0/12)*(y-1)*(1+y)*(-1+x)*(3*x+7*y2+6*z));
  dphi[22] = svector<3>((7.0/96)*y*(-7*y+sqrt21)*(y-1)*(1+y)*(z-1),
			(7.0/96)*(-1+x)*(-4*y2+sqrt21*y-1)*(7*y+sqrt21)*(z-1),
			(7.0/96)*(-1+x)*y*(-7*y+sqrt21)*(y-1)*(1+y));
  dphi[23] = svector<3>(-(7.0/96)*y*(-7*y+sqrt21)*(y-1)*(1+y)*(z-1),
			-(7.0/96)*(1+x)*(-4*y2+sqrt21*y-1)*(7*y+sqrt21)*(z-1),
			-(7.0/96)*(1+x)*y*(y-1)*(-7*y+sqrt21)*(1+y));
  dphi[24] = svector<3>((1.0/12)*(y-1)*(1+y)*(z-1)*(-7*y2-3*z+6*x),
			(1.0/6)*y*(z-1)*(1+x)*(3*x-3*z+4-14*y2),
			(1.0/12)*(y-1)*(1+y)*(1+x)*(3*x-7*y2-6*z));
  dphi[25] = svector<3>((7.0/96)*y*(y-1)*(7*y+sqrt21)*(1+y)*(z-1),
			-(7.0/96)*(1+x)*(4*y2+sqrt21*y+1)*(-7*y+sqrt21)*(z-1),
			(7.0/96)*(1+x)*y*(7*y+sqrt21)*(y-1)*(1+y));
  dphi[26] = svector<3>((7.0/96)*y*(7*y+sqrt21)*(y-1)*(1+y)*(1+z),
			-(7.0/96)*(-1+x)*(4*y2+sqrt21*y+1)*(-7*y+sqrt21)*(1+z),
			(7.0/96)*(-1+x)*y*(y-1)*(7*y+sqrt21)*(1+y));
  dphi[27] = svector<3>(-(1.0/12)*(y-1)*(1+y)*(1+z)*(7*y2-3*z+6*x),
			-(1.0/6)*y*(1+z)*(-1+x)*(3*x-4-3*z+14*y2),
			-(1.0/12)*(y-1)*(1+y)*(-1+x)*(3*x+7*y2-6*z));
  dphi[28] = svector<3>(-(7.0/96)*y*(-7*y+sqrt21)*(y-1)*(1+y)*(1+z),
			-(7.0/96)*(-1+x)*(-4*y2+sqrt21*y-1)*(7*y+sqrt21)*(1+z),
			-(7.0/96)*(-1+x)*y*(-7*y+sqrt21)*(y-1)*(1+y));
  dphi[29] = svector<3>((7.0/96)*y*(-7*y+sqrt21)*(y-1)*(1+y)*(1+z),
			(7.0/96)*(1+x)*(-4*y2+sqrt21*y-1)*(7*y+sqrt21)*(1+z),
			(7.0/96)*(1+x)*y*(y-1)*(-7*y+sqrt21)*(1+y));
  dphi[30] = svector<3>(-(1.0/12)*(y-1)*(1+y)*(1+z)*(-7*y2+3*z+6*x),
			-(1.0/6)*y*(1+z)*(1+x)*(3*x+3*z-14*y2+4),
			-(1.0/12)*(y-1)*(1+y)*(1+x)*(3*x-7*y2+6*z));
  dphi[31] = svector<3>(-(7.0/96)*y*(7*y+sqrt21)*(y-1)*(1+y)*(1+z),
			(7.0/96)*(1+x)*(4*y2+sqrt21*y+1)*(-7*y+sqrt21)*(1+z),
			-(7.0/96)*(1+x)*y*(7*y+sqrt21)*(y-1)*(1+y));
  dphi[32] = svector<3>((7.0/96)*(y-1)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			(7.0/96)*(-1+x)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			(7.0/96)*(-1+x)*(y-1)*(-4*z2+sqrt21*z-1)*(7*z+sqrt21));
  dphi[33] = svector<3>((1.0/12)*(1+z)*(y-1)*(z-1)*(3*y+7*z2+6*x),
			(1.0/12)*(-1+x)*(1+z)*(z-1)*(3*x+7*z2+6*y),
			(1.0/6)*z*(y-1)*(-1+x)*(3*x+3*y+14*z2-4));
  dphi[34] = svector<3>(-(7.0/96)*(y-1)*(7*z+sqrt21)*(1+z)*(z-1)*z,
			-(7.0/96)*(-1+x)*(7*z+sqrt21)*(1+z)*(z-1)*z,
			(7.0/96)*(-1+x)*(y-1)*(4*z2+sqrt21*z+1)*(-7*z+sqrt21));
  dphi[35] = svector<3>(-(7.0/96)*(y-1)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			-(7.0/96)*(1+x)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			-(7.0/96)*(1+x)*(y-1)*(-4*z2+sqrt21*z-1)*(7*z+sqrt21));
  dphi[36] = svector<3>((1.0/12)*(1+z)*(y-1)*(z-1)*(-3*y-7*z2+6*x),
			(1.0/12)*(1+x)*(1+z)*(z-1)*(3*x-7*z2-6*y),
			(1.0/6)*z*(y-1)*(1+x)*(3*x-3*y-14*z2+4));
  dphi[37] = svector<3>((7.0/96)*(y-1)*(7*z+sqrt21)*(1+z)*(z-1)*z,
			(7.0/96)*(1+x)*(1+z)*(z-1)*(7*z+sqrt21)*z,
			-(7.0/96)*(1+x)*(y-1)*(4*z2+sqrt21*z+1)*(-7*z+sqrt21));
  dphi[38] = svector<3>(-(7.0/96)*(1+y)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			-(7.0/96)*(-1+x)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			-(7.0/96)*(-1+x)*(1+y)*(-4*z2+sqrt21*z-1)*(7*z+sqrt21));
  dphi[39] = svector<3>(-(1.0/12)*(1+z)*(1+y)*(z-1)*(-3*y+7*z2+6*x),
			-(1.0/12)*(-1+x)*(1+z)*(z-1)*(3*x+7*z2-6*y),
			-(1.0/6)*z*(1+y)*(-1+x)*(3*x+14*z2-4-3*y));
  dphi[40] = svector<3>((7.0/96)*(1+y)*(1+z)*(z-1)*(7*z+sqrt21)*z,
			(7.0/96)*(-1+x)*(7*z+sqrt21)*(1+z)*(z-1)*z,
			-(7.0/96)*(-1+x)*(1+y)*(4*z2+sqrt21*z+1)*(-7*z+sqrt21));
  dphi[41] = svector<3>((7.0/96)*(1+y)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			(7.0/96)*(1+x)*(1+z)*(z-1)*(-7*z+sqrt21)*z,
			(7.0/96)*(1+x)*(1+y)*(-4*z2+sqrt21*z-1)*(7*z+sqrt21));
  dphi[42] = svector<3>(-(1.0/12)*(1+z)*(1+y)*(z-1)*(3*y-7*z2+6*x),
			-(1.0/12)*(1+x)*(1+z)*(z-1)*(3*x-7*z2+6*y),
			-(1.0/6)*z*(1+y)*(1+x)*(3*x+3*y-14*z2+4));
  dphi[43] = svector<3>(-(7.0/96)*(1+y)*(1+z)*(z-1)*(7*z+sqrt21)*z,
			-(7.0/96)*(1+x)*(1+z)*(z-1)*(7*z+sqrt21)*z,
			(7.0/96)*(1+x)*(1+y)*(4*z2+sqrt21*z+1)*(-7*z+sqrt21));
  dphi[44] = svector<3>(-x*(y-1)*(1+y)*(z-1),
			-y*(-1+x)*(1+x)*(z-1),
			-(1.0/2)*(y-1)*(1+y)*(-1+x)*(1+x));
  dphi[45] = svector<3>(x*(y-1)*(1+y)*(1+z),
			y*(-1+x)*(1+x)*(1+z),
			(1.0/2)*(y-1)*(1+y)*(-1+x)*(1+x));
  dphi[46] = svector<3>(-x*(z-1)*(1+z)*(y-1),
			-(1.0/2)*(z-1)*(1+z)*(-1+x)*(1+x),
			-z*(-1+x)*(1+x)*(y-1));
  dphi[47] = svector<3>(x*(z-1)*(1+z)*(1+y),
			(1.0/2)*(z-1)*(1+z)*(-1+x)*(1+x),
			z*(-1+x)*(1+x)*(1+y));
  dphi[48] = svector<3>(-(1.0/2)*(z-1)*(1+z)*(y-1)*(1+y),
			-y*(z-1)*(1+z)*(-1+x),
			-z*(y-1)*(1+y)*(-1+x));
  dphi[49] = svector<3>((1.0/2)*(z-1)*(1+z)*(y-1)*(1+y),
			y*(z-1)*(1+z)*(1+x),
			z*(y-1)*(1+y)*(1+x));     

  return;
}

void cube50Permutation(const int *nodes, int *dof)
{
  int temp;

  int ee[12][2] = {{0,1}, {3,2}, {5,4}, {7,6}, {2,0}, {1,3}, 
		   {6,4}, {5,7}, {0,4}, {1,5}, {2,6}, {3,7}};

  // permute the edge dof

  for(int i = 0; i < 12; i++)
    {
      if(nodes[ee[i][0]] > nodes[ee[i][1]])  // swap the edge dof
	{
	  temp = dof[8 + 3*i]; 
	  dof[8  + 3*i] = dof[10 + 3*i];
	  dof[10 + 3*i] = temp;
	}
    }

  return;
}
