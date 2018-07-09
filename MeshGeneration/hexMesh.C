// Set up a mesh with each hexahedra affine equivalent to the unit cube

#include<string>
#include<sstream>
#include<fstream>

#include"../LinearAlgebra/svector.C"
#include"gelement.C"
#include"../FiniteElement/basisFns.C"

#define HEXMESH_NCHI  8
#define HEXMESH_NCHIB 4

typedef GeometricElement<3, 8> gelement;
typedef BoundaryGelement<3, 4> bgelement;

class hexMesh
{
public:
  hexMesh(int);
  hexMesh(std::string);
  hexMesh(int, int*);

  ~hexMesh() { return; }
  
  int nNode()  const {return nodes;}
  int nElem()  const {return nelem;}
  int nBelem() const {return nbelem;}
  
  int locate(const svector<3>&, gelement&, svector<3>&) const;

  void info();
  void print();

  void print(std::ostream& dump);
  
  class geiterator {
  public:
    geiterator();
    geiterator(int l, hexMesh * g);
    void operator++();
    void operator++(int);
    geiterator& operator=(const geiterator & it);
    bool operator==(const geiterator & it);
    bool operator!=(const geiterator & it);
    gelement operator*();
    
  private:
    int loc;
    hexMesh * gm;
  };
  
  geiterator gebegin();
  geiterator geend();
  
  class bgeiterator {
    public:
    bgeiterator();
    bgeiterator(int l, hexMesh * g);
    void operator++();
    void operator++(int);
    bgeiterator& operator=(const bgeiterator & it);
    bool operator==(const bgeiterator & it);
    bool operator!=(const bgeiterator & it);
    bgelement operator*();

  private:
    int loc;
    hexMesh * gm;
  };
  
  bgeiterator bgebegin();
  bgeiterator bgeend();

  LagrangeBasisFunction<3> *gBasis;

  gelement  intToGe(int) const;
  bgelement intToBge(int) const;
  
private:
  int nelem;
  int nodes;
  int nbelem;
  
  int nex, ney, nez, bcType[6];

  double xmin, xmax, ymin, ymax, zmin, zmax;
};

hexMesh::hexMesh(std::string ss) : gBasis(Cube8)
{
  int nx;

  if(!ss.empty())
    {
      std::istringstream ist(ss);
     
     if( !(ist >> nx) )
       {
	 std::cerr << "hexMesh(string): Invalid integer: " << ss
		   << ", terminating" << std::endl;
	 throw("Invalid input");
       }
     
     if(nx <= 0) 
       {
	 std::cerr << "hexMesh(string): Invalid integer: " 
		   << nx << ", terminating" << std::endl;
	 throw("Invalid input argument");
       }
    }
 else
   {
     nx = 30;
   }

  nelem = nx*nx*nx; 
  nodes = ((nx+1)*(nx+1)*(nx+1));
  nbelem= (6*nx*nx);

  nex = (nx);
  ney = (nx); 
  nez = (nx); 

  xmin = (-1); xmax = (1); 
  ymin = (-1); ymax = (1);
  zmin = (-1); zmax = (1);

  for(int i = 0; i < 6; i++) bcType[i] = 0;

  bcType[3] = 1;   // Test Neumann BC

  info();
}

hexMesh::hexMesh(int nx) 
  : gBasis(Cube8),
    nelem(nx*nx*nx), nodes((nx+1)*(nx+1)*(nx+1)),
    nbelem(6*nx*nx), 
    nex(nx), ney(nx), nez(nx),
    xmin(-1), xmax(1), ymin(-1), ymax(1), zmin(-1), zmax(1)
{
  for(int i = 0; i < 6; i++) bcType[i] = 0;

  info();
}

hexMesh::hexMesh(int nx, int bct[6]) 
  : gBasis(Cube8),
    nelem(nx*nx*nx), nodes((nx+1)*(nx+1)*(nx+1)),
    nbelem(6*nx*nx), 
    nex(nx), ney(nx), nez(nx),
    xmin(-1), xmax(1), ymin(-1), ymax(1), zmin(-1), zmax(1)
{
  for(int i = 0; i < 6; i++) bcType[i] = bct[i];

  info();
}

void hexMesh::info()
{
  std::cout << "Initializing hexMesh(): hexagonal mesh" << std::endl;
  std::cout << "Boundary codes bcType[zmin,zmax,ymin,ymax,xmin,xmax] = ["
	    << bcType[0] << "," << bcType[1] << ","
	    << bcType[2] << "," << bcType[3] << "," 
	    << bcType[4] << "," << bcType[5] << "]" 
	    << std::endl;
  std::cout << "Number of elements in each x-direction = ......... " << nex 
	    << std::endl;
  std::cout << "Number of elements in each y-direction = ......... " << ney 
	    << std::endl;
  std::cout << "Number of elements in each z-direction = ......... " << nez
	    << std::endl;
  std::cout << "[xmin,xmax] = ............................... [" 
	    << xmin << "," << xmax << "]" << std::endl;
  std::cout << "[ymin,ymax] = ............................... [" 
	    << ymin << "," << ymax << "]" << std::endl;
  std::cout << "[zmin,zmax] = ............................... [" 
	    << zmin << "," << zmax << "]" << std::endl;

  return;
}

int hexMesh::locate(const svector<3>& xx, gelement& ee, svector<3>& xi) const
{
  int ie, je, ke;
  
  if((xx[0] < xmin) || (xx[0] > xmax) || 
     (xx[1] < ymin) || (xx[1] > ymax) ||
     (xx[2] < zmin) || (xx[2] > zmax))
    {
      //  std::cerr << "Point out of range: (" << xx[0] << "," << xx[1]
      //   << ") not in [xmin,ymin] x [xmax,ymax]" << std::endl;
      
      throw("Point out of range");
    }
     
  double dx = (xx[0]-xmin) / (xmax-xmin);
  double dy = (xx[1]-ymin) / (ymax-ymin);
  double dz = (xx[2]-zmin) / (zmax-zmin);

  ie = (int) (dx * (nex - 1.0e-14));  // if xx[0] = xmax evaluates to nex-1
  je = (int) (dy * (ney - 1.0e-14));  // if xx[1] = ymax evaluates to ney-1
  ke = (int) (dz * (nez - 1.0e-14));  // if xx[2] = zmax evaluates to nez-1

  xi = svector<3>(2.0*(dx*nex-ie) - 1.0, 
	       2.0*(dy*ney-je) - 1.0,
	       2.0*(dz*nez-ke) - 1.0);

  ee = intToGe(ke*nex*ney + je*nex + ie);

  return(ke*nex*ney + je*nex + ie);
}

hexMesh::geiterator::geiterator() {
    loc = 0;
    gm = NULL;
}

hexMesh::geiterator::geiterator(int l, hexMesh * g) {
    loc = l;
    gm = g;
}

void hexMesh::geiterator::operator++() {
    (*this)++;
}

void hexMesh::geiterator::operator++(int) {
    loc++;
}

hexMesh::geiterator& hexMesh::geiterator::operator=(const geiterator & it) {
    loc = it.loc;
    gm = it.gm;
    return *this;
}

bool hexMesh::geiterator::operator==(const geiterator & it) {
    if (gm != it.gm || loc != it.loc) {
	return false;
    }
    return true;
}

bool hexMesh::geiterator::operator!=(const geiterator & it) {
    return !(*this == it);
}

gelement hexMesh::geiterator::operator*() {

  gelement ee = gm->intToGe(loc);
    
  return ee;
}


hexMesh::geiterator hexMesh::gebegin() {
    return geiterator(0, this);
}

hexMesh::geiterator hexMesh::geend() {
    return geiterator(nelem, this);
}

hexMesh::bgeiterator::bgeiterator() {
    loc = 0;
    gm = NULL;
}

hexMesh::bgeiterator::bgeiterator(int l, hexMesh * g) {
    loc = l;
    gm = g;
}

void hexMesh::bgeiterator::operator++() {
    (*this)++;
}

void hexMesh::bgeiterator::operator++(int) {
    loc++;
}

hexMesh::bgeiterator& hexMesh::bgeiterator::operator=(const bgeiterator & it) {
    loc = it.loc;
    gm = it.gm;
    return *this;
}

bool hexMesh::bgeiterator::operator==(const bgeiterator & it) {
    if (loc != it.loc || gm != it.gm) {
	return false;
    }
    return true;
}

bool hexMesh::bgeiterator::operator!=(const bgeiterator & it) {
    return !(*this == it);
}

bgelement hexMesh::bgeiterator::operator*() {

  bgelement be = gm->intToBge(loc);

  return be;

}

hexMesh::bgeiterator hexMesh::bgebegin() {
    return bgeiterator(0,this);
}

hexMesh::bgeiterator hexMesh::bgeend() {
    return bgeiterator(nbelem,this);
}

gelement hexMesh::intToGe(int n) const
{
  int ie, je, ke, nnx, nny, nnz;
  double hx,hy,hz, xemin[3],xemax[3];
  gelement ee;

  // *********** Set up for Hex 8 Elements **************
  // *********** Set up for Hex 8 Elements **************
  // *********** Set up for Hex 8 Elements **************

  if((n < 0) || (n >= nelem))
    {
      std::cerr << "hexMesh::element() Element number " << n 
		<< " out of range, nelem = " << nelem << std::endl;
      
      throw("Element number out of range");
    }

  ke =  n / (nex * ney);
  je = (n - ke*nex*ney) / nex;
  ie =  n - ke*nex*ney - je * nex;

  nnx = nex + 1;
  nny = ney + 1;
  nnz = nez + 1;

  ee.nodes[0] = ( ke   *nny +  je   )*nnx +  ie;
  ee.nodes[1] = ( ke   *nny +  je   )*nnx + (ie+1);
  ee.nodes[2] = ( ke   *nny + (je+1))*nnx +  ie;
  ee.nodes[3] = ( ke   *nny + (je+1))*nnx + (ie+1);
  ee.nodes[4] = ((ke+1)*nny +  je   )*nnx +  ie;
  ee.nodes[5] = ((ke+1)*nny +  je   )*nnx + (ie+1);
  ee.nodes[6] = ((ke+1)*nny + (je+1))*nnx +  ie;
  ee.nodes[7] = ((ke+1)*nny + (je+1))*nnx + (ie+1);

 // coordinates

  hx = (xmax-xmin) / nex;
  hy = (ymax-ymin) / ney;
  hz = (zmax-zmin) / nez;
  
  xemin[0] = xmin + ie*hx; xemax[0] = xemin[0] + hx;
  xemin[1] = ymin + je*hy; xemax[1] = xemin[1] + hy;
  xemin[2] = zmin + ke*hz; xemax[2] = xemin[2] + hz;

  ee.xx[0][0] = xemin[0]; ee.xx[0][1] = xemin[1]; ee.xx[0][2] = xemin[2];
  ee.xx[1][0] = xemax[0]; ee.xx[1][1] = xemin[1]; ee.xx[1][2] = xemin[2];
  ee.xx[2][0] = xemin[0]; ee.xx[2][1] = xemax[1]; ee.xx[2][2] = xemin[2];
  ee.xx[3][0] = xemax[0]; ee.xx[3][1] = xemax[1]; ee.xx[3][2] = xemin[2];
  ee.xx[4][0] = xemin[0]; ee.xx[4][1] = xemin[1]; ee.xx[4][2] = xemax[2];
  ee.xx[5][0] = xemax[0]; ee.xx[5][1] = xemin[1]; ee.xx[5][2] = xemax[2];
  ee.xx[6][0] = xemin[0]; ee.xx[6][1] = xemax[1]; ee.xx[6][2] = xemax[2];
  ee.xx[7][0] = xemax[0]; ee.xx[7][1] = xemax[1]; ee.xx[7][2] = xemax[2];

  // set the transformation used for the geometry

  ee.master = gBasis->master;
  ee.basis  = gBasis->basis;

  return(ee);
}

bgelement hexMesh::intToBge(int n) const
{
  int fid;
  gelement ee;
  bgelement be;

  // *********** Set up for Hex 8 Elements **************
  // *********** Set up for Hex 8 Elements **************
  // *********** Set up for Hex 8 Elements **************

 if((n < 0) || (n >= nbelem))
   {
     std::cerr << "Boundary element number " << n << " out of range, nelem = " 
	       << nelem << std::endl;

     throw("Boundary element number out of range");
   }

// bcType[6] is an int array.
//
// Bottom: type = bcType[0]    z=0
// Top:    type = bcType[1]    z=1
// Front:  type = bcType[2]    y=0
// Back:   type = bcType[3]    y=1
// Letf:   type = bcType[4]    x=0
// Right:  type = bcType[5]    x=1
//

 if(n < nex*ney)  // bottom  
  {
    fid = 0;

    ee = intToGe(n);
  }
 else if(n < 2*(nex*ney)) // top
   {
    fid = 1;

    ee = intToGe( (n - nex*ney) + nex*ney*(nez-1) );
  }
 else if(n < 2*nex*ney + nex*nez )  // front
  {
    fid = 2;

    int nxz = n - 2*(nex*ney);

    int ke = nxz / nex;
    int ie = nxz - ke * nex;

    ee = intToGe( ke*(nex*ney) + ie );
  }
 else if(n < 2*(nex*ney + nex*nez) ) // back
   {
    fid = 3;

    int nxz = n - 2*nex*ney - nex*nez;

    int ke = nxz / nex;
    int ie = nxz - ke * nex;

    ee = intToGe( ke*(nex*ney) + nex*(ney-1) + ie );
   }
 else if(n < 2*(nex*ney + nex*nez) + ney*nez)  // lhs
  {
    fid = 4;

    int nyz = n - 2*(nex*ney + nex*nez);

    int ke = nyz / ney;
    int je = nyz - ke * ney;

    ee = intToGe( ke*(nex*ney) + je*nex );
  }
 else if(n < 2*(nex*ney + nex*nez + ney*nez)) // rhs
   {
    fid = 5;

    int nyz = n - 2*(nex*ney + nex*nez) - ney*nez;

    int ke = nyz / ney;
    int je = nyz - ke * ney;

    ee = intToGe( ke*(nex*ney) + je*nex + (nex-1) );
   }
 else
   {
     std::cerr << "hexMesh::intToBge() Boundary element number " << n 
	       << " out of range, nbelem = " << nbelem << std::endl;

     throw("Boundary element number out of range");
   }

 // extract the face data from the element

 be.type = bcType[fid];

 for(int i = 0; i < 4; i++)
   {
     be.xx[i] = ee.xx[ (ee.master)->subs(2,fid,i) ];

     be.nodes[i] = ee.nodes[ (ee.master)->subs(2,fid,i) ];
   }

 // for pure Dirichlet boundary data adjust the
 // zeroth boundary element in case pressure needs specification

 if( (bcType[0] == 0) && (bcType[1] == 0) &&
     (bcType[2] == 0) && (bcType[3] == 0) &&
     (bcType[4] == 0) && (bcType[5] == 0) && (n==0) ) be.type += 10;

 be.master = gBasis->boundary->master;
 be.basis  = gBasis->boundary->basis;

 return(be);
}

void hexMesh::print()
{
  int i;

  gelement ge;
  bgelement be;

  std::ofstream dump("meshOut");

  dump.setf(std::ios::fixed, std::ios::floatfield);

  dump << "Number of elements          = " << nelem << std::endl;
  dump << "Number of boundary elements = " << nbelem << std::endl;
  dump << "Number of nodes             = " << nodes << std::endl;

  for(i = 0; i < nelem; i++)
    {
      dump << "Element " << i << std::endl;

      ge = intToGe(i);

      // ge.print(dump);
    }

  dump << std::endl;

  for(i = 0; i < nbelem; i++)
    {
      dump << "Boundary Element " << i << ",";

      be = intToBge(i);

      // be.print(dump);
    }

  dump << std::endl;

  dump.close();
}
