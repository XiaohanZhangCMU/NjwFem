// Set up a mesh with each quadralateral "affine" equivalent to the unit square

#include<string>
#include<sstream>
#include<fstream>

#include<vector>
#include"../LinearAlgebra/svector.C"
#include"gelement.C"
#include"../FiniteElement/basisFns.C"

#define QUADMESH_NCHI  4
#define QUADMESH_NCHIB 2

//#define MYQUAD

typedef GeometricElement<2, 4> gelement;
typedef BoundaryGelement<2, 2> bgelement;

class quadMesh
{
public:
  quadMesh(int);
  quadMesh(int, int);
  quadMesh(int, int*);
  quadMesh(int, int, int*);

  quadMesh(std::string ss);

  int nNode()  const {return nodes;}
  int nElem()  const {return nelem;}
  int nBelem() const {return nbelem;}
  int getNx() const {return nex;}
  int getNumElmsInLayer() {return this->getNx();} //return number of 1 layer of elements
  int getNy() const {return ney;}
  int getLayer() const {return layer;}
  double getXmax() const{return xmax;}
  double getXmin() const{return xmin;}

  double getYmax() const{return ymax;}
  double getYmin() const{return ymin;}

  int getBndyElm() const{return nbelem;}
  void locate(const svector<2>&, gelement&, svector<2>&) const;
  void setBoundaryElems(int); 
  int getNodeByElm(int);
  double setLayerElems(int, double, double, vector<double>& );
  vector<double> params(); 
  vector<double> center_elem(int);
  
  void info();
  void print();
  void print(std::ostream&);

  class geiterator 
  {
  public:
    geiterator();
    geiterator(int l, quadMesh * g);
    void operator++();
    void operator++(int);
    geiterator& operator=(const geiterator & it);
    bool operator==(const geiterator & it);
    bool operator!=(const geiterator & it);
    gelement operator*();  
  private:
    int loc;
    quadMesh * gm;
  };
  
  geiterator gebegin();
  geiterator geend();
  
  class bgeiterator 
  {
  public:
    bgeiterator();
    bgeiterator(int l, quadMesh * g);
    void operator++();
    void operator++(int);
    bgeiterator& operator=(const bgeiterator & it);
    bool operator==(const bgeiterator & it);
    bool operator!=(const bgeiterator & it);
    bgelement operator*();
  private:
    int loc;
    quadMesh * gm;
  };
  
  bgeiterator bgebegin();
  bgeiterator bgeend();

  LagrangeBasisFunction<2> *gBasis;

  gelement  intToGe(int) const;
  bgelement intToBge(int) const;

private:
  int nelem;
  int nodes;
  int nbelem;
  
  int nex, ney;
  int bandth;
  int layer;
  double xmin, xmax, ymin, ymax, xcc, ycc;
  int bcType[4];
};

quadMesh::quadMesh(std::string ss) : gBasis(Square4)
{
  int nx;
  int ny;

  if(!ss.empty())
    {
      std::stringstream ist(ss);
     
     if( !(ist >> nx) )
       {
	 std::cerr << "quadMesh(string): Invalid integer: " << ss
		   << ", terminating" << std::endl;
	 throw("Invalid input");
       }
     
     if(nx <= 0) 
       {
	 std::cerr << "quadMesh(string): Invalid integer: " 
		   << nx << ", terminating" << std::endl;
	 throw("Invalid input argument");
       }
    }
 else
   {
     nx = 551; /////////////g//////////////////////////////8;
     ny = 131;
     layer = 1;
   }


  nelem = nx*ny; 
  nodes = ((nx+1)*(ny+1));
  nbelem= 2*nx+2*ny;

  nex = (nx);
  ney = (ny);

  //xmin = (-0.5); xmax = (0.5); 
  //ymin = (-0.5); ymax = (0.5);

  for(int i = 0; i < 4; i++) bcType[i] = 0;

  // bcType[2] = 1;   // Test Neumann BC

//  bcType[1] = 1;   // Test Neumann BC
//  std::cout << "quadMesh(): BC on RHS for wave equation test" << std::endl;

//   bcType[0] = bcType[2] = 1;   // Test Neumann BC
//   std::cout 
//	<< "quadMesh(): Neumann BC on Bottom and Top for wave equation test" 
//      	<< std::endl;

  // Mesh and BC for Examples in Chapters 8 and 9 of FEbook
  // using Example No. 1.

  //std::cout 
  // << "quadMesh(): Mesh and BC for Examples in Chapters 8 and 9 of FEbook" 
  // << std::endl;

  // xmin = (-0.5); xmax = (0.5); 
  // ymin = (-0.5); ymax = (0.5);

  // std::cout[botom,rhs,top,lhs] = ["
  //	    << bcType[0] << bcType[1] << bcType[2] << bcType[3]
  
  // //quasi-
    xmin = -100; xmax = 100; xcc = 0;
    ymin = -80.5; ymax = 80.5; ycc = 0;
 
  // //xiaohan
   bcType[0] = 1;
   bcType[1] = 1;
   bcType[2] = 1;
   bcType[3] = 0;
    
   
   // //dyn-
   // xmin = 0; xmax = 5;
   // ymin = 0; ymax = 1;
 
   // //xiaohan
   // bcType[0] = 0;
   // bcType[1] = 0;
   // bcType[2] = 0;
   // bcType[3] = 0;
   

  info();
}

vector<double> quadMesh::center_elem(int tag)
{
  svector<2> xx;
  gelement ee = intToGe(tag);
  
  for(int i = 0;i<Nchi;i++)
    {
      xx += ee.xx[i];
    }
  double bar_x = xx[0]/((double)Nchi);
  double bar_y = xx[1]/((double)Nchi);
  double vals[] = {bar_x, bar_y};
  vector<double> v(vals, vals + sizeof(vals) / sizeof(double) );
  return v;
}

vector<double> quadMesh::params()
{
  double myVals[] = {xmin, ymin, xmax, ymax, xcc, ycc};
  vector<double> v(myVals, myVals+sizeof(myVals)/sizeof(double));
  return v;
}

void quadMesh:: setBoundaryElems(int rb)
{
  return;
}

double quadMesh:: setLayerElems(int layer, double layer_top, double layer_bot, vector<double>& imagNodes )
{
  int n = nex;
  //double height_layer = ymax - ymin;
  double length_layer = xmax-xmin;
  //  double quadElmHeight = height_layer/((double)layer);
  double dx = length_layer/((double)n);
  double quadElmLength = length_layer/((double)n);
  
  for(int i = 0;i<n;i++)
    {
      if(i==0)
        imagNodes.push_back( xmin+quadElmLength/2.0);
      else
        imagNodes.push_back(imagNodes[i-1]+quadElmLength);
    }

  return dx;
}

int quadMesh::getNodeByElm(int tag)
{
  //return -1 if not contains a 1d node, +n if contains the index of 1d node of imageNodes[]                        

  int n = nex;                                                                                                   
  svector<2> xx;
  gelement ee = intToGe(tag);

  for(int i = 0;i<Nchi;i++)
    {
      xx += ee.xx[i];
    }
  double bar_x = xx[0]/((double)Nchi);
  double bar_y = xx[1]/((double)Nchi);

    // std::cout<<"layer = "<<layer<<std::endl;                                                                       
  double height_layer = ymax - ymin;
  double length_layer = xmax-xmin;
  double quadElmHeight = height_layer/((double)ney);
  //  double dx = length_layer/((double)n);
  double quadElmLength = length_layer/((double)n);

  //  std::cout<<"quadelmhei = " <<quadElmHeight <<" bar_x = "<<bar_x << " bar_y = "<<bar_y<<" Nchi = "<<Nchi<<" layer*quadElmHeight/2.0 = "<<layer*quadElmHeight/2.0<<std::endl;                                     

  if((bar_y)>(layer*quadElmHeight/2.0)+ycc || (bar_y)<(-1*layer*quadElmHeight/2.0)+ycc )
    return -1;
  else
    {
      double temp = (bar_x-xmin)/quadElmLength;
      int index = (int)(temp);
      return index;
    }
}


quadMesh::quadMesh(int nx)
  : gBasis(Square4),
    nelem(nx*nx), nodes((nx+1)*(nx+1)),
    nbelem(4*nx), nex(nx), ney(nx), 
    xmin(-1), xmax(1), ymin(-1), ymax(1)
{
  for(int i = 0; i < 4; i++) bcType[i] = 0;

  bcType[3] = 1;   // Test Neumann BC

  info();
}

quadMesh::quadMesh(int nx, int bct[4]) 
  : gBasis(Square4),
    nelem(nx*nx), nodes((nx+1)*(nx+1)),
    nbelem(4*nx), nex(nx), ney(nx), 
    xmin(-1), xmax(1), ymin(-1), ymax(1)
{
  for(int i = 0; i < 4; i++) bcType[i] = bct[i];

  info();
}

quadMesh::quadMesh(int nx, int ny) 
  : gBasis(Square4),
    nelem(nx*ny), nodes((nx+1)*(ny+1)),
    nbelem(2*(nx+ny)), nex(nx), ney(ny),
    xmin(-1), xmax(1), ymin(-1), ymax(1)
{
  for(int i = 0; i < 4; i++) bcType[i] = 0;

  bcType[3] = 1;   // Test Neumann BC

  info();
}

quadMesh::quadMesh(int nx, int ny, int bct[4]) 
  : gBasis(Square4),
    nelem(nx*ny), nodes((nx+1)*(ny+1)),
    nbelem(2*(nx+ny)), nex(nx), ney(ny),
    xmin(-1), xmax(1), ymin(-1), ymax(1)
{
  for(int i = 0; i < 4; i++) bcType[i] = bct[i];

  info();
}

void quadMesh::info()
{
  std::cout << "quadMesh::info(): quadrilateral mesh" << std::endl;
  std::cout << "Boundary codes bcType[botom,rhs,top,lhs] = ["
	    << bcType[0] << "," << bcType[1] << ","
	    << bcType[2] << "," << bcType[3] << "]" 
	    << std::endl;
  std::cout << "Number of elements in each x-direction = ......... " << nex 
	    << std::endl;
  std::cout << "Number of elements in each y-direction = ......... " << ney 
	    << std::endl;
  std::cout << "[xmin,xmax] = ............................... [" 
	    << xmin << "," << xmax << "]" << std::endl;
  std::cout << "[ymin,ymax] = ............................... [" 
	    << ymin << "," << ymax << "]" << std::endl;
  std::cout << "Boundary elements number = ["<<this->nbelem<<"]" << std::endl;
  return;
}

void quadMesh::locate(const svector<2>& xx, gelement& ee, svector<2>& xi) const
{
  int ie, je;

  if((xx[0] < xmin) || (xx[0] > xmax) || (xx[1] < ymin) || (xx[1] > ymax))
    {
      //  std::cerr << "Point out of range: (" << xx[0] << "," << xx[1]
      //   << ") not in [0,1]^2" << std::endl;

      throw("Point out of range");
    }

  double dx = (xx[0]-xmin) / (xmax-xmin);
  double dy = (xx[1]-ymin) / (ymax-ymin);

  ie = (int) (dx * (nex - 1.0e-14));  // if xx[0] = xmax evaluates to nex-1
  je = (int) (dy * (ney - 1.0e-14));  // if xx[1] = ymax evaluates to ney-1

  xi = svector<2>(2.0*(xx[0]*nex-ie) - 1.0, 2.0*(xx[1]*ney-je) - 1.0);

  ee = intToGe(je*nex+ie);

  return;
}

quadMesh::geiterator::geiterator() {
    loc = 0;
    gm = NULL;
}

quadMesh::geiterator::geiterator(int l, quadMesh * g) {
    loc = l;
    gm = g;
}

void quadMesh::geiterator::operator++() {
    (*this)++;
}

void quadMesh::geiterator::operator++(int) {
    loc++;
}

quadMesh::geiterator& quadMesh::geiterator::operator=(const geiterator & it) {
    loc = it.loc;
    gm = it.gm;
    return *this;
}

bool quadMesh::geiterator::operator==(const geiterator & it) {
    if (gm != it.gm || loc != it.loc) {
	return false;
    }
    return true;
}

bool quadMesh::geiterator::operator!=(const geiterator & it) {
    return !(*this == it);
}

gelement quadMesh::geiterator::operator*() {

  gelement ee = gm->intToGe(loc);
    
  return ee;
}


quadMesh::geiterator quadMesh::gebegin() {
    return geiterator(0, this);
}

quadMesh::geiterator quadMesh::geend() {
    return geiterator(nelem, this);
}

quadMesh::bgeiterator::bgeiterator() {
    loc = 0;
    gm = NULL;
}

quadMesh::bgeiterator::bgeiterator(int l, quadMesh * g) {
    loc = l;
    gm = g;
}

void quadMesh::bgeiterator::operator++() {
    (*this)++;
}

void quadMesh::bgeiterator::operator++(int) {
    loc++;
}

quadMesh::bgeiterator& quadMesh::bgeiterator::operator=(const bgeiterator & it) {
    loc = it.loc;
    gm = it.gm;
    return *this;
}

bool quadMesh::bgeiterator::operator==(const bgeiterator & it) {
    if (loc != it.loc || gm != it.gm) {
	return false;
    }
    return true;
}

bool quadMesh::bgeiterator::operator!=(const bgeiterator & it) {
    return !(*this == it);
}

bgelement quadMesh::bgeiterator::operator*() {

  bgelement be = gm->intToBge(loc);

  return be;

}

quadMesh::bgeiterator quadMesh::bgebegin() {
    return bgeiterator(0,this);
}

quadMesh::bgeiterator quadMesh::bgeend() {
    return bgeiterator(nbelem,this);
}

gelement quadMesh::intToGe(int n) const
{
int ie, je;
double hx,hy, x0, y0;
gelement ee;

//*********** Set up for Quad 4 Elements **************
//*********** Set up for Quad 4 Elements **************
//*********** Set up for Quad 4 Elements **************

 if((n < 0) || (n >= nelem))
   {
     std::cerr << "Element number " << n << " out of range, nelem = " 
	       << nelem << std::endl;

     throw("Element number out of range");
   }

 je = n / nex;
 ie = n - je * nex;

 ee.nodes[0] = (nex+1) * je + ie;
 ee.nodes[1] = ee.nodes[0] + 1;
 ee.nodes[3] = ee.nodes[0] + (nex+1);
 ee.nodes[2] = ee.nodes[3] + 1;
 
 hx = (xmax-xmin) / nex;
 hy = (ymax-ymin) / ney;

 x0 = xmin + ie * hx;
 y0 = ymin + je * hy;

 ee.xx[0] = svector<2>(x0,    y0);
 ee.xx[1] = svector<2>(x0+hx, y0);
 ee.xx[2] = svector<2>(x0+hx, y0+hy);
 ee.xx[3] = svector<2>(x0,    y0+hy);

 ee.tag = n;
 ee.master = gBasis->master;
 ee.basis  = gBasis->basis;
 
 return(ee);
}

/*
#ifdef MYQUAD
gelement quadMesh::intToGe(int n) const
{
  int ie, je;
  double x0, y0;
  gelement ee;


  if((n < 0) || (n >= nelem))
    {
      std::cerr << "Element number " << n << " out of range, nelem = " 
		<< nelem << std::endl;

      throw("Element number out of range");
    }

  if(bandth%2 != 1)
    std::cout<<"the bandth should be odd!"<<std::endl;

  double hx = (xmax - xmin)/nex;
  double hy_bandth = hx;
  double hy = (ymax-ymin-bandth*hy_bandth)/(ney-bandth);
  double ymin_bandth = -0.5 * hy_bandth - (bandth-1)/2 * hy_bandth;
  //double ymax_bandth = -1.0 * ymin_bandth;
  int nymin_bandth = (ney-bandth)/2;
  int nymax_bandth = nymin_bandth + (bandth-1);

  je = n / nex;
  ie = n - je * nex;

  ee.nodes[0] = (nex+1) * je + ie;
  ee.nodes[1] = ee.nodes[0] + 1;
  ee.nodes[3] = ee.nodes[0] + (nex+1);
  ee.nodes[2] = ee.nodes[3] + 1;
 
  x0 = xmin + ie * hx;



  if(je<=nymax_bandth && je>=nymin_bandth)
    {
      hy = hy_bandth;
      y0 = ymin_bandth+(je-nymin_bandth)*hy_bandth;
    }
  else if(je<nymin_bandth)
    {
      y0 = ymin + je * hy;
    }
  else
    {
      y0 = ymax - (ney-je) * hy;
    }

  ee.xx[0] = svector<2>(x0,    y0);
  ee.xx[1] = svector<2>(x0+hx, y0);
  ee.xx[2] = svector<2>(x0+hx, y0+hy);
  ee.xx[3] = svector<2>(x0,    y0+hy);

  ee.master = gBasis->master;
  ee.basis  = gBasis->basis;

  return(ee);
}

bgelement quadMesh::intToBge(int n) const
{
  double hx,hy;
  bgelement ee;

  if((n < 0) || (n >= nbelem))
    {
      std::cerr << "Boundary element number " << n << " out of range, nelem = " 
		<< nelem << std::endl;

      throw("Boundary element number out of range");
    }

  hx = (xmax-xmin) / nex;
  hy = (ymax-ymin) / ney;

// bcType[4] sets the boundary code
//
// Bottom boundary: type = bcType[0] 
// Right hand side: type = bcType[1]
// Top boundary:    type = bcType[2]
// Left hand side:  type = bcType[3]

 if(n < nex)  // bottom
  {
    ee.type = bcType[0];

    ee.nodes[0] = n;
    ee.nodes[1] = ee.nodes[0] + 1;

    ee.xx[0] = svector<2>(xmin + n * hx, ymin);
    ee.xx[1] = ee.xx[0] + svector<2>(hx  , 0.0);
  }
 else if(n < nex+ney) // rhs
   {
    ee.type = bcType[1];

    ee.nodes[0] = ((n-nex)+1) * (nex+1) - 1;
    ee.nodes[1] = ee.nodes[0] + (nex+1);

    ee.xx[0] = svector<2>(xmax, xmin + (n-nex)*hy);
    ee.xx[1] = ee.xx[0] + svector<2>(0.0, hy);
  }
 else if(n < 2*nex+ney)  // top
  {
    ee.type = bcType[2];

    ee.nodes[1] = (nex+1)*ney + (n-nex-ney);
    ee.nodes[0] = ee.nodes[1] + 1;

    ee.xx[1] = svector<2>(xmin + (n-nex-ney) * hx, ymax);
    ee.xx[0] = ee.xx[1] + svector<2>(hx  , 0.0);
  }
 else if(n < 2*nex + 2*ney) // lhs
   {
    ee.type = bcType[3];

    ee.nodes[1] = (n-2*nex-ney) * (nex+1);
    ee.nodes[0] = ee.nodes[1] +   (nex+1);

    ee.xx[1] = svector<2>(xmin, ymin + (n-2*nex-ney)*hy);
    ee.xx[0] = ee.xx[1] + svector<2>(0.0, hy);
  }
 else
   {
     std::cerr 
       << "Boundary element number " << n << " out of range, nbelem = " 
       << nbelem << std::endl;

     throw("Boundary element number out of range");
   }

 if( (bcType[0] == 0) && (bcType[1] == 0) &&
     (bcType[2] == 0) && (bcType[3] == 0) && (n==0) ) ee.type += 10;

 ee.master = gBasis->boundary->master;
 ee.basis  = gBasis->boundary->basis;

 return(ee);
}

#endif
*/


bgelement quadMesh::intToBge(int n) const
{
  double hx,hy;
  bgelement ee;

  //*********** Set up for Quad 4 Elements **************
  //*********** Set up for Quad 4 Elements **************
  //*********** Set up for Quad 4 Elements **************

  if((n < 0) || (n >= nbelem))
    {
      std::cerr << "Boundary element number " << n << " out of range, nelem = " 
		<< nelem << std::endl;

      throw("Boundary element number out of range");
    }

  hx = (xmax-xmin) / nex;
  hy = (ymax-ymin) / ney;

  //  std::cout<<"hx = " << hx << "hy = " << hy<<std::endl;
  // bcType[4] sets the boundary code
  //
  // Bottom boundary: type = bcType[0] 
  // Right hand side: type = bcType[1]
  // Top boundary:    type = bcType[2]
  // Left hand side:  type = bcType[3]

 if(n < nex)  // bottom
  {
    ee.type = bcType[0];

    ee.nodes[0] = n;
    ee.nodes[1] = ee.nodes[0] + 1;

    ee.xx[0] = svector<2>(xmin + n * hx, ymin);
    ee.xx[1] = ee.xx[0] + svector<2>(hx  , 0.0);
  }
 else if(n < nex+ney) // rhs
   {
    ee.type = bcType[1];

    ee.nodes[0] = ((n-nex)+1) * (nex+1) - 1;
    ee.nodes[1] = ee.nodes[0] + (nex+1);

    ee.xx[0] = svector<2>(xmax, ymin + (n-nex)*hy);
    ee.xx[1] = ee.xx[0] + svector<2>(0.0, hy);

    // std::cout<<"xmax = "<<xmax<<" xmin = "<<xmin<<std::endl;
    // std::cout<<"n = "<<n<<std::endl;
    // std::cout<<"ee.nodes[0]"<<ee.nodes[0]<<std::endl;
    // std::cout<<"ee.nodes[1]"<<ee.nodes[1]<<std::endl;
    // std::cout<<"ee.xx[0]"<<ee.xx[0]<<std::endl;
    // std::cout<<"ee.xx[1]"<<ee.xx[1]<<std::endl;
  }
 else if(n < 2*nex+ney)  // top
  {
    ee.type = bcType[2];

    ee.nodes[1] = (nex+1)*ney + (n-nex-ney);
    ee.nodes[0] = ee.nodes[1] + 1;

    ee.xx[1] = svector<2>(xmin + (n-nex-ney) * hx, ymax);
    ee.xx[0] = ee.xx[1] + svector<2>(hx  , 0.0);
  }
 else if(n < 2*nex + 2*ney) // lhs
   {
    ee.type = bcType[3];

    ee.nodes[1] = (n-2*nex-ney) * (nex+1);
    ee.nodes[0] = ee.nodes[1] +   (nex+1);

    ee.xx[1] = svector<2>(xmin, ymin + (n-2*nex-ney)*hy);
    ee.xx[0] = ee.xx[1] + svector<2>(0.0, hy);
  }
 else
   {
     std::cerr 
       << "Boundary element number " << n << " out of range, nbelem = " 
       << nbelem << std::endl;

     throw("Boundary element number out of range");
   }

 if( (bcType[0] == 0) && (bcType[1] == 0) &&
     (bcType[2] == 0) && (bcType[3] == 0) && (n==0) ) ee.type += 10;

 ee.master = gBasis->boundary->master;
 ee.basis  = gBasis->boundary->basis;
 ee.tag = n;
 return(ee);
}



void quadMesh::print()
{
  int i;

  gelement ge;
  bgelement be;

  std::ofstream dump("output");

  dump.setf(std::ios::fixed, std::ios::floatfield);

  dump << "Number of elements = " << nelem << std::endl;
  dump << "Number of nodes    = " << nodes << std::endl;
  dump << "Number of bandth = " << bandth << std::endl;
  
  for(i = 0; i < nelem; i++)
    {
      dump << "Element " << i << std::endl;

      ge = intToGe(i);

      // ge.print(dump);
    }

  dump << std::endl << "Number of boundary elements = " << nbelem << std::endl;

  for(i = 0; i < nbelem; i++)
    {
      dump << "Boundary Element " << i << ",";

      be = intToBge(i);

      // be.print(dump);
    }

  dump.close();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
#ifdef MYQUAD
quadMesh::quadMesh(std::string ss) : gBasis(Square4)
{
  int nx;
  int ny;

  if(!ss.empty())
    {
      std::stringstream ist(ss);
     
     if( !(ist >> nx) )
       {
	 std::cerr << "quadMesh(string): Invalid integer: " << ss
		   << ", terminating" << std::endl;
	 throw("Invalid input");
       }
     
     if(nx <= 0) 
       {
	 std::cerr << "quadMesh(string): Invalid integer: " 
		   << nx << ", terminating" << std::endl;
	 throw("Invalid input argument");
       }
     ny = nx;
    }
 else
   {
     int choice;
     std::cout<<"please specify size of domain in x direction:(1, 3, 30, 60, 90, 120, 150)"<<std::endl;
     std::cin>>choice;
     
     switch (choice)
       {

       case 1:{
	 std::cout<<"choice == 1 for domain setting"<<std::endl;
	 nx = 1;
	 ny = 1;
	 bandth = 1;
	 xmin = -1.0;xmax=1.0;
	 ymin = -1.0;ymax = 1.0;
	 break;
       }
       case 3:{
	 std::cout<<"choice == 3 for domain setting"<<std::endl;
	 nx = 3;
	 ny = 3;
	 bandth = 1;
	 xmin = -3.0;xmax=3.0;
	 ymin = -3.0;ymax = 3.0;
	 break;
       }
       case 30:{
	 std::cout<<"choice == 30 for domain setting"<<std::endl;
	 nx =151; 
	 ny = 151;
	 bandth = 21;
     
	 xmin = -150; xmax = 150;
	 ymin = -150; ymax = 150;
	 break;
       }

       case 60:{
	 std::cout<<"choice == 60 for domain setting"<<std::endl;
	 nx =301; 
	 ny = 41;
	 bandth = 19;
     
	 xmin = -30.0; xmax = 30.0;
	 ymin = -22.5; ymax = 22.5;
	 break;
       }
       case 90:{
	 std::cout<<"choice == 90 for domain setting"<<std::endl;
	 nx =451; 
	 ny = 61;
	 bandth = 27;
     
	 xmin = -45.0; xmax = 45.0;
	 ymin = -33.75; ymax = 33.75;
	 break;
       }
       case 120:{
	 std::cout<<"choice == 120 for domain setting"<<std::endl;
	 nx = 601; 
	 ny = 81;
	 bandth = 35;
     
	 xmin = -60.0; xmax = 60.0;
	 ymin = -45; ymax = 45;
	 break;
       }
       case 150:{
	 std::cout<<"choice == 150 for domain setting"<<std::endl;
	 nx = 751; 
	 ny = 101;
	 bandth = 43;
     
	 xmin = -75.0; xmax = 75.0;
	 ymin = -56.25; ymax = 56.25;
	 break;
       }
       default:{
	 std::cout<<"have you forgot to input domain size? or you have input number other than the choices"<<std::endl;
	 throw("invalid input");	 
       }
     }
   }

  nelem = nx*ny; 
  nodes = ((nx+1)*(nx+1));
  nbelem= 2*(nx + ny);

  nex = (nx);
  ney = (ny);

  //xmin = (-0.5); xmax = (0.5); 
  //ymin = (-0.5); ymax = (0.5);

  for(int i = 0; i < 4; i++) bcType[i] = 0;

  // bcType[2] = 1;   // Test Neumann BC

//  bcType[1] = 1;   // Test Neumann BC
//  std::cout << "quadMesh(): BC on RHS for wave equation test" << std::endl;

//   bcType[0] = bcType[2] = 1;   // Test Neumann BC
//   std::cout 
//	<< "quadMesh(): Neumann BC on Bottom and Top for wave equation test" 
//      	<< std::endl;

  // Mesh and BC for Examples in Chapters 8 and 9 of FEbook
  // using Example No. 1.

  //std::cout 
  // << "quadMesh(): Mesh and BC for Examples in Chapters 8 and 9 of FEbook" 
  // << std::endl;

  // xmin = (-0.5); xmax = (0.5); 
  // ymin = (-0.5); ymax = (0.5);

  // std::cout[botom,rhs,top,lhs] = ["
  //	    << bcType[0] << bcType[1] << bcType[2] << bcType[3]

  //xiaohan
  bcType[0] = 0;
  bcType[3] = 1;
  bcType[1] = 1;
  bcType[2] = 1;
  info();
}
#endif

*/
