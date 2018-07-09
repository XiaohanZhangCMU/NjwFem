// Set up a mesh with each quadralateral "affine" equivalent to the unit square

#include<string>
#include<sstream>
#include<fstream>

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
  int getNx() const {return nlayer;} //get no of layer elements
  int getNy() const {return ney;}
  int getLayer() const {return layer;}
  double getXmax() const{return xmax;}
  double getXmin() const{return xmin;}

  double getYmax() const{return ymax;}
  double getYmin() const{return ymin;}

  int getBndyElm() const{return nbelem;}
  void locate(const svector<2>&, gelement&, svector<2>&) const;

  void info();
  void print();
  void print(std::ostream&);
  void input(std::string meshFile);
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
  int nlayer;
  int layer;

  int nex, ney;
  int bandth;
  double xmin, xmax, ymin, ymax;
  int bcType[4];

  map<int,int*>* Elements;
  
  map<int,double*>* Nodes;
  
  map<int,int*>* Belems;
  
  vector<int>* Layers;

};
quadMesh::quadMesh(std::string meshFile) : gBasis(Square4)
{
  std::cout<<"Initializing quadMesh("<<meshFile<<"):"
	   <<"(read elements, nodes, belements)" << std::endl;

  input(meshFile);

  std::cout << "Number of vertices read in = ..................... " << nodes
	    << std::endl;
  std::cout << "Number of boundary segments = .................... " << nbelem 
	    << std::endl;
  std::cout << "Number of triangles = ............................ " << nelem 
	    << std::endl;

  return;
}

void quadMesh::input(std::string meshFile)
{
 std::ifstream in(meshFile.c_str());
  
  if( !in.is_open() ) 
    {
      std::cerr << "can't open file:" << meshFile  << std::endl;

      throw("Can't open file");
    }
  
  Elements = new map<int,int*>();
  map<int,int*>::iterator it_ele;
  
  Nodes = new map<int,double*>();
  map<int,double*>::iterator it_node;
  
  Belems = new map<int,int*>();
  map<int,int*>::iterator it_bele;
  
  Layers = new vector<int>();
  vector<int>::iterator it_layer;

  xmin = 1e50; xmax = -1e50;
  ymin = 1e50; ymax = -1e50;
 
  nelem = 0; 
  nodes = 0;
  nbelem = 0;
  nlayer = 0;

  bool readElem = false;
  bool readNode = false;
  bool readBelem = false;
  bool readLayer = false;


  std::string line;

  while(std::getline(in,line))
    {
      std::stringstream linestream(line);

      if(line.find("Elements") != string::npos && line.find("Boundary")==string::npos)
	{
	  readElem = true;
	  readNode = false;
	  readBelem = false;
	  readLayer = false;
	  std::cout<<"mark1"<<std::endl;
	  continue;
	}
      if(line.find("Nodes") != string::npos)
	{
	  readElem = false;
	  readNode = true;
	  readBelem = false;
	  readLayer = false;
	  std::cout<<"mark2"<<std::endl;
	  continue;
	}
      if(line.find("Boundary Elements") != string::npos)
	{
	  readElem = false;
	  readNode = false;
	  readBelem = true;
	  readLayer = false;
	  std::cout<<"mark3"<<std::endl;
	  continue;
	}
      if(line.find("Layer") != string::npos)
	{
	  readElem = false;
	  readNode = false;
	  readBelem = false;
	  readLayer = true;
	  std::cout<<"mark33"<<std::endl;
	  continue;
	}
      
      if(readElem)
	{
	  int tmp;
	  linestream>>tmp;
	  int* temp = new int[4];
	  for (int i = 0;i<4;i++)
	    {
	      linestream>>temp[i];
	      temp[i]--;
	    }
	  
	  (*Elements)[nelem] = temp;
	  nelem++;
	  std::cout<<"mark4"<<std::endl;
	}
      if(readNode)
	{
	  int tmp;
	  linestream>>tmp;
	  double* temp = new double[2];
	  for (int i = 0;i<2;i++)
	    {
	      if(i == 0)
		{
		  linestream>>temp[i];
		  if(temp[i] >= xmax)
		    {
		      xmax = temp[i];
		    }
		  if(temp[i] <=xmin)
		    {
		      xmin = temp[i];
		    }
		}
	      if(i == 1)
		{
		  linestream>>temp[i];
		  if(temp[i] >= ymax)
		    {
		      ymax = temp[i];
		    }
		  if(temp[i] <=ymin)
		    {
		      ymin = temp[i];
		    }
		}
	    }
	  
	  (*Nodes)[nodes] = temp;
	  nodes++;
	  std::cout<<"mark5"<<std::endl;
	}
      if(readBelem)
	{
	  int tmp;
	  linestream>>tmp;
	  int* temp = new int[3];
	  for(int i = 0;i<3;i++)
	    {
	      linestream>>temp[i];
	      if(i>0)
		temp[i]--;
	    }
	  (*Belems)[nbelem] = temp;
	  nbelem++;
	  std::cout<<"mark6"<<std::endl;
	}
      if(readLayer)
	{
	  int tmp;
	  linestream >>tmp;
	  tmp--;
	  (*Layers).push_back(tmp);
	  nlayer++;
	}
    }
  in.close();

  for(it_ele = (*Elements).begin();it_ele!=(*Elements).end();++it_ele)
    {
      std::cout<<it_ele->first<<" ";
      for(int i = 0;i<4;i++)
	std::cout<<(it_ele->second)[i]<<"  ";
      std::cout<<std::endl;
    }


  for(it_node = (*Nodes).begin();it_node!=(*Nodes).end();++it_node)
    {
      std::cout<<it_node->first<<" ";
      for(int i = 0;i<2;i++)
	std::cout<<(it_node->second)[i]<<"  ";
      std::cout<<std::endl;
    }


  for(it_bele = (*Belems).begin();it_bele!=(*Belems).end();++it_bele)
    {
      std::cout<<it_bele->first<<" ";
      for(int i = 0;i<3;i++)
	std::cout<<(it_bele->second)[i]<<"  ";
      std::cout<<std::endl;
    }

 for(it_layer = (*Layers).begin();it_layer!=(*Layers).end();++it_layer)
    {
      std::cout<<*it_layer;
      std::cout<<std::endl;
    }

  //xiaohan
  bcType[0] = 1;
  bcType[1] = 1;
  bcType[2] = 1;
  bcType[3] = 1;
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
  assert(n>=0 && n<nelem);
  
  gelement ee;
  
  for(int i = 0;i<4;i++)
    {
      ee.nodes[i] = (*Elements)[n][i];
      ee.xx[i] = svector<2>((*Nodes)[ee.nodes[i]][0], (*Nodes)[ee.nodes[i]][1]);
    }
  
  ee.master = gBasis->master;
  ee.basis = gBasis->basis;
  
  return(ee);
}

bgelement quadMesh::intToBge(int n) const
{
  bgelement ee;
  
  assert(n>=0 && n<nbelem);
  
  for(int i = 0;i<2;i++)
    {
      ee.nodes[i] = (*Belems)[n][i+1];
      ee.xx[i] = svector<2>((*Nodes)[ee.nodes[i]][0], (*Nodes)[ee.nodes[i]][1]);
    }

  ee.type = bcType[(*Belems)[n][0]];

  ee.master = gBasis->boundary->master;
  ee.basis  = gBasis->boundary->basis;

  if( (bcType[0] == 0) && (bcType[1] == 0) &&
      (bcType[2] == 0) && (bcType[3] == 0) && (n==0) ) ee.type += 10;

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
