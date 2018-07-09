
// Set up triangular from output from Triangle

#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include"intArray.h"

#include"../LinearAlgebra/svector.C"
#include"gelement.C"
#include"../FiniteElement/basisFns.C"

typedef GeometricElement<2, 3> gelement;
typedef BoundaryGelement<2, 2> bgelement;

class triangleMesh
{
public:
  triangleMesh(std::string);                        // read mesh.poly file
  triangleMesh(const triangleMesh &mm) 
  {
    std::cerr << "triangleMesh(): don't copy a triangleMesh\n";
    throw("copying a triangleMesh");
  }
    
  triangleMesh() 
  {
    std::cerr << "triangleMesh(): Default constructor called, nyi\n";
    throw("bad triangleMesh constructor");
  }

  ~triangleMesh() {return;}
  
  int nNode()  const {return nodes;}
  int nElem()  const {return nelem;}
  int nBelem() const {return nbelem;}

  void print();

  bool checkBdyOrientation(svector<2> cc);
    
  //  void locate(const svector<2>&, gelement&, svector<2>&) const;

  class geiterator : public std::vector<intArray<3> >::iterator
  {
  public:
    geiterator(const std::vector<intArray<3> >::iterator &it, triangleMesh *g) 
      : std::vector<intArray<3> >::iterator(it), gm(g) {}

    gelement operator*();
    
  private:
    triangleMesh * gm;
  };

  geiterator gebegin();
  geiterator geend();

  class bgeiterator : public std::vector<intArray<3> >::iterator
  {
  public:
    bgeiterator(const std::vector<intArray<3> >::iterator &it, triangleMesh *g)
      : std::vector<intArray<3> >::iterator(it), gm(g) {}

    bgelement operator*();
    
  private:
    triangleMesh * gm;
  };
  
  bgeiterator bgebegin();
  bgeiterator bgeend();
  
  LagrangeBasisFunction<2> *gBasis;
  
private:
  int nelem;
  int nodes;
  int nbelem;

  std::vector<svector<2> > vertices;       // vertex coordinates

  std::vector<intArray<3> > segments;      // boundary edges + colors

  std::vector<intArray<3> > triangles;     // triangles

  void input(std::string meshFile);         // auxillary function to read data
};

triangleMesh::triangleMesh(std::string meshFile) : gBasis(Triangle3)
{
  std::cout << "Initializing triangleMesh(" << meshFile << "):" 
	    << " (reads .node, .poly, and .ele files)" << std::endl;

  input(meshFile);

  std::cout << "Number of vertices read in = ..................... " << nodes
	    << std::endl;
  std::cout << "Number of boundary segments = .................... " << nbelem 
	    << std::endl;
  std::cout << "Number of triangles = ............................ " << nelem 
	    << std::endl;

  return;
}

gelement triangleMesh::geiterator::operator*() 
{
  gelement ee;
  int j, *tt = std::vector<intArray<3> >::iterator::operator*();
  
  for(j = 0; j < 3; j++) 
    {
      ee.xx[j] = gm->vertices[tt[j]];
      ee.nodes[j] = tt[j];
    }
  
  ee.master = gm->gBasis->master;
  ee.basis  = gm->gBasis->basis;
  
  return ee;
}

triangleMesh::geiterator triangleMesh::gebegin() 
{
  return geiterator(triangles.begin(), this);
}

triangleMesh::geiterator triangleMesh::geend() 
{
  return geiterator(triangles.end(), this);
}

bgelement triangleMesh::bgeiterator::operator*() 
{
  bgelement be;
  int j, *tt = std::vector<intArray<3> >::iterator::operator*();

  be.type = tt[2];

  for(j = 0; j < 2; j++) 
  {
      be.xx[j] = gm->vertices[tt[j]];
      
      be.nodes[j] = tt[j];
  }  

  be.master = gm->gBasis->boundary->master;
  be.basis  = gm->gBasis->boundary->basis;

  return be;
}

triangleMesh::bgeiterator triangleMesh::bgebegin() 
{
  return bgeiterator(segments.begin(), this);
}

triangleMesh::bgeiterator triangleMesh::bgeend() 
{
  return bgeiterator(segments.end(), this);
}

void triangleMesh::print()
{
  int i;

  std::ofstream dump("outputGmesh");

  dump.setf(std::ios::fixed, std::ios::floatfield);

  dump << "Number of gelements = " << nelem << std::endl;
  dump << "Number of nodes     = " << nodes << std::endl;

  i = 0;

 for (geiterator geit = gebegin(); 
      geit != geend(); geit++ )
   {
     gelement ge = *geit;

     dump << "gElement " << i++ << ": nodes[] = ";

     for(int j = 0; j < 3; j++) dump << ge.nodes[j] << "  ";

      dump << "\n\t (" 
	   << ge.xx[0] << "," << ge.xx[1] << "," << ge.xx[2] 
	   << ")\n"; 
    }

  dump << "\nNumber of boundary gelements = " << nbelem << std::endl;

  i = 0;

  for (bgeiterator beit = bgebegin(); 
      beit != bgeend(); beit++ )
   {
     bgelement be = *beit;

     dump << "Boundary gElement " << i++ << ": nodes[] = ";

     for(int j = 0; j < 2; j++) dump << be.nodes[j] << "  ";

      dump << "\n\t (" 
	  << be.xx[0] << "," 
	  << be.xx[1] << "), type = " << be.type << std::endl;
    }

  dump.close();
}

void triangleMesh::input(std::string meshFile)
{
  int nAttr, nBmkr, color, dim, index, i0,i1, i2;
  double xx,yy;

  std::ifstream in((meshFile + ".node").c_str());

  if( !in.is_open() ) 
    {
      std::cerr << "can't open file:" << meshFile + ".node" << std::endl;

      throw("Can't open file");
    }

  // Read in the vertices

  in >> nodes >> dim >> nAttr >> nBmkr;

  assert(dim == 2);

  assert(nodes > 0);

  vertices.resize(nodes);

  for(int i = 0; i < nodes; i++)
    {
      in.ignore(256,'\n');

      in >> index;

      if(index != i) 
	std::cerr << "strToMesh(): Warning: data in .poly may be bad\n";

      in >> xx >> yy;

      vertices[i] = svector<2>(xx,yy);

      // ignore the rest
    }

  in.close();

  // Read in segments from the .poly file

  in.open((meshFile + ".poly").c_str());

  if( !in.is_open() ) 
    {
      std::cerr << "can't open file:" << meshFile + ".poly" << std::endl;

      throw("Can't open file");
    }

  in >> i0 >> dim >> nAttr >> nBmkr;

  assert(dim == 2);

  assert(i0 == 0);    // i0 = 0 means vertices in .node file
    
  in.ignore(256,'\n');
  
  in >> nbelem >> nBmkr;

  color = 0;      // default color if nBmkr = 0 (Dirichlet bc)

  segments.resize(nbelem);

  for(int i = 0; i < nbelem; i++)
    {
      in.ignore(256,'\n');

      in >> index;

      if(index != i) 
	std::cerr << "Warning: data in file may be bad" << std::endl;

      in >> i0 >> i1;

      assert((-1 < i0) && (i0 < nodes));
      assert((-1 < i1) && (i1 < nodes));

      segments[i][0] = i0;
      segments[i][1] = i1;

      if(nBmkr > 0) in >> color;

      segments[i][2] = color-1;   // triangle messes up zero labels
    }

  in.close();

  // Read in the triangles

  in.open( (meshFile + ".ele").c_str() );

  if( !in.is_open() ) 
    {
      std::cerr << "can't open file:" << meshFile + ".ele" << std::endl;

      throw("Can't open file");
    }

  in >> nelem >> i0 >> nAttr;

  triangles.resize(nelem);

  for(int i = 0; i < nelem; i++)
    {
      in.ignore(256,'\n');

      in >> index;

      if(index != i) 
	std::cerr << "Warning: data in file may be bad" << std::endl;

      in >> i0 >> i1 >> i2;

      assert((-1 < i0) && (i0 < nodes));
      assert((-1 < i1) && (i1 < nodes));
      assert((-1 < i2) && (i2 < nodes));

      triangles[i][0] = i0;
      triangles[i][1] = i1;
      triangles[i][2] = i2;

      // ignore the rest
    }

  in.close();

  checkBdyOrientation(svector<2>(0.375,0.25));

  return;
}

bool triangleMesh::checkBdyOrientation(svector<2> cc)
{
  // loop over boundary elements and check that they
  // go counterclockwise around the point cc

  // for this to work the domain must be star shaped wrt cc

  bool test = true;
  double cross;
  svector<2> v0, v1;

  for(int i = 0; i < nbelem; i++)
    {
      v0 = vertices[ segments[i][0]  ] - cc;
      v1 = vertices[ segments[i][1] ] - cc;

      cross = v0[0]*v1[1] - v0[1]*v1[0];

      if(cross < 0)   // fix thes
	{
	  test = false;

	  int temp = segments[i][0];

	  segments[i][0] = segments[i][1];
	  segments[i][1] = temp;
	}   
    }

  if(!test)
    std::cerr << "triangleMesh::checkBdyOrientation" << cc << ": "
	      << "boundary orientation test failed; fixup attempted\n";

  return(test);
}

