// Set up a mesh with each simplex affine equivalent to the unit tetrahedra

#include<string>
#include<fstream>
#include<map>
#include<vector>
#include"intArray.h"

#include"../LinearAlgebra/svector.C"
#include"gelement.C"
#include"../FiniteElement/basisFns.C"

typedef GeometricElement<3, 4> gelement;
typedef BoundaryGelement<3, 3> bgelement;

class dir3Mesh
{
public:
  dir3Mesh(std::string);                        // read mesh.poly file

  dir3Mesh(const dir3Mesh &mm) 
  {
    std::cerr << "dir3Mesh(): Don't copy me, deep copy nyi\n";
    throw("don't copy a dir3Mesh");
  }

  dir3Mesh() 
  {
    std::cerr << "dir3Mesh(): Default constructor called, nyi\n";
    throw("bad dir3Mesh constructor");
  }
  
  ~dir3Mesh() { return; }
  
  int nNode()  const {return nodes;}
  int nElem()  const {return nelem;}
  int nBelem() const {return nbelem;}
  
  void print(std::ostream& dump);
  
  class geiterator : public std::vector<intArray<4> >::iterator
  {
  public:
    geiterator(const std::vector<intArray<4> >::iterator &it, dir3Mesh *g) 
      : std::vector<intArray<4> >::iterator(it), gm(g) {}

    gelement operator*();
    
  private:
    dir3Mesh * gm;
  };

  geiterator gebegin();
  geiterator geend();

  class bgeiterator : public std::vector<intArray<4> >::iterator
  {
  public:
    bgeiterator(const std::vector<intArray<4> >::iterator &it, dir3Mesh *g) 
      : std::vector<intArray<4> >::iterator(it), gm(g) {}

    bgelement operator*();
    
  private:
    dir3Mesh *gm;
  };

  bgeiterator bgebegin();
  bgeiterator bgeend();

  LagrangeBasisFunction<3> *gBasis;

private:
  int nelem;
  int nodes;
  int nbelem;
  
  std::map<int, svector<3> > vertices;    // vertices[id] = vertex coordinates

  std::vector<intArray<4> > triangles;   // boundary triangles + id

  std::vector<intArray<4> > tets;        // tetrahedra

  void input(std::string meshFile); // auxillary function to read .tplc file
};

dir3Mesh::dir3Mesh(std::string meshFile) : gBasis(Tet4)
{
  std::cout << "Initializing dir3Mesh(" << meshFile << "):" 
	    << " (reads .tplc file)" << std::endl;
  
  input(meshFile);
  
  std::cout << "Number of vertices read in = ..................... " << nodes
	    << std::endl;
  std::cout << "Number of boundary triangles = ................... " << nbelem 
	    << std::endl;
  std::cout << "Number of tetrahedra = ........................... " << nelem 
	    << std::endl;
  
  return;
}

gelement dir3Mesh::geiterator::operator*() 
{
  gelement ee;
  int j, *tt = std::vector<intArray<4> >::iterator::operator*();
  
  ee.type = 0;

  for(j = 0; j < 4; j++) 
    {
      ee.xx[j] = gm->vertices[tt[j]];
      ee.nodes[j] = tt[j];
    }

  ee.master = gm->gBasis->master;
  ee.basis  = gm->gBasis->basis;

  return ee;
}

dir3Mesh::geiterator dir3Mesh::gebegin() 
{
  return geiterator(tets.begin(), this);
}

dir3Mesh::geiterator dir3Mesh::geend() 
{
  return geiterator(tets.end(), this);
}

bgelement dir3Mesh::bgeiterator::operator*() 
{
  int j, *tt = std::vector<intArray<4> >::iterator::operator*();
  bgelement be;

  be.type = tt[3];

  for(j = 0; j < 3; j++) 
  {
      be.xx[j] = gm->vertices[tt[j]];
      
      be.nodes[j] = tt[j];
  }  

  be.master = gm->gBasis->boundary->master;
  be.basis  = gm->gBasis->boundary->basis;

  return be;
}

dir3Mesh::bgeiterator dir3Mesh::bgebegin() 
{
  return bgeiterator(triangles.begin(),this);
}

dir3Mesh::bgeiterator dir3Mesh::bgeend() 
{
  return bgeiterator(triangles.end(),this);
}

void dir3Mesh::input(std::string meshFile)
{
  int i, j, dim, id, i0,i1, i2, i3, nFace, nSegment, np, ntri, faceId;
  double xx, yy, zz;

  std::ifstream in((meshFile).c_str());

  if( !in.is_open() ) in.open((meshFile + ".tplc").c_str());

  if( !in.is_open() ) 
    {
      std::cerr << "dir3Mesh::input(): can't open .tplc file:" 
		<< meshFile << std::endl;

      throw("Can't open file");
    }

  in.ignore(256,'\n');    // first line is TPLC
  in.ignore(256,' ');     // Second line is DIMENSION dim
  in >> dim;
  assert(dim == 3);

  in.ignore(256,'\n');  
  in.ignore(256,' ');    // Next line is POINTS nodes
  in >> nodes;
  assert(nodes > 0);

  for(i = 0; i < nodes; i++)
    {
      in.ignore(256,'\n');

      in >> id;

      in >> xx >> yy >> zz;

      vertices[id] = svector<3>(xx,yy,zz);
    }

  // Read and discard segments

  in.ignore(256,'\n');
  in.ignore(256,' ');       // next line is SEGMENTS nSegment
  in >> nSegment;
  assert(nSegment >= 0);
  
  for(i = 0; i < nSegment; i++)
    {
    in.ignore(256,'\n');
   
    in >> id >> np;

    for(j = 0; j < np; j++) in >> i0; 
    }

  // Read the triangulated faces

  in.ignore(256,'\n');
  in.ignore(256,' ');       // next line is FACES nFace
  in >> nFace;
  assert(nFace >= 0);

  nbelem = 0;

  for(i = 0; i < nFace; i++)
    {
    in.ignore(256,'\n');

    in >> faceId >> ntri;

    assert(ntri >= 0);

    triangles.resize(nbelem + ntri);

    for(j = 0; j < ntri; j++)
      {
      in.ignore(256,'\n');

      in >> i0 >> i1 >> i2;

      assert(vertices.find(i0) != vertices.end());
      assert(vertices.find(i1) != vertices.end());
      assert(vertices.find(i2) != vertices.end());

      triangles[nbelem][0] = i0; 
      triangles[nbelem][1] = i1; 
      triangles[nbelem][2] = i2; 

      triangles[nbelem++][3] = faceId; 
      }
    }

  // Read the tetrahedra

  in.ignore(256,'\n');
  in.ignore(256,' ');       // next line is TETS nTet
  in >> nelem;
  assert(nelem >= 0);

  tets.resize(nelem);

  for(i = 0; i < nelem; i++)
    {
    in.ignore(256,'\n');

    in >> i0 >> i1 >> i2 >> i3;

    assert(vertices.find(i0) != vertices.end());
    assert(vertices.find(i1) != vertices.end());
    assert(vertices.find(i2) != vertices.end());
    assert(vertices.find(i3) != vertices.end());

    tets[i][0] = i0; 
    tets[i][1] = i1; 
    tets[i][2] = i2; 
    tets[i][3] = i3; 
    }

  in.close();

  return;
}
