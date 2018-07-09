// Set up mesh with elements affine equivalent to the unit tetrahedra

#include<string>
#include<fstream>
#include "allh.h"

#include"../LinearAlgebra/svector.C"
#include"gelement.C"
#include"../FiniteElement/basisFns.C"

typedef GeometricElement<3, 4> gelement;
typedef BoundaryGelement<3, 3> bgelement;

class dir3Mesh
{
public:
  dir3Mesh(std::string); 

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
    
  ~dir3Mesh() { delete C; return; }
    
  void print(std::ostream& dump);  // FIXME: implement this
  void print(std::string fileName);  // FIXME: implement this
    
  //  void locate(const dvector&, gelement&, dvector&) const;

  class geiterator : public Mesh::siterator
  {
  public:
    geiterator(const Mesh::siterator &it, dir3Mesh *g) 
      : Mesh::siterator(it), gm(g) {}

    gelement operator*();
    
  private:
    dir3Mesh *gm;
  };
        
  geiterator gebegin();
  geiterator geend();
    
  class bgeiterator {
  public:
    bgeiterator();
    bgeiterator(const Mesh::siterator & s, PLC * C, int face);
    void operator++();
    void operator++(int);
    bgeiterator& operator=(const bgeiterator & it);
    bool operator==(const bgeiterator & it);
    bool operator!=(const bgeiterator & it);
    bgelement operator*();
    
  private:
    PLC * C;
    int face;
    Mesh::siterator sit;
  };
  
  bgeiterator bgebegin();
  bgeiterator bgeend();
    
  LagrangeBasisFunction<3> *gBasis;

private:
  PLC * C;
};

dir3Mesh::dir3Mesh(std::string meshFile) : gBasis(Tet4)
{
  std::cout << "Initializing dir3Mesh(" << meshFile << "):" 
	    << " (reads .plc file)\n";
  
  C = new PLC();
  // FIXME
  // set up configuration here
  // to supress unneeded output
  C->ReadAndRefine(meshFile);

  return;
}

gelement dir3Mesh::geiterator::operator*() 
{
  gelement ee;
  Simplex *s = Mesh::siterator::operator*();

  for(int j = 0; j < 4; j++) 
    {
      ee.xx[j] = svector<3>(gm->C->P->V[s->v[j]].xyz);	
      ee.nodes[j] = s->v[j];
    }
    
  ee.master = gm->gBasis->master;
  ee.basis  = gm->gBasis->basis;

  return ee;
}

dir3Mesh::geiterator dir3Mesh::gebegin() 
{
  return geiterator(C->M->sbegin(), this);
}

dir3Mesh::geiterator dir3Mesh::geend() 
{
  return geiterator(C->M->send(), this);
}

dir3Mesh::bgeiterator::bgeiterator() {
  C = NULL;
  face = -1;
}

dir3Mesh::bgeiterator::bgeiterator(const Mesh::siterator &s, 
				   PLC *complex, int f) 
{
  sit = s;
  C = complex;
  face = f;
}

void dir3Mesh::bgeiterator::operator++() {
  (*this)++;
}

void dir3Mesh::bgeiterator::operator++(int) 
{
  int found = 0;
  int reset = 0;
  while (found == 0) 
    {
      while (face < (int)C->M->children.size() && 
	     (face < 0 || C->M->children[face]->dim != 2) ) 
	{
	  face++;
	  reset = 1;
	    
	}
      if (face >= (int)C->M->children.size()) 
	{
	  return;
	}
      Mesh * m = C->M->children[face];
	
      if (reset == 1) {
	sit = m->sbegin();
      } else {
	sit++;
      }

      found = 1;
      while(sit != m->send() && 
	    m->Exterior(*sit) )
	{
	  sit++;
	}
      if (sit==m->send()) 
	{
	  // done with this face... move to the next one
	  found = 0;
	  reset=1;
	  face++;
	}
      if (face >= (int)C->M->children.size()) 
	{
	  return;
	}
    }
}

dir3Mesh::bgeiterator& dir3Mesh::bgeiterator::operator=(const bgeiterator & it)
{
  C = it.C;
  face = it.face;
  sit = it.sit;
  return *this;
}

bool dir3Mesh::bgeiterator::operator==(const bgeiterator & it) {
  if (C != it.C) {
    return false;
  }
  if (C == NULL) {
    return true;
  }
    
  // .end() iterator
  if (face == it.face && face >= (int)C->M->children.size()) {
    return true;
  }

  if (face == it.face) {
    if (sit == it.sit) {
      return true;
    }
  }
  return false;
}

bool dir3Mesh::bgeiterator::operator!=(const bgeiterator & it) {
  return !(*this == it);
}

bgelement dir3Mesh::bgeiterator::operator*() {

  bgelement be;
  Simplex * s = *sit;
  
  for(int j = 0; j < 3; j++) 
    {
      be.xx[j] = svector<3>(C->P->V[s->v[j]].xyz);	
      be.nodes[j] = s->v[j];
    }
  
  be.type = C->M->children[face]->id;
  
  be.master = &Triangle;
  be.basis  = &tri3;

  return be;
}

dir3Mesh::bgeiterator dir3Mesh::bgebegin() 
{
  Mesh::siterator it;
  bgeiterator bgeit(it,C,-1);
  bgeit++;
  return bgeit;
}

dir3Mesh::bgeiterator dir3Mesh::bgeend() {
  return bgeiterator(C->M->send(), C, C->M->children.size());
}
