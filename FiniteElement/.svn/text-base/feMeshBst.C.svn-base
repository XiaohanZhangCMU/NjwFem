#include"bst.C"
#include"felement.C"

#define FeMeshVectorGrouping

class feMesh : public gmMesh
{
public:
  feMesh(vector<basisFunction<Ndim> *> &bases, 
		 vector<int> &mults, string gmeshArg);

  feMesh(basisFunction<Ndim> *basis, int mult, string gmeshArg);
  
  ~feMesh() { bstfree(bst); }
  
  int band();
  int ndof() {return(ndofm);}
  
  class feiterator : public gmMesh::geiterator
  {
  public:
    feiterator(const geiterator & it, feMesh * f);
    felement operator*();
    
  private:
    feMesh * fm;
  };
  
  feiterator febegin();
  feiterator feend();
  
  class bfeiterator : public gmMesh::bgeiterator
  {
  public:
    bfeiterator(const bgeiterator & it, feMesh * f);
    bfelement operator*();
    
  private:
    feMesh * fm;
  };
  
  bfeiterator bfebegin();
  bfeiterator bfeend();

  felement geToFe(const gelement &);

  void print();  
private:
  mutable bnode *bst;
  
  int ndofm;
  
  int key2dof(int * key);
  
  void mkFeMesh();

  bfelement bgeToBfe(const bgelement &);

  vector<int> multiplicities;
  vector<basisFunction<Ndim> *  >  basisFns;
  vector<basisFunction<Ndim-1> *> bBasisFns;

  friend ostream& operator<<(ostream& out, feMesh &dd);
};

feMesh::feMesh(basisFunction<Ndim> *basis, int mult, string gmeshArg)
  : gmMesh(gmeshArg), bst(NULL)
{

  basisFns.resize(1); basisFns[0] = basis;

  multiplicities.resize(1); multiplicities[0] = mult;

  mkFeMesh();
}

feMesh::feMesh(vector<basisFunction<Ndim> *> &bases, 
	       vector<int> &mults, string gmeshArg)
  : gmMesh(gmeshArg), bst(NULL), multiplicities(mults), basisFns(bases)
{
  assert(basisFns.size() == multiplicities.size());

  mkFeMesh();
}

void feMesh::mkFeMesh()
{
  unsigned int j;
  int k, nn, dim, ss, key[Nkey];

  assert(Nkey >= 4*(Ndim-1));  // **** Square & Cube **** FIX THIS **********
  assert(Nkey >= Ndim+1);      // **** Simplex ********** FIX THIS **********

  bBasisFns.resize(basisFns.size());

  masterCell<Ndim-1> *bmaster = basisFns[0]->boundary->master;

  for(j = 0; j < basisFns.size(); j++) 
    {
    bBasisFns[j] = basisFns[j]->boundary;

    assert(bmaster == bBasisFns[j]->master);
    }

  masterCell<Ndim> *master = basisFns[0]->master;

  for (gmMesh::geiterator geit = gebegin(); 
       geit != geend(); geit++ )
    {
    gelement ee = *geit;

    assert(ee.master == master);

    for(dim = 0; dim < Ndim+1; dim++) // allocate dof for subcells[dim]
      {
      nn = 0;   // number of dof allocated for subsimplices of dimension dim

      for(j = 0; j < basisFns.size(); j++) 
	nn += basisFns[j]->subDof[dim] * multiplicities[j];
      
      for(ss = 0; ss < ee.master->nsubs[dim]; ss++)  // sub-cells
	{
	for(k = 0; k < Nkey; k++) key[k] = -1;
	  
	for(k = 0; k < ee.master->subSize[dim]; k++) 
	  key[k] = ee.nodes[ ee.master->subs(dim,ss,k) ];
	
	bstadd(&bst, bstsort(key), nn);
	}
      }
    }

 ndofm = bstset(bst);        // Assign dof numbers

 //  cout << "Number of degrees of freedom = ................... " << ndofm
 //       << endl;
 
 //  cout << "Leaving constructor" << endl;
}

int feMesh::band()
{
  int bb = 0, dmin, dmax;

  for (feMesh::feiterator feit = febegin(); 
      feit != feend(); feit++ )
    {
    felement ee = *feit;

    dmin = dmax = ee.dof[0];
   
    for(int j = 0; j < ee.ndof(); j++)
      {
      if(dmin > ee.dof[j]) dmin = ee.dof[j];
      if(dmax < ee.dof[j]) dmax = ee.dof[j];
      }
    
    if(bb < (dmax-dmin+1)) bb = dmax-dmin+1;
    }
 
  return(bb);
}

int feMesh::key2dof(int * key) {
    return bstdat(&bst, bstsort(key));
}

felement feMesh::geToFe(const gelement &ge)
{
  unsigned int i;
  int dd, dim, ii, j, k, m, *dof0, ss, key[Nkey], idx[Ndim+1];

  felement fe(ge);
    
  fe.ndofm = 0;
  fe.functions = &basisFns;
  fe.multiplicities = &multiplicities;

  for(i = 0; i < basisFns.size(); i++)
    fe.ndofm += multiplicities[i] * basisFns[i]->ndof;

  fe.resize(fe.ndofm);

  // set dof
  
  for(int dim = 0; dim < Ndim+1; dim++) idx[dim] = 0;

  ii = 0;

  for(i = 0; i < basisFns.size();   i++)   // each basis function
  for(m = 0; m < multiplicities[i]; m++)   // each copy
    {
    dof0 = fe.dof + ii;  // pointer for permutation

    for(dim = 0; dim < Ndim+1; dim++)
      {
      for(ss = 0; ss < ge.master->nsubs[dim]; ss++)  // sub-cells
	{
	for(k = 0; k < Nkey; k++) key[k] = -1;

	for(k = 0; k < ge.master->subSize[dim]; k++) 
	  key[k] = ge.nodes[ ge.master->subs(dim,ss,k) ];

	dd = key2dof(key) + idx[dim];

	for(j = 0; j < basisFns[i]->subDof[dim]; j++) 
	  fe.dof[ii++] = dd + j;
	}

      idx[dim] += basisFns[i]->subDof[dim];
      }

    basisFns[i]->permuteDof(ge.nodes, dof0);
    }

  // Local ordering has multiple copies of each element chained
  // (this is convenient for the permuteDof())

  // The following groups the variables

  #ifdef FeMeshVectorGrouping

  assert(ii == fe.ndofm);

  ii = 0;
  dd = 0;

  vector<int> dof(fe.ndofm);

  for(i = 0; i < basisFns.size();         i++)   // each basis function
    {
    for(j = 0; j < basisFns[i]->ndof; j++)  
    for(m = 0; m < multiplicities[i]; m++)   // each copy
      {
	dof[ii++] = fe.dof[dd + m * basisFns[i]->ndof + j];
      }
    
    dd += basisFns[i]->ndof * multiplicities[i];
    }

  for(j = 0; j < fe.ndofm; j++) fe.dof[j] = dof[j];

  #endif

  return(fe);
}

bfelement feMesh::bgeToBfe(const bgelement &ge)
{
  unsigned int i;
  int dd, dim, ii, j, k, m, *dof0, ss, key[Nkey], idx[Ndim];

  bfelement be(ge);

  be.ndofm = 0;
  be.functions = &bBasisFns;
  be.multiplicities = &multiplicities;

  for(i = 0; i < bBasisFns.size();i++)
    be.ndofm += multiplicities[i] * bBasisFns[i]->ndof;

  be.resize(be.ndofm);

  // set the dof

  for(int dim = 0; dim < Ndim; dim++) idx[dim] = 0;

  ii = 0;

  for(i = 0; i < bBasisFns.size();  i++)   // each basis function
  for(m = 0; m < multiplicities[i]; m++)   // each copy
    {
    dof0 = be.dof + ii;  // pointer for permutation

    for(dim = 0; dim < Ndim; dim++)
      {
      for(ss = 0; ss < ge.master->nsubs[dim]; ss++)  // sub-cells
	{
	for(k = 0; k < Nkey; k++) key[k] = -1;

	for(k = 0; k < ge.master->subSize[dim]; k++) 
	  key[k] = ge.nodes[ ge.master->subs(dim,ss,k) ];

	dd = key2dof(key) + idx[dim];
	
	for(j = 0; j < bBasisFns[i]->subDof[dim]; j++) 
	  be.dof[ii++] = dd + j;
	}

      idx[dim] += bBasisFns[i]->subDof[dim];
      }

    bBasisFns[i]->permuteDof(ge.nodes, dof0);
    }

  // Local ordering has multiple copies of each element chained
  // (this is convenient for the permuteDof())

  // The following groups the variables

  #ifdef FeMeshVectorGrouping

  assert(ii == be.ndofm);

  ii = 0;
  dd = 0;

  vector<int> dof(be.ndofm);

  for(i = 0; i < bBasisFns.size(); i++)   // each basis function
    {
    for(j = 0; j < bBasisFns[i]->ndof; j++)  
    for(m = 0; m < multiplicities[i];  m++)   // each copy
      {
	dof[ii++] = be.dof[dd + m * bBasisFns[i]->ndof + j];
      }
    
    dd += bBasisFns[i]->ndof * multiplicities[i];
    }

  for(j = 0; j < be.ndofm; j++) be.dof[j] = dof[j];

  #endif

  return(be);
}

ostream& operator<<(ostream& out, feMesh &dd)
{
  out << "Finite Element Mesh:" << endl;

  return(out);
}

// Fe iterator

feMesh::feiterator::feiterator(const feMesh::geiterator & it, feMesh * f) 
  : gmMesh::geiterator(it)
{ 
    fm = f;
}
    
felement feMesh::feiterator::operator*() 
{
  return fm->geToFe(this->geiterator::operator*());
}

feMesh::feiterator feMesh::febegin() 
{
    return feiterator(gebegin(),this);
}

feMesh::feiterator feMesh::feend() 
{
    return feiterator(geend(),this);
}

// be iterator

feMesh::bfeiterator::bfeiterator(const feMesh::bgeiterator & it, feMesh * f) 
  : gmMesh::bgeiterator(it)
{ 
    fm = f;
}

bfelement feMesh::bfeiterator::operator*() 
{
  return fm->bgeToBfe(this->bgeiterator::operator*());
}

feMesh::bfeiterator feMesh::bfebegin() {
    return bfeiterator(bgebegin(),this);
}

feMesh::bfeiterator feMesh::bfeend() {
    return bfeiterator(bgeend(),this);
}

void feMesh::print()
{
  unsigned int i;
  int j, k;

  ofstream dump("output");

  dump.setf(ios::fixed,ios::floatfield);

  dump << "Number of elements  = " << nelem  << endl;
  dump << "Number of dof       = " << ndof() << endl;
  dump << "Number of basis fns = " << multiplicities.size() 
       << ", Multiplicities";

  for(i = 0; i < multiplicities.size(); i++) dump << ", " << multiplicities[i];

  dump << endl << endl;

  int ne = 0;

  for(feiterator feit = febegin(); feit != feend(); feit++ )
    {
    felement ee = *feit;

    dump << "Element " << ne++ << endl;

    int dd = 0;
   
    for(i = 0; i < multiplicities.size(); i++)
      {
      int mm = multiplicities[i];	     

      int nd = basisFns[i]->ndof;

      dump << "  Variable " << i << " has " << nd << " dof";

     if(mm > 1)
       {
       #ifdef FeMeshVectorGrouping
	 dump << " with copies interlaced" << endl;
       #endif
       }
     else
       dump << endl;

     for(j = 0; j < mm; j++)
       {
       dump << "\tCopy " << j << "'s dof";

       #ifdef FeMeshVectorGrouping
         for(k = 0; k < nd; k++) dump << ", " << ee.dof[dd + mm*k + j];
       #else
         for(k = 0; k < nd; k++) dump << ", " << ee.dof[dd+ j*nd + k];
       #endif

       dump << endl;
       }

     dd += nd * mm;
     }
     
     dump << endl;
   }

  dump << endl << "Number of boundary elements = " << nbelem << endl;

  ne = 0;

  for(feMesh::bfeiterator bfeit = bfebegin(); 
       bfeit != bfeend(); bfeit++) 
    {
    bfelement be = *bfeit;     

    dump << "Boundary Element " << ne++ << endl;

    int dd = 0;
   
    for(i = 0; i < multiplicities.size(); i++)
     {
     int mm = multiplicities[i];	     

     int nd = bBasisFns[i]->ndof;

     dump << "  Variable " << i << " has " << nd << " dof";

     if(nd == 0)
       {
       dump << endl;
       continue;
       }

     if(mm > 1)
       {
       #ifdef FeMeshVectorGrouping
	 dump << " with copies interlaced" << endl;
       #endif
       }
     else
       dump << endl;

     for(j = 0; j < mm; j++)
       {
       dump << "\tCopy " << j << "'s dof";

       #ifdef FeMeshVectorGrouping
         for(k = 0; k < nd; k++) dump << ", " << be.dof[dd + mm*k + j];
       #else
         for(k = 0; k < nd; k++) dump << ", " << be.dof[dd+ j*nd + k];
       #endif

       dump << endl;
       }

     dd += nd * mm;
     }
     
    dump << endl;
    }

  dump.close();
}
