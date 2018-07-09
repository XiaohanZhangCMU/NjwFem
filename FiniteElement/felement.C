#include<vector>

class felement : public gelement {
public:
  felement();
  felement(const felement &fe);
  felement(const gelement &ge, int nd);
  ~felement();

  void resize(int);
  int ndof() { return(ndofm);}
  vector<vector<svector<Ndim> > > pts();      // interpolation points         

  int* dof;
  vector<int> *multiplicities;                // multiplicities
  vector<basisFunction<Ndim> *> *functions;   // functions
                                          
  felement& operator=(const felement &fe);
private:
  int ndofm;

  friend class feMesh;
};

class bfelement : public bgelement {
public:
  bfelement(const bgelement & bge, int nd);
  bfelement(const bfelement &fe);
  ~bfelement();  
  
  void resize(int);
  int ndof() { return(ndofm); }
  vector<vector<svector<Ndim> > > pts();      // interpolation points         
  vector<vector<svector<Ndim> > > normals();  // normals at interpolation pts
  vector<vector<svector<Ndim> > > areas();    // areas   at interpolation pts

  // Note: should eliminate areas(), 
  // should have normals(&vector<vector<double> > = NULL) (is NULL O.K?)

  int *dof;
  vector<int> *multiplicities;                // multiplicities
  vector<basisFunction<Ndim-1> *> *functions; // functions

private: 
  int ndofm; 

  friend class feMesh;
};

// ******************** felement ******************************
// ******************** felement ******************************
// ******************** felement ******************************

felement::felement()
  : gelement(), dof(NULL), ndofm(0)
{}

felement::felement(const gelement &ge, int nd = 0) 
  : gelement(ge), dof(NULL), ndofm(nd)
{
  if(nd != 0) resize(nd);
}

felement::felement(const felement &fe) 
  : gelement(fe), dof(NULL), multiplicities(fe.multiplicities), 
    functions(fe.functions), ndofm(fe.ndofm)
{
  if(fe.dof != NULL) 
    {
      dof = new int[ndofm];

      assert(dof != NULL);

      for(int i = 0; i < ndofm; i++) dof[i] = fe.dof[i];
    }
}

felement& felement::operator=(const felement &fe) 
{
  if(this == &fe) return(*this);
    
  gelement::operator=(fe);

  multiplicities = fe.multiplicities;

  functions = fe.functions;

  ndofm = fe.ndofm;

  assert(ndofm >= 0);

  if(ndofm != 0) 
    {
      int *newdof = new int[ndofm];

      assert(newdof != NULL);

      for(int i = 0; i < ndofm; i++) newdof[i] = fe.dof[i];

      if(dof != NULL) delete [] dof;

      dof = newdof;
    }
  else
    dof = NULL;

  return(*this);
}

felement::~felement() {if(dof != NULL) delete [] dof;}

void felement::resize(int nd) 
{
  assert(nd >= 0);

  if(dof != NULL) delete [] dof;
    
  ndofm = nd;

  if(ndofm != 0)
    {
      dof = new int[ndofm];

      assert(dof != NULL);
    }
  else
    dof = NULL;

  return;
}

vector<vector<svector<Ndim> > > felement::pts()  // interpolation points
{
  int npts;

  vector<vector<svector<Ndim> > > pp(functions->size());

  for(unsigned int i = 0; i < pp.size(); i++)
    {
      npts = (*functions)[i]->npts;

      pp[i].resize(npts);

      for(int j = 0; j < npts; j++) 
	pp[i][j] = map((*functions)[i]->pts[j]);
    }
      
  return(pp);
}

// ******************** bfelement ******************************
// ******************** bfelement ******************************
// ******************** bfelement ******************************

bfelement::bfelement(const bgelement &ge, int nd = 0) 
  : bgelement(ge), dof(NULL), ndofm(nd)
{
  if(nd != 0) resize(nd);
}

bfelement::bfelement(const bfelement &fe) 
  : bgelement(fe), dof(NULL), functions(fe.functions), ndofm(fe.ndofm)
{
  if(fe.dof != NULL) 
    {
      dof = new int[ndofm];

      assert(dof != NULL);

      for(int i = 0; i < ndofm; i++) dof[i] = fe.dof[i];
    }
}

bfelement::~bfelement() { if(dof != NULL) delete [] dof; }

void bfelement::resize(int nd) 
{
  assert(nd >= 0);

  if(dof != NULL) delete [] dof;
    
  ndofm = nd;

  if(ndofm != 0)
    {
      dof = new int[ndofm];

      assert(dof != NULL);
    }
  else
    dof = NULL;

  return;
}

vector<vector<svector<Ndim> > > bfelement::pts()  // interpolation points
{
  int npts;

  vector<vector<svector<Ndim> > > pp(functions->size());

  for(unsigned int i = 0; i < pp.size(); i++)
    {
      npts = (*functions)[i]->npts;

      pp[i].resize(npts);

      for(int j = 0; j < npts; j++) 
	pp[i][j] = map((*functions)[i]->pts[j]);
    }
      
  return(pp);
}

vector<vector<svector<Ndim> > > bfelement::normals() 
{
  int npts;                 // normals at interpolation pts
  double det;
  sbmatrix<Ndim> dxds;

  vector<vector<svector<Ndim> > > nn(functions->size());

  for(unsigned int i = 0; i < nn.size(); i++)
    {
      npts = (*functions)[i]->npts;

      nn[i].resize(npts);

      for(int j = 0; j < npts; j++) 
	{
	  map((*functions)[i]->pts[j], dxds);

	  nn[i][j] = dxds.normal(det);
	}
    }
      
  return(nn);
}

vector<vector<svector<Ndim> > > bfelement::areas() 
{
  int npts;                 // area vectors at interpolation pts
  double det, marea;
  sbmatrix<Ndim> dxds;

  vector<vector<svector<Ndim> > > aa(functions->size());

  for(unsigned int i = 0; i < aa.size(); i++)
    {
      npts = (*functions)[i]->npts;

      marea = (*functions)[i]->master->volume;

      aa[i].resize(npts);

      for(int j = 0; j < npts; j++) 
	{
	  map((*functions)[i]->pts[j], dxds);

	  aa[i][j] = dxds.normal(det);
	  aa[i][j] = (marea * det) * aa[i][j];
	}
    }
      
  return(aa);
}
