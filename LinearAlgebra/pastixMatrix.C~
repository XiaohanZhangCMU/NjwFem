#include<iostream>
#include<cstdlib>
#include<cmath>
#include<set>
#include<vector>
#include<complex>

extern "C"
{
#include "pastix.h"
}

const int Nproc = 48;

class pastixMV
{
public:
  class vector;
  class matrix;

  pastixMV(void *) {}
  pastixMV(int &argc, char ***argv)  {}

  ~pastixMV() 
  {
    std::cout << "~pastixMV(): matrix class with Nproc = " 
	      << Nproc << std::endl;
  }

  class matrix
  {
  public:
    matrix(FeMdata mdata);
    ~matrix();
    
    void zero();
    void assemble(int, int*, double*);
    void assemble(int, int*, double*, pastixMV::vector&, double*);
    void assemblepc(int, int*, double*);    

    void times(vector &x, vector &y);   // y = A x 

    void solve(vector&, vector&);    
    void reSolve(vector&, vector&);    
  private:

    // Sparse Storage of matrix A:
    // If jcol[j]-1 <= k < jcol[j+1]-1  
    // Then a[k] is the entry in A[irow[k]-1, j]  (Fortran indexing) 

    // Note: 
    // 1) jcol[] has size n+1
    // 2) number of non-zero entries = jcol[n]-1
    // 3) a[], irow[], and have size jcol[n]-1

    int n;
    double *a;      // stores non-zero matrix entries
    double *dd;     // scaling factor for rows

    int *irow;     // row indicies INCREMENTED by 1 for Fortran
    int *jcol;     // index to nonzero entries stored in irow[] and a[]

    pastix_data_t  *pastix_data;
    pastix_int_t   *perm;
    pastix_int_t   *invp;
    pastix_int_t    iparm[IPARM_SIZE];
    double          dparm[DPARM_SIZE];  
  };

  class vector
  {
  public:
    vector(int size) : n(size) 
    {
      v = new double [size]; 
      
      if(v == 0)
	{
	  std::cerr << "vector(int size): couldn't allocate memory" 
		    << std::endl;
	  throw("out of memory");
	}
    }
    
    ~vector() {delete [] v; }
    
    vector& zero();
    const int length() { return(n); }
    vector& operator=(const vector&);
    
    double norm();
    double normInfinity();
    vector& operator+=(vector&);
    vector& operator-=(vector&);
    vector& operator*(double );
    vector& operator*=(double );
    vector& add(vector& u, double s);
    vector& add(vector& u, vector& d, double s);

    double operator[](int i) { return(v[i]); }
    
    void assemble(int, int*, double*);
    void extract (int, int*, double*) const;
    void setvalues(int, int*, double*);
    
    const double* buffer() {return(v);}
  private:
    int n;
    double *v;
    
    friend void matrix::solve(vector &, vector&);
    friend void matrix::reSolve(vector &, vector&);
    friend void matrix::times(vector &x, vector &y);   // y = A x 
    friend void matrix::assemble(int, int*, double *a, vector&, double*);
  };
};

pastixMV::vector& pastixMV::vector::zero()  
{ 
  for(int i = 0; i < n; i++) v[i] = 0.0; 

  return(*this);
}

void pastixMV::vector::assemble(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "pastixMV::vector::assemble(): Index out of range: index[" 
		<< i << "] = " 
		<< index[i] << ", n = " << n << std::endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) v[index[i]] += fe[i];
}

void pastixMV::vector::setvalues(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "pastixMV::vector::setvalues(): Index out of range: index[" 
		<< i << "] = " 
		<< index[i] << ", n = " << n << std::endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) v[index[i]] = fe[i];
}

void pastixMV::vector::extract(int nvar, int index[], double fe[]) const
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "pastixMV::vector::extract(): Index out of range: index[" 
		<< i << "] = " 
		<< index[i] << ", n = " << n << std::endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) fe[i] = v[index[i]];
}

pastixMV::vector& pastixMV::vector::operator=(const vector &a)
{
  if(this != &a)
    {
    for(int i = 0; i < n; i++) v[i] = a.v[i];
    }

  return(*this);
}


pastixMV::vector& pastixMV::vector::operator+=(vector &a)
{
  for(int i = 0; i < n; i++) v[i] += a.v[i];

  return(*this);
}

pastixMV::vector& pastixMV::vector::operator-=(vector &a)
{
  for(int i = 0; i < n; i++) v[i] -= a.v[i];

  return(*this);
}


pastixMV::vector& pastixMV::vector::add(vector &a, double s)
{
  for(int i = 0; i < n; i++) v[i] += a.v[i]*s;

  return(*this);
}


pastixMV::vector& pastixMV::vector::add(vector &a, vector &b, double s)
{
  for(int i = 0; i < n; i++) v[i] = a.v[i]+ b.v[i]*s;

  return(*this);
}

pastixMV::vector& pastixMV::vector::operator*=(double a)
{
  for(int i = 0; i < n; i++) v[i] *= a;

  return(*this);
}

double pastixMV::vector::normInfinity()
{
  double vmax = fabs(v[0]);

  for(int i = 1; i < n; i++)
    if(fabs(v[i]) > vmax) vmax = fabs(v[i]);

  return(vmax);
}

pastixMV::matrix::matrix(FeMdata mdata) 
  : n(mdata.nvar), pastix_data(NULL)
{
  // Set up the indexing arrays for the sparse representation
  // of a matrix used by pastix etc.

  feMesh *mm = mdata.mm;

  assert(n == mm->ndof());

  ::vector<set<int> > cols(n);
  
  for (feMesh::feiterator feit = mm->febegin(); feit != mm->feend(); feit++)
   { 
   felement ee = *feit;

   // dd will have dof sorted (improves efficiency)

   set<int> dd(ee.dof, ee.dof+ee.ndof());

   for(int i = 0; i < ee.ndof(); i++) 
     cols[ee.dof[i]].insert(dd.begin(), dd.end());
   }

  // determine the number of non-zero's

  if( !(jcol = new int[n+1]) ) assert(false);
  
  jcol[0] = 1;     // fortran indexing

  for(int j = 0; j < n; j++) jcol[j+1] = jcol[j] + cols[j].size();

  if( !(irow = new int[jcol[n]-1]) ) assert(false);

  // sparse storage requires row indicies to be monotone assending
  // this code assumes that the sets cols[j] are sorted assending


  for(int j = 0; j < n; j++) 
    {
    set<int>::iterator cj = cols[j].begin();

    for(int i = 0; i < (int) cols[j].size(); i++) 
      irow[jcol[j]-1 + i] = *(cj++) + 1;
    }

  if( !(dd = new double [n]      ) ) assert(false);
  if( !(a  = new double [jcol[n]-1]) ) assert(false);

  if( !(perm = new pastix_int_t [n]) ) assert(false);
  if( !(invp = new pastix_int_t [n]) ) assert(false);

  return;
}

pastixMV::matrix::~matrix() 
{
  delete [] perm;
  delete [] invp;

  delete [] a;
  delete [] dd;

  delete [] irow;
  delete [] jcol;

  iparm[IPARM_START_TASK] = API_TASK_CLEAN;
  iparm[IPARM_END_TASK]   = API_TASK_CLEAN;

  pastix(&pastix_data,0, n, jcol,irow,a, perm,invp, NULL, 1, iparm,dparm);

}

void pastixMV::matrix::assemblepc(int nvar, int index[], double *a)
{
  return;    // preconditioner assemble for itterative schemes
}

void pastixMV::matrix::assemble(int nvar, int index[], double *ae)
{
  int i,j;
  int k;

#ifdef DeBug
  for(i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
	std::cerr << "pastixMV::matrix::assemble(): Index out of range: index[" 
		  << i << "] = " 
		  << index[i] << ", n = " << n << std::endl;

	throw("Index out of range");
      }
    }
#endif
  
  for(j = 0; j < nvar; j++)
  for(i = 0; i < nvar; i++)
    {
      int left  = jcol[ index[j]   ]-1;   // find index[i] in irow[left,right)
      int right = jcol[ index[j]+1 ]-1;
      
      for(k = (left+right)/2; k != left; k = (left+right)/2)
	{
	  if(irow[k]-1 <= index[i]) left  = k;
	  else                      right = k;
	}
      
      if(irow[k]-1 == index[i]) a[k] += ae[i*nvar+j];
      else                      assert(false);
    }

  return;
}
  
void pastixMV::matrix::assemble(int nvar, int index[], double *ae, 
			  pastixMV::vector& f, double *fe)
{
  int i,j;
  int k;

#ifdef DeBug
  for(i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "pastixMV::matrix::assemble(): Index out of range: index[" 
		<< i << "] = " 
		<< index[i] << ", n = " << n << std::endl;

      throw("Index out of range");
      }
    }
#endif

  for(j = 0; j < nvar; j++)
    {
    f.v[index[j]] += fe[j];

    for(i = 0; i < nvar; i++)
      {
      int left  = jcol[ index[j]   ]-1;   // find index[i] in irow[left,right)
      int right = jcol[ index[j]+1 ]-1;
      
      for(k = (left+right)/2; k != left; k = (left+right)/2)
	{
	  if(irow[k]-1 <= index[i]) left  = k;
	  else                      right = k;
	}
      
      if(irow[k]-1 == index[i]) a[k] += ae[i*nvar+j];
      else                      assert(false);
      }
    }

  return;
}

void pastixMV::matrix::zero()
{
  for(int i = 0; i < jcol[n]-1; i++) a[i] = 0.0;

  return;
}

void pastixMV::matrix::times(pastixMV::vector &x, pastixMV::vector &y) 
{
  int j;
  int k;

  assert(x.length() == n);
  assert(y.length() == n);

  y.zero();

  for(j = 0; j < n; j++)
    {  
    for(k = jcol[j]-1; k < jcol[j+1]-1; k++) y.v[irow[k]-1] += a[k] * x.v[j];
    }

  return;
}

void pastixMV::matrix::solve(pastixMV::vector& u, pastixMV::vector& f) 
{
  assert(u.length() == n);
  assert(f.length() == n);

  // scale the rows

  int i, j, k;

  for(i = 0; i < n; i++) dd[i] = 1.0;

  for(j = 0; j < n; j++) 
    {
    for(k = jcol[j]-1; k < jcol[j+1]-1; k++) 
      {
	double absak = fabs(a[k]);

	i = irow[k]-1;

	if(absak > dd[i]) dd[i] = absak;
      }
    }

  for(j = 0; j < n; j++) 
    {
    u.v[j] = f.v[j] / dd[j];

    for(k = jcol[j]-1; k < jcol[j+1]-1; k++) a[k] /= dd[irow[k]-1];
    }

  // First time through do a complete solve
  // Subsequently, assume only the matrix ENTRIES and rhs and have changed,
  // If only the rhs has changed, use reSolve()

  if(pastix_data == NULL)
    {
    // set parameters

    iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    iparm[IPARM_START_TASK]       = API_TASK_INIT;
    iparm[IPARM_END_TASK]         = API_TASK_INIT;

    pastix(&pastix_data,0, n, jcol,irow,a, perm,invp, u.v, 1, iparm,dparm);

    iparm[IPARM_VERBOSE]       = API_VERBOSE_NOT;
    iparm[IPARM_THREAD_NBR]    = Nproc;
    iparm[IPARM_SYM]           = API_SYM_NO;
    iparm[IPARM_FACTORIZATION] = API_FACT_LU;

    iparm[IPARM_START_TASK]    = API_TASK_ORDERING;
    }
  else
    iparm[IPARM_START_TASK]    = API_TASK_NUMFACT;

  // solve the system

  iparm[IPARM_END_TASK]        = API_TASK_SOLVE;

  pastix(&pastix_data,0, n, jcol,irow,a, perm,invp, u.v, 1, iparm,dparm);

  if(iparm[IPARM_ERROR_NUMBER] != NO_ERR)
    cerr << "pastixMV::matrix::reSolve): pastix() error code = "
	 << iparm[IPARM_ERROR_NUMBER] << endl;

 return;
}

void pastixMV::matrix::reSolve(pastixMV::vector& u, pastixMV::vector& f) 
{ 
  assert(u.length() == n);
  assert(f.length() == n);

  for(int i = 0; i < n; i++) u.v[i] = f.v[i] / dd[i];

  // solve the system

  iparm[IPARM_START_TASK] = API_TASK_SOLVE;
  iparm[IPARM_END_TASK]   = API_TASK_SOLVE;

  pastix(&pastix_data,0, n, jcol,irow,a, perm,invp, u.v, 1, iparm,dparm);

  if(iparm[IPARM_ERROR_NUMBER] != NO_ERR)
    cerr << "pastixMV::matrix::reSolve): pastix() error code = "
	 << iparm[IPARM_ERROR_NUMBER] << endl;

 return;
}
