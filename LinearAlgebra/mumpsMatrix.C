#include<iostream>
#include<cstdlib>
#include<cmath>
#include<set>
#include<vector>

#include "mpi.h"
#include "dmumps_c.h"

#define DeBug

class mumpsMV
{
public:
  class vector;
  class matrix;

  mumpsMV(void *) {}
  mumpsMV(int &argc, char ***argv)
  {
    int ierr, myid;

    ierr = MPI_Init(&argc, argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    std::cout<<"mumpsMV(): ierr= "<<ierr<<std::endl;
  }

  ~mumpsMV() 
  {
    int ierr;
    ierr = MPI_Finalize();

    std::cout << "~mumpsMV(): ierr= " <<ierr<< std::endl;
  }

  class matrix
  {
  public:
    matrix(FeMdata mdata);
    ~matrix();
    
    void zero();
    void assemble(int, int*, double*);
    void assemble(int, int*, double*, mumpsMV::vector&, double*);
    void assemblepc(int, int*, double*);    

    void times(vector &x, vector &y);   // y = A x 

    void solve(vector&, vector&);    
    void reSolve(vector&, vector&);    
  private:

    // Sparse Storage of matrix A:
    // 1) a[k] is entry in A[irow[k]-1, jidx[k]-1]  (Fortran indexing) 
    // 2) jidx[k]-1 = j for jcol[j] <= k < jcol[j+1]

    // Note: 
    // 1) jcol[] has size n+1
    // 2) number of non-zero entries = jcol[n]
    // 3) a[], irow[], and jidx[] have size jcol[n]

    int n;
    double *a;      // stores non-zero matrix entries
    double *dd;     // scaling factor for rows

    int *irow;     // row indicies INCREMENTED by 1 for Fortran
    int *jcol;     // index to nonzero entries stored in irow[] and a[]

    int *jidx;     // column indicies INCREMENTED by 1 for Fortran

  DMUMPS_STRUC_C id;
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
    
    //double norm();
    double normInfinity();
    vector& operator+=(vector&);
    vector& operator-=(vector&);
    vector& operator*(double );
    //xiaohan
    vector& operator*=(double scal);

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

mumpsMV::vector& mumpsMV::vector::zero()  
{ 
  for(int i = 0; i < n; i++) v[i] = 0.0; 

  return(*this);
}

void mumpsMV::vector::assemble(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "mumpsMV::vector::assemble(): Index out of range: index[" 
		<< i << "] = " 
		<< index[i] << ", n = " << n << std::endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) v[index[i]] += fe[i];
}

void mumpsMV::vector::setvalues(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "mumpsMV::vector::setvalues(): Index out of range: index[" 
		<< i << "] = " 
		<< index[i] << ", n = " << n << std::endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) v[index[i]] = fe[i];
}

void mumpsMV::vector::extract(int nvar, int index[], double fe[]) const
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "mumpsMV::vector::extract(): Index out of range: index[" 
		<< i << "] = " 
		<< index[i] << ", n = " << n << std::endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) fe[i] = v[index[i]];
}

mumpsMV::vector& mumpsMV::vector::operator=(const vector &a)
{
  if(this != &a)
    {
    for(int i = 0; i < n; i++) v[i] = a.v[i];
    }

  return(*this);
}


mumpsMV::vector& mumpsMV::vector::operator+=(vector &a)
{
  for(int i = 0; i < n; i++) v[i] += a.v[i];

  return(*this);
}

mumpsMV::vector& mumpsMV::vector::operator-=(vector &a)
{
  for(int i = 0; i < n; i++) v[i] -= a.v[i];

  return(*this);
}

mumpsMV::vector& mumpsMV::vector::operator*=(double scal)
{
  for(int i = 0;i<n;i++) v[i] *= scal;

  return(*this);
}
double mumpsMV::vector::normInfinity()
{
  double vmax = fabs(v[0]);

  for(int i = 1; i < n; i++)
    if(fabs(v[i]) > vmax) vmax = fabs(v[i]);

  return(vmax);
}

mumpsMV::matrix::matrix(FeMdata mdata) : n(mdata.nvar)
{
  // Set up the indexing arrays for the sparse representation
  // of a matrix used by mumps etc.

  feMesh *mm = mdata.mm;

  assert(n == mm->ndof());

  ::vector<set<int> > cols(n);
  
  for (feMesh::feiterator feit = mm->febegin(); feit != mm->feend(); feit++ )
    {
    felement ee = *feit;

    // dd will have dof sorted (improves efficiency)

    set<int> dd(ee.dof, ee.dof+ee.ndof());

    for(int i = 0; i < ee.ndof(); i++) 
      cols[ee.dof[i]].insert(dd.begin(), dd.end());
    }  

  // determine the number of non-zero's

  if( !(jcol = new int[n+1]) ) assert(false);
  
  jcol[0] = 0;

  for(int j = 0; j < n; j++) jcol[j+1] = jcol[j] + cols[j].size();

  if( !(irow = new int[jcol[n]]) ) assert(false);
  if( !(jidx = new int[jcol[n]]) ) assert(false);

  // sparse storage requires row indicies to be monotone assending
  // this code assumes that the sets cols[j] are sorted assending

  for(int j = 0; j < n; j++) 
    {
    set<int>::iterator cj = cols[j].begin();

    for(int i = 0; i < (int) cols[j].size(); i++) 
      {
      irow[jcol[j] + i] = *(cj++) + 1;
      jidx[jcol[j] + i] = j       + 1;
      }
    }

  if( !(dd = new double [n]      ) ) assert(false);
  if( !(a  = new double [jcol[n]]) ) assert(false);
      
  id.par=1; 
  id.sym=0;
  id.job=-1;
  id.comm_fortran = -987654;  // define USE_COMM_WORLD -987654

  dmumps_c(&id);

  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0) 
    {
    id.n  = n; 
    id.nz = jcol[n]; 
    id.irn= irow; 
    id.jcn= jidx;
    id.a = a;
    }

 // pivoting comsumes more memory - increase icntl[14-1] from default 20%

  if(mdata.pivot) id.cntl[14-1] = 50;

  return;
}

mumpsMV::matrix::~matrix() 
{
  delete [] a;
  delete [] dd;

  delete [] irow;
  delete [] jcol;
  delete [] jidx;

  id.job = -2;
 
  dmumps_c(&id); 
}

void mumpsMV::matrix::assemblepc(int nvar, int index[], double *a)
{
  return;    // preconditioner assemble for itterative schemes
}

void mumpsMV::matrix::assemble(int nvar, int index[], double *ae)
{
  int i,j;
  int k;

#ifdef DeBug
  for(i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
	std::cerr << "mumpsMV::matrix::assemble(): Index out of range: index[" 
		  << i << "] = " 
		  << index[i] << ", n = " << n << std::endl;

	throw("Index out of range");
      }
    }
#endif
  
  for(j = 0; j < nvar; j++)
  for(i = 0; i < nvar; i++)
    {
      int left  = jcol[ index[j]   ];   // find index[i] in irow[left,right)
      int right = jcol[ index[j]+1 ];
      
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
  
void mumpsMV::matrix::assemble(int nvar, int index[], double *ae, 
			  mumpsMV::vector& f, double *fe)
{
  int i,j;
  int k;

#ifdef DeBug
  for(i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      std::cerr << "mumpsMV::matrix::assemble(): Index out of range: index[" 
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
      int left  = jcol[ index[j]   ];   // find index[i] in irow[left,right)
      int right = jcol[ index[j]+1 ];
      
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

// void mumpsMV::matrix::setvalue(int row, int col, double val)
// {

// }

void mumpsMV::matrix::zero()
{
  for(int i = 0; i < jcol[n]; i++) a[i] = 0.0;

  return;
}

void mumpsMV::matrix::times(mumpsMV::vector &x, mumpsMV::vector &y) 
{
  int j;
  int k;

  assert(x.length() == n);
  assert(y.length() == n);

  y.zero();

  for(j = 0; j < n; j++)
    {  
    for(k = jcol[j]; k < jcol[j+1]; k++) y.v[irow[k]-1] += a[k] * x.v[j];
    }

  return;
}

void mumpsMV::matrix::solve(mumpsMV::vector& u, mumpsMV::vector& f) 
{
  assert(u.length() == n);
  assert(f.length() == n);

  // scale the rows

  int i, j, k;

  for(i = 0; i < n; i++) dd[i] = 1.0;

  for(j = 0; j < n; j++) 
    {
    for(k = jcol[j]; k < jcol[j+1]; k++) 
      {
	double absak = fabs(a[k]);

	i = irow[k]-1;

	if(absak > dd[i]) dd[i] = absak;
      }
    }

  for(j = 0; j < n; j++) 
    {
    u.v[j] = f.v[j] / dd[j];

    for(k = jcol[j]; k < jcol[j+1]; k++) a[k] /= dd[irow[k]-1];
    }

  #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

  id.ICNTL(1)=-1; 
  id.ICNTL(2)=-1; 
  id.ICNTL(3)=-1; 
  id.ICNTL(4)=0;

  id.job=6;
  id.rhs = u.v;    // right hand side

  dmumps_c(&id);

  if(id.info[0] == 0) return;

  if(id.info[0] < 0)
    {
      std::cerr << "mumpsMV::matrix::Solve(): Error in linear solver"
		<< ", id.info[0] = " << id.info[0] 
		<< ", id.info[1] = " << id.info[1] 
		<< ", aborting" 
		<< std::endl;

      assert(false);
    }

  #ifdef DeBug
  std::cerr << "mumpsMV::matrix::Solve(): Warning in linear solver,"
	    << " id.info[0] = " << id.info[0] 
	    <<  std::endl;
  #endif

 return;
}

void mumpsMV::matrix::reSolve(mumpsMV::vector& u, mumpsMV::vector& f) 
{ 
  assert(u.length() == n);
  assert(f.length() == n);

  for(int i = 0; i < n; i++) u.v[i] = f.v[i] / dd[i];

  id.job = 3;      // solve step
  id.rhs = u.v;    // right hand side

  dmumps_c(&id);

  if(id.info[0] == 0) return;

  if(id.info[0] < 0)
    {
      std::cerr << "mumpsMV::matrix::Solve(): Error in linear solver,"
		<< " id.info[0] = " << id.info[0] << ", aborting" 
		<< std::endl;

      assert(false);
    }

  #ifdef DeBug
  std::cerr << "mumpsMV::matrix::Solve(): Warning in linear solver,"
	    << " id.info[0] = " << id.info[0] 
	    <<  std::endl;
  #endif

 return;
}
