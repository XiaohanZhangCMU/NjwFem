#include<iostream>
#include<cstdlib>
#include<cmath>

// If Nproc > 1 p-thread the banded Gaussian elimination routine solve()
// Set Nproc to the number of processors (2 for dual, 4 for quad etc.)

#define Nproc 1

class bandMV
{
public:
  class vector;
  class matrix;

  bandMV(void *) {}
  bandMV(int &argc, char ***argv) {}

  ~bandMV() {
    #if Nproc > 1
      cout << "~bandMV() p-thread banded matrix class with Nproc = "
	   << Nproc << endl;
    #else
      cout << "~bandMV() serial banded matrix class" << endl;
    #endif
  }

  class matrix
  {
  public:
    matrix(FeMdata mdata);
    ~matrix();
    
    void zero();
    void assemble(int, int*, double*);
    void assemble(int, int*, double*, bandMV::vector&, double*);
    void assemblepc(int, int*, double*);    

    void times(vector &x, vector &y);   // y = A x 

    void solve(vector&, vector&);    
    void reSolve(vector&, vector&);    
  private:
    int n, band;
    static int imin(int,int), imax(int,int);

    double *A;

    #if Nproc > 1

      static void* inner(void *);

      typedef struct
      {
	int i;
	int ii;
	int jp;
	int bband;
	int jkmax;
	int jpstep;
	double *f;
	double *a;
      } Idata;
    #endif    
  };

  class vector
  {
  public:
    vector(int size) : n(size) 
	{
	v = new double [size]; 

	if(v == 0)
	  {
	  cerr << "vector(int size): couldn't allocate memory" << endl;
	  throw("out of memory");
    	  }
	}
    ~vector() {delete [] v; }
    
    vector& zero();
    const int length() { return(n); }
    vector& operator=(const vector&);

    double normInfinity();
    vector& operator+=(vector&);
    vector& operator-=(vector&);
    
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

bandMV::vector& bandMV::vector::zero()  
{ 
  for(int i = 0; i < n; i++) v[i] = 0.0; 

  return(*this);
}

void bandMV::vector::assemble(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "bandMV::vector::assemble(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) v[index[i]] += fe[i];
}

void bandMV::vector::setvalues(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "bandMV::vector::setvalues(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) v[index[i]] = fe[i];
}

void bandMV::vector::extract(int nvar, int index[], double fe[]) const
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "bandMV::vector::extract(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) fe[i] = v[index[i]];
}

bandMV::vector& bandMV::vector::operator=(const vector &a)
{
  if(this != &a)
    {
    for(int i = 0; i < n; i++) v[i] = a.v[i];
    }

  return(*this);
}


bandMV::vector& bandMV::vector::operator+=(vector &a)
{
  for(int i = 0; i < n; i++) v[i] += a.v[i];

  return(*this);
}

bandMV::vector& bandMV::vector::operator-=(vector &a)
{
  for(int i = 0; i < n; i++) v[i] -= a.v[i];

  return(*this);
}

double bandMV::vector::normInfinity()
{
  double vmax = fabs(v[0]);

  for(int i = 1; i < n; i++)
    if(fabs(v[i]) > vmax) vmax = fabs(v[i]);

  return(vmax);
}

bandMV::matrix::matrix(FeMdata mdata)
  : n(mdata.nvar), band(mdata.band)
{
  int size = n*(2*band-1), nbyte = size*8, check = nbyte/8;

  if(size != check)
    {
      cout << "bandMV::matrix::matrix(): Warning 32 bit calculations of "
	   << "memory may fail" << endl;
    }

  A = new double[n*(2*band-1)];

  if(A == 0)
    {
    cerr << "Couldn't allocate banded matrix: n = " << n
	 << ", band = " << band << endl;
    throw("Couldn't allocate memory for matrix");
    }
}

bandMV::matrix::~matrix() {delete [] A;}

void bandMV::matrix::assemblepc(int nvar, int index[], double *a)
{
  return;    // preconditioner assemble for itterative schemes
}

void bandMV::matrix::assemble(int nvar, int index[], double *a)
{
  int i,j, ii;

#ifdef DeBug
  for(i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "bandMV::matrix::assemble(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }

    for(j = i+1; j < nvar; j++)
      if( abs(index[i]-index[j]) >= band )
	{
	cerr << "bandMV::matrix::assemble(): Index out of range: index[" 
	     << i << "] - index[" 
	     << j << "] = " << index[i]-index[j] << ", band = "
	     << band << endl;
	  
	throw("Index out of range");
	}
    }
#endif
  
  for(i = 0; i < nvar; i++)
    {
    ii = index[i]*(2*band-1) - index[i];
      
    for(j = 0; j < nvar; j++) A[ii+index[j]] += a[i*nvar+j];
    }
}

void bandMV::matrix::assemble(int nvar, int index[], double *ae, 
			  bandMV::vector& f, double *fe)
{
  int i,j, ii;

#ifdef DeBug
  for(i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "bandMV::matrix::assemble(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }

    for(j = i+1; j < nvar; j++)
      if( abs(index[i]-index[j]) >= band )
	{
	cerr << "bandMV::matrix::assemble(): Index out of range: index[" 
	     << i << "] - index[" 
	     << j << "] = " << index[i]-index[j] << ", band = "
	     << band << endl;
	  
	throw("Index out of range");
	}
    }
#endif
  
  for(i = 0; i < nvar; i++)
    {
    ii = index[i]*(2*band-1) - index[i];
      
    for(j = 0; j < nvar; j++) A[ii+index[j]] += ae[i*nvar+j];

    f.v[index[i]] += fe[i];
    }
}

void bandMV::matrix::zero()
{
  for(long i = 0; i < n*(2*band-1); i++) 
    {
      A[i] = 0.0;
    }
}

void bandMV::matrix::times(bandMV::vector &x, bandMV::vector &y) 
{
  int i,j,ii, jmin, jmax;    /*  Variables labeled 0,1,2, ... ,n-1 */ 
  double temp;               /* A[i,j] = A[i*(2*band-1) + j-i ]     */

  for(i = 0; i < n; i++)		
    {
    ii = i * (2*band-1);

    jmax = imin(n,i+band);
    jmin = imax(0,i-band+1);

    temp = 0.0;

    for(j = jmin; j < jmax; j++) temp += A[ii + j-i] * x.v[j];

    y.v[i] = temp;
    }

  return;
}

int bandMV::matrix::imin(int i, int j)
  {
  if(i < j) return(i);
  else      return(j);
  }

int bandMV::matrix::imax(int i, int j)
  {
  if(i > j) return(i);
  else      return(j);
  }

#if Nproc > 1  // use p-threaded solve()

#include <errno.h>
#include<pthread.h>

void bandMV::matrix::solve(bandMV::vector& uu, bandMV::vector& ff) 
{ 
 int j,jp;      /* Split the elimination over several processors */
 double temp;   /* Variables labeled 0,1,2, ... ,n-1  */
                /* A[i,j] = a[i*(2*band-1) + j-i ]     */

 double *u = uu.v, *f = ff.v;

 volatile Idata tdata[Nproc];
 int i,ii,jkmax,jpstep = band / Nproc + 1;

 int check = 0;
 void *retval;
 pthread_t threads[Nproc];

 errno = 0;

 for(i = 0; i < n-1; i++)		/* Gauss Eliminate */
   {
   ii = i * (2*band-1);
   jkmax = imin(n,i+band);

   if( A[ii] == 0.0 ) goto error;

   j = 0;

   for(jp = i+1; jp < jkmax; jp += jpstep) 
     {
       tdata[j].i  = i;
       tdata[j].ii = ii;
       tdata[j].jp = jp;
       tdata[j].bband = band;
       tdata[j].jkmax = jkmax;
       tdata[j].jpstep = jpstep;
       tdata[j].f = f;
       tdata[j].a = A;

       // inner((void *) &(tdata[j]));

       check |= pthread_create(threads+j, NULL, inner, (void *) &(tdata[j]));
       j++;
     }

   if(check)
     {
       cerr << "Can't make threads in bandMatrixp::solve() -- terminating\n";
       perror("bandMV::matrix::solve(): error in pthread_create()");       

       assert(false);
     }

   j = 0;

   for(jp = i+1; jp < jkmax; jp += jpstep) 
     check |= pthread_join(threads[j++], &retval);

   if(check)
     {
       cerr << "pthread_join() failed in solve() -- terminating\n";
       perror("bandMV::matrix::solve(): error in pthread_join()");

       assert(false);
     }
   }

 ii = (n-1) * (2*band-1);
 if( A[ii] == 0.0 ) goto error;

 for(i = n-1; i >= 0; i--)	/* back substitute	*/
   {
   ii = i * (2*band-1);
   jkmax = imin(n,i+band);

   temp = f[i];

   for(j = i+1; j < jkmax; j++)
     {
     temp -= A[ii + j-i] * u[j];
     }
  
   u[i] = temp / A[ii];
   }

 /* if(errno) perror("bandMV::matrix::solve(): at completion of solve()"); */

 return;

 error:
  cerr << "bandMV::matrix::solve(): zero pivot -- terminating" << endl;
  assert(false);

 exit(1);
}

void* bandMV::matrix::inner(void *vtd)
{
  int i,ii,j,jj,k,jp;
  double temp,*a,*f;
  Idata *td = (Idata *) vtd;

  i  = td->i;
  ii = td->ii;
  jp = td->jp;
  f  = td->f;
  a  = td->a;

  for(j = jp; j < imin(jp+td->jpstep,td->jkmax); j++)
   {
   jj = j * (2*td->bband-1);
     
   temp = a[jj + i-j] / a[ii];
   if(temp == 0.0) continue;

   f[j] -= temp * f[i];

   for(k = i+1; k < td->jkmax; k++) a[jj+k-j] -= temp * a[ii+k-i]; 

   a[jj + i-j] = temp;
   }

 return(NULL);
}

#else  // use serial solve

#include<omp.h>

void bandMV::matrix::solve(bandMV::vector& u, bandMV::vector& f) 
{ 
  int i,j,k,ii,jj,jkmax;     /*  Variables labeled 0,1,2, ... ,n-1 */ 
  double temp;               /* A[i,j] = A[i*(2*band-1) + j-i ]     */

  for(i = 0; i < n-1; i++)		/* Gauss Eliminate */
    {
    ii = i * (2*band-1);
    jkmax = imin(n,i+band);

    if( A[ii] == 0.0 ) goto error;

    // #pragma omp parallel for num_threads(4) private(jj,k,temp) 

    for(j = i+1; j < jkmax; j++)
      {
      jj = j * (2*band-1);
     
      temp = A[jj + i-j] / A[ii];
      if(temp == 0.0) continue;

      f.v[j] -= temp * f.v[i];

      for(k = i+1; k < jkmax; k++) A[jj+k-j] -= temp * A[ii+k-i]; 

      A[jj + i-j] = temp;
      }
    }

  ii = (n-1) * (2*band-1);
  if( A[ii] == 0.0 ) goto error;

  for(i = n-1; i >= 0; i--)	/* back substitute	*/
    {
    ii = i * (2*band-1);
    jkmax = imin(n,i+band);

    temp = f.v[i];

    for(j = i+1; j < jkmax; j++)
      {
      temp -= A[ii + j-i] * u.v[j];
      }
  
    u.v[i] = temp / A[ii];
    }

  return;

 error:
  cerr << "Zero pivot in banded solve() -- Execution Terminates" << endl;
  throw("Zero pivot in banded solve()");
}
#endif  // serial solve

void bandMV::matrix::reSolve(bandMV::vector& u, bandMV::vector& f) 
{ 
  int i,j,ii,jj,jkmax;     /*  Variables labeled 0,1,2, ... ,n-1 */ 
  double temp;               /* A[i,j] = A[i*(2*band-1) + j-i ]     */

  for(i = 0; i < n-1; i++)		/* forward substitute */
    {
    ii = i * (2*band-1);
    jkmax = imin(n,i+band);

    for(j = i+1; j < jkmax; j++)
      {
      jj = j * (2*band-1);
     
      temp = A[jj + i-j];

      f.v[j] -= temp * f.v[i];
      }
    }

  for(i = n-1; i >= 0; i--)	/* back substitute	*/
    {
    ii = i * (2*band-1);
    jkmax = imin(n,i+band);

    temp = f.v[i];

    for(j = i+1; j < jkmax; j++)
      {
      temp -= A[ii + j-i] * u.v[j];
      }
  
    u.v[i] = temp / A[ii];
    }

  return;
}
