#include<fstream>
#include<iostream>
#include<cstdlib>
#include<cmath>

// #define MKL

#ifdef MKL
#include "mkl.h"
//#extern "C" void dgbsv(int *N,int *lband,int *uband,int *NRHS,double *AB,int *ldab,int *ipiv,double *y,int *ldb,int *info);

//#extern "C" void dgbtrs(char *c, int *N,int *lband,int *uband,int *NRHS,double *AB,int *ldab,int *ipiv,double *y,int *ldb,int *info);
#else
extern "C" void dgbsv_(int *N,int *lband,int *uband,int *NRHS,double *AB,int *ldab,int *ipiv,double *y,int *ldb,int *info);

extern "C" void dgbtrs_(char *c, int *N,int *lband,int *uband,int *NRHS,double *AB,int *ldab,int *ipiv,double *y,int *ldb,int *info);
#endif

// LaPack Banded Solver

#define DeBug

class lapackMV
{
public:
class vector;
class matrix;

lapackMV(void *) {}
lapackMV(int &argc, char ***argv) {}

~lapackMV() 
{
  #ifdef MKL
  cout << "~lapackMV() banded matrix class using dgbsv(), mkl version" << endl;
  #else
  cout << "~lapackMV() banded matrix class using dgbsv_()" << endl;
  #endif
}

class matrix
{
public:
matrix(FeMdata mdata);
~matrix();

void zero();
void assemble(int, int*, double*);
void assemble(int, int*, double*, lapackMV::vector&, double*);
void assemblepc(int, int*, double*);    

void times(vector &x, vector &y);   // y = A x 

void solve(vector&, vector&);    
void reSolve(vector&, vector&);    

void print(ostream& dump); // print matrix to ofstream
void print(const char[]);  // print sparse matrix to file in Matlab format.
private:
int n, ldab;
int band, *ipiv;

static int imin(int,int), imax(int,int);

double *A, *dd;    // dd is scale factor for rows
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

double norm();
double normInfinity();
vector& operator+=(vector&);
vector& operator-=(vector&);
vector& operator*(double );

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

lapackMV::vector& lapackMV::vector::zero()  
{ 
for(int i = 0; i < n; i++) v[i] = 0.0; 

return(*this);
}

void lapackMV::vector::assemble(int nvar, int index[], double fe[])
{
#ifdef DeBug
for(int i = 0; i < nvar; i++)
{
if((index[i] < 0) || (index[i] >= n))
{
cerr << "lapackMV::vector::assemble(): Index out of range: index[" 
   << i << "] = " 
   << index[i] << ", n = " << n << endl;

throw("Index out of range");
}
}
#endif

for(int i = 0; i < nvar; i++) v[index[i]] += fe[i];
}

void lapackMV::vector::setvalues(int nvar, int index[], double fe[])
{
#ifdef DeBug
for(int i = 0; i < nvar; i++)
{
if((index[i] < 0) || (index[i] >= n))
{
cerr << "lapackMV::vector::setvalues(): Index out of range: index[" 
   << i << "] = " 
   << index[i] << ", n = " << n << endl;

throw("Index out of range");
}
}
#endif

for(int i = 0; i < nvar; i++) v[index[i]] = fe[i];
}

void lapackMV::vector::extract(int nvar, int index[], double fe[]) const
{
#ifdef DeBug
for(int i = 0; i < nvar; i++)
{
if((index[i] < 0) || (index[i] >= n))
{
cerr << "lapackMV::vector::extract(): Index out of range: index[" 
   << i << "] = " 
   << index[i] << ", n = " << n << endl;

throw("Index out of range");
}
}
#endif

for(int i = 0; i < nvar; i++) fe[i] = v[index[i]];
}

lapackMV::vector& lapackMV::vector::operator=(const vector &a)
{
if(this != &a)
{
for(int i = 0; i < n; i++) v[i] = a.v[i];
}

return(*this);
}


lapackMV::vector& lapackMV::vector::operator+=(vector &a)
{
for(int i = 0; i < n; i++) v[i] += a.v[i];

return(*this);
}

lapackMV::vector& lapackMV::vector::operator-=(vector &a)
{
for(int i = 0; i < n; i++) v[i] -= a.v[i];

return(*this);
}

lapackMV::vector& lapackMV::vector::operator*(double c)
{
for(int i = 0; i < n; i++) v[i] = c * v[i];

return(*this);
}

double lapackMV::vector::normInfinity()
{
double vmax = fabs(v[0]);

for(int i = 1; i < n; i++)
if(fabs(v[i]) > vmax) vmax = fabs(v[i]);

return(vmax);
}

double lapackMV::vector::norm()
{
double vnorm = 0.0;

for(int i = 0; i < n; i++) vnorm += v[i]*v[i];

return(sqrt(vnorm));
}

lapackMV::matrix::matrix(FeMdata mdata)
: n(mdata.nvar) ,band(mdata.band)
{
ldab=3*band-2;

int size = n*ldab, nbyte = size*8, check = nbyte/8;

if(size != check)
{
cout << "lapackMV::matrix::matrix(): Warning 32 bit calculations of "
   << "memory may fail" << endl;
}

ipiv = new int[n];

if(ipiv == 0)
{
cerr << "lapackMV::matrix:matrix(): Couldn't allocate memory" << endl;
throw("out of memory");
}

dd = new double[n];    

if(dd == 0)
{
cerr << "lapackMV::matrix:matrix(): Couldn't allocate memory" << endl;
throw("out of memory");
}

A = new double[size];

if(A == 0)
{
cerr << "lapackMV::matrix:matrix(): "
 << "Couldn't allocate banded matrix: n = " << n
 << ", band = " << band << endl;

throw("Couldn't allocate memory for matrix");
}
}

lapackMV::matrix::~matrix() 
{
delete [] A; 
delete [] dd; 
delete [] ipiv;
}

void lapackMV::matrix::assemblepc(int nvar, int index[], double *a)
{
return;    // preconditioner assemble for itterative schemes
}

void lapackMV::matrix::assemble(int nvar, int index[], double *a)
{
int i,j;

#ifdef DeBug
for(i = 0; i < nvar; i++)
{
if((index[i] < 0) || (index[i] >= n))
{
cerr << "lapackMV::matrix::assemble(): Index out of range: index[" 
   << i << "] = " 
   << index[i] << ", n = " << n << endl;

throw("Index out of range");
}

for(j = i+1; j < nvar; j++)
if( abs(index[i]-index[j]) >= band )
{
cerr << "lapackMV::matrix::assemble(): Index out of range: index[" 
     << i << "] - index[" 
     << j << "] = " << index[i]-index[j] << ", band = "
     << band << endl;
  
throw("Index out of range");
}
}
#endif

for(i = 0; i < nvar; i++)
{

for(j = 0; j < nvar; j++)
A[index[j]*ldab+2*(band-1)-index[j]+index[i]] += a[i*nvar+j];
}
}

void lapackMV::matrix::assemble(int nvar, int index[], double *ae, 
		  lapackMV::vector& f, double *fe)
{
int i,j;

#ifdef DeBug
for(i = 0; i < nvar; i++)
{
if((index[i] < 0) || (index[i] >= n))
{
cerr << "lapackMV::matrix::assemble(): Index out of range: index[" 
   << i << "] = " 
   << index[i] << ", n = " << n << endl;

throw("Index out of range");
}

for(j = i+1; j < nvar; j++)
if( abs(index[i]-index[j]) >= band )
{
cerr << "lapackMV::matrix::assemble(): Index out of range: index[" 
     << i << "] - index[" 
     << j << "] = " << index[i]-index[j] << ", band = "
     << band << endl;
  
throw("Index out of range");
}
}
#endif

for(i = 0; i < nvar; i++)
{

for(j = 0; j < nvar; j++){
A[index[j]*ldab+2*(band-1)-index[j]+index[i]] += ae[i*nvar+j];
}

f.v[index[i]] += fe[i];
}
}

void lapackMV::matrix::zero()
{
for(long i = 0; i < n*ldab; i++) A[i] = 0.0;
}

void lapackMV::matrix::times(lapackMV::vector &x, lapackMV::vector &y) 
{
int i,j;    /*  Variables labeled 0,1,2, ... ,n-1 */ 
double temp;

for(i = 0; i < n; i++)		
{

temp = 0.0;

int jmax = imin(n,i+band);
int jmin = imax(0,i-band+1);

for(j = jmin; j < jmax; j++) temp += A[j*ldab+2*(band-1)+i-j] * x.v[j];

y.v[i] = temp;
}

return;
}

int lapackMV::matrix::imin(int i, int j)
{
if(i < j) return(i);
else      return(j);
}

int lapackMV::matrix::imax(int i, int j)
{
if(i > j) return(i);
else      return(j);
}

void lapackMV::matrix::solve(lapackMV::vector& u, lapackMV::vector& f) 
{ 
// solve A u = f

int nrhs=1,ldb=n,info = 0;
int kl=band-1,ku=band-1;

for(int i=0;i<n;i++) u.v[i]=f.v[i];

// before doing the solve, get rid of those 1e50 diagonal entries

for(int i = 0; i < n; i++) 
{
double diag=A[i*ldab+2*(band-1)];

dd[i] = 1.0;

if(fabs(diag)>1.e40)
{
dd[i] = diag;

u.v[i]=u.v[i]/diag;

for(int j=max(0,i-band+1);j<=min(n-1,i+band-1);j++)
A[j*ldab+2*(band-1)+i-j] = A[j*ldab+2*(band-1)+i-j]/diag;
}
}

#ifdef MKL
  dgbsv(&n,&kl,&ku,&nrhs,A,&ldab,ipiv,u.v,&ldb,&info);
#else
  dgbsv_(&n,&kl,&ku,&nrhs,A,&ldab,ipiv,u.v,&ldb,&info);
#endif
}

void lapackMV::matrix::reSolve(lapackMV::vector &u, lapackMV::vector &f) 
{ 
  // scale the rows of the rhs
  
  for(int i=0; i < n; i++) u.v[i] = f.v[i] / dd[i];

  // solve the system

  char c = 'N';

  int nrhs=1,ldb=n,info = 0;
  int kl=band-1,ku=band-1;

  #ifdef MKL
  dgbtrs(&c, &n,&kl,&ku,&nrhs,A,&ldab,ipiv,u.v,&ldb,&info);
  #else
  dgbtrs_(&c, &n,&kl,&ku,&nrhs,A,&ldab,ipiv,u.v,&ldb,&info);
  #endif
}


void lapackMV::matrix::print(ostream& dump){
  // print matrix to ofstream
      dump.precision(17);
  for(int i=0;i<n;i++){

    for(int j=max(0,i-band+1);j<=min(n-1,i+band-1);j++){
      if(A[j*ldab+2*(band-1)-j+i] != 0.0)
        dump << i+1 << " " << j+1 << " " << A[j*ldab+2*(band-1)-j+i] << endl;
    }
  }
}

void lapackMV::matrix::print(const char *filename){
  // print sparse matrix to file in Matlab format.
  ofstream dump(filename);
  dump.precision(17);
  dump << "zzz=[ ..." << endl;
  for(int i=0;i<n;i++){

    for(int j=max(0,i-band+1);j<=min(n-1,i+band-1);j++){
      if(A[j*ldab+2*(band-1)-j+i] != 0.0)
        dump << i+1 << " " << j+1 << " " << A[j*ldab+2*(band-1)-j+i] << endl;
    }
  }
  dump << "];" << endl;
  dump << "Mat=sparse(zzz(:,1),zzz(:,2),zzz(:,3));" << endl;
}
