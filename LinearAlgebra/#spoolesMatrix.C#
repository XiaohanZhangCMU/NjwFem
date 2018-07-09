// Spools interface
// Spools interface
// Spools interface

#include<iostream>
#include<cstdlib>
#include<cmath>

#define Nthread 1

#if (Nthread > 1)
  #include "BridgeMT.h"
#else
  #include "Bridge.h"
#endif

// To compile spools 2.2 files to run from a C++ application

// A) USING THE TAR FILE MODIFIED BY njwSpooles.2.2.tar.gz

// 1) tar -xvf njwSpooles.2.2.tar.gz   (this creates the Spooles directory)
// 2) cd Spooles
// 3) njwMake

// The spooles serial and MT library is made

// B) USING THE TAR FILE FROM NETLIB

// 1) Download and gunzip and tar -xvf the files into a Spooles directory

// WARNING: tar -xvf does not create a directory, it dumps into the current

// 2) Edit ./Make.inc and change the compiler to CC = g++
// 3) Change the cast (void *) to (void **) in the files
//    ./FrontMtx/src/instance.c, lines 346 and 390
//    ./SemiImplMtx/src/init.c,  lines 609 and 639
// 4) From the spooles directory "make lib" will make spooles.a
// 5) Make the following correction to ./LinSol/srcST/solve.c
//      The scope of if(msglvl > 1) { ... } on line 104 should
// 	be extended to contain the if(cputotal > 0.0) { ... }
//	statment which starts on line 109.
// 6) From the ./LinSol/srcST directory, "make Bridge.a"

// 7) To compile an application code let $(SPOOLES_DIR) be set, then

//  g++ -O files.C -I$(SPOOLES_DIR)/LinSol $(SPOOLES_DIR)/LinSol/srcST/Bridge.a $(SPOOLES_DIR)/spooles.a 

// WARNING: The ordering of the .a files is important!
// WARNING: The ordering of the .a files is important!
// WARNING: The ordering of the .a files is important!

// To use the pthreads additionally 

// 8) From the ./MT/src "make spoolesMT.a"
// 9) Make the following correction to ./LinSol/srcMT/solve.c
//      The scope of if(msglvl > 1) { ... } on line 105 should
// 	be extended to contain the if(cputotal > 0.0) { ... }
//	statment which starts on line 109.
//10) From the ./LinSol/srcMT directory, "make BridgeMT.a"

//11) To compile, set $(SPOOLES_DIR), then

// g++ -O files.C -I$(SPOOLES_DIR)/LinSol ${SPOOLES_DIR}/LinSol/srcMT/BridgeMT.a $(SPOOLES_DIR)/MT/src/spoolesMT.a $(SPOOLES_DIR)/spooles.a -lpthread

class spoolesMV
{
public:
  class vector;
  class matrix;
  class input;

  spoolesMV(void *) {}
  spoolesMV(int &argc, char ***argv) {}

  ~spoolesMV() {
    #if Nthread > 1
      cout << "~spoolesMV() p-thread matrix class with Nthread = "
	   << Nthread << endl;
    #else
      cout << "~spoolesMV() serial matrix class" << endl;
    #endif
  }

  class matrix
  {
  public:
    matrix(FeMdata mdata);
    ~matrix();
    
    void zero();
    void assemble(int, int*, double*);
    void assemble(int, int*, double*, spoolesMV::vector&, double*);
    void assemblepc(int, int*, double*);    

    void times(vector &x, vector &y);   // y = A x 

    void solve(vector&, vector&);    
    void reSolve(vector&, vector&);  

    const int size() const {return (n);}

    //xiaohan
    void setvalue(int,int,double);
  private:
    int n, permuteflag;
    bool pivotFlag;
    InpMtx *A;

  #if (Nthread > 1)
    BridgeMT *bridge;
  #else
    Bridge *bridge;
  #endif
  };

  class vector
  {
  public:
    vector(int);
    ~vector();
    
    vector& zero();
    const int length() { return(n); }
    vector& operator=(const vector&);

    double normInfinity();
    vector& operator+=(vector&);
    vector& operator-=(vector&);
    
    double operator[](int i) { return(entries[i]); }
    //xiaohan
    vector& operator*=(double scal);

    void assemble(int, int*, double*);
    void extract (int, int*, double*) const;
    void setvalues(int, int*, double*);
    void setvalue(int, double);

    const double* buffer() {return(entries);}
  private:
    int n;
    double *entries;
    DenseMtx *v;

    friend void matrix::solve(vector &, vector&);
    friend void matrix::reSolve(vector &, vector&);
    friend void matrix::times(vector &x, vector &y);   // y = A x 
    friend void matrix::assemble(int, int*, double *a, vector&, double*);
  };
};

spoolesMV::vector::vector(int size) 
  : n(size)
{
  v = DenseMtx_new() ;

  DenseMtx_init(v, SPOOLES_REAL, 0, 0, n, 1, 1, n);  // column major

  entries = DenseMtx_entries(v);

  DenseMtx_zero(v);
}

spoolesMV::vector::~vector() 
{
  DenseMtx_free(v);
} 

spoolesMV::vector& spoolesMV::vector::zero()  
{ 
  DenseMtx_zero(v);

  return(*this);
}

void spoolesMV::vector::assemble(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "spoolesMV::vector::assemble(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) entries[index[i]] += fe[i];
}

void spoolesMV::vector::setvalue(int row, double fe)
{
 #ifdef DeBug
  if(row <0 || row >=n)
    {
      cerr << "spoolesMV::vector::setvalue(): Index out of range: index " << std::endl; 
      throw("Index out of range");
    }
  #endif
  entries[row] = fe;
}

void spoolesMV::vector::setvalues(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "spoolesMV::vector::setvalues(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) entries[index[i]] = fe[i];
}

void spoolesMV::vector::extract(int nvar, int index[], double fe[]) const
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "spoolesMV::vector::extract(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(int i = 0; i < nvar; i++) fe[i] = entries[index[i]];
}

spoolesMV::vector& spoolesMV::vector::operator=(const vector &a)
{
  if(this != &a)
    {
    for(int i = 0; i < n; i++) entries[i] = a.entries[i];
    }

  return(*this);
}
//xiaohan
spoolesMV::vector& spoolesMV::vector::operator*=(double scal)
{
  for(int i = 0;i<n;i++) entries[i] *=scal;
  return(*this);
}
spoolesMV::vector& spoolesMV::vector::operator+=(vector &a)
{
  for(int i = 0; i < n; i++) entries[i] += a.entries[i];

  return(*this);
}

spoolesMV::vector& spoolesMV::vector::operator-=(vector &a)
{
  for(int i = 0; i < n; i++) entries[i] -= a.entries[i];

  return(*this);
}

double spoolesMV::vector::normInfinity()
{
  return( DenseMtx_maxabs(v) );
}

spoolesMV::matrix::matrix(FeMdata mdata)
  : n(mdata.nvar), pivotFlag(mdata.pivot)
{
  permuteflag = 1;

  #if (Nthread > 1)
    bridge = BridgeMT_new();
  #else
    bridge = Bridge_new();
  #endif

  A = InpMtx_new();
 
  InpMtx_init(A, INPMTX_BY_ROWS, SPOOLES_REAL, 0, 0);
}

spoolesMV::matrix::~matrix()
{
  #if (Nthread > 1)
    BridgeMT_free(bridge);
  #else
    Bridge_free(bridge);
  #endif

  InpMtx_free(A);
}

void spoolesMV::matrix::assemblepc(int nvar, int index[], double *ae)
{
  return;    // preconditioner assemble for itterative schemes
}

void spoolesMV::matrix::assemble(int nvar, int index[], double *ae)
{

#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "spoolesMV::matrix::assemble(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif

  InpMtx_inputRealMatrix(A, nvar, nvar, nvar, 1, index, index, ae);  
}

void spoolesMV::matrix::assemble(int nvar, int index[], double *ae, 
			  spoolesMV::vector& f, double *fe)
{
  int i;

#ifdef DeBug
  for(i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "spoolesMV::matrix::assemble(): Index out of range: index[" 
	   << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  for(i = 0; i < nvar; i++) f.entries[index[i]] += fe[i];

  InpMtx_inputRealMatrix(A, nvar, nvar, nvar, 1, index, index, ae);  
}

void spoolesMV::matrix::zero()
{
  InpMtx_clearData(A);
}

void spoolesMV::matrix::times(spoolesMV::vector &x, spoolesMV::vector &y) 
{
  double a0 = 0.0, a1 = 1.0;

  InpMtx_nonsym_gmmm(A, &a0, y.v, &a1, x.v);

  return;
}

#if (Nthread > 1)

void spoolesMV::matrix::solve(spoolesMV::vector& u, spoolesMV::vector& f) 
{ 
  int err;

  if(pivotFlag) bridge->pivotingflag = SPOOLES_PIVOTING;

  if( BridgeMT_setMatrixParams(bridge, n, SPOOLES_REAL, SPOOLES_NONSYMMETRIC)
      != 1) goto error;

  if( BridgeMT_setup(bridge, A) != 1 ) goto error;

  if( BridgeMT_factorSetup(bridge, Nthread, 0, 0.0) != 1 ) goto error;

  if( BridgeMT_factor(bridge, A, permuteflag, &err) != 1 ) goto error;

  if( BridgeMT_solveSetup(bridge) != 1 ) goto error;

  if( BridgeMT_solve(bridge, permuteflag, u.v, f.v) != 1 ) goto error;

  return;

 error:
  cerr << "spoolesMV::matrix::solve(): Error -- Execution Terminates" << endl;
  throw("spoolesMV::matrix::solve() error");
}

void spoolesMV::matrix::reSolve(spoolesMV::vector& u, spoolesMV::vector& f) 
{ 
  if( BridgeMT_solve(bridge, permuteflag, u.v, f.v) != 1 ) goto error;

  return;

 error:
  cerr << "spoolesMV::matrix::solve(): Error -- Execution Terminates" << endl;
  throw("spoolesMV::matrix::solve() error");
}

#else   // serial solve

void spoolesMV::matrix::solve(spoolesMV::vector& u, spoolesMV::vector& f) 
{ 
  int err;

  if(pivotFlag) bridge->pivotingflag = SPOOLES_PIVOTING;

  if( Bridge_setMatrixParams(bridge, n, SPOOLES_REAL, SPOOLES_NONSYMMETRIC)
      != 1) goto error;

  if( Bridge_setup(bridge, A) != 1 ) goto error;

  if( Bridge_factor(bridge, A, permuteflag, &err) != 1 ) goto error;

  if( Bridge_solve(bridge, permuteflag, u.v, f.v) != 1 ) goto error;

  return;

 error:
  cerr << "spoolesMV::matrix::solve(): Error -- Execution Terminates" << endl;
  throw("spoolesMV::matrix::solve() error");
}

void spoolesMV::matrix::reSolve(spoolesMV::vector& u, spoolesMV::vector& f) 
{ 
  if( Bridge_solve(bridge, permuteflag, u.v, f.v) != 1 ) goto error;

  return;

 error:
  cerr << "spoolesMV::matrix::solve(): Error -- Execution Terminates" << endl;
  throw("spoolesMV::matrix::solve() error");
}
//xiaohan
void spoolesMV::matrix::setvalue(int row, int col, double value)
{
  
}
#endif
