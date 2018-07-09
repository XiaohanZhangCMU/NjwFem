#include<iostream>
#include<cstdlib>
#include<cmath>

#include "petscksp.h"

class petscMV
{
public:
  class vector;
  class matrix;

  petscMV(int &argc, char ***argv);
  ~petscMV();

  class matrix
  {
  public:
    matrix(FeMdata mdata);
    ~matrix();

    void zero();
    const int size() const {return(n);}
    
    void assemble(int, int*, double*);
    void assemble(int, int*, double*, vector&, double*);
    void assemblepc(int, int*, double*);    

    void times(vector &x, vector &y);   // y = A x 

    void solve(vector&, vector&);
    void reSolve(vector&, vector&);

  private:
    double rtol;
    bool pcflag;
    int ierr, n;
    Mat a, apc;
    KSP ksp;

    Vec dd, ii;   // vectors used to scale the matrix
  };

  class vector	
  {
  public:	
    vector(int) ;
    ~vector();

    void zero();

    double* buffer();
    const int length() const {return(n);}

    double dot(vector &);
    double norm() {return(sqrt(dot(*this)));}

    // vector operator+(const vector&) const;
    // vector operator-() const;
    // vector operator-(const vector&) const;

    vector& operator+=(vector&);
    vector& operator-=(vector&);
    double normInfinity();
    
    double& operator[](int);        // optional, extract() gets subvector
    void print();

    vector& operator=(vector &);
    
    void assemble(int, int *, double *); 
    void extract (int, int *, double *); 
    void setvalues (int, int *, double *); 

  private:
    bool restored, assembled;
    int ierr, n;
    PetscScalar *varray;
    Vec v;

    friend void matrix::solve(vector&, vector&);
    friend void matrix::reSolve(vector&, vector&);
    friend void matrix::times(vector &x, vector &y);   // y = A x 
    friend void matrix::assemble(int, int*, double *a, vector&, double*);
  };

};

petscMV::petscMV(int &argc, char ***argv)
{
  PetscInitialize (&argc, argv, PETSC_NULL, PETSC_NULL) ;
}

petscMV::~petscMV()
{
  PetscFinalize();

  cout << "petscMV used" << endl;
  cout << "Tolerences: rtol from matrix constructor"
       << ", atol = default, PC = default" << endl;
}

petscMV::vector::vector(int nvar): n(nvar)
{ 
  restored  = true;
  assembled = true;

  ierr = VecCreateMPI(PETSC_COMM_SELF, PETSC_DECIDE, nvar, &v);
}

petscMV::vector::~vector() 
{
  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

  ierr = VecDestroy(&v);
}

void petscMV::vector::zero()  
{ 
  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

  ierr = VecZeroEntries(v);

  return;
}

double* petscMV::vector::buffer()  
{ 
  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

  if(restored)
    {  
      ierr = VecGetArray(v, &varray);
      restored = false;
    }

  return(varray);
}

double petscMV::vector::dot(petscMV::vector &w)
{ 
  PetscScalar vdotw;

  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  if(!w.restored)  
    {
      ierr = VecRestoreArray(w.v, &(w.varray));
      w.restored = true;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

  if(!(w.assembled))  
    {
      ierr = VecAssemblyBegin(w.v) ;
      ierr = VecAssemblyEnd(w.v) ;
      w.assembled = true;
    }

  VecDot(this->v, w.v, &vdotw);

  return(vdotw);
}

void petscMV::vector::assemble(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "Index out of range: index[" << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif

  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  ierr =  VecSetValues(v, nvar, index, fe, ADD_VALUES) ; 

  assembled = false;
}

void petscMV::vector::setvalues(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "Index out of range: index[" << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif

  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  ierr =  VecSetValues(v, nvar, index, fe, INSERT_VALUES) ; 

  assembled = false;
}

void petscMV::vector::extract(int nvar, int index[], double fe[])
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "Index out of range: index[" << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif

  if(restored)
    {  
      ierr = VecGetArray(v, &varray);
      restored = false;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

  for(int i = 0; i < nvar; i++) fe[i] = varray[index[i]];
}

petscMV::vector& petscMV::vector::operator=(vector &a)
{
  if(this != &a)
    {
      if(!restored)  
	{
	  ierr = VecRestoreArray(v, &varray);
	  restored = true;
	}
      
      if(!(a.restored)) 
	{
	  ierr = VecRestoreArray(a.v, &(a.varray));
	  a.restored = true;
	}
      
      if(!assembled)  
	{
	  ierr = VecAssemblyBegin(v) ;
	  ierr = VecAssemblyEnd(v) ;
	  assembled = true;
	}

      if(!(a.assembled))  
	{
	  ierr = VecAssemblyBegin(a.v) ;
	  ierr = VecAssemblyEnd(a.v) ;
	  a.assembled = true;
	}
      
      ierr = VecCopy(a.v, v);
    }
  
  return(*this);
}

double& petscMV::vector::operator[](int i)
{
  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }
  
  if(restored)
    {  
      ierr = VecGetArray(v, &varray);
      restored = false;
    }

  return(varray[i]);
}

petscMV::vector& petscMV::vector::operator+=(vector &a)
{
  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }
  
  VecAXPY(v, 1.0, a.v);

  return(*this);
}

petscMV::vector& petscMV::vector::operator-=(vector &a)
{
  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

  VecAXPY(v, -1.0, a.v);

  return(*this);
}

double petscMV::vector::normInfinity()
{
  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

#ifdef PETSC23
  return( VecNorm(v, NORM_INFINITY) );
#endif


#ifdef PETSC30
  PetscReal value;

  VecNorm(v, NORM_INFINITY, &value);
  return value;
#endif

#ifdef PETSC40
  PetscReal value;

  VecNorm(v, NORM_INFINITY, &value);
  return value;
#endif
}

void petscMV::vector::print()
{
  if(!restored)  
    {
      ierr = VecRestoreArray(v, &varray);
      restored = true;
    }

  if(!assembled)  
    {
      ierr = VecAssemblyBegin(v) ;
      ierr = VecAssemblyEnd(v) ;
      assembled = true;
    }

  VecView(v, PETSC_VIEWER_STDOUT_SELF);

  return;
}

petscMV::matrix::matrix(FeMdata mdata)
  :  rtol(mdata.tol), pcflag(mdata.pivot), n(mdata.nvar)
{
  ierr = MatCreateAIJ (PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
		   n, n, mdata.nzero, PETSC_NULL, 0, PETSC_NULL, &a);

 ierr = KSPCreate(PETSC_COMM_SELF, &ksp) ;
    
 VecCreateMPI(PETSC_COMM_SELF, PETSC_DECIDE, n, &dd);
 VecCreateMPI(PETSC_COMM_SELF, PETSC_DECIDE, n, &ii);
  
 VecSet(ii, 1.0);   // vectors used to scale the matrix (and pc)

 if(pcflag)
   {
     ierr = MatCreateAIJ (PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
		      n, n, mdata.nzero, PETSC_NULL, 0, PETSC_NULL, &apc);

     ierr = KSPSetOperators(ksp, a, apc, DIFFERENT_NONZERO_PATTERN);
   }
 else
   {
     ierr = KSPSetOperators(ksp, a, a, DIFFERENT_NONZERO_PATTERN);
   }

 PC pc;  
 ierr = KSPGetPC(ksp, &pc);

 // ierr = PCSetType(pc, PCNONE);
 // ierr = PCSetType(pc, PCJACOBI); 
 // ierr = PCSetType(pc, PCILU);

 // KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,
 //                          PetscReal dtol,PetscInt maxits) 
 // KSPDefaultConverged() reaches convergence when
 //     rnorm < MAX (rtol * rnorm_0, abstol);  
 // rnom = 2-norm residual value
 // rnorm_0 is the two norm of the right hand side.

 ierr = KSPSetTolerances(ksp, rtol, PETSC_DEFAULT,
  			 PETSC_DEFAULT, PETSC_DEFAULT);
}

petscMV::matrix::~matrix() 
{
  ierr = MatDestroy(&a);
  ierr = KSPDestroy(&ksp);

  if(pcflag) ierr = MatDestroy(&apc);
}

void petscMV::matrix::assemble(int nvar, int index[], double *ae)
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "Index out of range: index[" << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  MatSetValues(a, nvar, index, nvar, index, ae, ADD_VALUES);

  if(pcflag)
    MatSetValues(apc, nvar, index, nvar, index, ae, ADD_VALUES);
}

void petscMV::matrix::assemblepc(int nvar, int index[], double *ae)
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
    if((index[i] < 0) || (index[i] >= n))
      {
      cerr << "Index out of range: index[" << i << "] = " 
	   << index[i] << ", n = " << n << endl;

      throw("Index out of range");
      }
    }
#endif
  
  if(!pcflag) return;   // ignoring this allows same code to work in both modes

  MatSetValues(apc, nvar, index, nvar, index, ae, ADD_VALUES);
}

void petscMV::matrix::assemble(int nvar, int index[], double *ae, 
			  petscMV::vector& f, double *fe)
{
#ifdef DeBug
  for(int i = 0; i < nvar; i++)
    {
      if((index[i] < 0) || (index[i] >= n))
	{
	  cerr << "Index out of range: index[" << i << "] = " 
	       << index[i] << ", n = " << n << endl;
	  
	  throw("Index out of range");
	}
    }
#endif
      
  if(!(f.restored))  
    {
      ierr = VecRestoreArray(f.v, &(f.varray));
      f.restored = true;
    }

  ierr = VecSetValues(f.v, nvar, index, fe, ADD_VALUES) ; 
  
  ierr = MatSetValues(a, nvar, index, nvar, index, ae, ADD_VALUES);
  
  if(pcflag)
    MatSetValues(apc, nvar, index, nvar, index, ae, ADD_VALUES);

  return;
}
  
void petscMV::matrix::zero()
{
  ierr = MatZeroEntries(a);

  if(pcflag)
    {
      ierr = MatZeroEntries(apc);
      
      ierr = KSPSetOperators(ksp, a, apc, DIFFERENT_NONZERO_PATTERN);
    }
  else
    {
      ierr = KSPSetOperators(ksp, a, a, DIFFERENT_NONZERO_PATTERN);
    }

  return;
}


void petscMV::matrix::times(petscMV::vector &x, petscMV::vector &y) 
{
  if(!(x.restored))  
    {
      ierr = VecRestoreArray(x.v, &(x.varray));
      x.restored = true;
    }

  if(!(y.restored))  
    {
      ierr = VecRestoreArray(y.v, &(y.varray));
      y.restored = true;
    }

  ierr = MatMult(a, x.v, y.v);

  return;
}

void petscMV::matrix::solve(petscMV::vector &u, petscMV::vector &f) 
{ 
  if(!(f.restored))  
    {
      ierr = VecRestoreArray(f.v, &(f.varray));
      f.restored = true;
    }

  if(!(u.restored))  
    {
      ierr = VecRestoreArray(u.v, &(u.varray));
      u.restored = true;
    }

  if(pcflag) ierr = MatAssemblyBegin(apc, MAT_FINAL_ASSEMBLY) ;
  ierr = MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY) ;
  ierr = VecAssemblyBegin(f.v) ;
  ierr = VecAssemblyBegin(u.v) ;
  ierr = VecAssemblyEnd(f.v) ;
  ierr = VecAssemblyEnd(u.v) ;
  ierr = MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY) ;
  if(pcflag) ierr = MatAssemblyEnd(apc, MAT_FINAL_ASSEMBLY) ;

  f.assembled = true;
  u.assembled = true;

  // scale the rows of the matrix

  if(pcflag) MatGetRowMaxAbs(apc, dd, PETSC_NULL);
  else       MatGetRowMaxAbs(a  , dd, PETSC_NULL);

  // if(pcflag) MatGetRowMax(apc, dd); // petsc 2.3.0 gives max abs
  // else       MatGetRowMax(a  , dd); 
  
  VecReciprocal(dd);
  
  VecPointwiseMult(f.v, f.v, dd);   // scale f
  
  MatDiagonalScale(a,   dd , ii);  
  if(pcflag) MatDiagonalScale(apc, dd , ii);

  // solve the system

  ierr = KSPSolve(ksp, f.v, u.v);
}

void petscMV::matrix::reSolve(petscMV::vector &u, petscMV::vector &f) 
{ 
  if(!(f.restored))  
    {
      ierr = VecRestoreArray(f.v, &(f.varray));
      f.restored = true;
    }

  if(!(u.restored))  
    {
      ierr = VecRestoreArray(u.v, &(u.varray));
      u.restored = true;
    }

  ierr = VecAssemblyBegin(f.v) ;
  ierr = VecAssemblyBegin(u.v) ;
  ierr = VecAssemblyEnd(f.v) ;
  ierr = VecAssemblyEnd(u.v) ;

  f.assembled = true;
  u.assembled = true;

  // scale the rows of the matrix
  
  VecPointwiseMult(f.v, f.v, dd);   // scale f
  
  // solve the system

  ierr = KSPSolve(ksp, f.v, u.v);

  return;
}

