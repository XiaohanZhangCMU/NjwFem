#if 0
Need to be done: MatSetValuesBlocked, instead ee.dof() but ee.node()
Pass in negative index,automatically eliminate dirichlet b.c
Make preallocation of matrix and vectors specific. increase perfomance by a factor of 50.
#endif

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
    matrix(FeMdata mdata,int n, int m, int mype, int npe, MPI_Comm comm);
    ~matrix();

    void zero();
    const int size() const {return(n);}
    
    void assemble(int, int*, double*);
    void assemble(int, int*, double*, vector&, double*);
    void assemblepc(int, int*, double*);    

    void times(vector &x, vector &y);   // y = A x 

    void solve(vector&, vector&);
    void reSolve(vector&, vector&);       
    void getInfo(void);
    void print();
  private:
    MPI_Comm    comm;
    double      rtol;
    bool        pcflag;
    PetscInt    ierr, n, m, mype, npe;
    Mat         a, apc;
    KSP         ksp;
    PC          pc;
    PCType      type;
    Vec         dd, ii;   // vectors used to scale the matrix    
    MatInfo     info;
    bool        assembled;
  };

  class vector	
  {
  public:	
    vector(int n, int m, int mype, int npe, MPI_Comm comm) ;
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
    void getInfo(void);
  private:
    MPI_Comm comm;
    bool restored, assembled;
    int ierr, n, m, mype, npe;
;
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

  cout << "MPI petscMV used (Krylov-CG)" << endl;
  cout << "Tolerences: rtol from matrix constructor"
       << ", atol = default, PC = default" << endl;
}

petscMV::vector::vector(int n, int m, int mype, int npe, MPI_Comm comm)
{ 
  int Ndim = NDIM;

  this->n = n;
  this->m = m;
  this->mype = mype;
  this->npe = npe;
  this->comm = comm;

  restored  = true;
  assembled = true;
  
  ierr = VecCreate(PETSC_COMM_WORLD,&v);
  ierr = VecSetSizes(v,m,n);
  ierr = VecSetBlockSize(v,Ndim);
  ierr = VecSetFromOptions(v);
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
#if 0
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

  VecScatter     ctx;
  Vec            v_all;

  VecScatterCreateToAll(v,&ctx,&v_all);

  VecScatterBegin(ctx,v,v_all,INSERT_VALUES,SCATTER_FORWARD);

  VecScatterEnd(ctx,v,v_all,INSERT_VALUES,SCATTER_FORWARD);

  ierr = VecGetArray(v_all, &varray);
  // if(!assembled)  
  //   {
  //     ierr = VecAssemblyBegin(v) ;
  //     ierr = VecAssemblyEnd(v) ;
  //     assembled = true;
  //   }

  for(int i = 0; i < nvar; i++) fe[i] = varray[index[i]];
  // for(int i = 0; i < nvar; i++) std::cout<<"; mype = "<<mype<<"; index["<<i<<"]="<<index[i];
  // std::cout<<std::endl;

  VecRestoreArray(v_all,&varray);

  VecScatterDestroy(&ctx); 
  VecDestroy(&v_all); 
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
  
  ierr = VecAssemblyBegin(v) ;
  ierr = VecAssemblyEnd(v) ;
  assembled = true;
  
  VecView(v, PETSC_VIEWER_STDOUT_WORLD);

  return;
}

petscMV::matrix::matrix(FeMdata mdata, int n, int m, int mype, int npe, MPI_Comm comm)
  :  rtol(mdata.tol), pcflag(mdata.pivot)
{
  int Ndim = NDIM;
  this->n = n;
  this->m = m;
  this->mype = mype;
  this->npe = npe;
  this->comm = comm;

  ierr = MatCreate(PETSC_COMM_WORLD,&a); 
  ierr = MatSetSizes(a,m,m,n,n); 
  ierr = MatSetBlockSize(a,Ndim); 
  ierr = MatSetType(a,MATAIJ); 
  ierr = MatSeqAIJSetPreallocation(a,PETSC_DEFAULT,NULL); 
  ierr = MatMPIAIJSetPreallocation(a,PETSC_DEFAULT,NULL,PETSC_DEFAULT,NULL); 

  //KSP
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
  ierr = KSPSetType(ksp, KSPCG);
  ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE);
  ierr = KSPGetPC(ksp, &pc);
  ierr = PCSetType(pc, PCGAMG); /* default */
  ierr = KSPSetFromOptions(ksp);
  ierr = PCGetType(pc, &type);
  ierr = KSPSetOperators(ksp, a, a, SAME_NONZERO_PATTERN); 
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

  ierr = MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY) ;
  ierr = VecAssemblyBegin(f.v) ;
  ierr = VecAssemblyBegin(u.v) ;
  ierr = VecAssemblyEnd(f.v) ;
  ierr = VecAssemblyEnd(u.v) ;
  ierr = MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY) ;
  
  f.assembled = true;
  u.assembled = true;

  // solve the system
  ierr = KSPSetUp(ksp); 
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
  ierr = KSPSetUp(ksp); 
  ierr = KSPSolve(ksp, f.v, u.v);

  return;
}

void petscMV::matrix:: getInfo(void)
{ 
  MatGetInfo(a,MAT_LOCAL,&info);
  PetscInt localRow, localCol, globalRow, globalCol;
  MatGetSize(a,&globalRow,&globalCol);
  MatGetLocalSize(a,&localRow,&localCol);
  PetscPrintf(comm,"------------------------------------\n");
  PetscSynchronizedPrintf(comm,"mype = %d, npe = %d\n",mype,npe);
  PetscSynchronizedPrintf(comm,"info.mallocs = %d\n", info.mallocs);
  PetscSynchronizedPrintf(comm,"info.nz_allocated = %d\n",info.nz_allocated);
  PetscSynchronizedPrintf(comm,"block_size = %d\n",info.block_size);
  PetscSynchronizedPrintf(comm,"local  row,col = %d, %d \n",localRow, localCol);
  PetscSynchronizedPrintf(comm,"global row,col = %d, %d \n",globalRow, globalCol);
  PetscSynchronizedPrintf(comm,"------------------------------------\n");
  PetscSynchronizedFlush(comm);
}



void petscMV::vector:: getInfo(void)
{ 
  PetscInt localsz, globalsz;
  VecGetLocalSize(v,&localsz);
  VecGetSize(v,&globalsz);
  PetscPrintf(comm,"------------------------------------\n");
  PetscSynchronizedPrintf(comm,"mype = %d, npe = %d\n",mype,npe);
  PetscSynchronizedPrintf(comm,"local = %d, global = %d\n",localsz,globalsz);  
  PetscSynchronizedPrintf(comm,"------------------------------------\n");
  PetscSynchronizedFlush(comm);
}

void petscMV::matrix:: print()
{
  ierr = MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY) ;
  ierr = MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY) ;
  
  MatView(a,PETSC_VIEWER_STDOUT_WORLD);
}

