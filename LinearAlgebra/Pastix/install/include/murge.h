#ifndef PREC_DOUBLE
#define PREC_DOUBLE
#endif
#ifndef INTSIZE32
#define INTSIZE32
#endif

/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
**
** This file is part of the PaStiX parallel sparse matrix package.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/* 
   Title: Murge
     Murge is close to the Merge.
  
   Description:

     Murge is a common interface definition to multiple solver.

     It has been initiated by *HIPS* and *PaStiX* solvers developpers
     in january 2009.
 
  Integers:

     Depending of your compilation options, *INTS* and *INTL* can be
     32 or 64 bits.
     
     In C user will be abble to use *INTS* types.  In Fortran user
     will have to use the correct type, corresponding to the
     compilation options.
   
     If user doesn't define any compilation option, *INTS* and *INTL*
     will both be int.
   
     If user defines *-DINTSIZE32*, *INTS* and *INTL* will both be 32
     bits integers.
   
     If user defines *-DINTSIZE64*, *INTS* will be 32 bits integers
     whereas *INTL* will be INT64 integers.
     
     If user defines *-DINTSSIZE64*, both *INTS* and *INTL* will be 64
     bit integers.
     
     We recommand not to use *-DINTSSIZE64*, as *INTS* integers are
     supposed to keep being small integers, and there is no need to
     make them longs.

  Coefficients:

     The coefficients of the matrices and vectors is of type *COEF*.
   
     *COEF* can be simple or double precision and real or complex
     depending on compilation options.
   
     If user defines nothing, *COEF* will be real simple precision.
     
     If user defines *-DPREC_DOUBLE*, *COEF* will be double
     precision.
     
     Defining nothing or *-DPREC_SIMPLE* will result in setting *COEF*
     precision to simple.
   
     If user defines *-DTYPE_COMPLEX*, *COEF* will be complex, otherwise,
     if user defines nothing or *-DTYPE_REAL* it will be REAL.

  Floating parameters:
    
     Floating parameters are of type *REAL* which can be simple or 
     double precision.
     
     As for Coefficients, *REAL* precision also depends on 
     *-DPREC_DOUBLE* compilation option.
     
     Defining nothing or *-DPREC_SIMPLE* will result in setting *REAL*
     precision to simple.

  Nodes and Unknowns:

     In murge exists node and unknown, a node correspond to a node 
     in the physical problem, which can be composed of 
     multiple degrees of freedom.
     
     Thus, the number of unknowns correspond to the number of node 
     times the degree of freedom.
     

  Authors:

     HIPS and PaStiX developpers teams :

       Mathieu Faverge   - faverge@labri.fr 
       Jérémie Gaidamour - gaidamou@labri.fr 
       Pascal  Hénon     - henon@labri.fr 
       Xavier  Lacoste   - lacoste@labri.fr 
       Pierre  Ramet     - ramet@labri.fr     
       
*/
#ifndef MURGE_H
#define MURGE_H

#define MURGE_INTERFACE_VERSION 0

#if !(defined INTL) || !(defined INTS)

#ifdef INTSIZE32
#define INTS           int32_t
#define INTL           int32_t
#elif (defined INTSIZE64)
#define INTS           int32_t
#define INTL           int64_t
#elif (defined INTSSIZE64)
#define INTS           int64_t
#define INTL           int64_t
#elif (defined LONG)
#define INTS           int
#define INTL           long
#else
#define INTS           int
#define INTL           int
#endif

#endif


#ifdef    PREC_DOUBLE
#define REAL double 
#ifdef    TYPE_COMPLEX
#define COEF double complex
#else  /* TYPE_COMPLEX */
#define COEF double
#endif /* TYPE_COMPLEX */
#else  /* PREC_DOUBLE  */
#define REAL float
#ifdef    TYPE_COMPLEX
#define COEF float complex
#else  /* TYPE_COMPLEX */
#define COEF float
#endif /* TYPE_COMPLEX */
#endif /* PREC_DOUBLE  */

#ifdef MPI_VERSION
#ifdef    PREC_DOUBLE
#define MURGE_MPI_REAL MPI_DOUBLE
#ifdef    TYPE_COMPLEX
#define MURGE_MPI_COEF MPI_DOUBLE_COMPLEX
#else  /* TYPE_COMPLEX */
#define MURGE_MPI_COEF MPI_DOUBLE
#endif /* TYPE_COMPLEX */
#else  /* PREC_DOUBLE  */
#define MURGE_MPI_REAL MPI_FLOAT
#ifdef    TYPE_COMPLEX
#define MURGE_MPI_COEF MPI_COMPLEX
#else  /* TYPE_COMPLEX */
#define MURGE_MPI_COEF MPI_FLOAT
#endif /* TYPE_COMPLEX */
#endif /* PREC_DOUBLE  */
#endif /* MPI_VERSION */
/* 
   Group: Solver setup functions
*/
/*
  Function: MURGE_Initialize

  Allocate the instance arrays which will keeps intern data for all
  solver instances.

  If user is creating several threads calling the solver, this function 
  has to be called before creating threads to insure that the execution
  is thread safe.

  Parameters: 
    idnbr - Maximum number of solver instances that will be
            launched.

  Returns: 
    MURGE_SUCCESS - If function runned successfully.
    MURGE_ERR_ALLOCATE - If for some reason, allocation was not
                     successfull.


  Fortran interface:
  >
  > SUBROUTINE MURGE_INITIALIZE(IDNBR, IERROR)
  >   INTS, INTENT(IN)  :: IDNBR
  >   INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_INITIALIZE
*/
INTS MURGE_Initialize(INTS idnbr);

/*
  Function: MURGE_SetDefaultOptions

  Sets default options, for solver instance number *id*.

  The default option set correspond to *stratnum* strategy ID, 
  depending on the solver.

  Needs <MURGE_Initialize> to be called before 
  to allocate solver instances array.
  
  Parameters: 
    id       - Solver instance identification number.
    stratnum - Strategy for the default option Set.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_Initialize> was not called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range 
                          or *stratnum* is not valid.
    MURGE_ERR_ALLOCATE  - If couldn't create solver instance.


  Fortran interface:
  >
  > SUBROUTINE MURGE_SETDEFAULTOPTIONS(ID, STRATNUM, IERROR)
  >   INTS, INTENT(IN)  :: ID, STRATNUM
  >   INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_SETDEFAULTOPTIONS
 */
INTS MURGE_SetDefaultOptions(INTS id, INTS stratnum);

/*
  Function: MURGE_SetOptionINT

  Sets integer option, indicated by *number*, to *value* for the
  solver instance number *id*.

  Needs <MURGE_SetDefaultOptions> to be called before to initiate
  solver instance data.

  Parameters: 
    id     - Solver instance identification number.  
    number - Identification of the integer parameter.  
    value  - value to set the parameter to.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
                          *number* or *value* are not valid. 

  Fortran interface:
  >
  > SUBROUTINE MURGE_SETOPTIONINT(ID, NUMBER, VALUE, IERROR)
  >   INTS, INTENT(IN)    :: ID, NUMBER, VALUE
  >   INTS, INTENT(OUT)   :: IERROR
  > END SUBROUTINE MURGE_SETOPTIONINT
*/
INTS MURGE_SetOptionINT(INTS id, INTS number, INTS value);

/*
  Function: MURGE_SetOptionREAL

  Sets real option, indicated by *number*, to *value* for the
  solver instance number *id*.

  Needs <MURGE_SetDefaultOptions> to be called before to initiate
  solver instance data.

  Parameters: 
    id     - Solver instance identification number.  
    number - Identification of the integer parameter.  
    value  - value to set the parameter to.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
                          *number* or *value* are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_SETOPTIONREAL(ID, NUMBER, VALUE, IERROR)
  >   INTS, INTENT(IN)    :: ID, NUMBER
  >   REAL, INTENT(IN)    :: VALUE
  >   INTS, INTENT(OUT)   :: IERROR
  > END SUBROUTINE MURGE_SETOPTIONREAL
*/
INTS MURGE_SetOptionREAL(INTS id, INTS number, REAL value);

#ifdef MPI_VERSION 
/*
  Function: MURGE_SetCommunicator

  Sets MPI communicator for the given solver instance.

  Needs <MURGE_SetDefaultOptions> to be called before to initiate
  solver instance data.

  Musn't be called before <MURGE_SAVE>, <MURGE_LOAD>, 
  <MURGE_GetLocalNodeNbr> nor <MURGE_GetLocalUnknownNbr> 
  because the solver as to be runned with the same MPI 
  communicator all along.

  If this function is not called, MPI communicator will be
  *MPI_COMM_WORLD*.

  This function may not exist if the solver 
  has been compiled without MPI.
  
  Parameters: 
    id      - Solver instance identification number.  
    mpicomm - MPI communicator to affect the solver to.  

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before or if it is called after 
		          the solver starts its computing tasks.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_SETCOMMUNICATOR(ID, MPICOMM, IERROR)
  >   INTS,            INTENT(IN)  :: ID
  >   INTEGER,         INTENT(IN)  :: MPICOMM
  >   INTS,            INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_SETCOMMUNICATOR
*/
INTS MURGE_SetCommunicator(INTS id, MPI_Comm mpicom);
#endif 


/* 
   Group: Graph setup function
*/
/*
  Function: MURGE_GraphBegin

  Begin building the adjency graph for renumbering and all 
  preprocessing.

  Needs <MURGE_SetDefaultOptions> to be called before to initiate
  solver instance data.
  
  Allocate temporary structures needed to build the graph.
  
  Parameters: 
    id      - Solver instance identification number.  
    N       - Global number of nodes in the graph.
    edgenbr - Number of edges which will be added in the graph by proc.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before.
    MURGE_ERR_ALLOCATE  - If Allocate didn't worked.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
                          *N* or *edgenbr* are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GRAPHBEGIN(ID, N, EDGENBR, IERROR)
  >   INTS,      INTENT(IN)  :: ID, N
  >   INTL,      INTENT(IN)  :: EDGENBR
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GRAPHBEGIN
*/ 
INTS MURGE_GraphBegin(INTS id, INTS N, INTL edgenbr);


/*
  Function: MURGE_GraphEdge

  Adds an edge to the graph user is currently building.

  Needs <MURGE_GraphBegin> to be called before to initiate
  building sequence.

  This function depends on integer parameter *MURGE_BASEVAL*.
  
  Parameters: 
    id      - Solver instance identification number.  
    ROW     - First vertex of the edge.
    COL     - Second vertex of the edge.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_GraphBegin> was not 
                          called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
                          *I* or *J* are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GRAPHEDGE(ID, ROW, COL, IERROR)
  >   INTS,      INTENT(IN)  :: ID, ROW, COL
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GRAPHEDGE
*/ 
INTS MURGE_GraphEdge(INTS id, INTS COL, INTS ROW);

/*
  Function: MURGE_GraphEnd

  End the graph building.
  
  Needs <MURGE_GraphBegin> to be called before.

  Parameters: 
    id      - Solver instance identification number.  

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_GraphBegin> was not 
                          called before or there are missing edges.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GRAPHEND(ID, IERROR)
  >   INTS,      INTENT(IN)  :: ID
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GRAPHEND
*/ 
INTS MURGE_GraphEnd(INTS id);

/*
  Function: MURGE_GraphGlobalCSR

  Build an adjency graph from a Compress Sparse Row matrix pattern.
  
  Needs <MURGE_SetDefaultOptions> to be called before.

  This function depends on integer parameter *MURGE_BASEVAL*.

  Parameters: 
    id      - Solver instance identification number.
    N       - Global number of columns
    rowptr  - Index of the first element of each row in *COLS* array.
    COLS    - Global column numbers array.
    root    - Root processor : this processor enter the global data 
              (-1 for all processors).

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or 
                          CSR is not correct.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GRAPHGLOBALCSR(ID, N, ROWPTR, COLS, ROOT, IERROR)
  >   INTS, DIMENSION(0), INTENT(IN)  :: ID, N, ROOT
  >   INTL, DIMENSION(0), INTENT(IN)  :: ROWPTR
  >   INTS, DIMENSION(0), INTENT(IN)  :: COLS
  >   INTS, DIMENSION(0), INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GRAPHGLOBALCSR
*/ 
INTS MURGE_GraphGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS, INTS root);

/*
  Function: MURGE_GraphGlobalCSC

  Build an adjency graph from a Compress Sparse Column matrix pattern.
  
  Needs <MURGE_SetDefaultOptions> to be called before.

  This function depends on integer parameter *MURGE_BASEVAL*.

  Parameters: 
    id      - Solver instance identification number.
    N       - Global number of columns
    COLPTR  - Index of the first element of each column in *ROWS* array.
    ROWS    - Global row number array.
    root    - Root processor : this processor enter the global data 
              (-1 for all processors).  


  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or 
                          CSC is not correct.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GRAPHGLOBALCSC(ID, N, COLPTR, ROWS, ROOT, IERROR)
  >   INTS,               INTENT(IN)  :: ID, N, ROOT
  >   INTL, DIMENSION(0), INTENT(IN)  :: COLPTR
  >   INTL, DIMENSION(0), INTENT(IN)  :: ROWS
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GRAPHGLOBALCSC
*/ 
INTS MURGE_GraphGlobalCSC(INTS id, INTS N, INTL *colptr, INTS *ROWS, INTS root);

/*
  Function: MURGE_GraphGlobalIJV

  Build an adjency graph from a Compress Sparse Column matrix pattern.
  
  Needs <MURGE_SetDefaultOptions> to be called before.

  This function depends on integer parameter *MURGE_BASEVAL*.

  Parameters: 
    id      - Solver instance identification number.
    N       - Global number of unknowns.
    NNZ     - Global number of non zeros.
    ROW     - Global column number array.
edges.
    COL     - Global row number array.
    root    - Root processor : this processor enter the global data 
              (-1 for all processors).  


  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or 
                          graph IJ is not correct.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GRAPHGLOBALIJV(ID, N, NNZ, ROW, COL, ROOT, IERROR)
  >   INTS,               INTENT(IN)  :: ID, N, ROOT
  >   INTL,               INTENT(IN)  :: NNZ
  >   INTS, DIMENSION(0), INTENT(IN)  :: ROW
  >   INTS, DIMENSION(0), INTENT(IN)  :: COL
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GRAPHGLOBALIJV
*/ 
INTS MURGE_GraphGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROW, INTS *COL, INTS root);

/* 
   Group: IO functions

   Allows to save and load solver state after preprocessing.
*/
/*
  Function: MURGE_Save
  
  Runs preprocessing step, if not done yet, and save the result to disk,
  into *directory*, so that it can be resume using <MURGE_Load>.

  Needs the graph to be built.

  Parameters: 
    id        - Solver instance identification number.  
    directory - Path to the directory where to save the solver step.

  In Fortran, *STR_LEN* is the length of the string directory.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If graph hasn't been built.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_IO        - If file(s) couldn't be writen.

  Fortran interface:
  >
  > SUBROUTINE MURGE_SAVE(ID, DIRECTORY, STR_LEN, IERROR)
  >   INTS,             INTENT(IN)  :: ID, STR_LEN
  >   CHARACTER(len=*), INTENT(IN)  :: DIRECTORY
  >   INTS,             INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_SAVE
*/
INTS MURGE_Save(INTS id, char* directory);

/*
  Function: MURGE_Load
  
  Loads preprocessing result from disk, into *directory*, 
  where it had been saved by <MURGE_Save>.

  If preprocessing data was already computed or loaded, it will 
  be overwriten.
  
  Needs <MURGE_SetDefaultOptions> to be called before to initiate
  solver instance data.
  
  Parameters: 
    id        - Solver instance identification number.  
    directory - Path to the directory where to load the solver 
                preprocessing data.

  In Fortran, *STR_LEN* is the length of the string directory.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
                          called before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_IO        - If file(s) couldn't be read.

  Fortran interface:
  >
  > SUBROUTINE MURGE_LOAD(ID, DIRECTORY, STR_LEN, IERROR)
  >   INTS,             INTENT(IN)  :: ID, STR_LEN
  >   CHARACTER(len=*), INTENT(IN)  :: DIRECTORY
  >   INTS,             INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_LOAD
*/
INTS MURGE_Load(INTS id, char* directory);


/* 
   Group: Getting new distribution
   
*/
/*
  Function: MURGE_GetLocalNodeNbr
  
  Computes preprocessing step, if not done, and the number of
  Nodes in the new ditribution of the matrix.

  Needs the graph to be built.

  Parameters:
    id        - Solver instance identification number.  
    nodenbr   - *INTS* where to store number of nodes.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If graph hasn't been built.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
                          or *nodenbr* is *NULL* (can occur in C).
    
  Fortran interface:
  >
  > SUBROUTINE MURGE_GETLOCALNODENBR(ID, NODENBR, IERROR)
  >   INTS, INTENT(IN)  :: ID
  >   INTS, INTENT(OUT) :: NODENBR
  >   INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETLOCALNODENBR 
*/
INTS MURGE_GetLocalNodeNbr(INTS id, INTS *nodenbr);

/*
  Function: MURGE_GetLocalNodeList
  
  Computes the local node list, corresponding to 
  the new distribution, after preprocessing.
  
  *nodelist* array has to be allocated before calling
  this function.

  As it's result determines the size of *nodelist*
  array, <MURGE_GetLocalNodeNbr> should be run before it.
    
  Parameters:
    id        - Solver instance identification number.  
    nodelist  - Array where to store the list of local nodes.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - if <MURGE_GetLocalNodeNbr> has not been called 
                          before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
                          or *nodelist* is *NULL* (can occur in C).
    
  Fortran interface:
  >
  > SUBROUTINE MURGE_GETLOCALNODELIST(ID, NODELIST, IERROR)
  >   INTS,               INTENT(IN)  :: ID      
  >   ! Warning : 0 is not the size of the array.
  >   ! Writing DIMENSION(:) does not work with
  >   ! the C function call (fortran send the array size?)
  >   INTS, DIMENSION(0), INTENT(OUT) :: NODELIST
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETLOCALNODELIST
*/
INTS MURGE_GetLocalNodeList(INTS id, INTS *nodelist);

/*
  Function: MURGE_GetLocalUnknownNbr
  
  Computes preprocessing step, if not done, and the number of
  Unkowns in the new ditribution of the matrix.

  Needs the graph to be built.

  Parameters:
    id        - Solver instance identification number.  
    unkownnbr     - *INTS* where to store number of unkowns.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If graph hasn't been built.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
                          or *unkownnbr* is *NULL* (can occur in C).
    
  Fortran interface:
  >
  > SUBROUTINE MURGE_GETLOCALUNKOWNNBR(ID, UNKOWNNBR, IERROR)
  >   INTS, INTENT(IN)  :: ID
  >   INTS, INTENT(OUT) :: UNKOWNNBR
  >   INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETLOCALUNKOWNNBR 
*/
INTS MURGE_GetLocalUnknownNbr(INTS id, INTS *unkownnbr);

/*
  Function: MURGE_GetLocalUnknownList
  
  Computes the local unkown list, corresponding to 
  the new distribution, after preprocessing.
  
  *unkownlist* array has to be allocated before calling
  this function.

  As it's result determines the size of *unkownlist*
  array, <MURGE_GetLocalUnkownNbr> should be run before it.
    
  Parameters:
    id        - Solver instance identification number.  
    unkownlist  - Array where to store the list of local unkowns.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - if <MURGE_GetLocalUnkownNbr> has not been called 
                          before.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
                          or *unkownlist* is *NULL* (can occur in C).
    
  Fortran interface:
  >
  > SUBROUTINE MURGE_GETLOCALUNKOWNLIST(ID, UNKOWNLIST, IERROR)
  >   INTS,               INTENT(IN)  :: ID      
  >   INTS, DIMENSION(0), INTENT(OUT) :: UNKOWNLIST
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETLOCALUNKOWNLIST
*/
INTS MURGE_GetLocalUnknownList(INTS id, INTS *unkownlist);

/* 
   Group: Filling the matrix   
*/

/*
  Function: MURGE_AssemblyBegin

  Begin Filling up sequence for the matrix, will allocate
  temporary structures used to build the matrix.

  It will perform preprocessing if it has not been done yet.

  It needs graph to be built.
  
  Several assembly sequences can be performed on the same matrix.

  *mode* shouldn't be *MURGE_ASSEMBLY_RESPECT* if neither 
  <MURGE_GetLocalNodeList> nor <MURGE_GetLocalUnknownList> has been called.

  Parameters: 
    id      - Solver instance identification number.  
    coefnbr - Number of coeficients he will had.
    op      - Operation to perform for coefficient which appear
              several tim (see <MURGE_ASSEMBLY_OP>).
    op2     - Operation to perform when a coefficient is set by
              two different processors (see <MURGE_ASSEMBLY_OP>).
    mode    - Indicates if user ensure he will respect solvers distribution 
              (see <MURGE_ASSEMBLY_MODE>).
    sym     - Indicates if user will give coefficient in a symmetric way 
              (ie: only triangullar part) or not. 

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If graph hasn't been built before.
    MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
                          *op*, *mode*, *sym*, or *coefnbr* are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_ASSEMBLYBEGIN(ID, COEFNBR, OP, OP2, MODE, SYM, IERROR)
  >   INTS,      INTENT(IN)  :: ID, OP, OP2, MODE, SYM
  >   INTL,      INTENT(IN)  :: COEFNBR
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_ASSEMBLYBEGIN
*/ 
INTS MURGE_AssemblyBegin(INTS id, INTL coefnbr, INTS op, INTS op2, INTS mode, INTS sym);

/*
  Function: MURGE_AssemblySetValue

  Set a coefficient value in the matrix.
  
  Needs <MURGE_AssemblyBegin> to be called before.

  Parameters: 
    id      - Solver instance identification number.  
    ROW     - Global row number of the coefficient.
    COL     - Global column number of the coefficient.
    value   - value of the coefficient.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If we are not in an assembly section.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
                          *ROW* or *COL* are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_ASSEMBLYSETVALUE(ID, ROW, COL, VALUE, IERROR)
  >   INTS,      INTENT(IN)  :: ID, ROW, COL
  >   COEF,      INTENT(IN)  :: VALUE
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_ASSEMBLYSETVALUE
*/ 
INTS MURGE_AssemblySetValue(INTS id, INTS ROW, INTS COL, COEF value);

/*
  Function: MURGE_AssemblySetNodeValues

  Set coefficients value for a node in the matrix.
  
  Needs <MURGE_AssemblyBegin> to be called before.

  Parameters: 
    id      - Solver instance identification number.  
    ROW     - Global row number of the coefficient.
    COL     - Global column number of the coefficient.
    values  - value of the coefficient.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If we are not in an assembly section.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
                          *ROW* or *COL* are not valid, or, in C, if 
			  *values* is NULL.

  Fortran interface:
  >
  > SUBROUTINE MURGE_ASSEMBLYSETNODEVALUES(ID, ROW, COL, VALUES, IERROR)
  >   INTS,               INTENT(IN)  :: ID, ROW, COL
  >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_ASSEMBLYSETNODEVALUES
*/ 
INTS MURGE_AssemblySetNodeValues(INTS id, INTS ROW, INTS COL, COEF *values);

/*
  Function: MURGE_AssemblySetBlockValues

  Set coefficients value for a dens block in the matrix.
  
  Needs <MURGE_AssemblyBegin> to be called before.

  Parameters: 
    id      - Solver instance identification number.  
    nROW    - Number of rows in the dense matrix.
    ROWlist - List of global row numbers.
    nCOL    - Number of columns in the dense matrix.
    COLlist - List of global column numbers.
    values  - Values array, by column (Fortran style)
    
  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If we are not in an assembly section.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
                          *nI* or *nJ* are not valid, or, in C, if 
			  *Ilist*, *Jlist* or *values* is NULL.

  Fortran interface:
  >
  > SUBROUTINE MURGE_ASSEMBLYSETBLOCKVALUES(ID, NROW, ROWLIST, &
  >                                 & NCOL, COLLIST, VALUES, IERROR)
  >   INTS,               INTENT(IN)  :: ID, NROW, NCOL
  >   INTS, DIMENSION(0), INTENT(IN)  :: ROWLIST
  >   INTL, DIMENSION(0), INTENT(IN)  :: COLLIST
  >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_ASSEMBLYSETBLOCKVALUES
*/ 
INTS MURGE_AssemblySetBlockValues(INTS id, INTS nROW, INTS *ROWlist, INTS nCOL, 
				  INTS *COLlist, COEF *values);
/*
  Function: MURGE_AssemblyEnd

  End Filling up sequence for the matrix.

  Needs <MURGE_AssemblyBegin> to be called before.

  Parameters: 
    id      - Solver instance identification number.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_ORDER     - If we are not in an assembly section.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_ASSEMBLYEND(ID, IERROR)
  >   INTS,      INTENT(IN)  :: ID
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_ASSEMBLYEND
*/ 
INTS MURGE_AssemblyEnd(INTS id);

/*
  Function: MURGE_MatrixReset

  Reset the matrix structure.

  Parameters: 
    id      - Solver instance identification number.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_MATRIXRESET(ID, IERROR)
  >   INTS,      INTENT(IN)  :: ID
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_MATRIXRESET
*/ 
INTS MURGE_MatrixReset(INTS id);

/*
  Function: MURGE_MatrixGlobalCSR

  Add the given global Compress Sparse Row matrix to the matrix.

  Parameters: 
    id      - Solver instance identification number.
    N       - Number of columns.
    rowptr  - Index of the first element of each row in *COLS* and 
              *values* array.
    COLS    - Column number array.
    values  - values array.
    root    - Root processor for MPI communications.
    op      - Operation to perform if a coefficient appear twice
              (see <MURGE_ASSEMBLY_OP>).
    sym     - Indicates if user will give coefficient in a symmetric way 
              (ie: only triangullar part) or not.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, 
                          if *root*, *op* or the CSR are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_MATRIXGLOBALCSR(ID, N, ROWPTR, COLS, VALUES, &
  >                                & ROOT, OP, SYM,  IERROR)
  >   INTS,               INTENT(IN)  :: ID, N, ROOT, OP, SYM
  >   INTL, DIMENSION(0), INTENT(IN)  :: ROWPTR
  >   INTS, DIMENSION(0), INTENT(IN)  :: COLS
  >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_MATRIXGLOBALCSR
*/ 
INTS MURGE_MatrixGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS, 
			   COEF *values, INTS root, INTS op, INTS sym);

/*
  Function: MURGE_MatrixGlobalCSC

  Add the given global Compress Sparse Column matrix to the matrix.

  Parameters: 
    id      - Solver instance identification number.
    N       - Number of columns.
    colptr  - Index of the first element of each column in *ROWS* and 
              *values* array.
    ROWS    - Row number array.
    values  - values array.
    root    - Root processor for MPI communications.
    op      - Operation to perform if a coefficient appear twice
              (see <MURGE_ASSEMBLY_OP>).
    sym     - Indicates if user will give coefficient in a symmetric way 
              (ie: only triangullar part) or not.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, 
                          if *root*, *op* or the CSC are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_MATRIXGLOBALCSC(ID, N, COLPTR, ROWS, &
  >                              & VALUES, ROOT, OP, SYM, IERROR)
  >   INTS,               INTENT(IN)  :: ID, N, ROOT, OP, SYM
  >   INTL, DIMENSION(0), INTENT(IN)  :: COLPTR
  >   INTS, DIMENSION(0), INTENT(IN)  :: ROWS
  >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_MATRIXGLOBALCSC
*/ 
INTS MURGE_MatrixGlobalCSC(INTS id, INTS N, INTL *COLPTR, INTS *ROWS, 
			   COEF *values, INTS root, INTS op, INTS sym);

/*
  Function: MURGE_MatrixGlobalIJV

  Add the given global Compress Sparse Column matrix to the matrix.

  Parameters: 
    id      - Solver instance identification number.
    N       - Number of edges.
    NNZ     - Number of non zeros.
    ROWS    - Global row number array.
    COLS    - Global column number array.
    values  - values array.
    root    - Root processor for MPI communications.
    op      - Operation to perform if a coefficient appear twice
              (see <MURGE_ASSEMBLY_OP>).
    sym     - Indicates if user will give coefficient in a symmetric way 
              (ie: only triangullar part) or not.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, 
                          if *root*, *op*, *ROWS* or *COLS* are not valid.

  Fortran interface:
  >
  > SUBROUTINE MURGE_MATRIXGLOBALIJV(ID, N, NNZ, ROWS, COLS, VALUES, &
  >                                & ROOT, OP, SYM, IERROR)
  >   INTS,               INTENT(IN)  :: ID, ROOT, OP, SYM, N
  >   INTL,               INTENT(IN)  :: NNZ
  >   INTS, DIMENSION(0), INTENT(IN)  :: ROWS, COLS
  >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_MATRIXGLOBALIJV
*/ 
INTS MURGE_MatrixGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROWS, INTS *COLS, 
			   COEF *values, INTS root, INTS op, INTS sym);

/* 
   Group: Filling the right-hand-side member
   
*/
/*
  Function: MURGE_SetGlobalRHS

  Set the right-hand-side member in global mode.

  Parameters: 
    id      - Solver instance identification number.
    b       - Array of size global column number which correspond to the 
              right-hand-side member.
    op      - Operation to perform if a coefficient appear twice
              (see <MURGE_ASSEMBLY_OP>).
    root    - Indicates which processor has the right-hand-side member,
              -1 for all.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_SETGLOBALRHS(ID, B, OP, ROOT, IERROR)
  >   INTS,               INTENT(IN)  :: ID, OP, ROOT
  >   COEF, DIMENSION(0), INTENT(IN)  :: B
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_SETGLOBALRHS
*/
INTS MURGE_SetGlobalRHS(INTS id, COEF *b, INTS root, INTS op); 

/*
  Function: MURGE_SetLocalRHS

  Set the right-hand-side member in local mode.

  Parameters: 
    id      - Solver instance identification number.
    b       - Array of size local column number which correspond to the 
              right-hand-side member.
    op      - Operation to perform if a coefficient appear twice
              (see <MURGE_ASSEMBLY_OP>).
    op2     - Operation to perform when a coefficient is set by
              two different processors (see <MURGE_ASSEMBLY_OP>).

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_SETLOCALRHS(ID, B, OP, OP2, IERROR)
  >   INTS,               INTENT(IN)  :: ID, OP, OP2
  >   COEF, DIMENSION(0), INTENT(IN)  :: B
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_SETLOCALRHS
*/
INTS MURGE_SetLocalRHS(INTS id, COEF *b, INTS op, INTS op2); 

/*
  Function: MURGE_SetRHS

  Set the right-hand-side member, giving the list of 
  coefficient that we set.

  *mode* shouldn't be *MURGE_ASSEMBLY_RESPECT* if neither 
  <MURGE_GetLocalNodeList> nor <MURGE_GetLocalUnknownList> has been called.

  Parameters: 
    id       - Solver instance identification number.
    n        - Number of coefficients to set.
    coefsidx - List of global index of the coefficients to set.
    B        - Array of coefficients values. 
    op       - Operation to perform if a coefficient appear twice
               (see <MURGE_ASSEMBLY_OP>).
    op2     - Operation to perform when a coefficient is set by
              two different processors (see <MURGE_ASSEMBLY_OP>).
    mode     - Indicates if user ensure he will respect solvers distribution 
               (see <MURGE_ASSEMBLY_MODE>).

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
                          if *mode* or *op* are not valid, or, in C, 
			  if *coeflist* or *b* are NULL.

  Fortran interface:
  >
  > SUBROUTINE MURGE_SETRHS(ID, N, COEFSIDX, B, OP, OP2, MODE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, N, OP, OP2, MODE
  >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
  >   COEF, DIMENSION(0), INTENT(IN)  :: B
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_SETRHS
*/
INTS MURGE_SetRHS(INTS id, INTS n, INTS *coefsidx, COEF *b, 
		  INTS op, INTS op2, INTS mode); 

/*
  Function: MURGE_RHSReset

  Reset the right-hand-side.

  Parameters: 
    id      - Solver instance identification number.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_RHSRESET(ID, IERROR)
  >   INTS,      INTENT(IN)  :: ID
  >   INTS,      INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_RHSRESET
*/ 
INTS MURGE_RHSReset(INTS id);

/* 
   Group: Getting the solution
   
*/
/*
  Function: MURGE_GetGlobalSolution

  Perform Factorization and Solve, if needed, 
  and then fill the global solution in *x*.

  Parameters: 
    id      - Solver instance identification number.
    x       - Array of size global column number*dof which will contain
              the solution
    root    - Indicates which processor will have the solution 
              at the end of the call, -1 for all.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GETGLOBALSOLUTION(ID, X, ROOT, IERROR)
  >   INTS,               INTENT(IN)  :: ID, ROOT
  >   COEF, DIMENSION(0), INTENT(OUT) :: X
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETGLOBALSOLUTION
*/
INTS MURGE_GetGlobalSolution(INTS id, COEF *x, INTS root); 

/*
  Function: MURGE_GetLocalSolution

  Perform Factorization and Solve, if needed, 
  and then fill the local solution in *x*.

  Parameters: 
    id      - Solver instance identification number.
    x       - Array of size local column number*dof which will contain
              the solution

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GETLOCALSOLUTION(ID, X, IERROR)
  >   INTS,               INTENT(IN)  :: ID
  >   COEF, DIMENSION(0), INTENT(OUT) :: X
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETLOCALSOLUTION
*/
INTS MURGE_GetLocalSolution(INTS id, COEF *x); 

/*
  Function: MURGE_GetSolution

  Perform Factorization and Solve, if needed, 
  and then fill the solution in *x* followin the given
  index list.

  Parameters: 
    id       - Solver instance identification number.
    n        - Number of coefficients user wants to get.
    coefsidx - List of the coefficients user wants to get.
    x        - Array of size dof*n which will contain
               the solution.
    mode     - Indicates if the user is sure to respect the distribution.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GETSOLUTION(ID, N, COEFSIDX, X, MODE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, MODE, N
  >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
  >   COEF, DIMENSION(0), INTENT(OUT) :: X
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETSOLUTION
*/
INTS MURGE_GetSolution(INTS id,  INTS n, INTS *coefsidx, COEF *x, INTS mode); 

/*
  Group: Scaling
*/

/*
  Function: MURGE_GetGlobalNorm

  Compute the global norm array following a norm rule.

  Must be performed after assembly step.

  Parameters: 
    id      - Solver instance identification number.
    norm    - Array of size global column number*dof which will contain
              the norm values
    root    - Indicates which processor will have the norm array
              at the end of the call, -1 for all.
    rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_ORDER     - If the assembly has not been performed.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GETGLOBALNORM(ID, NORM, ROOT, RULE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, ROOT, RULE
  >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETGLOBALNORM
*/
INTS MURGE_GetGlobalNorm(INTS id, REAL *norm, INTS root, INTS rule); 

/*
  Function: MURGE_GetLocalNorm

  Compute the local norm array following a norm rule.
  
  Must be performed after assembly step.

  Parameters: 
    id      - Solver instance identification number.
    norm    - Array of size local column number*dof which will contain
              the solution
    rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_ORDER     - If the assembly has not been performed.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GETLOCALNORM(ID, NORM, RULE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, RULE
  >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETLOCALNORM
*/
INTS MURGE_GetLocalNorm(INTS id, REAL *norm, INTS rule); 

/*
  Function: MURGE_GetNorm

  Compute the indicated part of the norm array 
  following a norm rule.
  
  Must be performed after assembly step.
 

  Parameters: 
    id       - Solver instance identification number.
    n        - Number of coefficients user wants to get norm of.
    coefsidx - List of the coefficients user wants to get norm of.
    norm     - Array of size dof*n which will contain
               the solution.
    rule     - Rule to follow to build norm array, see <MURGE_NORM_RULES>
    mode     - Indicates if the user is sure to respect the distribution.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_ORDER     - If the assembly has not been performed.

  Fortran interface:
  >
  > SUBROUTINE MURGE_GETNORM(ID, N, COEFSIDX, NORM, RULE, MODE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, MODE, N, RULE
  >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
  >   COEF, DIMENSION(0), INTENT(OUT) :: NORM
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETNORM
*/
INTS MURGE_GetNorm(INTS id,  INTS n, INTS *coefsidx, REAL *norm, INTS rule, INTS mode); 

/* 
   Function: MURGE_ApplyGlobalScaling
   
   Apply scaling to local unknowns.

   Must be performed after assembly step.

   Parameters:
     id      - Solver instance identification number.
     scal    - Scaling user wants to apply.
     root    - Indicates which processor that posses the scaling array, 
               -1 for all.
     sc_mode - Indicate if the scaling is applied on rows or on columns.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_ORDER     - If the assembly has not been performed.

  Fortran interface:
  >
  > SUBROUTINE MURGE_APPLYGLOBALSCALING(ID, SCAL, ROOT, SC_MODE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, ROOT, SC_MODE
  >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_APPLYGLOBALSCALING

*/
INTS MURGE_ApplyGlobalScaling(INTS id, REAL *scal,  INTS root, INTS sc_mode);

/*
  Function: MURGE_ApplyLocalScaling

  Apply the local scaling array on the matrix.
  
  Must be performed after assembly step.

  Parameters: 
    id      - Solver instance identification number.
    scal    - Array of size local column number*dof which will contain
              the solution.
    sc_mode - Indicate if the scaling is applied on rows or on columns.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_ORDER     - If the assembly has not been performed.

  Fortran interface:
  >
  > SUBROUTINE MURGE_APPLYLOCALSCALING(ID, SCAL, SC_MODE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, SC_MODE
  >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_APPLYLOCALSCALING
*/
INTS MURGE_ApplyLocalScaling(INTS id, REAL *scal, INTS sc_mode);

/*
  Function: MURGE_ApplyScaling

  Apply the scaling array on the indicated part of the matrix
  
  Must be performed after assembly step.
 

  Parameters: 
    id       - Solver instance identification number.
    n        - Number of coefficients user wants to scale.
    coefsidx - List of the coefficients user wants to scale.
    scal     - Array of size dof*n which will contain
               the solution.
    sc_mode  - Indicate if the scaling is applied on rows or on columns.
    mode     - Indicates if the user is sure to respect the distribution.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
    MURGE_ERR_ORDER     - If the assembly has not been performed.

  Fortran interface:
  >
  > SUBROUTINE MURGE_APPLYSCALING(ID, N, COEFSIDX, SCAL, SC_MODE, MODE, IERROR)
  >   INTS,               INTENT(IN)  :: ID, SC_MODE, MODE, N
  >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
  >   COEF, DIMENSION(0), INTENT(OUT) :: SCAL
  >   INTS,               INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_APPLYSCALING
*/
INTS MURGE_ApplyScaling(INTS id,  INTS n, INTS *coefsidx, REAL *scal, INTS sc_mode, INTS mode);

/* 
   Group: Cleaning up this mess
   
*/
/*
  Function: MURGE_Clean

  Clean the given instance of the solver structure's.

  Parameters: 
    id      - Solver instance identification number.

  Returns: 
    MURGE_SUCCESS       - If function runned successfully.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.

  Fortran interface:
  >
  > SUBROUTINE MURGE_CLEAN(ID, IERROR)
  >   INTS, INTENT(IN)  :: ID
  >   INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_CLEAN
*/
INTS MURGE_Clean(INTS id);

/*
  Function: MURGE_Finalize
  
  Clean all not cleaned instances and instances ID array.
  
  Returns: 
    MURGE_SUCCESS       - If function runned successfully.

  Fortran interface:
  >
  > SUBROUTINE MURGE_FINALIZE(IERROR)
  >   INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_FINALIZE
*/
INTS MURGE_Finalize();

/* 
  Function: MURGE_GetSolver
   
  Return the solver ID Murge was compiled with.
   
  Parameters:
    solver - Integer to store solver ID.
   
  Returns:
    MURGE_SUCCESS - If execution ended normaly.
   
  Fortran interface:
  >
  > SUBROUTINE MURGE_GETSOLVER(SOLVER, IERROR)
  >   INTS,  INTENT(OUT) :: SOLVER
  >   INTS,  INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETSOLVER   
*/
INTS MURGE_GetSolver(INTS * solver);

/*
  Group: Getting Murge's metrics
 */
/*
  Function: MURGE_GetMetricINT

  Get an integer metric from MURGE.

  See <MURGE_IINFOS> and the solver documentation to 
  get available metric list.
  
  Parameters:
    id     - Solver instance identification number.
    metric - Wanted metric identification number.
    value  - Integer which will contain the value of the metric.

  Returns:
    MURGE_SUCCESS - If execution ended normaly.
    MURGE_ERR_ORDER     - If metric is not available in the current 
                          solver state.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
                          *metric* or *value* are not valid.

    
  Fortran interface:
  > SUBROUTINE MURGE_GETINFOINT(ID, INFO, VALUE, IERROR)
  >  INTS, INTENT(IN)  :: ID, INFO
  >  INTL, INTENT(OUT) :: VALUE
  >  INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETINFOINT
 */
INTS MURGE_GetInfoINT(INTS id,  INTS metric, INTL * value);

/*
  Function: MURGE_GetInfoREAL

  Get a real metric value from MURGE.

  See <MURGE_RINFOS> and the solver documentation to 
  get available metric list.
  
  Parameters:
    id     - Solver instance identification number.
    metric - Wanted metric identification number.
    value  - Real which will contain the value of the metric.

  Returns:
    MURGE_SUCCESS - If execution ended normaly.
    MURGE_ERR_ORDER     - If metric is not available in the current 
                          solver state.
    MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
                          *metric* or *value* are not valid.

    
  Fortran interface:
  > SUBROUTINE MURGE_GETINFOREAL(ID, INFO, VALUE, IERROR)
  >  INTS, INTENT(IN)  :: ID, INFO
  >  REAL, INTENT(OUT) :: VALUE
  >  INTS, INTENT(OUT) :: IERROR
  > END SUBROUTINE MURGE_GETINFOREAL
 */
INTS MURGE_GetInfoREAL(INTS id,  INTS metric, REAL * value);

/*
  Function: MURGE_PrintError
  
  Print the error message corresponding to ierror
  Parameters: 
    error_number  - Error identification number.

  Returns: 
    MURGE_ERR_PARAMETER - If ierror does not match an error number
    MURGE_SUCCESS       - If function runned successfully.

  Fortran interface:
  >
  > SUBROUTINE MURGE_PRINTERROR(ERROR_NUMBER, IERROR)
  >   INTS, INTENT(IN)  :: IERROR
  >   INTS, INTENT(OUT) :: ERROR_NUMBER
  > END SUBROUTINE MURGE_PRINTERROR
*/
INTS MURGE_PrintError(INTS error_number);

/*
  Function: MURGE_ExitOnError
  
  Print the error message corresponding to ierror.
  If the ierr is not MURGE_SUCCESS then the program is stopped.

  Parameters: 
    ierror         - Error identification number.

  Returns: 
    MURGE_SUCCESS   - If function runned successfully; stop the program otherwise.

  Fortran interface:
  >
  > SUBROUTINE MURGE_EXITONERROR(ERROR_NUMBER, IERROR)
  >   INTS, INTENT(IN)  :: IERROR
  >   INTS, INTENT(OUT) :: ERROR_NUMBER
  > END SUBROUTINE MURGE_EXITONERROR
*/
INTS MURGE_ExitOnError(INTS error_number);

/*
  Group: Murge's constants
 */
/*
  Enum: MURGE_RETURNS
  
  Murge error return values.
  
  Contains:
    MURGE_SUCCESS             - If function runs correctly.
    MURGE_ERR_ALLOCATE        - If couldn't allocate.
    MURGE_ERR_IO              - If there was an input or output error.
    MURGE_ERR_PARAMETER       - If one parameter is not correct.
    MURGE_ERR_ORDER           - If function were run in wrong order.
    MURGE_ERR_SOLVER          - Internal solver error.
    MURGE_ERR_NOT_IMPLEMENTED - Not yet implemented.
*/
enum MURGE_RETURNS {
  MURGE_SUCCESS             = 1,
  MURGE_ERR_ALLOCATE        = 2,
  MURGE_ERR_IO              = 3,
  MURGE_ERR_PARAMETER       = 4,
  MURGE_ERR_ORDER           = 5,
  MURGE_ERR_SOLVER          = 6,
  MURGE_ERR_NOT_IMPLEMENTED = 7
};
  

/*
  Enum: MURGE_IPARAM

  Murge integer parameters identifiers.
  
  Solvers may implement is own list of parameters.

  MURGE_IPARAM_BASEVAL - Numbering style , 0 for C, 1 for fortran. 
  MURGE_IPARAM_DOF     - Number of degrees of freedom.
*/
enum MURGE_IPARAM {
  MURGE_IPARAM_BASEVAL = 1024,
  MURGE_IPARAM_DOF,
  MURGE_IPARAM_SYM
};

/*
  Enum: MURGE_RPARAM

  Murge real parameters identifiers.
  
  Solvers may implement is own list of parameters.

  Contains: 
    MURGE_RPARAM_EPSILON_ERROR - Wanted norm error at the end of solve.
*/
enum MURGE_RPARAM {
  MURGE_RPARAM_EPSILON_ERROR = 1024
};
/*
  Enum: MURGE_IINFO

  Murge integer metrics identifiers.
  
  Solvers may implement is own list of parameters.
  
  Contains:
    MURGE_IINFOS_NNZ - Number of non zeros in factorized matrix.
*/
enum MURGE_IINFO {
  MURGE_IINFO_NNZ = 1024
};

/*
  Enum: MURGE_RINFO

  Murge real metrics identifiers.
  
  Solvers may implement is own list of parameters.

  Contains: 
    MURGE_RINFO_FACT_TIME  - Factorization time.
    MURGE_RINFO_SOLVE_TIME - Solving time.
*/
enum MURGE_RINFO {
  MURGE_RINFO_FACT_TIME = 1024,
  MURGE_RINFO_SOLVE_TIME
};

/*
  Enum: MURGE_ASSEMBLY_MODE

    Indicates if user can ensure that the information he is giving respects
    the solver distribution.
  
    MURGE_ASSEMBLY_RESPECT - User ensure he respects distribution 
                             during assembly.
			     See solver documentation.
    MURGE_ASSEMBLY_FOOL    - User is not sure he will respect ditribution
                             during assembly
*/
enum MURGE_ASSEMBLY_MODE {
  MURGE_ASSEMBLY_RESPECT,
  MURGE_ASSEMBLY_FOOL
};

/*
  Enum: MURGE_ASSEMBLY_OP
   
    Operations possible when a coefficient appear twice.

    MURGE_ASSEMBLY_ADD - Coefficients will be added during assembly.
    MURGE_ASSEMBLY_OVW - Coefficients will be overwriten during assembly.
    MURGE_ASSEMBLY_MAX - Maximum value will be used for assembly.
    MURGE_ASSEMBLY_MIN - Minimum value will be used for assembly.
*/
enum MURGE_ASSEMBLY_OP {
  MURGE_ASSEMBLY_ADD,
  MURGE_ASSEMBLY_OVW,
  MURGE_ASSEMBLY_MAX,
  MURGE_ASSEMBLY_MIN
};

/*
  Enum: MURGE_SOLVER
  
  Solver ID for murge compliant solvers.
  
  Contains: 
    MURGE_SOLVER_HIPS   - HIPS hybrid solver.
    MURGE_SOLVER_PASTIX - PaStiX direct solver.
*/
enum MURGE_SOLVER {
  MURGE_SOLVER_HIPS,
  MURGE_SOLVER_PASTIX
};

/*
  Enum: MURGE_BOOLEAN

  Boolean for murge parameters

  Contains: 
    MURGE_BOOLEAN_FALSE - False value  
    MURGE_BOOLEAN_TRUE  - True value
*/

enum MURGE_BOOLEAN {
  MURGE_BOOLEAN_FALSE,
  MURGE_BOOLEAN_TRUE
};

/* 
   Enum: MURGE_NORM_RULES
   
   Flags for Murge's norm rules

   Contains:
     MURGE_NORM_MAX_COL  - Get maximum column value  (absolute value).
     MURGE_NORM_MAX_ROW  - Get maximum row value     (absolute value).
     MURGE_NORM_2_COL    - Get the norm 2 of columns.
     MURGE_NORM_2_ROW    - Get the norm 2 of rows.
*/
enum MURGE_NORM_RULES {
  MURGE_NORM_MAX_COL,
  MURGE_NORM_MAX_ROW,
  MURGE_NORM_2_COL,
  MURGE_NORM_2_ROW
};

/* 
   Enum: MURGE_SCAL_MODES
   
   Flags for Murge's scaling rules

   Contains:
     MURGE_SCAL_COL  - Perform scaling on columns
     MURGE_SCAL_ROW  - Perform scaling on rows.
*/
enum MURGE_SCAL_MODES {
  MURGE_SCAL_COL,
  MURGE_SCAL_ROW
};

#endif /* MURGE_H */
/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
**
** This file is part of the PaStiX parallel sparse matrix package.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/******************************************************************************
 * Title: PaStiX specific addons to Murge                                     *
 *                                                                            *
 * Functions declaration for Murge function defined only in PaStiX.           *
 *                                                                            *
 * Authors:                                                                   *
 *   Xavier Lacoste - xavier.lacoste@inria.fr                                 *
 *                                                                            *
 * More informations can be found in <murge.h> documentation.                 *
 *                                                                            *
 ******************************************************************************/


/******************************************************************************
 * Function: MURGE_Analyze                                                    *
 *                                                                            *
 * Perform matrix analyze.                                                    *
 *                                                                            *
 * Follows several steps :                                                    *
 *   - Compute a new ordering of the unknows                                  *
 *   - Compute the symbolic factorisation of the matrix                       *
 *   - Distribute column blocks and computation on processors                 *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Analyze(INTS id);

/******************************************************************************
 * Function: MURGE_Factorize                                                  *
 *                                                                            *
 * Perform matrix factorization.                                              *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Factorize(INTS id);

/******************************************************************************
 * Function: MURGE_SetOrdering                                                *
 *                                                                            *
 * Set a personal ordering to perform factorization.                          *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   permutation - Permutation to perform factorization.                      *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_SetOrdering(INTS   id,
                       INTS * permutation);

/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeNbr                                     *
 *                                                                            *
 * Set local node number if users only perform product and doesn't need       *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   n  - Number of local nodes.                                              *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ProductSetLocalNodeNbr (INTS id, INTS n);

/******************************************************************************
 * Function: MURGE_ProductSetGlobalNodeNbr                                    *
 *                                                                            *
 * Set global node number if users only perform product and doesn't need      *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   N  - Number of global nodes.                                             *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ProductSetGlobalNodeNbr (INTS id, INTS N);

/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeList                                    *
 *                                                                            *
 * Set local node list if users only perform product and doesn't need         *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id  - Solver instance identification number.                             *
 *   l2g - Local to global node numbers.                                      *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ProductSetLocalNodeList (INTS id, INTS * l2g);

/******************************************************************************
 * Function: MURGE_GetLocalProduct                                            *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <MURGE_SetLocalRHS> or              *
 * <MURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   x  - Array in which the local part of the product will be stored.        *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_GetLocalProduct (INTS id, COEF *x);

/******************************************************************************
 * Function: MURGE_GetGlobalProduct                                           *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <MURGE_SetLocalRHS> or              *
 * <MURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id   - Solver instance identification number.                            *
 *   x    - Array in which the product will be stored.                        *
 *   root - Rank of the process which will own the product at end of call,    *
 *          use -1 for all processes.                                         *
 * Returns:                                                                   *
 *   MURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_GetGlobalProduct (INTS id, COEF *x, INTS root);

/******************************************************************************
 * Function: MURGE_ForceNoFacto                                               *
 *                                                                            *
 * Prevent Murge from running factorisation even if matrix has changed.       *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 * Returns:                                                                   *
 *   MURGE_SUCCESS                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ForceNoFacto(INTS id);

/******************************************************************************
 * Function: MURGE_SetLocalNodeList                                           *
 *                                                                            *
 * Deprecated, need to be checked                                             *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_SetLocalNodeList(INTS id, INTS n, INTS * list);

/******************************************************************************
 * Function: MURGE_AssemblySetSequence                                        *
 *                                                                            *
 * Create a sequence of entries to build a matrix and store it for being      *
 * reused.                                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id      - Solver instance identification number.                         *
 *   coefnbr - Number of entries.                                             *
 *   ROWs    - List of rows in the sequence.                                  *
 *   COLs    - List of columns in the sequence.                               *
 *   op      - Operation to perform for coefficient which appear              *
 *             several tim (see <MURGE_ASSEMBLY_OP>).                         *
 *   op2     - Operation to perform when a coefficient is set by              *
 *             two different processors (see <MURGE_ASSEMBLY_OP>).            *
 *   mode    - Indicates if user ensure he will respect solvers distribution  *
 *             (see <MURGE_ASSEMBLY_MODE>).                                   *
 *   nodes   - 0 entries are entered value by value,                          *
 *             1 entries are entries node by node.                            *
 *   id_seq  - Sequence ID.                                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.   *
 ******************************************************************************/
int MURGE_AssemblySetSequence (INTS id, INTL coefnbr, INTS * ROWs, INTS * COLs,
                               INTS op, INTS op2, INTS mode, INTS nodes,
                               INTS * id_seq);

/******************************************************************************
 * MURGE_AssemblySetSequence                                                  *
 *                                                                            *
 * Assembly the matrix using a stored sequence.                               *
 *                                                                            *
 * Parameters:                                                                *
 *   id      - Solver instance identification number.                         *
 *   id_seq  - Sequence ID.                                                   *
 *   values  - Values to insert in the CSC.                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *id_seq* or *values* are not valid.                *
 ******************************************************************************/
INTS MURGE_AssemblyUseSequence(INTS id, INTS id_seq, COEF * values);

/******************************************************************************
 * Function: MURGE_AssemblyDeleteSequence                                     *
 *                                                                            *
 * Destroy an assembly sequence                                               *
 *                                                                            *
 *   id      - Solver instance identification number.                         *
 *   id_seq  - Sequence ID.                                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *id_seq* is not valid.                             *
 ******************************************************************************/
INTS MURGE_AssemblyDeleteSequence(INTS id, INTS id_seq);
