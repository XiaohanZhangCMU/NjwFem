! 
! **
! ** This file is part of the PaStiX parallel sparse matrix package.
! **
! ** This software is governed by the CeCILL-C license under French law
! ** and abiding by the rules of distribution of free software. You can
! ** use, modify and/or redistribute the software under the terms of the
! ** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
! ** URL: "http://www.cecill.info".
! ** 
! ** As a counterpart to the access to the source code and rights to copy,
! ** modify and redistribute granted by the license, users are provided
! ** only with a limited warranty and the software's author, the holder of
! ** the economic rights, and the successive licensors have only limited
! ** liability.
! ** 
! ** In this respect, the user's attention is drawn to the risks associated
! ** with loading, using, modifying and/or developing or reproducing the
! ** software by the user in light of its specific status of free software,
! ** that may mean that it is complicated to manipulate, and that also
! ** therefore means that it is reserved for developers and experienced
! ** professionals having in-depth computer knowledge. Users are therefore
! ** encouraged to load and test the software's suitability as regards
! ** their requirements in conditions enabling the security of their
! ** systems and/or data to be ensured and, more generally, to use and
! ** operate it in the same conditions as regards security.
! ** 
! ** The fact that you are presently reading this means that you have had
! ** knowledge of the CeCILL-C license and that you accept its terms.
! 
! 
!    Title: Murge
!      Murge is close to the Merge.
!   
!    Description:
! 
!      Murge is a common interface definition to multiple solver.
! 
!      It has been initiated by *HIPS* and *PaStiX* solvers developpers
!      in january 2009.
!  
!   Integers:
! 
!      Depending of your compilation options, *INTS* and *INTL* can be
!      32 or 64 bits.
!      
!      In C user will be abble to use *INTS* types.  In Fortran user
!      will have to use the correct type, corresponding to the
!      compilation options.
!    
!      If user doesn't define any compilation option, *INTS* and *INTL*
!      will both be int.
!    
!      If user defines *-DINTSIZE32*, *INTS* and *INTL* will both be 32
!      bits integers.
!    
!      If user defines *-DINTSIZE64*, *INTS* will be 32 bits integers
!      whereas *INTL* will be INT64 integers.
!      
!      If user defines *-DINTSSIZE64*, both *INTS* and *INTL* will be 64
!      bit integers.
!      
!      We recommand not to use *-DINTSSIZE64*, as *INTS* integers are
!      supposed to keep being small integers, and there is no need to
!      make them longs.
! 
!   Coefficients:
! 
!      The coefficients of the matrices and vectors is of type *COEF*.
!    
!      *COEF* can be simple or double precision and real or complex
!      depending on compilation options.
!    
!      If user defines nothing, *COEF* will be real simple precision.
!      
!      If user defines *-DPREC_DOUBLE*, *COEF* will be double
!      precision.
!      
!      Defining nothing or *-DPREC_SIMPLE* will result in setting *COEF*
!      precision to simple.
!    
!      If user defines *-DTYPE_COMPLEX*, *COEF* will be complex, otherwise,
!      if user defines nothing or *-DTYPE_REAL* it will be REAL.
! 
!   Floating parameters:
!     
!      Floating parameters are of type *REAL* which can be simple or 
!      double precision.
!      
!      As for Coefficients, *REAL* precision also depends on 
!      *-DPREC_DOUBLE* compilation option.
!      
!      Defining nothing or *-DPREC_SIMPLE* will result in setting *REAL*
!      precision to simple.
! 
!   Nodes and Unknowns:
! 
!      In murge exists node and unknown, a node correspond to a node 
!      in the physical problem, which can be composed of 
!      multiple degrees of freedom.
!      
!      Thus, the number of unknowns correspond to the number of node 
!      times the degree of freedom.
!      
! 
!   Authors:
! 
!      HIPS and PaStiX developpers teams :
! 
!        Mathieu Faverge   - faverge@labri.fr 
!        Jérémie Gaidamour - gaidamou@labri.fr 
!        Pascal  Hénon     - henon@labri.fr 
!        Xavier  Lacoste   - lacoste@labri.fr 
!        Pierre  Ramet     - ramet@labri.fr     
!        
! 
! 
!    Group: Solver setup functions
! 
! 
!   Function: MURGE_Initialize
! 
!   Allocate the instance arrays which will keeps intern data for all
!   solver instances.
! 
!   If user is creating several threads calling the solver, this function 
!   has to be called before creating threads to insure that the execution
!   is thread safe.
! 
!   Parameters: 
!     idnbr - Maximum number of solver instances that will be
!             launched.
! 
!   Returns: 
!     MURGE_SUCCESS - If function runned successfully.
!     MURGE_ERR_ALLOCATE - If for some reason, allocation was not
!                      successfull.
! 
! 
!   Fortran interface:
!   >
INTERFACE
   SUBROUTINE MURGE_INITIALIZE(IDNBR, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: IDNBR
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_INITIALIZE
   ! 
   ! 
   !   Function: MURGE_SetDefaultOptions
   ! 
   !   Sets default options, for solver instance number *id*.
   ! 
   !   The default option set correspond to *stratnum* strategy ID, 
   !   depending on the solver.
   ! 
   !   Needs <MURGE_Initialize> to be called before 
   !   to allocate solver instances array.
   !   
   !   Parameters: 
   !     id       - Solver instance identification number.
   !     stratnum - Strategy for the default option Set.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_Initialize> was not called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range 
   !                           or *stratnum* is not valid.
   !     MURGE_ERR_ALLOCATE  - If couldn't create solver instance.
   ! 
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SETDEFAULTOPTIONS(ID, STRATNUM, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: ID, STRATNUM
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETDEFAULTOPTIONS
   !  
   ! 
   !   Function: MURGE_SetOptionINT
   ! 
   !   Sets integer option, indicated by *number*, to *value* for the
   !   solver instance number *id*.
   ! 
   !   Needs <MURGE_SetDefaultOptions> to be called before to initiate
   !   solver instance data.
   ! 
   !   Parameters: 
   !     id     - Solver instance identification number.  
   !     number - Identification of the integer parameter.  
   !     value  - value to set the parameter to.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
   !                           *number* or *value* are not valid. 
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SETOPTIONINT(ID, NUMBER, VALUE, IERROR)
     INTEGER(KIND=4), INTENT(IN)    :: ID, NUMBER, VALUE
     INTEGER(KIND=4), INTENT(OUT)   :: IERROR
   END SUBROUTINE MURGE_SETOPTIONINT
   ! 
   ! 
   !   Function: MURGE_SetOptionREAL
   ! 
   !   Sets real option, indicated by *number*, to *value* for the
   !   solver instance number *id*.
   ! 
   !   Needs <MURGE_SetDefaultOptions> to be called before to initiate
   !   solver instance data.
   ! 
   !   Parameters: 
   !     id     - Solver instance identification number.  
   !     number - Identification of the integer parameter.  
   !     value  - value to set the parameter to.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
   !                           *number* or *value* are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SETOPTIONREAL(ID, NUMBER, VALUE, IERROR)
     INTEGER(KIND=4), INTENT(IN)    :: ID, NUMBER
     REAL(KIND=8),    INTENT(IN)    :: VALUE
     INTEGER(KIND=4), INTENT(OUT)   :: IERROR
   END SUBROUTINE MURGE_SETOPTIONREAL
   ! 
   ! 
   !   Function: MURGE_SetCommunicator
   ! 
   !   Sets MPI communicator for the given solver instance.
   ! 
   !   Needs <MURGE_SetDefaultOptions> to be called before to initiate
   !   solver instance data.
   ! 
   !   Musn't be called before <MURGE_SAVE>, <MURGE_LOAD>, 
   !   <MURGE_GetLocalNodeNbr> nor <MURGE_GetLocalUnknownNbr> 
   !   because the solver as to be runned with the same MPI 
   !   communicator all along.
   ! 
   !   If this function is not called, MPI communicator will be
   !   *MPI_COMM_WORLD*.
   ! 
   !   This function may not exist if the solver 
   !   has been compiled without MPI.
   !   
   !   Parameters: 
   !     id      - Solver instance identification number.  
   !     mpicomm - MPI communicator to affect the solver to.  
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before or if it is called after 
   ! 		          the solver starts its computing tasks.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SETCOMMUNICATOR(ID, MPICOMM, IERROR)
     INTEGER(KIND=4),            INTENT(IN)  :: ID
     INTEGER,         INTENT(IN)  :: MPICOMM
     INTEGER(KIND=4),            INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETCOMMUNICATOR
   ! 
   ! 
   !    Group: Graph setup function
   ! 
   ! 
   !   Function: MURGE_GraphBegin
   ! 
   !   Begin building the adjency graph for renumbering and all 
   !   preprocessing.
   ! 
   !   Needs <MURGE_SetDefaultOptions> to be called before to initiate
   !   solver instance data.
   !   
   !   Allocate temporary structures needed to build the graph.
   !   
   !   Parameters: 
   !     id      - Solver instance identification number.  
   !     N       - Global number of nodes in the graph.
   !     edgenbr - Number of edges which will be added in the graph by proc.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before.
   !     MURGE_ERR_ALLOCATE  - If Allocate didn't worked.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
   !                           *N* or *edgenbr* are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GRAPHBEGIN(ID, N, EDGENBR, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID, N
     INTEGER(KIND=4),      INTENT(IN)  :: EDGENBR
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GRAPHBEGIN
   ! 
   ! 
   !   Function: MURGE_GraphEdge
   ! 
   !   Adds an edge to the graph user is currently building.
   ! 
   !   Needs <MURGE_GraphBegin> to be called before to initiate
   !   building sequence.
   ! 
   !   This function depends on integer parameter *MURGE_BASEVAL*.
   !   
   !   Parameters: 
   !     id      - Solver instance identification number.  
   !     ROW     - First vertex of the edge.
   !     COL     - Second vertex of the edge.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_GraphBegin> was not 
   !                           called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
   !                           *I* or *J* are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GRAPHEDGE(ID, ROW, COL, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID, ROW, COL
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GRAPHEDGE
   ! 
   ! 
   !   Function: MURGE_GraphEnd
   ! 
   !   End the graph building.
   !   
   !   Needs <MURGE_GraphBegin> to be called before.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.  
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_GraphBegin> was not 
   !                           called before or there are missing edges.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GRAPHEND(ID, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GRAPHEND
   ! 
   ! 
   !   Function: MURGE_GraphGlobalCSR
   ! 
   !   Build an adjency graph from a Compress Sparse Row matrix pattern.
   !   
   !   Needs <MURGE_SetDefaultOptions> to be called before.
   ! 
   !   This function depends on integer parameter *MURGE_BASEVAL*.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     N       - Global number of columns
   !     rowptr  - Index of the first element of each row in *COLS* array.
   !     COLS    - Global column numbers array.
   !     root    - Root processor : this processor enter the global data 
   !               (-1 for all processors).
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or 
   !                           CSR is not correct.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GRAPHGLOBALCSR(ID, N, ROWPTR, COLS, ROOT, IERROR)
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ID, N, ROOT
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ROWPTR
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COLS
     INTEGER(KIND=4), DIMENSION(0), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GRAPHGLOBALCSR
   ! 
   ! 
   !   Function: MURGE_GraphGlobalCSC
   ! 
   !   Build an adjency graph from a Compress Sparse Column matrix pattern.
   !   
   !   Needs <MURGE_SetDefaultOptions> to be called before.
   ! 
   !   This function depends on integer parameter *MURGE_BASEVAL*.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     N       - Global number of columns
   !     COLPTR  - Index of the first element of each column in *ROWS* array.
   !     ROWS    - Global row number array.
   !     root    - Root processor : this processor enter the global data 
   !               (-1 for all processors).  
   ! 
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or 
   !                           CSC is not correct.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GRAPHGLOBALCSC(ID, N, COLPTR, ROWS, ROOT, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, N, ROOT
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COLPTR
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ROWS
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GRAPHGLOBALCSC
   ! 
   ! 
   !   Function: MURGE_GraphGlobalIJV
   ! 
   !   Build an adjency graph from a Compress Sparse Column matrix pattern.
   !   
   !   Needs <MURGE_SetDefaultOptions> to be called before.
   ! 
   !   This function depends on integer parameter *MURGE_BASEVAL*.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     N       - Global number of unknowns.
   !     NNZ     - Global number of non zeros.
   !     ROW     - Global column number array.
   ! edges.
   !     COL     - Global row number array.
   !     root    - Root processor : this processor enter the global data 
   !               (-1 for all processors).  
   ! 
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or 
   !                           graph IJ is not correct.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GRAPHGLOBALIJV(ID, N, NNZ, ROW, COL, ROOT, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, N, ROOT
     INTEGER(KIND=4),               INTENT(IN)  :: NNZ
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ROW
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COL
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GRAPHGLOBALIJV
   ! 
   ! 
   !    Group: IO functions
   ! 
   !    Allows to save and load solver state after preprocessing.
   ! 
   ! 
   !   Function: MURGE_Save
   !   
   !   Runs preprocessing step, if not done yet, and save the result to disk,
   !   into *directory*, so that it can be resume using <MURGE_Load>.
   ! 
   !   Needs the graph to be built.
   ! 
   !   Parameters: 
   !     id        - Solver instance identification number.  
   !     directory - Path to the directory where to save the solver step.
   ! 
   !   In Fortran, *STR_LEN* is the length of the string directory.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If graph hasn't been built.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_IO        - If file(s) couldn't be writen.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SAVE(ID, DIRECTORY, STR_LEN, IERROR)
     INTEGER(KIND=4),             INTENT(IN)  :: ID, STR_LEN
     CHARACTER(len=*),            INTENT(IN)  :: DIRECTORY
     INTEGER(KIND=4),             INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SAVE
   ! 
   ! 
   !   Function: MURGE_Load
   !   
   !   Loads preprocessing result from disk, into *directory*, 
   !   where it had been saved by <MURGE_Save>.
   ! 
   !   If preprocessing data was already computed or loaded, it will 
   !   be overwriten.
   !   
   !   Needs <MURGE_SetDefaultOptions> to be called before to initiate
   !   solver instance data.
   !   
   !   Parameters: 
   !     id        - Solver instance identification number.  
   !     directory - Path to the directory where to load the solver 
   !                 preprocessing data.
   ! 
   !   In Fortran, *STR_LEN* is the length of the string directory.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If <MURGE_SetDefaultOptions> was not 
   !                           called before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_IO        - If file(s) couldn't be read.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_LOAD(ID, DIRECTORY, STR_LEN, IERROR)
     INTEGER(KIND=4),             INTENT(IN)  :: ID, STR_LEN
     CHARACTER(len=*),            INTENT(IN)  :: DIRECTORY
     INTEGER(KIND=4),             INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_LOAD
   ! 
   ! 
   !    Group: Getting new distribution
   !    
   ! 
   ! 
   !   Function: MURGE_GetLocalNodeNbr
   !   
   !   Computes preprocessing step, if not done, and the number of
   !   Nodes in the new ditribution of the matrix.
   ! 
   !   Needs the graph to be built.
   ! 
   !   Parameters:
   !     id        - Solver instance identification number.  
   !     nodenbr   - *INTS* where to store number of nodes.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If graph hasn't been built.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
   !                           or *nodenbr* is *NULL* (can occur in C).
   !     
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETLOCALNODENBR(ID, NODENBR, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: ID
     INTEGER(KIND=4), INTENT(OUT) :: NODENBR
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETLOCALNODENBR 
   ! 
   ! 
   !   Function: MURGE_GetLocalNodeList
   !   
   !   Computes the local node list, corresponding to 
   !   the new distribution, after preprocessing.
   !   
   !   *nodelist* array has to be allocated before calling
   !   this function.
   ! 
   !   As it's result determines the size of *nodelist*
   !   array, <MURGE_GetLocalNodeNbr> should be run before it.
   !     
   !   Parameters:
   !     id        - Solver instance identification number.  
   !     nodelist  - Array where to store the list of local nodes.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - if <MURGE_GetLocalNodeNbr> has not been called 
   !                           before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
   !                           or *nodelist* is *NULL* (can occur in C).
   !     
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETLOCALNODELIST(ID, NODELIST, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID      
     ! Warning : 0 is not the size of the array.
     ! Writing DIMENSION(:) does not work with
     ! the C function call (fortran send the array size?)
     INTEGER(KIND=4), DIMENSION(0), INTENT(OUT) :: NODELIST
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETLOCALNODELIST
   ! 
   ! 
   !   Function: MURGE_GetLocalUnknownNbr
   !   
   !   Computes preprocessing step, if not done, and the number of
   !   Unkowns in the new ditribution of the matrix.
   ! 
   !   Needs the graph to be built.
   ! 
   !   Parameters:
   !     id        - Solver instance identification number.  
   !     unkownnbr     - *INTS* where to store number of unkowns.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If graph hasn't been built.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
   !                           or *unkownnbr* is *NULL* (can occur in C).
   !     
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETLOCALUNKOWNNBR(ID, UNKOWNNBR, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: ID
     INTEGER(KIND=4), INTENT(OUT) :: UNKOWNNBR
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETLOCALUNKOWNNBR 
   ! 
   ! 
   !   Function: MURGE_GetLocalUnknownList
   !   
   !   Computes the local unkown list, corresponding to 
   !   the new distribution, after preprocessing.
   !   
   !   *unkownlist* array has to be allocated before calling
   !   this function.
   ! 
   !   As it's result determines the size of *unkownlist*
   !   array, <MURGE_GetLocalUnkownNbr> should be run before it.
   !     
   !   Parameters:
   !     id        - Solver instance identification number.  
   !     unkownlist  - Array where to store the list of local unkowns.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - if <MURGE_GetLocalUnkownNbr> has not been called 
   !                           before.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
   !                           or *unkownlist* is *NULL* (can occur in C).
   !     
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETLOCALUNKOWNLIST(ID, UNKOWNLIST, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID      
     INTEGER(KIND=4), DIMENSION(0), INTENT(OUT) :: UNKOWNLIST
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETLOCALUNKOWNLIST
   ! 
   ! 
   !    Group: Filling the matrix   
   ! 
   ! 
   !   Function: MURGE_AssemblyBegin
   ! 
   !   Begin Filling up sequence for the matrix, will allocate
   !   temporary structures used to build the matrix.
   ! 
   !   It will perform preprocessing if it has not been done yet.
   ! 
   !   It needs graph to be built.
   !   
   !   Several assembly sequences can be performed on the same matrix.
   ! 
   !   *mode* shouldn't be *MURGE_ASSEMBLY_RESPECT* if neither 
   !   <MURGE_GetLocalNodeList> nor <MURGE_GetLocalUnknownList> has been called.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.  
   !     coefnbr - Number of coeficients he will had.
   !     op      - Operation to perform for coefficient which appear
   !               several tim (see <MURGE_ASSEMBLY_OP>).
   !     op2     - Operation to perform when a coefficient is set by
   !               two different processors (see <MURGE_ASSEMBLY_OP>).
   !     mode    - Indicates if user ensure he will respect solvers distribution 
   !               (see <MURGE_ASSEMBLY_MODE>).
   !     sym     - Indicates if user will give coefficient in a symmetric way 
   !               (ie: only triangullar part) or not. 
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If graph hasn't been built before.
   !     MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
   !                           *op*, *mode*, *sym*, or *coefnbr* are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_ASSEMBLYBEGIN(ID, COEFNBR, OP, OP2, MODE, SYM, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID, OP, OP2, MODE, SYM
     INTEGER(KIND=4),      INTENT(IN)  :: COEFNBR
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYBEGIN
   ! 
   ! 
   !   Function: MURGE_AssemblySetValue
   ! 
   !   Set a coefficient value in the matrix.
   !   
   !   Needs <MURGE_AssemblyBegin> to be called before.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.  
   !     ROW     - Global row number of the coefficient.
   !     COL     - Global column number of the coefficient.
   !     value   - value of the coefficient.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If we are not in an assembly section.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
   !                           *ROW* or *COL* are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_ASSEMBLYSETVALUE(ID, ROW, COL, VALUE, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID, ROW, COL
     REAL(KIND=8),         INTENT(IN)  :: VALUE
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYSETVALUE
   ! 
   ! 
   !   Function: MURGE_AssemblySetNodeValues
   ! 
   !   Set coefficients value for a node in the matrix.
   !   
   !   Needs <MURGE_AssemblyBegin> to be called before.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.  
   !     ROW     - Global row number of the coefficient.
   !     COL     - Global column number of the coefficient.
   !     values  - value of the coefficient.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If we are not in an assembly section.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
   !                           *ROW* or *COL* are not valid, or, in C, if 
   ! 			  *values* is NULL.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_ASSEMBLYSETNODEVALUES(ID, ROW, COL, VALUES, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, ROW, COL
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: VALUES
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYSETNODEVALUES
   ! 
   ! 
   !   Function: MURGE_AssemblySetBlockValues
   ! 
   !   Set coefficients value for a dens block in the matrix.
   !   
   !   Needs <MURGE_AssemblyBegin> to be called before.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.  
   !     nROW    - Number of rows in the dense matrix.
   !     ROWlist - List of global row numbers.
   !     nCOL    - Number of columns in the dense matrix.
   !     COLlist - List of global column numbers.
   !     values  - Values array, by column (Fortran style)
   !     
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If we are not in an assembly section.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
   !                           *nI* or *nJ* are not valid, or, in C, if 
   ! 			  *Ilist*, *Jlist* or *values* is NULL.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_ASSEMBLYSETBLOCKVALUES(ID, NROW, ROWLIST, &
                                   & NCOL, COLLIST, VALUES, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, NROW, NCOL
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ROWLIST
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COLLIST
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: VALUES
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYSETBLOCKVALUES
   ! 
   ! 
   !   Function: MURGE_AssemblyEnd
   ! 
   !   End Filling up sequence for the matrix.
   ! 
   !   Needs <MURGE_AssemblyBegin> to be called before.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_ORDER     - If we are not in an assembly section.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_ASSEMBLYEND(ID, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYEND
   ! 
   ! 
   !   Function: MURGE_MatrixReset
   ! 
   !   Reset the matrix structure.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_MATRIXRESET(ID, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_MATRIXRESET
   ! 
   ! 
   !   Function: MURGE_MatrixGlobalCSR
   ! 
   !   Add the given global Compress Sparse Row matrix to the matrix.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     N       - Number of columns.
   !     rowptr  - Index of the first element of each row in *COLS* and 
   !               *values* array.
   !     COLS    - Column number array.
   !     values  - values array.
   !     root    - Root processor for MPI communications.
   !     op      - Operation to perform if a coefficient appear twice
   !               (see <MURGE_ASSEMBLY_OP>).
   !     sym     - Indicates if user will give coefficient in a symmetric way 
   !               (ie: only triangullar part) or not.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, 
   !                           if *root*, *op* or the CSR are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_MATRIXGLOBALCSR(ID, N, ROWPTR, COLS, VALUES, &
                                  & ROOT, OP, SYM,  IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, N, ROOT, OP, SYM
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ROWPTR
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COLS
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: VALUES
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_MATRIXGLOBALCSR
   ! 
   ! 
   !   Function: MURGE_MatrixGlobalCSC
   ! 
   !   Add the given global Compress Sparse Column matrix to the matrix.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     N       - Number of columns.
   !     colptr  - Index of the first element of each column in *ROWS* and 
   !               *values* array.
   !     ROWS    - Row number array.
   !     values  - values array.
   !     root    - Root processor for MPI communications.
   !     op      - Operation to perform if a coefficient appear twice
   !               (see <MURGE_ASSEMBLY_OP>).
   !     sym     - Indicates if user will give coefficient in a symmetric way 
   !               (ie: only triangullar part) or not.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, 
   !                           if *root*, *op* or the CSC are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_MATRIXGLOBALCSC(ID, N, COLPTR, ROWS, &
                                & VALUES, ROOT, OP, SYM, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, N, ROOT, OP, SYM
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COLPTR
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ROWS
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: VALUES
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_MATRIXGLOBALCSC
   ! 
   ! 
   !   Function: MURGE_MatrixGlobalIJV
   ! 
   !   Add the given global Compress Sparse Column matrix to the matrix.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     N       - Number of edges.
   !     NNZ     - Number of non zeros.
   !     ROWS    - Global row number array.
   !     COLS    - Global column number array.
   !     values  - values array.
   !     root    - Root processor for MPI communications.
   !     op      - Operation to perform if a coefficient appear twice
   !               (see <MURGE_ASSEMBLY_OP>).
   !     sym     - Indicates if user will give coefficient in a symmetric way 
   !               (ie: only triangullar part) or not.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, 
   !                           if *root*, *op*, *ROWS* or *COLS* are not valid.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_MATRIXGLOBALIJV(ID, N, NNZ, ROWS, COLS, VALUES, &
                                  & ROOT, OP, SYM, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, ROOT, OP, SYM, N
     INTEGER(KIND=4),               INTENT(IN)  :: NNZ
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: ROWS, COLS
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: VALUES
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_MATRIXGLOBALIJV
   ! 
   ! 
   !    Group: Filling the right-hand-side member
   !    
   ! 
   ! 
   !   Function: MURGE_SetGlobalRHS
   ! 
   !   Set the right-hand-side member in global mode.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     b       - Array of size global column number which correspond to the 
   !               right-hand-side member.
   !     op      - Operation to perform if a coefficient appear twice
   !               (see <MURGE_ASSEMBLY_OP>).
   !     root    - Indicates which processor has the right-hand-side member,
   !               -1 for all.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SETGLOBALRHS(ID, B, OP, ROOT, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, OP, ROOT
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: B
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETGLOBALRHS
   ! 
   ! 
   !   Function: MURGE_SetLocalRHS
   ! 
   !   Set the right-hand-side member in local mode.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     b       - Array of size local column number which correspond to the 
   !               right-hand-side member.
   !     op      - Operation to perform if a coefficient appear twice
   !               (see <MURGE_ASSEMBLY_OP>).
   !     op2     - Operation to perform when a coefficient is set by
   !               two different processors (see <MURGE_ASSEMBLY_OP>).
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SETLOCALRHS(ID, B, OP, OP2, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, OP, OP2
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: B
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETLOCALRHS
   ! 
   ! 
   !   Function: MURGE_SetRHS
   ! 
   !   Set the right-hand-side member, giving the list of 
   !   coefficient that we set.
   ! 
   !   *mode* shouldn't be *MURGE_ASSEMBLY_RESPECT* if neither 
   !   <MURGE_GetLocalNodeList> nor <MURGE_GetLocalUnknownList> has been called.
   ! 
   !   Parameters: 
   !     id       - Solver instance identification number.
   !     n        - Number of coefficients to set.
   !     coefsidx - List of global index of the coefficients to set.
   !     B        - Array of coefficients values. 
   !     op       - Operation to perform if a coefficient appear twice
   !                (see <MURGE_ASSEMBLY_OP>).
   !     op2     - Operation to perform when a coefficient is set by
   !               two different processors (see <MURGE_ASSEMBLY_OP>).
   !     mode     - Indicates if user ensure he will respect solvers distribution 
   !                (see <MURGE_ASSEMBLY_MODE>).
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
   !                           if *mode* or *op* are not valid, or, in C, 
   ! 			  if *coeflist* or *b* are NULL.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_SETRHS(ID, N, COEFSIDX, B, OP, OP2, MODE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, N, OP, OP2, MODE
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COEFSIDX
     REAL(KIND=8),    DIMENSION(0), INTENT(IN)  :: B
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETRHS
   ! 
   ! 
   !   Function: MURGE_RHSReset
   ! 
   !   Reset the right-hand-side.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_RHSRESET(ID, IERROR)
     INTEGER(KIND=4),      INTENT(IN)  :: ID
     INTEGER(KIND=4),      INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_RHSRESET
   ! 
   ! 
   !    Group: Getting the solution
   !    
   ! 
   ! 
   !   Function: MURGE_GetGlobalSolution
   ! 
   !   Perform Factorization and Solve, if needed, 
   !   and then fill the global solution in *x*.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     x       - Array of size global column number*dof which will contain
   !               the solution
   !     root    - Indicates which processor will have the solution 
   !               at the end of the call, -1 for all.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETGLOBALSOLUTION(ID, X, ROOT, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, ROOT
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: X
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETGLOBALSOLUTION
   ! 
   ! 
   !   Function: MURGE_GetLocalSolution
   ! 
   !   Perform Factorization and Solve, if needed, 
   !   and then fill the local solution in *x*.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     x       - Array of size local column number*dof which will contain
   !               the solution
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETLOCALSOLUTION(ID, X, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: X
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETLOCALSOLUTION
   ! 
   ! 
   !   Function: MURGE_GetSolution
   ! 
   !   Perform Factorization and Solve, if needed, 
   !   and then fill the solution in *x* followin the given
   !   index list.
   ! 
   !   Parameters: 
   !     id       - Solver instance identification number.
   !     n        - Number of coefficients user wants to get.
   !     coefsidx - List of the coefficients user wants to get.
   !     x        - Array of size dof*n which will contain
   !                the solution.
   !     mode     - Indicates if the user is sure to respect the distribution.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETSOLUTION(ID, N, COEFSIDX, X, MODE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, MODE, N
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COEFSIDX
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: X
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETSOLUTION
   ! 
   ! 
   !   Group: Scaling
   ! 
   ! 
   !   Function: MURGE_GetGlobalNorm
   ! 
   !   Compute the global norm array following a norm rule.
   ! 
   !   Must be performed after assembly step.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     norm    - Array of size global column number*dof which will contain
   !               the norm values
   !     root    - Indicates which processor will have the norm array
   !               at the end of the call, -1 for all.
   !     rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_ORDER     - If the assembly has not been performed.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETGLOBALNORM(ID, NORM, ROOT, RULE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, ROOT, RULE
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: NORM
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETGLOBALNORM
   ! 
   ! 
   !   Function: MURGE_GetLocalNorm
   ! 
   !   Compute the local norm array following a norm rule.
   !   
   !   Must be performed after assembly step.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     norm    - Array of size local column number*dof which will contain
   !               the solution
   !     rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_ORDER     - If the assembly has not been performed.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETLOCALNORM(ID, NORM, RULE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, RULE
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: NORM
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETLOCALNORM
   ! 
   ! 
   !   Function: MURGE_GetNorm
   ! 
   !   Compute the indicated part of the norm array 
   !   following a norm rule.
   !   
   !   Must be performed after assembly step.
   !  
   ! 
   !   Parameters: 
   !     id       - Solver instance identification number.
   !     n        - Number of coefficients user wants to get norm of.
   !     coefsidx - List of the coefficients user wants to get norm of.
   !     norm     - Array of size dof*n which will contain
   !                the solution.
   !     rule     - Rule to follow to build norm array, see <MURGE_NORM_RULES>
   !     mode     - Indicates if the user is sure to respect the distribution.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_ORDER     - If the assembly has not been performed.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETNORM(ID, N, COEFSIDX, NORM, RULE, MODE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, MODE, N, RULE
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COEFSIDX
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: NORM
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETNORM
   ! 
   ! 
   !    Function: MURGE_ApplyGlobalScaling
   !    
   !    Apply scaling to local unknowns.
   ! 
   !    Must be performed after assembly step.
   ! 
   !    Parameters:
   !      id      - Solver instance identification number.
   !      scal    - Scaling user wants to apply.
   !      root    - Indicates which processor that posses the scaling array, 
   !                -1 for all.
   !      sc_mode - Indicate if the scaling is applied on rows or on columns.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_ORDER     - If the assembly has not been performed.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_APPLYGLOBALSCALING(ID, SCAL, ROOT, SC_MODE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, ROOT, SC_MODE
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: SCAL
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_APPLYGLOBALSCALING
   ! 
   ! 
   ! 
   !   Function: MURGE_ApplyLocalScaling
   ! 
   !   Apply the local scaling array on the matrix.
   !   
   !   Must be performed after assembly step.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   !     scal    - Array of size local column number*dof which will contain
   !               the solution.
   !     sc_mode - Indicate if the scaling is applied on rows or on columns.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_ORDER     - If the assembly has not been performed.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_APPLYLOCALSCALING(ID, SCAL, SC_MODE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, SC_MODE
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: SCAL
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_APPLYLOCALSCALING
   ! 
   ! 
   !   Function: MURGE_ApplyScaling
   ! 
   !   Apply the scaling array on the indicated part of the matrix
   !   
   !   Must be performed after assembly step.
   !  
   ! 
   !   Parameters: 
   !     id       - Solver instance identification number.
   !     n        - Number of coefficients user wants to scale.
   !     coefsidx - List of the coefficients user wants to scale.
   !     scal     - Array of size dof*n which will contain
   !                the solution.
   !     sc_mode  - Indicate if the scaling is applied on rows or on columns.
   !     mode     - Indicates if the user is sure to respect the distribution.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   !     MURGE_ERR_ORDER     - If the assembly has not been performed.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_APPLYSCALING(ID, N, COEFSIDX, SCAL, SC_MODE, MODE, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, SC_MODE, MODE, N
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: COEFSIDX
     REAL(KIND=8),    DIMENSION(0), INTENT(OUT) :: SCAL
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_APPLYSCALING
   ! 
   ! 
   !    Group: Cleaning up this mess
   !    
   ! 
   ! 
   !   Function: MURGE_Clean
   ! 
   !   Clean the given instance of the solver structure's.
   ! 
   !   Parameters: 
   !     id      - Solver instance identification number.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_CLEAN(ID, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: ID
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_CLEAN
   ! 
   ! 
   !   Function: MURGE_Finalize
   !   
   !   Clean all not cleaned instances and instances ID array.
   !   
   !   Returns: 
   !     MURGE_SUCCESS       - If function runned successfully.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_FINALIZE(IERROR)
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_FINALIZE
   ! 
   ! 
   !   Function: MURGE_GetSolver
   !    
   !   Return the solver ID Murge was compiled with.
   !    
   !   Parameters:
   !     solver - Integer to store solver ID.
   !    
   !   Returns:
   !     MURGE_SUCCESS - If execution ended normaly.
   !    
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_GETSOLVER(SOLVER, IERROR)
     INTEGER(KIND=4),  INTENT(OUT) :: SOLVER
     INTEGER(KIND=4),  INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETSOLVER   
   ! 
   ! 
   !   Group: Getting Murge's metrics
   !  
   ! 
   !   Function: MURGE_GetMetricINT
   ! 
   !   Get an integer metric from MURGE.
   ! 
   !   See <MURGE_IINFOS> and the solver documentation to 
   !   get available metric list.
   !   
   !   Parameters:
   !     id     - Solver instance identification number.
   !     metric - Wanted metric identification number.
   !     value  - Integer which will contain the value of the metric.
   ! 
   !   Returns:
   !     MURGE_SUCCESS - If execution ended normaly.
   !     MURGE_ERR_ORDER     - If metric is not available in the current 
   !                           solver state.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
   !                           *metric* or *value* are not valid.
   ! 
   !     
   !   Fortran interface:
   SUBROUTINE MURGE_GETINFOINT(ID, INFO, VALUE, IERROR)
    INTEGER(KIND=4), INTENT(IN)  :: ID, INFO
    INTEGER(KIND=4), INTENT(OUT) :: VALUE
    INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETINFOINT
   !  
   ! 
   !   Function: MURGE_GetInfoREAL
   ! 
   !   Get a real metric value from MURGE.
   ! 
   !   See <MURGE_RINFOS> and the solver documentation to 
   !   get available metric list.
   !   
   !   Parameters:
   !     id     - Solver instance identification number.
   !     metric - Wanted metric identification number.
   !     value  - Real which will contain the value of the metric.
   ! 
   !   Returns:
   !     MURGE_SUCCESS - If execution ended normaly.
   !     MURGE_ERR_ORDER     - If metric is not available in the current 
   !                           solver state.
   !     MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
   !                           *metric* or *value* are not valid.
   ! 
   !     
   !   Fortran interface:
   SUBROUTINE MURGE_GETINFOREAL(ID, INFO, VALUE, IERROR)
    INTEGER(KIND=4), INTENT(IN)  :: ID, INFO
    REAL(KIND=8),    INTENT(OUT) :: VALUE
    INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_GETINFOREAL
   !  
   ! 
   !   Function: MURGE_PrintError
   !   
   !   Print the error message corresponding to ierror
   !   Parameters: 
   !     error_number  - Error identification number.
   ! 
   !   Returns: 
   !     MURGE_ERR_PARAMETER - If ierror does not match an error number
   !     MURGE_SUCCESS       - If function runned successfully.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_PRINTERROR(ERROR_NUMBER, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: IERROR
     INTEGER(KIND=4), INTENT(OUT) :: ERROR_NUMBER
   END SUBROUTINE MURGE_PRINTERROR
   ! 
   ! 
   !   Function: MURGE_ExitOnError
   !   
   !   Print the error message corresponding to ierror.
   !   If the ierr is not MURGE_SUCCESS then the program is stopped.
   ! 
   !   Parameters: 
   !     ierror         - Error identification number.
   ! 
   !   Returns: 
   !     MURGE_SUCCESS   - If function runned successfully; stop the program otherwise.
   ! 
   !   Fortran interface:
   !   >
   SUBROUTINE MURGE_EXITONERROR(ERROR_NUMBER, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: IERROR
     INTEGER(KIND=4), INTENT(OUT) :: ERROR_NUMBER
   END SUBROUTINE MURGE_EXITONERROR
   ! 
   ! 
!******************************************************************************
!* File: murge_pastix.inc                                                     *
!*                                                                            *
!* Fortran function declaration for Murge function used only in PaStiX.       *
!*                                                                            *
!* Authors:                                                                   *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                 *
!*                                                                            *
!******************************************************************************


!******************************************************************************
!* Function: MURGE_ANALYZE                                                    *
!*                                                                            *
!* Perform matrix analyse:                                                    *
!*   - Compute a new ordering of the unknows                                  *
!*   - Compute the symbolic factorisation of the matrix                       *
!*   - Distribute column blocks and computation on processors                 *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned succesfuly.                     *
!*   MURGE_ERR_ORDER     - If function the graph is not built.                *
!*   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
!*                                                                            *
!******************************************************************************

   SUBROUTINE MURGE_ANALYZE(ID, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: ID
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ANALYZE

!******************************************************************************
!* Function: MURGE_FACTORIZE                                                  *
!*                                                                            *
!* Perform matrix factorisation.                                              *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned succesfuly.                     *
!*   MURGE_ERR_ORDER     - If function the graph is not built.                *
!*   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
!*                                                                            *
!******************************************************************************

   SUBROUTINE MURGE_FACTORIZE(ID, IERROR)
     INTEGER(KIND=4), INTENT(IN)  :: ID
     INTEGER(KIND=4), INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_FACTORIZE


!******************************************************************************
!* Function: MURGE_SETORDERING                                                *
!*                                                                            *
!* Set a personal ordering to perform factorization.                          *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*   permutation - Permutation to perform factorization.                      *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned succesfuly.                     *
!*   MURGE_ERR_ORDER     - If graph has not been built.                       *
!*   MURGE_ERR_PARAMETER - If permutation is not a valid memory area.         *
!*                                                                            *
!******************************************************************************

  SUBROUTINE MURGE_SETORDERING(ID, PERMUTATION, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: PERMUTATION
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_SETORDERING

!******************************************************************************
!* Function: MURGE_GETLOCALPRODUCT                                            *
!*                                                                            *
!* Perform the product A * X.                                                 *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*   x  - Array in which the local part of the product will be stored.        *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************

  SUBROUTINE MURGE_GETLOCALPRODUCT(ID, X, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     REAL(KIND=8), DIMENSION(0), INTENT(OUT) :: X
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_GETLOCALPRODUCT

!******************************************************************************
!* Function: MURGE_GETGLOBALPRODUCT                                           *
!*                                                                            *
!* Perform the product A * X.                                                 *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   x    - Array in which the product will be stored.                        *
!*   root - Rank of the process which will own the product at end of call,    *
!*          use -1 for all processes.                                         *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************
  SUBROUTINE MURGE_GETGLOBALPRODUCT(ID, X, root, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     REAL(KIND=8), DIMENSION(0), INTENT(OUT) :: X
     INTEGER(KIND=4),               INTENT(IN)  :: root
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_GETGLOBALPRODUCT

!******************************************************************************
!* Function: MURGE_PRODUCTSETLOCALNODENBR                                     *
!*                                                                            *
!*                                                                            *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   n    - Number of local nodes.                                            *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************
  SUBROUTINE MURGE_PRODUCTSETLOCALNODENBR(ID, n, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     INTEGER(KIND=4),               INTENT(IN)  :: n
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_PRODUCTSETLOCALNODENBR

!******************************************************************************
!* Function: MURGE_PRODUCTSETGLOBALNODENBR                                    *
!*                                                                            *
!*                                                                            *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   N    - Number of global nodes.                                           *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************
  SUBROUTINE MURGE_PRODUCTSETGLOBALNODENBR(ID, N, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     INTEGER(KIND=4),               INTENT(IN)  :: N
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_PRODUCTSETGLOBALNODENBR

!******************************************************************************
!* Function: MURGE_PRODUCTSETLOCALNODELIST                                    *
!*                                                                            *
!*                                                                            *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   l2g  - List of local nodes.                                              *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************

  SUBROUTINE MURGE_PRODUCTSETLOCALNODELIST(ID, l2g, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     INTEGER(KIND=4), DIMENSION(0), INTENT(IN)  :: l2g
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_PRODUCTSETLOCALNODELIST

!******************************************************************************
!* Function: MURGE_FORCENOFACTO                                               *
!*                                                                            *
!* Prevent Murge from running factorisation even if matrix has changed.       *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!* Returns:                                                                   *
!*   MURGE_SUCCESS                                                            *
!*                                                                            *
!******************************************************************************!
   SUBROUTINE MURGE_FORCENOFACTO(ID, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_FORCENOFACTO

   SUBROUTINE MURGE_SETLOCALNODELIST(ID, N, LIST, IERROR)
     INTEGER(KIND=4),               INTENT(IN)  :: ID, N
     INTEGER(KIND=4), DIMENSION(N), INTENT(IN)  :: LIST
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETLOCALNODELIST

!******************************************************************************
!* Function: MURGE_AssemblySetSequence                                        *
!*                                                                            *
!* Create a sequence of entries to build a matrix and store it for being      *
!* reused.                                                                    *
!*                                                                            *
!* Parameters:                                                                *
!*   id      - Solver instance identification number.                         *
!*   coefnbr - Number of entries.                                             *
!*   ROWs    - List of rows in the sequence.                                  *
!*   COLs    - List of columns in the sequence.                               *
!*   op      - Operation to perform for coefficient which appear              *
!*             several tim (see <MURGE_ASSEMBLY_OP>).                         *
!*   op2     - Operation to perform when a coefficient is set by              *
!*             two different processors (see <MURGE_ASSEMBLY_OP>).            *
!*   mode    - Indicates if user ensure he will respect solvers distribution  *
!*             (see <MURGE_ASSEMBLY_MODE>).                                   *
!*   nodes   - 0 entries are entered value by value,                          *
!*             1 entries are entries node by node.                            *
!*   id_seq  - Sequence ID.                                                   *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned successfully.                   *
!*   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
!*   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
!*   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
!*                         *op*, *mode*, *sym*, or *coefnbr* are not valid.   *
!******************************************************************************/
   SUBROUTINE MURGE_ASSEMBLYSETSEQUENCE (ID, NNZ, ROWs, COLs, OP, OP2, MODE, NODES, ID_SEQ, IERROR)
     INTEGER(KIND=4),                     INTENT(IN)  :: ID, OP, OP2, MODE, NODES
     INTEGER(KIND=4),                     INTENT(IN)  :: NNZ
     INTEGER(KIND=4), DIMENSION(NNZ),     INTENT(IN)  :: ROWs, COLs
     INTEGER(KIND=4),                     INTENT(OUT) :: ID_SEQ, IERROR
   END SUBROUTINE MURGE_ASSEMBLYSETSEQUENCE
          

!******************************************************************************
!* MURGE_AssemblySetSequence                                                  *
!*                                                                            *
!* Assembly the matrix using a stored sequence.                               *
!*                                                                            *
!* Parameters:                                                                *
!*   id      - Solver instance identification number.                         *
!*   id_seq  - Sequence ID.                                                   *
!*   values  - Values to insert in the CSC.                                   *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned successfully.                   *
!*   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
!*   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
!*   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
!*                         *id_seq* or *values* are not valid.                *
!******************************************************************************/
   SUBROUTINE MURGE_ASSEMBLYUSESEQUENCE(ID, ID_SEQ, VALUES, IERROR);
     INTEGER(KIND=4),               INTENT(IN)  :: ID, ID_SEQ 
     REAL(KIND=8), DIMENSION(0), INTENT(IN)  :: VALUES
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYUSESEQUENCE
      
!******************************************************************************
!* Function: MURGE_AssemblyDeleteSequence                                     *
!*                                                                            *
!* Destroy an assembly sequence                                               *
!*                                                                            *
!*   id      - Solver instance identification number.                         *
!*   id_seq  - Sequence ID.                                                   *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned successfully.                   *
!*   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
!*   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
!*   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
!*                         *id_seq* is not valid.                             *
!******************************************************************************/
   SUBROUTINE MURGE_ASSEMBLYDELETESEQUENCE(ID, ID_SEQ, IERROR);
     INTEGER(KIND=4),               INTENT(IN)  :: ID, ID_SEQ 
     INTEGER(KIND=4),               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYDELETESEQUENCE
END INTERFACE

!
!   Enum: KINDS
!
!   Solvers type kinds
!
!   Contains:
!     MURGE_COEF_KIND - Kind to use for COEF
!     MURGE_REAL_KIND - Kind to use for REAL
!     MURGE_INTS_KIND - Kind to use for INTS
!     MURGE_INTL_KIND - Kind to use for INTL
!
INTEGER, PARAMETER :: MURGE_COEF_KIND                = 8
INTEGER, PARAMETER :: MURGE_REAL_KIND                = 8
INTEGER, PARAMETER :: MURGE_INTS_KIND                = 4
INTEGER, PARAMETER :: MURGE_INTL_KIND                = 4
!   Group: Murge's constants
!  
! 
!   Enum: MURGE_RETURNS
!   
!   Murge error return values.
!   
!   Contains:
!     MURGE_SUCCESS             - If function runs correctly.
!     MURGE_ERR_ALLOCATE        - If couldn't allocate.
!     MURGE_ERR_IO              - If there was an input or output error.
!     MURGE_ERR_PARAMETER       - If one parameter is not correct.
!     MURGE_ERR_ORDER           - If function were run in wrong order.
!     MURGE_ERR_SOLVER          - Internal solver error.
!     MURGE_ERR_NOT_IMPLEMENTED - Not yet implemented.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_SUCCESS                  = 1
INTEGER(KIND=4), PARAMETER :: MURGE_ERR_ALLOCATE             = 2
INTEGER(KIND=4), PARAMETER :: MURGE_ERR_IO                   = 3
INTEGER(KIND=4), PARAMETER :: MURGE_ERR_PARAMETER            = 4
INTEGER(KIND=4), PARAMETER :: MURGE_ERR_ORDER                = 5
INTEGER(KIND=4), PARAMETER :: MURGE_ERR_SOLVER               = 6
INTEGER(KIND=4), PARAMETER :: MURGE_ERR_NOT_IMPLEMENTED      = 7
! 
!   Enum: MURGE_IPARAM
! 
!   Murge integer parameters identifiers.
!   
!   Solvers may implement is own list of parameters.
! 
!   MURGE_IPARAM_BASEVAL - Numbering style , 0 for C, 1 for fortran. 
!   MURGE_IPARAM_DOF     - Number of degrees of freedom.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_IPARAM_BASEVAL           = 1024
INTEGER(KIND=4), PARAMETER :: MURGE_IPARAM_DOF               = 1025
INTEGER(KIND=4), PARAMETER :: MURGE_IPARAM_SYM               = 1026
! 
!   Enum: MURGE_RPARAM
! 
!   Murge real parameters identifiers.
!   
!   Solvers may implement is own list of parameters.
! 
!   Contains: 
!     MURGE_RPARAM_EPSILON_ERROR - Wanted norm error at the end of solve.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_RPARAM_EPSILON_ERROR     = 1024
! 
!   Enum: MURGE_IINFO
! 
!   Murge integer metrics identifiers.
!   
!   Solvers may implement is own list of parameters.
!   
!   Contains:
!     MURGE_IINFOS_NNZ - Number of non zeros in factorized matrix.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_IINFO_NNZ                = 1024
! 
!   Enum: MURGE_RINFO
! 
!   Murge real metrics identifiers.
!   
!   Solvers may implement is own list of parameters.
! 
!   Contains: 
!     MURGE_RINFO_FACT_TIME  - Factorization time.
!     MURGE_RINFO_SOLVE_TIME - Solving time.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_RINFO_FACT_TIME          = 1024
INTEGER(KIND=4), PARAMETER :: MURGE_RINFO_SOLVE_TIME         = 1025
! 
!   Enum: MURGE_ASSEMBLY_MODE
! 
!     Indicates if user can ensure that the information he is giving respects
!     the solver distribution.
!   
!     MURGE_ASSEMBLY_RESPECT - User ensure he respects distribution 
!                              during assembly.
! 			     See solver documentation.
!     MURGE_ASSEMBLY_FOOL    - User is not sure he will respect ditribution
!                              during assembly
! 
INTEGER(KIND=4), PARAMETER :: MURGE_ASSEMBLY_RESPECT         = 0
INTEGER(KIND=4), PARAMETER :: MURGE_ASSEMBLY_FOOL            = 1
! 
!   Enum: MURGE_ASSEMBLY_OP
!    
!     Operations possible when a coefficient appear twice.
! 
!     MURGE_ASSEMBLY_ADD - Coefficients will be added during assembly.
!     MURGE_ASSEMBLY_OVW - Coefficients will be overwriten during assembly.
!     MURGE_ASSEMBLY_MAX - Maximum value will be used for assembly.
!     MURGE_ASSEMBLY_MIN - Minimum value will be used for assembly.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_ASSEMBLY_ADD             = 0
INTEGER(KIND=4), PARAMETER :: MURGE_ASSEMBLY_OVW             = 1
INTEGER(KIND=4), PARAMETER :: MURGE_ASSEMBLY_MAX             = 2
INTEGER(KIND=4), PARAMETER :: MURGE_ASSEMBLY_MIN             = 3
! 
!   Enum: MURGE_SOLVER
!   
!   Solver ID for murge compliant solvers.
!   
!   Contains: 
!     MURGE_SOLVER_HIPS   - HIPS hybrid solver.
!     MURGE_SOLVER_PASTIX - PaStiX direct solver.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_SOLVER_HIPS              = 0
INTEGER(KIND=4), PARAMETER :: MURGE_SOLVER_PASTIX            = 1
! 
!   Enum: MURGE_BOOLEAN
! 
!   Boolean for murge parameters
! 
!   Contains: 
!     MURGE_BOOLEAN_FALSE - False value  
!     MURGE_BOOLEAN_TRUE  - True value
! 
INTEGER(KIND=4), PARAMETER :: MURGE_BOOLEAN_FALSE            = 0
INTEGER(KIND=4), PARAMETER :: MURGE_BOOLEAN_TRUE             = 1
! 
!    Enum: MURGE_NORM_RULES
!    
!    Flags for Murge's norm rules
! 
!    Contains:
!      MURGE_NORM_MAX_COL  - Get maximum column value  (absolute value).
!      MURGE_NORM_MAX_ROW  - Get maximum row value     (absolute value).
!      MURGE_NORM_2_COL    - Get the norm 2 of columns.
!      MURGE_NORM_2_ROW    - Get the norm 2 of rows.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_NORM_MAX_COL             = 0
INTEGER(KIND=4), PARAMETER :: MURGE_NORM_MAX_ROW             = 1
INTEGER(KIND=4), PARAMETER :: MURGE_NORM_2_COL               = 2
INTEGER(KIND=4), PARAMETER :: MURGE_NORM_2_ROW               = 3
! 
!    Enum: MURGE_SCAL_MODES
!    
!    Flags for Murge's scaling rules
! 
!    Contains:
!      MURGE_SCAL_COL  - Perform scaling on columns
!      MURGE_SCAL_ROW  - Perform scaling on rows.
! 
INTEGER(KIND=4), PARAMETER :: MURGE_SCAL_COL                 = 0
INTEGER(KIND=4), PARAMETER :: MURGE_SCAL_ROW                 = 1
!  Copyright 2008 BORDEAUX I UNIVERSITY & INRIA ! **
! ** This file is part of the PaStiX parallel sparse matrix package.
! **
! ** This software is governed by the CeCILL-C license under French law
! ** and abiding by the rules of distribution of free software. You can
! ** use, modify and/or redistribute the software under the terms of the
! ** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
! ** URL: "http://www.cecill.info".
! ** 
! ** As a counterpart to the access to the source code and rights to copy,
! ** modify and redistribute granted by the license, users are provided
! ** only with a limited warranty and the software's author, the holder of
! ** the economic rights, and the successive licensors have only limited
! ** liability.
! ** 
! ** In this respect, the user's attention is drawn to the risks associated
! ** with loading, using, modifying and/or developing or reproducing the
! ** software by the user in light of its specific status of free software,
! ** that may mean that it is complicated to manipulate, and that also
! ** therefore means that it is reserved for developers and experienced
! ** professionals having in-depth computer knowledge. Users are therefore
! ** encouraged to load and test the software's suitability as regards
! ** their requirements in conditions enabling the security of their
! ** systems and/or data to be ensured and, more generally, to use and
! ** operate it in the same conditions as regards security.
! ** 
! ** The fact that you are presently reading this means that you have had
! ** knowledge of the CeCILL-C license and that you accept its terms.
! !   header: api.h
! 
!   Header file containing constants used in PaStiX and provided to users.
! 
!   Authors:
!     Mathieu Faverge - faverge@labri.fr
!     Xavier   Lacoste - lacoste@labri.fr
!     Pierre Ramet     - ramet@labri.fr
! 
! !    enum: IPARM_ACCESS
! 
!    Integer parameters tabular accessors
! 
!    IPARM_MODIFY_PARAMETER      - Indicate if parameters have been set by user             Default: API_YES             IN
!    IPARM_START_TASK            - Indicate the first step to execute (see PaStiX steps)    Default: API_TASK_ORDERING   IN
!    IPARM_END_TASK              - Indicate the last step to execute (see PaStiX steps)     Default: API_TASK_CLEAN      IN
!    IPARM_VERBOSE               - Verbose mode (see Verbose modes)                         Default: API_VERBOSE_NO      IN
!    IPARM_DOF_NBR               - Degree of freedom per node                               Default: 1                   IN
!    IPARM_ITERMAX               - Maximum iteration number for refinement                  Default: 250                 IN
!    IPARM_MATRIX_VERIFICATION   - Check the input matrix                                   Default: API_NO              IN
!    IPARM_MC64                  - MC64 operation <pastix.h> IGNORE                         Default: 0                   IN
!    IPARM_ONLY_RAFF             - Refinement only                                          Default: API_NO              IN
!    IPARM_CSCD_CORRECT          - Indicate if the cscd has been redistributed after blend  Default: API_NO              IN
!    IPARM_NBITER                - Number of iterations performed in refinement       Default: -                   OUT
!    IPARM_TRACEFMT              - Trace format (see Trace modes)                           Default: API_TRACE_PICL      IN
!    IPARM_GRAPHDIST             - Specify if the given graph is distributed or not         Default: API_YES             IN
!    IPARM_AMALGAMATION_LEVEL    - Amalgamation level                                       Default: 5                   IN
!    IPARM_ORDERING              - Choose ordering                                          Default: API_ORDER_SCOTCH    IN
!    IPARM_DEFAULT_ORDERING      - Use default ordering parameters with \scotch{} or \metis{} Default: API_YES             IN
!    IPARM_ORDERING_SWITCH_LEVEL - Ordering switch level    (see \scotch{} User's Guide)    Default: 120                 IN
!    IPARM_ORDERING_CMIN         - Ordering cmin parameter  (see \scotch{} User's Guide)    Default: 0                   IN
!    IPARM_ORDERING_CMAX         - Ordering cmax parameter  (see \scotch{} User's Guide)    Default: 100000              IN
!    IPARM_ORDERING_FRAT         - Ordering frat parameter  (see \scotch{} User's Guide)    Default: 8                   IN
!    IPARM_STATIC_PIVOTING       - Static pivoting                                          Default: -                   OUT
!    IPARM_METIS_PFACTOR         - \metis{} pfactor                                         Default: 0                   IN
!    IPARM_NNZEROS               - Number of nonzero entries in the factorized matrix       Default: -                   OUT
!    IPARM_ALLOCATED_TERMS       - Maximum memory allocated for matrix terms                Default: -                   OUT
!    IPARM_BASEVAL               - Baseval used for the matrix                              Default: 0                   IN
!    IPARM_MIN_BLOCKSIZE         - Minimum block size                                       Default: 60                  IN
!    IPARM_MAX_BLOCKSIZE         - Maximum block size                                       Default: 120                 IN
!    IPARM_SCHUR                 - Schur mode                                               Default: API_NO              IN
!    IPARM_ISOLATE_ZEROS         - Isolate null diagonal terms at the end of the matrix     Default: API_NO              IN
!    IPARM_RHSD_CHECK            - Set to API_NO to avoid RHS redistribution                Default: API_YES             IN
!    IPARM_FACTORIZATION         - Factorization mode (see Factorization modes)             Default: API_FACT_LDLT       IN
!    IPARM_NNZEROS_BLOCK_LOCAL   - Number of nonzero entries in the local block factorized matrix Default: -                   OUT
!    IPARM_CPU_BY_NODE           - Number of CPUs per SMP node                              Default: 0                   IN
!    IPARM_BINDTHRD              - Thread binding mode (see Thread binding modes)           Default: API_BIND_AUTO       IN
!    IPARM_THREAD_NBR            - Number of threads per MPI process                        Default: 1                   IN
!    IPARM_DISTRIBUTION_LEVEL    - Distribution level IGNORE                                Default:                     IN
!    IPARM_LEVEL_OF_FILL         - Level of fill for incomplete factorization               Default: 1                   IN
!    IPARM_IO_STRATEGY           - IO strategy (see Checkpoints modes)                      Default: API_IO_NO           IN
!    IPARM_RHS_MAKING            - Right-hand-side making (see Right-hand-side modes)      Default: API_RHS_B           IN
!    IPARM_REFINEMENT            - Refinement type (see Refinement modes)                   Default: API_RAF_GMRES       IN
!    IPARM_SYM                   - Symmetric matrix mode (see Symmetric modes)              Default: API_SYM_YES         IN
!    IPARM_INCOMPLETE            - Incomplete factorization                                 Default: API_NO              IN
!    IPARM_ABS                   - ABS level (Automatic Blocksize Splitting)                Default: 1                   IN
!    IPARM_ESP                   - ESP (Enhanced Sparse Parallelism)                        Default: API_NO              IN
!    IPARM_GMRES_IM              - GMRES restart parameter                                  Default: 25                  IN
!    IPARM_FREE_CSCUSER          - Free user CSC                                            Default: API_CSC_PRESERVE    IN
!    IPARM_FREE_CSCPASTIX        - Free internal CSC (Use only without call to Refin. step) Default: API_CSC_PRESERVE    IN
!    IPARM_OOC_LIMIT             - Out of core memory limit (Mo)                            Default: 2000                IN
!    IPARM_OOC_THREAD            - Out of core thread number IGNORE                         Default: 1                   IN
!    IPARM_OOC_ID                - Out of core run ID        IGNORE                         Default: -                   OUT
!    IPARM_NB_SMP_NODE_USED      - Number of SMP node used   IGNORE                         Default:                     IN
!    IPARM_THREAD_COMM_MODE      - Threaded communication mode (see Communication modes)    Default: API_THREAD_MULT     IN
!    IPARM_NB_THREAD_COMM        - Number of thread(s) for communication                    Default: 1                   IN
!    IPARM_FILL_MATRIX           - Initialize matrix coefficients (for test only)  IGNORE   Default:                     IN
!    IPARM_INERTIA               - Return the inertia (symmetric matrix without pivoting)   Default: -                   OUT
!    IPARM_ESP_NBTASKS           - Return the number of tasks generated by ESP              Default: -                   OUT
!    IPARM_ESP_THRESHOLD         - Minimal block sizee to switch in ESP mode (128 * 128)    Default: 16384               IN
!    IPARM_DOF_COST              - Degree of freedom for cost computation (If different from IPARM_DOF_NBR) Default: 0                    IN
!    IPARM_MURGE_REFINEMENT      - Enable refinement in MURGE                               Default: API_YES             IN
!    IPARM_STARPU                - Use StarPU runtime                                       Default: API_NO              IN
!    IPARM_AUTOSPLIT_COMM        - Automaticaly split communicator to have one MPI task by node             Default: API_NO               IN
!    IPARM_FLOAT                 - Indicate the floating point type  IGNORE                 Default: -                   INOUT
!    IPARM_PID                   - Pid of the first process (used for naming the log directory) Default: -1                  OUT
!    IPARM_ERROR_NUMBER          - Return value                                             Default: -                   OUT
!    IPARM_CUDA_NBR              - Number of cuda devices                                   Default: 0                   IN
!    IPARM_TRANSPOSE_SOLVE       - Use transposed matrix during solve                       Default: API_NO              IN
!    IPARM_SIZE                  - Iparm Size                IGNORE                         Default:                     IN
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MODIFY_PARAMETER         = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_START_TASK               = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_END_TASK                 = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_VERBOSE                  = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DOF_NBR                  = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ITERMAX                  = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MATRIX_VERIFICATION      = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MC64                     = 7
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ONLY_RAFF                = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_CSCD_CORRECT             = 9
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NBITER                   = 10
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_TRACEFMT                 = 11
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_GRAPHDIST                = 12
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_AMALGAMATION_LEVEL       = 13
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING                 = 14
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DEFAULT_ORDERING         = 15
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_SWITCH_LEVEL    = 16
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_CMIN            = 17
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_CMAX            = 18
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ORDERING_FRAT            = 19
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_STATIC_PIVOTING          = 20
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_METIS_PFACTOR            = 21
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NNZEROS                  = 22
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ALLOCATED_TERMS          = 23
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_BASEVAL                  = 24
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MIN_BLOCKSIZE            = 25
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MAX_BLOCKSIZE            = 26
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_SCHUR                    = 27
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ISOLATE_ZEROS            = 28
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_RHSD_CHECK               = 29
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FACTORIZATION            = 30
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NNZEROS_BLOCK_LOCAL      = 31
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_CPU_BY_NODE              = 32
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_BINDTHRD                 = 33
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_THREAD_NBR               = 34
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DISTRIBUTION_LEVEL       = 35
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_LEVEL_OF_FILL            = 36
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_IO_STRATEGY              = 37
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_RHS_MAKING               = 38
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_REFINEMENT               = 39
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_SYM                      = 40
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_INCOMPLETE               = 41
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ABS                      = 42
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ESP                      = 43
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_GMRES_IM                 = 44
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FREE_CSCUSER             = 45
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FREE_CSCPASTIX           = 46
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_OOC_LIMIT                = 47
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_OOC_THREAD               = 48
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_OOC_ID                   = 49
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NB_SMP_NODE_USED         = 50
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_THREAD_COMM_MODE         = 51
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_NB_THREAD_COMM           = 52
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FILL_MATRIX              = 53
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_INERTIA                  = 54
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ESP_NBTASKS              = 55
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ESP_THRESHOLD            = 56
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_DOF_COST                 = 57
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_MURGE_REFINEMENT         = 58
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_STARPU                   = 59
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_AUTOSPLIT_COMM           = 60
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_FLOAT                    = 61
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_PID                      = 62
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_ERROR_NUMBER             = 63
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_CUDA_NBR                 = 64
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_TRANSPOSE_SOLVE          = 65
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IPARM_SIZE                     = 128
! !    Enum: DPARM_ACCESS
! 
!    Floating point parameters tabular accossors
! 
!    DPARM_FILL_IN            - Fill-in                                           Default: -                OUT
!    DPARM_MEM_MAX            - Maximum memory (-DMEMORY_USAGE)                   Default: -                OUT
!    DPARM_EPSILON_REFINEMENT - Epsilon for refinement                            Default: 1e^{-12}         IN
!    DPARM_RELATIVE_ERROR     - Relative backward error                           Default: -                OUT
!    DPARM_EPSILON_MAGN_CTRL  - Epsilon for magnitude control                     Default: 1e^{-31}         IN
!    DPARM_ANALYZE_TIME       - Time for Analyse step (wallclock)                 Default: -                OUT
!    DPARM_PRED_FACT_TIME     - Predicted factorization time                      Default: -                OUT
!    DPARM_FACT_TIME          - Time for Numerical Factorization step (wallclock) Default: -                OUT
!    DPARM_SOLV_TIME          - Time for Solve step (wallclock)                   Default: -                OUT
!    DPARM_FACT_FLOPS         - Numerical Factorization flops (rate!)             Default: -                OUT
!    DPARM_SOLV_FLOPS         - Solve flops (rate!)                               Default: -                OUT
!    DPARM_RAFF_TIME          - Time for Refinement step (wallclock)              Default: -                OUT
!    DPARM_SIZE               - Dparm Size         IGNORE                         Default: -                IN
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_FILL_IN                  = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_MEM_MAX                  = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_EPSILON_REFINEMENT       = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_RELATIVE_ERROR           = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_EPSILON_MAGN_CTRL        = 10
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_ANALYZE_TIME             = 18
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_PRED_FACT_TIME           = 19
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_FACT_TIME                = 20
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_SOLV_TIME                = 21
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_FACT_FLOPS               = 22
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_SOLV_FLOPS               = 23
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_RAFF_TIME                = 24
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: DPARM_SIZE                     = 64
! !   Enum: API_TASK
! 
!   PaStiX step modes (index IPARM_START_TASK and IPARM_END_TASK)
! 
!   API_TASK_INIT       - Set default parameters
!   API_TASK_ORDERING   - Ordering
!   API_TASK_SYMBFACT   - Symbolic factorization
!   API_TASK_ANALYSE    - Tasks mapping and scheduling
!   API_TASK_NUMFACT    - Numerical factorization
!   API_TASK_SOLVE      - Numerical solve
!   API_TASK_REFINE     - Numerical refinement
!   API_TASK_CLEAN      - Clean

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_INIT                  = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_ORDERING              = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SYMBFACT              = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_ANALYSE               = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_NUMFACT               = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SOLVE                 = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_REFINE                = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_CLEAN                 = 7
! !   Enum: API_TASK_OLD
! 
!   API_TASK_SCOTCH     - Ordering
!   API_TASK_FAX        - Symbolic factorization
!   API_TASK_BLEND      - Tasks mapping and scheduling
!   API_TASK_SOPALIN    - Numerical factorization
!   API_TASK_UPDOWN     - Numerical solve
!   API_TASK_REFINEMENT - Numerical refinement

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SCOTCH                = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_FAX                   = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_BLEND                 = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_SOPALIN               = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_UPDOWN                = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TASK_REFINEMENT            = 6
! !    Enum: API_VERBOSE
! 
!    Verbose modes (index IPARM_VERBOSE)
! 
!    API_VERBOSE_NOT        - Silent mode, no messages
!    API_VERBOSE_NO         - Some messages
!    API_VERBOSE_YES        - Many messages
!    API_VERBOSE_CHATTERBOX - Like a gossip
!    API_VERBOSE_UNBEARABLE - Really talking too much...

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_NOT                = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_NO                 = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_YES                = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_CHATTERBOX         = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_VERBOSE_UNBEARABLE         = 4
! !   Enum: API_IO
! 
!   Check-points modes (index IPARM_IO)
! 
!   API_IO_NO         - No output or input
!   API_IO_LOAD       - Load ordering during ordering step and symbol matrix instead of symbolic factorisation.
!   API_IO_SAVE       - Save ordering during ordering step and symbol matrix instead of symbolic factorisation.
!   API_IO_LOAD_GRAPH - Load graph during ordering step.
!   API_IO_SAVE_GRAPH - Save graph during ordering step.
!   API_IO_LOAD_CSC   - Load CSC(d) during ordering step.
!   API_IO_SAVE_CSC   - Save CSC(d) during ordering step.

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_NO                      = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_LOAD                    = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_SAVE                    = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_LOAD_GRAPH              = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_SAVE_GRAPH              = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_LOAD_CSC                = 16
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_IO_SAVE_CSC                = 32
! !   Enum: API_RHS
! 
!   Right-hand-side modes (index IPARM_RHS)
! 
!   API_RHS_B - User's right hand side
!   API_RHS_1 - $ \forall i, X_i = 1 $
!   API_RHS_I - $ \forall i, X_i = i $
! 

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_B                      = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_1                      = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_I                      = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RHS_0                      = 3
! !   Enum: API_RAF
! 
!   Refinement modes (index IPARM_REFINEMENT)
! 
!   API_RAF_GMRES - GMRES
!   API_RAF_GRAD  - Conjugate Gradient ($LL^t$ or $LDL^t$ factorization)
!   API_RAF_PIVOT - Iterative Refinement (only for $LU$ factorization)
! 

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RAF_GMRES                  = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RAF_GRAD                   = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_RAF_PIVOT                  = 1
! !   Enum: API_FACT
! 
!   Factorization modes (index IPARM_FACTORISATION)
! 
!   API_FACT_LLT  - $LL^t$ Factorization
!   API_FACT_LDLT - $LDL^t$ Factorization
!   API_FACT_LU   - $LU$ Factorization
!   API_FACT_LDLH - $LDL^h$ hermitian factorization

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LLT                   = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LDLT                  = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LU                    = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_FACT_LDLH                  = 3
! !   Enum: API_SYM
! 
!   Symmetric modes (index IPARM_SYM)
! 
!   API_SYM_YES - Symmetric matrix
!   API_SYM_NO  - Nonsymmetric matrix
!   API_SYM_HER - Hermitian
! 

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_SYM_YES                    = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_SYM_NO                     = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_SYM_HER                    = 2
! !   Enum: API_ERASE_CSC
! 
!   CSC Management modes (index IPARM_FREE_CSCUSER and IPARM_FREE_CSCPASTIX)
! 
!   API_CSC_PRESERVE - Do not free the CSC
!   API_CSC_FREE     - Free the CSC when no longer needed

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_CSC_PRESERVE               = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_CSC_FREE                   = 1
! !   Enum: API_THREAD_MODE
! 
!   Comunication modes (index IPARM_NB_THREAD_COMM)
! 
!   API_THREAD_MULTIPLE      - All threads communicate.
!   API_THREAD_FUNNELED      - One thread perform all the MPI Calls.
!   API_THREAD_COMM_ONE      - One dedicated communication thread will receive messages.
!   API_THREAD_COMM_DEFINED  - Then number of threads receiving the messages is given by IPARM_NB_THREAD_COMM.
!   API_THREAD_COMM_NBPROC   - One communication thread per computation thread will receive messages.

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_MULTIPLE            = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_FUNNELED            = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_COMM_ONE            = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_COMM_DEFINED        = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_THREAD_COMM_NBPROC         = 16
! !   Enum: API_BIND_MODE
! 
!   Thread-binding modes (index IPARM_BINTHRD)
! 
!   API_BIND_NO   - Do not bind thread
!   API_BIND_AUTO - Default binding
!   API_BIND_TAB  - Use vector given by pastix_setBind

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_BIND_NO                    = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_BIND_AUTO                  = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_BIND_TAB                   = 2
! !   Enum: API_BOOLEAN
! 
!   Boolean modes (All boolean except IPARM_SYM)
! 
!   API_NO  - No
!   API_YES - Yes

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_NO                         = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_YES                        = 1
! !   Enum: API_TRACEFMT
! 
!   Trace modes (index IPARM_TRACEFMT)
! 
!   API_TRACE_PICL       - Use PICL trace format
!   API_TRACE_PAJE       - Use Paje trace format
!   API_TRACE_HUMREAD    - Use human-readable text trace format
!   API_TRACE_UNFORMATED - Unformated trace format

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_PICL                 = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_PAJE                 = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_HUMREAD              = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_TRACE_UNFORMATED           = 3
! !   Enum: API_ORDER
! 
!   Ordering modes (index IPARM_ORDERING)
! 
!   API_ORDER_SCOTCH   - Use \scotch{} ordering
!   API_ORDER_METIS    - Use \metis{} ordering
!   API_ORDER_PERSONAL - Apply user's permutation
!   API_ORDER_LOAD     - Load ordering from disk

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_SCOTCH               = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_METIS                = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_PERSONAL             = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_ORDER_LOAD                 = 3
! !   Enum: API_FLOAT
! 
!   Ordering modes (index IPARM_ORDERING)
! 
!   API_REALSINGLE    - Use \scotch{} ordering
!   API_REALDOUBLE    - Use \metis{} ordering
!   API_COMPLEXSINGLE - Apply user's permutation
!   API_COMPLEXDOUBLE - Load ordering from disk

INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_REALSINGLE                 = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_REALDOUBLE                 = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_COMPLEXSINGLE              = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: API_COMPLEXDOUBLE              = 3
! !   Enum: MODULES
! 
!   Module Identification number.
! 
!   If an error occurs, error value is set to
!   MODULE + EER_NUMBER.
! 
!   User can catch error by computing iparm[IPARM_ERROR_NUMBER]%100.
! 
!   MODULE can be catch by computing iparm[IPARM_ERROR_NUMBER] - iparm[IPARM_ERROR_NUMBER]%100.
! 
!   MOD_UNKNOWN - Unknown module
!   MOD_SOPALIN - Numerical factorisation module
!   MOD_BLEND   - Analysing module
!   MOD_SCOTCH  - Scotch module
!   MOD_FAX     - Symbolic factorisation module
!   MOD_ORDER   - Order module
!   MOD_COMMON  - Common module
!   MOD_SI      -
!   MOD_GRAPH   - Graph module
!   MOD_SYMBOL  - Symbol structure module
!   MOD_KASS    - Kass module
!   MOD_BUBBLE  - Bubble
!   MOD_MURGE   - Murge
! 
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_UNKNOWN                    = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SOPALIN                    = 100
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_BLEND                      = 200
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SCOTCH                     = 300
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_FAX                        = 400
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_ORDER                      = 500
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_COMMON                     = 600
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SI                         = 700
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_GRAPH                      = 800
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_SYMBOL                     = 900
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_KASS                       = 1000
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_BUBBLE                     = 1100
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MOD_MURGE                      = 1200
!  Enum: ERR_NUMBERS! 
!    Error Numbers
! 
!    NO_ERR             - No error
!    UNKNOWN_ERR        - Unknown error
!    ALLOC_ERR          - Allocation error
!    ASSERT_ERR         - Error in one assertion
!    NOTIMPLEMENTED_ERR - Not implemented feature
!    OUTOFMEMORY_ERR    - Not enough memory (OOC)
!    THREAD_ERR         - Error with threads
!    INTERNAL_ERR       - Internal error
!    BADPARAMETER_ERR   - Bad parameters given
!    FILE_ERR           - Error in In/Out operations
!    BAD_DEFINE_ERROR   - Error with defines during compilation
!    INTEGER_TYPE_ERR   - Error with integer types
!    IO_ERR             - Error with input/output
!    MATRIX_ERR         - Wrongly defined matrix
!    FLOAT_TYPE_ERR     - Wrong type of floating point values
!    STEP_ORDER_ERR     - Error in ordering
!    MPI_ERR            - Error with MPI calls
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: NO_ERR                         = 0
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: UNKNOWN_ERR                    = 1
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: ALLOC_ERR                      = 2
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: ASSERT_ERR                     = 3
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: NOTIMPLEMENTED_ERR             = 4
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: OUTOFMEMORY_ERR                = 5
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: THREAD_ERR                     = 6
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: INTERNAL_ERR                   = 7
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: BADPARAMETER_ERR               = 8
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: FILE_ERR                       = 9
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: BAD_DEFINE_ERR                 = 10
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: INTEGER_TYPE_ERR               = 11
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: IO_ERR                         = 12
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MATRIX_ERR                     = 13
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: FLOAT_TYPE_ERR                 = 14
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: STEP_ORDER_ERR                 = 15
INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: MPI_ERR                        = 16

