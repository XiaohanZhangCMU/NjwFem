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
  File: cscd_utils.h

  Several operations on CSCD.

 */

#ifndef CSCD_UTILS_H
#define CSCD_UTILS_H

#ifndef _GLIBCXX_HAVE_COMPLEX_H
#  define _GLIBCXX_HAVE_COMPLEX_H 0
#endif

#ifdef   __cplusplus
#  if (_GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
#    define  COMPLEX  std::complex<float>
#    define  DCOMPLEX std::complex<double>
#  endif
#else /* not __cplusplus */
#  if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__)
#    define  COMPLEX float complex
#    define  DCOMPLEX double complex
#  endif
#endif /* not __cplusplus */

/*
  Enum: CSCD_OPERATIONS

  Operation when adding CSCD

  CSCD_ADD  - Add coefficient values.
  CSCD_KEEP - Keep value from first CSCD.
  CSCD_MAX  - Keep maximum of first and second CSCD.
  CSCD_MIN  - Keep minimum of first and second CSCD.
  CSCD_OVW  - Overwrite with second CSCD value.
*/
enum CSCD_OPERATIONS {
  CSCD_ADD,
  CSCD_KEEP,
  CSCD_MAX,
  CSCD_MIN,
  CSCD_OVW
};
typedef enum CSCD_OPERATIONS CSCD_OPERATIONS_t;

/*
  Enum: CSC_DISPATCH_OP

  Operation when dispatching the csc into a cscd

  CSC_DISP_SIMPLE - Reparts linearly the columns over the proc.
  CSC_DISP_CYCLIC - Reparts cyclicly the columns over the proc.

*/
enum CSC_DISPATCH_OP {
  CSC_DISP_SIMPLE,
  CSC_DISP_CYCLIC
};
typedef enum CSC_DISPATCH_OP CSCD_DISPATCH_OP_t;

/* Section: Functions */
#if (defined FLOAT)
/*
   Function: csc_dispatch

   Distribute a CSC to a CSCD

   Parameters:
      gN                - global number of columns
      gcolptr           - global starting index of each column in grows ans gavals.
      grows             - global rows of each element.
      gavals            - global values of each element.
      gperm             - global permutation tabular.
      ginvp             - global reverse permutation tabular.
      lN                - local number of columns (output).
      lcolptr           - starting index of each local column (output).
      lrowptr           - row number of each local element (output).
      lavals            - values of each local element (output).
      lrhs              - local part of the right hand side (output).
      lperm             - local part of the permutation tabular (output).
      loc2glob          - global numbers of local columns (before permutation).
      dispatch          - choose how to dispatch the csc
      pastix_comm       - PaStiX MPI communicator.
*/
void csc_dispatch(INT  gN, INT *  gcolptr, INT *  grow, FLOAT *  gavals,
                  FLOAT *  grhs, INT *  gperm, INT *  ginvp,
                  INT *lN, INT ** lcolptr, INT ** lrow, FLOAT ** lavals,
                  FLOAT ** lrhs, INT ** lperm,
                  INT **loc2glob, int dispatch, MPI_Comm pastix_comm);
#endif
void s_csc_dispatch(INT  gN, INT *  gcolptr, INT *  grow, float *  gavals,
                    float *  grhs, INT *  gperm, INT *  ginvp,
                    INT *lN, INT ** lcolptr, INT ** lrow, float ** lavals,
                    float ** lrhs, INT ** lperm,
                    INT **loc2glob, int dispatch, MPI_Comm pastix_comm);
void d_csc_dispatch(INT  gN, INT *  gcolptr, INT *  grow, double *  gavals,
                    double *  grhs, INT *  gperm, INT *  ginvp,
                    INT *lN, INT ** lcolptr, INT ** lrow, double ** lavals,
                    double ** lrhs, INT ** lperm,
                    INT **loc2glob, int dispatch, MPI_Comm pastix_comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
void c_csc_dispatch(INT  gN, INT *  gcolptr, INT *  grow, COMPLEX *  gavals,
                    COMPLEX *  grhs, INT *  gperm, INT *  ginvp,
                    INT *lN, INT ** lcolptr, INT ** lrow, COMPLEX ** lavals,
                    COMPLEX ** lrhs, INT ** lperm,
                    INT **loc2glob, int dispatch, MPI_Comm pastix_comm);
void z_csc_dispatch(INT  gN, INT *  gcolptr, INT *  grow, DCOMPLEX *  gavals,
                    DCOMPLEX *  grhs, INT *  gperm, INT *  ginvp,
                    INT *lN, INT ** lcolptr, INT ** lrow, DCOMPLEX ** lavals,
                    DCOMPLEX ** lrhs, INT ** lperm,
                    INT **loc2glob, int dispatch, MPI_Comm pastix_comm);
#endif
#if (defined FLOAT)
/*
  Function: csc_cyclic_distribution

  Distribute the CSC cyclicaly.

  Parameters:
    column      - column number to distribute
    columnnbr   - Number of colmuns.
    pastix_comm - PaStiX MPI communicator

  Return:
    owner of the column (column%commSize)
*/
INT csc_cyclic_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
#endif
INT Scsc_cyclic_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
INT Dcsc_cyclic_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
INT Ccsc_cyclic_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
INT Zcsc_cyclic_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
#endif
#if (defined FLOAT)
/*
  Function: csc_simple_distribution

  Distribute the CSC.
  First columns are for first proc and so on.

  Parameters:
    column      - column number to distribute
    columnnbr   - Number of colmuns.
    pastix_comm - PaStiX MPI communicator

  Return:
    owner of the column (column/commSize)
*/
INT csc_simple_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);

#endif
INT s_csc_simple_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
INT d_csc_simple_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
INT c_csc_simple_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
INT z_csc_simple_distribution(INT column, INT columnnbr, MPI_Comm pastix_comm);
#endif

#if (defined FLOAT)
/*
  Function: cscd_symgraph

  Check if the CSCD graph is symetric.

  Parameters:
    n           - Number of local columns
    ia          - Starting index of each columns in *ja* and *a*
    ja          - Row of each element.
    a           - Values of each element.
    newn        - New number of local columns
    newia       - Starting index of each columns in *newja* and *newa*
    newja       - Row of each element.
    newa        - Values of each element.
    l2g         - global number of each local column.
    malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int cscd_symgraph(INT      n, INT *     ia, INT *     ja, FLOAT *     a,
                  INT * newn, INT ** newia, INT ** newja, FLOAT ** newa,
                  INT *     l2g,  MPI_Comm comm);
#endif
int s_cscd_symgraph(INT      n, INT *     ia, INT *     ja, float *     a,
                    INT * newn, INT ** newia, INT ** newja, float ** newa,
                    INT *     l2g,  MPI_Comm comm);
int d_cscd_symgraph(INT      n, INT *     ia, INT *     ja, double *     a,
                    INT * newn, INT ** newia, INT ** newja, double ** newa,
                    INT *     l2g,  MPI_Comm comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
int c_cscd_symgraph(INT      n, INT *     ia, INT *     ja, COMPLEX *     a,
                    INT * newn, INT ** newia, INT ** newja, COMPLEX ** newa,
                    INT *     l2g,  MPI_Comm comm);
int z_cscd_symgraph(INT      n, INT *     ia, INT *     ja, DCOMPLEX *     a,
                    INT * newn, INT ** newia, INT ** newja, DCOMPLEX ** newa,
                    INT *     l2g,  MPI_Comm comm);
#endif
#if (defined FLOAT)
/*
  Function: cscd_addlocal

  Add second cscd to first cscd into third cscd (unallocated)

  Parameters:
    n           - First cscd size
    ia          - First cscd starting index of each column in *ja* and *a*
    ja          - Row of each element in first CSCD
    a           - value of each cscd in first CSCD (can be NULL)
    l2g         - local 2 global column numbers for first cscd
    addn        - CSCD to add size
    addia       - CSCD to add starting index of each column in *addja* and *adda*
    addja       - Row of each element in second CSCD
    adda        - value of each cscd in second CSCD (can be NULL -> add 0)
    addl2g      - local 2 global column numbers for second cscd
    newn        - new cscd size (same as first)
    newia       - CSCD to add starting index of each column in *newja* and *newwa*
    newja       - Row of each element in third CSCD
    newa        - value of each cscd in third CSCD
    OP          - Operation to manage common CSCD coefficients.
    dof         - Number of degrees of freedom.
*/

int cscd_addlocal(INT   n   , INT *  ia   , INT *  ja   , FLOAT *  a   , INT * l2g,
      INT   addn, INT *  addia, INT *  addja, FLOAT *  adda, INT * addl2g,
      INT * newn, INT ** newia, INT ** newja, FLOAT ** newa, CSCD_OPERATIONS_t OP, int dof);
#endif
int s_cscd_addlocal(INT   n   , INT *  ia   , INT *  ja   , float *  a   , INT * l2g,
      INT   addn, INT *  addia, INT *  addja, float *  adda, INT * addl2g,
      INT * newn, INT ** newia, INT ** newja, float ** newa, CSCD_OPERATIONS_t OP, int dof);

int d_cscd_addlocal(INT   n   , INT *  ia   , INT *  ja   , double *  a   , INT * l2g,
      INT   addn, INT *  addia, INT *  addja, double *  adda, INT * addl2g,
      INT * newn, INT ** newia, INT ** newja, double ** newa, CSCD_OPERATIONS_t OP, int dof);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
int c_cscd_addlocal(INT   n   , INT *  ia   , INT *  ja,
                    COMPLEX *  a   , INT * l2g,
                    INT   addn, INT *  addia, INT *  addja,
                    COMPLEX *  adda, INT * addl2g,
                    INT * newn, INT ** newia, INT ** newja,
                    COMPLEX ** newa, CSCD_OPERATIONS_t OP, int dof);

int z_cscd_addlocal(INT   n   , INT *  ia   , INT *  ja,
                    DCOMPLEX *  a   , INT * l2g,
                    INT   addn, INT *  addia, INT *  addja,
                    DCOMPLEX *  adda, INT * addl2g,
                    INT * newn, INT ** newia, INT ** newja,
                    DCOMPLEX ** newa, CSCD_OPERATIONS_t OP, int dof);
#endif

#if (defined FLOAT)
/**
    Function: csc2cscd

    Transform a csc to a cscd.
    Allocate the CSCD.
    If grhs == NULL forget right hand side part.
    If gperm == NULL forget permutation and reverse permutation part.

    Parameters:
      gN       - global number of columns
      gcolptr  - global starting index of each column in grows ans gavals.
      grows    - global rows of each element.
      gavals   - global values of each element.
      gperm    - global permutation tabular.
      ginvp    - global reverse permutation tabular.
      lN       - local number of columns.
      lcolptr  - starting index of each local column.
      lrowptr  - row number of each local element.
      lavals   - values of each local element.
      lrhs     - local part of the right hand side (output).
      lperm    - local part of the permutation tabular (output).
      linvp    - local part of the reverse permutation tabular (output).
      loc2glob - global numbers of local columns (before permutation).
*/
void  csc2cscd(INT gN, INT *  gcolptr, INT *  grow, FLOAT *  gavals,
               FLOAT *  grhs, INT *  gperm, INT *  ginvp,
               INT lN, INT ** lcolptr, INT ** lrow, FLOAT ** lavals,
               FLOAT ** lrhs, INT ** lperm, INT ** linvp,
               INT *loc2glob);
#endif
void  s_csc2cscd(INT gN, INT *  gcolptr, INT *  grow, float *  gavals,
                 float *  grhs, INT *  gperm, INT *  ginvp,
                 INT lN, INT ** lcolptr, INT ** lrow, float ** lavals,
                 float ** lrhs, INT ** lperm, INT ** linvp,
                 INT *loc2glob);
void  d_csc2cscd(INT gN, INT *  gcolptr, INT *  grow, double *  gavals,
                 double *  grhs, INT *  gperm, INT *  ginvp,
                 INT lN, INT ** lcolptr, INT ** lrow, double ** lavals,
                 double ** lrhs, INT ** lperm, INT ** linvp,
                 INT *loc2glob);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
void  c_csc2cscd(INT gN, INT *  gcolptr, INT *  grow, COMPLEX *  gavals,
                 COMPLEX *  grhs, INT *  gperm, INT *  ginvp,
                 INT lN, INT ** lcolptr, INT ** lrow, COMPLEX ** lavals,
                 COMPLEX ** lrhs, INT ** lperm, INT ** linvp,
                 INT *loc2glob);
void  z_csc2cscd(INT gN, INT *  gcolptr, INT *  grow, DCOMPLEX *  gavals,
                 DCOMPLEX *  grhs, INT *  gperm, INT *  ginvp,
                 INT lN, INT ** lcolptr, INT ** lrow, DCOMPLEX ** lavals,
                 DCOMPLEX ** lrhs, INT ** lperm, INT ** linvp,
                 INT *loc2glob);
#endif

#if (defined FLOAT)
/**
    Function: cscd2csc

    Transform a cscd to a csc.
    colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.

    Parameters:
       lN          - number of local column.
       lcolptr     - starting index of each local column in row and avals.
       lrow        _ row number of each local element.
       lavals      - values of each local element.
       lrhs        - local part of the right hand side.
       lperm       - local part of the permutation tabular.
       linvp       - local part of the reverse permutation tabular.
       gN          - global number of columns (output).
       gcolptr     - starting index of each column in row2 and avals2 (output).
       grow        - row number of each element (output).
       gavals      - values of each element (output).
       grhs        - global right hand side (output).
       gperm       - global permutation tabular (output).
       ginvp       - global reverse permutation tabular (output).
       loc2glob    - global number of each local column.
       pastix_comm - PaStiX MPI communicator.

*/

void  cscd2csc(INT  lN, INT *  lcolptr, INT * lrow, FLOAT * lavals,
               FLOAT * lrhs, INT * lperm, INT * linvp,
               INT *gN, INT ** gcolptr, INT **grow, FLOAT **gavals,
               FLOAT **grhs, INT **gperm, INT **ginvp,
               INT *loc2glob, MPI_Comm pastix_comm);
#endif
void  s_cscd2csc(INT  lN, INT *  lcolptr, INT * lrow, float * lavals,
                 float * lrhs, INT * lperm, INT * linvp,
                 INT *gN, INT ** gcolptr, INT **grow, float **gavals,
                 float **grhs, INT **gperm, INT **ginvp,
                 INT *loc2glob, MPI_Comm pastix_comm);
void  d_cscd2csc(INT  lN, INT *  lcolptr, INT * lrow, double * lavals,
                 double * lrhs, INT * lperm, INT * linvp,
                 INT *gN, INT ** gcolptr, INT **grow, double **gavals,
                 double **grhs, INT **gperm, INT **ginvp,
                 INT *loc2glob, MPI_Comm pastix_comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
void  c_cscd2csc(INT  lN, INT *  lcolptr, INT * lrow, COMPLEX * lavals,
                 COMPLEX * lrhs, INT * lperm, INT * linvp,
                 INT *gN, INT ** gcolptr, INT **grow, COMPLEX **gavals,
                 COMPLEX **grhs, INT **gperm, INT **ginvp,
                 INT *loc2glob, MPI_Comm pastix_comm);
void  z_cscd2csc(INT  lN, INT *  lcolptr, INT * lrow, DCOMPLEX * lavals,
                 DCOMPLEX * lrhs, INT * lperm, INT * linvp,
                 INT *gN, INT ** gcolptr, INT **grow, DCOMPLEX **gavals,
                 DCOMPLEX **grhs, INT **gperm, INT **ginvp,
                 INT *loc2glob, MPI_Comm pastix_comm);
#endif

#if (defined FLOAT)
/*
 * Function: cscd_redispatch
 *
 * Redistribute the first cscd into a new one using *dl2g*.
 *
 * - gather all new loc2globs on all processors.
 * - allocate *dia*, *dja* and *da*.
 * - Create new CSC for each processor and send it.
 * - Merge all new CSC to the new local CSC with <cscd_addlocal_int>.
 *
 * If communicator size is one, check that n = dn and
 * l2g = dl2g and simply create a copy of the first cscd.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   rhs         - right-hand-side member corresponding to the first CSCD
 *                 (can be NULL)
 *   nrhs        - number of right-hand-side.
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - New cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in new CSCD
 *   da          - value of each cscd in new CSCD
 *   rhs         - right-hand-side member corresponding to the new CSCD
 *   dl2g        - local 2 global column numbers for new cscd
 *   comm        - MPI communicator
 *
 * Returns:
 *   EXIT_SUCCESS - If all goes well
 *   EXIT_FAILURE - If commsize = 1 and *n* != *dn* or *l2g* != *dl2g*.
 */
int cscd_redispatch(INT   n, INT *   ia, INT *   ja, FLOAT *   a,
                    FLOAT *  rhs,  INT nrhs, INT *   l2g,
                    INT  dn, INT ** dia, INT ** dja, FLOAT ** da,
                    FLOAT ** drhs,  INT *  dl2g,
                    MPI_Comm comm);
#endif
int s_cscd_redispatch(INT   n, INT *   ia, INT *   ja, float *   a,
                      float *  rhs,  INT nrhs,   INT *   l2g,
                      INT  dn, INT ** dia, INT ** dja, float ** da,
                      float ** drhs,  INT *  dl2g,
                      MPI_Comm comm);
int d_cscd_redispatch(INT   n, INT *   ia, INT *   ja, double *   a,
                      double *  rhs,  INT nrhs,   INT *   l2g,
                      INT  dn, INT ** dia, INT ** dja, double ** da,
                      double ** drhs,  INT *  dl2g,
                      MPI_Comm comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
int c_cscd_redispatch(INT   n, INT *   ia, INT *   ja, COMPLEX *   a,
                      COMPLEX *  rhs,  INT nrhs,   INT *   l2g,
                      INT  dn, INT ** dia, INT ** dja, COMPLEX ** da,
                      COMPLEX ** drhs,  INT *  dl2g,
                      MPI_Comm comm);
int z_cscd_redispatch(INT   n, INT *   ia, INT *   ja, DCOMPLEX *   a,
                      DCOMPLEX *  rhs,  INT nrhs,   INT *   l2g,
                      INT  dn, INT ** dia, INT ** dja, DCOMPLEX ** da,
                      DCOMPLEX ** drhs,  INT *  dl2g,
                      MPI_Comm comm);
#endif

#if (defined FLOAT)
/*
   Function: cscd_save

   save a distributed csc to disk.
   files are called $(filename) and $(filename)$(RANK)
   if filename is NULL then filename = cscd_matrix.

   file filename contains the number of processors/files
   on first line. Then each line contain the name of each file
   (here $(filename)$(RANK)).



   Parameters:
     n           - Number of local columns
     ia          - First cscd starting index of each column in *ja* and *a*
     ja          - Row of each element in first CSCD
     a           - value of each cscd in first CSCD (can be NULL)
     rhs         - Right hand side.
     l2g         - local 2 global column numbers for first cscd
     dof         - Number of degrees of freedom
     filename    - name of the files.
     comm        - MPI communicator

*/
int cscd_save(INT n, INT *ia, INT *ja, FLOAT * a, FLOAT * rhs, INT* l2g,
              int dof, const char * filename, MPI_Comm comm);
#endif
int s_cscd_save(INT n, INT *ia, INT *ja, float * a, float * rhs, INT* l2g,
                int dof, const char * filename, MPI_Comm comm);
int d_cscd_save(INT n, INT *ia, INT *ja, double * a, double * rhs, INT* l2g,
                int dof, const char * filename, MPI_Comm comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
int c_cscd_save(INT n, INT *ia, INT *ja, COMPLEX * a, COMPLEX * rhs,
                INT* l2g,  int dof,  const char * filename, MPI_Comm comm);
int z_cscd_save(INT n, INT *ia, INT *ja, DCOMPLEX * a, DCOMPLEX * rhs,
                INT* l2g,  int dof,  const char * filename, MPI_Comm comm);
#endif
#if (defined FLOAT)
/*
   Function: cscd_load

   Loads a distributed csc from disk.
   if filename is NULL then filename = cscd_matrix.

   Parameters:
     n           - Number of local columns
     ia          - First cscd starting index of each column in *ja* and *a*
     ja          - Row of each element in first CSCD
     a           - value of each cscd in first CSCD (can be NULL)
     rhs         - Right hand side.
     l2g         - local 2 global column numbers for first cscd
     filename    - name of the files.
     comm        - MPI communicator

*/
int cscd_load(INT *n, INT ** ia, INT ** ja, FLOAT ** a, FLOAT ** rhs, INT ** l2g,
              const char * filename, MPI_Comm mpi_comm);
#endif
int s_cscd_load(INT *n, INT ** ia, INT ** ja, float ** a, float ** rhs, INT ** l2g,
                const char * filename, MPI_Comm mpi_comm);
int d_cscd_load(INT *n, INT ** ia, INT ** ja, double ** a, double ** rhs,
                INT ** l2g, const char * filename, MPI_Comm mpi_comm);
#if (defined _COMPLEX_H || defined _H_COMPLEX || defined __COMPLEX__ || _GLIBCXX_HAVE_COMPLEX_H == 1 || defined __STD_COMPLEX)
int c_cscd_load(INT *n, INT ** ia, INT ** ja, COMPLEX ** a,
                COMPLEX ** rhs, INT ** l2g, const char * filename,
                MPI_Comm mpi_comm);
int z_cscd_load(INT *n, INT ** ia, INT ** ja, DCOMPLEX ** a,
                DCOMPLEX ** rhs, INT ** l2g, const char * filename,
                MPI_Comm mpi_comm);
#endif
#undef COMPLEX
#undef DCOMPLEX
#endif /* CSCD_UTILS_H */
