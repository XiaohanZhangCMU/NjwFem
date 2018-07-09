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
  Function: add_two_floats

  Adds two integers.

  Parameters :
    a - first integer
    b - second integer

  Returns: a + b
 */
FLOAT add_two_floats(FLOAT a, FLOAT b);

/*
  Function: keep_first

  Returns first integer.

  Parameters :
    a - first integer
    b - second integer

  Returns: a
 */
FLOAT keep_first(FLOAT a, FLOAT b);

/*
  Function: keep_last

  Returns last integer.

  Parameters :
    a - first integer
    b - second integer

  Returns: b
 */
FLOAT keep_last(FLOAT a, FLOAT b);

/*
  Function: get_max

  Returns maximum value from two integers.

  Parameters :
    a - first integer
    b - second integer

  Returns: MAX(a,b)
 */
FLOAT get_max(FLOAT a, FLOAT b);

/*
  Function: get_min

  Returns minimum value from two integers.

  Parameters :
    a - first integer
    b - second integer

  Returns: MIN(a,b)
 */
FLOAT get_min(FLOAT a, FLOAT b);


int cscd_addlocal_int(INT   n   , INT *  ia   , INT *  ja   , FLOAT *  a   , INT * l2g,
		      INT   addn, INT *  addia, INT *  addja, FLOAT *  adda, INT * addl2g,
		      INT * newn, INT ** newia, INT ** newja, FLOAT ** newa, 
		      FLOAT (*add_fct)(FLOAT , FLOAT), int dof, int malloc_flag);

/*
 * Function: cscd_redispatch_int
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
 *   rhs         - right-hand-side member corresponding to the first CSCD (can be NULL)
 *   nrhs        - number of right-hand-side.
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - New cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in new CSCD
 *   da          - value of each cscd in new CSCD
 *   rhs         - right-hand-side member corresponding to the new CSCD
 *   dl2g        - local 2 global column numbers for new cscd
 *   malloc_flag - Internal (API_YES) or external (API_NO) malloc use.
 *   comm        - MPI communicator
 *
 * Returns:
 *   EXIT_SUCCESS - If all goes well
 *   EXIT_FAILURE - If commsize = 1 and *n* != *dn* or *l2g* != *dl2g*.
 */
int cscd_redispatch_int(INT   n, INT *   ia, INT *   ja, FLOAT *   a, FLOAT *  rhs,  INT nrhs, INT *   l2g, 
			INT  dn, INT ** dia, INT ** dja, FLOAT ** da, FLOAT ** drhs, INT *  dl2g,
			int  malloc_flag, MPI_Comm comm);

int cscd_symgraph_int(INT      n, INT *      ia, INT *        ja, FLOAT *     a,
		      INT * newn, INT **  newia, INT **    newja, FLOAT ** newa,
		      INT *  l2g, MPI_Comm comm, int malloc_flag);


/*
  Function: cscd_build_g2l

  Construct global to local tabular containing local number of global columns
  if one column is local, and -owner if column is not local.

  For i in 0, gN
     g2l[i] = i local number if i is local
     g2l[i] = -p if p is the owner of the column i
  
  Parameters:
    n        - Number of local columns
    colptr   - Starting index of each columns in *ja*
    rows     - Row of each element.
    values   - Value of each element.
    l2g      - global number of each local column.
    correct  - Flag indicating if we can correct the symmetry.
    dof      - Number of degrees of freedom.
    comm     - MPI communicator
 */
int cscd_build_g2l(INT       ncol,
		   INT      *loc2glob,
		   MPI_Comm  comm,
		   INT      *gN,
		   INT     **g2l);
/* 
   Function: cscd_noDiag
   
   Removes diagonal elements from a CSCD.
   *ja* and *a* can be reallocated to 
   ia[n]-1 elements after this call.
   
   Parameters:
     n           - Number of local columns
     ia          - First cscd starting index of each column in *ja* and *a*
     ja          - Row of each element in first CSCD
     a           - value of each cscd in first CSCD (can be NULL)
     l2g         - local 2 global column numbers for first cscd

*/
int cscd_noDiag(INT n, INT *ia, INT *ja, FLOAT * a, INT * l2g);


/*
  Function: cscd_checksym

  Check if the CSCD graph is symetric.

  Parameters:
    n   - Number of local columns
    ia  - Starting index of each columns in *ja*
    ja  - Row of each element.
    l2g - global number of each local column.
    correct  - Flag indicating if we can correct the symmetry.
    alloc    - indicate if allocation on CSC uses internal malloc. 
    dof      - Number of degrees of freedom.
    comm     - MPI communicator
*/
int cscd_checksym(INT      n, 
		  INT     *colptr, 
		  INT    **rows, 
		  FLOAT  **values, 
		  INT     *l2g, 
		  int      correct,
		  int      alloc,
		  int      dof,
		  MPI_Comm comm);


/** 
    Function: cscd2csc_int        
    
    Transform a cscd to a csc.
    colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.

    External function, allocation are not of the internal type.

    Parameters:
       lN          - number of local column.
       lcolptr     - starting index of each local column in row and avals.
       lrow        _ row number of each local element.
       lavals      - values of each local element.
       lrhs        - local part of the right hand side.
       lperm       - local part of the permutation tabular.
       linvp       - Means nothing, to suppress.
       gN          - global number of columns (output).
       gcolptr     - starting index of each column in row2 and avals2 (output).
       grow        - row number of each element (output).
       gavals      - values of each element (output).
       grhs        - global right hand side (output).
       gperm       - global permutation tabular (output).
       ginvp       - global reverse permutation tabular (output).
       loc2glob    - global number of each local column.
       pastix_comm - PaStiX MPI communicator.
       intern_flag - Decide if malloc will use internal or external macros.

*/
void  cscd2csc_int(INT  lN, INT *  lcolptr, INT * lrow, FLOAT * lavals, 
		   FLOAT * lrhs, INT * lperm, INT * linvp, 
		   INT *gN, INT ** gcolptr, INT **grow, FLOAT **gavals, 
		   FLOAT **grhs, INT **gperm, INT **ginvp, 
		   INT *loc2glob, MPI_Comm pastix_comm, int intern_flag);

/*
  Function: cscd_redispatch_scotch

  Renumber the columns to have first columns on first proc for Scotch

  Parameters:
    n           - Number of local columns
    ia          - First cscd starting index of each column in *ja* and *a*
    ja          - Row of each element in first CSCD
    a           - value of each cscd in first CSCD (can be NULL)
    l2g         - local 2 global column numbers for first cscd
    dn          - Number of local columns
    dia         - First cscd starting index of each column in *ja* and *a*
    dja         - Row of each element in first CSCD
    da          - value of each cscd in first CSCD (can be NULL)
    l2g         - local 2 global column numbers for first cscd
    comm        - MPI communicator

  Returns:
    EXIT_SUCCESS if already well distributed, 2 if redistributed
    
 */
int cscd_redispatch_scotch(INT   n, INT *   ia, INT *   ja, FLOAT *   a, INT *   l2g, 
			   INT *dn, INT ** dia, INT ** dja, FLOAT ** da, INT ** dl2g,
			   MPI_Comm comm);
