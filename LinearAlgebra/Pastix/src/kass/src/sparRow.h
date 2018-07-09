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
/************************************************************/
/**                                                        **/
/**   NAME       : sparRow.h                               **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 17 May 2005     **/
/**                                                        **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/

#ifndef SPARSE_ROW_
#define SPARSE_ROW_

typedef struct SparRow *csptr;
typedef struct SparRow {
  /*--------------------------------------------- 
    | C-style CSR format - used internally
    | for all matrices in CSR format 
    |---------------------------------------------*/

  INT      n;
  INT     *nnzrow; /* length of each row                               */
  double **ma;     /* pointer-to-pointer to store nonzero entries      */
  INT    **ja;     /* pointer-to-pointer to store column indices       */
  INT      inarow; /* This flag means the matrix has been allocated as 
		      a single block of memory; it must be desallocated 
		      in consequence                                   */
  INT     *jatab;  /* Used if inarow == 1 to store the the matrix in 
		      two contigue block of memory                     */
  double  *matab;
} SparMat;


INT initCS (csptr amat, INT len);
INT cleanCS(csptr amat);
INT CSnnz  (csptr mat);
INT CS_Perm(csptr mat, INT *perm);

#endif
