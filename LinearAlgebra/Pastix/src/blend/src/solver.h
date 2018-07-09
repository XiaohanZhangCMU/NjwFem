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
/**   NAME       : solver.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the solver matrix.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     28 oct 1998     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to     06 jun 2002     **/
/**                                                        **/
/************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

/*
**  The type and structure definitions.
*/

#define COMP_1D                     0
#define DIAG                        1
#define E1                          2
#define E2                          3
#define DRUNK                       4

typedef struct Task_ {
  INT                       taskid;               /*+ COMP_1D DIAG E1 E2                                     +*/
  INT                       prionum;              /*+ Priority value for the factorization                   +*/
  INT                       prionum2;             /*+ Priority value for solve steps                         +*/
  INT                       cblknum;              /*+ Attached column block                                  +*/
  INT                       bloknum;              /*+ Attached block                                         +*/
  INT volatile              ftgtcnt;              /*+ Number of fan-in targets                               +*/
  INT volatile              ctrbcnt;              /*+ Total number of contributions                          +*/
  volatile BlockTarget *    btagptr;              /*+ si non local, pointeur sur la structure (NB reception) +*/
  INT                       indnum;               /*+ For E2 (COMP_1D), index of ftgt (>0) else if local = -taskdest
                                                      For DIAG and E1 , index of btag (>0) if there is a
						      local one it must be the first of the chain of local E1   +*/
  INT                       tasknext;             /*+ chainage des E1 ou E2, si fin = -1 => liberer les btagptr +*/
  INT                       taskmstr;             /*+ Master task for E1 or E2 tasks                         +*/
                                                  /*+ Index of DIAG (or first E1) task for E1                +*/
                                                  /*+ Index of E1 (or first E2) task for E2                  +*/
#if (defined PASTIX_DYNSCHED) || (defined TRACE_SOPALIN)
  INT                       threadid;             /*+ Index of the bubble which contains the task +*/
  INT                       cand;		  /*+ Thread candidate in static version          +*/
#endif
#ifdef TRACE_SOPALIN
  INT                       fcandnum;             /*+ First thread candidate                      +*/
  INT                       lcandnum;		  /*+ Last thread candidate                       +*/
  INT                       id;                   /*+ Global cblknum of the attached column block +*/
#endif
} Task;

/*+ Solver column block structure. +*/

typedef struct SolverCblk_  {
  INT                       fcolnum;              /*+ First column index                     +*/
  INT                       lcolnum;              /*+ Last column index (inclusive)          +*/
  INT                       bloknum;              /*+ First block in column (diagonal)       +*/
  INT                       stride;               /*+ Column block stride                    +*/
  INT			    color;		  /*+ Color of column block (PICL trace)     +*/
  INT                       procdiag;             /*+ Processor owner of diagonal block      +*/
  INT                       cblkdiag;             /*+ Column block owner of diagonal block   +*/
  FLOAT * restrict          coeftab;              /*+ Coefficients access vector             +*/
  FLOAT * restrict          ucoeftab;             /*+ Coefficients access vector             +*/
} SolverCblk; 

/*+ Solver block structure. +*/

typedef struct SolverBlok_ {
  INT                       frownum;              /*+ First row index            +*/
  INT                       lrownum;              /*+ Last row index (inclusive) +*/
  INT                       cblknum;              /*+ Facing column block        +*/
  INT                       levfval;              /*+ Level-of-fill value        +*/
  INT                       coefind;              /*+ Index in coeftab           +*/
} SolverBlok;

/*+ Solver matrix structure. +*/

typedef struct SolverMatrix_ {
  INT                       baseval;              /*+ Base value for numberings                 +*/
  INT                       nodenbr;              /*+ Number of nodes in matrix                 +*/
  INT                       cblknbr;              /*+ Number of column blocks                   +*/
  INT                       bloknbr;              /*+ Number of blocks                          +*/
  SolverCblk * restrict     cblktab;              /*+ Array of solver column blocks             +*/
  SolverBlok * restrict     bloktab;              /*+ Array of solver blocks                    +*/
  INT                       coefnbr;              /*+ Number of coefficients                    +*/

  INT                       ftgtnbr;              /*+ Number of fanintargets                    +*/
  INT                       ftgtcnt;              /*+ Number of fanintargets to receive         +*/
  FanInTarget * restrict    ftgttab;              /*+ Fanintarget access vector                 +*/

  INT                       coefmax;              /*+ Working block max size (cblk coeff 1D)    +*/
  INT                       bpftmax;              /*+ Maximum of block size for btag to receive +*/
  INT                       cpftmax;              /*+ Maximum of block size for ftgt to receive +*/
  INT                       nbftmax;              /*+ Maximum block number in ftgt              +*/
  INT                       arftmax;              /*+ Maximum block area in ftgt                +*/

  INT                       clustnum;             /*+ current processor number                  +*/
  INT                       clustnbr;             /*+ number of processors                      +*/
  INT                       procnbr;              /*+ Number of physical processor used         +*/
  INT                       thrdnbr;              /*+ Number of local computation threads       +*/
  INT                       bublnbr;              /*+ Number of local computation threads       +*/
  BubbleTree  * restrict    btree;                /*+ Bubbles tree                              +*/

  BlockTarget * restrict    btagtab;              /*+ Blocktarget access vector                 +*/
  INT                       btagnbr;              /*+ Number of Blocktargets                    +*/
  INT                       btgsnbr;              /*+ Number of Blocktargets to send            +*/
  INT                       btgrnbr;              /*+ Number of Blocktargets to recv            +*/
  BlockCoeff  * restrict    bcoftab;              /*+ BlockCoeff access vector                  +*/
  INT                       bcofnbr;

  INT                       indnbr;
  INT * restrict            indtab;
  Task * restrict           tasktab;              /*+ Task access vector                        +*/
  INT                       tasknbr;              /*+ Number of Tasks                           +*/
  INT **                    ttsktab;              /*+ Task access vector by thread              +*/
  INT *                     ttsknbr;              /*+ Number of tasks by thread                 +*/

  INT *                     proc2clust;           /*+ proc -> cluster                           +*/
  INT                       gridldim;             /*+ Dimensions of the virtual processors      +*/
  INT                       gridcdim;             /*+ grid if dense end block                   +*/
  UpDownVector              updovct;              /*+ UpDown vector                             +*/
} SolverMatrix;

#endif /* SOLVER_H */
