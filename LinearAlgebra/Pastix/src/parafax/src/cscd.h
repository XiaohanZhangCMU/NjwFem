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
#ifndef CSCD_H
#define CSCD_H

#define CSCD_IA(x)         (x)->csc->ia
#define CSCD_JA(x)         (x)->csc->ja
#define CSCD_LOC2GLB(x)    (x)->loc2glb
#define CSCD_AVAL(x)       (x)->csc->aval
#define CSCD_RHS(x)        (x)->csc->rhs
#define CSCD_VERTLOC(x)    (x)->csc->vertlocnbr
#define CSCD_EDGELOC(x)    (x)->csc->edgelocnbr


typedef struct CSC_ {
  INT               vertlocnbr;  /*+ Number of local verteces, correspond to ia local vector size */
  INT               edgelocnbr;  /*+ Number of local edges, correpond to ja local vector size +*/
  INT             * ia;          /*+ Table des sommets locaux +*/
  INT             * ja;          /*+ Table des sommets au bout de l'arete :) (num global) +*/
  FLOAT           * aval;        /*+ Table de poids des aretes locales +*/
  FLOAT           * rhs;         /*+ Right Hand Side +*/
} CSC;

typedef struct CSCD_ {
  CSC             * csc;      
#ifdef DISTRIBUTED
  SCOTCH_Dgraph   * dgrafptr;    
#endif
  INT             * loc2glb;     /*+ Table of global indices for local vertices +*/
  MPI_Comm          mpi_comm;
  INT               j_size;
  INT               j_cst;
} CSCD;

/* Alloue, la structure, les champs et effectue les recopies memoires */
double second              ();

int    csc_Init            (CSC  *, INT,     INT,   INT *, INT *,   FLOAT *, FLOAT *);
int    csc_Free            (CSC  *);
int    cscd_Init           (CSCD *, INT *,   INT,   INT,   INT *,   INT *,   FLOAT *, FLOAT *, MPI_Comm);
int    cscd_Free           (CSCD *);

int    cscd_Save           (CSCD  *, const char *);    /* Enregistre la CSCd                                 */
int    cscd_Load           (CSCD  *, const char *, MPI_Comm ); 

int    cscd_Sym            (CSCD  *, INT);             /* Symetrise la CSCd                                  */
int    cscd_Addedge        (CSCD  *, INT, INT, FLOAT); /* Ajoute une arete a une CSCd                        */
int    cscd_Convert_dgraph (CSCD  *, INT);             /* Prepare et effectue l'appel a dgraph_build         */
int    cscd_Fusion         (CSCD  *);                  /* Concentre la CSCd sur le noeud maitre 0, pas utile */
int    cscd_Explode        (CSCD  *, MPI_Comm);        /* Divise la CSC sur le comm MPI                      */
int    cscd_Remvertices    (CSCD  *, CSCD *, INT *, INT);
int    cscd_Addvertices    (CSCD  *, INT  *, INT *, INT *, FLOAT *, FLOAT *, INT, INT);
int    cscd_Addvalues      (CSCD  *, CSCD *);

INT    cscd_getGloballocal (CSCD  *, INT);
INT    cscd_indIfLoc       (CSCD  *, INT);          /* donne l'indice local si le sommet est local sinon -1 */
int    cscd_add_diag       (CSCD  *);
int    cscd_NoDiag         (CSCD * cscd_graph) ;
int    cscd_sort           (CSCD  *);
void   cscd_quicksort      (INT   *, FLOAT *, INT,   INT);
int    cscd_partition      (INT   *, FLOAT *, INT,   INT,   INT);


/*
  Function: cscd_Redispatch

  Renumber the columns to have first columns on first proc for Scotch

  Parameters:
    cscd_graph - the CSC

  Returns:
    EXIT_SUCCESS if already well distributed, 2 if redistributed
    
 */
int cscd_Redispatch(CSCD * cscd_graph);

#endif /* CSCD_H */
