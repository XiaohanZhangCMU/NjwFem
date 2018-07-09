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
#ifndef SOPALIN_INIT_H
#define SOPALIN_INIT_H

/* Type of thread asking for allocation (must be a power of 2) */
#define INIT_COMPUTE   1
#define INIT_SEND      2
#define INIT_RECV      4
#define INIT_OOC       8

/************************************************/
/*     Structure pour le backup des donn�es     */
/************************************************/

typedef struct Backup_ {
  INT              cpftmax;         /* Double by LU version             */
  INT              arftmax;         /* Double by LU version             */
  INT              nbftmax;
  INT *            task_ctrbcnt;    /* no inital value                  */
  INT *            task_ftgtcnt;    /* no inital value                  */
  INT *            fanin_ctrbnbr;   /* change updown information        */
  INT *            fanin_prionum;   /* both used for tag and pack       */
  INT *            btagtaskcnt;     /* btag taskcnt                     */
  INT *            bcofsendcnt;     /* bcof sendcnt                     */
  INT *            symbol_cblknum;  /* sopalin add negative information */
  INT              symbol_nodenbr;  /* ???                              */
} Backup;

typedef struct BackupSolve_ {
  INT *            fanin_ctrbnbr;   /* change updown information        */
  INT *            symbol_cblknum;  /* sopalin add negative information */
} BackupSolve_t;


/* Allocate and initialize/Free globale data for solver */
void sopalin_init     (Sopalin_Data_t *sopalin_data, SolverMatrix *m, SopalinParam *sopaparam, int fact);
void sopalin_clean    (Sopalin_Data_t *sopalin_data, int step);

/* Allocate and initialize/Free thread data for solver */
void sopalin_init_smp (Sopalin_Data_t *sopalin_data, INT me, int fact, int thrdtype);
void sopalin_clean_smp(Sopalin_Data_t *sopalin_data, INT me);

/* Restore/backup des donn�es modifi�es pendant l'encha�nement facto/solve */
void sopalin_backup (SolverMatrix *datacode, Backup *b);
void sopalin_restore(SolverMatrix *datacode, Backup *b);

/* Restore/backup des donn�es modifi�es pendant le solve */
void solve_backup (SolverMatrix *datacode, BackupSolve_t *b);
void solve_restore(SolverMatrix *datacode, BackupSolve_t *b);

#endif /* SOPALIN_INIT_H */
