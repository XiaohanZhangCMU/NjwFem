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
/**   NAME       : updown.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the UpDown step  .                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     28 oct 1998     **/
/**                                                        **/
/************************************************************/

#ifndef UPDOWN_H
#define UPDOWN_H

/*+ UpDown block structure. +*/

typedef struct UpDownCblk_  {
  INT                       sm2xind;              /*+ Index in the rhs local vector of the unknowns corresponding to the diag blok +*/
  INT *                     browproctab;          /*+ Brow                               +*/
  INT *                     browcblktab;          /*+ Brow                               +*/
  INT                       browprocnbr;          /*+ Brow size                          +*/
  INT                       msgnbr;               /*+ Number of messages                 +*/
  INT volatile              msgcnt;               /*+ Number of messages                 +*/
  INT                       ctrbnbr;              /*+ Number of contributions            +*/
  INT volatile              ctrbcnt;              /*+ Number of contributions            +*/
} UpDownCblk;


/*+ UpDown vector structure. +*/

typedef struct UpDownVector_ {
  UpDownCblk *              cblktab;              /*+ Array of solver column blocks      +*/
  FLOAT *                   sm2xtab;              /*+ Unknown vector                     +*/
  INT                       sm2xmax;              /*+ Maximum of coefficients per unknown vector +*/
  INT                       sm2xsze;              /*+ Size of sm2xtab                    +*/
  INT                       sm2xnbr;              /*+ Number of sm2x                     +*/
  INT *                     gcblk2list;           /*+ Global cblknum -> index in listptr +*/
  INT                       gcblk2listnbr;        /*+ Size of gcblk2list                 +*/
  INT *                     listptr;              /*+ Index in list                      +*/
  INT                       listptrnbr;           /*+ Size of listptr                    +*/
  INT *                     listcblk;             /*+ List of cblk in a same row         +*/
  INT *                     listblok;             /*+ List of blok in a same row         +*/
  INT                       listnbr;              /*+ Size of list                       +*/
  INT *                     loc2glob;             /*+ Local cblknum -> global cblknum    +*/
  INT                       loc2globnbr;          /*+ Size of loc2glob                   +*/
  INT *                     lblk2gcblk;           /*+ Local blok -> global facing cblk   +*/
  INT                       gcblknbr;             /*+ total number of cblk               +*/
  INT                       gnodenbr;             /*+ total number of nodes              +*/
  INT                       downmsgnbr;           /*+ Nb messages receive during down    +*/
  INT                       upmsgnbr;             /*+ Nb messages receive during up      +*/
} UpDownVector;

#endif /* UPDOWN_H */
