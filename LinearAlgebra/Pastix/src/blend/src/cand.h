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
/**   NAME       : cand.h                                  **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the candidate group of a cblk       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 sep 1998     **/
/**                                 to     08 oct 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/
#define D1 0
#define D2 1
#define DENSE 3

/*+ Processor candidate group to own a column blok      +*/
typedef struct Cand_{
  INT treelevel;    /*+ Level of the cblk in the elimination tree (deepness from the root) +*/
  double costlevel; /*+ Cost from root to node +*/
  INT fcandnum;     /*+ first processor number of this candidate group  +*/
  INT lcandnum;     /*+ last processor number of this candidate group   +*/
  INT fccandnum;    /*+ first cluster number of the cluster candidate group +*/
  INT lccandnum;    /*+ last cluster number of the cluster candidate group +*/
  INT distrib;      /*+ type of the distribution +*/
  INT cluster;      /*+ TRUE if cand are clusters (bubble number) +*/
#if defined(TRACE_SOPALIN) || defined(PASTIX_DYNSCHED)
  INT cand;         /*+ TRUE if cand are clusters (bubble number) +*/
#endif
} Cand;

