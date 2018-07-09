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
/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : splifunc.h                              +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Compute optimal split                   +*/
/*+                                                        +*/
/*+   DATES      : # Version 0.0  : from : 01 jan 1999     +*/
/*+                                 to     20 jan 1999     +*/
/*+                                                        +*/
/*+********************************************************+*/

#ifndef SPLITFUNC_H
#define static
#endif
INT           splitSeqCblk2D    (INT, INT, const SymbolMatrix *, const ExtraSymbolMatrix *, const Dof *, const BlendCtrl *,
			         INT (*P)(INT , INT, const SymbolMatrix *, const ExtraSymbolMatrix *, const Dof *, INT, const BlendCtrl *));
static void   virtualSplit      (INT, INT, const SymbolBlok *, INT *, SymbolBlok *);
static double  cblkCost          (INT, const SymbolBlok *, const Dof *);
static void   build_cblk        (INT, const SymbolMatrix *, const ExtraSymbolMatrix *, SymbolBlok *);
static INT    cblkNbr           (INT,  const SymbolMatrix *, const ExtraSymbolMatrix *);
INT           P1D               (INT, INT, const SymbolMatrix *, const ExtraSymbolMatrix *const extrasymbptr, const Dof *, INT, const BlendCtrl *);
INT           P2D               (INT, INT, const SymbolMatrix *, const ExtraSymbolMatrix *const extrasymbptr, const Dof *, INT, const BlendCtrl *);
#undef static


