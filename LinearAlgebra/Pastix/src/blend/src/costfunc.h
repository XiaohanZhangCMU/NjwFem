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
/*+   NAME       : costfunc.h                              +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Cost compute functions                  +*/
/*+                                                        +*/
/*+   DATES      : # Version 0.0  : from : 27 sep 1998     +*/
/*+                                 to     03 oct 1998     +*/
/*+                                                        +*/
/*+********************************************************+*/

#ifndef COSTFUNC_H
#define static
#endif


void            costMatrixBuild       (CostMatrix *, const SymbolMatrix *, const Dof *);
void            costMatrixCorrect     (CostMatrix *, const SymbolMatrix *, Cand * candtab,  const Dof *);
double          subtreeUpdateCost     (INT, CostMatrix *, const EliminTree *);
double          subtreeUpdateCostLocal(INT, const BlendCtrl *, const SymbolMatrix *, const SimuCtrl *, const Dof *, INT);
double          cblkComputeCost       (INT, CostMatrix *, const SymbolMatrix *, const Dof *);
double          cblkComputeCost2D     (INT, CostMatrix *, const SymbolMatrix *, const Dof *);
				 		   
/** 2D **/			 	   
double          DIAGCost              (INT);
double          E1Cost                (INT, INT);
double          E2Cost                (INT, INT, INT);
				 		   
static double   computeCost           (INT, INT);
static double   contribCompCost       (INT, INT, INT);
static double   contribAddCost        (INT, INT);
double          costFtgtSend          (INT, INT, FanInTarget *, BlendCtrl *,  const Dof *);
				 		   
double          costFtgtAdd           (FanInTarget *, const Dof *);
double          cblkMaxCost           (INT, const CostMatrix *);
double          totalCost             (INT, const CostMatrix *);
void            printSolverInfo       (FILE *, const SolverMatrix *, const SymbolMatrix *, const Dof * const dofptr);
double          memorySpaceCost       (const SolverMatrix *);
static double   solverSpaceCost       (const SolverMatrix *);
static double   symbolSpaceCost       (const SymbolMatrix *);

#undef static




