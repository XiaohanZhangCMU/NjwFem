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
/**   NAME       : solverMatrixGen.h                       **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Genere the local solver matrix          **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 Oct 1998     **/
/**                                 to     15 Oct 1998     **/
/**                                                        **/
/************************************************************/
#ifndef SOLVERMATRIXGEN_H
#define SOLVERMATRIXGEN_H
INT *               solverMatrixGen (const INT, SolverMatrix *, const SymbolMatrix *, const SimuCtrl *, const BlendCtrl *, const Dof *);
void                allSolverMatrixSave(const char *, const SymbolMatrix *, const SimuCtrl *, const BlendCtrl *, const Dof *);
#endif /* SOLVERMATRIXGEN_H */