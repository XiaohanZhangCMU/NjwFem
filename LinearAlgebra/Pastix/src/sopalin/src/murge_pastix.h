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
/******************************************************************************
 * Title: PaStiX specific addons to Murge                                     *
 *                                                                            *
 * Functions declaration for Murge function defined only in PaStiX.           *
 *                                                                            *
 * Authors:                                                                   *
 *   Xavier Lacoste - xavier.lacoste@inria.fr                                 *
 *                                                                            *
 * More informations can be found in <murge.h> documentation.                 *
 *                                                                            *
 ******************************************************************************/


/******************************************************************************
 * Function: MURGE_Analyze                                                    *
 *                                                                            *
 * Perform matrix analyze.                                                    *
 *                                                                            *
 * Follows several steps :                                                    *
 *   - Compute a new ordering of the unknows                                  *
 *   - Compute the symbolic factorisation of the matrix                       *
 *   - Distribute column blocks and computation on processors                 *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Analyze(INTS id);

/******************************************************************************
 * Function: MURGE_Factorize                                                  *
 *                                                                            *
 * Perform matrix factorization.                                              *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Factorize(INTS id);

/******************************************************************************
 * Function: MURGE_SetOrdering                                                *
 *                                                                            *
 * Set a personal ordering to perform factorization.                          *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   permutation - Permutation to perform factorization.                      *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_SetOrdering(INTS   id,
                       INTS * permutation);

/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeNbr                                     *
 *                                                                            *
 * Set local node number if users only perform product and doesn't need       *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   n  - Number of local nodes.                                              *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ProductSetLocalNodeNbr (INTS id, INTS n);

/******************************************************************************
 * Function: MURGE_ProductSetGlobalNodeNbr                                    *
 *                                                                            *
 * Set global node number if users only perform product and doesn't need      *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   N  - Number of global nodes.                                             *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ProductSetGlobalNodeNbr (INTS id, INTS N);

/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeList                                    *
 *                                                                            *
 * Set local node list if users only perform product and doesn't need         *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id  - Solver instance identification number.                             *
 *   l2g - Local to global node numbers.                                      *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ProductSetLocalNodeList (INTS id, INTS * l2g);

/******************************************************************************
 * Function: MURGE_GetLocalProduct                                            *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <MURGE_SetLocalRHS> or              *
 * <MURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   x  - Array in which the local part of the product will be stored.        *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_GetLocalProduct (INTS id, COEF *x);

/******************************************************************************
 * Function: MURGE_GetGlobalProduct                                           *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <MURGE_SetLocalRHS> or              *
 * <MURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id   - Solver instance identification number.                            *
 *   x    - Array in which the product will be stored.                        *
 *   root - Rank of the process which will own the product at end of call,    *
 *          use -1 for all processes.                                         *
 * Returns:                                                                   *
 *   MURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_GetGlobalProduct (INTS id, COEF *x, INTS root);

/******************************************************************************
 * Function: MURGE_ForceNoFacto                                               *
 *                                                                            *
 * Prevent Murge from running factorisation even if matrix has changed.       *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 * Returns:                                                                   *
 *   MURGE_SUCCESS                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ForceNoFacto(INTS id);

/******************************************************************************
 * Function: MURGE_SetLocalNodeList                                           *
 *                                                                            *
 * Deprecated, need to be checked                                             *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_SetLocalNodeList(INTS id, INTS n, INTS * list);

/******************************************************************************
 * Function: MURGE_AssemblySetSequence                                        *
 *                                                                            *
 * Create a sequence of entries to build a matrix and store it for being      *
 * reused.                                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id      - Solver instance identification number.                         *
 *   coefnbr - Number of entries.                                             *
 *   ROWs    - List of rows in the sequence.                                  *
 *   COLs    - List of columns in the sequence.                               *
 *   op      - Operation to perform for coefficient which appear              *
 *             several tim (see <MURGE_ASSEMBLY_OP>).                         *
 *   op2     - Operation to perform when a coefficient is set by              *
 *             two different processors (see <MURGE_ASSEMBLY_OP>).            *
 *   mode    - Indicates if user ensure he will respect solvers distribution  *
 *             (see <MURGE_ASSEMBLY_MODE>).                                   *
 *   nodes   - 0 entries are entered value by value,                          *
 *             1 entries are entries node by node.                            *
 *   id_seq  - Sequence ID.                                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.   *
 ******************************************************************************/
int MURGE_AssemblySetSequence (INTS id, INTL coefnbr, INTS * ROWs, INTS * COLs,
                               INTS op, INTS op2, INTS mode, INTS nodes,
                               INTS * id_seq);

/******************************************************************************
 * MURGE_AssemblySetSequence                                                  *
 *                                                                            *
 * Assembly the matrix using a stored sequence.                               *
 *                                                                            *
 * Parameters:                                                                *
 *   id      - Solver instance identification number.                         *
 *   id_seq  - Sequence ID.                                                   *
 *   values  - Values to insert in the CSC.                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *id_seq* or *values* are not valid.                *
 ******************************************************************************/
INTS MURGE_AssemblyUseSequence(INTS id, INTS id_seq, COEF * values);

/******************************************************************************
 * Function: MURGE_AssemblyDeleteSequence                                     *
 *                                                                            *
 * Destroy an assembly sequence                                               *
 *                                                                            *
 *   id      - Solver instance identification number.                         *
 *   id_seq  - Sequence ID.                                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *id_seq* is not valid.                             *
 ******************************************************************************/
INTS MURGE_AssemblyDeleteSequence(INTS id, INTS id_seq);
