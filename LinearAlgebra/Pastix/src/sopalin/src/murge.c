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
/*
 * File: murge.c
 *
 * This file implements <Murge> interface.
 *
 * About: Authors
 *   Mathieu Faverge - faverge@labri.fr
 *   Xavier Lacoste  - xavier.lacoste@inria.fr
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#ifndef FORCE_NOSMP
#  include <pthread.h>
#endif
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else /* not FORCE_NOMPI */
#  include <mpi.h>
#endif /* not FORCE_NOMPI */
#ifdef WITH_SEM_BARRIER
#  include <semaphore.h>
#endif
#include "common_pastix.h"
#include "tools.h"
#include "sopalin_define.h"

#ifdef WITH_SCOTCH
#  ifdef    DISTRIBUTED
#    include "ptscotch.h"
#  else
#include "scotch.h"
#endif /* DISTRIBUTED */
#endif /* WITH_SCOTCH */

#include "ftgt.h"
#include "symbol.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "order.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"
#include "pastixstr.h"
#include "pastix.h"
#include "cscd_utils.h"
#include "cscd_utils_intern.h"
#include "csc_intern_compute.h"
#include "csc_intern_updown.h"
#include "sopalin_init.h"
#include "sopalin_compute.h"

#include "murge.h"

#include "murge_pastix.h"
#include "murge_defines.h"
#ifdef DISTRIBUTED
/******************************************************************************/
/***                           Section: Structures                          ***/
/******************************************************************************/

/*
 * Structure: ijv_
 *
 * Structure to represente coefficients.
 *
 * Contains:
 *   i     - row
 *   j     - column
 *   v     - pointer to the value array (can be several degree of freedom)
 *   owner - process which own the coefficient.
 */
struct ijv_ {
  INT    i;
  INT    j;
  FLOAT* v;
#ifdef MURGE_FOLLOW_INPUT_ORDER
  INT    idx;
#endif
  int    owner;
};

/*
 * Typedef: ijv_t
 *
 * Alias to structure <ijv_>.
 */
typedef struct ijv_ ijv_t;

/*
 * struct: murge_seq_t
 *
 * Structure used to store assembly sequence.
 *
 * Contains:
 *   indexes              - Order sequences of indexes to be set.
 *   coefnbr              - Number of entries in the sequence.
 *   recv_nbr             - Number of entries to receive from each processor.
 *   recv_indexes         - Indexes of entries which will be received.
 *   mode                 - Local entries or communicating mode
 *   fusion_local_entries - Operation to perform when a coefficient appear twice
 *   fusion_dist_entries  - Operation to perform when a coefficient appear
 *                            twice, given by two processors.
 *   ijv_size             - size of the required array to store
 *                            not local entries.
 *   nodes                - 0 entries are entered value by value,
 *                          1 entries are entries node by node.
 *   next                 - Next entry in the list of sequences.
 */
typedef struct murge_seq_ murge_seq_t;
struct murge_seq_ {
  INTL        * indexes;
  INTL          coefnbr;
  INTL        * recv_nbr;
  INTL       ** recv_indexes;
  INTS          mode;
  FLOAT       (*fusion_local_entries)(FLOAT , FLOAT);
  FLOAT       (*fusion_dist_entries)(FLOAT , FLOAT);
  INTL          ijv_size;
  INTS          nodes;
  murge_seq_t * next;
  INT           ID;
};


/*
 * struct: murge_data_t
 *
 * Structure used to store murge data
 *
 * Contains:
 *   pastix_data - Pointer to the <pastix_data_t> associated to
 *                 the solver instance
 *   n           - Number of local column indicated by murge user
 *   N           - Number of global column indicated by murge user
 *   colptr      - Colptr in murge's user CSCd
 *   rows        - Rows in murge's user CSCd
 *   values      - Values in murge's user CSCd
 *   l2g         - Local to global column number in murge's user CSCd
 *   g2l         - Global to local column number in murge's user CSCd
 *   perm        - Permtab for murge's user
 *   b           - Right-hand-side member(s) given by murge's user
 *   nrhs        - Number of right-hand-side member(s) given by murge's user
 *   tmpv        - Temporary values array
 *   tmpv_node   - Temporary values array (node entries)
 *   tmpijv      - Temporary ijv structure array
 *   cnt         - Iterator for number of entered edges
 *   edgenbr     - Number of edges
 *   state       - State of the solver
 *   mode        - Local entries or communicating mode
 *   op          - Operation to perform when a coefficient appear twice
 *   op2         - Operation to perform when a coefficient appear twice,
 *                 given by two processors.
 *   sym         - Indicate if we have to check that the matrix is symmetric
 *   sequences   - List of saved sequences.
 */
struct murge_data_{
  pastix_data_t   *pastix_data;
  INT              n;
  INT              N;
  INT             *colptr;
  INT             *rows;
  FLOAT           *values;
  INT             *l2g;
  INT             *g2l;
  INT             *perm;
#ifdef CENTRALISED
  INT             *invp;
#endif
  FLOAT           *b;
  INT              nrhs;
  FLOAT           *tmpv;
  FLOAT           *tmpv_node;
  ijv_t           *tmpijv;
  ijv_t           *tmpijv_node;
  INT              tmpijv_size;
  INT              tmpijv_size_node;
#ifdef MURGE_THREADSAFE
  pthread_mutex_t  mutex_tmpmatrix;
#endif
  INT              cnt;
  INT              cnt_zero;
  INT              cnt_node;
  INT              edgenbr;
  INT              coefnbr;
  INT              nodenbr;
  int              state;
  int              mode;
  int              op;
  int              op2;
  int              sym;
  int              dynamic;
  murge_seq_t     *sequences;
  INT              seq_ID;
  int              ndump;
};

/*
 * Typedef: murge_data_t
 *
 * alias to structure <murge_data_>.
 */
typedef struct murge_data_ murge_data_t;


/******************************************************************************/
/***                           Section: Global variables                    ***/
/******************************************************************************/


/*
 * Variables: Global variables
 *
 * idnbr   - Number of solvers instances.
 * solvers - Murge solver instances array (<murge_data_t>).
 */
INTS           idnbr   = 0;
murge_data_t **solvers = NULL;


/******************************************************************************/
/***                           Section: Functions                           ***/
/******************************************************************************/


/*******************************************************************************
 * Group: Auxilary functions
 */

#ifdef CENTRALISED
#define ALLOC_INVP                                              \
  if (NULL == solvers[id]->invp)                                \
    MURGE_MEMALLOC(solvers[id]->invp, solvers[id]->n, INT);
#define INVP solvers[id]->invp
#else
#define ALLOC_INVP {};
#define INVP NULL
#endif
/*
 * Function: check_preprocessing
 *
 * Checks if preprocessing (blend) has been called.
 *
 * If it hasn't, it will allocate permutation tabular
 * and call preprocessing step.
 *
 * After calling preprocessing, it will set local number of column
 * and local to global column number tabular to their new values.
 *
 * Colptr and rows will be destroyed because it is obsolete,
 * and state will be set to indicate that preprocessing has been performed.
 *
 * Parameters:
 *   id - Solver instance ID we want to check
 *
 * Returns:
 *   MURGE_ERR_ALLOCATE - If any allocation error occurs.
 */
static inline
int check_preprocessing(int id) {
  pastix_data_t   *pastix_data = solvers[id]->pastix_data;
  INT             *iparm       = pastix_data->iparm;

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_BLEND_OK)))
    {
      /* Si il n'est pas fait on effectue le pretraitement */
      if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_SYMB_OK)))
        iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
      else
        iparm[IPARM_START_TASK]  = API_TASK_BLEND;
      iparm[IPARM_END_TASK]    = API_TASK_BLEND;
      pastix_welcome_print(solvers[id]->pastix_data,
                           solvers[id]->colptr,
                           solvers[id]->n);

      if (solvers[id]->perm == NULL)
        MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, INT);

      ALLOC_INVP;

      if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES)
        {
          INT err;
          if (NO_ERR !=
              ((err =
                pastix_checkMatrix_int(pastix_data->pastix_comm,
                                       iparm[IPARM_VERBOSE],
                                       iparm[IPARM_SYM],
                                       API_YES,
                                       solvers[id]->n,
                                       &(solvers[id]->colptr),
                                       &(solvers[id]->rows),
                                       NULL,
                                       &(solvers[id]->l2g),
                                       iparm[IPARM_DOF_NBR],
                                       API_YES))))
            {
              errorPrint("pastix_checkMatrix : err %ld\n", (long)err);
              return MURGE_ERR_PARAMETER;
            }
        }

      dpastix(&(pastix_data),
              pastix_data->pastix_comm,
              solvers[id]->n,
              solvers[id]->colptr,
              solvers[id]->rows,
              solvers[id]->values,
              solvers[id]->l2g,
              solvers[id]->perm,
              INVP,
              solvers[id]->b,
              solvers[id]->nrhs,
              iparm,
              pastix_data->dparm);
#ifdef MURGE_INSERT_DIRECTLY
      /* We build the new CSC which will receive the coefficients */
      {
        INT newn;
        INT * newl2g;
        INT * newcolptr;
        INT * newrows;

        newn = pastix_getLocalNodeNbr(&(pastix_data));
        memFree_null(solvers[id]->perm);
        MURGE_MEMALLOC(solvers[id]->perm, newn, INT);
        MURGE_MEMALLOC(newl2g,            newn, INT);
        pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                               newl2g);

        cscd_redispatch_int(solvers[id]->n,  solvers[id]->colptr,
                            solvers[id]->rows,
                            NULL, NULL, 0,  solvers[id]->l2g,
                            newn,           &newcolptr,           &newrows,
                            NULL, NULL, newl2g, API_YES,
                            solvers[id]->pastix_data->pastix_comm);

        memFree_null(solvers[id]->l2g);
        memFree_null(solvers[id]->colptr);
        memFree_null(solvers[id]->rows);

        solvers[id]->n      = newn;
        solvers[id]->N      = -1;
        solvers[id]->l2g    = newl2g;
        solvers[id]->colptr = newcolptr;
        solvers[id]->rows   = newrows;
        solvers[id]->values = NULL;
      }
#else /* not MURGE_INSERT_DIRECTLY */
      /* On corrige n et l2g */

      solvers[id]->n = pastix_getLocalNodeNbr(&(pastix_data));
      memFree_null(solvers[id]->l2g);
      memFree_null(solvers[id]->perm);
      memFree_null(solvers[id]->colptr);
      memFree_null(solvers[id]->rows);
      MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, INT);
      MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, INT);
      pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                             solvers[id]->l2g);
      solvers[id]->N=-1;

#endif /* not MURGE_INSERT_DIRECTLY */
      /*
       * Building global to local column number array
       */
      cscd_build_g2l(solvers[id]->n,
                     solvers[id]->l2g,
                     solvers[id]->pastix_data->pastix_comm,
                     &solvers[id]->N,
                     &solvers[id]->g2l);

      MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
      MURGE_STATE_TRUE(solvers[id]->state, MURGE_SYMB_OK);
    }

  return MURGE_SUCCESS;
}

static inline
int extend_node_list(INT   *  edgenbr_recv_node,
                     INT      procnum,
                     INT   *  tmpijvsize_node,
                     ijv_t ** tmpijv_node,
                     ijv_t ** ijvptr,
                     int      dof,
                     FLOAT ** tmpvalues,
                     FLOAT ** valuesptr){
  INT nodenbr = edgenbr_recv_node[procnum];
  INT newsize;
  /* Extend the nodelist if needed */
  if (*tmpijvsize_node < nodenbr+1)
    {
      INT iterator;
      ijv_t * tmp  = NULL;
      FLOAT * tmpv = NULL;
      newsize = *tmpijvsize_node + *tmpijvsize_node/2 + 1;

      MURGE_MEMALLOC(tmp, newsize, ijv_t);
      memset(tmp, 0, newsize*sizeof(ijv_t));
      if (nodenbr>0)
        {
          memcpy(tmp,
                 *tmpijv_node,
                 nodenbr*sizeof(ijv_t));

          memFree_null(*tmpijv_node);
        }
      *tmpijv_node = tmp;
      *ijvptr = tmp;
      MURGE_MEMALLOC(tmpv, newsize*dof*dof,
                     FLOAT);
      memset(tmpv, 0, (newsize*dof*dof) *
             sizeof(FLOAT));

      if (nodenbr>0)
        {
          memcpy(tmpv, *tmpvalues,
                 nodenbr*dof*dof*sizeof(FLOAT));
          memFree_null(*tmpvalues);
        }
      *tmpvalues = tmpv;
      *valuesptr = tmpv;

      *tmpijvsize_node = newsize;
      for (iterator = 0; iterator < *tmpijvsize_node; iterator++)
        {
          (*tmpijv_node)[iterator].v = &(tmpv[dof*dof*iterator]);
        }
    }
  return MURGE_SUCCESS;
}

/*
 * Function: cmp_ijv
 *
 * Compare to <ijv_t> structures on their column value (j) then, if equal,
 * their row value (i).
 *
 * Used for qsort.
 *
 * Parameters:
 *   p1 - pointer to the first element.
 *   p2 - pointer to the second element.
 *
 * Returns:
 *   A positive number - if p1.j > p2.j or p1.j == p2.j and pi.i > pi.j
 *   0                 - if columns and rows are aquals.
 *   A negative number - otherwise.
 */
static int cmp_ijv(const void *p1, const void *p2)
{
  /* Les arguments de cette fonction sont des "pointeurs de
   pointeurs sur des caractères", mais les arguments de
   strcmp(3) sont des "pointeurs sur des caractères", d’où
   le forçage de type et l’utilisation de l’astérisque */
  ijv_t ijv1 = *((ijv_t *)p1);
  ijv_t ijv2 = *((ijv_t *)p2);


  if (ijv1.j != ijv2.j)
    return ijv1.j - ijv2.j;
  else
    return ijv1.i - ijv2.i;
}

/*
 * Function: cmp_ijv_own
 *
 * Compare to <ijv_t> structures on their owner, column value (j) then, if equal,
 * their row value (i).
 *
 * Used for qsort.
 *
 * Parameters:
 *   p1 - pointer to the first element.
 *   p2 - pointer to the second element.
 *
 * Returns:
 *   A positive number - if p1.owner > p2.owner or
 *                       p1.j > p2.j or p1.j == p2.j and pi.i > pi.j
 *   0                 - if columns and rows are aquals.
 *   A negative number - otherwise.
 */
static int cmp_ijv_own(const void *p1, const void *p2)
{
  /* Les arguments de cette fonction sont des "pointeurs de
   pointeurs sur des caractères", mais les arguments de
   strcmp(3) sont des "pointeurs sur des caractères", d’où
   le forçage de type et l’utilisation de l’astérisque */
  ijv_t ijv1 = *((ijv_t *)p1);
  ijv_t ijv2 = *((ijv_t *)p2);

  if (ijv1.owner != ijv2.owner)
    return ijv1.owner - ijv2.owner;

  if (ijv1.j != ijv2.j)
    return ijv1.j - ijv2.j;

  return ijv1.i - ijv2.i;
}


/*
 * Function: MurgeTmpijvOwnSort
 *
 * Sort ijv_t structure, sorting also associated values.
 * Sort following owner, then j and then i.
 *
 * Parameters:
 *   pbase       - Array of two pointers, one  to the first element
 *                 of the array to sort, and the second to dofnbr.
 *   total_elems - Number of element in the array.
 *
 * Returns:
 *   Nothing
 */
#define INTSORTNAME            MurgeTmpijvOwnSort
#define INTSORTSIZE(x)         (sizeof (ijv_t))
#define INTSORTNTAB            1
#define OWNER(p)               ((ijv_t *) (p))->owner
#define IJV_T_J(p)             ((ijv_t *) (p))->j
#define IJV_T_I(p)             ((ijv_t *) (p))->i
#define IJV_T_INDEX(p)         ((ijv_t *) (p))->idx
#define INTSORTSWAP(p,q)       do {                             \
    ijv_t   t;                                                  \
    FLOAT   f;                                                  \
    int     dof = *((int *)(*(pbase+1)));                       \
    INTS    iter;                                               \
    void   *ptr;                                                \
    /* swap integers */                                         \
    t = *((ijv_t *) (p));                                       \
    *((ijv_t *) (p)) = *((ijv_t *) (q));                        \
    *((ijv_t *) (q)) = t;                                       \
    /* swap values */                                           \
    for (iter = 0; iter < dof*dof; iter++)                      \
      {                                                         \
        f = ((ijv_t *) (p))->v[iter];                           \
        ((ijv_t *) (p))->v[iter] = ((ijv_t *) (q))->v[iter];    \
        ((ijv_t *) (q))->v[iter] = f;                           \
      }                                                         \
    ptr = ((ijv_t *) (p))->v;                                   \
    ((ijv_t *) (p))->v = ((ijv_t *) (q))->v;                    \
    ((ijv_t *) (q))->v = ptr;                                   \
  } while (0)
#ifdef MURGE_FOLLOW_INPUT_ORDER
#define INTSORTCMP(p,q)  ((OWNER(p) < OWNER(q) ||                       \
                                               (OWNER(p) == OWNER(q) && \
                                                  (IJV_T_J(p) < IJV_T_J(q) || \
                                                                           (IJV_T_J(p) == IJV_T_J(q) && \
                                                                              (IJV_T_I(p) < IJV_T_I(q) || \
                                                                                                       (IJV_T_I(p) == IJV_T_I(q) && \
                                                                                                          IJV_T_INDEX(p) < IJV_T_INDEX(q))))))))
#else
#define INTSORTCMP(p,q)  ((OWNER(p) < OWNER(q) ||               \
                           (OWNER(p) == OWNER(q) &&             \
                            (IJV_T_J(p) < IJV_T_J(q) ||         \
                             (IJV_T_J(p) == IJV_T_J(q) &&       \
                              (IJV_T_I(p) < IJV_T_I(q)))))))
#endif
#include "../../common/src/common_sort3.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#undef OWNER
#undef IJV_T_J
#undef IJV_T_I

/*
 * Function: MurgeTmpijvSort
 *
 * Sort ijv_t structure, sorting also associated values.
 * Sort following j and then i.
 *
 * Parameters:
 *   pbase       - Array of two pointers, one  to the first element
 *                 of the array to sort, and the second to dofnbr.
 *   total_elems - Number of element in the array.
 *
 * Returns:
 *   Nothing
 */
#define INTSORTNAME            MurgeTmpijvSort
#define INTSORTSIZE(x)         (sizeof (ijv_t))
#define INTSORTNTAB            1
#define OWNER(p)               ((ijv_t *) (p))->owner
#define IJV_T_J(p)             ((ijv_t *) (p))->j
#define IJV_T_I(p)             ((ijv_t *) (p))->i
#define INTSORTSWAP(p,q)       do {                             \
    ijv_t   t;                                                  \
    FLOAT   f;                                                  \
    int     dof = *((int *)(*(pbase+1)));                       \
    INTS    iter;                                               \
    void   *ptr;                                                \
    /* swap integers */                                         \
    t = *((ijv_t *) (p));                                       \
    *((ijv_t *) (p)) = *((ijv_t *) (q));                        \
    *((ijv_t *) (q)) = t;                                       \
    /* swap values */                                           \
    for (iter = 0; iter < dof*dof; iter++)                      \
      {                                                         \
        f = ((ijv_t *) (p))->v[iter];                           \
        ((ijv_t *) (p))->v[iter] = ((ijv_t *) (q))->v[iter];    \
        ((ijv_t *) (q))->v[iter] = f;                           \
      }                                                         \
    ptr = ((ijv_t *) (p))->v;                                   \
    ((ijv_t *) (p))->v = ((ijv_t *) (q))->v;                    \
    ((ijv_t *) (q))->v = ptr;                                   \
  } while (0)
#define INTSORTCMP(p,q)  ((IJV_T_J(p) < IJV_T_J(q) ||                   \
                                                   (IJV_T_J(p) == IJV_T_J(q) && \
                                                      IJV_T_I(p) < IJV_T_I(q))))
#include "../../common/src/common_sort3.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#undef OWNER
#undef IJV_T_J
#undef IJV_T_I
#endif /* DISTRIBUTED */
/*******************************************************************************
 * Group: Solver setup functions
 */

/*
 * Function: MURGE_GetSolver
 *
 * returns MURGE_SOLVER_PASTIX
 */
INTS MURGE_GetSolver(INTS * solver_id)
{
#ifdef DISTRIBUTED
  *solver_id = MURGE_SOLVER_PASTIX;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
#endif /* DISTRIBUTED */
  return MURGE_SUCCESS;
}

/*
 * Function: MURGE_Initialize
 *
 * Allocate the instance arrays which will keeps intern data for all
 * solver instances.
 *
 * If user is creating several threads calling the solver, this function
 * has to be called before creating threads to insure solver is thread safe.
 *
 * Parameters:
 *   idnbr - Maximum number of solver instances that will be
 *           launched.
 *
 * Returns:
 *   MURGE_SUCCESS      - If function runned successfully.
 *   MURGE_ERR_ALLOCATE - If for some reason, allocation was not
 *                        successfull.
 */
INTS MURGE_Initialize(INTS id_nbr){
#ifdef DISTRIBUTED
  INTS i;

  print_debug(DBG_MURGE, ">> Murge_Initialize\n");

  if (sizeof(COEF) != sizeof(FLOAT))
    {
      errorPrint("Incompatible coefficient type\n");
      return MURGE_ERR_PARAMETER;
    }

  if ( (solvers != NULL) )
    {
      errorPrint("MURGE_Initialize has been already called");
      return MURGE_ERR_ORDER;
    }

  idnbr = id_nbr;

  solvers = (murge_data_t**)malloc(idnbr*sizeof(murge_data_t*));

  for (i=0; i< idnbr; i++)
    {
      solvers[i] = NULL;
      solvers[i] = (murge_data_t*)malloc(sizeof(murge_data_t));

      solvers[i]->n           = 0;
      solvers[i]->N           = 0;
      solvers[i]->colptr      = NULL;
      solvers[i]->rows        = NULL;
      solvers[i]->values      = NULL;
      solvers[i]->l2g         = NULL;
      solvers[i]->g2l         = NULL;
      solvers[i]->perm        = NULL;
#ifdef CENTRALISED
      solvers[i]->invp        = NULL;
#endif
      solvers[i]->b           = NULL;
      solvers[i]->nrhs        = 1;
      solvers[i]->state       = MURGE_INIT_OK;
      solvers[i]->pastix_data = NULL;
      solvers[i]->tmpv        = NULL;
      solvers[i]->tmpv_node   = NULL;
      solvers[i]->tmpijv      = NULL;
      solvers[i]->tmpijv_node = NULL;
      solvers[i]->sequences   = NULL;
      solvers[i]->seq_ID      = 0;
      solvers[i]->ndump       = 0;

      pastix_task_init(&(solvers[i]->pastix_data), MPI_COMM_WORLD, NULL, NULL);
    }

  print_debug(DBG_MURGE, "<< Murge_Initialize\n");
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_SetDefaultOptions
 *
 * Create a solver instance if not created yet.
 *
 * Sets default options, for solver instance number *id*.
 *
 * The default option set correspond to *stratnum* strategy ID,
 * depending on the solver.
 *
 * Needs <MURGE_Initialize> to be called before
 * to allocate solver instances array.
 *
 * Parameters:
 *   id       - Solver instance identification number.
 *   stratnum - Strategy for the default option Set.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_Initialize> was not called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *stratnum* is not valid.
 *   MURGE_ERR_ALLOCATE  - If couldn't create solver instance.
 */
INTS MURGE_SetDefaultOptions(INTS id, INTS stratnum)
{
#ifdef DISTRIBUTED
  print_debug(DBG_MURGE, ">> MURGE_SetDefaultOptions\n");
  CHECK_SOLVER_ID(id);

  solvers[id]->pastix_data->iparm = (INT*)malloc(IPARM_SIZE*sizeof(INT));
  solvers[id]->pastix_data->dparm = (double*)malloc(DPARM_SIZE*sizeof(double));

  pastix_initParam(solvers[id]->pastix_data->iparm,
                   solvers[id]->pastix_data->dparm);

  solvers[id]->pastix_data->iparm[IPARM_PID] = solvers[id]->pastix_data->pastix_id;

#ifdef MURGE_THREADSAFE
  pthread_mutex_init(&solvers[id]->mutex_tmpmatrix, NULL);
#endif
  solvers[id]->pastix_data->iparm[IPARM_RHS_MAKING] = API_RHS_B;
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_RHS_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/*
 * Function: MURGE_SetOptionINT
 *
 * Sets integer option, indicated by *number*, to *value* for the
 * solver instance number *id*.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   number - Identification of the integer parameter.
 *   value  - value to set the parameter to.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *                         *number* or *value* are not valid.
 *
 */
INTS MURGE_SetOptionINT (INTS id, INTS number, INTS value){
#ifdef DISTRIBUTED
  INT murge_param[64];
  INT * iparm = NULL;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  iparm = solvers[id]->pastix_data->iparm;

  murge_param[MURGE_IPARAM_BASEVAL       - 1024] =  IPARM_BASEVAL;
  murge_param[MURGE_IPARAM_DOF           - 1024] =  IPARM_DOF_NBR;
  murge_param[MURGE_IPARAM_SYM           - 1024] =  IPARM_SYM;

  if (number >= 1024)
    {
      number = murge_param[number-1024];
      if (number == IPARM_SYM)
        {
          if (value == MURGE_BOOLEAN_TRUE)
            {
              value = API_SYM_YES;
              if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
                iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
            }
          else
            {
              if (value == MURGE_BOOLEAN_FALSE)
                {
                  value = API_SYM_NO;
                  if (iparm[IPARM_FACTORIZATION] != API_FACT_LU)
                    iparm[IPARM_FACTORIZATION] = API_FACT_LU;
                }
              else
                {
                  errorPrint("Invalid value");
                  return MURGE_ERR_PARAMETER;
                }
            }

        }
    }

  if (!( number < IPARM_SIZE ))
    {
      errorPrint("number is too big");
      return MURGE_ERR_PARAMETER;
    }

  if (number < 0)
    {
      errorPrint("number is negative");
      return MURGE_ERR_PARAMETER;
    }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  iparm[number] = (INT)value;
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 Function: MURGE_SetOptionREAL

 Sets real option, indicated by *number*, to *value* for the
 solver instance number *id*.

 Needs <MURGE_SetDefaultOption> to be called before to initiate
 solver instance data.

 Parameters:
 id     - Solver instance identification number.
 number - Identification of the integer parameter.
 value  - value to set the parameter to.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 called before.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *number* or *value* are not valid.

 */
INTS MURGE_SetOptionREAL(INTS id, INTS number, REAL value)
{
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  INT murge_param[64];

  murge_param[MURGE_RPARAM_EPSILON_ERROR- 1024] =  DPARM_EPSILON_REFINEMENT;

  if (number >= 1024)
    {
      number = murge_param[number-1024];
    }

  if (!( number < DPARM_SIZE ))
    {
      errorPrint("number is too big");
      return MURGE_ERR_PARAMETER;
    }

  if (number < 0)
    {
      errorPrint("number is negative");
      return MURGE_ERR_PARAMETER;
    }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  solvers[id]->pastix_data->dparm[number] = (double)value;
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_SetCommunicator
 *
 * Sets MPI communicator for the given solver instance.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Musn't be called before <MURGE_SAVE>, <MURGE_LOAD>,
 * <MURGE_GetLocalNodeNbr> nor <MURGE_GetLocalUnknownNbr>
 * because the solver as to be runned with the same MPI
 * communicator all along.
 *
 * If this function is not called, MPI communicator will be
 * *MPI_COMM_WORLD*.
 *
 * This function may not exist if the solver
 * has been compiled without MPI.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   mpicomm - MPI communicator to affect the solver to.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before or if it is called after
 *                         the solver starts its computing tasks.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *                         *number* or *value* are not valid.
 */
INTS MURGE_SetCommunicator(INTS id, MPI_Comm mpicomm)
{
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  solvers[id]->pastix_data->pastix_comm = mpicomm;
  solvers[id]->pastix_data->inter_node_comm = mpicomm;
  solvers[id]->pastix_data->intra_node_comm = MPI_COMM_SELF;
  MPI_Comm_size((solvers[id]->pastix_data)->inter_node_comm,
                &((solvers[id]->pastix_data)->inter_node_procnbr));
  MPI_Comm_rank((solvers[id]->pastix_data)->inter_node_comm,
                &((solvers[id]->pastix_data)->inter_node_procnum));
  MPI_Comm_size((solvers[id]->pastix_data)->intra_node_comm,
                &((solvers[id]->pastix_data)->intra_node_procnbr));
  MPI_Comm_rank((solvers[id]->pastix_data)->intra_node_comm,
                &((solvers[id]->pastix_data)->intra_node_procnum));

  MPI_Comm_size(mpicomm, &(solvers[id]->pastix_data->procnbr));
  MPI_Comm_rank(mpicomm, &(solvers[id]->pastix_data->procnum));
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/*******************************************************************************
 * Group: I/O functions
 */

/*
 * Function: MURGE_Save
 *
 * Runs preprocessing step, if not done yet, and save the result to disk,
 * into *directory*, so that it can be resume using <MURGE_Load>.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   directory - Path to the directory where to save the solver step.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 *   MURGE_ERR_IO        - If file(s) couldn't be writen.
 */
INTS MURGE_Save(INTS id, char* directory)
{
#ifdef DISTRIBUTED
  char * dest    = NULL;
  char * src     = NULL;
  int    procnum;
  FILE * stream  = NULL;

  CHECK_SOLVER_ID(id);
  procnum = (int)solvers[id]->pastix_data->procnum;

  if (solvers[id]->pastix_data->iparm == NULL)
    {
      errorPrint("You need to call MURGE_SetDefaultOptions before");
      return MURGE_ERR_ORDER;
    }
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("You need to set graph before");
      return MURGE_ERR_ORDER;
    }

  solvers[id]->pastix_data->iparm[IPARM_IO_STRATEGY] = API_IO_SAVE;
  solvers[id]->pastix_data->iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
  solvers[id]->pastix_data->iparm[IPARM_END_TASK]    = API_TASK_SYMBFACT;

  if (NULL == solvers[id]->perm)
    MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, INT);
  pastix_welcome_print(solvers[id]->pastix_data,
                       solvers[id]->colptr,
                       solvers[id]->n);
  dpastix(&(solvers[id]->pastix_data),
          solvers[id]->pastix_data->pastix_comm,
          solvers[id]->n,
          solvers[id]->colptr,
          solvers[id]->rows,
          solvers[id]->values,
          solvers[id]->l2g,
          solvers[id]->perm,
          NULL,
          solvers[id]->b,
          solvers[id]->nrhs,
          solvers[id]->pastix_data->iparm,
          solvers[id]->pastix_data->dparm);

  MURGE_MEMALLOC(dest, (strlen(directory)+20), char);
  MURGE_MEMALLOC(src, 20, char);

  if (procnum == 0)
    {
      sprintf(dest,"%s/ordergen", directory);
      RENAME("ordergen", dest);
      sprintf(dest,"%s/symbgen", directory);
      RENAME("symbgen", dest);

      sprintf(dest,"%s/murge.save", directory);
      PASTIX_FOPEN(stream, dest, "w");
      fclose(stream);

    }

  memFree_null(src);
  memFree_null(dest);

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_Load
 *
 * Loads preprocessing result from disk, into *directory*,
 * where it had been saved by <MURGE_Save>.
 *
 * If preprocessing data was already computed or loaded, it will
 * be overwriten.
 *
 * Needs <MURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   directory - Path to the directory where to load the solver
 *               preprocessing data.
 *
 * In Fortran, *STR_LEN* is the length of the string directory.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 *   MURGE_ERR_IO        - If file(s) couldn't be read.
 */
INTS MURGE_Load(INTS id, char* directory)
{
#ifdef DISTRIBUTED
  char * src     = NULL;
  int    procnum;

  CHECK_SOLVER_ID(id);
  procnum = (int)solvers[id]->pastix_data->procnum;

  if (solvers[id]->pastix_data->iparm == NULL)
    {
      errorPrint("You need to call MURGE_SetDefaultOptions before");
      return MURGE_ERR_ORDER;
    }

  solvers[id]->pastix_data->iparm[IPARM_IO_STRATEGY] = API_IO_LOAD;
  solvers[id]->pastix_data->iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
  solvers[id]->pastix_data->iparm[IPARM_END_TASK]    = API_TASK_BLEND;

  pastix_welcome_print(solvers[id]->pastix_data,
                       solvers[id]->colptr,
                       solvers[id]->n);

  if (NULL == solvers[id]->perm)
    MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, INT);
  MURGE_MEMALLOC(src, (strlen(directory)+20), char);
  if (procnum == 0)
    {
      sprintf(src,"%s/ordergen", directory);
      LINK(src, "ordername");
      sprintf(src,"%s/symbgen", directory);
      LINK(src, "symbname");
    }
  MPI_Barrier(solvers[id]->pastix_data->pastix_comm);
  memFree_null(src);

  dpastix(&(solvers[id]->pastix_data),
          solvers[id]->pastix_data->pastix_comm,
          solvers[id]->n,
          solvers[id]->colptr,
          solvers[id]->rows,
          solvers[id]->values,
          solvers[id]->l2g,
          solvers[id]->perm,
          NULL,
          solvers[id]->b,
          solvers[id]->nrhs,
          solvers[id]->pastix_data->iparm,
          solvers[id]->pastix_data->dparm);
  solvers[id]->N = solvers[id]->pastix_data->ordemesh.rangtab[solvers[id]->pastix_data->ordemesh.cblknbr];
  if (procnum == 0)
    {
      UNLINK("ordername");
      UNLINK("symbname");
    }


  MURGE_STATE_TRUE(solvers[id]->state, MURGE_SYMB_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_GRAPH_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/*******************************************************************************
 * Group: Getting solver's distribution
 */

/*
 * Function: MURGE_GetLocalNodeNbr
 *
 * Computes preprocessing step, if not done, and the number of
 * Nodes in the new ditribution of the matrix.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   nodenbr   - *INTS* where to store number of nodes.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *nodenbr* is *NULL* (can occur in C).
*/
INTS MURGE_GetLocalNodeNbr    (INTS id, INTS *nodenbr){
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  *nodenbr = (INTS)pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODENBR_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_GetLocalNodeList
 *
 * Computes the local node list, corresponding to
 * the new distribution, after preprocessing.
 *
 * *nodelist* array has to be allocated before calling
 * this function.
 *
 * As it's result determines the size of *nodelist*
 * array, <MURGE_GetLocalNodeNbr> should be run before it.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   nodelist  - Array where to store the list of local nodes.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - if <MURGE_GetLocalNodeNbr> has not been called
 *                         before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *nodelist* is *NULL* (can occur in C).
 */
INTS MURGE_GetLocalNodeList   (INTS id, INTS *nodelist)
{
#ifdef DISTRIBUTED
  int ret;
  INT nodenbr = 0;
  INT i;
  INT * intern_nodelist;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODENBR_OK)))
    {
      errorPrint("You need to call MURGE_GetLocalNodeNbr before");
      return MURGE_ERR_ORDER;
    }
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODELST_OK);

  if (sizeof(INT) != sizeof(INTS))
    {
      nodenbr = pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));
      MURGE_MEMALLOC(intern_nodelist, nodenbr, INT);
    }
  else
    {
      intern_nodelist = (INT*)nodelist;
    }

  if (EXIT_SUCCESS != ( ret =
                        pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                                               intern_nodelist)))
    return MURGE_ERR_SOLVER;

  if (sizeof(INT) != sizeof(INTS))
    {
      for (i = 0; i < nodenbr; i++)
        {
          nodelist[i] = (INTS) intern_nodelist[i];
        }
      memFree_null(intern_nodelist);
    }
  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
    {
      for (i = 0; i < nodenbr; i++)
        nodelist[i] -= 1;
    }

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_GetLocalUnkownNbr
 *
 * Computes preprocessing step, if not done, and the number of
 * Unkowns in the new ditribution of the matrix.
 *
 * Parameters:
 *   id            - Solver instance identification number.
 *   unkownnbr     - *INTS* where to store number of unkowns.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *unkownnbr* is *NULL* (can occur in C).
 */
INTS MURGE_GetLocalUnknownNbr (INTS id, INTS *unkownnbr)
{
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  *unkownnbr = (INTS)pastix_getLocalUnknownNbr(&(solvers[id]->pastix_data));
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODENBR_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_GetLocalUnkownList
 *
 * Computes the local unkown list, corresponding to
 * the new distribution, after preprocessing.
 *
 * *unkownlist* array has to be allocated before calling
 * this function.
 *
 * As it's result determines the size of *unkownlist*
 * array, <MURGE_GetLocalUnkownNbr> should be run before it.
 *
 * Parameters:
 *   id          - Solver instance identification number.
 *   unkownlist  - Array where to store the list of local unkowns.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - if <MURGE_GetLocalUnkownNbr> has not been called
 *                         before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *unkownlist* is *NULL* (can occur in C).
 */
INTS MURGE_GetLocalUnknownList(INTS id, INTS *unkownlist){
#ifdef DISTRIBUTED
  int ret;
  INT nodenbr = 0;
  INT i;
  INT * intern_nodelist = NULL;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODENBR_OK)))
    {
      errorPrint("You need to call MURGE_GetLocalNodeNbr before");
      return MURGE_ERR_ORDER;
    }

  nodenbr = pastix_getLocalUnknownNbr(&(solvers[id]->pastix_data));
  if (sizeof(INT) != sizeof(INTS))
    {
      MURGE_MEMALLOC(intern_nodelist, nodenbr, INT);
    }
  else
    {
      intern_nodelist = (INT*)unkownlist;
    }

  if (EXIT_SUCCESS != ( ret =
                        pastix_getLocalUnknownLst(&(solvers[id]->pastix_data),
                                                  intern_nodelist)))
    return MURGE_ERR_SOLVER;

  if (sizeof(INT) != sizeof(INTS))
    {
      for (i = 0; i < nodenbr; i++)
        {
          unkownlist[i] = (INTS) intern_nodelist[i];
        }
      memFree_null(intern_nodelist);
    }
  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
    {
      for (i = 0; i < nodenbr; i++)
        unkownlist[i] -= 1;
    }
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}



/*******************************************************************************
 * Group: Graph setup functions
 */

/*
 * Function: MURGE_GraphBegin
 *
 * - Allocate temporary structure which will contain graph entries.
 * - Set number of unkowns in the graph.
 * - Set the number of entries that are expected in this building session.
 * - Reset the number of entries for this build session.
 * - Set all states except MURGE_GRAPH_BUILD to FALSE (graph, values, blend,
 * nodelst, nodenbr, facto)
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   N       - Number of unkowns.
 *   edgenbr - Number of edges in this building session.
 *             If edgenbr is negative, PaStiX will perform dynamic
 *             reallocation of the array, with the first allocation of
 *             size -edgenbr.
 *
 * Returns:
 *   MURGE_ERR_ORDER     - MURGE_GraphBegin has already been called, or if
 *                         *solvers* or *solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - If *id* is not in correct range.
 *   MURGE_SUCCESS       - Otherwise.
 */
INTS MURGE_GraphBegin(INTS id, INTS N, INTL edgenbr)
{
#ifdef DISTRIBUTED
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD))
    {
      errorPrint("MURGE_GraphBegin has been called before");
      return MURGE_ERR_ORDER;
    }

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (edgenbr < 0) {
    edgenbr = -edgenbr;
    solvers[id]->dynamic = API_YES;
  }
  else {
    solvers[id]->dynamic = API_NO;
  }
  MURGE_MEMALLOC(solvers[id]->tmpijv, edgenbr, ijv_t);
  memset(solvers[id]->tmpijv, 0, edgenbr*sizeof(ijv_t));

  solvers[id]->N       = N;
  solvers[id]->edgenbr = edgenbr;
  solvers[id]->cnt     = 0;

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);

  MURGE_STATE_TRUE(solvers[id]->state,  MURGE_GRAPH_BUILD);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_GraphEdge
 *
 * - Check that the number of entries has not been reach for
 * this session.
 * - Increments ROW and COL if baseval is set to 0.
 * - Checks that ROW and COL ranges are corrects.
 * - Adds an entry to the temporary ijv structure.
 *
 * Parameters:
 *   id  - Solver instance identification number.
 *   ROW - Row of the entry.
 *   COL - Column of the entry.
 *
 * Return:
 *   MURGE_ERR_ORDER     - if we are not in a graph building session, or if
 *                         two many edges have been entered, or if
 *                         *solvers* or *solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - *ROW* or *COL* are out of range or if *id* is not
 *                         in correct range.
 *   MURGE_SUCCESS       - Otherwise
 */
INTS MURGE_GraphEdge (INTS id, INTS ROW, INTS COL)
{
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD)))
    {
      errorPrint("Need to call MURGE_GraphBegin first");
      return MURGE_ERR_ORDER;
    }
  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
    {
      COL += 1;
      ROW += 1;
    }

  if (ROW < 1 || COL < 1 || ROW > solvers[id]->N || COL > solvers[id]->N)
    {
      errorPrint("ROW or COL is out of range");
      return MURGE_ERR_PARAMETER;
    }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (solvers[id]->cnt + 1 > solvers[id]->edgenbr)
    {
      if (solvers[id]->dynamic == API_NO) {
        errorPrint("Too many edges in graph description (%ld > %ld)",
                   solvers[id]->cnt + 1, solvers[id]->edgenbr);
#ifdef MURGE_THREADSAFE
        pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif

        return MURGE_ERR_ORDER;
      }
      else {
        INT old_edgnbr = solvers[id]->edgenbr;
        solvers[id]->edgenbr += solvers[id]->edgenbr/2 + 1;
        MURGE_REALLOC(solvers[id]->tmpijv, solvers[id]->edgenbr, ijv_t);
        memset(&(solvers[id]->tmpijv[old_edgnbr]),
               0, (solvers[id]->edgenbr-old_edgnbr)*sizeof(ijv_t));
      }
    }

  solvers[id]->tmpijv[solvers[id]->cnt].i = ROW;
  solvers[id]->tmpijv[solvers[id]->cnt].j = COL;
  solvers[id]->cnt++;
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_GraphEnd
 *
 * - Sort temporary IJV structure with cols as key.
 * - Distribute columns onto processors.
 * (first column on first proc and so on...)
 * - Build a distributed CSC that will be given to PaStiX.
 *
 * TODO:
 * - In the case of a triangular matrix, count each extra-diagonal twice.
 * - Use initial distribution to compute column distribution,
 * in order to reduce communications.
 *
 * Parameters :
 *   id  - Solver instance identification number.
 *
 * Returns:
 *   MURGE_ERR_ORDER     - if we are not in a graph building session, or if
 *                         all edges have not been entered, or if
 *                         *solvers* or *solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - *ROW* or *COL* are out of range or if *id* is not
 *                         in correct range.
 *   MURGE_SUCCESS       - Otherwise
 */
INTS MURGE_GraphEnd  (INTS id)
{
#ifdef DISTRIBUTED
  INT              currentcol;
  INT              ncol = 0;
  INT              iter;
  INT             *sizecols      = NULL;
  INT             *sizecols_recv = NULL;
  INTS            *coldist       = NULL;
  INT              totaledgenbr;
  INT              moyedgenbr;
  INT             *localedges    = NULL;
  INT              procnum;
  INT             *tosend        = NULL;
  INT             *torecv        = NULL;
  INT              newedgenbr;
  MPI_Request     *requests      = NULL;
  MPI_Request     *requests2     = NULL;
  MPI_Status       status;
  INT              index;
  ijv_t           *tmpijv        = NULL;
  INT              baseval;
  INT              iter2;
  pastix_data_t   *pastix_data = solvers[id]->pastix_data;
#ifdef MURGE_TIME
  Clock            clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

#ifdef CENTRALISED
  pastix_data->iparm[IPARM_GRAPHDIST] = API_NO;
#endif

  /*
   * Checking that the function is called at the right time
   */
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD)))
    {
      errorPrint("Need to call MURGE_GraphBegin first");
      return MURGE_ERR_ORDER;
    }

  if (solvers[id]->dynamic == API_NO &&
      solvers[id]->cnt < solvers[id]->edgenbr)
    {
      errorPrint("Missing edges entry, expected %ld, entered %ld",
                 (long)solvers[id]->edgenbr,
                 (long)solvers[id]->cnt);
      return MURGE_ERR_ORDER;
    }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (solvers[id]->cnt > solvers[id]->edgenbr)
    {
      errorPrint("Too many edges entry, expected %ld, entered %ld",
                 (long)solvers[id]->edgenbr,
                 (long)solvers[id]->cnt);
#ifdef MURGE_THREADSAFE
      pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
      return MURGE_ERR_ORDER;
    }
  /*
   * Sort the temporary structure following column numbers
   */
  qsort(solvers[id]->tmpijv,
        solvers[id]->cnt,
        sizeof(ijv_t),
        cmp_ijv);

  /*
   * Remove doubles
   */

  index = 0;
  if (solvers[id]->cnt >0)
    {
      for (iter = 0; iter < solvers[id]->cnt; iter++)
        {
          if (solvers[id]->tmpijv[iter].i != solvers[id]->tmpijv[index].i||
              solvers[id]->tmpijv[iter].j != solvers[id]->tmpijv[index].j)
            {
              index++;
              solvers[id]->tmpijv[index].i = solvers[id]->tmpijv[iter].i;
              solvers[id]->tmpijv[index].j = solvers[id]->tmpijv[iter].j;
            }
        }
      solvers[id]->cnt= index+1;
    }

  /*
   * Count the number of element in each column.
   */

  /* decision de la repartition sur les procs */
  MURGE_MEMALLOC(sizecols, solvers[id]->N, INT);

  for (iter = 0; iter < solvers[id]->N; iter ++)
    sizecols[iter] = 0;

  /* on compte combien on en a sur chaque colonne */
  for (iter = 0; iter < solvers[id]->cnt; iter ++)
    {
      /******************************************************************
       * TODO: Dans le cas ou la matrice donnée ne contient que la partie
       *       triangulaire, il faudra  compter 2 fois les elements
       *       non diagonaux.
       ******************************************************************/
      sizecols[solvers[id]->tmpijv[iter].j - 1]++;
    }


  MURGE_MEMALLOC(sizecols_recv, solvers[id]->N, INT);

  for (iter = 0; iter < solvers[id]->N; iter ++)
    sizecols_recv[iter] = 0;

  MPI_Allreduce (sizecols, sizecols_recv,
                 solvers[id]->N, COMM_INT,
                 MPI_SUM, pastix_data->pastix_comm);


  totaledgenbr = 0;
  for (iter = 0; iter < solvers[id]->N; iter ++)
    totaledgenbr += sizecols_recv[iter];

  moyedgenbr = totaledgenbr/pastix_data->procnbr;

  MURGE_MEMALLOC(coldist, solvers[id]->N, INTS);

  /*
   * Distribute the column.
   */

  /* TODO: prendre en compte la distribution initiale */

  for (iter = 0; iter < solvers[id]->N; iter++)
    coldist[iter] = pastix_data->procnbr - 1;

  procnum    = 0;
  iter       = 0;

  MURGE_MEMALLOC(localedges, pastix_data->procnbr , INT);

  while (iter < solvers[id]->N)
    {
      localedges[procnum] = 0;
      while ((iter < solvers[id]->N) &&
             ((localedges[procnum] < moyedgenbr)||
              (procnum == pastix_data->procnbr -1) ))
        {
          coldist[iter] = procnum;
          if (procnum == pastix_data->procnum)
            {
              ncol++;
            }
          /* columns belonging to procnum */
          localedges[procnum] +=  sizecols_recv[iter];
          iter ++;
        }

      procnum++;

    }

  MURGE_MEMALLOC(solvers[id]->l2g, ncol, INT);
  ncol = 0;
  iter = 0;
  while (iter < solvers[id]->N)
    {
      procnum = coldist[iter];
      if (procnum == pastix_data->procnum)
        {
          solvers[id]->l2g[ncol] = iter + 1;
          ncol++;
        }
      iter ++;
    }

  MURGE_MEMALLOC(tosend, pastix_data->procnbr, INT);
  for (iter = 0; iter < pastix_data->procnbr; iter++)
    tosend[iter] = 0;
  for (iter = 0; iter < solvers[id]->N; iter++)
    {
      tosend[coldist[iter]] += sizecols[iter];
    }
#ifdef CENTRALISED
  /* Pour les appels avec graphe centralisé on envoi tout */
  for (iter = 0; iter < pastix_data->procnbr; iter++)
    {
      tosend[iter] = solvers[id]->cnt;
      localedges[iter] = totaledgenbr;
    }
  ncol = solvers[id]->N;
#endif /* CENTRALISED */
  MURGE_MEMALLOC(requests, pastix_data->procnbr, MPI_Request);

  /* envoi du nombre d'arretes a echanger a chaque proc */
  for (procnum = 0; procnum < pastix_data->procnbr; procnum++)
    {
      if (procnum != pastix_data->procnum)
        {
          MPI_Isend(&tosend[procnum],
                    1, COMM_INT, procnum,
                    TAG_SIZE, pastix_data->pastix_comm,
                    &requests[procnum]);
        }
    }

  MURGE_MEMALLOC(torecv, pastix_data->procnbr, INT);

  /*reception du nombre d'arretes a recevoir de chaque proc */
  newedgenbr = tosend[pastix_data->procnum];
  for (procnum = 0; procnum < pastix_data->procnbr; procnum++)
    {
      if (procnum != pastix_data->procnum)
        {
          MPI_Recv(&torecv[procnum],
                   1, COMM_INT, procnum,
                   TAG_SIZE, pastix_data->pastix_comm,
                   &status);
          newedgenbr += torecv[procnum];
        }
      else
        {
          torecv[procnum] = tosend[procnum];
        }
    }

  for (procnum = 0; procnum < pastix_data->procnbr; procnum++)
    {
      if (procnum != pastix_data->procnum)
        {
          MPI_Wait(&requests[procnum],
                   &status);
        }
    }


  MURGE_MEMALLOC(tmpijv, localedges[pastix_data->procnum], ijv_t);

  MURGE_MEMALLOC(requests2, pastix_data->procnbr, MPI_Request);

  index = 0;
  for (procnum = 0; procnum < pastix_data->procnbr; procnum++)
    {
      if (procnum != pastix_data->procnum)
        {
#ifdef CENTRALISED
          /* On envoi tout */
          index = 0;
#endif
          MPI_Isend(&(solvers[id]->tmpijv[index]),
                    tosend[procnum]*sizeof(ijv_t), MPI_BYTE, procnum,
                    TAG_IJV, pastix_data->pastix_comm,
                    &requests[procnum]);
        }
      index += tosend[procnum];

    }
  index = 0;
  for (procnum = 0; procnum < pastix_data->procnbr; procnum++)
    {
      if (procnum != pastix_data->procnum)
        {
          MPI_Recv(&(tmpijv[index]),
                   torecv[procnum]*sizeof(ijv_t), MPI_BYTE, procnum,
                   TAG_IJV, pastix_data->pastix_comm, &status);
        }
      else
        {
          if (torecv[procnum] != 0)
            {
              iter = 0;
#ifndef CENTRALISED
              while (iter < solvers[id]->edgenbr &&
                     coldist[solvers[id]->tmpijv[iter].j-1] != procnum)
                iter++;
#endif
              if (iter == solvers[id]->edgenbr)
                {
#ifdef MURGE_THREADSAFE
                  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
                  return MURGE_ERR_SOLVER;
                }
              memcpy(&(tmpijv[index]),
                     &(solvers[id]->tmpijv[iter]),
                     torecv[procnum]*sizeof(ijv_t));
            }
        }
      index += torecv[procnum];
    }
  for (procnum = 0; procnum < pastix_data->procnbr; procnum++)
    {
      if (procnum != pastix_data->procnum)
        {
          MPI_Wait(&requests[procnum],
                   &status);
        }
    }
  memFree_null(solvers[id]->tmpijv);
  memFree_null(requests);
  memFree_null(requests2);
  /* on retrie le nouveau tmpijv */
  qsort(tmpijv,
        localedges[pastix_data->procnum],
        sizeof(ijv_t),
        cmp_ijv);
  /*
   * Remove doubles
   */
  index = 0;
  if (localedges[pastix_data->procnum] >0)
    {
      for (iter = 0; iter < localedges[pastix_data->procnum]; iter++)
        {
          if (tmpijv[iter].i != tmpijv[index].i||
              tmpijv[iter].j != tmpijv[index].j)
            {
              index++;
              tmpijv[index].i = tmpijv[iter].i;
              tmpijv[index].j = tmpijv[iter].j;
            }
        }
      localedges[pastix_data->procnum]= index+1;
    }

  solvers[id]->n = ncol;

  MURGE_MEMALLOC(solvers[id]->colptr, ncol+1, INT);

  baseval=1; /* Attention on base a 1 */
  iter2=baseval;
  for (iter=0; iter<(ncol); iter++)
    {
      solvers[id]->colptr[iter] = iter2;
      if (iter2 - baseval < localedges[solvers[id]->pastix_data->procnum])
        {
          currentcol = tmpijv[iter2-baseval].j;
          while (((iter2-baseval) < localedges[pastix_data->procnum])
                 && (tmpijv[iter2-baseval].j  == currentcol))
            {
              iter2++;
            }
        }
    }
  solvers[id]->colptr[ncol] = iter2;

  MURGE_MEMALLOC(solvers[id]->rows,
                 localedges[pastix_data->procnum],
                 INT);

  for (iter=0; iter<localedges[pastix_data->procnum]; iter++)
    {
      solvers[id]->rows[iter] = tmpijv[iter].i;
    }

  memFree_null(tmpijv);
  memFree_null(tosend);
  memFree_null(torecv);
  memFree_null(coldist);
  memFree_null(localedges);
  memFree_null(sizecols);
  memFree_null(sizecols_recv);
  CLOCK_STOP;
#ifdef MURGE_TIME
  fprintf(stdout, " > MURGE_GraphEnd computed in %.3lg seconds (id=%d)\n", (double)CLOCK_GET, id);
#endif
  MURGE_DUMP_GRAPH;

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_GRAPH_BUILD);

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 * Function: MURGE_GraphGlobalCSR
 *
 * Enter the adjency graph in a Column Sparse Row form.
 *
 *
 * If the matrix is symmetric, calls <MURGE_GraphGlobalCSC>
 * else uses <MURGE_GraphBegin>, <MURGE_GraphEdge>,
 * <MURGE_GraphEnd> sequence.
 *
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of rows in the CSR.
 *   rowptr - Indexes of each row in COLS array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS MURGE_GraphGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS, INTS root)
{
#ifdef DISTRIBUTED
  INT  iter;
  INT  iter2;
  INTS ret;
  int  baseval;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) )
    {
      errorPrint("Do not call MURGE_GraphBegin before");
      return MURGE_ERR_ORDER;
    }

  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  /* Si on a un graph symetrique autant faire un MURGE_GraphGlobalCSC */
  if (solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES)
    return MURGE_GraphGlobalCSC(id, N, rowptr, COLS, root);


  if (solvers[id]->pastix_data->procnum == root || root == -1)
    {
      if (MURGE_SUCCESS != (ret = MURGE_GraphBegin(id, N, rowptr[N]- baseval)))
        return ret;

      for (iter = 0; iter < N; iter++)
        {
          for (iter2 = rowptr[iter]; iter2 < rowptr[iter+1]; iter2++)
            {
              ret =
                MURGE_GraphEdge(id,
                                iter + baseval,
                                COLS[iter2 - baseval]);
              if (MURGE_SUCCESS != ret)
                return ret;
            }
        }
      if (MURGE_SUCCESS != (ret = MURGE_GraphEnd(id)))
        return ret;
    }
  else
    {
      if (MURGE_SUCCESS != (ret = MURGE_GraphBegin(id, N, 0)))
        return ret;
      if (MURGE_SUCCESS != (ret = MURGE_GraphEnd(id)))
        return ret;
    }

  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 * Function: MURGE_GraphGlobalCSC
 *
 * Distribute the CSC on the processors and use it for PaStiX calls.
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of columns in the CSR.
 *   colptr - Indexes of each columns in ROWS array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS MURGE_GraphGlobalCSC(INTS id, INTS N, INTL *colptr, INTS *ROWS, INTS root)
{
#ifdef DISTRIBUTED
  INT          globaledgenbr;
  INT          averageedgenbr;
  INT          localedgenbr;
  INT          firstcol;
  INT          lastcol         = 0;
  INT          iter;
  INT          procnum;
  INT          firstlast[2];
  INTS        *tmpj            = NULL;
  MPI_Request *requests_fl     = NULL;
  MPI_Request *requests_colptr = NULL;
  MPI_Request *requests_rows   = NULL;
  MPI_Status   status;
  int          baseval;             /* User baseval               */

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) )
    {
      errorPrint("Do not call MURGE_GraphBegin before");
      return MURGE_ERR_ORDER;
    }

  solvers[id]->N = N;
  baseval        = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];

  if (root == solvers[id]->pastix_data->procnum || root == -1)
    {
      globaledgenbr   = colptr[N] - baseval;
      averageedgenbr  = globaledgenbr / solvers[id]->pastix_data->procnbr;
      firstcol        = 0;
      if (root != -1)
        {
          MURGE_MEMALLOC(requests_fl,
                         solvers[id]->pastix_data->procnbr,
                         MPI_Request);
          MURGE_MEMALLOC(requests_colptr,
                         solvers[id]->pastix_data->procnbr,
                         MPI_Request);
          MURGE_MEMALLOC(requests_rows,
                         solvers[id]->pastix_data->procnbr,
                         MPI_Request);
        }
      /* Pour chaque prrocesseur,
       on attribue au processeur un certain nombre de colonnes et donc
       d'arrêtes

       Si le processeur est local on construit le loc2glob
       On construit le colptr
       On copie rows

       Sinon, on envoi les numéros de première et dernière colonnes et
       les morceaux de colptr et rows correspondants.

       */
      for (procnum = 0; procnum <  solvers[id]->pastix_data->procnbr; procnum++)
        {
          while (lastcol < N - 1&&
                 colptr[lastcol+1] - colptr[firstcol]< averageedgenbr)
            lastcol++;
          localedgenbr =  colptr[lastcol+1] - colptr[firstcol];

          if (procnum == solvers[id]->pastix_data->procnum)
            {
              solvers[id]->n = lastcol-firstcol+1;

              MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, INT);

              for (iter = 0; iter < solvers[id]->n; iter++)
                {
                  solvers[id]->l2g[iter] = firstcol+iter+1;
                }

              MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, INT);

              for (iter = 0; iter < solvers[id]->n+1; iter++)
                {
                  solvers[id]->colptr[iter] =
                    colptr[firstcol+iter] - colptr[firstcol]+1;
                }

              MURGE_MEMALLOC(solvers[id]->rows, localedgenbr, INT);

              for (iter = 0; iter < localedgenbr; iter++)
                solvers[id]->rows[iter] = ROWS[colptr[firstcol]+iter-1];

            }
          else
            {
              if (root != -1)
                {
                  firstlast[0] = firstcol;
                  firstlast[1] = lastcol;


                  MPI_Isend(firstlast,
                            2, COMM_INT, procnum,
                            TAG_FL, solvers[id]->pastix_data->pastix_comm,
                            &requests_fl[procnum]);

                  MPI_Isend(&colptr[firstcol],
                            lastcol-firstcol+2,
                            COMM_INT, procnum,
                            TAG_COL, solvers[id]->pastix_data->pastix_comm,
                            &requests_colptr[procnum]);

                  MPI_Isend(&ROWS[colptr[firstcol]],
                            localedgenbr*sizeof(INTS),
                            MPI_BYTE, procnum,
                            TAG_ROW, solvers[id]->pastix_data->pastix_comm,
                            &requests_rows[procnum]);

                }
            }
          firstcol = lastcol + 1;
        }
      if (root != -1)
        {
          for (procnum = 0;
               procnum < solvers[id]->pastix_data->procnbr;
               procnum++)
            {
              if (procnum != solvers[id]->pastix_data->procnum)
                {
                  MPI_Wait(&requests_fl[procnum], &status);
                  MPI_Wait(&requests_colptr[procnum], &status);
                  MPI_Wait(&requests_rows[procnum], &status);
                }
            }
        }
      memFree_null(requests_rows);
      memFree_null(requests_colptr);
      memFree_null(requests_fl);
    }
  else
    {
      /* Si on est pas le processeur racine

       On recoit les numeros de première et dernière colonnes
       On en déduit le loca2glob
       On recoit les parties locales de colptr et rows
       On construit le colptr local et rows local.
       */
      MPI_Recv(firstlast,
               2, COMM_INT, root,
               TAG_FL, solvers[id]->pastix_data->pastix_comm,
               &status);
      firstcol = firstlast[0];
      lastcol  = firstlast[1];

      solvers[id]->n = lastcol-firstcol+1;

      MURGE_MEMALLOC(solvers[id]->l2g, lastcol-firstcol+1, INT);

      for (iter = 0; iter < lastcol-firstcol+1; iter++)
        {
          solvers[id]->l2g[iter] = firstcol+iter;
        }

      MURGE_MEMALLOC(solvers[id]->colptr, lastcol-firstcol+2, INT);

      MPI_Recv(solvers[id]->colptr,
               lastcol-firstcol+2, COMM_INT, root,
               TAG_COL, solvers[id]->pastix_data->pastix_comm,
               &status);


      for (iter = 0; iter < lastcol-firstcol+2; iter++)
        {
          solvers[id]->colptr[lastcol-firstcol+1 - iter] -=
            solvers[id]->colptr[0];
        }

      localedgenbr = solvers[id]->colptr[lastcol-firstcol+1]-1;

      MURGE_MEMALLOC(tmpj, localedgenbr, INTS);

      MPI_Recv(tmpj,
               localedgenbr*sizeof(INTS), MPI_BYTE, root,
               TAG_ROW, solvers[id]->pastix_data->pastix_comm,
               &status);

      if (sizeof(INTS) == sizeof(INT))
        {
          solvers[id]->rows = (INT*)tmpj;
        }
      else
        {
          MURGE_MEMALLOC(solvers[id]->rows, localedgenbr, INT);

          for (iter = 0; iter < localedgenbr; iter++)
            solvers[id]->rows[iter] = (INT)tmpj[iter];
          memFree_null(tmpj);
        }
    }

  MURGE_DUMP_GRAPH;

  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
/*
 * Function: MURGE_GraphGlobalIJV
 *
 * Distribute the graph on the processors, compress the columns
 * array and use the built CSCd to call PaStiX.
 *
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of columns in the CSR.
 *   NNZ    - Number of non-zeros in the matrix.
 *   ROWS   - Rows array.
 *   COLS   - Columns array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS MURGE_GraphGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROWS,
                          INTS *COLS, INTS root)
{
#ifdef DISTRIBUTED
  INT       *localn     = NULL;   /* Number of local column on each proc  */
  INT       *localedges = NULL;   /* Number of local edges on each proc   */
  INT        lnnz;                /* Local number of edges                */
  INT        sizes[2];            /* Array to send n and nnz to each proc */
  INTS      *tmprows  = NULL;     /* Temporary local rows tabular         */
  INTS      *tmpcols  = NULL;     /* Temporary local columns tabular      */
  INT       *sizecols = NULL;     /* Number of rows in each column        */
  INT        totaledgenbr;        /* Total number of edges                */
  INT        avredgenbr;          /* Average number of edges              */
  int        baseval_int     = 1; /* Internal baseval, always 1           */
  int        baseval;             /* User baseval                         */
  INT        iter, iter2, index;  /* Iterators                            */
  INT       *coldist  = NULL;     /* Owner of each column                 */
  int        procnum;             /* Processor number iterator            */
  MPI_Status status;              /* MPI status */
  INTS       currentcol;
  void      *sortptr[2];

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) )
    {
      errorPrint("Do not call MURGE_GraphBegin before");
      return MURGE_ERR_ORDER;
    }


  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];


  if ((solvers[id]->pastix_data->procnum == root) || (root == -1))
    {
      solvers[id]->N = N;
    }
  if (root != -1)
    {
      MPI_Bcast(&solvers[id]->N, 1, COMM_INT, root,
                solvers[id]->pastix_data->pastix_comm);
      MPI_Bcast(&NNZ, sizeof(INTL), MPI_BYTE, root,
                solvers[id]->pastix_data->pastix_comm);
    }
  if ((solvers[id]->pastix_data->procnum == root) || (root == -1))
    {
      /* Sort Col and Rows with first key, col */
      sortptr[0] = COLS;
      sortptr[1] = ROWS;
#ifdef INTSSIZE64
      qsort2IntAsc(sortptr, NNZ);
#else
      qsort2SmallIntAsc(sortptr, NNZ);
#endif
      /* decide how to distribute the graph */
      MURGE_MEMALLOC(sizecols, solvers[id]->N, INT);
      memset(sizecols, 0, solvers[id]->N*sizeof(INT));
      totaledgenbr = 0;
      /* Count how long is each column */
      for (iter = 0; iter < NNZ; iter ++)
        {
          /* TODO: Dans le cas ou la matrice donnée ne contient que la partie
           triangulaire, il faudra  compter 2 fois les elements
           non diagonaux.
           */
          sizecols[COLS[iter] - 1]++;
          totaledgenbr++;
        }

      avredgenbr = totaledgenbr/solvers[id]->pastix_data->procnbr;

      MURGE_MEMALLOC(coldist, solvers[id]->N, INT);

      for (iter = 0; iter < solvers[id]->N; iter++)
        coldist[iter] = solvers[id]->pastix_data->procnbr - 1;

      procnum    = 0;
      iter       = 0;

      MURGE_MEMALLOC(localedges, solvers[id]->pastix_data->procnbr, INT);
      MURGE_MEMALLOC(localn,     solvers[id]->pastix_data->procnbr, INT);
      memset(localedges, 0, solvers[id]->pastix_data->procnbr*sizeof(INT));
      memset(localn, 0, solvers[id]->pastix_data->procnbr*sizeof(INT));

      while (iter < solvers[id]->N)
        {
          localedges[procnum] = 0;
          while (iter < solvers[id]->N &&
                 (localedges[procnum] < avredgenbr||
                  (procnum == solvers[id]->pastix_data->procnbr -1)))
            {
              coldist[iter] = procnum;
              localn[procnum]++;
              localedges[procnum] +=  sizecols[iter];
              iter ++;
            }
          procnum++;
        }

      memFree_null(coldist);

      /* Send data to each processor */

      for (index = 0, procnum = 0;
           procnum < solvers[id]->pastix_data->procnbr;
           procnum++)
        {
          if (procnum != solvers[id]->pastix_data->procnum)
            {
              if (root != -1)
                {
                  sizes[0] = localn[procnum];
                  sizes[1] = localedges[procnum];

                  /* envoi du nombre de non zeros */
                  MPI_Send(sizes,
                           2, COMM_INT, procnum,
                           TAG_SIZE, solvers[id]->pastix_data->pastix_comm);
                  /* envoi des lignes */
                  MPI_Send(&(ROWS[index]),
                           localedges[procnum]*sizeof(INTS), MPI_BYTE, procnum,
                           TAG_ROW, solvers[id]->pastix_data->pastix_comm);
                  /* envoi des colonnes */
                  MPI_Send(&(COLS[index]),
                           localedges[procnum]*sizeof(INTS), MPI_BYTE, procnum,
                           TAG_COL, solvers[id]->pastix_data->pastix_comm);
                }
            }
          else
            {
              tmprows = &(ROWS[index]);
              tmpcols = &(COLS[index]);
            }
          index += localedges[procnum];
        }

      solvers[id]->n = localn[solvers[id]->pastix_data->procnum];
      lnnz = localedges[solvers[id]->pastix_data->procnum];
      memFree_null(localn);
      memFree_null(localedges);
    }
  else
    {
      MPI_Recv(sizes,
               2, COMM_INT, root,
               TAG_SIZE, solvers[id]->pastix_data->pastix_comm,
               &status);
      solvers[id]->n = sizes[0];
      lnnz           = sizes[1];

      MURGE_MEMALLOC(tmprows, lnnz, INTS);
      MPI_Recv(tmprows,
               lnnz*sizeof(INTS), MPI_BYTE, root,
               TAG_ROW, solvers[id]->pastix_data->pastix_comm,
               &status);
      MURGE_MEMALLOC(tmpcols, lnnz, INTS);
      MPI_Recv(tmpcols,
               lnnz*sizeof(INTS), MPI_BYTE, root,
               TAG_COL, solvers[id]->pastix_data->pastix_comm,
               &status);
    }

  MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, INT);
  MURGE_MEMALLOC(solvers[id]->l2g,    solvers[id]->n  , INT);
  /* convert tmpcols/tmprows to CSCd */
  iter2=baseval_int;
  for (iter=0; iter<(solvers[id]->n); iter++)
    {
      solvers[id]->colptr[iter] = iter2;
      solvers[id]->l2g[iter]    = tmpcols[iter2-baseval_int]
        - baseval + baseval_int;
      currentcol = tmpcols[iter2-baseval_int];
      while (((iter2-baseval) < lnnz) &&
             ((tmpcols[iter2-baseval_int]) == (currentcol)))
        {
          iter2++;
        }
    }

  if ((solvers[id]->pastix_data->procnum != root) && (root != -1))
    memFree_null(tmpcols);

  solvers[id]->colptr[solvers[id]->n] = iter2;

  if (iter2 != lnnz+baseval)
    {
      errorPrint("Mauvais nombre d'arrête");
      return MURGE_ERR_PARAMETER;
    }


  MURGE_MEMALLOC(solvers[id]->rows, lnnz, INT);

  for (iter=0; iter<lnnz; iter++)
    {
      solvers[id]->rows[iter] = tmprows[iter];
    }
  if (!((solvers[id]->pastix_data->procnum == root) || (root == -1)))
    memFree_null(tmprows);

  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

INTS MURGE_SetOrdering(INTS id, INTS * permutation)
{
#ifdef DISTRIBUTED
  INTS i;
  print_debug(DBG_MURGE, ">> MURGE_SetOrdering\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("Build graph before calling MURGE_SetOrdering");
      return MURGE_ERR_ORDER;
    }
  if (solvers[id]->l2g == NULL)
    {
      errorPrint("Local to global array is not set");
      return MURGE_ERR_PARAMETER;
    }
  if (permutation == NULL)
    {
      errorPrint("NULL parameter");
      return MURGE_ERR_PARAMETER;
    }
  MURGE_MEMALLOC(solvers[id]->perm, solvers[id]->n, INT);

  for (i = 0; i < solvers[id]->n; i++)
    solvers[id]->perm[i] = permutation[solvers[id]->l2g[i]-1];

  solvers[id]->pastix_data->iparm[IPARM_ORDERING] = API_ORDER_PERSONAL;
  solvers[id]->pastix_data->iparm[IPARM_LEVEL_OF_FILL] = -1;

  print_debug(DBG_MURGE, "<< MURGE_SetOrdering\n");
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*******************************************************************************
 * Group: Matrix assembly functions
 */

/* Function: MURGE_AssemblySetSequence
 *
 * Create a sequence of entries to build a matrix and store it for being reused.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   coefnbr - Number of entries.
 *   ROWs    - List of rows in the sequence.
 *   COLs    - List of columns in the sequence.
 *   op      - Operation to perform for coefficient which appear
 *             several tim (see <MURGE_ASSEMBLY_OP>).
 *   op2     - Operation to perform when a coefficient is set by
 *             two different processors (see <MURGE_ASSEMBLY_OP>).
 *   mode    - Indicates if user ensure he will respect solvers distribution
 *             (see <MURGE_ASSEMBLY_MODE>).
 *   nodes   - 0 entries are entered value by value,
 *             1 entries are entries node by node.
 *   id_seq  - Sequence ID.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.
 */
int MURGE_AssemblySetSequence (INTS id, INTL coefnbr, INTS * ROWs, INTS * COLs,
                               INTS op, INTS op2, INTS mode, INTS nodes,
                               INTS * id_seq)
{
#ifdef DISTRIBUTED
  murge_seq_t * sequence     = NULL;
  INTL         iter;
  ijv_t      **send_ijv      = NULL;
  INTL        *send_ijv_size = NULL;
  INTL        *send_nbr      = NULL;
  INTS         dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("Graph has to be built before");
      return MURGE_ERR_ORDER;
    }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);


  MURGE_MEMALLOC(sequence, 1, murge_seq_t);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (solvers[id]->sequences == NULL)
    {
      solvers[id]->sequences = sequence;
    }
  else
    {
      murge_seq_t *last_sequence = solvers[id]->sequences;
      while(sequence->next != NULL)
        {
          last_sequence = sequence->next;
        }
      last_sequence->next = sequence;
    }

  sequence->next         = NULL;
  sequence->indexes      = NULL;
  sequence->recv_indexes = NULL;
  sequence->recv_nbr     = NULL;
  CHOOSE_FUNC(sequence->fusion_local_entries, op);
  CHOOSE_FUNC(sequence->fusion_dist_entries, op);
  sequence->mode  = mode;
  sequence->nodes = (nodes==MURGE_BOOLEAN_FALSE)?API_NO:API_YES;
  sequence->ID = solvers[id]->seq_ID++;
  *id_seq = sequence->ID;

  sequence->coefnbr = coefnbr;
  MURGE_MEMALLOC(sequence->indexes, coefnbr, INTL);
  if (sequence->mode == MURGE_ASSEMBLY_FOOL)
    {
      MURGE_MEMALLOC(send_nbr, solvers[id]->pastix_data->procnbr, INTL);
      MURGE_MEMALLOC(sequence->recv_nbr,
                     solvers[id]->pastix_data->procnbr,
                     INTL);
      MURGE_MEMALLOC(sequence->recv_indexes,
                     solvers[id]->pastix_data->procnbr,
                     INTL*);
      MURGE_MEMALLOC(send_ijv, solvers[id]->pastix_data->procnbr, ijv_t*);
      MURGE_MEMALLOC(send_ijv_size, solvers[id]->pastix_data->procnbr, INTL);

      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        {
          send_nbr[iter] = 0;
          send_ijv_size[iter] = 1+coefnbr/10;
          MURGE_MEMALLOC(send_ijv[iter], send_ijv_size[iter], ijv_t);
          sequence->recv_indexes[iter] = NULL;
        }
    }

  if (solvers[id]->colptr == NULL)
    {
      /* Need to build a CSC */
      ijv_t * mat_ijv = NULL;
      INTS inc = (solvers[id]->pastix_data->iparm[IPARM_BASEVAL]==0)?1:0;

      MURGE_MEMALLOC(mat_ijv, coefnbr, ijv_t);
      /* Copy ROWs and COLs into IJV structso that we can sort it */
      for (iter = 0; iter < coefnbr; iter++)
        {
          INTS node_col;
          INTS node_row;
          INTS node_col_loc;
          if ((COLs[iter]+inc < 1) || (ROWs[iter]+inc < 1) ||
              (ROWs[iter]+inc > ((solvers[id]->N))) ||
              (COLs[iter]+inc > ((solvers[id]->N))))
            {
              errorPrint("COLs[%ld] (%ld) or ROWs[%ld] (%ld) is out of range [1-%ld]",
                         (long)iter, (long)COLs[iter]+inc, (long)iter, (long)ROWs[iter]+inc, (long)(solvers[id]->N));
              return MURGE_ERR_PARAMETER;
            }

          /* If sequence entries are given coefficient by coefficient we have to
           compute the coordinate of the corresponding node.
           Else we already have nodes coordinate.
           */
          if (dof > 1 && sequence->nodes == API_NO) {
            node_col = (COLs[iter]+inc-1 - (COLs[iter]+inc-1)%dof)/dof + 1;
            node_row = (ROWs[iter]+inc-1 - (ROWs[iter]+inc-1)%dof)/dof + 1;
          }
          else {
            node_row = ROWs[iter]+inc;
            node_col = COLs[iter]+inc;
          }
          /* if (node_row == 36084 && node_col == 34866) */
          /*   fprintf(stdout, "ITER %d\n", iter); */

          /* Local column number of the node */
          node_col_loc = solvers[id]->g2l[node_col-1];
          /* if (node_row == 36084 && node_col == 34866) */
          /*   fprintf(stdout, "node_col_loc %d\n", node_col_loc); */
          mat_ijv[iter].i = node_row;
          mat_ijv[iter].j = node_col;

          /* Compute the owner of the entry */
          if ( node_col_loc > 0 )
            {
              mat_ijv[iter].owner =  solvers[id]->pastix_data->procnum;
            }
          else
            {
              if (sequence->mode == MURGE_ASSEMBLY_RESPECT)
                {
                  errorPrint("COLs[%ld] (%ld) is not local",
                             (long)iter, (long)COLs[iter]);
                  return MURGE_ERR_PARAMETER;
                }
              mat_ijv[iter].owner = -node_col_loc;
            }
        }

      if (sequence->mode != MURGE_ASSEMBLY_RESPECT)
        {
          INT     index;
          int   * coefnbr_recv;
          int   * coefnbr_send;
          int   * send_count;
          int   * recv_count;
          int   * send_disp;
          int   * recv_disp;
          ijv_t * mat_ijv_recv;
          INT     total_coefnbr_recv;

          /* Sort using owner, column then row */
          qsort(mat_ijv,
                coefnbr,
                sizeof(ijv_t),
                cmp_ijv_own);

          /*
           * Remove doubles
           */
          index = 0;
          for (iter = 0; iter < coefnbr; iter++)
            {
              if (mat_ijv[iter].i != mat_ijv[index].i||
                  mat_ijv[iter].j != mat_ijv[index].j)
                {
                  index++;
                  mat_ijv[index].i     = mat_ijv[iter].i;
                  mat_ijv[index].j     = mat_ijv[iter].j;
                  mat_ijv[index].owner = mat_ijv[iter].owner;
                }
            }
          coefnbr= index+1;


          /* count the number of entries to send to each processor */
          MURGE_MEMALLOC(coefnbr_send, solvers[id]->pastix_data->procnbr, int);
          MURGE_MEMALLOC(coefnbr_recv, solvers[id]->pastix_data->procnbr, int);
          memset(coefnbr_send, 0, solvers[id]->pastix_data->procnbr*sizeof(int));
          memset(coefnbr_recv, 0, solvers[id]->pastix_data->procnbr*sizeof(int));

          for (iter = 0; iter < coefnbr; iter++)
            {
              coefnbr_send[mat_ijv[iter].owner]++;
            }

          MPI_Alltoall(coefnbr_send, 1, MPI_INT,
                       coefnbr_recv, 1, MPI_INT,
                       solvers[id]->pastix_data->pastix_comm);

          /* Exchange the entries */
          total_coefnbr_recv = 0;
          for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
            total_coefnbr_recv += coefnbr_recv[iter];

          MURGE_MEMALLOC(send_disp, solvers[id]->pastix_data->procnbr, int);
          MURGE_MEMALLOC(send_count, solvers[id]->pastix_data->procnbr, int);
          send_disp[0] = 0;

          for (iter = 1; iter < solvers[id]->pastix_data->procnbr; iter++)
            {
              send_disp[iter] = send_disp[iter-1] +
                coefnbr_send[iter-1]*sizeof(ijv_t);
            }

          for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
            {
              send_count[iter] = coefnbr_send[iter]*sizeof(ijv_t);
            }

          MURGE_MEMALLOC(recv_disp, solvers[id]->pastix_data->procnbr, int);
          MURGE_MEMALLOC(recv_count, solvers[id]->pastix_data->procnbr, int);
          recv_disp[0] = 0;

          for (iter = 1; iter < solvers[id]->pastix_data->procnbr; iter++)
            {
              recv_disp[iter] = recv_disp[iter-1] +
                coefnbr_recv[iter-1]*sizeof(ijv_t);
            }

          for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
            {
              recv_count[iter] = coefnbr_recv[iter]*sizeof(ijv_t);
            }

          MURGE_MEMALLOC(mat_ijv_recv, total_coefnbr_recv, ijv_t);

          MPI_Alltoallv(mat_ijv, send_count, send_disp, MPI_BYTE,
                        mat_ijv_recv, recv_count, recv_disp, MPI_BYTE,
                        solvers[id]->pastix_data->pastix_comm);

          memFree_null(mat_ijv);
          mat_ijv = mat_ijv_recv;
          coefnbr = total_coefnbr_recv;

          memFree_null(coefnbr_send);
          memFree_null(coefnbr_recv);
          memFree_null(send_count);
          memFree_null(recv_count);
          memFree_null(send_disp);
          memFree_null(recv_disp);
        }

      {
        /* build the CSC */
        INTL index;
        INTL iter2;
        INTS baseval;
        INTS currentcol;

        /* Sort elements on columns then rows */
        qsort(mat_ijv,
              coefnbr,
              sizeof(ijv_t),
              cmp_ijv);

        /*
         * Remove doubles
         */
        index = 0;
        for (iter = 0; iter < coefnbr; iter++)
          {
            if (mat_ijv[iter].i != mat_ijv[index].i||
                mat_ijv[iter].j != mat_ijv[index].j)
              {
                index++;
                mat_ijv[index].i = mat_ijv[iter].i;
                mat_ijv[index].j = mat_ijv[iter].j;
              }
          }
        coefnbr= index+1;

        MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, INT);

        baseval=1; /* Attention on base a 1 */
        iter2=baseval;
        for (iter=0; iter<(solvers[id]->n); iter++)
          {
            solvers[id]->colptr[iter] = iter2;
            if (iter2 - baseval < coefnbr)
              {
                currentcol = solvers[id]->l2g[iter];
                while (((iter2-baseval) < coefnbr)
                       && (mat_ijv[iter2-baseval].j  == currentcol))
                  {
                    iter2++;
                  }
              }
            else
              {
                /* Should not go there */
              }
          }
        solvers[id]->colptr[solvers[id]->n] = iter2;

        MURGE_MEMALLOC(solvers[id]->rows,
                       coefnbr,
                       INT);

        for (iter=0; iter < coefnbr; iter++)
          solvers[id]->rows[iter] = mat_ijv[iter].i;

        memFree_null(mat_ijv);
      }
    }

  coefnbr = sequence->coefnbr;

  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD))
    {
      /* TODO: CHECK That all entries exists in CSC and if not insert it*/
    }

  for (iter = 0; iter < coefnbr; iter++)
    {
      INTL iter2;
      INTS inc = (solvers[id]->pastix_data->iparm[IPARM_BASEVAL]==0)?1:0;
      INTS col = COLs[iter]+inc; /* 1 based */
      INTS row = ROWs[iter]+inc; /* 1 based */
      INTS node_col; /* 0 based */
      INTS node_row; /* 0 based */
      INTS in_node_col; /* 0 based */
      INTS in_node_row; /* 0 based */
      INTS node_col_loc; /* 1 based */

      if (dof > 1 && sequence->nodes == API_NO)
        {
          node_col     = (col-1 - (col-1)%dof)/dof;
          in_node_col  = (col-1)%dof;
          node_row     = (row-1 - (row-1)%dof)/dof;
          in_node_row  = (row-1)%dof;
        }
      else
        {
          node_col     = col-1;
          in_node_col  = 0;
          node_row     = row-1;
          in_node_row  = 0;
        }

      node_col_loc = solvers[id]->g2l[node_col];
      if ( node_col_loc > 0 )
        {
          node_col_loc--;
          /* Entry is local */
          for (iter2 = solvers[id]->colptr[node_col_loc]-1;
               iter2 < solvers[id]->colptr[node_col_loc+1]-1;
               iter2++)
            {
              if (solvers[id]->rows[iter2] == row)
                break;
            }
          if (solvers[id]->colptr[node_col_loc+1]-1 == iter2)
            {
              /* Entry not found in CSC */
              errorPrint("ROW (%ld:%ld) not found in COL (%d:%d) %d",
                         (long)row, node_row, (long)col, node_col,node_col_loc);
              return MURGE_ERR_PARAMETER;
            }
          else
            {
              if (iter2*dof*dof + in_node_col*dof+in_node_row >
                  dof*dof*(solvers[id]->colptr[solvers[id]->n]-1))
                {
                  return MURGE_ERR_PARAMETER;
                }

              sequence->indexes[iter] = iter2*dof*dof + in_node_col*dof+in_node_row;
            }
        }
      else
        {
          /* Entry not local */
          if (sequence->mode == MURGE_ASSEMBLY_RESPECT)
            {
              errorPrint("COL (%d) is not local (row %d, owner %d)",
                         (long)col+1, (long)row+1, -node_col_loc);
              return MURGE_ERR_PARAMETER;
            }
          else
            {
              int owner = -node_col_loc;
              sequence->indexes[iter] = -(owner+1);

              /* => send buffer */
              if (send_nbr[owner] == send_ijv_size[owner])
                {
                  send_ijv_size[owner] = 1.5*send_ijv_size[owner] + 1;
                  MURGE_REALLOC(send_ijv[owner], send_ijv_size[owner], ijv_t);
                }
              send_ijv[owner][send_nbr[owner]].i = row;
              send_ijv[owner][send_nbr[owner]].j = col;
              send_nbr[owner]++;
            }
        }
    }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT)
    {
      MPI_Request * requests;
      ijv_t       * recv_ijv;
      int           size;
      int           lastsize;
      INT           iter_coef;

      MPI_Alltoall(send_nbr,           1, MPI_INTL,
                   sequence->recv_nbr, 1, MPI_INTL,
                   solvers[id]->pastix_data->pastix_comm);

      MURGE_MEMALLOC(requests, solvers[id]->pastix_data->procnbr, MPI_Request);
      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        {
          if (send_nbr[iter] > 0)
            MPI_Isend(send_ijv[iter], send_nbr[iter]*sizeof(ijv_t), MPI_BYTE,
                      iter, TAG_IJV, solvers[id]->pastix_data->pastix_comm,
                      &(requests[iter]));
        }

      lastsize = sequence->recv_nbr[0];
      MURGE_MEMALLOC(recv_ijv, lastsize, ijv_t);
      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        {
          MPI_Status status;
          size = sequence->recv_nbr[iter];
          MURGE_MEMALLOC(sequence->recv_indexes[iter], size, INTL);
          if (lastsize < size)
            {
              MURGE_REALLOC(recv_ijv, size, ijv_t);
              lastsize = size;
            }
          if (size > 0)
            MPI_Recv(recv_ijv, size*sizeof(ijv_t), MPI_BYTE,
                     iter, TAG_IJV, solvers[id]->pastix_data->pastix_comm,
                     &status);

          for (iter_coef = 0; iter_coef < size; iter_coef++)
            {
              INTL iter2;
              INTS col = recv_ijv[iter_coef].j;
              INTS row = recv_ijv[iter_coef].i;
              INTS node_col;
              INTS node_row;
              INTS in_node_col;
              INTS in_node_row;
              INTS node_col_loc;


              if (dof > 1 && sequence->nodes == API_NO)
                {
                  node_row = (row-1-(row-1)%dof)/dof;
                  node_col = (col-1-(col-1)%dof)/dof;
                  in_node_row = (row-1)%dof;
                  in_node_col = (col-1)%col;
                }
              else
                {
                  node_row = row-1;
                  node_col = col-1;
                  in_node_row = 0;
                  in_node_col = 0;
                }
              node_col_loc = solvers[id]->g2l[node_col];

              if ( node_col_loc > 0 )
                {
                  /* Entry is local */
                  for (iter2 = solvers[id]->colptr[node_col_loc-1]-1;
                       iter2 < solvers[id]->colptr[node_col_loc]-1;
                       iter2++)
                    {
                      if (solvers[id]->rows[iter2] == row)
                        break;
                    }
                  if (solvers[id]->colptr[node_col_loc]-1 == iter2)
                    {
                      /* Entry not found in CSC */
                      errorPrint("ROW (%ld) not found in COL (%d) %d",
                                 (long)row, (long)col, node_col_loc);

                      return MURGE_ERR_PARAMETER;
                    }
                  else
                    {
                      sequence->recv_indexes[iter][iter_coef] =
                        iter2*dof*dof + in_node_col*dof+in_node_row;
                    }
                }
              else
                {
                  /* Entry not local */
                  errorPrint("%s:%d COL (%d) is not local (row %d, owner %d)",
                             __FILE__, __LINE__, (long)col, row, -node_col_loc);
                  return MURGE_ERR_PARAMETER;
                }
            }
        }

      memFree_null(recv_ijv);

      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        {
          MPI_Status status;

          if (send_nbr[iter] > 0)
            MPI_Wait(&(requests[iter]), &status);
          memFree_null(send_ijv[iter]);
        }
      memFree_null(send_ijv);
      memFree_null(send_ijv_size);
      memFree_null(requests);
      memFree_null(send_nbr);
    }

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

#ifdef DISTRIBUTED
static inline
int try_complete_reception(INTS    id,
                           murge_seq_t * sequence,
                           int     src,
                           int     send_nbr,
                           int     blocksize,
                           FLOAT * recv_data,
                           INT   * recv_cnt)
{
  MPI_Status TCR_status;
  int        TCR_flag = 1;

  while (TCR_flag)
    {
      /* Do we have something to receive ? */
      INT        TCR_iter;
      MPI_Iprobe(src, TAG_VAL,
                 solvers[id]->pastix_data->pastix_comm,
                 &TCR_flag, &TCR_status);
      if (TCR_flag)
        {
          /* Receive and add it */
          int     TCR_src = TCR_status.MPI_SOURCE;
          int     TCR_tag = TCR_status.MPI_TAG;
          FLOAT * TCR_vals_ptr;
          MPI_Recv(recv_data, send_nbr*blocksize, COMM_FLOAT,
                   TCR_src, TCR_tag,
                   solvers[id]->pastix_data->pastix_comm, &TCR_status);
          TCR_vals_ptr=recv_data;
          for (TCR_iter = 0; TCR_iter < send_nbr; TCR_iter++)
            {
              INT     TCR_index =
                sequence->recv_indexes[TCR_src][recv_cnt[TCR_src]];
              FLOAT * TCR_node  = &(solvers[id]->values[TCR_index]);
              INT     TCR_iterdof;
              recv_cnt[TCR_src]++;
              for (TCR_iterdof = 0; TCR_iterdof < blocksize;
                   TCR_iterdof++, TCR_node++, TCR_vals_ptr++)
                *TCR_node= sequence->fusion_local_entries(*TCR_node,
                                                          *TCR_vals_ptr);
            }
        }
    }
  return MURGE_SUCCESS;
}
#endif
/*
 * MURGE_AssemblyUseSequence
 *
 * Assembly the matrix using a stored sequence.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   id_seq  - Sequence ID.
 *   values  - Values to insert in the CSC.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *id_seq* or *values* are not valid.
 */
INTS MURGE_AssemblyUseSequence(INTS id, INTS id_seq, COEF * values)
{
#ifdef DISTRIBUTED
  murge_seq_t * sequence;
  INTL          iter;
  INTS          dof;
  FLOAT       * send_array = NULL;
  INT         * send_cnt   = NULL;
  MPI_Request * send_reqs  = NULL;
  MPI_Request * send_reqs2 = NULL;
  FLOAT       * recv_data  = NULL;
  INT         * recv_cnt   = NULL;
  int           send_nbr   = -1;
  int           blocksize;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  sequence = solvers[id]->sequences;
  while(sequence != NULL && sequence->ID != id_seq)
    sequence = sequence->next;

  if (sequence == NULL)
    {
      errorPrint("Sequence %d not found", id_seq);
      sequence = solvers[id]->sequences;
      return MURGE_ERR_PARAMETER;
    }

  if (values == NULL)
    {
      errorPrint("NULL value Pointer");
      return MURGE_ERR_PARAMETER;
    }

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("Graph has to be built before");
      return MURGE_ERR_ORDER;
    }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (solvers[id]->values == NULL)
    {
      MURGE_MEMALLOC(solvers[id]->values,
                     (solvers[id]->colptr[solvers[id]->n]-1)*dof*dof,
                     FLOAT);
      memset(solvers[id]->values, 0,
             (solvers[id]->colptr[solvers[id]->n]-1)*dof*dof*sizeof(FLOAT));
    }

  if (sequence->nodes == API_YES)
    {
      blocksize = dof*dof;
    }
  else
    {
      blocksize = 1;
    }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT)
    {
      INT coefnbr, recv_coefnbr;
      coefnbr = sequence->coefnbr;
      MPI_Allreduce(&coefnbr, &recv_coefnbr, 1, COMM_INT, MPI_SUM,
                    solvers[id]->pastix_data->pastix_comm);
      send_nbr = recv_coefnbr/100 + 1;

      MURGE_MEMALLOC(send_reqs, solvers[id]->pastix_data->procnbr, MPI_Request);
      MURGE_MEMALLOC(send_reqs2, solvers[id]->pastix_data->procnbr, MPI_Request);
      MURGE_MEMALLOC(send_array,
                     blocksize*send_nbr*solvers[id]->pastix_data->procnbr, FLOAT);
      MURGE_MEMALLOC(recv_data, blocksize*send_nbr, FLOAT);
      MURGE_MEMALLOC(recv_cnt,  solvers[id]->pastix_data->procnbr, INT);
      MURGE_MEMALLOC(send_cnt, solvers[id]->pastix_data->procnbr, INT);
      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        {
          recv_cnt[iter] = 0;
          send_cnt[iter] = 0;
        }
    }

  for (iter = 0; iter < sequence->coefnbr; iter++)
    {
      INTL index;
      FLOAT*node;
      INTS iterdof;

      if (sequence->mode != MURGE_ASSEMBLY_RESPECT)
        {
          try_complete_reception(id, sequence,
                                 MPI_ANY_SOURCE, send_nbr, blocksize,
                                 recv_data, recv_cnt);
        }

      index = sequence->indexes[iter];

      if (index>=0)
        {
          if (index > dof*dof*(solvers[id]->colptr[solvers[id]->n]-1))
            {
              return -1;
            }
          node = &(solvers[id]->values[index]);
          for (iterdof = 0; iterdof < blocksize; iterdof++, node++, values++)
            *node= sequence->fusion_local_entries(*node, *values);
        }
      else
        {
          if (sequence->mode == MURGE_ASSEMBLY_RESPECT)
            {
              errorPrint("Non local entry incompatible with MURGE_ASSEMBLY_RESPECT");
              return MURGE_ERR_PARAMETER;
            }

          if (send_cnt[-index-1] == send_nbr)
            {
              MPI_Status status;
              MPI_Wait(&(send_reqs[-index-1]), &status);
              send_cnt[-index-1] = 0;
            }

          /* Prepare to send, if send_buff is full send */
          node = &(send_array[blocksize*(send_cnt[-index-1]+send_nbr*(-index-1))]);
          for (iterdof = 0; iterdof < blocksize; iterdof++, node++, values++)
            *node= *values;

          send_cnt[-index-1]++;

          if (send_cnt[-index-1] == send_nbr)
            {
              MPI_Isend(&(send_array[send_nbr*(-index-1)]),
                        send_cnt[-index-1]*blocksize, COMM_FLOAT,
                        -index-1, TAG_VAL,
                        solvers[id]->pastix_data->pastix_comm,
                        &(send_reqs[-index-1]));
            }
        }
    }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT)
    {
      INT * done, done_cnt;
      MURGE_MEMALLOC(done, solvers[id]->pastix_data->procnbr, INT);
      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        done[iter] = API_NO;
      done_cnt = 0;
      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        {
          /* Wait last sends */
          if (send_cnt[iter] == send_nbr)
            {
              MPI_Status status;
              MPI_Wait(&(send_reqs[iter]), &status);
              send_cnt[iter] = 0;
            }

          /* Send last entries */
          MPI_Isend(&(send_cnt[iter]), 1, COMM_INT,
                    iter, TAG_SIZE,
                    solvers[id]->pastix_data->pastix_comm,
                    &(send_reqs[iter]));
          MPI_Isend(&(send_array[blocksize*send_nbr*(iter)]),
                    send_cnt[iter]*blocksize, COMM_FLOAT,
                    iter, TAG_VAL2,
                    solvers[id]->pastix_data->pastix_comm,
                    &(send_reqs2[iter]));
        }

      while (done_cnt < solvers[id]->pastix_data->procnbr)
        {
          for (iter =0; iter < solvers[id]->pastix_data->procnbr; iter++)
            {
              if (done[iter] == API_NO)
                {
                  try_complete_reception(id, sequence,
                                         iter, send_nbr, blocksize,
                                         recv_data, recv_cnt);

                  /* recv last count /entries */
                  MPI_Status status;
                  INT        cnt;
                  INT        iter2;
                  FLOAT     *myvals_ptr;
                  int       flag;
                  /* Receive and add it */
                  MPI_Iprobe(iter, TAG_SIZE,
                             solvers[id]->pastix_data->pastix_comm, &flag, &status);
                  if (flag)
                    {
                      MPI_Recv(&cnt, 1, COMM_INT,
                               iter, TAG_SIZE,
                               solvers[id]->pastix_data->pastix_comm, &status);
                      MPI_Recv(recv_data, cnt*blocksize, COMM_FLOAT,
                               iter, TAG_VAL2,
                               solvers[id]->pastix_data->pastix_comm, &status);
                      myvals_ptr = recv_data;

                      for (iter2 = 0; iter2 < cnt; iter2++)
                        {
                          INTS iterdof;
                          FLOAT*node;
                          INT index;

                          index = sequence->recv_indexes[iter][recv_cnt[iter]];
                          recv_cnt[status.MPI_SOURCE]++;

                          node  = &(solvers[id]->values[index]);
                          for (iterdof = 0; iterdof < dof*dof;
                               iterdof++, node++, myvals_ptr++)
                            *node= sequence->fusion_local_entries(*node, *myvals_ptr);
                        }
                      done[iter] = API_YES;
                      done_cnt++;
                    }
                }
            }
        }
      for (iter =0; iter < solvers[id]->pastix_data->procnbr; iter++)
        {
          MPI_Status status;
          MPI_Wait(&(send_reqs[iter]), &(status));
          MPI_Wait(&(send_reqs2[iter]), &(status));
        }

      memFree_null(done);
      memFree_null(send_reqs);
      memFree_null(send_reqs2);
      memFree_null(send_array);
      memFree_null(recv_data);
      memFree_null(recv_cnt);
      memFree_null(send_cnt);
    }
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_DUMP_MATRIX;
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 * Function: MURGE_AssemblyDeleteSequence
 *
 * Destroy an assembly sequence
 *
 *   id      - Solver instance identification number.
 *   id_seq  - Sequence ID.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *id_seq* is not valid.
 */
INTS MURGE_AssemblyDeleteSequence(INTS id, INTS id_seq)
{
#ifdef DISTRIBUTED
  murge_seq_t * sequence;
  murge_seq_t * psequence;
  INT iter;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  psequence = NULL;
  sequence  = solvers[id]->sequences;
  while(sequence != NULL && sequence->ID != id_seq)
    {
      psequence = sequence;
      sequence  = sequence->next;
    }

  if (sequence == NULL)
    {
      errorPrint("Sequence %d not found", id_seq);
      return MURGE_ERR_PARAMETER;
    }

  if (psequence != NULL)
    {
      psequence->next = sequence->next;
    }
  else
    {
      solvers[id]->sequences = sequence->next;
    }
  memFree_null(sequence->indexes);
  memFree_null(sequence->recv_nbr);
  if (NULL != sequence->recv_indexes)
    {
      for (iter = 0; iter < solvers[id]->pastix_data->procnbr; iter++)
        memFree_null(sequence->recv_indexes[iter]);
      memFree_null(sequence->recv_indexes);
    }
  memFree_null(sequence);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
/*
 * Function: MURGE_AssemblyBegin
 *
 * Check that preprocessing has been performed, if not performs it.
 *
 * Allocate ijv structure which will be used to store I,J,v[dof*dof].
 *
 * Parameters:
 *   op      - Operation to perform for coefficient which appear
 *             several time (see <MURGE_ASSEMBLY_OP>).
 *   op2     - Operation to perform when a coefficient is set by
 *             two different processors (see <MURGE_ASSEMBLY_OP>).
 *   mode    - Indicates if user ensure he will respect solvers distribution
 *             (see <MURGE_ASSEMBLY_MODE>).
 *   sym     - Indicates if user will give coefficient in a symmetric way
 *             (ie: only triangullar part) or not.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.
 */
INTS MURGE_AssemblyBegin(INTS id, INTL coefnbr, INTS op,
                         INTS op2, INTS mode, INTS sym)
{
#ifdef DISTRIBUTED
  INT iter;
  int dof;

  print_debug(DBG_MURGE, ">> MURGE_AssemblyBegin\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("Graph has to be built before");
      return MURGE_ERR_ORDER;
    }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (coefnbr < 0) {
    solvers[id]->dynamic = API_YES;
    coefnbr = -coefnbr;
  }
  else {
#ifdef MURGE_INSERT_DIRECTLY
    if (solvers[id]->colptr != NULL)
      {
        solvers[id]->dynamic = API_YES;
        coefnbr = 1+coefnbr/1000;
      }
    else
#endif /* MURGE_INSERT_DIRECTLY */
      {
        solvers[id]->dynamic = API_NO;
      }
  }
#ifdef MURGE_INSERT_DIRECTLY
  /* allocate values */
  if (solvers[id]->colptr != NULL && solvers[id]->values == NULL)
    {
      MURGE_MEMALLOC(solvers[id]->values, dof*dof*(solvers[id]->colptr[solvers[id]->n]-1), FLOAT);
      memset(solvers[id]->values, 0, (dof*dof*(solvers[id]->colptr[solvers[id]->n]-1))*sizeof(FLOAT));
    }
#endif /* MURGE_INSERT_DIRECTLY */
  solvers[id]->coefnbr     = coefnbr;
  solvers[id]->nodenbr     = coefnbr/(dof*dof);
  solvers[id]->tmpijv_node = NULL;
  solvers[id]->tmpijv      = NULL;
  solvers[id]->tmpv        = NULL;
  solvers[id]->tmpv_node   = NULL;
  solvers[id]->tmpijv_size = 0;
  solvers[id]->tmpijv_size_node = 0;

  if (dof == 1)
    {
      /* We will receive only coef by coef value */

      MURGE_MEMALLOC(solvers[id]->tmpv, coefnbr, FLOAT);
      memset(solvers[id]->tmpv, 0, coefnbr*sizeof(FLOAT));

      MURGE_MEMALLOC(solvers[id]->tmpijv,  coefnbr, ijv_t);
      memset(solvers[id]->tmpijv, 0, coefnbr*sizeof(ijv_t));

      solvers[id]->tmpijv_size = coefnbr;
      for (iter = 0; iter < coefnbr; iter++)
        {
          solvers[id]->tmpijv[iter].v  = &(solvers[id]->tmpv[iter]);
        }
    }
  else
    {
      /* we expect to receive node value arrays and
       only few isolated values */
      MURGE_MEMALLOC(solvers[id]->tmpv_node,
                     solvers[id]->nodenbr*dof*dof, FLOAT);
      memset(solvers[id]->tmpv_node, 0,
             solvers[id]->nodenbr*dof*dof*sizeof(FLOAT));

      MURGE_MEMALLOC(solvers[id]->tmpijv_node, solvers[id]->nodenbr, ijv_t);
      memset(solvers[id]->tmpijv_node, 0, solvers[id]->nodenbr*sizeof(ijv_t));
      solvers[id]->tmpijv_size_node = solvers[id]->nodenbr;
      for (iter = 0; iter < solvers[id]->nodenbr; iter++)
        {
          solvers[id]->tmpijv_node[iter].v  =
            &(solvers[id]->tmpv_node[iter*dof*dof]);
        }


      /* and for the isolated values */
      solvers[id]->tmpijv_size = coefnbr-solvers[id]->nodenbr*dof*dof+1;

      MURGE_MEMALLOC(solvers[id]->tmpv,
                     (solvers[id]->tmpijv_size), FLOAT);
      memset(solvers[id]->tmpv, 0,
             (solvers[id]->tmpijv_size)*sizeof(FLOAT));

      MURGE_MEMALLOC(solvers[id]->tmpijv,
                     (solvers[id]->tmpijv_size), ijv_t);
      memset(solvers[id]->tmpijv, 0,
             (solvers[id]->tmpijv_size)*sizeof(ijv_t));


      for (iter = 0; iter < solvers[id]->tmpijv_size; iter++)
        {
          solvers[id]->tmpijv[iter].v  = &(solvers[id]->tmpv[iter]);
        }
    }
  solvers[id]->edgenbr  = coefnbr;
  solvers[id]->cnt      = 0;
  solvers[id]->cnt_zero = 0;
  solvers[id]->cnt_node = 0;
  solvers[id]->mode     = mode;
  solvers[id]->op       = op;
  solvers[id]->op2      = op2;
  solvers[id]->sym      = sym;


#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif

  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_FACTO_OK);

  MURGE_STATE_TRUE(solvers[id]->state,  MURGE_MATR_BUILD);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 * Function:MURGE_AssemblySetValue
 *
 * Check that we are in an assembly section.
 *
 * Check that the number of coefficient entered will not
 * overpass the number of coefficient waited.
 *
 * Adds ROW, COL and value into the solvers[id]->tmpijv structure.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   ROW     - Global row number of the coefficient.
 *   COL     - Global column number of the coefficient.
 *   value   - value of the coefficient.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If we are not in an assembly section.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 */
INTS MURGE_AssemblySetValue     (INTS id, INTS ROW, INTS COL, COEF value)
{
#ifdef DISTRIBUTED
  int dof;
  FLOAT (*func)(FLOAT , FLOAT);

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHOOSE_FUNC(func, solvers[id]->op);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];


  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD)))
    {
      errorPrint("Need to call MURGE_AssemblyBegin first");
      return MURGE_ERR_ORDER;
    }

  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
    {
      COL += 1;
      ROW += 1;
    }

  if ((COL < 1) || (ROW < 1) ||
      (ROW > (((solvers[id]->N)*dof))) ||
      (COL > (((solvers[id]->N)*dof))))
    {
      errorPrint("COL (%ld) or ROW (%ld) is out of range [1-%ld]",
                 (long)COL, (long)ROW, (long)(solvers[id]->N*dof));
      return MURGE_ERR_PARAMETER;
    }

#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
  if (COL == MURGE_TRACE_COL && MURGE_TRACE_ROW == ROW)
    fprintf(stdout, "Setting A(%d,%d) <- %g\n",
            MURGE_TRACE_ROW, MURGE_TRACE_COL, value);
#endif
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD) && value == 0.0)
    {
      solvers[id]->cnt_zero ++;
      return MURGE_SUCCESS;
    }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&(solvers[id]->mutex_tmpmatrix));
#endif

  if ((solvers[id]->dynamic == API_NO) &&
      (solvers[id]->cnt + solvers[id]->cnt_zero +
       (solvers[id]->cnt_node )*dof*dof + 1 > solvers[id]->coefnbr))
    {
      errorPrint("Too many coef added in matrix building session (%ld > %ld)",
                 (long)(solvers[id]->cnt + solvers[id]->cnt_zero +
                        (solvers[id]->cnt_node )*dof*dof + 1),
                 (long)solvers[id]->coefnbr);
#ifdef MURGE_THREADSAFE
      pthread_mutex_unlock(&(solvers[id]->mutex_tmpmatrix));
#endif
      return MURGE_ERR_ORDER;
    }


  if (solvers[id]->cnt >= solvers[id]->tmpijv_size)
    {

      int iter;
      solvers[id]->tmpijv_size += solvers[id]->tmpijv_size/2 + 1;
      MURGE_REALLOC(solvers[id]->tmpv,
                    solvers[id]->tmpijv_size,
                    FLOAT);
      MURGE_REALLOC(solvers[id]->tmpijv,
                    solvers[id]->tmpijv_size,
                    ijv_t);
      for (iter = 0; iter < solvers[id]->tmpijv_size; iter++)
        {
          solvers[id]->tmpijv[iter].v = &(solvers[id]->tmpv[iter]);
        }

    }

  {
    INT node_col_glob;
    INT node_row_glob;
    INT node_col_loc;


    node_col_glob = (COL-1 - (COL-1)%dof)/dof;
    node_row_glob = (ROW-1 - (ROW-1)%dof)/dof;
    node_col_loc = solvers[id]->g2l[node_col_glob];
#ifdef MURGE_INSERT_DIRECTLY
    if ( solvers[id]->colptr != NULL && node_col_loc > 0 )
      {
        INT node_idx;
        int coef_idx;

        node_col_loc--;
        /* The column is local we add it into the local CSC */
        for (node_idx = solvers[id]->colptr[node_col_loc]-1;
             node_idx < solvers[id]->colptr[node_col_loc+1]-1;
             node_idx++)
          if (solvers[id]->rows[node_idx]-1 == node_row_glob) break;

        if (node_idx == solvers[id]->colptr[node_col_loc+1]-1)
          {
            if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD))
              {
                /* we will add it later */
                solvers[id]->tmpijv[solvers[id]->cnt].owner  =
                  solvers[id]->pastix_data->procnum;
                solvers[id]->tmpijv[solvers[id]->cnt].i      = ROW;
                solvers[id]->tmpijv[solvers[id]->cnt].j      = COL;
#ifdef MURGE_FOLLOW_INPUT_ORDER
                solvers[id]->tmpijv[solvers[id]->cnt].idx  = solvers[id]->cnt;
#endif
                solvers[id]->tmpijv[solvers[id]->cnt].v    =
                  &solvers[id]->tmpv[solvers[id]->cnt];

                solvers[id]->tmpijv[solvers[id]->cnt].v[0] = value;
                solvers[id]->cnt++;
              }
            else
              {
                errorPrint("ROW (%ld) not found in COL (%d)",
                           (long)ROW, (long)COL);
#ifdef MURGE_THREADSAFE
                pthread_mutex_unlock(&(solvers[id]->mutex_tmpmatrix));
#endif
                return MURGE_ERR_PARAMETER;
              }
          }
        else
          {
            /* we found the correct node */
            coef_idx = node_idx*dof*dof + ((COL-1)%dof)*dof + (ROW-1)%dof;
            solvers[id]->values[coef_idx] =
              func(solvers[id]->values[coef_idx], value);
          }
      }
    else
#endif /* MURGE_INSERT_DIRECTLY */
      {
        /* The column has to be sent to the correct CPU */

        if (node_col_loc < 0)
          {
            if (solvers[id]->mode == MURGE_ASSEMBLY_RESPECT)
              {
                errorPrint("Column %d is not local", COL);
                return MURGE_ERR_PARAMETER;
              }
            solvers[id]->tmpijv[solvers[id]->cnt].owner  = -node_col_loc;
          }
        else
          {
            solvers[id]->tmpijv[solvers[id]->cnt].owner  =
              solvers[id]->pastix_data->procnum;

          }
        solvers[id]->tmpijv[solvers[id]->cnt].i      = ROW;
        solvers[id]->tmpijv[solvers[id]->cnt].j      = COL;
#ifdef MURGE_FOLLOW_INPUT_ORDER
        solvers[id]->tmpijv[solvers[id]->cnt].idx  = solvers[id]->cnt;
#endif
        solvers[id]->tmpijv[solvers[id]->cnt].v    =
          &solvers[id]->tmpv[solvers[id]->cnt];

        solvers[id]->tmpijv[solvers[id]->cnt].v[0] = value;
        solvers[id]->cnt++;
      }
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&(solvers[id]->mutex_tmpmatrix));
#endif
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/*
 * Function:MURGE_AssemblySetNodeValues
 *
 * Check that we are in an assembly section.
 *
 * Check that the number of coefficient entered will not
 * overpass the number of coefficient waited.
 *
 * Adds ROW, COL and value into the solvers[id]->tmpijv structure.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   ROW     - Global row number of the coefficient.
 *   COL     - Global column number of the coefficient.
 *   values  - value of the coefficient.(dof^2 matrix)
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If we are not in an assembly section.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 */
INTS MURGE_AssemblySetNodeValues (INTS id, INTS ROW, INTS COL, COEF *values)
{
#ifdef DISTRIBUTED
  int dof;
  FLOAT (*func)(FLOAT , FLOAT);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (dof == 1)
    return MURGE_AssemblySetValue(id, ROW, COL, *values);

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHOOSE_FUNC(func, solvers[id]->op);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD)))
    {
      errorPrint("Need to call MURGE_AssemblyBegin first");
      return MURGE_ERR_ORDER;
    }
#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&solvers[id]->mutex_tmpmatrix);
#endif
  if (solvers[id]->dynamic == API_NO &&
      ( (solvers[id]->cnt_node+1)*dof*dof +
        solvers[id]->cnt + solvers[id]->cnt_zero > solvers[id]->coefnbr))
    {
      errorPrint("Too many coef added in matrix building session (%ld > %ld)",
                 (long)((solvers[id]->cnt_node+1)*dof*dof + solvers[id]->cnt),
                 (long)(solvers[id]->coefnbr));
#ifdef MURGE_THREADSAFE
      pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
      return MURGE_ERR_ORDER;
    }

  if (solvers[id]->cnt_node >= solvers[id]->tmpijv_size_node) {
    INT iterator;
    solvers[id]->tmpijv_size_node += solvers[id]->tmpijv_size_node/2 + 1;
    MURGE_REALLOC(solvers[id]->tmpijv_node,
                  solvers[id]->tmpijv_size_node, ijv_t);
    MURGE_REALLOC(solvers[id]->tmpv_node,
                  solvers[id]->tmpijv_size_node*dof*dof, FLOAT);
    for (iterator = 0; iterator < solvers[id]->tmpijv_size_node; iterator++)
      solvers[id]->tmpijv_node[iterator].v  =
        &(solvers[id]->tmpv_node[dof*dof*iterator]);
  }

  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
    {
      ROW += 1;
      COL += 1;
    }


  if ((COL < 1) || (ROW < 1) ||
      (ROW > ((solvers[id]->N))) ||
      (COL > ((solvers[id]->N))))
    {
      errorPrint("COL (%ld) or ROW (%ld) is out of range [1-%ld]",
                 (long)COL, (long)ROW, (long)(solvers[id]->N));
      return MURGE_ERR_PARAMETER;
    }
  {
    INT node_col_loc;

    node_col_loc = solvers[id]->g2l[COL-1];
#ifdef MURGE_INSERT_DIRECTLY
    if ( solvers[id]->colptr != NULL && node_col_loc > 0 )
      {
        INT node_idx;
        int coef_idx;

        node_col_loc--;
        /* The column is local we add it into the local CSC */
        for (node_idx = solvers[id]->colptr[node_col_loc]-1;
             node_idx < solvers[id]->colptr[node_col_loc+1]-1;
             node_idx++)
          if (solvers[id]->rows[node_idx] == ROW) break;

        if (node_idx == solvers[id]->colptr[node_col_loc+1]-1)
          {
            if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_ONLY_PROD))
              {
                solvers[id]->tmpijv_node[solvers[id]->cnt_node].owner =
                  solvers[id]->pastix_data->procnum;
                solvers[id]->tmpijv_node[solvers[id]->cnt_node].i = ROW;
                solvers[id]->tmpijv_node[solvers[id]->cnt_node].j = COL;
#ifdef MURGE_FOLLOW_INPUT_ORDER
                solvers[id]->tmpijv_node[solvers[id]->cnt_node].idx  = solvers[id]->cnt_node;
#endif
                memcpy(solvers[id]->tmpijv_node[solvers[id]->cnt_node].v,
                       values,
                       dof*dof*sizeof(COEF));
                solvers[id]->cnt_node++;
              }
            else
              {
                errorPrint("ROW (%ld) not found in COL (%d)",
                           (long)ROW, (long)COL);
#ifdef MURGE_THREADSAFE
                pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
                return MURGE_ERR_PARAMETER;
              }
          }
        /* we found the correct node */
        for ( coef_idx = 0;
              coef_idx < dof*dof;
              coef_idx++)
          solvers[id]->values[node_idx*dof*dof + coef_idx] =
            func(solvers[id]->values[node_idx*dof*dof+coef_idx],
                 values[coef_idx]);
      }
    else
#endif /* MURGE_INSERT_DIRECTLY */
      {

        if (node_col_loc < 0)
          {
            if (solvers[id]->mode == MURGE_ASSEMBLY_RESPECT)
              {
                errorPrint("Column %d is not local", COL);
                return MURGE_ERR_PARAMETER;
              }
            solvers[id]->tmpijv[solvers[id]->cnt].owner  = -node_col_loc;
          }
        else
          {
            solvers[id]->tmpijv[solvers[id]->cnt].owner  =
              solvers[id]->pastix_data->procnum;

          }
        solvers[id]->tmpijv_node[solvers[id]->cnt_node].i = ROW;
        solvers[id]->tmpijv_node[solvers[id]->cnt_node].j = COL;
#ifdef MURGE_FOLLOW_INPUT_ORDER
        solvers[id]->tmpijv_node[solvers[id]->cnt_node].idx  =
          solvers[id]->cnt_node;
#endif
        memcpy(solvers[id]->tmpijv_node[solvers[id]->cnt_node].v,
               values,
               dof*dof*sizeof(COEF));
        solvers[id]->cnt_node++;
      }
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

INTS MURGE_AssemblySetBlockValues(INTS id, INTS nROW, INTS *ROWlist,
                                  INTS nCOL, INTS *COLlist, COEF *values)
{
#ifdef DISTRIBUTED
  INT  iter;
  INT  iter2;
  INTS ret;
  INT  iterv;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD)))
    {
      errorPrint("Need to call MURGE_GraphBegin first");
      return MURGE_ERR_ORDER;
    }

  iterv = 0;
  for (iter = 0; iter < nROW; iter ++)
    {
      for (iter2 = 0; iter2 < nCOL; iter2++)
        {
          if (solvers[id]->pastix_data->iparm[IPARM_DOF_NBR] == 1)
            {
              ret = MURGE_AssemblySetValue(id,
                                           ROWlist[iter],
                                           COLlist[iter2],
                                           values[iterv]);
              if (MURGE_SUCCESS != ret)
                return ret;
              iterv++;
            }
          else
            {
              if (solvers[id]->pastix_data->iparm[IPARM_DOF_NBR] < 1)
                return MURGE_ERR_PARAMETER;

              ret = MURGE_AssemblySetNodeValues(id,
                                                ROWlist[iter],
                                                COLlist[iter2] ,
                                                &(values[iterv]));
              if (MURGE_SUCCESS != ret)
                return ret;

              iterv+=solvers[id]->pastix_data->iparm[IPARM_DOF_NBR]*
                solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
            }
        }
    }

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 Function: MURGE_AssemblyEnd

 We have on each proc a part of the matrix in
 two structure, one containing nodes to add
 to the CSCd the other containing simple values.

 We send all data to his owner:
 - We sort our data structures (IJV structures)
 using the "owner" attribute.
 - We send non local data to other processors.

 We merge all data in the node structure.
 - We receive Data and merge node structure with simple
 values one.
 - We look for each coef in node structure, if present we modify the node, if
 not, we search in the CSCd and directly modify it. Else we construct
 a new node and add it.

 We Add this structure to the local CSCd.

 */

INTS MURGE_AssemblyEnd(INTS id)
{
#ifdef DISTRIBUTED
  INT          size;
  INT          index;
  INT          procnum;
  INT          iter, iter2, iter3;
  INT          startidx           = 0;
  INT          startidx2          = 0;
  INT          firstlast[2];
  INT         *edgenbr            = NULL;
  INT         *edgenbr_node       = NULL;
  INT         *edgenbr_recv       = NULL;
  INT         *edgenbr_recv_node  = NULL;
  INT          tmpijvsize         = 0;
  ijv_t       *tmpijv             = NULL;
  INT          tmpijvsize_node    = 0;
  ijv_t       *tmpijv_node        = NULL;
  ijv_t       *ijvptr             = NULL;
  MPI_Status   status;
  MPI_Request *requests_size      = NULL;
  MPI_Request *requests_ijv       = NULL;
  MPI_Request *requests_val       = NULL;
  MPI_Request *requests_size_node = NULL;
  MPI_Request *requests_ijv_node  = NULL;
  MPI_Request *requests_val_node  = NULL;
  int          baseval;
  int          dof;
  INT         *tmpcolptr          = NULL;
  INT         *tmprows            = NULL;
  FLOAT       *tmpvalues          = NULL;
  FLOAT       *valuesptr          = NULL;
  INT         *tmpcolptr2         = NULL;
  INT         *tmprows2           = NULL;
  FLOAT       *tmpvalues2         = NULL;
  FLOAT       *tmpvalues3         = NULL;
  FLOAT       *tmpvalue_mdof      = NULL;
  void        *sortptr[2];
  FLOAT        tmpvalue;
  FLOAT (*func)(FLOAT , FLOAT);
#ifdef CENTRALISED
  INT         *total_nodelist;
#endif
#ifdef MURGE_TIME
  Clock       clock;
  double      time1;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  print_debug(DBG_MURGE, ">> MURGE_AssemblyEnd\n");
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD)))
    {
      errorPrint("Need to call MURGE_GraphBegin first");
      return MURGE_ERR_ORDER;
    }


  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (solvers[id]->dynamic == API_NO &&
      solvers[id]->cnt_node*dof*dof + solvers[id]->cnt + solvers[id]->cnt_zero != solvers[id]->edgenbr)
    {
      errorPrint("Wrong number of entries  (%ld != %ld) ",
                 (long)solvers[id]->cnt_node*dof*dof + solvers[id]->cnt + solvers[id]->cnt_zero ,
                 (long)solvers[id]->edgenbr);
      return MURGE_ERR_ORDER;
    }
  if (solvers[id]->cnt_zero != 0) {
    if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
      {
        fprintf(stdout,
                "%ld (%.2g %%) zero entries were skipped on proc %ld\n",
                (long)solvers[id]->cnt_zero,
                (double)((double)100.0*((double)solvers[id]->cnt_zero)/
                         ((double)(solvers[id]->cnt_node*dof*dof +
                                   solvers[id]->cnt + solvers[id]->cnt_zero))),
                (long)solvers[id]->pastix_data->procnum);
      }
    else
      {
        if (solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
          {
            INT nz_glob;
            INT zeros_glob;
            INT nz = solvers[id]->colptr[solvers[id]->n]-1;
            MPI_Reduce( &(solvers[id]->cnt_zero), &zeros_glob,
                        1, COMM_INT,
                        MPI_SUM, 0, solvers[id]->pastix_data->pastix_comm);
            MPI_Reduce( &nz, &nz_glob,
                        1, COMM_INT,
                        MPI_SUM, 0, solvers[id]->pastix_data->pastix_comm);
            if (solvers[id]->pastix_data->procnum == 0)
              {
                fprintf(stdout,
                        "%ld zero entries were skipped"
                        " (from %ld (%.3g%%))\n",
                        (long int)zeros_glob,
                        (long int)nz_glob,
                        100.0*((double)(zeros_glob)/
                               ((double)(nz_glob))));
              }

          }
      }
  }
#ifdef CENTRALISED
  MURGE_MEMALLOC(total_nodelist, solvers[id]->N, INT);
  for (iter = 0; iter < solvers[id]->N; iter++)
    total_nodelist[iter] = iter+1;
#endif

  CHOOSE_FUNC(func, solvers[id]->op);

  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof     = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];


  /*
   * Distributing tmpijv structure onto processors.
   *
   * tmpijv contains the (i,j,v), value by value, added on
   * local processor, but this (i,j,v) can be not local.
   * We fill the owner field of the ijv_t structure to be abble
   * to sort the array following the owner and then (i,j).
   *
   */
  for (iter = 0; iter < solvers[id]->cnt; iter++)
    {
      index = (solvers[id]->tmpijv[iter].j-1 -
               (solvers[id]->tmpijv[iter].j-1)%dof)/dof;
      if ( solvers[id]->g2l[index] > 0)
        {
          solvers[id]->tmpijv[iter].owner = solvers[id]->pastix_data->procnum;
        }
      else
        {
          solvers[id]->tmpijv[iter].owner =
            -solvers[id]->g2l[index];
        }
    }

  /*
   * The same is done on tmpijv_node which contains node
   * (i,j,v[dof*dof]) structure.
   *
   */
  for (iter = 0; iter < solvers[id]->cnt_node; iter++)
    {
      index = solvers[id]->tmpijv_node[iter].j-1;
      if ( solvers[id]->g2l[index] > 0)
        {
          solvers[id]->tmpijv_node[iter].owner =
            solvers[id]->pastix_data->procnum;
        }
      else
        {
          solvers[id]->tmpijv_node[iter].owner =
            -solvers[id]->g2l[index];
        }
    }

  CLOCK_STOP;
#ifdef MURGE_TIME
  time1 = CLOCK_GET;
#endif

  /*****************************************************************************
   * The arrays are sorted following the owner and the column
   */

#ifdef CENTRALISED

  sortptr[0] = solvers[id]->tmpijv;
  dof = 1;
  sortptr[1] = &dof;
  MurgeTmpijvSort(sortptr, solvers[id]->cnt);
  sortptr[0] = solvers[id]->tmpijv_node;
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  sortptr[1] = &dof;
  MurgeTmpijvSort(sortptr, solvers[id]->cnt_node);
#else
  sortptr[0] = solvers[id]->tmpijv;
  dof = 1;
  sortptr[1] = &dof;
  if (solvers[id]->mode != MURGE_ASSEMBLY_RESPECT)
    MurgeTmpijvOwnSort(sortptr, solvers[id]->cnt);
  else
    MurgeTmpijvSort(sortptr, solvers[id]->cnt);
  sortptr[0] = solvers[id]->tmpijv_node;
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  sortptr[1] = &dof;
  if (solvers[id]->mode != MURGE_ASSEMBLY_RESPECT)
    MurgeTmpijvOwnSort(sortptr, solvers[id]->cnt_node);
  else
    MurgeTmpijvSort(sortptr, solvers[id]->cnt_node);
#endif

#ifdef MURGE_TIME
  CLOCK_STOP;
  fprintf(stdout, " > > MURGE_AssemblyEnd tmpijv and tmpijv_node sorted in %.3lg seconds (id=%d)\n", CLOCK_GET - time1, id);
  time1 = CLOCK_GET;
#endif

  for (iter = 0; iter < solvers[id]->cnt_node; iter++)
    {
      ASSERT(solvers[id]->tmpijv_node[iter].v ==
             &(solvers[id]->tmpv_node[iter*dof*dof]),
             MOD_MURGE);
    }

  for (iter = 0; iter < solvers[id]->cnt; iter++)
    {
      ASSERT(solvers[id]->tmpijv[iter].v ==
             &(solvers[id]->tmpv[iter]),
             MOD_MURGE);
    }

#ifdef MURGE_TIME
  CLOCK_STOP;
  fprintf(stdout, " > > MURGE_AssemblyEnd double entries removed in %.3lg seconds (id=%d)\n", CLOCK_GET - time1, id);
  time1 = CLOCK_GET;
#endif

  /*****************************************************************************
   * Fusion double entries
   *
   * Depend on the fact that tmpijv and tmpv are sorted on
   * tmpijv[iter].j and tmpijv[iter].j
   *
   */
  /* For node entries */
  if (dof > 1)
    {

      MURGE_MEMALLOC(tmpvalue_mdof, dof*dof, FLOAT);

      index = 0;
      for (iter2 = 0; iter2 < dof*dof; iter2++)
        tmpvalue_mdof[iter2] = 0.0;
      for (iter = 0; iter < solvers[id]->cnt_node; iter++)
        {
          if (solvers[id]->tmpijv_node[index].i ==
              solvers[id]->tmpijv_node[iter].i &&
              solvers[id]->tmpijv_node[index].j ==
              solvers[id]->tmpijv_node[iter].j)
            {
              for (iter2 = 0; iter2 < dof*dof; iter2++)
                tmpvalue_mdof[iter2] =
                  func(tmpvalue_mdof[iter2],
                       solvers[id]->tmpijv_node[iter].v[iter2]);


            }
          else
            {
              for (iter2 = 0; iter2 < dof*dof; iter2++)
                solvers[id]->tmpv_node[index*dof*dof+iter2] =
                  tmpvalue_mdof[iter2];

              solvers[id]->tmpijv_node[index].v =
                &(solvers[id]->tmpv_node[index*dof*dof]);

              index++;
              solvers[id]->tmpijv_node[index].i =
                solvers[id]->tmpijv_node[iter].i;
              solvers[id]->tmpijv_node[index].j =
                solvers[id]->tmpijv_node[iter].j;
              solvers[id]->tmpijv_node[index].owner =
                solvers[id]->tmpijv_node[iter].owner;

              for (iter2 = 0; iter2 < dof*dof; iter2++)
                tmpvalue_mdof[iter2] =
                  solvers[id]->tmpv_node[iter*dof*dof+iter2];
            }
        }
      if (solvers[id]->cnt_node >0 )
        {
          for (iter2 = 0; iter2 < dof*dof; iter2++)
            solvers[id]->tmpv_node[index*dof*dof+iter2] = tmpvalue_mdof[iter2];

          solvers[id]->tmpijv_node[index].v =
            &(solvers[id]->tmpv_node[index*dof*dof]);
          index++;

          solvers[id]->cnt_node = index;
        }
      memFree_null(tmpvalue_mdof);
    }

  /* For simple entries */
  if (solvers[id]->cnt > 0)
    {
      index = 0;
      tmpvalue = solvers[id]->tmpijv[0].v[0];
      for (iter = 1; iter < solvers[id]->cnt; iter++)
        {
          /* If we find same (i,j) we add value to tmpvalue
           Else, we store tmpvalue and we change reference index and
           start a new tmpvalue
           */
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
          if (solvers[id]->tmpijv[index].i == MURGE_TRACE_ROW &&
              solvers[id]->tmpijv[index].j == MURGE_TRACE_COL)
            fprintf(stdout, "%d,%d %g + %g = %g (%d)\n",
                    MURGE_TRACE_ROW, MURGE_TRACE_COL, tmpvalue, solvers[id]->tmpijv[iter].v[0],
                    func(tmpvalue, solvers[id]->tmpijv[iter].v[0]), solvers[id]->tmpijv[iter].owner);
#endif
          if (solvers[id]->tmpijv[index].i == solvers[id]->tmpijv[iter].i &&
              solvers[id]->tmpijv[index].j == solvers[id]->tmpijv[iter].j)
            {
              tmpvalue = func(tmpvalue, solvers[id]->tmpijv[iter].v[0]);
            }
          else
            {
              solvers[id]->tmpv[index] = tmpvalue;
              solvers[id]->tmpijv[index].v = &(solvers[id]->tmpv[index]);
              index++;
              solvers[id]->tmpijv[index].i = solvers[id]->tmpijv[iter].i;
              solvers[id]->tmpijv[index].j = solvers[id]->tmpijv[iter].j;
              solvers[id]->tmpijv[index].owner = solvers[id]->tmpijv[iter].owner;
              tmpvalue = solvers[id]->tmpijv[iter].v[0];
            }
        }
      if (solvers[id]->cnt >0 )
        {
          solvers[id]->tmpv[index] = tmpvalue;
          solvers[id]->tmpijv[index].v = &(solvers[id]->tmpv[index]);
          index++;
          solvers[id]->cnt = index;
        }
    }
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
  for (iter = 1; iter < solvers[id]->cnt; iter++)
    {
      if (solvers[id]->tmpijv[iter].i == MURGE_TRACE_ROW &&
          solvers[id]->tmpijv[iter].j == MURGE_TRACE_COL)
        fprintf(stdout, "A(%d,%d) = %g (%d, %d)\n", MURGE_TRACE_ROW, MURGE_TRACE_COL,
                solvers[id]->tmpijv[iter].v[0], iter,  solvers[id]->tmpijv[iter].owner);
    }
#endif
#ifdef MURGE_TIME
  CLOCK_STOP;
  fprintf(stdout, " > > MURGE_AssemblyEnd double entries removed in %.3lg seconds (id=%d)\n", CLOCK_GET - time1, id);
  time1 = CLOCK_GET;
#endif

  /*****************************************************************************
   * If there is node stucture (dof > 1) we send locally added nodes
   * to its owner.
   *
   * We start computing the index and size of the part of tmpijv_node/tmpv
   * to send and we send it.
   *
   * For local data we just store starting index and number of
   * nodes which is sufficient.
   *
   */
  if (dof > 1)
    {
      MURGE_MEMALLOC(requests_val_node,  solvers[id]->pastix_data->procnbr, MPI_Request);
      MURGE_MEMALLOC(requests_ijv_node,  solvers[id]->pastix_data->procnbr, MPI_Request);
      MURGE_MEMALLOC(requests_size_node, solvers[id]->pastix_data->procnbr, MPI_Request);
      MURGE_MEMALLOC(edgenbr_node, solvers[id]->pastix_data->procnbr, INT);

      firstlast[0] = 0;
      firstlast[1] = 0;

      for (procnum = 0; procnum <solvers[id]->pastix_data->procnbr; procnum++)
        {
          /* On envoi a chacun les infos ijv par noeud*/
          while(firstlast[1]< solvers[id]->cnt_node &&
                solvers[id]->tmpijv_node[firstlast[1]].owner ==
                procnum)
            firstlast[1]++;

#ifdef CENTRALISED
          /* On envoi tout */
          firstlast[0] = 0;
          firstlast[1] = solvers[id]->cnt_node;
#endif

          edgenbr_node[procnum] = firstlast[1] - firstlast[0];

          if (procnum != solvers[id]->pastix_data->procnum)
            {
#ifdef MURGE_DEBUG_COMM
              fprintf(stdout, "%d sends %d nodes to %d\n", solvers[id]->pastix_data->procnum, edgenbr_node[procnum], procnum);
#endif
              MPI_Isend(&edgenbr_node[procnum],
                        1, COMM_INT, procnum,
                        TAG_SIZE, solvers[id]->pastix_data->pastix_comm,
                        &(requests_size_node[procnum]));

              if (edgenbr_node[procnum] > 0)
                {
                  MPI_Isend(&(solvers[id]->tmpijv_node[firstlast[0]]),
                            edgenbr_node[procnum]*sizeof(ijv_t),
                            MPI_BYTE, procnum,
                            TAG_IJV, solvers[id]->pastix_data->pastix_comm,
                            &(requests_ijv_node[procnum]));

                  /* tmpv is in the same order as tmpijv, sorted by owner */
                  MPI_Isend(solvers[id]->tmpijv_node[firstlast[0]].v,
                            dof*dof*edgenbr_node[procnum]*sizeof(FLOAT),
                            MPI_BYTE, procnum,
                            TAG_VAL, solvers[id]->pastix_data->pastix_comm,
                            &(requests_val_node[procnum]));
                }
            }
          else
            {
              startidx = firstlast[0];
            }
          firstlast[0] = firstlast[1];
        }
    }

  /*
   * In all cases we can have data to send in
   * tmpijv (value by value structure)
   *
   * As for tmpijv_node, we send it.
   *
   */
  firstlast[0] = 0;
  firstlast[1] = 0;
  MURGE_MEMALLOC(requests_val,  solvers[id]->pastix_data->procnbr, MPI_Request);
  MURGE_MEMALLOC(requests_ijv,  solvers[id]->pastix_data->procnbr, MPI_Request);
  MURGE_MEMALLOC(requests_size, solvers[id]->pastix_data->procnbr, MPI_Request);
  PRINT_DEALLOC(requests_size, __FILE__, __LINE__);
  MURGE_MEMALLOC(edgenbr,       solvers[id]->pastix_data->procnbr, INT);

  for (procnum = 0; procnum <solvers[id]->pastix_data->procnbr; procnum++)
    {
      /* On envoi a chacun les infos ijv coef par coef*/
      while(firstlast[1]< solvers[id]->cnt &&
            solvers[id]->tmpijv[firstlast[1]].owner ==
            procnum)
        firstlast[1]++;

#ifdef CENTRALISED
      /* On envoi tout */
      firstlast[0] = 0;
      firstlast[1] = solvers[id]->cnt;
#endif

      edgenbr[procnum] = firstlast[1] - firstlast[0];

      if (procnum != solvers[id]->pastix_data->procnum)
        {
#ifdef MURGE_DEBUG_COMM
          fprintf(stdout, "%d sends %d coefs to %d\n", solvers[id]->pastix_data->procnum, edgenbr[procnum], procnum);
#endif
          MPI_Isend(&edgenbr[procnum],
                    1, COMM_INT, procnum,
                    TAG_SIZE2, solvers[id]->pastix_data->pastix_comm,
                    &(requests_size[procnum]));

          if (edgenbr[procnum] > 0)
            {
              MPI_Isend(&(solvers[id]->tmpijv[firstlast[0]]),
                        edgenbr[procnum]*sizeof(ijv_t),
                        MPI_BYTE, procnum,
                        TAG_IJV2, solvers[id]->pastix_data->pastix_comm,
                        &(requests_ijv[procnum]));
              /* tmpv is in the same order as tmpijv, sorted by owner */
              MPI_Isend(solvers[id]->tmpijv[firstlast[0]].v,
                        edgenbr[procnum]*sizeof(FLOAT),
                        MPI_BYTE, procnum,
                        TAG_VAL2, solvers[id]->pastix_data->pastix_comm,
                        &(requests_val[procnum]));
            }
        }
      else
        {
          startidx2 = firstlast[0];
        }
      firstlast[0] = firstlast[1];
    }
  PRINT_DEALLOC(requests_size, __FILE__, __LINE__);
#ifdef MURGE_TIME
  CLOCK_STOP;
  fprintf(stdout, " > > MURGE_AssemblyEnd send data in %.3lg seconds (id=%d)\n", CLOCK_GET - time1, id);
  time1 = CLOCK_GET;
#endif
  MURGE_MEMALLOC(edgenbr_recv, solvers[id]->pastix_data->procnbr, INT);
  MURGE_MEMALLOC(edgenbr_recv_node, solvers[id]->pastix_data->procnbr, INT);

  PRINT_DEALLOC(requests_size, __FILE__, __LINE__);
  /*****************************************************************************
   * For each node, we receive node and values structures, we merge it into
   * the node structure and then we build a CSCd that we add to the existing
   * CSCd if it exists.
   *
   */
  for (procnum = 0; procnum <solvers[id]->pastix_data->procnbr; procnum++)
    {

      /*
       * Receiving node structure if it exists
       */
      edgenbr_recv_node[procnum] = 0;
      if (dof > 1)
        {
          if (procnum != solvers[id]->pastix_data->procnum)
            {
              MPI_Recv(&edgenbr_recv_node[procnum],
                       1, COMM_INT, procnum,
                       TAG_SIZE, solvers[id]->pastix_data->pastix_comm,
                       &status);
#ifdef MURGE_DEBUG_COMM
              fprintf(stdout, "%d received %d coefs from %d\n", solvers[id]->pastix_data->procnum, edgenbr_recv_node[procnum], procnum);
#endif
            }
          else
            {
              edgenbr_recv_node[procnum] = edgenbr_node[procnum];
            }

          if (edgenbr_recv_node[procnum] > 0)
            {
              /* We keep the same buffer if it fits in last one */
              if (tmpijvsize_node < edgenbr_recv_node[procnum])
                {
                  if (NULL != tmpijv_node)
                    memFree_null(tmpijv_node);
                  if (edgenbr_recv_node[procnum] > 0)
                    MURGE_MEMALLOC(tmpijv_node,
                                   edgenbr_recv_node[procnum],
                                   ijv_t);


                  if (NULL != tmpvalues)
                    memFree_null(tmpvalues);
                  MURGE_MEMALLOC(tmpvalues, edgenbr_recv_node[procnum]*dof*dof,
                                 FLOAT);

                  tmpijvsize_node = edgenbr_recv_node[procnum];
                }

              if (procnum != solvers[id]->pastix_data->procnum)
                {
                  MPI_Recv(tmpijv_node,
                           edgenbr_recv_node[procnum]*sizeof(ijv_t),
                           MPI_BYTE, procnum, TAG_IJV,
                           solvers[id]->pastix_data->pastix_comm,
                           &status);

                  MPI_Recv(tmpvalues,
                           dof*dof*edgenbr_recv_node[procnum]*sizeof(FLOAT),
                           MPI_BYTE, procnum, TAG_VAL,
                           solvers[id]->pastix_data->pastix_comm,
                           &status);
                }
              else
                {
                  memcpy(tmpijv_node,
                         &(solvers[id]->tmpijv_node[startidx]),
                         edgenbr_recv_node[procnum]*sizeof(ijv_t));
                  memcpy(tmpvalues,
                         solvers[id]->tmpijv_node[startidx].v,
                         dof*dof*edgenbr_recv_node[procnum]*sizeof(FLOAT));
                }

              for (iter = 0; iter < edgenbr_recv_node[procnum]; iter++)
                tmpijv_node[iter].v = &(tmpvalues[iter*dof*dof]);

              ijvptr    = tmpijv_node;
              valuesptr = tmpvalues;
            }
          else
            {
              ijvptr    = NULL;
              valuesptr = NULL;
            }
        }


      if (procnum != solvers[id]->pastix_data->procnum)
        {
          /*
           * Receiving value by value structure
           */
          MPI_Recv(&edgenbr_recv[procnum],
                   1, COMM_INT, procnum,
                   TAG_SIZE2, solvers[id]->pastix_data->pastix_comm,
                   &status);
        }
      else
        {
          edgenbr_recv[procnum] = edgenbr[procnum];
        }

#ifdef MURGE_DEBUG_COMM
      fprintf(stdout, "%d received %d nodes from %d\n", solvers[id]->pastix_data->procnum, edgenbr_recv[procnum], procnum);
#endif

      if (edgenbr_recv[procnum] > 0)
        {
          /* We keep the same buffer if it fits in last one */

          if (tmpijvsize < edgenbr_recv[procnum])
            {
              if (NULL != tmpijv)
                memFree_null(tmpijv);
              MURGE_MEMALLOC(tmpijv, edgenbr_recv[procnum], ijv_t);


              if (NULL != tmpvalues3)
                memFree_null(tmpvalues3);
              MURGE_MEMALLOC(tmpvalues3, edgenbr_recv[procnum],
                             FLOAT);

              tmpijvsize = edgenbr_recv[procnum];
            }

          if (procnum != solvers[id]->pastix_data->procnum &&
              edgenbr_recv[procnum] != 0)
            {
              MPI_Recv(tmpijv,
                       edgenbr_recv[procnum]*sizeof(ijv_t),
                       MPI_BYTE, procnum, TAG_IJV2,
                       solvers[id]->pastix_data->pastix_comm,
                       &status);
              MPI_Recv(tmpvalues3,
                       edgenbr_recv[procnum]*sizeof(FLOAT),
                       MPI_BYTE, procnum, TAG_VAL2,
                       solvers[id]->pastix_data->pastix_comm,
                       &status);
#ifdef MURGE_DEBUG_COMM
              fprintf(stdout, "%d received %d from %d\n", solvers[id]->pastix_data->procnum, edgenbr_recv_node[procnum], procnum);
#endif
            }
          else
            {
              memcpy(tmpijv,
                     &(solvers[id]->tmpijv[startidx2]),
                     edgenbr_recv[procnum]*sizeof(ijv_t));
              memcpy(tmpvalues3,
                     solvers[id]->tmpijv[startidx2].v,
                     edgenbr_recv[procnum]*sizeof(FLOAT));
            }

          for (iter = 0; iter < edgenbr_recv[procnum]; iter++)
            tmpijv[iter].v = &(tmpvalues3[iter]);

          /*
           * Merging ijv value by value structure into node structure.
           */
          if (dof > 1)
            {
              for (iter = 0; iter < edgenbr_recv[procnum]; iter++)
                {
                  /*
                   Search in ijvptr if the node containing
                   (tmpijv[iter].i, tmpijv[iter].j) exists

                   If not, search in CSCd and use the CSCd one as a
                   base for a new ijv node.

                   If not, set a new node.
                   */

                  /* compress row number, 0 based */
                  INT c_row = (tmpijv[iter].i - 1 - (tmpijv[iter].i-1)%dof)/dof;
                  /* compress column number, 0 based */
                  INT c_col = (tmpijv[iter].j - 1 - (tmpijv[iter].j-1)%dof)/dof;
                  /* index of the value to set in the node*/
                  INT intra_nd_idx = (tmpijv[iter].i-1)%dof +
                    (tmpijv[iter].j-1)%dof*dof;

                  for (iter2 = 0; iter2  < edgenbr_recv_node[procnum]; iter2++)
                    if (tmpijv_node[iter2].i-1 == c_row &&
                        tmpijv_node[iter2].j-1 == c_col)
                      break;

                  if (iter2 <  edgenbr_recv_node[procnum])
                    {
                      /* found in tmpijv, at index iter2 */
                      tmpijv_node[iter2].v[intra_nd_idx] =
                        func( tmpijv_node[iter2].v[intra_nd_idx],
                              tmpijv[iter].v[0]);
                    }
                  else
                    {
                      if (solvers[id]->colptr != NULL)
                        {
                          for (iter2 = 0; iter2 < solvers[id]->n; iter2++)
                            if (solvers[id]->l2g[iter2] - 1 == c_col)
                              break;

                          if (iter2 != solvers[id]->n)
                            {

                              /* column found, column iter2 */
                              for (iter3 = solvers[id]->colptr[iter2]-1;
                                   iter3 < solvers[id]->colptr[iter2+1]-1;
                                   iter3++)
                                {
                                  if (solvers[id]->rows[iter3] - 1== c_row)
                                    break;
                                }

                              if (iter3 !=  solvers[id]->colptr[iter2+1]-1)
                                {
                                  /* line found, line index is iter3
                                   *
                                   * Add a new coef to the tmpijv_node
                                   * based on the values of the CSCd
                                   */

                                  /* extend tmpijv_node and valuesptr if needed */
                                  EXTEND_NODE_LIST;

                                  memcpy(tmpijv_node[edgenbr_recv_node[procnum]].v,
                                         &(solvers[id]->values[iter3*dof*dof]),
                                         dof*dof*sizeof(FLOAT));

                                  /* found in tmpijv */
                                  solvers[id]->values[iter3*dof*dof+ intra_nd_idx] =
                                    func(solvers[id]->values[iter3*dof*dof+
                                                             intra_nd_idx],
                                         tmpijv[iter].v[0]);
                                }
                              else /* iter3 !=  solvers[id]->colptr[iter2+1] -1 */
                                {
                                  COPY_ELEMENT_TO_NODE;
                                }
                            } /* iter2 != solvers[id]->n */
                          else
                            {
                              errorPrint("Shouldn't go there\n");
                              return MURGE_ERR_PARAMETER;
                            }
                        } /* solvers[id]->colptr != NULL */
                      else
                        {
                          /* Add a coefficient to the node structure */
                          COPY_ELEMENT_TO_NODE;
                        }
                    }
                }
            }
          else
            {
              ijvptr = tmpijv;
              valuesptr = tmpvalues3;
              edgenbr_recv_node[procnum] = edgenbr_recv[procnum];
            }
        }
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
      for (iter = 1; iter < edgenbr_recv_node[procnum]; iter++)
        {
          if (ijvptr[iter].i == MURGE_TRACE_ROW &&
              ijvptr[iter].j == MURGE_TRACE_COL)
            fprintf(stdout, ".A(%d,%d) = %g\n", MURGE_TRACE_ROW, MURGE_TRACE_COL, ijvptr[iter].v[0]);
        }
#endif
      /*
       * We build the local CSCd and add it to the
       * existing one if it exists.
       */
#ifdef CENTRALISED
      solvers[id]->n = solvers[id]->N;
      memFree_null(solvers[id]->l2g);
      MURGE_MEMALLOC(solvers[id]->l2g,solvers[id]->N, INT);
      for (iter = 0; iter < solvers[id]->N; iter++)
        solvers[id]->l2g[iter] = iter+1;
#endif

      if (edgenbr_recv_node[procnum] > 0)
        {
          /* On ajoute a la cscd locale */

          /* Building colptr */
          MURGE_MEMALLOC(tmpcolptr, solvers[id]->n+1, INT);

          for (iter=0; iter<(solvers[id]->n); iter++)
            tmpcolptr[iter] = 1;

          iter2 = 1;
          for (iter=0; iter<(solvers[id]->n); iter++)
            {
              tmpcolptr[iter] = iter2;
              while ((iter2 < edgenbr_recv_node[procnum] + 1) &&
                     (ijvptr[iter2-1].j == solvers[id]->l2g[iter]))
                {
                  iter2++;
                }
            }

          tmpcolptr[solvers[id]->n] = iter2;

          size = tmpcolptr[solvers[id]->n]-1;
          if (size > 0) /* size < 0 Shouldn't happen */
            MURGE_MEMALLOC(tmprows, size, INT);


          for (iter=0; iter<size; iter++)
            {
              tmprows[iter] = ijvptr[iter].i;
            }


          if (solvers[id]->colptr != NULL)
            {
              /* if it is not the first add */
              cscd_addlocal_int(solvers[id]->n,
                                solvers[id]->colptr,
                                solvers[id]->rows,
                                solvers[id]->values,
                                solvers[id]->l2g,
                                solvers[id]->n, tmpcolptr, tmprows, valuesptr,
                                solvers[id]->l2g,
                                &solvers[id]->n, &tmpcolptr2, &tmprows2,
                                &tmpvalues2, func, dof, API_YES);
              memFree_null(solvers[id]->colptr);
              memFree_null(solvers[id]->rows);
              memFree_null(solvers[id]->values);

              memFree_null(tmpcolptr);
              memFree_null(tmprows);
              solvers[id]->colptr = tmpcolptr2;
              solvers[id]->rows   = tmprows2;
              solvers[id]->values = tmpvalues2;
            }
          else
            {
              /* For first add, we juste copy the CSCd */
              solvers[id]->colptr = tmpcolptr;
              solvers[id]->rows   = tmprows;

              size = solvers[id]->colptr[solvers[id]->n]-1;
              MURGE_MEMALLOC(solvers[id]->values, size*dof*dof, FLOAT);
              memcpy(solvers[id]->values,
                     valuesptr,
                     size*dof*dof*sizeof(FLOAT));
              tmpcolptr = NULL;
              tmprows   = NULL;
              valuesptr = NULL;
            }
        }
    }

#ifdef MURGE_TIME
  CLOCK_STOP;
  fprintf(stdout, " > > MURGE_AssemblyEnd receved and merged data in %.3lg seconds (id=%d)\n", CLOCK_GET - time1, id);
  time1 = CLOCK_GET;
#endif
  memFree_null(tmpvalues3);
  memFree_null(edgenbr_recv_node);
  if (NULL != tmpijv)
    memFree_null(tmpijv);
  if (NULL != tmpvalues)
    memFree_null(tmpvalues);

  memFree_null(edgenbr_recv);

  /* Just wait for send to be finished */
  for (procnum = 0; procnum <solvers[id]->pastix_data->procnbr; procnum++)
    {
      if (procnum != solvers[id]->pastix_data->procnum)
        {
          if (dof > 1)
            {
              MPI_Wait(&(requests_size_node[procnum]),
                       &status);
              if (edgenbr_node[procnum] > 0)
                {
                  MPI_Wait(&(requests_ijv_node[procnum]),
                           &status);
                  memFree_null(tmpijv_node);

                  MPI_Wait(&(requests_val_node[procnum]),
                           &status);
                }
            }
          MPI_Wait(&(requests_size[procnum]),
                   &status);
          if (edgenbr[procnum] > 0)
            {
              MPI_Wait(&(requests_ijv[procnum]),
                       &status);
              memFree_null(tmpijv);

              MPI_Wait(&(requests_val[procnum]),
                       &status);
            }

        }
    }
  PRINT_DEALLOC(requests_size, __FILE__, __LINE__);
#ifdef MURGE_TIME
  CLOCK_STOP;
  fprintf(stdout, " > > MURGE_AssemblyEnd terminating send in %.3lg seconds (id=%d)\n", CLOCK_GET - time1, id);
  time1 = CLOCK_GET;
#endif

  if (NULL != tmpijv_node)
    memFree_null(tmpijv_node);

  if (solvers[id]->colptr == NULL )
    {
      MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, INT);
      for (iter = 0; iter < solvers[id]->n+1; iter++)
        solvers[id]->colptr[iter]=1;
    }

#ifdef CENTRALISED
  memFree_null(total_nodelist);
#endif

  MURGE_DUMP_MATRIX;

  memFree_null(requests_size);
  memFree_null(requests_ijv);
  memFree_null(requests_val);
  memFree_null(edgenbr);
  if (dof > 1) {
    memFree_null(requests_size_node);
    memFree_null(requests_ijv_node);
    memFree_null(requests_val_node);
    memFree_null(edgenbr_node);
    memFree_null(solvers[id]->tmpijv_node);
    memFree_null(solvers[id]->tmpv_node);
  }
  memFree_null(solvers[id]->tmpijv);
  memFree_null(solvers[id]->tmpv);
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
  fprintf(stdout, "..A(%d,%d) = %g\n", MURGE_TRACE_ROW, MURGE_TRACE_COL,
          solvers[id]->values[solvers[id]->colptr[MURGE_TRACE_COL-1]-1+MURGE_TRACE_ROW-1]);
#endif

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&solvers[id]->mutex_tmpmatrix);
#endif

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);

  MURGE_STATE_FALSE(solvers[id]->state, MURGE_MATR_BUILD);
#ifdef MURGE_TIME
  CLOCK_STOP;
  fprintf(stdout, " > MURGE_AssemblyEnd computed in %.3lg seconds (id=%d)\n", (double)CLOCK_GET, id);
#endif

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

INTS MURGE_MatrixReset(INTS id){
#ifdef DISTRIBUTED
  INT nbcoef;
  INT dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (solvers[id]->values != NULL)
    {
      nbcoef = solvers[id]->colptr[solvers[id]->n]-1;
      dof    = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

      memset(solvers[id]->values, 0, nbcoef*dof*dof*sizeof(FLOAT));
    }
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

INTS MURGE_MatrixGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS,
                           COEF *values, INTS root, INTS op, INTS sym)
{
#ifdef DISTRIBUTED
  INT  dof;
  INT  coefnbr;
  INT  iter;
  INT  iter2;
  INTS ret;
  int baseval;

  CHECK_SOLVER_ID(id);
  if (solvers[id]->pastix_data->procnum == 0)
    {
      errorPrintW("MURGE_MatrixGlobalCSR is not optimal with PaStiX, try using MURGE_MatrixGlobalCSC instead");
    }
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("Graph has to be built before");
      return MURGE_ERR_ORDER;
    }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);
  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof     = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  /* Si on a une matrice symetrique autant faire un MURGE_MatrixGlobalCSC */
  if (solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES)
    return  MURGE_MatrixGlobalCSC(id, N, rowptr, COLS, values, root, op, sym);

  coefnbr = rowptr[N] - baseval;

  if (solvers[id]->pastix_data->procnum == root ||
      (root == -1 && solvers[id]->pastix_data->procnum == 0))
    {
      ret = MURGE_AssemblyBegin(id, coefnbr,  op, op,MURGE_ASSEMBLY_FOOL,
                                solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;
      for (iter = 0; iter < N; iter++)
        {
          for (iter2 = rowptr[iter]; iter2 < rowptr[iter+1]; iter2++)
            {
              if (dof == 1)
                {

                  ret = MURGE_AssemblySetValue(id,
                                               iter+baseval,
                                               COLS[iter2-baseval],
                                               values[iter2-baseval]);
                  if (MURGE_SUCCESS != ret)
                    return ret;
                }
              else
                {
                  ret = MURGE_AssemblySetNodeValues(id,
                                                    iter+baseval,
                                                    COLS[iter2-baseval],
                                                    &(values[(iter2-baseval)*
                                                             dof*dof]));
                  if (MURGE_SUCCESS != ret)
                    return ret;
                }
            }
        }
    }
  else
    {
      ret = MURGE_AssemblyBegin(id, 0, op, op,
                                MURGE_ASSEMBLY_FOOL,
                                solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;
    }

  if (MURGE_SUCCESS != (ret = MURGE_AssemblyEnd(id)))
    return ret;

  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
/*
 Function: MURGE_MatrixGlobalCSC

 Give a CSC on one processor to PaStiX.
 */
INTS MURGE_MatrixGlobalCSC(INTS id, INTS N, INTL *COLPTR, INTS *ROWS,
                           COEF *values, INTS root, INTS op, INTS sym)
{
#ifdef DISTRIBUTED
  INT          *l2g = NULL;
  INT           procnum;
  INT           localn;
  INT          *tmpcolptr;
  INT          *tmprows;
  FLOAT        *tmpvalues;
  INT          *tmpcolptr2;
  INT          *tmprows2;
  INT           tmpn;
  MPI_Status    status;
  INT           iter;
  int           dof;
  FLOAT (*func)(FLOAT , FLOAT);

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);


  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("Graph has to be built before");
      return MURGE_ERR_ORDER;
    }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

  CHOOSE_FUNC(func, op);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  /* Si tout le monde est racine */
  if (root == -1 || root == solvers[id]->pastix_data->procnum)
    {

      if (sizeof(INTS) != sizeof(INT))
        {
          MURGE_MEMALLOC(tmprows2, (COLPTR[N]-1), INT);
          for (iter = 0; iter <  COLPTR[N]-1; iter++)
            {
              tmprows2[iter] = (INT)ROWS[iter];
            }
        }
      else
        {
          tmprows2 = (INT*)ROWS;
        }

      if (sizeof(INTL) != sizeof(INT))
        {
          MURGE_MEMALLOC(tmpcolptr2, N+1, INT);
          for (iter = 0; iter <  N+1; iter++)
            {
              tmpcolptr2[iter] = (INT)COLPTR[iter];
            }
        }
      else
        {
          tmpcolptr2 = (INT*)COLPTR;
        }


      if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
        {
          tmprows2 --;
          values   --;
        }

      MURGE_MEMALLOC(l2g, N, INT);

      for (iter = 0; iter < N; iter++)
        l2g[iter] = iter+1;

      /* colptr must be allocated for cscd_addlocal_int */
      if (NULL == solvers[id]->colptr)
        {
          MURGE_MEMALLOC(solvers[id]->colptr, solvers[id]->n+1, INT);
          for (iter = 0; iter < solvers[id]->n+1; iter++)
            {
              solvers[id]->colptr[iter] = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
            }
        }

      cscd_addlocal_int(solvers[id]->n,
                        solvers[id]->colptr,
                        solvers[id]->rows,
                        solvers[id]->values,
                        solvers[id]->l2g,
                        N, tmpcolptr2, tmprows2, (FLOAT*)values, l2g,
                        &tmpn, &tmpcolptr, &tmprows, &tmpvalues, func,
                        dof,  API_YES);

      if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
        {
          tmprows2 ++;
          values   ++;
        }

      if (root == -1 && sizeof(INTS) != sizeof(INT))
        {
          memFree_null(tmprows2);
        }
      if (root == -1 && sizeof(INTL) != sizeof(INT))
        {
          memFree_null(tmpcolptr2);
        }
      memFree_null(solvers[id]->colptr);
      memFree_null(solvers[id]->rows);
      memFree_null(solvers[id]->values);

      solvers[id]->colptr = tmpcolptr;
      solvers[id]->rows   = tmprows;
      solvers[id]->values = tmpvalues;

    }
  if (root != -1)
    {
      /* si on est le processeur racine */
      if (root == solvers[id]->pastix_data->procnum)
        {

          /* Pour chaque processeur, on calcule
           la CSCD a ajouter puis on l'envoi.
           */
          for (procnum = 0; procnum < solvers[id]->pastix_data->procnbr; procnum++)
            {
              if (procnum != solvers[id]->pastix_data->procnum)
                {
                  MPI_Recv(&localn, 1,  COMM_INT, procnum, TAG_SIZE,
                           solvers[id]->pastix_data->pastix_comm, &status);
                  MPI_Recv(l2g, localn, COMM_INT, procnum, TAG_L2G,
                           solvers[id]->pastix_data->pastix_comm, &status);

                  MURGE_MEMALLOC(tmpcolptr, localn + 1, INT);

                  for (iter = 0; iter < localn+1; iter++)
                    {
                      tmpcolptr[iter] = 1;
                    }
                  if (sizeof(INTS) != sizeof(INT))
                    {
                      MURGE_MEMALLOC(tmprows2, (COLPTR[N]-1), INT);
                      for (iter = 0; iter <  COLPTR[N]-1; iter++)
                        {
                          tmprows2[iter] = (INT)ROWS[iter];
                        }
                    }
                  else
                    {
                      tmprows2 = (INT*)ROWS;
                    }

                  if (sizeof(INTL) != sizeof(INT))
                    {
                      MURGE_MEMALLOC(tmpcolptr2, N+1, INT);
                      for (iter = 0; iter <  N+1; iter++)
                        {
                          tmpcolptr2[iter] = (INT)COLPTR[iter];
                        }
                    }
                  else
                    {
                      tmpcolptr2 = (INT*)ROWS;
                    }

                  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
                    {
                      tmprows2 --;
                      values   --;
                    }
                  cscd_addlocal_int(localn,
                                    tmpcolptr,
                                    NULL,
                                    NULL,
                                    l2g,
                                    N, tmpcolptr2, tmprows2, (FLOAT*)values, l2g,
                                    &localn,
                                    &tmpcolptr,
                                    &tmprows,
                                    &tmpvalues, func, dof, API_YES);
                  if (solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0)
                    {
                      tmprows2 ++;
                      values   ++;
                    }

                  /* On envoi tmpcolptr, tmprows, tmpvalues */
                  MPI_Send(tmpcolptr,
                           localn+1, COMM_INT, procnum,
                           TAG_COL, solvers[id]->pastix_data->pastix_comm);
                  memFree_null(tmpcolptr);

                  MPI_Send(tmprows,
                           tmpcolptr[localn - 1], COMM_INT, procnum,
                           TAG_ROW, solvers[id]->pastix_data->pastix_comm);
                  memFree_null(tmprows);

                  MPI_Send(tmpvalues,
                           tmpcolptr[localn - 1], COMM_FLOAT, procnum,
                           TAG_VAL, solvers[id]->pastix_data->pastix_comm);
                  memFree_null(tmpvalues);

                  if (sizeof(INTS) != sizeof(INT))
                    {
                      memFree_null(tmprows2);
                    }
                  if (sizeof(INTL) != sizeof(INT))
                    {
                      memFree_null(tmpcolptr2);
                    }
                }
              else
                {
                  /* La CSCd local a déjà été traitée */
                }

            }
          if (sizeof(INTS) != sizeof(INT))
            {
              memFree_null(tmprows2);
            }
          if (sizeof(INTL) != sizeof(INT))
            {
              memFree_null(tmpcolptr2);
            }
        }
      else
        {
          /* Si on est pas la racine, on recoit de la racine la CSCd a ajouter
           et on l'ajoute
           */

          MPI_Send(&solvers[id]->n,
                   1, COMM_INT, root,
                   TAG_SIZE, solvers[id]->pastix_data->pastix_comm);
          localn = solvers[id]->n;
          MPI_Send(solvers[id]->l2g,
                   solvers[id]->n,
                   COMM_INT, root,
                   TAG_L2G, solvers[id]->pastix_data->pastix_comm);

          MPI_Recv(solvers[id]->colptr,
                   localn+1, COMM_INT, root,
                   TAG_COL, solvers[id]->pastix_data->pastix_comm, &status);

          MPI_Recv(solvers[id]->rows,
                   tmpcolptr[localn - 1], COMM_INT, root,
                   TAG_ROW, solvers[id]->pastix_data->pastix_comm, &status);

          MPI_Recv(solvers[id]->values,
                   tmpcolptr[localn - 1], COMM_FLOAT, root,
                   TAG_VAL, solvers[id]->pastix_data->pastix_comm, &status);


        }
    }

  MURGE_DUMP_MATRIX;


  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 Function: MURGE_MatrixGlobalIJV

 Add the given global Compress Sparse Column matrix to the matrix.

 Parameters:
 id      - Solver instance identification number.
 N       - Number of edges.
 NNZ     - Number of non zeros.
 ROWS    - Global row number array.
 COLS    - Global column number array.
 values  - values array.
 root    - Root processor for MPI communications.
 op      - Operation to perform if a coefficient appear twice
 (see <MURGE_ASSEMBLY_OP>).
 sym     - Indicates if user will give coefficient in a symmetric way
 (ie: only triangullar part) or not.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range,
 if *root*, *op*, *ROWS* or *COLS* are not valid.

 Fortran interface:
 >
 > SUBROUTINE MURGE_MATRIXGLOBALIJV(ID, N, NNZ, ROWS, COLS, VALUES, &
 >                                & ROOT, OP, SYM, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, OP, SYM, N
 >   INTL,               INTENT(IN)  :: NNZ
 >   INTS, DIMENSION(0), INTENT(IN)  :: ROWS, COLS
 >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_MATRIXGLOBALIJV
 */
INTS MURGE_MatrixGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROWS, INTS *COLS,
                           COEF *values, INTS root, INTS op, INTS sym)
{
#ifdef DISTRIBUTED
  int        baseval;             /* User baseval                         */
  INT        iter;                /* Iterators                            */
  INTS       ret;                 /* Return value                         */
  int        dof;                 /* Number of degree of freedom          */

  CHECK_SOLVER_ID(id);
  if (solvers[id]->pastix_data->procnum == 0)
    errorPrintW("MURGE_MatrixGlobalIJV is not optimal with PaStiX, try using MURGE_MatrixGlobalCSC instead");
  CHECK_SOLVER_PARAM(id);


  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_BUILD) )
    {
      errorPrint("Do not call MURGE_GraphBegin before");
      return MURGE_ERR_ORDER;
    }

  if ( MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_MATR_BUILD) )
    {
      errorPrint("Do not call MURGE_AssemblyBegin before");
      return MURGE_ERR_ORDER;
    }
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_GRAPH_OK)))
    {
      errorPrint("Graph has to be built before");
      return MURGE_ERR_ORDER;
    }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);
  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof     = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (solvers[id]->pastix_data->procnum == root ||
      (root == -1 && solvers[id]->pastix_data->procnum == 0))
    {
      ret = MURGE_AssemblyBegin(id, NNZ, op, op, MURGE_ASSEMBLY_FOOL,
                                solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;

      for (iter = 0; iter < NNZ; iter++)
        {
          if (dof == 1)
            {

              ret = MURGE_AssemblySetValue(id,
                                           ROWS[iter],
                                           COLS[iter],
                                           values[iter]);
              if (MURGE_SUCCESS != ret)
                return ret;
            }
          else
            {
              ret = MURGE_AssemblySetNodeValues(id,
                                                ROWS[iter],
                                                COLS[iter],
                                                &(values[iter*dof*dof]));
              if (MURGE_SUCCESS != ret)
                return ret;
            }
        }
    }
  else
    {
      ret = MURGE_AssemblyBegin(id, 0, op, op,
                                MURGE_ASSEMBLY_FOOL,
                                solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;
    }

  if (MURGE_SUCCESS != (ret = MURGE_AssemblyEnd(id)))
    return ret;

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_VALUES_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/*******************************************************************************
 * Group: Filling the right-hand-side member
 */


INTS MURGE_SetGlobalRHS(INTS id, COEF *b, INTS root, INTS op)
{
#ifdef DISTRIBUTED
  INT        iter;
  INT        procnum;
  INT        localn;
  INT        lastn = 0;
  INT       *l2g   = NULL;
  FLOAT     *tmpb  = NULL;
  FLOAT    (*func)(FLOAT , FLOAT);
  MPI_Status status;
  int        dof;
  int        iterdof;
  COEF      *b_recv;
  int        allocated = API_NO;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

  CHOOSE_FUNC(func, op);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (root != -1) {
    /* Broadcast and then run algorithm for root = -1 */
    /* TODO : fix mistake and remove BCast            */
    MURGE_MEMALLOC(b_recv, solvers[id]->N*dof, FLOAT);
    if (root ==  solvers[id]->pastix_data->procnum) {
      memcpy(b_recv, b, solvers[id]->N*dof*sizeof(FLOAT));
    }
    MPI_Bcast(b_recv, solvers[id]->N*dof, COMM_FLOAT, root,
              solvers[id]->pastix_data->pastix_comm);
    root = -1;
    allocated = API_YES;
  }
  else {
    b_recv = b;
  }

  if (root == -1)
    {
      /* Si tous le monde à la racine */
      if (NULL == solvers[id]->b)
        {
          MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof, FLOAT);
          memset(solvers[id]->b, 0, solvers[id]->n*dof*sizeof(FLOAT));
        }
      for (iter = 0; iter < solvers[id]->n; iter++)
        {
          for (iterdof = 0; iterdof < dof; iterdof++)
            {

              solvers[id]->b[iter*dof +iterdof] =
                func(solvers[id]->b[iter*dof +iterdof],
                     b_recv[(solvers[id]->l2g[iter]-1)*dof+iterdof]);
            }
        }
    }
  else
    {
      /* Sinon, on recupère sur la racine tous les loc2globs
       On construit et on envoi tous les right-hand-side locaux
       */
      if (root == solvers[id]->pastix_data->procnum)
        {
          lastn = 0;
          for (procnum = 0;
               procnum < solvers[id]->pastix_data->procnbr;
               procnum++)
            {
              if (procnum != solvers[id]->pastix_data->procnum)
                {

                  MPI_Recv(&localn, 1, COMM_INT, procnum, TAG_SIZE,
                           solvers[id]->pastix_data->pastix_comm, &status);

                  if (lastn < localn)
                    {
                      if (NULL != l2g)
                        memFree_null(l2g);
                      MURGE_MEMALLOC(l2g, localn, INT);
                      if (tmpb != NULL)
                        memFree_null(tmpb);
                      MURGE_MEMALLOC(tmpb, localn*dof, FLOAT);
                      lastn = localn;
                    }
                  MPI_Recv(l2g, localn, COMM_INT, procnum, TAG_L2G,
                           solvers[id]->pastix_data->pastix_comm, &status);

                  for (iter = 0; iter < localn; iter++)
                    {
                      for (iterdof = 0; iterdof < dof; iterdof++)
                        {
                          tmpb[iter*dof+iterdof]= b[(l2g[iter]-1)*dof+iterdof];
                        }
                    }
                  MPI_Send(tmpb,
                           localn*dof, COMM_FLOAT, procnum,
                           TAG_VAL, solvers[id]->pastix_data->pastix_comm);
                }
            }
          if (NULL != tmpb)
            memFree_null(tmpb);
          if (NULL != l2g)
            memFree_null(l2g);

          /* Le processeur racine construit son bout de RHS */
          if (NULL == solvers[id]->b)
            {
              MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof, FLOAT);
              memset(solvers[id]->b, 0, solvers[id]->n*dof*sizeof(FLOAT));
            }

          for (iter = 0; iter < solvers[id]->n; iter++)
            {
              for (iterdof = 0; iterdof < dof; iterdof++)
                {
                  solvers[id]->b[iter*dof+iterdof] =
                    func(solvers[id]->b[iter*dof+iterdof],
                         b[(solvers[id]->l2g[iter]-1)*dof+iterdof]);
                }
            }
        }
      else
        {
          /* Sur les procs non racine on recoit simplement le RHS a ajouter */
          MPI_Send(&solvers[id]->n,
                   1, COMM_INT, root,
                   TAG_SIZE, solvers[id]->pastix_data->pastix_comm);
          MPI_Send(solvers[id]->l2g,
                   solvers[id]->n,
                   COMM_INT, root,
                   TAG_L2G, solvers[id]->pastix_data->pastix_comm);

          MURGE_MEMALLOC(tmpb, solvers[id]->n, FLOAT);
          if (NULL == solvers[id]->b)
            {
              solvers[id]->b = tmpb;
            }

          MPI_Recv(tmpb,
                   solvers[id]->n*dof, COMM_FLOAT, root,
                   TAG_VAL, solvers[id]->pastix_data->pastix_comm, &status);

          if (tmpb != solvers[id]->b)
            {
              for (iter = 0; iter < solvers[id]->n*dof; iter++)
                {
                  solvers[id]->b[iter] =
                    func(solvers[id]->b[iter],
                         tmpb[iter]);
                }
              memFree_null(tmpb);
            }
        }
    }
  if (allocated == API_YES) memFree_null(b_recv);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_RHS_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

INTS MURGE_SetLocalRHS (INTS id, COEF *b, INTS op, INTS op2)
{
#ifdef DISTRIBUTED
  INT        iter;
  FLOAT    (*func)(FLOAT , FLOAT) = NULL;
  int        dof;
#ifdef CENTRALISED
  INT        nodenbr;
  INT       *intern_nodelist;
  FLOAT     *tmpb;
#endif

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  CHOOSE_FUNC(func, op);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

#ifdef CENTRALISED
  nodenbr = pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));
  MURGE_MEMALLOC(intern_nodelist, nodenbr, INT);
  if (NO_ERR != ( pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                                         intern_nodelist)))
    return MURGE_ERR_SOLVER;

  if (NULL == solvers[id]->b)
    {
      MURGE_MEMALLOC(solvers[id]->b, solvers[id]->N*dof, FLOAT);
      memset(solvers[id]->b, 0, solvers[id]->N*dof*sizeof(FLOAT));
    }
  MURGE_MEMALLOC(tmpb, solvers[id]->N*dof /* solvers[id]->n*dof */, FLOAT);
  memset(tmpb, 0, solvers[id]->N*dof*sizeof(FLOAT));

  for (iter = 0; iter < nodenbr*dof; iter++)
    {
      tmpb[(intern_nodelist[(iter-iter%dof)/dof]-1)*dof+iter%dof] =
        func(tmpb[(intern_nodelist[(iter-iter%dof)/dof]-1)*dof+iter%dof],
             b[iter]);
    }
  MPI_Allreduce(tmpb, solvers[id]->b, solvers[id]->N*dof, COMM_FLOAT, MPI_SUM,
                solvers[id]->pastix_data->pastix_comm);

  memFree_null(tmpb);
  memFree_null(intern_nodelist);

#else /* CENTRALISED */
  if (NULL == solvers[id]->b)
    {
      MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof , FLOAT);
      memset(solvers[id]->b, 0, solvers[id]->n*dof*sizeof(FLOAT));
    }

  for (iter = 0; iter < solvers[id]->n*dof; iter++)
    {
      solvers[id]->b[iter] = func(solvers[id]->b[iter],
                                  b[iter]);
    }
#endif /* CENTRALISED */

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_RHS_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

INTS MURGE_SetRHS      (INTS id, INTS n, INTS *coefsidx, COEF *b, INTS op,
                        INTS op2, INTS mode)
{
#ifdef DISTRIBUTED
  INT             iter;
  FLOAT         (*func)(FLOAT , FLOAT) = NULL;
  INT             index;
  int             baseval;
  int             dof;
  int             iterdof;
  pastix_data_t * pastix_data = solvers[id]->pastix_data;
  INTS            procnbr     = pastix_data->procnbr;


  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);
  baseval = pastix_data->iparm[IPARM_BASEVAL];
  dof = pastix_data->iparm[IPARM_DOF_NBR];
  CHOOSE_FUNC(func, op);

  if (mode == MURGE_ASSEMBLY_FOOL)
    {
      INTS            coefs_rcv_size;
      INTS         *  coefnbr;
      INTS         ** coefs_idx;
      COEF         ** coefs_vals;
      MPI_Request  *  request_cfnbr;
      MPI_Request  *  request_cfidx;
      MPI_Request  *  request_cfvals;

      MURGE_MEMALLOC(coefnbr, procnbr, INTS);
      for (iter = 0; iter <procnbr; iter++)
        coefnbr[iter] = 0;

      /* Count the entries to send to each processor */
      for (iter = 0; iter <n; iter++)
        {
          INT procnum;
          if (solvers[id]->g2l[coefsidx[iter]-1] > 0)
            procnum = pastix_data->procnum;
          else
            procnum = -solvers[id]->g2l[coefsidx[iter]-1];

          coefnbr[procnum]++;
        }
      MURGE_MEMALLOC(coefs_idx,  procnbr, INTS*);
      MURGE_MEMALLOC(coefs_vals, procnbr, COEF*);

      for (iter = 0; iter < procnbr; iter++)
        {
          MURGE_MEMALLOC(coefs_idx[iter],  coefnbr[iter],     INTS);
          MURGE_MEMALLOC(coefs_vals[iter], coefnbr[iter]*dof, COEF);
          coefnbr[iter] = 0;
        }
      /* Prepare the arrays to send to each processors */
      for (iter = 0; iter <n; iter++)
        {
          INT procnum;
          if (solvers[id]->g2l[coefsidx[iter]-1] > 0)
            procnum = pastix_data->procnum;
          else
            procnum = -solvers[id]->g2l[coefsidx[iter]-1];

          coefs_idx[procnum][coefnbr[procnum]] = coefsidx[iter];
          for (iterdof = 0; iterdof < dof; iterdof++)
            {
              coefs_vals[procnum][coefnbr[procnum]*dof+iterdof] =
                b[iter*dof+iterdof];
            }

          coefnbr[procnum]++;
        }

      MURGE_MEMALLOC(request_cfnbr,  procnbr, MPI_Request);
      MURGE_MEMALLOC(request_cfidx,  procnbr, MPI_Request);
      MURGE_MEMALLOC(request_cfvals, procnbr, MPI_Request);

      /* Send the data to the processors */
      for (iter = 0; iter < procnbr; iter++)
        {
          MPI_Isend(&(coefnbr[iter]),    1,                 MPI_INTS,
                    iter, TAG_SIZE, pastix_data->pastix_comm,
                    &(request_cfnbr[iter]));
          if (coefnbr[iter] > 0) {
            MPI_Isend(coefs_idx[iter],     coefnbr[iter],     MPI_INTS,
                      iter, TAG_ROW, pastix_data->pastix_comm,
                      &(request_cfidx[iter]));
            MPI_Isend(coefs_vals[iter],    coefnbr[iter]*dof, MURGE_MPI_COEF,
                      iter, TAG_VAL, pastix_data->pastix_comm,
                      &(request_cfvals[iter]));
          }
        }

      /* receive the data and run MPI_SetRHS with MURGE_ASSEMBLY_RESPECT */
      coefs_rcv_size = 0;
      for (iter = 0; iter < procnbr; iter++)
        {
          INTS         coefnbr_rcv;
          INTS       * coefs_idx_rcv  = NULL;
          COEF       * coefs_vals_rcv = NULL;
          MPI_Status   status;

          MPI_Recv(&coefnbr_rcv, 1, MPI_INTS, iter, TAG_SIZE,
                   pastix_data->pastix_comm, &status);

          if (coefnbr_rcv > 0)
            {
              if (coefnbr_rcv > coefs_rcv_size)
                {
                  if (coefs_rcv_size != 0)
                    {
                      memFree_null(coefs_idx_rcv);
                      memFree_null(coefs_vals_rcv);
                    }
                  MURGE_MEMALLOC(coefs_idx_rcv,  coefnbr_rcv,         INTS);
                  MURGE_MEMALLOC(coefs_vals_rcv, coefnbr_rcv*dof,     COEF);
                  coefs_rcv_size = coefnbr_rcv;
                }
              MPI_Recv(coefs_idx_rcv, coefnbr_rcv,     MPI_INTS,
                       iter, TAG_ROW, pastix_data->pastix_comm, &status);
              MPI_Recv(coefs_vals_rcv,coefnbr_rcv*dof, MURGE_MPI_COEF,
                       iter, TAG_VAL, pastix_data->pastix_comm, &status);

              MURGE_SetRHS(id, coefnbr_rcv, coefs_idx_rcv, coefs_vals_rcv, op,
                           op2, MURGE_ASSEMBLY_RESPECT);
            }
          if (iter == procnbr-1 && coefs_rcv_size != 0)
            {
              memFree_null(coefs_idx_rcv);
              memFree_null(coefs_vals_rcv);
            }
        }

      /* Now we clean it all */
      for (iter = 0; iter < procnbr; iter++)
        {
          MPI_Status status;

          MPI_Wait(&(request_cfnbr[iter]), &status);
          MPI_Wait(&(request_cfidx[iter]), &status);
          MPI_Wait(&(request_cfvals[iter]), &status);
        }
      memFree_null(request_cfnbr);
      memFree_null(request_cfidx);
      memFree_null(request_cfvals);

      for (iter = 0; iter <procnbr; iter++)
        {
          memFree_null(coefs_idx[iter]);
          memFree_null(coefs_vals[iter]);
        }
      memFree_null(coefs_idx);
      memFree_null(coefs_vals);
      memFree_null(coefnbr);
    }
  else
    {
      if (!MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODELST_OK))
        {
          errorPrint("Assembly can be respected only if user asked for node list");
          return MURGE_ERR_ORDER;
        }
      for (iter = 0; iter < n; iter++)
        {
          index = coefsidx[iter]- baseval;
          index = solvers[id]->g2l[index] - baseval;

          for (iterdof = 0; iterdof < dof; iterdof++)
            {
              solvers[id]->b[index*dof+iterdof] =
                func(solvers[id]->b[index*dof+iterdof],
                     b[iter*dof+iterdof]);
            }
        }
    }

  MURGE_STATE_TRUE(solvers[id]->state, MURGE_RHS_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

INTS MURGE_RHSReset(INTS id){
#ifdef DISTRIBUTED
  INT iter, dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (NULL == solvers[id]->b)
    MURGE_MEMALLOC(solvers[id]->b, solvers[id]->n*dof, FLOAT);

  for (iter = 0; iter < solvers[id]->n*dof; iter++)
    solvers[id]->b[iter] =0;
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_RHS_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/*******************************************************************************
 * Group: Getting the solution
 */

INTS MURGE_GetGlobalSolution(INTS id, COEF *x, INTS root)
{
#ifdef DISTRIBUTED
  int    dof;
#ifndef CENTRALISED
  FLOAT *tmpx = NULL;
  INT    iter;
  int    iterdof;
#endif

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  CHECK_FACT(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if ((root == solvers[id]->pastix_data->procnum ||
       root == -1)
      && NULL == x)
    return MURGE_ERR_PARAMETER;
#ifdef CENTRALISED
  if ((root == solvers[id]->pastix_data->procnum ||
       root == -1))
    memcpy(x, solvers[id]->b, solvers[id]->N*dof*sizeof(FLOAT));
#else
  MURGE_MEMALLOC(tmpx, solvers[id]->N*dof, FLOAT);

  for (iter = 0; iter < solvers[id]->N*dof; iter ++)
    tmpx[iter] = 0;
  for (iter = 0; iter < solvers[id]->n; iter ++)
    {
      for (iterdof = 0; iterdof < dof; iterdof++)
        {
          tmpx[(solvers[id]->l2g[iter]-1)*dof+iterdof] = solvers[id]->b[iter*dof+iterdof];
        }
    }

  if (root == -1)
    {
      MPI_Allreduce(tmpx, x, solvers[id]->N*dof, COMM_FLOAT, COMM_SUM,
                    solvers[id]->pastix_data->pastix_comm);
    }
  else
    {
      MPI_Reduce(tmpx, x, solvers[id]->N*dof, COMM_FLOAT, COMM_SUM, root,
                 solvers[id]->pastix_data->pastix_comm);
    }

  memFree_null(tmpx);
#endif
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
INTS MURGE_GetLocalSolution (INTS id, COEF *x)
{
#ifdef DISTRIBUTED
  INT    iter;
  int    dof;
#ifdef CENTRALISED
  INT    nodenbr;
  INT   *intern_nodelist;
  int    iterdof;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  CHECK_FACT(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (NULL == x)
    {
      return MURGE_ERR_PARAMETER;
    }
#ifdef CENTRALISED
  nodenbr = pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));
  MURGE_MEMALLOC(intern_nodelist, nodenbr, INT);
  if (NO_ERR != ( pastix_getLocalNodeLst(&(solvers[id]->pastix_data),
                                         intern_nodelist)))
    return MURGE_ERR_SOLVER;

  for (iter = 0; iter < nodenbr; iter ++)
    {
      for (iterdof = 0; iterdof < dof; iterdof++)
        {
          x[iter*dof+iterdof] = solvers[id]->b[(intern_nodelist[iter]-1)*
                                               dof+iterdof];
        }
    }
  memFree_null(intern_nodelist);
#else
  for (iter = 0; iter < solvers[id]->n*dof; iter ++)
    {
      x[iter] = solvers[id]->b[iter];
    }
#endif
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
/* TODO : optimiser */
INTS MURGE_GetSolution      (INTS id, INTS n, INTS *coefsidx, COEF *x,
                             INTS mode)
{
#ifdef DISTRIBUTED
  INT    iter;
  COEF  *tmpx = NULL;
  int    dof;
  int    iterdof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  CHECK_FACT(id);

  MURGE_MEMALLOC(tmpx, solvers[id]->N, COEF);
  return MURGE_ERR_ALLOCATE;
  MURGE_GetGlobalSolution(id, tmpx, -1);
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  for (iter = 0; iter < n; iter ++)
    {
      for (iterdof = 0; iterdof < dof; iterdof++)
        {
          x[iter*dof+iterdof] = tmpx[coefsidx[iter]*dof+iterdof-
                                     solvers[id]->pastix_data->iparm[IPARM_BASEVAL]];
        }
    }
  memFree_null(tmpx);


  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*******************************************************************************
 * Group: Cleaning up this mess
 */

INTS MURGE_Clean(INTS id){
#ifdef DISTRIBUTED
  INT    * iparm;
  double * dparm;
  iparm = solvers[id]->pastix_data->iparm;
  dparm = solvers[id]->pastix_data->dparm;

  if (NULL != solvers[id]->colptr)
    memFree_null(solvers[id]->colptr);
  if (NULL != solvers[id]->rows)
    memFree_null(solvers[id]->rows);
  if (NULL != solvers[id]->values)
    memFree_null(solvers[id]->values);
  if (NULL != solvers[id]->b)
    memFree_null(solvers[id]->b);
  if (NULL != solvers[id]->l2g)
    memFree_null(solvers[id]->l2g);
  if (NULL != solvers[id]->g2l)
    memFree_null(solvers[id]->g2l);
  if (NULL != solvers[id]->tmpv)
    memFree_null(solvers[id]->tmpv);
  if (NULL != solvers[id]->perm)
    memFree_null(solvers[id]->perm);
#ifdef CENTRALISED
  if (NULL != solvers[id]->invp)
    memFree_null(solvers[id]->invp);
#endif

  while (NULL != solvers[id]->sequences)
    {
      MURGE_AssemblyDeleteSequence(id, solvers[id]->sequences->ID);
    }
  pastix_task_clean(&(solvers[id]->pastix_data),
                    solvers[id]->pastix_data->pastix_comm);

#ifdef MURGE_THREADSAFE
  pthread_mutex_destroy(&solvers[id]->mutex_tmpmatrix);
#endif
  free(iparm); iparm = NULL;
  free(dparm); dparm = NULL;
  free(solvers[id]); solvers[id] = NULL;

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


INTS MURGE_Finalize(){
#ifdef DISTRIBUTED
  INT i;

  for (i=0; i< idnbr; i++)
    {
      if (solvers[i] != NULL)
        MURGE_Clean(i);
    }

  free(solvers); solvers = NULL;

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

INTS MURGE_GetInfoINT(INTS id,  INTS metric, INTL * value)
{
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  INT murge_param[1];

  murge_param[MURGE_IINFO_NNZ  - 1024] =  IPARM_NNZEROS;

  if (metric >= 1024)
    {
      metric = murge_param[metric-1024];
    }

  if (!( metric < IPARM_SIZE ))
    {
      errorPrint("metric is too big");
      return MURGE_ERR_PARAMETER;
    }

  if (metric < 0)
    {
      errorPrint("metric is negative");
      return MURGE_ERR_PARAMETER;
    }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  *value = solvers[id]->pastix_data->iparm[metric];

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

INTS MURGE_GetInfoREAL(INTS id,  INTS metric, REAL * value)
{
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  INT murge_param[1];

  murge_param[MURGE_RPARAM_EPSILON_ERROR  - 1024] =  DPARM_RELATIVE_ERROR;

  if (metric >= 1024)
    {
      metric = murge_param[metric-1024];
    }

  if (!( metric < IPARM_SIZE ))
    {
      errorPrint("metric is too big");
      return MURGE_ERR_PARAMETER;
    }

  if (metric < 0)
    {
      errorPrint("metric is negative");
      return MURGE_ERR_PARAMETER;
    }

  *value = solvers[id]->pastix_data->dparm[metric];
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
/*
 Function: MURGE_PrintError

 Print the error message corresponding to ierror
 Parameters:
 error_number  - Error identification number.

 Returns:
 MURGE_ERR_PARAMETER - If ierror does not match an error number
 MURGE_SUCCESS       - If function runned successfully.

 Fortran interface:
 >
 > SUBROUTINE MURGE_PRINTERROR(ERROR_NUMBER, IERROR)
 >   INTS, INTENT(IN)  :: IERROR
 >   INTS, INTENT(OUT) :: ERROR_NUMBER
 > END SUBROUTINE MURGE_PRINTERROR
 */
INTS MURGE_PrintError(INTS error_number)
{
#ifdef DISTRIBUTED
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/*
 Function: MURGE_ExitOnError

 Print the error message corresponding to ierror.
 If the ierr is not MURGE_SUCCESS then the program is stopped.

 Parameters:
 ierror         - Error identification number.

 Returns:
 MURGE_SUCCESS   - If function runned successfully,
 stop the program otherwise.

 Fortran interface:
 >
 > SUBROUTINE MURGE_EXITONERROR(ERROR_NUMBER, IERROR)
 >   INTS, INTENT(IN)  :: IERROR
 >   INTS, INTENT(OUT) :: ERROR_NUMBER
 > END SUBROUTINE MURGE_EXITONERROR
 */
INTS MURGE_ExitOnError(INTS error_number)
{
#ifdef DISTRIBUTED
  if  (error_number == MURGE_SUCCESS)
    return MURGE_SUCCESS;
  else
    exit(1);
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/*
 Group: Scaling
 */

/*
 Function: MURGE_GetGlobalNorm

 Compute the global norm array following a norm rule.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 norm    - Array of size global column number*dof which will contain
 the norm values
 root    - Indicates which processor will have the norm array
 at the end of the call, -1 for all.
 rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETGLOBALNORM(ID, NORM, ROOT, RULE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, RULE
 >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETGLOBALNORM
 */
INTS MURGE_GetGlobalNorm(INTS id, REAL *norm, INTS root, INTS rule)
{
#ifdef DISTRIBUTED
  INTS itercol;       /* Each column*/
  INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;    /* Column number of the value */
  INT  value_idx;     /* Index of the value */
  REAL*local_norm = NULL;
  INTS dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  print_debug(DBG_MURGE, ">> MURGE_GetGlobalNorm\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_NORM_RULE(rule);


  MURGE_MEMALLOC(local_norm, solvers[id]->N*dof, REAL);

  for(itercol = 0; itercol <  solvers[id]->N*dof; itercol++)
    {
      local_norm[itercol] = 0;
    }
  for(itercol = 0; itercol <  solvers[id]->n; itercol++)
    {
      for (iterrow = solvers[id]->colptr[itercol]-1;
           iterrow < solvers[id]->colptr[itercol+1]-1;
           iterrow++)
        {
          for (iterdof_col = 0;
               iterdof_col < dof;
               iterdof_col++)
            {
              for (iterdof_row = 0;
                   iterdof_row < dof;
                   iterdof_row++)
                {
                  value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
                  switch(rule)
                    {
                    case MURGE_NORM_MAX_COL:
                      scal_idx = iterdof_col + (solvers[id]->l2g[itercol]-1) * dof;
                      local_norm[scal_idx] =
                        MAX(local_norm[scal_idx],
                            ABS_FLOAT(solvers[id]->values[value_idx]));
                      break;

                    case MURGE_NORM_MAX_ROW:
                      scal_idx = iterdof_row + (solvers[id]->rows[iterrow]-1) * dof;
                      local_norm[scal_idx] =
                        MAX(local_norm[scal_idx],
                            ABS_FLOAT(solvers[id]->values[value_idx]));
                      break;

                    case MURGE_NORM_2_COL:
                      scal_idx = iterdof_col + (solvers[id]->l2g[itercol]-1) * dof;
                      local_norm[scal_idx] += solvers[id]->values[value_idx]*
                        solvers[id]->values[value_idx];
                      break;

                    case MURGE_NORM_2_ROW:
                      scal_idx = iterdof_row + (solvers[id]->rows[iterrow]-1) * dof;
                      local_norm[scal_idx] += solvers[id]->values[value_idx]*
                        solvers[id]->values[value_idx];
                      break;

                    default:
                      errorPrint("Rule not implemented");
                      return MURGE_ERR_NOT_IMPLEMENTED;
                    }
                }
            }
        }
    }


  if (rule == MURGE_NORM_2_COL ||
      rule == MURGE_NORM_2_ROW)
    {
      fprintf(stderr, "Reduce on norm\n");
      MPI_Allreduce(local_norm,norm,solvers[id]->N*dof,
                    MURGE_MPI_REAL,
                    MPI_SUM,
                    solvers[id]->pastix_data->pastix_comm);
      for(itercol = 0; itercol <  solvers[id]->N*dof; itercol++)
        {
          local_norm[itercol] = (REAL)sqrt((double)local_norm[itercol]);
        }
    }
  else
    {
      fprintf(stderr, "Reduce on norm 2\n");
      MPI_Allreduce(local_norm,norm,solvers[id]->N*dof,
                    MURGE_MPI_REAL,
                    MPI_MAX,
                    solvers[id]->pastix_data->pastix_comm);
    }
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 Function: MURGE_GetLocalNorm

 Compute the local norm array following a norm rule.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 norm    - Array of size local column number*dof which will contain
 the solution
 rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETLOCALNORM(ID, NORM, RULE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, RULE
 >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETLOCALNORM
 */
INTS MURGE_GetLocalNorm(INTS id, REAL *norm, INTS rule){
#ifdef DISTRIBUTED
  INTS itercol;       /* Each column*/
  INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS column_dof;    /* Column number of the value */
  INT  value_idx;     /* Index of the value */
  INTS dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_NORM_RULE(rule);
  CHECK_L2G(id);

  if (rule == MURGE_NORM_MAX_ROW || rule == MURGE_NORM_2_ROW)
    {
      errorPrint("PaStiX uses column distribution, local norm can't be a norm on rows");
      return MURGE_ERR_PARAMETER;
    }

  for(itercol = 0; itercol <  solvers[id]->n*dof; itercol++)
    norm[itercol] = 0;
  for(itercol = 0; itercol <  solvers[id]->n; itercol++)
    {
      for (iterrow = solvers[id]->colptr[itercol]-1;
           iterrow < solvers[id]->colptr[itercol+1]-1;
           iterrow++)
        {
          for (iterdof_col = 0;
               iterdof_col < dof;
               iterdof_col++)
            {
              for (iterdof_row = 0;
                   iterdof_row < dof;
                   iterdof_row++)
                {
                  column_dof = iterdof_col + itercol * dof;
                  value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;

                  if (rule == MURGE_NORM_2_COL)
                    {
                      norm[column_dof] += solvers[id]->values[value_idx];
                    }
                  else
                    {
                      norm[column_dof] =
                        MAX(norm[column_dof],
                            ABS_FLOAT(solvers[id]->values[value_idx]));
                    }
                }
            }
        }
    }

  for(itercol = 0; itercol <  solvers[id]->n*dof; itercol++)
    norm[itercol] = sqrt(norm[itercol]);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}


/*
 Function: MURGE_GetNorm

 Compute the indicated part of the norm array
 following a norm rule.

 Must be performed after assembly step.


 Parameters:
 id       - Solver instance identification number.
 n        - Number of coefficients user wants to get norm of.
 coefsidx - List of the coefficients user wants to get norm of.
 norm     - Array of size dof*n which will contain
 the solution.
 rule     - Rule to follow to build norm array, see <MURGE_NORM_RULES>
 mode     - Indicates if the user is sure to respect the distribution.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETNORM(ID, N, COEFSIDX, NORM, RULE, MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, MODE, N, RULE
 >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
 >   COEF, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETNORM
 */
INTS MURGE_GetNorm(INTS id,  INTS n, INTS *coefsidx, REAL *norm, INTS rule, INTS mode){
#ifdef DISTRIBUTED
  errorPrint("Not yet implemented");
  return MURGE_ERR_NOT_IMPLEMENTED;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}


/*
 Function: MURGE_ApplyGlobalScaling

 Apply scaling to local unknowns.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 scal    - Scaling user wants to apply.
 sc_mode - Indicate if the scaling is applied on rows or on columns.
 root    - Indicates which processor that posses the scaling array,
 -1 for all.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYGLOBALSCALING(ID, SCAL, SC_MODE, ROOT, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, SC_MODE
 >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYGLOBALSCALING

 */
INTS MURGE_ApplyGlobalScaling(INTS id, REAL *scal, INTS root, INTS sc_mode){
#ifdef DISTRIBUTED
  INTS itercol;       /* Each column*/
  INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;      /* Scaling array index */
  INT  value_idx;     /* Index of the value */
  INTS dof;
  REAL *scaling = NULL;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  CHECK_L2G(id);

  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (root == -1)
    {
      scaling = scal;
    }
  else
    {
      if (root != (solvers[id]->pastix_data)->procnum)
        {
          MURGE_MEMALLOC(scaling, solvers[id]->N*dof, REAL);
        }
      else
        {
          scaling = scal;
        }
      MPI_Bcast( scaling, solvers[id]->N*dof,
                 MURGE_MPI_REAL,
                 root,
                 solvers[id]->pastix_data->pastix_comm);
    }

  for(itercol = 0; itercol <  solvers[id]->n; itercol++)
    {
      for (iterrow = solvers[id]->colptr[itercol]-1;
           iterrow < solvers[id]->colptr[itercol+1]-1;
           iterrow++)
        {
          for (iterdof_col = 0;
               iterdof_col < dof;
               iterdof_col++)
            {
              for (iterdof_row = 0;
                   iterdof_row < dof;
                   iterdof_row++)
                {
                  if (sc_mode == MURGE_SCAL_COL)
                    {
                      scal_idx = iterdof_col + (solvers[id]->l2g[itercol]-1) * dof;
                    }
                  else
                    {
                      scal_idx = iterdof_row + (solvers[id]->rows[iterrow]-1) * dof;
                    }
                  value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
                  solvers[id]->values[value_idx] =
                    solvers[id]->values[value_idx] /scaling[scal_idx];
                }
            }
        }
    }
  if (root != -1 && root != (solvers[id]->pastix_data)->procnum)
    memFree_null(scaling);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 Function: MURGE_ApplyLocalScaling

 Apply the local scaling array on the matrix.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 scal    - Array of size local column number*dof which will contain
 the solution.
 sc_mode - Indicate if the scaling is applied on rows or on columns.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYLOCALSCALING(ID, SCAL, SC_MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, SC_MODE
 >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYLOCALSCALING
 */
INTS MURGE_ApplyLocalScaling(INTS id, REAL *scal, INTS sc_mode){
#ifdef DISTRIBUTED
  INTS itercol;       /* Each column*/
  INT  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;      /* Index in scaling array */
  INT  value_idx;     /* Index of the value */
  INTS dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (sc_mode == MURGE_SCAL_ROW)
    {
      /*
       * Building global to local column number array
       */
      CHECK_L2G(id);
    }
  for(itercol = 0; itercol <  solvers[id]->n; itercol++)
    {
      for (iterrow = solvers[id]->colptr[itercol]-1;
           iterrow < solvers[id]->colptr[itercol+1]-1;
           iterrow++)
        {
          for (iterdof_col = 0;
               iterdof_col < dof;
               iterdof_col++)
            {
              for (iterdof_row = 0;
                   iterdof_row < dof;
                   iterdof_row++)
                {

                  if (sc_mode == MURGE_SCAL_COL)
                    {
                      scal_idx = iterdof_col + itercol * dof;
                    }
                  else
                    {
                      scal_idx = iterdof_row + (solvers[id]->g2l[solvers[id]->rows[iterrow]-1]-1)*dof;
                    }
                  value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
                  solvers[id]->values[value_idx] =
                    solvers[id]->values[value_idx] /scal[scal_idx];
                }
            }
        }
    }
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 Function: MURGE_ApplyScaling

 Apply the scaling array on the indicated part of the matrix

 Must be performed after assembly step.


 Parameters:
 id       - Solver instance identification number.
 n        - Number of coefficients user wants to scale.
 coefsidx - List of the coefficients user wants to scale.
 scal     - Array of size dof*n which will contain
 the solution.
 sc_mode  - Indicate if the scaling is applied on rows or on columns.
 mode     - Indicates if the user is sure to respect the distribution.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYSCALING(ID, N, COEFSIDX, SCAL, SC_MODE, MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, SC_MODE, MODE, N
 >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
 >   COEF, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYSCALING
 */
INTS MURGE_ApplyScaling(INTS id,  INTS n, INTS *coefsidx, REAL *scal,
                        INTS sc_mode, INTS mode){
#ifdef DISTRIBUTED
  INTS itercol;       /* Each column*/
  INTS iterdof_col;   /* each dof on column */
  REAL *scaling = NULL;
  INTS dof;
  INTS baseval;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (mode == MURGE_ASSEMBLY_RESPECT)
    {
      /*
       * Building global to local column number array
       */
      CHECK_L2G(id);

      MURGE_MEMALLOC(scaling, solvers[id]->n*dof, REAL);
      for (itercol = 0; itercol < solvers[id]->n*dof; itercol++)
        scaling[itercol] = 1.0;
      for (itercol = 0; itercol < n; itercol++)
        for (iterdof_col =0; iterdof_col < dof; iterdof_col++)
          scaling[(solvers[id]->g2l[coefsidx[itercol]-baseval]-1)*dof+iterdof_col] = scal[itercol*dof+iterdof_col];
      MURGE_ApplyLocalScaling(id, scaling, sc_mode);
      memFree_null(scaling);
    }
  else
    {
      REAL * scaling_recv = NULL;
      MURGE_MEMALLOC(scaling, solvers[id]->N*dof, REAL);
      for (itercol = 0; itercol < solvers[id]->N*dof; itercol++)
        scaling[itercol] = 0.0;
      for (itercol = 0; itercol < n; itercol++)
        for (iterdof_col =0; iterdof_col < dof; iterdof_col++)
          scaling[(coefsidx[itercol]-baseval)*dof+iterdof_col] = scal[itercol*dof+iterdof_col];

      MURGE_MEMALLOC(scaling_recv, solvers[id]->N*dof, REAL);
      MPI_Allreduce(scaling, scaling_recv,
                    solvers[id]->N*dof,
                    MURGE_MPI_REAL,
                    MPI_SUM,
                    solvers[id]->pastix_data->pastix_comm);
      memFree_null(scaling);
      for (itercol = 0; itercol < solvers[id]->N*dof; itercol++)
        if (scaling_recv[itercol] == 0.0)
          scaling_recv[itercol] = 1.0;

      for (itercol = 0; itercol < n; itercol++)
        {
          for (iterdof_col =0; iterdof_col < dof; iterdof_col++)
            {
              if (scaling[(coefsidx[itercol]-baseval)*dof+iterdof_col] != scal[itercol*dof+iterdof_col])
                {
                  errorPrint("Multiple entries for the same scaling entry");
                  return MURGE_ERR_PARAMETER;
                }
            }
        }


      MURGE_ApplyGlobalScaling(id, scaling_recv, sc_mode, -1);
      memFree_null(scaling_recv);
    }
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}





/******************************************************************************
 * Group: Specific PaStiX functions.                                          *
 ******************************************************************************/

/******************************************************************************
 * Function: MURGE_Analyze                                                    *
 *                                                                            *
 * Perform matrix analyze:                                                    *
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

INTS MURGE_Analyze(INTS id){
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

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

INTS MURGE_Factorize(INTS id){
#ifdef DISTRIBUTED
  pastix_data_t   *pastix_data = solvers[id]->pastix_data;
  INT             *iparm       = pastix_data->iparm;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);

  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK)))
    {
      errorPrint("Need to set values before.");
      return MURGE_ERR_ORDER;
    }
  if (iparm[IPARM_ONLY_RAFF] ==  API_YES)
    {
      errorPrint("MURGE_Factorize is not compatible with IPARM_ONLY_RAFF == API_YES\n");
      return MURGE_ERR_PARAMETER;
    }
  pastix_data->iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
  iparm[IPARM_END_TASK]                = API_TASK_NUMFACT;

  dpastix(&pastix_data,
          pastix_data->pastix_comm,
          solvers[id]->n,
          solvers[id]->colptr,
          solvers[id]->rows,
          solvers[id]->values,
          solvers[id]->l2g,
          solvers[id]->perm,
          NULL,
          NULL,
          solvers[id]->nrhs,
          pastix_data->iparm,
          pastix_data->dparm);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_FACTO_OK);

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

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
INTS MURGE_ForceNoFacto(INTS id)
{
#ifdef DISTRIBUTED
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_FACTO_OK);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeNbr                                     *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   n  - Number of local nodes.                                              *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetLocalNodeNbr (INTS id, INTS n) {
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  solvers[id]->n = n;
  MPI_Allreduce(&solvers[id]->n,
                &solvers[id]->N, 1, COMM_INT,
                MPI_SUM, solvers[id]->pastix_data->pastix_comm);
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}


/******************************************************************************
 * Function: MURGE_ProductSetGlobalNodeNbr                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   N  - Number of global nodes.                                             *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetGlobalNodeNbr (INTS id, INTS N) {
#ifdef DISTRIBUTED
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_NODELST_OK))
    {
      errorPrint("%s muste be called before MURGE_ProductSetLocalNodeList",
                 __FUNCTION__);
      return MURGE_ERR_ORDER;
    }
  solvers[id]->N = N;
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");

  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeList                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id  - Solver instance identification number.                             *
 *   l2g - Local to global node numbers.                                      *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetLocalNodeList (INTS id, INTS * l2g) {
#ifdef DISTRIBUTED
  INTS i;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, INT);
  for (i = 0; i < solvers[id]->n; i++) {
    solvers[id]->l2g[i] = l2g[i];
  }

  cscd_build_g2l(solvers[id]->n,
                 solvers[id]->l2g,
                 solvers[id]->pastix_data->pastix_comm,
                 &solvers[id]->N,
                 &solvers[id]->g2l);

  /* No need to work on the graph nor factorize
   when we only perform product */
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_ONLY_PROD);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_TRUE (solvers[id]->state, MURGE_NODELST_OK);
  MURGE_STATE_FALSE(solvers[id]->state, MURGE_VALUES_OK);


  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

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
 * Returns:                                                                   *
 *   MURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_GetLocalProduct (INTS id, COEF *x)
{
#ifdef DISTRIBUTED
  COEF * glob_prod;
  INTS ierr, iter;
  MURGE_MEMALLOC(glob_prod, solvers[id]->N, COEF);

  if (MURGE_SUCCESS != (ierr = MURGE_GetGlobalProduct(id, glob_prod, -1)))
    return ierr;

  for (iter  = 0; iter < solvers[id]->n; iter++)
    {
      x[iter] = glob_prod[solvers[id]->l2g[iter]];
    }
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}

#ifdef DISTRIBUTED
struct product_data_
{
  INTS thread_id;
  murge_data_t * solver;
  COEF * t_prod;
  INTS ret;
};
typedef struct product_data_ product_data_t;

void* product_thread(void * data)
{
  INTS first, last;
  INTS gfirst, glast;
  INTS nnz, nnz_per_thread;
  INTS itercol;
  INTS iterrows;
  INTS dof, baseval;
  INTS row;
  COEF * mat = NULL;
  product_data_t * pdata = (product_data_t *)data;

  baseval = pdata->solver->pastix_data->iparm[IPARM_BASEVAL];
  dof = pdata->solver->pastix_data->iparm[IPARM_DOF_NBR];
  pdata->ret = MURGE_SUCCESS;


  nnz = pdata->solver->colptr[pdata->solver->n]-baseval;
  nnz_per_thread = nnz/pdata->solver->pastix_data->iparm[IPARM_THREAD_NBR];
  gfirst = pdata->thread_id * nnz_per_thread;
  glast = (pdata->thread_id +1)* nnz_per_thread;

  if (pdata->thread_id == pdata->solver->pastix_data->iparm[IPARM_THREAD_NBR]-1)
    {
      glast = nnz;
    }
  for (itercol = 0; itercol < pdata->solver->n; itercol++)
    {
      if (pdata->solver->colptr[itercol+1]-baseval > gfirst ||
          pdata->solver->colptr[itercol]-baseval  < glast)
        {
          first = MAX(gfirst, pdata->solver->colptr[itercol]-baseval);
          last  = MIN(glast,  pdata->solver->colptr[itercol+1]-baseval);
          for (iterrows = first;
               iterrows < last;
               iterrows++)
            {
              row = pdata->solver->rows[iterrows]-baseval;
              mat = &(pdata->solver->values[iterrows*dof*dof]);
              SOPALIN_GEMV("N", dof, dof, 1.0, mat, dof,
                           &(pdata->solver->b[itercol*dof]), 1, 1.0,
                           &(pdata->t_prod[row*dof]), 1);
            }

        }
    }
  return 0;
}
#endif

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
INTS MURGE_GetGlobalProduct (INTS id, COEF *x, INTS root)
{
#ifdef DISTRIBUTED
  INTS iter, iter2;
  INTS dof, baseval;
  int ret;
  product_data_t * pdata;
  COEF * my_prod = NULL;
  pthread_t        *calltab = NULL;
  int threadnbr;
#ifdef MURGE_TIME
  Clock            clock;
#endif

  CLOCK_INIT;
  if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK)))
    {
      errorPrint("Need to set values before.");
      return MURGE_ERR_ORDER;
    }

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  baseval = solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof = solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES) {
    errorPrint("Product only available with unsymmetric matrices.");
    return MURGE_ERR_NOT_IMPLEMENTED;
  }
  threadnbr = solvers[id]->pastix_data->iparm[IPARM_THREAD_NBR];
  MURGE_MEMALLOC(pdata, threadnbr, product_data_t);
  MURGE_MEMALLOC(calltab, threadnbr, pthread_t);
#ifdef MURGE_PRODUCT_CHECK_ZEROS
  {
    INTS itercol, iterrows, row;
    COEF max, sum, norm, critere;
    INT cnt, cnt2, cnt_sum, cnt2_sum, nz_glob;
    COEF *values = solvers[id]->values;
    critere = solvers[id]->pastix_data->dparm[DPARM_EPSILON_MAGN_CTRL];
    if (critere < 0.0)
      {
        critere = -critere;
      }
    else
      {
        max = 0;
        for (itercol = 0; itercol < solvers[id]->n; itercol++)
          {
            for (iter = 0; iter < dof; iter++) {
              sum = 0;
              for (iterrows = solvers[id]->colptr[itercol]-baseval;
                   iterrows < solvers[id]->colptr[itercol+1]-baseval;
                   iterrows++)
                {
                  row = solvers[id]->rows[iterrows]-baseval;
                  for (iter2 = 0; iter2 < dof; iter2++)
                    {
                      INT idx = iterrows*dof*dof+iter*dof+iter2;

                      sum = sum + ABS_FLOAT(values[idx]);
                    }

                }
              max = MAX(max,sum);
            }
          }

        MPI_Allreduce(&max, &norm, 1, MURGE_MPI_COEF, MPI_MAX,
                      solvers[id]->pastix_data->pastix_comm);

        critere = norm*sqrt(critere);
      }
    cnt = 0;
    cnt2 = 0;
    for (itercol = 0; itercol < solvers[id]->n; itercol++)
      {
        for (iter = 0; iter < dof; iter++) {
          sum = 0;
          for (iterrows = solvers[id]->colptr[itercol]-baseval;
               iterrows < solvers[id]->colptr[itercol+1]-baseval;
               iterrows++)
            {
              row = solvers[id]->rows[iterrows]-baseval;
              for (iter2 = 0; iter2 < dof; iter2++)
                {
                  INT idx = iterrows*dof*dof+iter*dof+iter2;
                  if (ABS_FLOAT(values[idx]) < critere)
                    cnt = cnt + 1;
                  if (values[idx] ==  0.0)
                    cnt2 = cnt2 + 1;
                }
            }
          max = MAX(max,sum);
        }
      }
    cnt_sum = 0;
    MPI_Reduce(&cnt, &cnt_sum, 1, COMM_INT, MPI_SUM, 0,
               solvers[id]->pastix_data->pastix_comm);
    MPI_Reduce(&cnt2, &cnt2_sum, 1, COMM_INT, MPI_SUM, 0,
               solvers[id]->pastix_data->pastix_comm);
    cnt = solvers[id]->colptr[solvers[id]->n]-1;
    MPI_Reduce(&cnt, &nz_glob, 1, COMM_INT, MPI_SUM, 0,
               solvers[id]->pastix_data->pastix_comm);
    nz_glob = nz_glob *dof*dof;
    if ((solvers[id]->pastix_data)->procnum == 0)
      {
        fprintf(stdout, "%d zeros in matrix from %d (%.3lg %%) critere :"
                " %.20lg\n",
                cnt_sum, nz_glob,
                (100.0*(double)(cnt_sum)/(double)(nz_glob)), critere);
        fprintf(stdout, "%d real zeros in matrix from %d (%.3lg %%)\n",
                cnt2_sum, nz_glob,
                (100.0*(double)(cnt2_sum)/(double)(nz_glob)));
      }
  }
#endif

  for (iter = 0; iter < threadnbr; iter++)
    {
      pthread_attr_t attr;
      pthread_attr_init(&attr);

      pdata[iter].thread_id = iter;
      pdata[iter].solver = solvers[id];
      MURGE_MEMALLOC(pdata[iter].t_prod, solvers[id]->N*dof, COEF);
      memset(pdata[iter].t_prod, 0, solvers[id]->N*dof*sizeof(COEF));
      ret = pthread_create(&(calltab[iter]), &attr,
                           product_thread, (void *)&(pdata[iter]));
      if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }
  for (iter = 0; iter < threadnbr; iter++)
    {
      ret = pthread_join(calltab[iter],(void**)NULL);
      if (ret) {errorPrint("thread join."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    }
  memFree_null(calltab);
  if (pdata[0].ret != MURGE_SUCCESS) return pdata[0].ret;
  /* Threads reduction */
  for (iter = 1; iter< threadnbr; iter++)
    {
      if (pdata[iter].ret != MURGE_SUCCESS) return pdata[iter].ret;
      for (iter2 = 0; iter2 < solvers[id]->N*dof; iter2++)
        pdata[0].t_prod[iter2] += pdata[iter].t_prod[iter2];
      memFree_null(pdata[iter].t_prod);
    }
  my_prod = pdata[0].t_prod;
  memFree_null(pdata);

  if (root == -1)
    {
      MPI_Allreduce(my_prod, x, solvers[id]->N*dof, MURGE_MPI_COEF, MPI_SUM,
                    solvers[id]->pastix_data->pastix_comm);
    }
  else
    {
      MPI_Reduce(my_prod, x, solvers[id]->N*dof, MURGE_MPI_COEF, MPI_SUM,
                 root, solvers[id]->pastix_data->pastix_comm);
    }
  memFree_null(my_prod);
  CLOCK_STOP;
#ifdef MURGE_TIME
  fprintf(stdout, " > MURGE_GetGlobalProduct computed in %.3g seconds (id=%d)\n",
          (double)CLOCK_GET, id);
#endif

  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */

}

/*
 WARNING: NEEDS TO BE CHECKED !
 */
INTS MURGE_SetLocalNodeList   (INTS id, INTS nodenbr, INTS *nodelist)
{
#ifdef DISTRIBUTED
  INT i;
  solvers[id]->n = nodenbr;
  /* On detruit colptr et rows, la distribution a changé */
  memFree_null(solvers[id]->colptr);
  memFree_null(solvers[id]->rows);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_GRAPH_OK)
    MURGE_STATE_TRUE(solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(solvers[id]->state, MURGE_SYMB_OK);
  MURGE_MEMALLOC(solvers[id]->l2g, nodenbr, INT);
  for (i = 0; i < nodenbr; i++) {
    solvers[id]->l2g[i] = nodelist[i];
  }
  return MURGE_SUCCESS;
#else /* DISTRIBUTED */
  fprintf (stderr,
           "Murge interface needs to  compile PaStiX with -DDISTRIBUTED\n");
  return MURGE_ERR_NOT_IMPLEMENTED;
#endif /* DISTRIBUTED */
}
