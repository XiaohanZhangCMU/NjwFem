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
 * File: murge_defines.h
 *
 *  This file define defines, macros and external funtion definition
 *  used to build <Murge> interface.
 *
 * About: Authors
 *   Xavier Lacoste  - xavier.lacoste@inria.fr
 */
#ifndef MURGE_DEFINE_H
#define MURGE_DEFINE_H
#include "pastix_internal.h"

/******************************************************************************/
/***                           Section: External functions                  ***/
/******************************************************************************/

FLOAT add_two_floats(FLOAT a, FLOAT b);
FLOAT keep_last(FLOAT a, FLOAT b);
FLOAT get_max(FLOAT a, FLOAT b);
FLOAT get_min(FLOAT a, FLOAT b);
int cscd_addlocal_int(INT   n   , INT *  ia   , INT *  ja   ,
		      FLOAT *  a   , INT * l2g,
		      INT   addn, INT *  addia, INT *  addja,
		      FLOAT *  adda, INT * addl2g,
		      INT * newn, INT ** newia, INT ** newja, FLOAT ** newa,
		      FLOAT (*add_fct)(FLOAT , FLOAT),
		      int dof, int malloc_flag);


void pastix_task_init(pastix_data_t **pastix_data,
		      MPI_Comm        pastix_comm,
		      INT            *iparm,
		      double         *dparm);


void pastix_initParam(INT    *iparm,
		      double *dparm);


void pastix_task_clean(pastix_data_t **pastix_data,
		       MPI_Comm        pastix_comm);

void pastix_welcome_print(pastix_data_t *pastix_data,
			  INT           *colptr,
			  INT            ln);

int pastix_checkMatrix_int(MPI_Comm pastix_comm,
			   int      verb,
			   int      flagsym,
			   int      flagcor,
			   INT      n,
			   INT    **colptr,
			   INT    **row,
			   FLOAT  **avals,
			   INT    **loc2glob,
			   int      dof,
			   int      flagalloc);

/******************************************************************************/
/***                           Section: Defines                             ***/
/******************************************************************************/

/*
  Defines: MPI Tags

  Tags for MPI communication.

  TAG_SIZE - To send size of a buffer.
  TAG_ROW  - To send rows.
  TAG_COL  - To send columns.
  TAG_FL   - To send an interval (first last).
  TAG_IJV  - To send <ijv_t> array.
  TAG_L2G  - To send local to global column numbers.
  TAG_VAL  - To send values.
*/
#define TAG_SIZE  1
#define TAG_ROW   2
#define TAG_COL   3
#define TAG_FL    4
#define TAG_IJV   5
#define TAG_L2G   6
#define TAG_VAL   7
#define TAG_SIZE2 8
#define TAG_VAL2  9
#define TAG_IJV2  10

/*
  Defines: State masks

  Bit masks used to define state variable.

  MURGE_INIT_OK     - If initialisation step has been called.
  MURGE_GRAPH_OK    - If graph of non zeros has been built.
  MURGE_GRAPH_BUILD - If we are in a graph building session.
  MURGE_VALUES_OK   - If Values of the matrix have been set.
  MURGE_BLEND_OK    - If preprocessing has been performed.
  MURGE_MATR_BUILD  - If we are in a matrix building session.
  MURGE_FACTO_OK    - If Factorization has been computed.
  MURGE_NODENBR_OK  - If node number has been given to user.
  MURGE_NODELST_OK  - If node list has been given to user.
  MURGE_RHS_OK      - If Right hand side has been set by user.
  MURGE_SYMB_OK     - If Symbolic factorization has been performed.
  MURGE_ONLY_PROD   - If we only compute producte
*/
#define MURGE_INIT_OK     1
#define MURGE_GRAPH_OK    2
#define MURGE_GRAPH_BUILD 4
#define MURGE_VALUES_OK   8
#define MURGE_BLEND_OK    16
#define MURGE_MATR_BUILD  32
#define MURGE_FACTO_OK    64
#define MURGE_NODENBR_OK  128
#define MURGE_NODELST_OK  256
#define MURGE_RHS_OK      512
#define MURGE_SYMB_OK     1024
#define MURGE_ONLY_PROD   2048

/******************************************************************************/
/***                           Section: Macros                              ***/
/******************************************************************************/

/*
  Macro: MURGE_MEMALLOC

  Allocate a space of size *size* x sizeof(*type*)
  at the adress indicated by ptr.

  Parameters:
  ptr   - address where to allocate.
  size  - Number of elements to allocate.
  types - Type of the elements to allocate.

  Returns:
  MURGE_ERR_ALLOCATE - If allocation fails.
*/
#define MURGE_MEMALLOC(ptr, size, type)                                 \
  {                                                                     \
    if (NULL == ((ptr) = (type *) memAlloc((size) * sizeof(type))))     \
      {                                                                 \
        errorPrint("Memory allocation error");                          \
        return MURGE_ERR_ALLOCATE;                                      \
      }                                                                 \
    PRINT_ALLOC(ptr, ((size)*sizeof(type)), __FILE__, __LINE__);        \
  }
#define MURGE_REALLOC(ptr, size, type)                                  \
  {                                                                     \
    if (NULL == ((ptr) = (type *) memRealloc(ptr,                       \
                                             (size) * sizeof(type))))   \
      {                                                                 \
        errorPrint("Memory allocation error");                          \
        return MURGE_ERR_ALLOCATE;                                      \
      }                                                                 \
  }

/*
  Macro: MURGE_STATE_FALSE

  Set *states* bit corresponding to *mask* to 0.

  Parameters:
  state - state variable
  mask  - information we want to set.
*/
#define MURGE_STATE_FALSE(state, mask)          \
  {                                             \
    state |= ( mask );                          \
    state ^= ( mask );                          \
  }


/*
  Macro: MURGE_STATE_TRUE

  Set *states* bit corresponding to *mask* to 1.

  Parameters:
  state - state variable
  mask  - information we want to set.
*/
#define MURGE_STATE_TRUE(state, mask)           \
  {                                             \
    state |= (mask);                            \
  }

/*
  Macro: MURGE_STATE_ISTRUE

  Check if *states* bit corresponding to *mask* is set to 1.

  Parameters:
  state - state variable
  mask  - information we want to test.

  Returns:
  true - if *mask* bit is set to 1
  else - otherwise
*/
#define MURGE_STATE_ISTRUE(state, mask) ( state == ( (state) | (mask) ) )

/*
  Macro: CHECK_SOLVER_ID

  Checks if solvers structure has been correctly set for instance *id*

  It checks that *id* value is in correct range and that *solvers*
  and solvers[id] have been allocated.

  Parameters:
  id - Solver instance ID we want to check.

  Returns:
  MURGE_ERR_PARAMETER - If *id* is not in correct range.
  MURGE_ERR_ORDER     - If *solvers* or *solvers[id]* are not allocated.
*/
#ifdef CENTRALISED
#  define INIT_CENTRALISED			\
  solvers[id]->invp        = NULL;
#else
#  define INIT_CENTRALISED {}
#endif
#define CHECK_SOLVER_ID(id)						\
  {									\
    if ( (idnbr > 0) && ((id < 0) || (id >= idnbr)))			\
      {									\
	errorPrint("Id is not in solvers array range");			\
	return MURGE_ERR_PARAMETER;					\
      }									\
    if (solvers == NULL)						\
      {									\
	errorPrint("You need to call MURGE_Initialize before");		\
	return MURGE_ERR_ORDER;						\
      }									\
    if (solvers[id] == NULL)						\
      {									\
	MURGE_MEMALLOC(solvers[id], 1, murge_data_t);			\
									\
	solvers[id]->n           = 0;					\
	solvers[id]->N           = 0;					\
	solvers[id]->colptr      = NULL;				\
	solvers[id]->rows        = NULL;				\
	solvers[id]->values      = NULL;				\
	solvers[id]->l2g         = NULL;				\
	solvers[id]->perm        = NULL;				\
	INIT_CENTRALISED;						\
	solvers[id]->b           = NULL;				\
	solvers[id]->nrhs        = 1;					\
	solvers[id]->state       = MURGE_INIT_OK;			\
	solvers[id]->pastix_data = NULL;				\
	solvers[id]->tmpv        = NULL;				\
	solvers[id]->tmpijv      = NULL;				\
	solvers[id]->tmpijv_node = NULL;				\
									\
	pastix_task_init(&(solvers[id]->pastix_data), MPI_COMM_WORLD,	\
			 NULL, NULL);					\
									\
      }									\
  }


/*
  Macro: CHECK_SOLVER_PARAM

  Checks if parameters have been set once for solver instance *id*.

  Checks if *iparm* or *dparm* are allocated.

  Parameters:
  id - Solver instance ID we want to check.

  Returns:
  MURGE_ERR_ORDER - If *iparm* or *dparm* are not allocated.
*/
#define CHECK_SOLVER_PARAM(id)						\
  {									\
    if (solvers[id]->pastix_data->iparm == NULL ||			\
	solvers[id]->pastix_data->dparm == NULL)			\
      {									\
	errorPrint("You need to call MURGE_SetDefaultOptions before");	\
	return MURGE_ERR_ORDER;						\
      }									\
  }

/*
  Macro: CHECK_PREPROCESSING

  Checks if preprocessing (blend) has been called.

  If it hasn't, it will allocate permutation tabular
  and call preprocessing step.

  After calling preprocessing, it will set local number of column
  and local to global column number tabular to their new values.

  Colptr and rows will be destroyed because it is obsolete,
  and state will be set to indicate that preprocessing has been performed.

  Parameters:
  id - Solver instance ID we want to check

  Returns:
  MURGE_ERR_ALLOCATE - If any allocation error occurs.
*/

#define CHECK_PREPROCESSING(id)			\
  {						\
    int err = check_preprocessing(id);		\
    if (err != MURGE_SUCCESS)			\
      return err;				\
  }

/*
  Macro: CHECK_FACT

  Checks if matrix values have been set.

  Checks if factorization has been performed.

  If not, it will call for it and set state.

  Parameters:
  id - Solver instance ID we want to check

  Returns:
  MURGE_ERR_ORDER - If values or right-hand-side member have not been set.

*/
#define CHECK_FACT(id)							\
  {									\
    pastix_data_t   *pastix_data = solvers[id]->pastix_data;		\
    INT             *iparm       = pastix_data->iparm;			\
									\
    if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_VALUES_OK)))	\
      {									\
	errorPrint("Need to set values before.");			\
	return MURGE_ERR_ORDER;						\
      }									\
    if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_RHS_OK)))	\
      {									\
	errorPrint("Need to set right-hand-side member before.");	\
	return MURGE_ERR_ORDER;						\
      }									\
    if (!(MURGE_STATE_ISTRUE(solvers[id]->state, MURGE_FACTO_OK)))	\
      {									\
	/* Si elle n'est pas faite, on effectue la factorisation */	\
	FILL_INTERNAL_CSC;						\
	solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =		\
	  API_TASK_NUMFACT;						\
      }									\
    else								\
      {									\
	/* Sinon la resolution suffit */				\
	solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =		\
	  API_TASK_SOLVE;						\
      }									\
    if (iparm[IPARM_ONLY_RAFF] ==  API_YES)				\
      {									\
	/* On set le second membre */					\
	iparm[IPARM_END_TASK] = API_TASK_SOLVE;				\
	dpastix(&pastix_data,						\
		pastix_data->pastix_comm,				\
		solvers[id]->n,						\
		solvers[id]->colptr,					\
		solvers[id]->rows,					\
		solvers[id]->values,					\
		solvers[id]->l2g,					\
		solvers[id]->perm,					\
		NULL,							\
		solvers[id]->b,						\
		solvers[id]->nrhs,					\
		pastix_data->iparm,					\
		pastix_data->dparm);					\
      }									\
    if (iparm[IPARM_MURGE_REFINEMENT] == API_YES) {			\
      iparm[IPARM_END_TASK] = API_TASK_REFINE;				\
    } else {								\
      iparm[IPARM_END_TASK] = API_TASK_SOLVE;				\
    }									\
    dpastix(&pastix_data,						\
	    pastix_data->pastix_comm,					\
	    solvers[id]->n,						\
	    solvers[id]->colptr,					\
	    solvers[id]->rows,						\
	    solvers[id]->values,					\
	    solvers[id]->l2g,						\
	    solvers[id]->perm,						\
	    NULL,							\
	    solvers[id]->b,						\
	    solvers[id]->nrhs,						\
	    pastix_data->iparm,						\
	    pastix_data->dparm);					\
    MURGE_STATE_TRUE(solvers[id]->state, MURGE_FACTO_OK);		\
  }

/*
  Macro: CHECK_L2G

  Checks if local to global array has been allocated.

  If not, it will correct number of local columns and set local to global
  column array.

  Parameters:
  id - Solver instance ID we want to check

  Returns:
  MURGE_ERR_ALLOCATE - If any allocation error occurs.
  MURGE_ERR_SOLVER   - If local to global array setting fails.

*/
#define CHECK_L2G(id)							\
  {									\
    if (solvers[id]->l2g == NULL)					\
      {									\
	solvers[id]->n =						\
	  pastix_getLocalNodeNbr(&(solvers[id]->pastix_data));		\
	if (solvers[id]->n != 0)					\
	  {								\
	    MURGE_MEMALLOC(solvers[id]->l2g, solvers[id]->n, INT);	\
	    if (EXIT_SUCCESS !=						\
		(pastix_getLocalNodeLst(&(solvers[id]->pastix_data),	\
					solvers[id]->l2g)))		\
	      return MURGE_ERR_SOLVER;					\
                                                                        \
            fprintf(stdout, "%s:%d Building G2L\n",__FILE__,__LINE__);  \
            cscd_build_g2l(solvers[id]->n,                              \
                           solvers[id]->l2g,                            \
                           solvers[id]->pastix_data->pastix_comm,       \
                           &solvers[id]->N,                             \
                           &solvers[id]->g2l);                          \
                                                                        \
	  }								\
      }									\
  }
/*
  Macro: CHECK_SCALING_MODE

  Check that the given mode exists.

  Parameters:
  mode - The given scaling mode.
  Returns:
  MURGE_ERR_PARAMETER - If mode is not defined.
*/
#define CHECK_SCALING_MODE(mode)			\
  if (mode != MURGE_SCAL_COL && mode != MURGE_SCAL_ROW)	\
    {							\
      errorPrint("Invalid scaling mode");		\
      return MURGE_ERR_PARAMETER;			\
    }
/*
  Macro: CHECK_NORM_RULE

  Check that the given norm rule exists.

  Parameters:
  rule - The given norm rule.
  Returns:
  MURGE_ERR_PARAMETER - If rule is not defined.
*/
#define CHECK_NORM_RULE(rule)			\
  if (rule != MURGE_NORM_MAX_COL &&		\
      rule != MURGE_NORM_2_COL   &&		\
      rule != MURGE_NORM_MAX_ROW &&		\
      rule != MURGE_NORM_2_ROW)			\
    {						\
      errorPrint("Invalid  norm rule");		\
      return MURGE_ERR_PARAMETER;		\
    }
/*
  Macro: CHOOSE_FUNC

  Will set *func* to the good function, depending to *op*.

  Parameters:
  op    - Operation flag.
  func  - Function pointer to set.

  Returns:
  MURGE_ERR_PARAMETER - If *op* doesn't exists.
*/
#ifdef TYPE_COMPLEX
#  define CHOOSE_FUNC(func, op)					\
  {								\
    /* Choix de l'operation a effectuer pour les doublons */	\
    switch(op)							\
      {								\
      case MURGE_ASSEMBLY_ADD:					\
	func = &add_two_floats;					\
	break;							\
      case MURGE_ASSEMBLY_OVW:					\
	func = &keep_last;					\
	break;							\
      default:							\
	func = NULL;						\
	return MURGE_ERR_PARAMETER;				\
	break;							\
      }								\
  }
#else
#  define CHOOSE_FUNC(func, op)					\
  {								\
    /* Choix de l'operation a effectuer pour les doublons */	\
    switch(op)							\
      {								\
      case MURGE_ASSEMBLY_ADD:					\
	func = &add_two_floats;					\
	break;							\
      case MURGE_ASSEMBLY_OVW:					\
	func = &keep_last;					\
	break;							\
      case MURGE_ASSEMBLY_MAX:					\
	func = &get_max;					\
	break;							\
      case MURGE_ASSEMBLY_MIN:					\
	func = &get_min;					\
	break;							\
      default:							\
	func = NULL;						\
	return MURGE_ERR_PARAMETER;				\
	break;							\
      }								\
  }
#endif

/*
  Macro: EXTEND_NODE_LIST

  Extend the number of node that the node list can receive if
  needed.

*/
#define EXTEND_NODE_LIST do {				\
    int ret;						\
    if (MURGE_SUCCESS !=				\
	( ret = extend_node_list(edgenbr_recv_node,	\
				 procnum,		\
				 &tmpijvsize_node,	\
				 &tmpijv_node,		\
				 &ijvptr,		\
				 dof,			\
				 &tmpvalues,		\
				 &valuesptr)))		\
      return ret;					\
  } while(0);

/*
  Macro: COPY_ELEMENT_TO_NODE

  Add a new node containing tmpijv[iter] at the
  end of the arrays tmpijv_node/tmpvalues.

*/
#define COPY_ELEMENT_TO_NODE do {				\
    EXTEND_NODE_LIST;						\
								\
    {								\
      INT nodenbr = edgenbr_recv_node[procnum];			\
      tmpijv_node[nodenbr].i =					\
	(tmpijv[iter].i - 1 - (tmpijv[iter].i-1)%dof)/dof + 1;	\
      tmpijv_node[nodenbr].j =					\
	(tmpijv[iter].j - 1 - (tmpijv[iter].j-1)%dof)/dof + 1;	\
      tmpijv_node[nodenbr].v =					\
	&(tmpvalues[nodenbr*dof*dof]);				\
								\
      /* line not found */					\
      memset(tmpijv_node[nodenbr].v,				\
	     0 , dof*dof*sizeof(FLOAT));			\
      tmpijv_node[nodenbr].v[(tmpijv[iter].i-1)%dof +		\
			     (tmpijv[iter].j-1)%dof*dof] =	\
	tmpijv[iter].v[0];					\
      edgenbr_recv_node[procnum]++;				\
    }								\
  } while(0)


/*
  Macro: FILL_INTERNAL_CSC

  Fill PaStiX internal CSC
*/
#define FILL_INTERNAL_CSC {						\
    INT err;								\
    INT * iparm = solvers[id]->pastix_data->iparm;			\
    if (NO_ERR !=							\
	((err =								\
	  pastix_checkMatrix_int(solvers[id]->pastix_data->pastix_comm, \
				 iparm[IPARM_VERBOSE],			\
				 iparm[IPARM_SYM],			\
				 API_YES,				\
				 solvers[id]->n,			\
				 &(solvers[id]->colptr),		\
				 &(solvers[id]->rows),			\
				 &(solvers[id]->values),		\
				 &(solvers[id]->l2g),			\
				 iparm[IPARM_DOF_NBR],			\
				 API_YES))))				\
      {									\
	errorPrint("pastix_checkMatrix : err %ld\n", (long)err);	\
	return MURGE_ERR_PARAMETER;					\
      }									\
									\
    pastix_fillin_csc(solvers[id]->pastix_data,				\
		      solvers[id]->pastix_data->pastix_comm,		\
		      solvers[id]->n,					\
		      solvers[id]->colptr,				\
		      solvers[id]->rows,				\
		      solvers[id]->values,				\
		      NULL,						\
                      0,                                                \
		      solvers[id]->l2g);				\
    solvers[id]->pastix_data->cscInternFilled = API_YES;		\
  }

/*
  Macro: UNLINK

  Suppress a file from the disk.

  Parameters:
  file - The file to suppress
  Returns:
  MURGE_ERR_IO - If an error occur.
*/
#define UNLINK(file)					\
  if (0 != unlink(file))				\
    {							\
      errorPrint("could not unlink %s\n", file);	\
      return MURGE_ERR_IO;				\
    }							\

/*
  Macro: LINK

  Create a symbolink link.

  Parameters:
  src  - file to link.
  dest - link path.
  Returns:
  MURGE_ERR_IO - If an error occur.
*/
#define LINK(src, dest)							\
  if (0 != symlink(src,dest))						\
    {									\
      if (errno == EEXIST)						\
	{								\
	  UNLINK(dest);							\
	  if (0 != symlink(src,dest))					\
	    {								\
	      errorPrint("Could not link %s to %s\n", dest, src);	\
	      return MURGE_ERR_IO;					\
	    }								\
	}								\
      else								\
	{								\
	  errorPrint("Could not link %s to %s\n", dest, src);		\
	  return MURGE_ERR_IO;						\
	}								\
    }
/*
  Macro: RENAME

  Move a file on disk.

  Parameters:
  src  - File to move.
  dest - New path.
  Returns:
  MURGE_ERR_IO - If an error occur.
*/
#define RENAME(src, dest)					\
  if (0 != rename(src, dest))					\
    {								\
      errorPrint("couldnt rename %s into %s", src, dest);	\
      return MURGE_ERR_IO;					\
    }


/*
  Macros: Time macros

  CLOCK_INIT - Start a clok
  CLOCK_STOP - Save clock time
  CLOCK_GET  - Get saved time (double value)
*/
#ifdef MURGE_TIME
#  define CLOCK_INIT {clockInit(&clock);clockStart(&clock);}
#  define CLOCK_STOP {clockStop(&clock);}
#  define CLOCK_GET  clockVal(&clock)
#else
#  define CLOCK_INIT
#  define CLOCK_STOP
#  define CLOCK_GET  0
#endif

#ifdef MURGE_DUMP_CSC
#  define MURGE_DUMP_GRAPH do {						\
    char name[256];							\
    sprintf(name,"Graph_%ld_%ld_",                                      \
            (long)solvers[id]->pastix_data->pastix_comm,                \
            (long)solvers[id]->ndump);                                  \
    solvers[id]->ndump++;                                               \
    cscd_save(solvers[id]->n,						\
	      solvers[id]->colptr,					\
	      solvers[id]->rows,					\
	      NULL,							\
	      NULL,							\
	      solvers[id]->l2g,						\
	      solvers[id]->pastix_data->iparm[IPARM_DOF_NBR],		\
	      name,							\
	      solvers[id]->pastix_data->pastix_comm);			\
  } while(0)

#  define MURGE_DUMP_MATRIX do {                                        \
    char name[256];							\
    sprintf(name,"Matrix_%ld_%ld_",                                     \
            (long)solvers[id]->pastix_data->pastix_comm,                \
            (long)solvers[id]->ndump);                                  \
    solvers[id]->ndump++;                                               \
    cscd_save(solvers[id]->n,						\
	      solvers[id]->colptr,					\
	      solvers[id]->rows,					\
	      solvers[id]->values,					\
	      solvers[id]->b,						\
	      solvers[id]->l2g,						\
	      solvers[id]->pastix_data->iparm[IPARM_DOF_NBR],		\
	      name,							\
	      solvers[id]->pastix_data->pastix_comm);			\
  } while(0)

#else  /* not MURGE_DUMP_CSC */
#  define MURGE_DUMP_GRAPH do {} while (0)
#  define MURGE_DUMP_MATRIX do {} while (0)
#endif /* not MURGE_DUMP_CSC */

#define MPI_INTL (sizeof(long) == sizeof(INTL))?MPI_LONG:MPI_INT
#define MPI_INTS (sizeof(long) == sizeof(INTS))?MPI_LONG:MPI_INT

#endif
