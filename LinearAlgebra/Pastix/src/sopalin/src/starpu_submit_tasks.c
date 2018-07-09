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
#ifdef WITH_STARPU
#  ifndef FORCE_NO_CUDA
#    include <cuda.h>
#  endif
#  ifdef STARPU_USE_DEPRECATED_API
#    undef STARPU_USE_DEPRECATED_API
#  endif
#  include <starpu.h>
#  include <starpu_profiling.h>
#  ifdef FORCE_NOMPI
#    include "nompi.h"
#  else
#    include <mpi.h>
#  endif
#  include "common_pastix.h"
#  include "out.h"
#  include "sopalin_define.h"
#  include "sopalin_acces.h"
#  include "symbol.h"
#  include "ftgt.h"
#  include "csc.h"
#  include "updown.h"
#  include "queue.h"
#  include "bulles.h"
#  include "solver.h"
#  include "sopalin_thread.h"
#  include "sopalin_time.h"
#  include "sopalin3d.h"
#  include "starpu_kernels.h"
#  include "../../perf/src/perf.h"
#  include "starpu_submit_tasks.h"
#  include "starpu_updo.h"
#  include "sopalin_init.h"
#define dump_all API_CALL(dump_all)
void  dump_all                 (SolverMatrix *, CscMatrix * cscmtx, int);

#  ifdef STARPU_USE_CUDA
#    if ((!defined PREC_DOUBLE) || (!(defined __CUDA_ARCH__) || __CUDA_ARCH__ >= 130))
#      if !(defined PREC_DOUBLE && defined TYPE_COMPLEX && CUDA_SM_VERSION < 20)
#        ifndef FORCE_NO_CUDA
#          define STARPU_USE_CUDA_GEMM_FUNC
#        endif
#      endif
#    endif
#  endif

#  ifdef TYPE_COMPLEX
#    ifdef PREC_DOUBLE
#      define PREFIX  "Z"
#    else
#      define PREFIX  "C"
#    endif
#  else
#    ifdef PREC_DOUBLE
#      define PREFIX  "D"
#    else
#      define PREFIX  "S"
#    endif
#  endif

#  define CUDA_CALL(x) do                                       \
    {                                                           \
      if (cudaSuccess != x)                                     \
        {                                                       \
          errorPrint("%s (%s,%d)\n",#x, __FILE__,__LINE__);     \
          assert(0);                                            \
        }                                                       \
    } while(0)


#  define USE_TASK_DEP

#  if (STARPU_MAJOR_VERSION < 1)
#    error "PaStiX requires STARPU >= 1"
#  endif /* STARPU_MAJOR_VERSION */

#  define SPARSE_GEMM
#  define ARCH_CPU  0
#  define ARCH_CUDA 1


static size_t trf_size(struct starpu_task *task,
                       enum starpu_perf_archtype arch,
                       unsigned nimpl)
{
  starpu_trf_data_t * args         = (starpu_trf_data_t*)task->cl_arg;
  Sopalin_Data_t    * sopalin_data = args->sopalin_data;
  SolverMatrix      * datacode     = sopalin_data->datacode;
  INT                 stride       = STARPU_MATRIX_GET_LD(task->handles[0]);
  size_t              dima         = CBLK_COLNBR(args->cblknum);
  return OPS_PPF(dima) + OPS_TRSM(dima, stride);

}

static size_t gemm_size(struct starpu_task *task,
                        enum starpu_perf_archtype arch,
                        unsigned nimpl)
{
  starpu_gemm_data_t         * args         = (starpu_gemm_data_t*)task->cl_arg;
  Sopalin_Data_t             * sopalin_data = args->sopalin_data;
  INT                          stride       = STARPU_MATRIX_GET_LD(task->handles[0]);
  INT                          bloknum      = args->bloknum;
  SolverMatrix               * datacode     = sopalin_data->datacode;
  INT                          indblok      = SOLV_COEFIND(bloknum);
  INT                          dimi         = stride - indblok;
  INT                          dimj         = BLOK_ROWNBR(bloknum);
  INT                          dima         = STARPU_MATRIX_GET_NY(task->handles[0]);
  size_t size = 0;

  size += OPS_GEMM(dimi,dimj,dima);
  return size;
}

static struct starpu_perfmodel GEMM_model = {
  .type = STARPU_HISTORY_BASED,
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = gemm_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = gemm_size }
  },
  .symbol = PREFIX "GEMM"
};

#  define starpu_partition_data                  API_CALL(starpu_partition_data)


#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
/* LU */
static struct starpu_perfmodel GETRF_TRSM_model = {
  .type = STARPU_HISTORY_BASED,
  .symbol = PREFIX "GETRF_TRSM",
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = trf_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = trf_size }
  }
};

#      define getrfsp1d_cl             PASTIX_PREFIX_F(getrfsp1d_cl)
#      define getrfsp1d_gemm_cl        PASTIX_PREFIX_F(getrfsp1d_gemm_cl)
#      define getrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(getrfsp1d_sparse_gemm_cl)

struct starpu_codelet getrfsp1d_cl =
{
  .where = STARPU_CPU
#ifndef FORCE_NO_CUDA
#  ifdef WITH_MAGMABLAS
  |STARPU_CUDA
#  endif /* WITH_MAGMABLAS */
#endif /* not FORCE_NO_CUDA */
  ,
  .cpu_func = getrfsp1d_starpu_cpu,
#ifndef FORCE_NO_CUDA
#  ifdef WITH_MAGMABLAS
  .cuda_func = getrfsp1d_starpu_cuda,
#  endif /* WITH_MAGMABLAS */
#endif /* not FORCE_NO_CUDA */
  .model = &GETRF_TRSM_model,
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_RW}
};

struct starpu_codelet getrfsp1d_gemm_cl =
{
  .where = STARPU_CPU
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  |STARPU_CUDA
#      endif
  ,
  .cpu_func = getrfsp1d_gemm_starpu_cpu,
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  .cuda_func = getrfsp1d_gemm_starpu_cuda,
#      endif
  .model = &GEMM_model,
  .nbuffers = 5,
  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH}
};

struct starpu_codelet getrfsp1d_sparse_gemm_cl =
{
  .where = STARPU_CPU
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  |STARPU_CUDA
#      endif
  ,
  .cpu_func = getrfsp1d_sparse_gemm_starpu_cpu,
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  .cuda_func = getrfsp1d_sparse_gemm_starpu_cuda,
#      endif
  .model = &GEMM_model,
#      if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 5,
#      else
  .nbuffers = 6,
#      endif
  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH,
#      if !(defined STARPU_BLOCKTAB_SELFCOPY)
            STARPU_R
#      endif
  }
};
#      define cl_trf       getrfsp1d_cl
#      ifdef SPARSE_GEMM
#        define cl_gemm      getrfsp1d_sparse_gemm_cl
#      else
#        define cl_gemm      getrfsp1d_gemm_cl
#      endif
#    else /* SOPALIN_LU */
/* LLT */
static struct starpu_perfmodel POTRF_TRSM_model = {
  .type = STARPU_HISTORY_BASED,
  .symbol = PREFIX "POTRF_TRSM",
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = trf_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = trf_size }
  }

};
#      define potrfsp1d_cl             PASTIX_PREFIX_F(potrfsp1d_cl)
#      define potrfsp1d_gemm_cl        PASTIX_PREFIX_F(potrfsp1d_gemm_cl)
#      define potrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(potrfsp1d_sparse_gemm_cl)
struct starpu_codelet potrfsp1d_cl =
{
  .where = STARPU_CPU
#ifndef FORCE_NO_CUDA
#  ifdef WITH_MAGMABLAS
  |STARPU_CPU
#  endif
#endif /* not FORCE_NO_CUDA */
  ,
  .cpu_func = potrfsp1d_starpu_cpu,
#ifndef FORCE_NO_CUDA
#  ifdef WITH_MAGMABLAS
  .cuda_func = getrfsp1d_starpu_cuda,
#  endif /* WITH_MAGMABLAS */
#endif /* not FORCE_NO_CUDA */
  .model = &POTRF_TRSM_model,
  .nbuffers = 1,
  .modes = {
    STARPU_RW}
};
struct starpu_codelet potrfsp1d_gemm_cl =
{
  .where = STARPU_CPU,
  .cpu_func = potrfsp1d_gemm_starpu_cpu,
  .model = &GEMM_model,
  .nbuffers = 3,
  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH}
};
struct starpu_codelet potrfsp1d_sparse_gemm_cl =
{
  .where = STARPU_CPU
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  | STARPU_CUDA
#      endif
  ,
  .cpu_func = potrfsp1d_sparse_gemm_starpu_cpu,
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  .cuda_func = potrfsp1d_sparse_gemm_starpu_cuda,
#      endif
  .model = &GEMM_model,
#      if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 3,
#      else
  .nbuffers = 4,
#      endif

  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH,
#      if !(defined STARPU_BLOCKTAB_SELFCOPY)
            STARPU_R
#      endif
  }
};
#      define cl_trf       potrfsp1d_cl
#      ifdef SPARSE_GEMM
#        define cl_gemm      potrfsp1d_sparse_gemm_cl
#      else
#        define cl_gemm      potrfsp1d_gemm_cl
#      endif
#    endif /* SOPALIN_LU */
#  else  /* CHOL_SOPALIN */
/* LDLT */
static struct starpu_perfmodel HETRF_TRSM_model = {
  .type = STARPU_HISTORY_BASED,
  .symbol = PREFIX "HETRF_TRSM",
  .per_arch =
  {
    [STARPU_CPU_DEFAULT][0]  = { .size_base = trf_size },
    [STARPU_CUDA_DEFAULT][0] = { .size_base = trf_size }
  }

};
#    ifdef HERMITIAN
#      define hetrfsp1d_cl             PASTIX_PREFIX_F(hetrfsp1d_cl)
#      define hetrfsp1d_gemm_cl        PASTIX_PREFIX_F(hetrfsp1d_gemm_cl)
#      define hetrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(hetrfsp1d_sparse_gemm_cl)
#    else
#      define hetrfsp1d_cl             PASTIX_PREFIX_F(sytrfsp1d_cl)
#      define hetrfsp1d_gemm_cl        PASTIX_PREFIX_F(sytrfsp1d_gemm_cl)
#      define hetrfsp1d_sparse_gemm_cl PASTIX_PREFIX_F(sytrfsp1d_sparse_gemm_cl)
#    endif
struct starpu_codelet hetrfsp1d_cl =
{
  .where = STARPU_CPU,
  .cpu_func = hetrfsp1d_starpu_cpu,
  .model = &HETRF_TRSM_model,
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_SCRATCH}
};

struct starpu_codelet hetrfsp1d_gemm_cl =
{
  .where = STARPU_CPU,
  .cpu_func = hetrfsp1d_gemm_starpu_cpu,
  .model = &GEMM_model,
  .nbuffers = 3,
  .modes = {STARPU_R,
            STARPU_RW,
            STARPU_SCRATCH}
};

struct starpu_codelet hetrfsp1d_sparse_gemm_cl =
{
  .where = STARPU_CPU,
  .cpu_func = hetrfsp1d_gemm_starpu_cpu,
  .model = &GEMM_model,
#    if (defined STARPU_BLOCKTAB_SELFCOPY)
  .nbuffers = 3,
#    else
  .nbuffers = 4,
#    endif
  .modes = { STARPU_R,
             STARPU_RW,
             STARPU_SCRATCH,
#    if !(defined STARPU_BLOCKTAB_SELFCOPY)
             STARPU_R
#    endif
  }
};
#    define cl_trf       hetrfsp1d_cl
#    ifdef SPARSE_GEMM
#      define cl_gemm      hetrfsp1d_sparse_gemm_cl
#    else
#      define cl_gemm      hetrfsp1d_gemm_cl
#    endif
#  endif /* CHOL_SOPALIN */



/*
 Function: starpu_data_partition

 Initialize column blocks handlers.

 Parameters:
 sopalin_data - PaStiX global data structure.
 L_handle     - Handles for L column blocks.
 U_handle     - Handles for U column blocks.
 */
static void starpu_partition_data(Sopalin_Data_t * sopalin_data,
                                  starpu_data_handle_t * L_handle,
                                  starpu_data_handle_t * U_handle)
{
  SolverMatrix       * datacode         = sopalin_data->datacode;
  INT itercblk;

  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++)
    {
      starpu_matrix_data_register(&(L_handle[itercblk]), 0,
                                  (uintptr_t)SOLV_COEFTAB(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  CBLK_COLNBR(itercblk),
                                  sizeof(FLOAT));
#  if defined USE_TASK_DEP
      starpu_data_set_sequential_consistency_flag(L_handle[itercblk], 0);
#  endif

#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
      starpu_matrix_data_register(&(U_handle[itercblk]), 0,
                                  (uintptr_t)SOLV_UCOEFTAB(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  (uint32_t)SOLV_STRIDE(itercblk),
                                  CBLK_COLNBR(itercblk),
                                  sizeof(FLOAT));
#      if (defined USE_TASK_DEP)
      starpu_data_set_sequential_consistency_flag(U_handle[itercblk], 0);
#      endif
#    endif
#  endif
    }
}

/*
 * Function: starpu_init_smp
 *
 * Initialize thread data structure for factorization when using StarPU.
 */
#  define starpu_init_smp API_CALL(starpu_init_smp)
void*
starpu_init_smp (void * arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  int init;
  init = INIT_COMPUTE;
  if (THREAD_FUNNELED_OFF)
    {
      init = init | INIT_SEND;
      if (THREAD_COMM_OFF)
        {
          init = init | INIT_RECV;
        }
    }
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
      sopalin_init_smp(sopalin_data, argument->me, API_YES, init);
    }
  else
    {
      sopalin_init_smp(sopalin_data, argument->me, API_NO, init);
    }
  return NULL;
}

/*
 Struct: sopthread_data

 Structure conotaining the thread number and a pointer to data given to thread in parameters.
 */
typedef struct starpu_loop_data {
  int                    me;
  starpu_data_handle_t * L_handle;
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  starpu_data_handle_t * U_handle;
#  endif
  starpu_data_handle_t   WORK_handle;
  starpu_data_handle_t   blocktab_handle;
  starpu_trf_data_t    * trf_args;
  starpu_gemm_data_t   * gemm_args;
  struct starpu_task  **starpu_tasktab;
  Sopalin_Data_t       * sopalin_data;                 /*+ Data given to thread as argument +*/
  int                  ctx;
  int                  ctx_nbr;
  int                  first, last;
#  ifdef STARPU_CONTEXT
  pthread_barrier_t    * barrier;
#  endif
} starpu_loop_data_t;

/*
 * Function: starpu_submit_loop
 *
 * Submit the tasks.
 */
#  define starpu_submit_loop API_CALL(starpu_submit_loop)
void*
starpu_submit_loop (void * arg)
{
  starpu_loop_data_t  *starpu_loop_data  = (starpu_loop_data_t*)(arg);
  Sopalin_Data_t      *sopalin_data      = (Sopalin_Data_t *)(starpu_loop_data->sopalin_data);
  SolverMatrix        *datacode          = sopalin_data->datacode;
  starpu_trf_data_t   *trf_args          = starpu_loop_data->trf_args;
  struct starpu_task **starpu_tasktab    = starpu_loop_data->starpu_tasktab;
  starpu_gemm_data_t  *gemm_args         = starpu_loop_data->gemm_args;
  int                  me                = starpu_loop_data->me;
#  ifdef STARPU_CONTEXT
  int                  ctx_nbr           = starpu_loop_data->ctx_nbr;
  int                  first             = starpu_loop_data->first;
  int                  last              = starpu_loop_data->last;
#  endif /* STARPU_CONTEXT */
  INT itertask;
  INT n_cblks = 0, n_tasks = 0;
#  if (defined USE_TASK_DEP)
  INT deps_tab_size = 1024;
  struct starpu_task ** deps;
  MALLOC_INTERN(deps, deps_tab_size, struct starpu_task *);
#  endif
#  ifdef STARPU_CONTEXT
  starpu_set_sched_ctx(&(starpu_loop_data->ctx));
#  endif

  /* For all column blocks we add a diag+trsm task.
   For all bloc in column block, we add a gemm task.
   */
#  ifdef STARPU_CONTEXT
  if (me == 0)
    for (itertask=0;itertask<SYMB_BLOKNBR;itertask++)
      {
        starpu_tasktab[itertask] = starpu_task_create();
      }
  pthread_barrier_wait(starpu_loop_data->barrier);
#  endif

  for (itertask=0;itertask<SOLV_TASKNBR;itertask++)
    {
      INT itercblk = TASK_CBLKNUM(itertask);
      INT ret;
      INT iterbloc;
      INT handle_idx;
      struct starpu_task *task_diag;
#  ifdef STARPU_CONTEXT
      INT threadid = TASK_THREADID(itertask);
      if ( ( me <  ctx_nbr-1 && (threadid > last || threadid < first)) ||
           ( me == ctx_nbr-1 && threadid < first ) )
        continue;
#  endif
      n_cblks++;
      n_tasks++;
#  ifdef STARPU_CONTEXT
      task_diag = starpu_tasktab[SYMB_BLOKNUM(itercblk)];
#  else
      task_diag = starpu_task_create();
#  endif
      /* We compute diagonal factorisation and TRSM */
      trf_args[itercblk].cblknum      = itercblk;
      trf_args[itercblk].sopalin_data = sopalin_data;

      task_diag->cl = &cl_trf;
      task_diag->cl_arg = &(trf_args[itercblk]);
      task_diag->destroy = 0;
      handle_idx = 0;

      task_diag->handles[handle_idx++] = starpu_loop_data->L_handle[itercblk];
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
      /* LU */
      task_diag->handles[handle_idx++] = starpu_loop_data->U_handle[itercblk];
#    endif
#  else /* CHOL_SOPALIN */
      /* LDLT */
      task_diag->handles[handle_idx++] = starpu_loop_data->WORK_handle;
#  endif

      ASSERT(handle_idx == cl_trf.nbuffers, MOD_SOPALIN);

#  if (defined USE_TASK_DEP)
      {
        INT gcblk2list = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB( itercblk ));
        INT browk      = ( gcblk2list != -1 )?
          UPDOWN_LISTPTR( gcblk2list    ):-1;
        INT browk1     = ( gcblk2list != -1 )?
          UPDOWN_LISTPTR( gcblk2list + 1):-1;
        INT ndeps, iter;
        ndeps = browk1-browk;

        starpu_tasktab[SYMB_BLOKNUM(itercblk)] = task_diag;
        if (ndeps > deps_tab_size)
          {
            memFree_null(deps);
            deps_tab_size = ndeps;
            MALLOC_INTERN(deps, ndeps, struct starpu_task *);
          }
        for (iter = 0; iter < ndeps; iter++)
          {
            deps[iter] = starpu_tasktab[UPDOWN_LISTBLOK(browk+iter)];
          }
        if (ndeps != 0)
          starpu_task_declare_deps_array(task_diag, ndeps, deps);
      }
#  endif /* USE_TASK_DEP */

      ret = starpu_task_submit(task_diag);
      STARPU_ASSERT(!ret);

      for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
           iterbloc < SYMB_BLOKNUM(itercblk+1);
           iterbloc ++)
        {
          struct starpu_task * task_gemm;
          INT blocnbr;
          n_tasks++;
#  ifdef   STARPU_CONTEXT
          task_gemm = starpu_tasktab[iterbloc];
#  else
          task_gemm = starpu_task_create();
#  endif
          blocnbr = SYMB_BLOKNUM(itercblk+1) - iterbloc;
          /* We compute GEMM */
          gemm_args[iterbloc].cblknum      = itercblk;
          gemm_args[iterbloc].bloknum      = iterbloc;
          gemm_args[iterbloc].nblocs       = blocnbr;
          gemm_args[iterbloc].fcblknum     = SYMB_CBLKNUM(iterbloc);
          gemm_args[iterbloc].sopalin_data = sopalin_data;
#  if (defined STARPU_BLOCKTAB_SELFCOPY)
          gemm_args[iterbloc].d_blocktab   = d_blocktab;
#  endif

          task_gemm->cl = &cl_gemm;
          task_gemm->cl_arg = &(gemm_args[iterbloc]);
          task_gemm->destroy = 0;

          handle_idx = 0;

          task_gemm->handles[handle_idx++] = starpu_loop_data->L_handle[itercblk];
          task_gemm->handles[handle_idx++] = starpu_loop_data->L_handle[SYMB_CBLKNUM(iterbloc)];
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
          task_gemm->handles[handle_idx++] = starpu_loop_data->U_handle[itercblk];
          task_gemm->handles[handle_idx++] = starpu_loop_data->U_handle[SYMB_CBLKNUM(iterbloc)];
#    endif
#  endif /* CHOL_SOPALIN */
          task_gemm->handles[handle_idx++] = starpu_loop_data->WORK_handle;
#  ifndef STARPU_BLOCKTAB_SELFCOPY
          task_gemm->handles[handle_idx++] = starpu_loop_data->blocktab_handle;
#  endif /* STARPU_BLOCKTAB_SELFCOPY */

          ASSERT(handle_idx == cl_gemm.nbuffers, MOD_SOPALIN);

#  ifdef USE_TASK_DEP
          starpu_tasktab[iterbloc] = task_gemm;
          deps[0] = task_diag;
          starpu_task_declare_deps_array(task_gemm, 1, deps);
#  endif /* USE_TASK_DEP */
          ret = starpu_task_submit(task_gemm);
          STARPU_ASSERT(!ret);
        }
    }

  /* Wait for all tasks to be finished */
#ifdef STARPU_CONTEXT
  starpu_task_wait_for_all();
#endif
  fprintf(stdout, "%ld cblks %ld tasks in context %d\n",
          (long int)n_cblks, (long int)n_tasks, (int)me);
#  if ( defined USE_TASK_DEP)
  memFree_null(deps);
#  endif
  return NULL;
}

/*
 * Funciton starpu_clean_smp
 *
 * Clean thread data structures when using starpu.
 */
#  define starpu_clean_smp API_CALL(starpu_clean_smp)
void*
starpu_clean_smp (void * arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  sopalin_clean_smp ( sopalin_data, argument->me );
  return NULL;
}

/*
 Function: starpu_submit_tasks

 Submit tasks to perform the decomposition of the matrix.

 Parameters:
 sopalin_data - PaStiX global data structure.

 Returns:
 NO_ERR
 */
int
starpu_submit_tasks(Sopalin_Data_t  * sopalin_data) {
  SolverMatrix         * datacode         = sopalin_data->datacode;
  Thread_Data_t        * thread_data;
  starpu_loop_data_t   * starpu_loop_data = NULL;
  starpu_trf_data_t    * trf_args;
  starpu_gemm_data_t   * gemm_args;
  starpu_data_handle_t * L_handle;
  starpu_data_handle_t * SM2X_handles = NULL;
#  ifdef USE_TASK_DEP
  INT                    task_number;
  struct starpu_task  ** starpu_tasktab;
#  endif
  starpu_data_handle_t * U_handle = NULL;
  starpu_data_handle_t   WORK_handle;
  INT itertask;
  int * blocktab;
#  ifdef STARPU_BLOCKTAB_SELFCOPY
  int ** d_blocktab;
#  else
  starpu_data_handle_t   blocktab_handle;
#  endif

  struct starpu_conf     conf;
  int                    ret;
  int                    cuda_nbr = sopalin_data->sopar->iparm[IPARM_CUDA_NBR];

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR+cuda_nbr, starpu_init_smp, sopalin_data,
                        0, NULL, NULL,
                        0, NULL, NULL);

  starpu_conf_init(&conf);
  if (NULL != conf.sched_policy_name)
    {
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        fprintf(stdout, OUT_STARPU_TP, conf.sched_policy_name);
    }
  else
    {
      conf.sched_policy_name = "heft";
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        fprintf(stdout, OUT_STARPU_STP, conf.sched_policy_name);
    }
  conf.ncpus = SOLV_THRDNBR;
  conf.ncuda = cuda_nbr;
  conf.nopencl = 0;
  conf.nspus   = 0;
  if (0 != (ret = starpu_init(&conf)))
    {
      errorPrint("Error %d initializing StarPU\n", ret);
    }
#  ifdef STARPU_PROFILING
  if ((ret = starpu_profiling_status_set(STARPU_PROFILING_ENABLE) < 0))
    {
      errorPrint("Error %d in starpu_profiling_enable\n", ret);
    }
#  endif
  {
    /* build blocktab */
    int iterblock;
    MALLOC_INTERN(blocktab, SYMB_BLOKNBR*2, int);
    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
      fprintf(stdout, "sizeof blocktab : %d integers\n",
              (int)(2*SYMB_BLOKNBR));

    for (iterblock = 0; iterblock < SYMB_BLOKNBR; iterblock++)
      {
        blocktab[2*iterblock]   = SYMB_FROWNUM(iterblock);
        blocktab[2*iterblock+1] = SYMB_LROWNUM(iterblock);
      }
#  ifdef STARPU_BLOCKTAB_SELFCOPY
    {
      int ndevices;
      int device_id;
      CUDA_CALL(cudaGetDeviceCount(&ndevices));
      MALLOC_INTERN(d_blocktab, ndevices, int*);
      for (device_id = 0; device_id < ndevices; device_id++)
        {
          cudaSetDevice(device_id);
          CUDA_CALL(cudaMalloc((void*)&(d_blocktab[device_id]),
                               2*SYMB_BLOKNBR*sizeof(int)));
          CUDA_CALL(cudaMemcpy((void*)d_blocktab[device_id], blocktab,
                               2*SYMB_BLOKNBR*sizeof(int),
                               cudaMemcpyHostToDevice));
        }
    }
#  else
    starpu_vector_data_register(&blocktab_handle, 0,
                                (uintptr_t)blocktab, 2*SYMB_BLOKNBR,
                                sizeof(int));
    starpu_data_set_sequential_consistency_flag(blocktab_handle, 0);
#  endif
  }

#  ifdef PASTIX_DUMP_FACTO
  dump_all(datacode, sopalin_data->sopar->cscmtx,
           ((datacode->updovct.sm2xtab!=NULL)?
            (DUMP_CSC | DUMP_SOLV | DUMP_SMB):(DUMP_CSC | DUMP_SOLV)));
#  endif

  thread_data = sopalin_data->thread_data[0];
  sopalin_data->sopar->diagchange = 0;
  SOPALIN_CLOCK_INIT;
#  ifdef STARPU_USE_CUDA
  /* starpu_helper_cublas_init(); */
#  endif


#  ifdef CHOL_SOPALIN
  starpu_vector_data_register(&WORK_handle, -1, (uintptr_t)NULL, SOLV_COEFMAX,
                              sizeof(FLOAT));
#  else
  starpu_vector_data_register(&WORK_handle, -1, (uintptr_t)NULL, 2*SOLV_COEFMAX,
                              sizeof(FLOAT));
#  endif

  MALLOC_INTERN(trf_args,    SYMB_CBLKNBR, starpu_trf_data_t);

  MALLOC_INTERN(L_handle,    SYMB_CBLKNBR, starpu_data_handle_t);

#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  MALLOC_INTERN(U_handle,    SYMB_CBLKNBR, starpu_data_handle_t);
#    endif
#  endif

  MALLOC_INTERN(gemm_args,      SYMB_BLOKNBR, starpu_gemm_data_t);


  {
    int itercblk;
    int max_cblksize = 0;
    int max_cblkcolnbr = 0;
    for (itercblk = 0; itercblk < SYMB_CBLKNBR; itercblk++)
      {
        max_cblksize   = MAX(max_cblksize,
                             CBLK_COLNBR(itercblk)*SOLV_STRIDE(itercblk));
        max_cblkcolnbr = MAX(max_cblkcolnbr,
                             CBLK_COLNBR(itercblk));
      }

    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
      fprintf(stdout, "Maximum cblk size %d, maximu cblk colnbr %d\n",
              max_cblksize, max_cblkcolnbr);
  }

#  ifdef USE_TASK_DEP
  task_number = 0;
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    task_number = SYMB_BLOKNBR;
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    task_number += 2*SYMB_BLOKNBR+SYMB_CBLKNBR;
  MALLOC_INTERN(starpu_tasktab, task_number, struct starpu_task *);
#  endif
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  starpu_partition_data(sopalin_data, L_handle, U_handle);
#  else
  starpu_partition_data(sopalin_data, L_handle, NULL);
#  endif
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    {
      if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_SOLVE)
        if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
          errorPrintW("Raffinement not available with StarPU,"
                      " only performing solve\n");
      MALLOC_INTERN(SM2X_handles, SYMB_CBLKNBR, starpu_data_handle_t);
      starpu_register_sm2x(sopalin_data, SM2X_handles);
    }
  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- Time after data registration %lf s\n",
            SOPALIN_CLOCK_GET);

#  ifdef STARPU_CONTEXT
  {
    INT iter;
    pthread_barrier_t barrier;
    pthread_t        *calltab = NULL;
    int ** contexts;
    int * contexts_size;
    int ret;
    int contexts_nbr;
    int cpu_per_ctx;
    int gpu_per_ctx;
    int cuda_idx = 0;
    int cpu_idx  = cuda_nbr;

    if (cuda_nbr == 0)
      {
        contexts_nbr = SOLV_THRDNBR+1;
        cpu_per_ctx  = SOLV_THRDNBR/(contexts_nbr-1);
        gpu_per_ctx  = 0;
      }
    else
      {
        contexts_nbr = cuda_nbr + 1;
        cpu_per_ctx  = SOLV_THRDNBR/(contexts_nbr-1);
        gpu_per_ctx  = cuda_nbr/(contexts_nbr-1);
      }

    MALLOC_INTERN(contexts, contexts_nbr, int*);
    MALLOC_INTERN(contexts_size, contexts_nbr, int);

    for (iter = 0; iter < contexts_nbr-1; iter++)
      {
        int size;
        int iter2 = 0;
        int this_cpu_per_ctx = cpu_per_ctx;
        int this_gpu_per_ctx = gpu_per_ctx;
        if (iter < SOLV_THRDNBR%(contexts_nbr-1))
          this_cpu_per_ctx++;
        if (iter <  cuda_nbr%(contexts_nbr-1))
          this_gpu_per_ctx++;
        /* if (iter == contexts_nbr - 2) */
        /*  { */
        /*    this_cpu_per_ctx += SOLV_THRDNBR%(contexts_nbr-1); */
        /*    this_gpu_per_ctx += cuda_nbr%(contexts_nbr-1); */
        /*  } */
        size = this_gpu_per_ctx + this_cpu_per_ctx;
        MALLOC_INTERN(contexts[iter], size, int);
        contexts_size[iter] = 0;
        for (iter2 = 0; iter2 < this_gpu_per_ctx; iter2++, cuda_idx++)
          {
            contexts[iter][contexts_size[iter]++] = cuda_idx;
          }
        for (iter2 = 0; iter2 < this_cpu_per_ctx; iter2++, cpu_idx++)
          {
            contexts[iter][contexts_size[iter]++] = cpu_idx;
          }
      }

    /* Global context */
    contexts_size[contexts_nbr-1] = SOLV_THRDNBR+cuda_nbr;
    MALLOC_INTERN(contexts[contexts_nbr-1],      SOLV_THRDNBR+cuda_nbr, int);
    for (iter = 0; iter < SOLV_THRDNBR+cuda_nbr; iter++)
      contexts[contexts_nbr-1][iter] = iter;

    MALLOC_INTERN(calltab, contexts_nbr, pthread_t);
    MALLOC_INTERN(starpu_loop_data, contexts_nbr, starpu_loop_data_t);
    pthread_barrier_init(&barrier, NULL, contexts_nbr);
    for (iter = 0; iter < contexts_nbr; iter++)
      {
        char string[256];
        pthread_attr_t attr;
        int iter2;
        pthread_attr_init(&attr);
        starpu_loop_data[iter].me               = iter;
        starpu_loop_data[iter].trf_args         = trf_args;
        starpu_loop_data[iter].gemm_args        = gemm_args;
#    ifdef USE_TASK_DEP
        starpu_loop_data[iter].starpu_tasktab   = starpu_tasktab;
#    endif
        starpu_loop_data[iter].L_handle         = L_handle;
#    if (defined CHOL_SOPALIN && defined SOPALIN_LU)
        starpu_loop_data[iter].U_handle         = U_handle;
#    endif
        starpu_loop_data[iter].WORK_handle      = WORK_handle;
        starpu_loop_data[iter].blocktab_handle  = blocktab_handle;
        starpu_loop_data[iter].sopalin_data     = sopalin_data;
        starpu_loop_data[iter].barrier          = &barrier;
        starpu_loop_data[iter].ctx_nbr          = contexts_nbr;
        if (iter < contexts_nbr - 1)
          {
            if (cuda_nbr > 0)
              starpu_loop_data[iter].first            = contexts[iter][1] - cuda_nbr;
            else
              starpu_loop_data[iter].first            = contexts[iter][0] - cuda_nbr;

            starpu_loop_data[iter].last             = contexts[iter][contexts_size[iter]-1]- cuda_nbr;
          }
        else
          {
            starpu_loop_data[iter].first = SOLV_THRDNBR;
            starpu_loop_data[iter].last  = -1;
          }
        sprintf(string, "ctx_%d", iter);
        fprintf(stdout, "%s [", string);
        for (iter2 = 0; iter2 < contexts_size[iter]-1; iter2++)
          fprintf(stdout, "%d, ", contexts[iter][iter2]);
        fprintf(stdout, "%d]\n", contexts[iter][iter2]);

        starpu_loop_data[iter].ctx              = starpu_create_sched_ctx("heft", contexts[iter], contexts_size[iter], &(string[0]));

        ret = pthread_create(&calltab[iter],&attr,starpu_submit_loop,(void *)&(starpu_loop_data[iter]));
        if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
      }
    for (iter=0;iter<contexts_nbr;iter++)
      {
        /* On ne recupere pas le thread qd il a pas été lancé */
        ret = pthread_join(calltab[iter],(void**)NULL);
        memFree_null(contexts[iter]);
        if (ret) {errorPrint("thread join."); EXIT(MOD_SOPALIN,THREAD_ERR);}
      }
    memFree_null(contexts);
    memFree_null(calltab);
    pthread_barrier_destroy(&barrier);
    memFree_null(starpu_loop_data);
  }

  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    starpu_submit_updown(sopalin_data, L_handle, U_handle, SM2X_handles, starpu_tasktab);
  starpu_task_wait_for_all();

#  else /* not STARPU_CONTEXT */
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
      MALLOC_INTERN(starpu_loop_data, 1, starpu_loop_data_t);
      starpu_loop_data[0].me               = 0;
      starpu_loop_data[0].trf_args         = trf_args;
      starpu_loop_data[0].gemm_args        = gemm_args;
      starpu_loop_data[0].starpu_tasktab   = starpu_tasktab;
      starpu_loop_data[0].L_handle         = L_handle;
#    if (defined CHOL_SOPALIN && defined SOPALIN_LU)
      starpu_loop_data[0].U_handle         = U_handle;
#    endif
      starpu_loop_data[0].WORK_handle      = WORK_handle;
      starpu_loop_data[0].blocktab_handle  = blocktab_handle;
      starpu_loop_data[0].sopalin_data     = sopalin_data;
      starpu_loop_data[0].ctx_nbr          = 1;
      starpu_loop_data[0].first            = 0;
      starpu_loop_data[0].last             = SOLV_THRDNBR;
      starpu_submit_loop (starpu_loop_data);
    }
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    starpu_submit_updown(sopalin_data, L_handle, U_handle, SM2X_handles, starpu_tasktab);
  starpu_task_wait_for_all();
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    memFree_null(starpu_loop_data);
#  endif

  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- submission and wait for all %lf s\n",
            SOPALIN_CLOCK_GET);


  /* Unregister buffers and leave starpu */
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
#  if(defined STARPU_PROFILING && defined USE_TASK_DEP)
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        {
          int worker;
          /* Display the occupancy of all workers during the test */

          for (worker = 0; worker < starpu_worker_get_count(); worker++)
            {
              struct starpu_worker_profiling_info worker_info;
              int ret = starpu_worker_get_profiling_info(worker, &worker_info);
              STARPU_ASSERT(!ret);

              double total_time     = starpu_timing_timespec_to_us(&worker_info.total_time);
              double executing_time = starpu_timing_timespec_to_us(&worker_info.executing_time);
              double sleeping_time  = starpu_timing_timespec_to_us(&worker_info.sleeping_time);

              float executing_ratio = 100.0*executing_time/total_time;
              float sleeping_ratio  = 100.0*sleeping_time/total_time;

              char workername[128];

              double delay_sum[2]  = {0.0, 0.0};
              double length_sum[2] = {0.0, 0.0};
              unsigned int cnt[2]  = {0,   0};

              for (itertask=0;itertask<SOLV_TASKNBR;itertask++)
                {
                  INT itercblk = TASK_CBLKNUM(itertask);
                  INT iterbloc;
                  struct starpu_task_profiling_info *info;
                  info = starpu_tasktab[SYMB_BLOKNUM(itercblk)]->profiling_info;
                  if (info->workerid == worker)
                    {
                      /* How much time did it take before the task started ? */
                      delay_sum[0] += starpu_timing_timespec_delay_us(&info->submit_time,
                                                                      &info->start_time);

                      /* How long was the task execution ? */
                      length_sum[0] += starpu_timing_timespec_delay_us(&info->start_time,
                                                                       &info->end_time);
                      cnt[0]++;
                    }

                  for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
                       iterbloc < SYMB_BLOKNUM(itercblk+1);
                       iterbloc ++)
                    {
                      info = starpu_tasktab[iterbloc]->profiling_info;
                      if (info->workerid == worker)
                        {
                          /* How much time did it take before the task started ? */
                          delay_sum[1] += starpu_timing_timespec_delay_us(&info->submit_time,
                                                                          &info->start_time);

                          /* How long was the task execution ? */
                          length_sum[1] += starpu_timing_timespec_delay_us(&info->start_time,
                                                                           &info->end_time);
                          cnt[1]++;

                        }
                    }
                }
              starpu_worker_get_name(worker, workername, 128);
              fprintf(stdout, "Worker %s:\n", workername);
              if (cnt[0] != 0)
                {
                  fprintf(stdout, "Avg. delay on XXTRF : %2.2lf us, %d tasks\n",
                          (delay_sum[0])/cnt[0], cnt[0]);
                  fprintf(stdout, "Avg. length on XXTRF : %2.2lf us\n",
                          (length_sum[0])/cnt[0]);
                }
              if (cnt[1] != 0)
                {
                  fprintf(stdout, "Avg. delay on XXMM : %2.2lf us, %d tasks\n",
                          (delay_sum[1])/cnt[1], cnt[1]);
                  fprintf(stdout, "Avg. length on XXMM : %2.2lf us\n",
                          (length_sum[1])/cnt[1]);
                }


              fprintf(stdout, "\ttotal time : %.2lf ms\n", total_time*1e-3);
              fprintf(stdout, "\texec time  : %.2lf ms (%.2f %%)\n", executing_time*1e-3, executing_ratio);
              fprintf(stdout, "\tblocked time  : %.2lf ms (%.2f %%)\n", sleeping_time*1e-3, sleeping_ratio);
            }
        }
#  endif /* USE_TASK_DEP && STARPU_PROFILING*/
    }
  for (itertask=0;itertask<SOLV_TASKNBR;itertask++)
    {
      INT itercblk = TASK_CBLKNUM(itertask);
      INT iterbloc;
#  ifdef USE_TASK_DEP
      for (iterbloc = SYMB_BLOKNUM(itercblk);
           iterbloc < SYMB_BLOKNUM(itercblk+1);
           iterbloc ++)
        {
          INT first_task = 0;
          if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
            {
              starpu_task_destroy(starpu_tasktab[iterbloc]);
              first_task = SYMB_BLOKNBR;
            }
          if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
            {
              starpu_task_destroy(starpu_tasktab[first_task+iterbloc]);
              first_task += SYMB_BLOKNBR;
#ifndef CHOL_SOPALIN
              first_task += SYMB_CBLKNBR;
#endif /* not CHOL_SOPALIN */
              starpu_task_destroy(starpu_tasktab[first_task+iterbloc]);
            }
        }
#    ifndef CHOL_SOPALIN
      {
        if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
          {
            INT first_task = SYMB_BLOKNBR;
            if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
              first_task += SYMB_BLOKNBR;
            starpu_task_destroy(starpu_tasktab[first_task+itercblk]);
          }
      }
#    endif /* CHOL_SOPALIN */
#  endif
      starpu_data_unregister(L_handle[itercblk]);
      if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
        {
          starpu_data_unregister(SM2X_handles[itercblk]);
        }
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
      starpu_data_unregister(U_handle[itercblk]);
#    endif
#  endif
    }
  starpu_data_unregister(WORK_handle);
#  ifdef USE_TASK_DEP
  memFree_null(starpu_tasktab);
#  endif

#  ifdef STARPU_BLOCKTAB_SELFCOPY
  {
    int ndevices;
    int device_id;
    CUDA_CALL(cudaGetDeviceCount(&ndevices));
    for (device_id = 0; device_id < ndevices; device_id++)
      {
        cudaSetDevice(device_id);
        CUDA_CALL(cudaFree(d_blocktab[device_id]));
      }
    memFree_null(d_blocktab);
  }
#  else
  starpu_data_unregister(blocktab_handle);
#  endif /* STARPU_BLOCKTAB_SELFCOPY */
  memFree_null(blocktab);

  /* Reduction on pivot number */
  sopalin_data->sopar->diagchange = 0;
  {
    INT me;
    for (me = 0; me < SOLV_THRDNBR+cuda_nbr; me++)
      {
        sopalin_data->sopar->diagchange += sopalin_data->thread_data[me]->nbpivot;
      }
  }
#  ifdef STARPU_USE_CUDA
  /*starpu_helper_cublas_shutdown();*/
#  endif

  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- sopalin time %lf\n",
            SOPALIN_CLOCK_GET);
  sopalin_data->sopar->dparm[DPARM_FACT_TIME] = SOPALIN_CLOCK_GET;
#  ifdef PASTIX_DUMP_FACTO
  dump_all(datacode, sopalin_data->sopar->cscmtx,
           ((datacode->updovct.sm2xtab!=NULL)?
            (DUMP_CSC | DUMP_SOLV | DUMP_SMB):(DUMP_CSC | DUMP_SOLV)));
#  endif
  memFree_null(trf_args);
  memFree_null(gemm_args);
  memFree_null(L_handle);
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    memFree_null(SM2X_handles);
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  memFree_null(U_handle);
#    endif /* CHOL_SOPALIN */
#  endif /* SOPALIN_LU   */
  starpu_shutdown();
  sopalin_clean(sopalin_data, 1);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR, starpu_clean_smp, sopalin_data,
                        0, NULL, NULL,
                        0, NULL, NULL);
  return NO_ERR;
}
#else
/* ISO C forbids an empty source file */
#  include "not_empty.h"
NOT_EMPTY(starpu_submit_tasks)
#endif
