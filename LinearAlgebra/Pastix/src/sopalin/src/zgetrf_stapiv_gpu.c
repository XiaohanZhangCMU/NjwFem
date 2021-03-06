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
    -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

       @precisions normal z -> s d c

*/
#ifndef CHOL_SOPALIN
#  define CHOL_SOPALIN
#endif /* CHOL_SOPALIN */

#ifndef SOPALIN_LU
#  define SOPALIN_LU
#endif /* SOPALIN_LU */

#include "common_pastix.h"
#include "sopalin_define.h"
#define inline static inline
#include "magma.h"
#undef inline

// === Define what BLAS to use ============================================
#define min                         MIN
#define max                         MAX
#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
#  else /* not PREC_DOUBLE */
#    define magma_zgemm             magma_cgemm
#    define magma_ztrsm             magma_ctrsm
#    define magma_zgetmatrix        magma_cgetmatrix
#    define magma_zsetmatrix        magma_csetmatrix
#    define magma_zmalloc_host      magma_cmalloc_host
#    undef  MAGMA_Z_ONE
#    undef  MAGMA_Z_NEG_ONE
#    define MAGMA_Z_ONE             MAGMA_C_ONE
#    define MAGMA_Z_NEG_ONE         MAGMA_C_NEG_ONE
#  endif /* not PREC_DOUBLE */
#else /* not TYPE_COMPLEX */
#  ifdef PREC_DOUBLE
#    define magma_zgemm             magmablas_dgemm
#    define magma_ztrsm             magmablas_dtrsm
#    define magma_zgetmatrix        magma_dgetmatrix
#    define magma_zsetmatrix        magma_dsetmatrix
#    define magma_zmalloc_host      magma_dmalloc_host
#    undef  MAGMA_Z_ONE
#    undef  MAGMA_Z_NEG_ONE
#    define MAGMA_Z_ONE             MAGMA_D_ONE
#    define MAGMA_Z_NEG_ONE         MAGMA_D_NEG_ONE
#  else /* not PREC_DOUBLE */
#    define magma_zgemm             magmablas_sgemm
#    define magma_ztrsm             magmablas_strsm
#    define magma_zgetmatrix        magma_sgetmatrix
#    define magma_zsetmatrix        magma_ssetmatrix
#    define magma_zmalloc_host      magma_smalloc_host
#    undef  MAGMA_Z_ONE
#    undef  MAGMA_Z_NEG_ONE
#    define MAGMA_Z_ONE             MAGMA_S_ONE
#    define MAGMA_Z_NEG_ONE         MAGMA_S_NEG_ONE
#  endif /* not PREC_DOUBLE */
#endif /* not TYPE_COMPLEX */

#include "pastix_cuda_helper.h"
#include "zgetrf_stapiv_gpu.h"

#define cuDoubleComplex CU_FLOAT

#define PASTIX_getrf_block API_CALL(PASTIX_getrf_block)
void PASTIX_getrf_block ( FLOAT *A, INT m, INT n, INT lda, INT *npvt, double crit );

// === End defining what BLAS to use =======================================


magma_int_t
magma_zgetrf_stapiv_gpu(magma_int_t m, magma_int_t n,
                        cuDoubleComplex *dA, magma_int_t ldda,
                        double criteria, INT * npiv,
                        magma_int_t *info)
{
/*  -- MAGMA (version 1.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       May 2012

    Purpose
    =======
    ZGETRF_STAPIV_GPU computes an LU factorization of a general M-by-N
    matrix A without any pivoting.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    Arguments
    =========
    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) COMPLEX_16 array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDDA     (input) INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    INFO    (output) INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.
    =====================================================================    */

#define inA(i,j) (dA + (i)*nb + (j)*nb*ldda)

    cuDoubleComplex c_one     = MAGMA_Z_ONE;
    cuDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;

    magma_int_t iinfo = 0, nb;
    magma_int_t maxm, maxn, mindim;
    magma_int_t i, rows, cols, s, lddwork;
    cuDoubleComplex *work;

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,m))
        *info = -4;

    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    /* Function Body */
    mindim = min(m, n);
    nb     = 2*magma_get_zgetrf_nb(m);
    s      = mindim / nb;

    if (nb <= 1 || nb >= min(m,n)) {
        /* Use CPU code. */
        work = (cuDoubleComplex*)malloc(m * n * sizeof(cuDoubleComplex));
        magma_zgetmatrix( m, n, dA, ldda, work, m );
        PASTIX_getrf_block((FLOAT*)work, m, n, m,
                           npiv,
                           criteria);
        magma_zsetmatrix( m, n, work, m, dA, ldda );
        free(work);
    }
    else {
        /* Use hybrid blocked code. */
        maxm = ((m + 31)/32)*32;
        maxn = ((n + 31)/32)*32;

        lddwork = maxm;

        if (MAGMA_SUCCESS != magma_zmalloc_host( &work, maxm*nb )) {
            *info = MAGMA_ERR_HOST_ALLOC;
            return *info;
        }

        for( i=0; i<s; i++ )
          {
            // download i-th panel
            cols = maxm - i*nb;
            magma_zgetmatrix( m-i*nb, nb, inA(i,i), ldda, work, lddwork );

            // make sure that gpu queue is empty
            magma_device_sync();

            if ( i>0 ){
              magma_ztrsm( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit,
                           nb, n - (i+1)*nb,
                           c_one, inA(i-1,i-1), ldda,
                           inA(i-1,i+1), ldda );
              magma_zgemm( MagmaNoTrans, MagmaNoTrans,
                           m-i*nb, n-(i+1)*nb, nb,
                           c_neg_one, inA(i,  i-1), ldda, inA(i-1,i+1), ldda,
                           c_one,     inA(i,  i+1), ldda );
            }

            // do the cpu part
            rows = m - i*nb;
            PASTIX_getrf_block((FLOAT*)work, rows, nb, lddwork,
                               npiv,
                               criteria);

            if ( (*info == 0) && (iinfo > 0) )
              *info = iinfo + i*nb;

            // upload i-th panel
            magma_zsetmatrix( m-i*nb, nb, work, lddwork, inA(i, i), ldda );

            // do the small non-parallel computations
            if ( s > (i+1) ) {
              magma_ztrsm( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit,
                           nb, nb,
                           c_one, inA(i, i  ), ldda,
                           inA(i, i+1), ldda);
              magma_zgemm( MagmaNoTrans, MagmaNoTrans,
                           m-(i+1)*nb, nb, nb,
                           c_neg_one, inA(i+1, i  ), ldda, inA(i,   i+1), ldda,
                           c_one,     inA(i+1, i+1), ldda );
            }
            else {
              magma_ztrsm( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit,
                           nb, n-s*nb,
                           c_one, inA(i, i  ), ldda,
                           inA(i, i+1), ldda);
              magma_zgemm( MagmaNoTrans, MagmaNoTrans,
                           m-(i+1)*nb, n-(i+1)*nb, nb,
                           c_neg_one, inA(i+1, i  ), ldda, inA(i,   i+1), ldda,
                           c_one,     inA(i+1, i+1), ldda );
            }
          }

        magma_int_t nb0 = min(m - s*nb, n - s*nb);
        rows = m - s*nb;
        cols = maxm - s*nb;
        magma_zgetmatrix( rows, nb0, inA(s,s), ldda, work, lddwork );

        // make sure that gpu queue is empty
        magma_device_sync();

        // do the cpu part
        PASTIX_getrf_block((FLOAT*)work, rows, nb0, lddwork,
                           npiv,
                           criteria);

        if ( (*info == 0) && (iinfo > 0) )
            *info = iinfo + s*nb;

        // upload i-th panel
        magma_zsetmatrix( rows, nb0, work, lddwork, inA(s,s), ldda );

        magma_ztrsm( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit,
                     nb0, n-s*nb-nb0,
                     c_one, inA(s,s),     ldda,
                            inA(s,s)+nb0, ldda);

        magma_free_host( work );
    }

    return *info;
} /* magma_zgetrf_stapiv_gpu */

#undef inA
