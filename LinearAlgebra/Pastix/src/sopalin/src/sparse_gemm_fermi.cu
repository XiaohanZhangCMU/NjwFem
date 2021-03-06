/*
  -- MAGMA (version 1.1) --
  Univ. of Tennessee, Knoxville
  Univ. of California, Berkeley
  Univ. of Colorado, Denver
  November 2011


  @precisions normal z -> z c d s
       
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cuda.h>

#include <starpu.h> 
#include <starpu_cuda.h>

#include "sparse_gemm_fermi.h"

#ifdef PRECISION_z
#include "zgemm_fermi_define.h"
#endif
#ifdef PRECISION_c
#include "cgemm_fermi_define.h"
#endif
#ifdef PRECISION_d
#include "dgemm_fermi_define.h"
#endif
#ifdef PRECISION_s
#include "sgemm_fermi_define.h"
#endif

extern "C" void
GENERATE_SM_VERSION_NAME(gemm)( char TRANSA, char TRANSB, int m , int n , int k ,
				cuDoubleComplex alpha, const cuDoubleComplex *d_A, int lda,
                                                       const cuDoubleComplex *d_B, int ldb,
				cuDoubleComplex beta,        cuDoubleComplex *d_C, int ldc,
				int blocknbr, const int *blocktab, int fblocknbr, const int *fblocktab,
				cudaStream_t stream )
{
    /*  -- MAGMA (version 1.1) --
        Univ. of Tennessee, Knoxville
        Univ. of California, Berkeley
        Univ. of Colorado, Denver
        November 2011

        Purpose
        =======
        ZGEMM  performs one of the matrix-matrix operations

        C := alpha*op( A )*op( B ) + beta*C,

        where  op( X ) is one of

        op( X ) = X   or   op( X ) = X',

        alpha and beta are scalars, and A, B and C are matrices, with op( A )
        an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

        Parameters
        ==========
        TRANSA - CHARACTER*1.
        On entry, TRANSA specifies the form of op( A ) to be used in
        the matrix multiplication as follows:
        TRANSA = 'N' or 'n',  op( A ) = A.
        TRANSA = 'T' or 't',  op( A ) = A'.
        TRANSA = 'C' or 'c',  op( A ) = A'.
        Unchanged on exit.

        TRANSB - CHARACTER*1.
        On entry, TRANSB specifies the form of op( B ) to be used in
        the matrix multiplication as follows:
        TRANSB = 'N' or 'n',  op( B ) = B.
        TRANSB = 'T' or 't',  op( B ) = B'.
        TRANSB = 'C' or 'c',  op( B ) = B'.
        Unchanged on exit.

        M      - INTEGER.
        On entry,  M  specifies  the number  of rows  of the  matrix
        op( d_A )  and of the  matrix d_C.  M  must  be at least  zero.
        Unchanged on exit.

        N      - INTEGER.
        On entry,  N  specifies the number  of columns of the matrix
        op( d_B ) and the number of columns of the matrix d_C. N must be
        at least zero.
        Unchanged on exit.

        K      - INTEGER.
        On entry,  K  specifies  the number of columns of the matrix
        op( d_A ) and the number of rows of the matrix op( d_B ). K must
        be at least  zero.
        Unchanged on exit.

        ALPHA  - COMPLEX_16
        On entry, ALPHA specifies the scalar alpha.
        Unchanged on exit.

        d_A    - COMPLEX_16 array of DIMENSION ( LDA, ka ), where ka is
        k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
        Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
        part of the array d_A must contain the matrix d_A, otherwise
        the leading  k by m  part of the array d_A must contain  the
        matrix d_A.
        Unchanged on exit.

        LDA    - INTEGER.
        On entry, LDA specifies the first dimension of A as declared
        in the calling (sub) program. When  TRANSA = 'N' or 'n' then
        LDA must be at least  max( 1, m ), otherwise  LDA must be at
        least  max( 1, k ).
        Unchanged on exit.

        d_B    - COMPLEX_16 array of DIMENSION ( LDB, kb ), where kb is
        n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
        Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
        part of the array d_B must contain the matrix d_B, otherwise
        the leading  n by k  part of the array d_B must contain  the
        matrix d_B.
        Unchanged on exit.
 
        LDB    - INTEGER.
        On entry, LDB specifies the first dimension of d_B as declared
        in the calling (sub) program. When  TRANSB = 'N' or 'n' then
        LDB must be at least  max( 1, k ), otherwise  LDB must be at
        least  max( 1, n ).
        Unchanged on exit.

        BETA   - COMPLEX_16.
        On entry,  BETA  specifies the scalar  beta.  When  BETA  is
        supplied as zero then d_C need not be set on input.
        Unchanged on exit.

        d_C    - COMPLEX_16 array of DIMENSION ( LDC, n ).
        Before entry, the leading  m by n  part of the array  d_C must
        contain the matrix  d_C,  except when  beta  is zero, in which
        case d_C need not be set on entry.
        On exit, the array  d_C  is overwritten by the  m by n  matrix
        ( alpha*op( d_A )*op( d_B ) + beta*d_C ).

        LDC    - INTEGER.
        On entry, LDC specifies the first dimension of d_C as declared
        in  the  calling  (sub)  program.   LDC  must  be  at  least
        max( 1, m ).
        Unchanged on exit.
        =====================================================================    */
    if (m<=0 || n<=0 || k<=0)
        return;

    size_t offsetA = 0;
    size_t offsetB = 0;
#if defined(PRECISION_z) || defined(PRECISION_c) 
    int TransA = 2, TransB = 2;
#else
    int TransA = 1, TransB = 1;
#endif
    if (TRANSA == 'T' ||  TRANSA == 't')
        TransA = 1;
    else
        if (TRANSA == 'N' ||  TRANSA == 'n')
            TransA = 0;
    
    if (TRANSB == 'T' ||  TRANSB == 't')
        TransB = 1;
    else
        if (TRANSB == 'N' ||  TRANSB == 'n')
            TransB = 0;

    size_t sizeA = (size_t) lda * (size_t) (!TransA ? k : m);
    size_t sizeB = (size_t) ldb * (size_t) (!TransB ? n : k);

    /* TODO: Check with Jakub what is this */
    size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512);
#if 0
    if (sizeA>=CUBLAS_MAX_1DBUF_SIZE ||
        sizeB>=CUBLAS_MAX_1DBUF_SIZE )
        {
            cublasZgemm(TRANSA, TRANSB, m, n, k, alpha,
                        d_A, lda, d_B, ldb,
                        beta, d_C, ldc);
            return;
        }
#else
    if (sizeA>=CUBLAS_MAX_1DBUF_SIZE ||
        sizeB>=CUBLAS_MAX_1DBUF_SIZE )
        {
            fprintf(stderr, "ERROR: The matrix size is too big to use texture\n");
            return;
        }

#endif

#ifdef TEXTURE_1D
    // Set textures parameters
    tex_ref_A.normalized = false;
    tex_ref_A.filterMode = cudaFilterModePoint;
    tex_ref_A.addressMode[0] = cudaAddressModeClamp;
    
    tex_ref_B.normalized = false;
    tex_ref_B.filterMode = cudaFilterModePoint;
    tex_ref_B.addressMode[0] = cudaAddressModeClamp;
    
    // Bind A and B to texture references
    assert(cudaBindTexture(&offsetA, tex_ref_A, d_A, sizeA*sizeof(cuDoubleComplex)) 
           == cudaSuccess);
    assert(cudaBindTexture(&offsetB, tex_ref_B, d_B, sizeB*sizeof(cuDoubleComplex))
           == cudaSuccess);
#endif

    // Set up grids
    // Warning: works because DIM_X and DIM_Y are equals for every cases of one precision
    dim3 dimBlock(DIM_X, DIM_Y);

    offsetA = offsetA/sizeof(d_A[0]);
    offsetB = offsetB/sizeof(d_B[0]);
    
    if (TransA==0 && TransB ==0){
        dim3 dimGrid(m/BLK_M_nn + (m%BLK_M_nn != 0),
                     n/BLK_N_nn + (n%BLK_N_nn != 0));
        GENERATE_SM_VERSION_NAME(gemm_nn)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    } 
    else if (TransA==0 && TransB ==1){
        dim3 dimGrid(m/BLK_M_nt + (m%BLK_M_nt != 0),
                     n/BLK_N_nt + (n%BLK_N_nt != 0));
        GENERATE_SM_VERSION_NAME(gemm_nt)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    }
    else if (TransA==1 && TransB ==0){
        dim3 dimGrid(m/BLK_M_tn + (m%BLK_M_tn != 0),
                     n/BLK_N_tn + (n%BLK_N_tn != 0));
        GENERATE_SM_VERSION_NAME(gemm_tn)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    }
    else if (TransA==1 && TransB ==1){
        dim3 dimGrid(m/BLK_M_tt + (m%BLK_M_tt != 0),
                     n/BLK_N_tt + (n%BLK_N_tt != 0));
        GENERATE_SM_VERSION_NAME(gemm_tt)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    }
#if defined(PRECISION_z) || defined(PRECISION_c) 
    else if (TransA==0 && TransB ==2){
        dim3 dimGrid(m/BLK_M_nt + (m%BLK_M_nt != 0),
                     n/BLK_N_nt + (n%BLK_N_nt != 0));
        GENERATE_SM_VERSION_NAME(gemm_nc)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    } 
    else if (TransA==1 && TransB ==2){
        dim3 dimGrid(m/BLK_M_tt + (m%BLK_M_tt != 0),
                     n/BLK_N_tt + (n%BLK_N_tt != 0));
        GENERATE_SM_VERSION_NAME(gemm_tc)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    }
    else if (TransA==2 && TransB ==0){
        dim3 dimGrid(m/BLK_M_tn + (m%BLK_M_tn != 0),
                     n/BLK_N_tn + (n%BLK_N_tn != 0));
        GENERATE_SM_VERSION_NAME(gemm_cn)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    }
    else if (TransA==2 && TransB ==1){
        dim3 dimGrid(m/BLK_M_tt + (m%BLK_M_tt != 0),
                     n/BLK_N_tt + (n%BLK_N_tt != 0));
        GENERATE_SM_VERSION_NAME(gemm_ct)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    } 
    else if (TransA==2 && TransB ==2){
        dim3 dimGrid(m/BLK_M_tt + (m%BLK_M_tt != 0),
                     n/BLK_N_tt + (n%BLK_N_tt != 0));
        GENERATE_SM_VERSION_NAME(gemm_cc)<<< dimGrid, dimBlock, 0, stream >>>(m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc,
                                                                         (int)offsetA, (int)offsetB, blocknbr, blocktab, fblocknbr, fblocktab);
    }
#endif
    else {
      fprintf(stderr, "ERROR: in GEMM kernel");
      assert(0);
    }
    cudaUnbindTexture ( tex_ref_A ) ;
    cudaUnbindTexture ( tex_ref_B ) ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
