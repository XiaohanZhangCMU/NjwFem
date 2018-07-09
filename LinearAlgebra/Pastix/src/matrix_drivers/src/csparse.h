#define CS_VER 1    /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
  INT nzmax ;    /* maximum number of entries */
  INT m ;    /* number of rows */
  INT n ;    /* number of columns */
  INT *p ;    /* column poINTers (size n+1) or col indices (size nzmax) */
  INT *i ;    /* row indices, size nzmax */
  FLOAT *x ;    /* numerical values, size nzmax */
  INT nz ;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

/* keep all large entries */
INT cs_droptol (cs *A, FLOAT tol) ;
/* keep all nonzero entries */
INT cs_dropzeros (cs *A) ;
/* C = alpha*A + beta*B */
cs *cs_add (const cs *A, const cs *B, FLOAT alpha, FLOAT beta) ;
/* removes duplicate entries from A */
INT cs_dupl (cs *A) ;
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
INT cs_entry (cs *T, INT i, INT j, FLOAT x) ;
/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
INT cs_fkeep (cs *A, INT (*fkeep) (INT, INT, FLOAT, void *), void *other) ;
/* y = A*x+y */
INT cs_gaxpy (const cs *A, const FLOAT *x, FLOAT *y) ;
/* C = A*B */
cs *cs_multiply (const cs *A, const cs *B) ;
/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
FLOAT cs_norm (const cs *A) ;
/* C = A(P,Q) where P and Q are permutations of 0..m-1 and 0..n-1. */
cs *cs_permute (const cs *A, const INT *P, const INT *Q, INT values) ;
/* Pinv = P', or P = Pinv' */
INT *cs_pinv (const INT *P, INT n) ;
/* C = A' */
cs *cs_transpose (const cs *A, INT values) ;
/* C = compressed-column form of a triplet matrix T */
cs *cs_triplet (const cs *T) ;
/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
INT cs_scatter (const cs *A, INT j, FLOAT beta, INT *w, FLOAT *x, INT mark,
                cs *C, INT nz) ;
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
INT cs_cumsum (INT *p, INT *c, INT n) ;

/* utilities */
/* wrapper for malloc */
void *cs_malloc (INT n, size_t size) ;
/* wrapper for calloc */
void *cs_calloc (INT n, size_t size) ;
/* wrapper for free */
void *cs_free (void *p) ;
/* wrapper for realloc */
void *cs_realloc (void *p, INT n, size_t size, INT *ok) ;
/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *cs_spalloc (INT m, INT n, INT nzmax, INT values, INT triplet) ;
/* change the max # of entries sparse matrix */
INT cs_sprealloc (cs *A, INT nzmax) ;
/* free a sparse matrix */
cs *cs_spfree (cs *A) ;
/* free workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, INT ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (INT) size)
