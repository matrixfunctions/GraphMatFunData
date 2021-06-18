#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_native_jl_rho0_95(double *A, const size_t n, double *output) {
    size_t max_memalloc = 5;
    size_t max_memslots = 6;

    /* Declarations and initializations. */
    double coeff1, coeff2;
    const double ZERO = 0.0;
    const double ONE = 1.0;
    lapack_int *ipiv = malloc(n*sizeof(*ipiv));
    size_t j;

    /* Memory management. */
    double *master_mem = malloc(n*n*max_memalloc*sizeof(*master_mem));
    double *memslots[max_memslots]; /* As many slots as nodes */
    memslots[0] = A; /* Overwrite A */
    /* The other slots are pointers to allocated memory. */
    for (j=0; j<max_memalloc; j++)
        memslots[j+1] = master_mem+j*n*n;
    /* Computation order: A2 Ua2 V2 A4 Ua3 V3 A6 Ua U V Z X P */

    /* Computing A2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[1], n);

    /* Computing Ua2 with operation: lincomb */
    coeff1 = 8.64864e6;
    coeff2 = 277200.0;
    memcpy(memslots[2], memslots[1], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing V2 with operation: lincomb */
    coeff1 = 1.729728e7;
    coeff2 = 1.99584e6;
    memcpy(memslots[3], memslots[1], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing A4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[1], n,
            ZERO, memslots[4], n);

    /* Computing Ua3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1512.0;
    /* Smart lincomb recycle Ua2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[2], 1);

    /* Computing V3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 25200.0;
    /* Smart lincomb recycle V2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing A6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[4], n,
            ZERO, memslots[5], n);
    /* Deallocating A2 in slot 2 */
    /* Deallocating A4 in slot 5 */

    /* Computing Ua with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ua3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[2], 1);

    /* Computing U with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating Ua in slot 3 */
    /* Deallocating A in slot 1 */

    /* Computing V with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 56.0;
    /* Smart lincomb recycle V3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);
    /* Deallocating A6 in slot 6 */

    /* Computing Z with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.0;
    memcpy(memslots[0], memslots[1], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[3], 1,
             coeff2, memslots[0], 1);

    /* Computing X with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle V */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);
    /* Deallocating U in slot 2 */

    /* Computing P with operation: ldiv */
    /* Reusing memory of Z for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[0], n, ipiv);
    /* Reusing memory of X for solution. */
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, n,
               memslots[0], n, ipiv,
               memslots[3], n);

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(ipiv);
    free(master_mem);
}

