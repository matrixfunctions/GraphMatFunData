#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for polynomial evaluation. */
void dexp_native_jl_rho172_8(double *A, const size_t n, double *output) {
    size_t max_memalloc = 7;
    size_t max_memslots = 8;

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
    /* Computation order: C A2 Ua2 Va2 A4 Va3 Ua3 Ub2 Vb2 A6 Ua Va Ub3 Ub Uc U Vb3 Vb V Z X P S1 S2 S3 S4 S5 */

    /* Computing C with operation: lincomb */
    coeff1 = 0.03125;
    coeff2 = 0.0;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff1, memslots[0], 1);
    cblas_daxpby(n, coeff2, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing A2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[1], n);

    /* Computing Ua2 with operation: lincomb */
    coeff1 = 3.238237626624e16;
    coeff2 = 1.1873537964288e15;
    memcpy(memslots[2], memslots[1], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Va2 with operation: lincomb */
    coeff1 = 6.476475253248e16;
    coeff2 = 7.7717703038976e15;
    memcpy(memslots[3], memslots[1], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing A4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[1], n,
            ZERO, memslots[4], n);

    /* Computing Va3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.29060195264e14;
    /* Smart lincomb recycle Va2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing Ua3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.05594705216e13;
    /* Smart lincomb recycle Ua2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[2], 1);

    /* Computing Ub2 with operation: lincomb */
    coeff1 = 4.08408e7;
    coeff2 = 16380.0;
    memcpy(memslots[5], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[1], 1,
             coeff2, memslots[5], 1);

    /* Computing Vb2 with operation: lincomb */
    coeff1 = 1.32324192e9;
    coeff2 = 960960.0;
    memcpy(memslots[6], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[1], 1,
             coeff2, memslots[6], 1);

    /* Computing A6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[4], n,
            ZERO, memslots[7], n);
    /* Deallocating A2 in slot 2 */
    /* Deallocating A4 in slot 5 */

    /* Computing Ua with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.352212864e10;
    /* Smart lincomb recycle Ua3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[2], 1);

    /* Computing Va with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 6.704425728e11;
    /* Smart lincomb recycle Va3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing Ub3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ub2 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[5], 1);

    /* Computing Ub with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[7], n,
            ZERO, memslots[1], n);
    /* Deallocating Ub3 in slot 6 */

    /* Computing Uc with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ub */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);
    /* Deallocating Ua in slot 3 */

    /* Computing U with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[1], n,
            ZERO, memslots[2], n);
    /* Deallocating C in slot 1 */
    /* Deallocating Uc in slot 2 */

    /* Computing Vb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 182.0;
    /* Smart lincomb recycle Vb2 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[6], 1);

    /* Computing Vb with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[7], n,
            ZERO, memslots[0], n);
    /* Deallocating Vb3 in slot 7 */
    /* Deallocating A6 in slot 8 */

    /* Computing V with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Vb */
    cblas_daxpby(n*n, coeff2, memslots[3], 1,
             coeff1, memslots[0], 1);
    /* Deallocating Va in slot 4 */

    /* Computing Z with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.0;
    memcpy(memslots[1], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[1], 1);

    /* Computing X with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle V */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[0], 1);
    /* Deallocating U in slot 3 */

    /* Computing P with operation: ldiv */
    /* Reusing memory of Z for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[1], n, ipiv);
    /* Reusing memory of X for solution. */
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, n,
               memslots[1], n, ipiv,
               memslots[0], n);

    /* Computing S1 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating P in slot 1 */

    /* Computing S2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[1], n,
            ZERO, memslots[0], n);
    /* Deallocating S1 in slot 2 */

    /* Computing S3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating S2 in slot 1 */

    /* Computing S4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[1], n,
            ZERO, memslots[0], n);
    /* Deallocating S3 in slot 2 */

    /* Computing S5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating S4 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(ipiv);
    free(master_mem);
}

