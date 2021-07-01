#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dps(double *A, const size_t n, double *output) {
    size_t max_memalloc = 5;
    size_t max_memslots = 6;

    /* Declarations and initializations. */
    double coeff1, coeff2;
    const double ZERO = 0.0;
    const double ONE = 1.0;
    size_t j;

    /* Memory management. */
    double *master_mem = malloc(n*n*max_memalloc*sizeof(*master_mem));
    double *memslots[max_memslots]; /* As many slots as nodes */
    memslots[0] = A; /* Overwrite A */
    /* The other slots are pointers to allocated memory. */
    for (j=0; j<max_memalloc; j++)
        memslots[j+1] = master_mem+j*n*n;
    /* Computation order: Bb4_2 Bb5_2 T2k2 B2 T2k3 Bb4_3 Bb5_3 B3 Bb4 B4 Bb5 B5 T2k7 */

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = -0.0205078125;
    coeff2 = 0.01611328125;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 0.0625;
    coeff2 = -0.0390625;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[4], n);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.125;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.013092041015625;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.02734375;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[2], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[0], n,
            ZERO, memslots[5], n);
    /* Deallocating B2 in slot 5 */
    /* Deallocating A in slot 1 */

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0109100341796875;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[1], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[1], n,
            ZERO, memslots[0], n);
    /* Deallocating Bb4 in slot 2 */

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B4 in slot 1 */

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[2], n,
            ZERO, memslots[0], n);
    /* Deallocating B3 in slot 6 */
    /* Deallocating Bb5 in slot 3 */

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B5 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

