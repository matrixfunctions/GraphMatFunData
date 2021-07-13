#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dps_m5(double *A, const size_t n, double *output) {
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
    /* Computation order: Bb5_2 T2k2 Bb6_2 B2 T2k3 Bb6_3 Bb5_3 B3 Bb6_4 T2k4 Bb5_4 B4 Bb5 B5 Bb6 B6 T2k8 */

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = -0.013092041015625;
    coeff2 = 0.0109100341796875;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -0.0390625;
    coeff2 = 0.02734375;
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
             coeff1, memslots[2], 1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0205078125;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.009273529052734375;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[0], n,
            ZERO, memslots[5], n);
    /* Deallocating B2 in slot 5 */

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01611328125;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0625;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.008008956909179688;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[1], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[0], n,
            ZERO, memslots[4], n);
    /* Deallocating B3 in slot 6 */
    /* Deallocating A in slot 1 */

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0070078372955322266;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[1], n,
            ZERO, memslots[0], n);
    /* Deallocating Bb5 in slot 2 */

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B5 in slot 1 */

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[3], n,
            ZERO, memslots[0], n);
    /* Deallocating B4 in slot 5 */
    /* Deallocating Bb6 in slot 4 */

    /* Computing T2k8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B6 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[2], n*n*sizeof(*output));
    free(master_mem);
}

