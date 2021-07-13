#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_sastre_m6(double *A, const size_t n, double *output) {
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
    /* Computation order: Ba7_2 B2 Ba7_3 Bb3 B3 Ba7_4 Bb4 Ba4_3 Ba4 B4 Ba7_5 Ba6_3 Bb7_3 Bb5_3 B5 Ba6 Bb7_6 Bb6 B6 Bb7 B7 T2k9 */

    /* Computing Ba7_2 with operation: lincomb */
    coeff1 = 3.09646797193604;
    coeff2 = 0.772360321294401;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[2], n);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.1673139636901279;
    /* Smart lincomb recycle Ba7_2 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 0.001052151783051235;
    coeff2 = 0.0004675683454147702;
    memcpy(memslots[3], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[3], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[3], n,
            ZERO, memslots[4], n);
    /* Deallocating Bb3 in slot 4 */

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 7.922322450524197;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 0.05317514832355802;
    coeff2 = 1.0;
    memcpy(memslots[3], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[2], 1,
             coeff2, memslots[3], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 0.2868706220817633;
    coeff2 = -0.03289442879547955;
    memcpy(memslots[5], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[5], 1);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[5], 1);
    /* Deallocating B3 in slot 5 */

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[3], n,
            ZERO, memslots[4], n);
    /* Deallocating Ba4 in slot 6 */
    /* Deallocating Bb4 in slot 4 */

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B4 in slot 5 */

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 0.39689859154115;
    coeff2 = 0.02219811707032801;
    memcpy(memslots[3], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[3], 1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 0.3229486011362678;
    coeff2 = 0.08092036376147299;
    memcpy(memslots[4], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[4], 1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 0.002688394980266927;
    coeff2 = 0.0004675683454147702;
    /* Smart lincomb recycle A */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[0], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[0], n,
            ZERO, memslots[5], n);
    /* Deallocating Bb5_3 in slot 1 */

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.930814505527068;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 0.0277140002806296;
    coeff2 = 1.0;
    /* Smart lincomb recycle B2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B5 in slot 6 */

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[2], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba6 in slot 4 */
    /* Deallocating Bb6 in slot 3 */

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb7_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);
    /* Deallocating B6 in slot 1 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[4], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba7_5 in slot 2 */
    /* Deallocating Bb7 in slot 5 */

    /* Computing T2k9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle B7 */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Prepare output. */
    memcpy(output, memslots[0], n*n*sizeof(*output));
    free(master_mem);
}

