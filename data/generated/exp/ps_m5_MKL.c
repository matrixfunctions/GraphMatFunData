#include<mkl/mkl.h>
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
    /* Computation order: B_0_1 B_1_1 B_2_1 A2 B_0_2 B_2_2 B_1_2 A3 B_2_3 B_0_3 B_1_3 A4 P2 C1 P1 C0 P0 */

    /* Computing B_0_1 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing B_1_1 with operation: lincomb */
    coeff1 = 0.041666666666666664;
    coeff2 = 0.008333333333333333;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing B_2_1 with operation: lincomb */
    coeff1 = 2.48015873015873e-5;
    coeff2 = 2.7557319223985893e-6;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing A2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[4], n);

    /* Computing B_0_2 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_0_1 */
    cblas_daxpby(n*n, coeff1, memslots[4], 1,
             coeff2, memslots[1], 1);

    /* Computing B_2_2 with operation: lincomb */
    coeff1 = 2.755731922398589e-7;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_2_1 */
    cblas_daxpby(n*n, coeff1, memslots[4], 1,
             coeff2, memslots[3], 1);

    /* Computing B_1_2 with operation: lincomb */
    coeff1 = 0.001388888888888889;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_1_1 */
    cblas_daxpby(n*n, coeff1, memslots[4], 1,
             coeff2, memslots[2], 1);

    /* Computing A3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[4], n,
            ZERO, memslots[5], n);
    /* Deallocating A2 in slot 5 */

    /* Computing B_2_3 with operation: lincomb */
    coeff1 = 2.505210838544172e-8;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_2_2 */
    cblas_daxpby(n*n, coeff1, memslots[5], 1,
             coeff2, memslots[3], 1);

    /* Computing B_0_3 with operation: lincomb */
    coeff1 = 0.16666666666666666;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_0_2 */
    cblas_daxpby(n*n, coeff1, memslots[5], 1,
             coeff2, memslots[1], 1);

    /* Computing B_1_3 with operation: lincomb */
    coeff1 = 0.0001984126984126984;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_1_2 */
    cblas_daxpby(n*n, coeff1, memslots[5], 1,
             coeff2, memslots[2], 1);

    /* Computing A4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[5], n,
            ZERO, memslots[4], n);
    /* Deallocating A in slot 1 */
    /* Deallocating A3 in slot 6 */

    /* Computing P2 with operation: lincomb */
    coeff1 = 2.08767569878681e-9;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_2_3 */
    cblas_daxpby(n*n, coeff1, memslots[4], 1,
             coeff2, memslots[3], 1);

    /* Computing C1 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[4], n,
            ZERO, memslots[0], n);
    /* Deallocating P2 in slot 4 */

    /* Computing P1 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle C1 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B_1_3 in slot 3 */

    /* Computing C0 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[4], n,
            ZERO, memslots[2], n);
    /* Deallocating P1 in slot 1 */
    /* Deallocating A4 in slot 5 */

    /* Computing P0 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle C0 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B_0_3 in slot 2 */

    /* Prepare output. */
    memcpy(output, memslots[2], n*n*sizeof(*output));
    free(master_mem);
}

