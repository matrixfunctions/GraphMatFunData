#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_ps_m8(double *A, const size_t n, double *output) {
    size_t max_memalloc = 7;
    size_t max_memslots = 8;

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
    /* Computation order: B_3_1 B_0_1 B_4_1 B_1_1 B_2_1 A2 B_4_2 B_0_2 B_2_2 B_3_2 B_1_2 A3 B_2_3 B_3_3 B_0_3 B_4_3 B_1_3 A4 B_3_4 B_4_4 B_2_4 B_1_4 B_0_4 A5 P4 C3 P3 C2 P2 C1 P1 C0 P0 */

    /* Computing B_3_1 with operation: lincomb */
    coeff1 = 7.647163731819816e-13;
    coeff2 = 4.779477332387385e-14;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing B_0_1 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing B_4_1 with operation: lincomb */
    coeff1 = 4.110317623312165e-19;
    coeff2 = 1.9572941063391263e-20;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing B_1_1 with operation: lincomb */
    coeff1 = 0.008333333333333333;
    coeff2 = 0.001388888888888889;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing B_2_1 with operation: lincomb */
    coeff1 = 2.755731922398589e-7;
    coeff2 = 2.505210838544172e-8;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing A2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[6], n);

    /* Computing B_4_2 with operation: lincomb */
    coeff1 = 8.896791392450574e-22;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_4_1 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[3], 1);

    /* Computing B_0_2 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_0_1 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[2], 1);

    /* Computing B_2_2 with operation: lincomb */
    coeff1 = 2.08767569878681e-9;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_2_1 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[5], 1);

    /* Computing B_3_2 with operation: lincomb */
    coeff1 = 2.8114572543455206e-15;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_3_1 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[1], 1);

    /* Computing B_1_2 with operation: lincomb */
    coeff1 = 0.0001984126984126984;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_1_1 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[4], 1);

    /* Computing A3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[6], n,
            ZERO, memslots[7], n);
    /* Deallocating A2 in slot 7 */

    /* Computing B_2_3 with operation: lincomb */
    coeff1 = 1.6059043836821613e-10;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_2_2 */
    cblas_daxpby(n*n, coeff1, memslots[7], 1,
             coeff2, memslots[5], 1);

    /* Computing B_3_3 with operation: lincomb */
    coeff1 = 1.5619206968586225e-16;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_3_2 */
    cblas_daxpby(n*n, coeff1, memslots[7], 1,
             coeff2, memslots[1], 1);

    /* Computing B_0_3 with operation: lincomb */
    coeff1 = 0.16666666666666666;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_0_2 */
    cblas_daxpby(n*n, coeff1, memslots[7], 1,
             coeff2, memslots[2], 1);

    /* Computing B_4_3 with operation: lincomb */
    coeff1 = 3.868170170630684e-23;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_4_2 */
    cblas_daxpby(n*n, coeff1, memslots[7], 1,
             coeff2, memslots[3], 1);

    /* Computing B_1_3 with operation: lincomb */
    coeff1 = 2.48015873015873e-5;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_1_2 */
    cblas_daxpby(n*n, coeff1, memslots[7], 1,
             coeff2, memslots[4], 1);

    /* Computing A4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[7], n,
            ZERO, memslots[6], n);
    /* Deallocating A3 in slot 8 */

    /* Computing B_3_4 with operation: lincomb */
    coeff1 = 8.22063524662433e-18;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_3_3 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[1], 1);

    /* Computing B_4_4 with operation: lincomb */
    coeff1 = 1.6117375710961184e-24;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_4_3 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[3], 1);

    /* Computing B_2_4 with operation: lincomb */
    coeff1 = 1.1470745597729725e-11;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_2_3 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[5], 1);

    /* Computing B_1_4 with operation: lincomb */
    coeff1 = 2.7557319223985893e-6;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_1_3 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[4], 1);

    /* Computing B_0_4 with operation: lincomb */
    coeff1 = 0.041666666666666664;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_0_3 */
    cblas_daxpby(n*n, coeff1, memslots[6], 1,
             coeff2, memslots[2], 1);

    /* Computing A5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[6], n,
            ZERO, memslots[7], n);
    /* Deallocating A in slot 1 */
    /* Deallocating A4 in slot 7 */

    /* Computing P4 with operation: lincomb */
    coeff1 = 6.446950284384474e-26;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_4_4 */
    cblas_daxpby(n*n, coeff1, memslots[7], 1,
             coeff2, memslots[3], 1);

    /* Computing C3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[7], n,
            ZERO, memslots[0], n);
    /* Deallocating P4 in slot 4 */

    /* Computing P3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle C3 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B_3_4 in slot 2 */

    /* Computing C2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[7], n,
            ZERO, memslots[1], n);
    /* Deallocating P3 in slot 1 */

    /* Computing P2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle C2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B_2_4 in slot 6 */

    /* Computing C1 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[7], n,
            ZERO, memslots[0], n);
    /* Deallocating P2 in slot 2 */

    /* Computing P1 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle C1 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B_1_4 in slot 5 */

    /* Computing C0 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[7], n,
            ZERO, memslots[1], n);
    /* Deallocating P1 in slot 1 */
    /* Deallocating A5 in slot 8 */

    /* Computing P0 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle C0 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B_0_4 in slot 3 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(master_mem);
}

