#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_sastre_m8(double *A, const size_t n, double *output) {
    size_t max_memalloc = 8;
    size_t max_memslots = 9;

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
    /* Computation order: Ba9_2 T2k2 Ba8_2 B2 Ba8_3 T2k3 Ba9_3 Bb6_3 Ba7_3 B3 Bb6_4 Ba9_4 Bb7_4 T2k4 Ba8_4 Ba7_4 B4 Ba8_5 Ba9_5 T2k5 Bb6_5 Bb7_5 Ba7_5 B5 Ba8_6 Bb7_6 Ba7_6 Bb6 B6 Ba8_7 Ba7 Bb7 B7 Ba8 B8 Ba9 B9 T2k11 */

    /* Computing Ba9_2 with operation: lincomb */
    coeff1 = 0.008333333333333333;
    coeff2 = 0.001388888888888889;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Ba8_2 with operation: lincomb */
    coeff1 = 2.755731922398589e-7;
    coeff2 = 2.505210838544172e-8;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[4], n);

    /* Computing Ba8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.08767569878681e-9;
    /* Smart lincomb recycle Ba8_2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba9_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0001984126984126984;
    /* Smart lincomb recycle Ba9_2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = -1.023660713518307e-11;
    coeff2 = -4.508311519886735e-13;
    memcpy(memslots[5], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[5], 1);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = -5.893435534477677e-5;
    coeff2 = -3.013961104055248e-6;
    memcpy(memslots[6], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[6], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[4], n,
            ZERO, memslots[7], n);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.980157255925737e-14;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba9_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.48015873015873e-5;
    /* Smart lincomb recycle Ba9_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = -5.100472475630675e-7;
    coeff2 = -4.032817333361947e-8;
    /* Smart lincomb recycle B2 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[4], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.16666666666666666;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.30531132637709e-10;
    /* Smart lincomb recycle Ba8_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.502070379373464e-7;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[6], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[7], n,
            ZERO, memslots[8], n);
    /* Deallocating B3 in slot 8 */

    /* Computing Ba8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 7.55676813469492e-12;
    /* Smart lincomb recycle Ba8_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba9_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.7557319223985893e-6;
    /* Smart lincomb recycle Ba9_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[1], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.041666666666666664;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -9.210033748491798e-16;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -2.785084196756015e-9;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -6.770221628797445e-9;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[6], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[8], n,
            ZERO, memslots[7], n);
    /* Deallocating A in slot 1 */
    /* Deallocating B4 in slot 9 */

    /* Computing Ba8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.024189993755686e-13;
    /* Smart lincomb recycle Ba8_5 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -3.294026127901678e-10;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.227011356117036e-10;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -6.140022498994532e-17;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[5], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[7], n, memslots[5], n,
            ZERO, memslots[0], n);
    /* Deallocating Bb6 in slot 6 */

    /* Computing Ba8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.001023463999572971;
    /* Smart lincomb recycle Ba8_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb7_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);
    /* Deallocating B6 in slot 1 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[4], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba7 in slot 7 */
    /* Deallocating Bb7 in slot 5 */

    /* Computing Ba8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba8_7 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B7 in slot 1 */

    /* Computing B8 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[7], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba8 in slot 4 */

    /* Computing Ba9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba9_5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B8 in slot 1 */

    /* Computing B9 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[7], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba9 in slot 2 */
    /* Deallocating B5 in slot 8 */

    /* Computing T2k11 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B9 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[2], n*n*sizeof(*output));
    free(master_mem);
}

