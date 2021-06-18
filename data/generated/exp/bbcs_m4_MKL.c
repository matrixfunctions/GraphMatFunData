#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void zbbcs_m4(MKL_Complex16 *A, const size_t n, MKL_Complex16 *output) {
    size_t max_memalloc = 7;
    size_t max_memslots = 8;

    /* Declarations and initializations. */
    MKL_Complex16 coeff1, coeff2;
    const MKL_Complex16 ZERO = {.real = 0.0, .imag = 0.0};
    const MKL_Complex16 ONE = {.real = 1.0, .imag = 0.0};
    size_t j;

    /* Memory management. */
    MKL_Complex16 *master_mem = malloc(n*n*max_memalloc*sizeof(*master_mem));
    MKL_Complex16 *memslots[max_memslots]; /* As many slots as nodes */
    memslots[0] = A; /* Overwrite A */
    /* The other slots are pointers to allocated memory. */
    for (j=0; j<max_memalloc; j++)
        memslots[j+1] = master_mem+j*n*n;
    /* Computation order: Ba5_2 Bb5_2 T2k2 B2 T2k3 Ba5_3 Bb5_3 Ba4_3 Bb4_3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 Bb5 B5 T2k7 */

    /* Computing Ba5_2 with operation: lincomb */
    coeff1.real = 2.6958430691533257;
    coeff1.imag = 0.0;
    coeff2.real = 0.0;
    coeff2.imag = 0.05272871327381115;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[1], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[1], n+1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1.real = 2.6958430691533257;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -1.3591092616886926;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[2], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1.real = -6.267569853502023;
    coeff1.imag = 0.0;
    coeff2.real = 0.0;
    coeff2.imag = 2.521796947120981;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[3], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[3], n+1);

    /* Computing B2 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[0], n, memslots[0], n,
            &ZERO, memslots[4], n);

    /* Computing T2k3 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 0.05786296656487002;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle T2k2 */
    cblas_zaxpby(n*n, &coeff2, memslots[4], 1,
             &coeff1, memslots[3], 1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.09896214548845832;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Ba5_2 */
    cblas_zaxpby(n*n, &coeff2, memslots[4], 1,
             &coeff1, memslots[1], 1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.09896214548845832;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Bb5_2 */
    cblas_zaxpby(n*n, &coeff2, memslots[4], 1,
             &coeff1, memslots[2], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1.real = 0.0;
    coeff1.imag = 0.13340427306445612;
    coeff2.real = 0.020226020298183107;
    coeff2.imag = 0.0;
    memcpy(memslots[5], memslots[4], n*n*sizeof(*master_mem));
    cblas_zaxpby(n*n, &coeff1, memslots[0], 1,
             &coeff2, memslots[5], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1.real = 0.0;
    coeff1.imag = 0.13340427306445612;
    coeff2.real = 0.020226020298183107;
    coeff2.imag = 0.0;
    memcpy(memslots[6], memslots[4], n*n*sizeof(*master_mem));
    cblas_zaxpby(n*n, &coeff1, memslots[0], 1,
             &coeff2, memslots[6], 1);

    /* Computing B3 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[4], n, memslots[0], n,
            &ZERO, memslots[7], n);
    /* Deallocating B2 in slot 5 */
    /* Deallocating A in slot 1 */

    /* Computing Ba4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -0.00674638241111651;
    /* Smart lincomb recycle Ba4_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[7], 1,
             &coeff1, memslots[5], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -0.00674638241111651;
    /* Smart lincomb recycle Bb4_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[7], 1,
             &coeff1, memslots[6], 1);

    /* Computing B4 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[5], n, memslots[6], n,
            &ZERO, memslots[0], n);
    /* Deallocating Ba4 in slot 6 */
    /* Deallocating Bb4 in slot 7 */

    /* Computing Ba5_4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 0.0;
    coeff2.imag = 0.007295441446830946;
    /* Smart lincomb recycle Ba5_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[7], 1,
             &coeff1, memslots[1], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 1.0;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Ba5_4 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[1], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -0.0776668640807187;
    /* Smart lincomb recycle T2k3 */
    cblas_zaxpby(n*n, &coeff2, memslots[7], 1,
             &coeff1, memslots[3], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 0.0;
    coeff2.imag = 0.015964794632994668;
    /* Smart lincomb recycle Bb5_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[7], 1,
             &coeff1, memslots[2], 1);
    /* Deallocating B3 in slot 8 */

    /* Computing Bb5 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 1.0;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Bb5_4 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[2], 1);
    /* Deallocating B4 in slot 1 */

    /* Computing B5 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[1], n, memslots[2], n,
            &ZERO, memslots[0], n);
    /* Deallocating Ba5 in slot 2 */
    /* Deallocating Bb5 in slot 3 */

    /* Computing T2k7 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 1.0;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle T2k4 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[3], 1);
    /* Deallocating B5 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

