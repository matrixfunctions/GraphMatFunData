#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void zexp_bbcs_m5(MKL_Complex16 *A, const size_t n, MKL_Complex16 *output) {
    size_t max_memalloc = 6;
    size_t max_memslots = 7;

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
    /* Computation order: Ba6_2 Bb6_2 B2 Ba6_3 Bb6_3 T2k3 Ba5_3 B3 Bb6_4 Ba5_4 T2k4 Bb5_4 Ba6_4 B4 T2k5 Bb5 B5 Ba6_5 Bb6_5 Ba6 Bb6 B6 T2k8 */

    /* Computing Ba6_2 with operation: lincomb */
    coeff1.real = 0.3420232802536553;
    coeff1.imag = 0.0;
    coeff2.real = 0.0;
    coeff2.imag = -0.2851997796332415;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[1], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[1], n+1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1.real = 2.9237775839655367;
    coeff1.imag = 0.0;
    coeff2.real = 0.0;
    coeff2.imag = 1.4451330034748826;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[2], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[2], n+1);

    /* Computing B2 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[0], n, memslots[0], n,
            &ZERO, memslots[3], n);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 0.047347067331271094;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Ba6_2 */
    cblas_zaxpby(n*n, &coeff2, memslots[3], 1,
             &coeff1, memslots[1], 1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 0.1240818356655045;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Bb6_2 */
    cblas_zaxpby(n*n, &coeff2, memslots[3], 1,
             &coeff1, memslots[2], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1.real = 0.0;
    coeff1.imag = -0.6604084076077131;
    coeff2.real = -1.093022784715649;
    coeff2.imag = 0.0;
    memcpy(memslots[4], memslots[3], n*n*sizeof(*master_mem));
    cblas_zaxpby(n*n, &coeff1, memslots[0], 1,
             &coeff2, memslots[4], 1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1.real = 0.12;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -0.008774760968797039;
    memcpy(memslots[5], memslots[3], n*n*sizeof(*master_mem));
    cblas_zaxpby(n*n, &coeff1, memslots[0], 1,
             &coeff2, memslots[5], 1);

    /* Computing B3 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[3], n, memslots[0], n,
            &ZERO, memslots[6], n);
    /* Deallocating A in slot 1 */

    /* Computing Bb6_4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -0.019571570936427238;
    /* Smart lincomb recycle Bb6_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[6], 1,
             &coeff1, memslots[2], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0009784845352378095;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Ba5_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[6], 1,
             &coeff1, memslots[5], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 0.0;
    coeff2.imag = 0.2537715581771087;
    /* Smart lincomb recycle T2k3 */
    cblas_zaxpby(n*n, &coeff2, memslots[6], 1,
             &coeff1, memslots[4], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1.real = 0.0;
    coeff1.imag = -0.12395369585828313;
    coeff2.real = -0.011202694841085593;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle B2 */
    cblas_zaxpby(n*n, &coeff2, memslots[6], 1,
             &coeff1, memslots[3], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -0.022186600635366212;
    /* Smart lincomb recycle Ba6_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[6], 1,
             &coeff1, memslots[1], 1);

    /* Computing B4 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[6], n, memslots[6], n,
            &ZERO, memslots[0], n);
    /* Deallocating B3 in slot 7 */

    /* Computing T2k5 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 0.0005437426743473122;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle T2k4 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[4], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -1.2367240538259895e-5;
    /* Smart lincomb recycle Bb5_4 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[3], 1);

    /* Computing B5 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[5], n, memslots[3], n,
            &ZERO, memslots[6], n);
    /* Deallocating Ba5_4 in slot 6 */
    /* Deallocating Bb5 in slot 4 */

    /* Computing Ba6_5 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -9.74758985615379e-6;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Ba6_4 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[1], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 2.425253007433925e-5;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Bb6_4 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[2], 1);
    /* Deallocating B4 in slot 1 */

    /* Computing Ba6 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 1.0;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Ba6_5 */
    cblas_zaxpby(n*n, &coeff2, memslots[6], 1,
             &coeff1, memslots[1], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 1.0;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle Bb6_5 */
    cblas_zaxpby(n*n, &coeff2, memslots[6], 1,
             &coeff1, memslots[2], 1);
    /* Deallocating B5 in slot 7 */

    /* Computing B6 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[1], n, memslots[2], n,
            &ZERO, memslots[0], n);
    /* Deallocating Ba6 in slot 2 */
    /* Deallocating Bb6 in slot 3 */

    /* Computing T2k8 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 1.0;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle T2k5 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[4], 1);
    /* Deallocating B6 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[4], n*n*sizeof(*output));
    free(master_mem);
}

