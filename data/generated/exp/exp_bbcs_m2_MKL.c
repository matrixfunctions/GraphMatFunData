#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void zexp_bbcs_m2(MKL_Complex16 *A, const size_t n, MKL_Complex16 *output) {
    size_t max_memalloc = 3;
    size_t max_memslots = 4;

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
    /* Computation order: T2k2 B2 Bb3 B3 T2k3 T2k5 */

    /* Computing T2k2 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.0;
    coeff2.imag = -0.9999999999998107;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[1], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[1], n+1);

    /* Computing B2 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[0], n, memslots[0], n,
            &ZERO, memslots[2], n);

    /* Computing Bb3 with operation: lincomb */
    coeff1.real = 0.0;
    coeff1.imag = 0.16666657785001893;
    coeff2.real = 0.04166664890333649;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle A */
    cblas_zaxpby(n*n, &coeff2, memslots[2], 1,
             &coeff1, memslots[0], 1);

    /* Computing B3 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[2], n, memslots[0], n,
            &ZERO, memslots[3], n);
    /* Deallocating Bb3 in slot 1 */

    /* Computing T2k3 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = -0.4999999999999432;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle T2k2 */
    cblas_zaxpby(n*n, &coeff2, memslots[2], 1,
             &coeff1, memslots[1], 1);
    /* Deallocating B2 in slot 3 */

    /* Computing T2k5 with operation: lincomb */
    coeff1.real = 1.0;
    coeff1.imag = 0.0;
    coeff2.real = 1.0;
    coeff2.imag = 0.0;
    /* Smart lincomb recycle T2k3 */
    cblas_zaxpby(n*n, &coeff2, memslots[3], 1,
             &coeff1, memslots[1], 1);
    /* Deallocating B3 in slot 4 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(master_mem);
}

