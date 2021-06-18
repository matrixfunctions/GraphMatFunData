#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void zbbcs_m3(openblas_complex_double *A, const size_t n, openblas_complex_double *output) {
    size_t max_memalloc = 4;
    size_t max_memslots = 5;

    /* Declarations and initializations. */
    openblas_complex_double coeff1, coeff2;
    const openblas_complex_double ZERO = 0.0 + 0.0*I;
    const openblas_complex_double ONE = 1.0 + 0.0*I;
    size_t j;

    /* Memory management. */
    openblas_complex_double *master_mem = malloc(n*n*max_memalloc*sizeof(*master_mem));
    openblas_complex_double *memslots[max_memslots]; /* As many slots as nodes */
    memslots[0] = A; /* Overwrite A */
    /* The other slots are pointers to allocated memory. */
    for (j=0; j<max_memalloc; j++)
        memslots[j+1] = master_mem+j*n*n;
    /* Computation order: Bb4_2 T2k2 B2 Bb3 T2k3 Bb4_3 B3 Ba4 Bb4 B4 T2k6 */

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = 0.0 + 0.5496085391143601*I;
    coeff2 = 0.1620095284677366 + 0.0*I;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[1], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[1], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 1.0 + 0.0*I;
    coeff2 = -0.0 + -0.9999999999999923*I;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[2], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[2], n+1);

    /* Computing B2 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[0], n, memslots[0], n,
            &ZERO, memslots[3], n);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 0.10775 + 0.0*I;
    coeff2 = -0.0 + -0.026939068735988708*I;
    /* Smart lincomb recycle A */
    cblas_zaxpby(n*n, &coeff2, memslots[3], 1,
             &coeff1, memslots[0], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0 + 0.0*I;
    coeff2 = -0.13549409636220702 + 0.0*I;
    /* Smart lincomb recycle T2k2 */
    cblas_zaxpby(n*n, &coeff2, memslots[3], 1,
             &coeff1, memslots[2], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0 + 0.0*I;
    coeff2 = -0.0 + -0.014179818052118045*I;
    /* Smart lincomb recycle Bb4_2 */
    cblas_zaxpby(n*n, &coeff2, memslots[3], 1,
             &coeff1, memslots[1], 1);

    /* Computing B3 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[3], n, memslots[0], n,
            &ZERO, memslots[4], n);
    /* Deallocating Bb3 in slot 1 */

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 0.0 + 0.6632100444166243*I;
    coeff2 = 1.0 + 0.0*I;
    /* Smart lincomb recycle B2 */
    cblas_zaxpby(n*n, &coeff2, memslots[4], 1,
             &coeff1, memslots[3], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0 + 0.0*I;
    coeff2 = -0.034159539168921116 + 0.0*I;
    /* Smart lincomb recycle Bb4_3 */
    cblas_zaxpby(n*n, &coeff2, memslots[4], 1,
             &coeff1, memslots[1], 1);
    /* Deallocating B3 in slot 5 */

    /* Computing B4 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[3], n, memslots[1], n,
            &ZERO, memslots[0], n);
    /* Deallocating Ba4 in slot 4 */
    /* Deallocating Bb4 in slot 2 */

    /* Computing T2k6 with operation: lincomb */
    coeff1 = 1.0 + 0.0*I;
    coeff2 = 1.0 + 0.0*I;
    /* Smart lincomb recycle T2k3 */
    cblas_zaxpby(n*n, &coeff2, memslots[0], 1,
             &coeff1, memslots[2], 1);
    /* Deallocating B4 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[2], n*n*sizeof(*output));
    free(master_mem);
}

