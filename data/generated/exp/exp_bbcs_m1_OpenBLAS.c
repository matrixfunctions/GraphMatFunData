#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void zexp_bbcs_m1(openblas_complex_double *A, const size_t n, openblas_complex_double *output) {
    size_t max_memalloc = 2;
    size_t max_memslots = 3;

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
    /* Computation order: T2k2 B2 T2k4 */

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 1.0 + 0.0*I;
    coeff2 = -0.0 + -0.999999999976195*I;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_zscal(n*n, &coeff2, memslots[1], 1);
    cblas_zaxpby(n, &coeff1, &ONE, 0,
             &ONE, memslots[1], n+1);

    /* Computing B2 with operation: mult */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            &ONE, memslots[0], n, memslots[0], n,
            &ZERO, memslots[2], n);
    /* Deallocating A in slot 1 */

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0 + 0.0*I;
    coeff2 = -0.499999999992065 + 0.0*I;
    /* Smart lincomb recycle T2k2 */
    cblas_zaxpby(n*n, &coeff2, memslots[2], 1,
             &coeff1, memslots[1], 1);
    /* Deallocating B2 in slot 3 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(master_mem);
}

