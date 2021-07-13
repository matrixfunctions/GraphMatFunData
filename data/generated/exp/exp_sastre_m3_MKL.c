#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_sastre_m3(double *A, const size_t n, double *output) {
    size_t max_memalloc = 4;
    size_t max_memslots = 5;

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
    /* Computation order: T2k2 B2 T2k3 Bb3 Ba4_3 B3 Ba4 Bb4 B4 T2k4 T2k6 */

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[2], n);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 0.01992047682223989;
    coeff2 = 0.004980119205559973;
    memcpy(memslots[3], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[3], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 0.8765009801785554;
    coeff2 = 0.07665265321119147;
    /* Smart lincomb recycle A */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[0], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[3], n,
            ZERO, memslots[4], n);
    /* Deallocating Bb3 in slot 4 */

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[0], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 0.1225521150112075;
    coeff2 = 1.0;
    /* Smart lincomb recycle B2 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[2], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[2], n,
            ZERO, memslots[3], n);
    /* Deallocating Ba4 in slot 1 */
    /* Deallocating Bb4 in slot 3 */

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.974307204847627;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B3 in slot 5 */

    /* Computing T2k6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[3], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B4 in slot 4 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(master_mem);
}

