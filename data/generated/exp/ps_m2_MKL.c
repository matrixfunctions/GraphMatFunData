#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for polynomial evaluation. */
void dps_m2(double *A, const size_t n, double *output) {
    size_t max_memalloc = 3;
    size_t max_memslots = 4;

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
    /* Computation order: B_0_1 B_1_1 A2 P1 C0 P0 */

    /* Computing B_0_1 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing B_1_1 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 0.16666666666666666;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing A2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[3], n);
    /* Deallocating A in slot 1 */

    /* Computing P1 with operation: lincomb */
    coeff1 = 0.041666666666666664;
    coeff2 = 1.0;
    /* Smart lincomb recycle B_1_1 */
    cblas_daxpby(n*n, coeff1, memslots[3], 1,
             coeff2, memslots[2], 1);

    /* Computing C0 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[3], n,
            ZERO, memslots[0], n);
    /* Deallocating P1 in slot 3 */
    /* Deallocating A2 in slot 4 */

    /* Computing P0 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle C0 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B_0_1 in slot 2 */

    /* Prepare output. */
    memcpy(output, memslots[0], n*n*sizeof(*output));
    free(master_mem);
}

