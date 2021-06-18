#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dsid_m4(double *A, const size_t n, double *output) {
    size_t max_memalloc = 6;
    size_t max_memslots = 7;

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
    /* Computation order: T2k2 B2 T2k3 Bb3 B3 T2k4 Bb4 Ba4_3 Ba4 B4 T2k5 Ba5_3 Bb5_3 Ba5_4 Bb5_4 Ba5 Bb5 B5 T2k7 */

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
    coeff2 = 0.5918659857804601;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 0.002945531440279683;
    coeff2 = 0.0004018761610201036;
    memcpy(memslots[3], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[3], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[3], n,
            ZERO, memslots[4], n);
    /* Deallocating Bb3 in slot 4 */

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -3.2733920099600837;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 0.03230762888122312;
    coeff2 = 1.0;
    memcpy(memslots[3], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[2], 1,
             coeff2, memslots[3], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 0.4017568440673568;
    coeff2 = -0.008709066576837676;
    memcpy(memslots[5], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[5], 1);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[5], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[3], n,
            ZERO, memslots[6], n);
    /* Deallocating Ba4 in slot 6 */
    /* Deallocating Bb4 in slot 4 */

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 10.408017352313541;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 2.224209172496374;
    coeff2 = 0.2614927977298117;
    memcpy(memslots[3], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[3], 1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = -0.04130276365929783;
    coeff2 = 0.02338576034271299;
    /* Smart lincomb recycle A */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 3 */

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 5.768988513026145;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.023373194047115575;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B3 in slot 5 */

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B4 in slot 7 */

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[0], n,
            ZERO, memslots[2], n);
    /* Deallocating Ba5 in slot 4 */
    /* Deallocating Bb5 in slot 1 */

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B5 in slot 3 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(master_mem);
}

