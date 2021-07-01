#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void ddenman_beavers(double *A, const size_t n, double *output) {
    size_t max_memalloc = 2;
    size_t max_memslots = 3;

    /* Declarations and initializations. */
    double coeff1, coeff2;
    const double ONE = 1.0;
    lapack_int *ipiv = malloc(n*sizeof(*ipiv));
    size_t j;

    /* Memory management. */
    double *master_mem = malloc(n*n*max_memalloc*sizeof(*master_mem));
    double *memslots[max_memslots]; /* As many slots as nodes */
    memslots[0] = A; /* Overwrite A */
    /* The other slots are pointers to allocated memory. */
    for (j=0; j<max_memalloc; j++)
        memslots[j+1] = master_mem+j*n*n;
    /* Computation order: A_shift Xinv0 Y1 X1 Yinv1 X2 Xinv1 Y2 Yinv2 X3 */

    /* Computing A_shift with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff1, memslots[0], 1);
    cblas_daxpby(n, coeff2, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Xinv0 with operation: ldiv */
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[1], n, ipiv);
    LAPACKE_dgetri(LAPACK_COL_MAJOR, n, memslots[1], n, ipiv);

    /* Computing Y1 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 0.5;
    /* Smart lincomb recycle Xinv0 */
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing X1 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 0.5;
    /* Smart lincomb recycle A_shift */
    cblas_dscal(n*n, coeff1, memslots[0], 1);
    cblas_daxpby(n, coeff2, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Yinv1 with operation: ldiv */
    memcpy(memslots[2], memslots[1], n*n*sizeof(*master_mem));
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[2], n, ipiv);
    LAPACKE_dgetri(LAPACK_COL_MAJOR, n, memslots[2], n, ipiv);

    /* Computing X2 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 0.5;
    /* Smart lincomb recycle Yinv1 */
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[2], 1);

    /* Computing Xinv1 with operation: ldiv */
    /* Reusing memory of X1 for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[0], n, ipiv);
    LAPACKE_dgetri(LAPACK_COL_MAJOR, n, memslots[0], n, ipiv);

    /* Computing Y2 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 0.5;
    /* Smart lincomb recycle Y1 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[1], 1);
    /* Deallocating Xinv1 in slot 1 */

    /* Computing Yinv2 with operation: ldiv */
    /* Reusing memory of Y2 for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[1], n, ipiv);
    LAPACKE_dgetri(LAPACK_COL_MAJOR, n, memslots[1], n, ipiv);

    /* Computing X3 with operation: lincomb */
    coeff1 = 0.5;
    coeff2 = 0.5;
    /* Smart lincomb recycle X2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);
    /* Deallocating Yinv2 in slot 2 */

    /* Prepare output. */
    memcpy(output, memslots[2], n*n*sizeof(*output));
    free(ipiv);
    free(master_mem);
}

