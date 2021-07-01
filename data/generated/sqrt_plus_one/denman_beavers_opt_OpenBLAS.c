#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void ddenman_beavers_opt(double *A, const size_t n, double *output) {
    size_t max_memalloc = 7;
    size_t max_memslots = 8;

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
    /* Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7 */

    /* Computing Bb2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.011141376535527443;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5152503668289159;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.03845623901890617;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.7610413081657074;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing B2 with operation: ldiv */
    /* Reusing memory of Ba2 for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[4], n, ipiv);
    /* Reusing memory of Bb2 for solution. */
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, n,
               memslots[4], n, ipiv,
               memslots[1], n);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.10150239642416546;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.026744368550352622;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba3_2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.07127247848266706;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9541397430184381;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = 0.49421282325538585;
    coeff2 = 0.05943417585870563;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5017024742690791;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.12762852235945304;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0015101375211018255;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[6], 1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 0.28457285753903816;
    coeff2 = 0.11495283295141033;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0344930890633602;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.02575361014723645;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.004282806676350282;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 2 */

    /* Computing B3 with operation: ldiv */
    /* Reusing memory of Ba3 for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[4], n, ipiv);
    /* Reusing memory of Bb3 for solution. */
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, n,
               memslots[4], n, ipiv,
               memslots[0], n);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.016950835467818528;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.029043605652931636;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing B4 with operation: ldiv */
    /* Reusing memory of Ba4 for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[2], n, ipiv);
    /* Reusing memory of Bb4 for solution. */
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, n,
               memslots[2], n, ipiv,
               memslots[3], n);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5276435956563207;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[3], 1,
             coeff1, memslots[5], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.34606871085737384;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.14486471812413823;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[6], 1);
    /* Deallocating B3 in slot 1 */

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.006442438999085179;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[3], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01583234295751802;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[3], 1,
             coeff1, memslots[6], 1);
    /* Deallocating B4 in slot 4 */

    /* Computing B5 with operation: ldiv */
    /* Reusing memory of Ba5 for LU factors. */
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, memslots[5], n, ipiv);
    /* Reusing memory of Bb5 for solution. */
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, n,
               memslots[5], n, ipiv,
               memslots[6], n);

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.4670072325444539;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);
    /* Deallocating B5 in slot 7 */

    /* Prepare output. */
    memcpy(output, memslots[7], n*n*sizeof(*output));
    free(ipiv);
    free(master_mem);
}

