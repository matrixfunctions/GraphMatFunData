#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dsid_m5(double *A, const size_t n, double *output) {
    size_t max_memalloc = 7;
    size_t max_memslots = 8;

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
    /* Computation order: T2k2 B2 T2k3 Ba6_3 Ba5_3 Bb4_3 Bb6_3 B3 Bb4 Bb6_4 Ba5_4 T2k4 Bb5_4 Ba6_4 B4 Ba5 T2k5 Bb5 B5 Ba6_5 Bb6_5 Ba6 Bb6 B6 T2k6 T2k8 */

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
    coeff2 = 0.2863243726334417;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 0.3112216227982407;
    coeff2 = 0.0852839259083158;
    memcpy(memslots[3], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[3], 1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 0.9418613214806352;
    coeff2 = 0.06974348269544424;
    memcpy(memslots[4], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[4], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 5.374708803114821e-5;
    coeff2 = 4.50085273957301e-6;
    memcpy(memslots[5], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[5], 1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 0.6865706355662834;
    coeff2 = 0.1392249143769798;
    memcpy(memslots[6], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[6], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[2], n,
            ZERO, memslots[7], n);
    /* Deallocating A in slot 1 */

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.16165883444488e-6;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.03151382711608315;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[6], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.002005403977292901;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[4], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.08776036732867759;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = -0.007544837153586671;
    coeff2 = 0.002852960512714315;
    /* Smart lincomb recycle B2 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0292447258748138;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[7], n, memslots[5], n,
            ZERO, memslots[0], n);
    /* Deallocating B3 in slot 8 */
    /* Deallocating Bb4 in slot 6 */

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.18995526739487723;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[2], n,
            ZERO, memslots[5], n);
    /* Deallocating Ba5 in slot 5 */
    /* Deallocating Bb5 in slot 3 */

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.829773504500424;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 11.173624766438472;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[6], 1);
    /* Deallocating B4 in slot 1 */

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[6], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[6], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba6 in slot 4 */
    /* Deallocating Bb6 in slot 7 */

    /* Computing T2k6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.23337016308538;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B5 in slot 6 */

    /* Computing T2k8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B6 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(master_mem);
}

