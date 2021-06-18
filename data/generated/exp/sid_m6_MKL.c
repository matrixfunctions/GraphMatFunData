#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dsid_m6(double *A, const size_t n, double *output) {
    size_t max_memalloc = 8;
    size_t max_memslots = 9;

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
    /* Computation order: T2k2 B2 T2k3 Ba6_3 Bb7_3 Bb5_3 Ba7_3 B3 Bb6_4 Bb7_4 T2k4 Bb5_4 Ba6_4 Ba7_4 B4 T2k5 Bb5 Ba6_5 Bb6_5 Bb7_5 Ba7_5 B5 Ba6 Bb7_6 Ba7_6 Bb6 B6 Ba7 B7 T2k9 */

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
    coeff2 = 0.1969779342112314;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[2], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 2.865001388641538;
    coeff2 = 0.1952545843107103;
    memcpy(memslots[3], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[3], 1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 0.03655234395347475;
    coeff2 = 0.01606091400855144;
    memcpy(memslots[4], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[4], 1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 2.294895435403922e-5;
    coeff2 = 1.406952242413849e-6;
    memcpy(memslots[5], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[5], 1);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 8.29008575139441;
    coeff2 = 1.739158441630994;
    memcpy(memslots[6], memslots[2], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[6], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[2], n,
            ZERO, memslots[7], n);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 0.02721930992200371;
    coeff2 = 0.002547056607231984;
    /* Smart lincomb recycle B2 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0001758035313846159;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[4], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.03005000525808178;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 9.379681616092325e-8;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01430688980356062;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.1965098904519709;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[6], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[7], n,
            ZERO, memslots[8], n);
    /* Deallocating A in slot 1 */
    /* Deallocating B3 in slot 8 */

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.002243394407902074;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.172460202011541e-8;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.002024281516007681;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.001204349003694297;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0002919349464582001;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.02018492049443954;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[6], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[8], n, memslots[5], n,
            ZERO, memslots[0], n);
    /* Deallocating B4 in slot 9 */
    /* Deallocating Bb5 in slot 6 */

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 249.896909254999;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B5 in slot 1 */

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[2], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba6 in slot 4 */
    /* Deallocating Bb6 in slot 3 */

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[6], 1);
    /* Deallocating B6 in slot 1 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[4], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba7 in slot 7 */
    /* Deallocating Bb7_6 in slot 5 */

    /* Computing T2k9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B7 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[1], n*n*sizeof(*output));
    free(master_mem);
}

