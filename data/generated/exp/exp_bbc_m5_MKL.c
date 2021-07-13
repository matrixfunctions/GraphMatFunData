#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_bbc_m5(double *A, const size_t n, double *output) {
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
    /* Computation order: Ba6_2 Bb6_2 B2 Ba6_3 Bb6_3 T2k3 Ba5_3 B3 Bb6_4 Ba5_4 T2k4 Bb5_4 Ba6_4 B4 T2k5 Bb5 B5 Ba6_5 Bb6_5 Ba6 Bb6 B6 T2k8 */

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = -11.058071288535288;
    coeff2 = 1.6125176868819238;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -0.09043168323908106;
    coeff2 = -0.06764045190713819;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[3], n);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.12477411482493252;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[3], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.06759613017704597;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[3], 1,
             coeff1, memslots[2], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 0.3978497494996451;
    coeff2 = 1.3678377846041172;
    memcpy(memslots[4], memslots[3], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[4], 1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = -0.10036558103014462;
    coeff2 = -0.00802924648241157;
    memcpy(memslots[5], memslots[3], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[5], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[3], n,
            ZERO, memslots[6], n);
    /* Deallocating A in slot 1 */

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.029555257042931552;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.00089213849804573;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[5], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.49828962252538267;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = -0.09233646193671186;
    coeff2 = -0.016936493900208172;
    /* Smart lincomb recycle B2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.02257315581805103;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[1], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[6], n,
            ZERO, memslots[0], n);
    /* Deallocating B3 in slot 7 */

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0006378981945947233;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.400867981820361e-5;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[3], n,
            ZERO, memslots[6], n);
    /* Deallocating Ba5_4 in slot 6 */
    /* Deallocating Bb5 in slot 4 */

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.9579475957000978e-5;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.391802575160607e-5;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B4 in slot 1 */

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[1], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B5 in slot 7 */

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[2], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba6 in slot 2 */
    /* Deallocating Bb6 in slot 3 */

    /* Computing T2k8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);
    /* Deallocating B6 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[4], n*n*sizeof(*output));
    free(master_mem);
}

