#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_sid_m10(double *A, const size_t n, double *output) {
    size_t max_memalloc = 10;
    size_t max_memslots = 11;

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
    /* Computation order: Bb9_2 Bb2 Ba9_2 Ba2 B2 Bb9_3 Ba9_3 Ba3 B3 Ba9_4 Bb9_4 Ba8_3 Ba8_4 Ba4 B4 Ba8_5 Ba9_5 Bb9_5 Ba5 B5 Bb9_6 Ba8_6 Ba9_6 Bb7_4 Bb7_5 Bb7_6 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba8_7 Bb7 Ba7_3 Bb8_3 Bb8_4 Ba7_4 Bb8_5 Ba7_5 Ba7_6 Bb8_6 Ba7 Bb8_7 B7 Ba8 B8 Ba9 Bb9 B9 B10 B11 T2k13 */

    /* Computing Bb9_2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.125;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Bb2 with operation: lincomb */
    coeff1 = 0.0;
    coeff2 = 0.125;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Ba9_2 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.125;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = 0.0;
    coeff2 = 0.125;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[2], n,
            ZERO, memslots[5], n);
    /* Deallocating Ba2 in slot 5 */
    /* Deallocating Bb2 in slot 3 */

    /* Computing Bb9_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5029302610017967;
    /* Smart lincomb recycle Bb9_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba9_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5029302610017967;
    /* Smart lincomb recycle Ba9_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 0.125;
    coeff2 = 0.0;
    memcpy(memslots[2], memslots[5], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[2], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[5], n,
            ZERO, memslots[4], n);
    /* Deallocating Ba3 in slot 3 */

    /* Computing Ba9_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.07705596948494946;
    /* Smart lincomb recycle Ba9_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb9_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.07705596948494946;
    /* Smart lincomb recycle Bb9_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba8_3 with operation: lincomb */
    coeff1 = 1.3883617476073524;
    coeff2 = 2.991654767354374;
    memcpy(memslots[2], memslots[5], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[2], 1);

    /* Computing Ba8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.2857950268422422;
    /* Smart lincomb recycle Ba8_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 0.125;
    coeff2 = 0.0;
    memcpy(memslots[6], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[6], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[4], n,
            ZERO, memslots[7], n);
    /* Deallocating Ba4 in slot 7 */

    /* Computing Ba8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.03005135891320298;
    /* Smart lincomb recycle Ba8_4 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba9_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.004985549176118462;
    /* Smart lincomb recycle Ba9_4 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb9_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.004985549176118462;
    /* Smart lincomb recycle Bb9_4 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 0.125;
    coeff2 = 0.0;
    memcpy(memslots[6], memslots[7], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[6], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[7], n,
            ZERO, memslots[8], n);
    /* Deallocating Ba5 in slot 7 */

    /* Computing Bb9_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 6.263526066651383e-5;
    /* Smart lincomb recycle Bb9_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.002742336655922557;
    /* Smart lincomb recycle Ba8_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba9_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 6.263526066651383e-5;
    /* Smart lincomb recycle Ba9_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 0.03306559506631931;
    coeff2 = 0.002630043177655382;
    memcpy(memslots[6], memslots[4], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[5], 1,
             coeff2, memslots[6], 1);

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0002100333647757715;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.7537107416419e-5;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.6421017771235788e-7;
    coeff2 = 6.204734935438909e-8;
    memcpy(memslots[9], memslots[5], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[9], 1);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.957106114715868e-9;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.556371639324141e-10;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.556371639324141e-11;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[9], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[8], n, memslots[9], n,
            ZERO, memslots[10], n);
    /* Deallocating Bb6 in slot 10 */

    /* Computing Ba8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 61.75954247606858;
    /* Smart lincomb recycle Ba8_6 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb7_6 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[6], 1);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 0.7439481579382581;
    coeff2 = 0.4155284057336423;
    memcpy(memslots[9], memslots[5], n*n*sizeof(*master_mem));
    cblas_daxpby(n*n, coeff1, memslots[0], 1,
             coeff2, memslots[9], 1);

    /* Computing Bb8_3 with operation: lincomb */
    coeff1 = -3.2977952779222e-5;
    coeff2 = 0.008139086096860678;
    /* Smart lincomb recycle A */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 6 */

    /* Computing Bb8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.001121744731945438;
    /* Smart lincomb recycle Bb8_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.02479095151834799;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B3 in slot 5 */

    /* Computing Bb8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 9.027588625491207e-5;
    /* Smart lincomb recycle Bb8_4 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.001283057135586989;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B4 in slot 8 */

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.501669195497238e-5;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 8.572383602707347e-6;
    /* Smart lincomb recycle Bb8_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B5 in slot 9 */

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb8_6 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B6 in slot 11 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[9], n, memslots[6], n,
            ZERO, memslots[4], n);
    /* Deallocating Ba7 in slot 10 */
    /* Deallocating Bb7 in slot 7 */

    /* Computing Ba8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba8_7 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[2], 1);

    /* Computing B8 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[0], n,
            ZERO, memslots[5], n);
    /* Deallocating Ba8 in slot 3 */
    /* Deallocating Bb8_7 in slot 1 */

    /* Computing Ba9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Ba9_6 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle Bb9_6 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[1], 1);
    /* Deallocating B8 in slot 6 */

    /* Computing B9 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[1], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba9 in slot 4 */
    /* Deallocating Bb9 in slot 2 */

    /* Computing B10 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating B9 in slot 1 */

    /* Computing B11 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[1], n,
            ZERO, memslots[0], n);
    /* Deallocating B10 in slot 2 */

    /* Computing T2k13 with operation: lincomb */
    coeff1 = 0.0;
    coeff2 = 1.0;
    /* Smart lincomb recycle B7 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);
    /* Deallocating B11 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[4], n*n*sizeof(*output));
    free(master_mem);
}

