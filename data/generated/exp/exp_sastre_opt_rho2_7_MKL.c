#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_sastre_opt_rho2_7(double *A, const size_t n, double *output) {
    size_t max_memalloc = 12;
    size_t max_memslots = 13;

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
    /* Computation order: Ba3_2 Ba6_2 T2k2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Ba3 Ba6_3 T2k3 Bb4_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 Ba7_2 Ba7_3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Bb5 B5 T2k6 Bb6 Ba6_4 Ba6_5 Ba6 B6 T2k7 Bb7_2 Bb7_3 Bb7_4 Ba7_4 Bb7_5 Ba7_5 Bb7_6 Ba7_6 Ba7 Bb7 B7 T2k9 */

    /* Computing Ba3_2 with operation: lincomb */
    coeff1 = 0.00203407200644078;
    coeff2 = -0.004855370488778173;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = 0.015636244039123465;
    coeff2 = 0.36835598195824903;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 0.9720224032451862;
    coeff2 = 0.08866505998167876;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Bb2 with operation: lincomb */
    coeff1 = 0.0014414740602762641;
    coeff2 = 1.015909177999167;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = -0.00216500404580565;
    coeff2 = 0.09646644408880453;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = -0.029613578842783986;
    coeff2 = 0.377669985575162;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = 0.0014414740602762641;
    coeff2 = 1.015909177999167;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[7], n, memslots[4], n,
            ZERO, memslots[8], n);
    /* Deallocating Ba2 in slot 8 */
    /* Deallocating Bb2 in slot 5 */

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.014294207683808999;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0000472141119134;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.04415666392494207;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.00025013219458514967;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.05304182885989784;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[6], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = -0.0022569437385435;
    coeff2 = 0.0010373403739564698;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.018857479770905;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = -0.0167704319125257;
    coeff2 = 0.0018989172164601057;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 7.889567871549961e-5;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba7_2 with operation: lincomb */
    coeff1 = 3.1043463304552805;
    coeff2 = 0.7528521412951402;
    memcpy(memslots[9], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[9], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[9], n+1);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.159971928559395;
    /* Smart lincomb recycle Ba7_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = 0.020097851537399642;
    coeff2 = 0.10667353947597173;
    memcpy(memslots[10], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[10], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[10], n+1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.034594555518270974;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[10], 1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = -0.0068121857554928915;
    coeff2 = 0.0018059085648684517;
    memcpy(memslots[11], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[11], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[11], n+1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.00010342034415488765;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[11], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[11], n,
            ZERO, memslots[12], n);
    /* Deallocating Ba3 in slot 2 */
    /* Deallocating Bb3 in slot 12 */

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.1172277610690486;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.022330965215494;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[6], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[6], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba4 in slot 6 */
    /* Deallocating Bb4 in slot 7 */

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.07438766844651382;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[10], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.02327178284493455;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.08243100849197117;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.10144129673946056;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[10], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.021251100930175442;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.13164782415883697;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0009931016008313642;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -2.623385669729595e-5;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[7], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[7], n,
            ZERO, memslots[5], n);
    /* Deallocating Ba5 in slot 5 */
    /* Deallocating Bb5 in slot 8 */

    /* Computing T2k6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0031660405732008064;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.987202584803334;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[10], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.019415313402916037;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.054241015390228235;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0061285566950167;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[2], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[10], n,
            ZERO, memslots[4], n);
    /* Deallocating Ba6 in slot 3 */
    /* Deallocating Bb6 in slot 11 */

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.015815360491620606;
    /* Smart lincomb recycle T2k6 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb7_2 with operation: lincomb */
    coeff1 = 0.008831439952196395;
    coeff2 = 0.29016677501097143;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.07731361344655992;
    /* Smart lincomb recycle Bb7_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 9 */

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01972134822932684;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 7.905233218278997;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B3 in slot 13 */

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.2595231983751941;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0995871761977756;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B4 in slot 2 */

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.9401869208092082;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.048631444054980266;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B5 in slot 6 */

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.2841742967744319;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9656182105520791;
    /* Smart lincomb recycle Bb7_6 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B6 in slot 5 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[9], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba7 in slot 10 */
    /* Deallocating Bb7 in slot 1 */

    /* Computing T2k9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9769074674369014;
    /* Smart lincomb recycle T2k7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B7 in slot 2 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

