#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_bbc_opt_rho1_9(double *A, const size_t n, double *output) {
    size_t max_memalloc = 9;
    size_t max_memslots = 10;

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
    /* Computation order: Ba3_2 Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba3 Ba4_3 Bb4_3 Ba6_2 Ba6_3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb6_2 Bb6_3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Bb6_4 Ba5_4 Ba5 Bb6_5 T2k4 T2k5 Bb5_4 Ba6_4 Bb5 Ba6_5 B5 Ba6 Bb6 B6 T2k6 T2k8 */

    /* Computing Ba3_2 with operation: lincomb */
    coeff1 = -0.38075288501512605;
    coeff2 = 0.8627222226009553;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Bb2 with operation: lincomb */
    coeff1 = -0.4242607265172309;
    coeff2 = 1.1950771657897332;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = 5.9466829480949973e-5;
    coeff2 = -0.0034728192074554225;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = 5.9466829480949973e-5;
    coeff2 = -0.0034728192074554225;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = -0.4242607265172309;
    coeff2 = 1.1950771657897332;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[2], n,
            ZERO, memslots[6], n);
    /* Deallocating Ba2 in slot 6 */
    /* Deallocating Bb2 in slot 3 */

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.021223247723999996;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.26578830898107336;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.26578830898107336;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = -11.089873841786654;
    coeff2 = 1.5695850557582407;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.10777934093048946;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = 0.015468999442833017;
    coeff2 = -0.24257817849158173;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.004371767791483258;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 0.025702668821420974;
    coeff2 = 0.038656450236714296;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.2493161704221847;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = -0.00923583715392215;
    coeff2 = 0.35247666711675785;
    memcpy(memslots[8], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[8], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[8], n+1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.598102175541146;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[8], 1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -0.0688394655881371;
    coeff2 = -0.10086875671987898;
    memcpy(memslots[9], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[9], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[9], n+1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.03935278313404004;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = -0.021792993546751255;
    coeff2 = -0.3404755645503017;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.876565481553418;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 7 */

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[0], n,
            ZERO, memslots[6], n);
    /* Deallocating Ba3 in slot 2 */
    /* Deallocating Bb3 in slot 1 */

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9731527852522692;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9731527852522692;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[4], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[3], n, memslots[4], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba4 in slot 4 */
    /* Deallocating Bb4 in slot 5 */

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.05834755130517858;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[9], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.00033737278567153267;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.0818706572379708e-7;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.1857830431851822e-5;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.4923984017940168;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[8], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -4.8910943480715535e-5;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[8], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.010629231499670367;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0547023784540487;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B3 in slot 7 */

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -2.2644581317176334e-6;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.1843377663664171e-5;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);
    /* Deallocating B4 in slot 1 */

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[7], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba5 in slot 6 */
    /* Deallocating Bb5 in slot 8 */

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9049940153322246;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.144156658763187;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[9], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[9], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba6 in slot 3 */
    /* Deallocating Bb6 in slot 10 */

    /* Computing T2k6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.037592213904829044;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[8], 1);
    /* Deallocating B5 in slot 1 */

    /* Computing T2k8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.1421298417742292;
    /* Smart lincomb recycle T2k6 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[8], 1);
    /* Deallocating B6 in slot 2 */

    /* Prepare output. */
    memcpy(output, memslots[8], n*n*sizeof(*output));
    free(master_mem);
}

