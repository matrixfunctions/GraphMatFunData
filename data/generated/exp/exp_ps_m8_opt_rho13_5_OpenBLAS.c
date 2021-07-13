#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_ps_m8_opt_rho13_5(double *A, const size_t n, double *output) {
    size_t max_memalloc = 13;
    size_t max_memslots = 14;

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
    /* Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Ba4_2 Ba9_2 Bb9_2 Bb2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba7_3 Ba9_3 Ba3 B3 Ba7_4 Ba4 Bb4_3 Bb9_3 Bb8_4 Bb4 B4 Bb8_5 Ba6_3 Ba9_4 Ba9_5 T2k3 T2k4 T2k5 Bb9_4 Bb9_5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 Bb9_6 Ba6 Ba9_6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Bb9_7 Ba7 Ba9_7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba9_8 Bb9_8 Ba8 Bb8 T2k8 B8 Ba9 T2k9 Bb9 B9 T2k11 */

    /* Computing Ba3_2 with operation: lincomb */
    coeff1 = 0.09654809779701808;
    coeff2 = 0.4657361097626277;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = -2.47688825760804e-9;
    coeff2 = 0.0002423691784933533;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 0.03877713825381048;
    coeff2 = -9.60542079781634e-5;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Ba7_2 with operation: lincomb */
    coeff1 = 5.241383560289876e-5;
    coeff2 = 0.007243536894625013;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = 0.9062397883450928;
    coeff2 = 0.05669626301307128;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = -0.06178357797354583;
    coeff2 = 0.06738540542304026;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Ba9_2 with operation: lincomb */
    coeff1 = 0.5991807949797006;
    coeff2 = 0.15278813100979383;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb9_2 with operation: lincomb */
    coeff1 = 0.8838850843889047;
    coeff2 = 0.13220125487009166;
    memcpy(memslots[8], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[8], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[8], n+1);

    /* Computing Bb2 with operation: lincomb */
    coeff1 = 0.42031771622120523;
    coeff2 = 0.154626052690562;
    memcpy(memslots[9], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[9], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[9], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = 0.7648819021317048;
    coeff2 = 0.15551270544310056;
    memcpy(memslots[10], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[10], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[10], n+1);

    /* Computing Bb8_2 with operation: lincomb */
    coeff1 = 0.09171588796170577;
    coeff2 = 0.00978234693719379;
    memcpy(memslots[11], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[11], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[11], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = 0.42031771622120523;
    coeff2 = 0.154626052690562;
    memcpy(memslots[12], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[12], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[12], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[12], n, memslots[9], n,
            ZERO, memslots[13], n);
    /* Deallocating Ba2 in slot 13 */
    /* Deallocating Bb2 in slot 10 */

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.12302744324936216;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01590736865892108;
    /* Smart lincomb recycle Bb8_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[11], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.7329752440067754;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[6], 1);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01054251273835961;
    /* Smart lincomb recycle Ba7_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba9_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.3877431664280993;
    /* Smart lincomb recycle Ba9_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.19731022777550086;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[1], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[5], n,
            ZERO, memslots[9], n);
    /* Deallocating Ba3 in slot 2 */
    /* Deallocating Bb3 in slot 6 */

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.12590635459023522;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.5453304515439815;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.13872927228389606;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[10], 1);

    /* Computing Bb9_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.3359270098679286;
    /* Smart lincomb recycle Bb9_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[8], 1);

    /* Computing Bb8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.010821193755227018;
    /* Smart lincomb recycle Bb8_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[11], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.031072070853003377;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[10], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[10], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba4 in slot 7 */
    /* Deallocating Bb4 in slot 11 */

    /* Computing Bb8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.002975380215829149;
    /* Smart lincomb recycle Bb8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[11], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.008759615364757254;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba9_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.15662767864214414;
    /* Smart lincomb recycle Ba9_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba9_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.3206562914411985;
    /* Smart lincomb recycle Ba9_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[7], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0001512696450918944;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -7.050561377429653e-5;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -9.005499714502881e-5;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb9_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.13565211594408094;
    /* Smart lincomb recycle Bb9_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[8], 1);

    /* Computing Bb9_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.27797106138779576;
    /* Smart lincomb recycle Bb9_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[8], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.024971360914667298;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.20141416233600964;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.35076815913299725;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = -0.0037226968956481066;
    coeff2 = 0.1513971118834432;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.16951280499357504;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.7173684223758208;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.6151562256858113;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 0.4984122280039027;
    coeff2 = 0.3084467039511569;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.32464232547450916;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.11022586483948126;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.00562561032811143;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[6], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[6], n,
            ZERO, memslots[10], n);
    /* Deallocating Ba5 in slot 6 */
    /* Deallocating Bb5 in slot 7 */

    /* Computing Bb9_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.05519918645255441;
    /* Smart lincomb recycle Bb9_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[8], 1);

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9838034227912541;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba9_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.06358374798917776;
    /* Smart lincomb recycle Ba9_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9784008568588043;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.005577510585906947;
    /* Smart lincomb recycle Bb8_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[11], 1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -7.332920860083383e-8;
    coeff2 = 3.6275534419311994e-10;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.599720660809858e-11;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 5.386947932212966e-11;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.696432230206371e-12;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.7604767541951467e-14;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[5], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[5], n,
            ZERO, memslots[6], n);
    /* Deallocating Ba6 in slot 3 */
    /* Deallocating Bb6 in slot 6 */

    /* Computing Bb9_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.39661780513350103;
    /* Smart lincomb recycle Bb9_6 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[8], 1);

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.756683002701541e-7;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba9_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.4636979253428603;
    /* Smart lincomb recycle Ba9_6 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.335908979618945;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -7.0435154500906185e-6;
    /* Smart lincomb recycle Bb8_6 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[11], 1);

    /* Computing Ba8_2 with operation: lincomb */
    coeff1 = 0.9900217594371604;
    coeff2 = -0.02132723543710295;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Ba8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0895038182070465;
    /* Smart lincomb recycle Ba8_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.4332397850086141;
    /* Smart lincomb recycle Ba8_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.6230790774477177;
    /* Smart lincomb recycle Ba8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.836798518267107;
    /* Smart lincomb recycle Ba8_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 6.715240076479665e-7;
    /* Smart lincomb recycle Ba8_6 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb7_2 with operation: lincomb */
    coeff1 = -0.004732216406789163;
    coeff2 = 1.2338220978385536e-5;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.805471807800756e-5;
    /* Smart lincomb recycle Bb7_2 */
    cblas_daxpby(n*n, coeff2, memslots[13], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 14 */

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 6.965400723989983e-6;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[9], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B3 in slot 10 */

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.4212591399083807e-6;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B4 in slot 2 */

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.2164967481538e-7;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[10], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B5 in slot 11 */

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0000000002008125;
    /* Smart lincomb recycle Bb7_6 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B6 in slot 7 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba7 in slot 5 */
    /* Deallocating Bb7 in slot 1 */

    /* Computing Ba9_8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.005883392744500455;
    /* Smart lincomb recycle Ba9_7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb9_8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.007845590149896397;
    /* Smart lincomb recycle Bb9_7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[8], 1);

    /* Computing Ba8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.03410244422083362;
    /* Smart lincomb recycle Ba8_7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0199299719677057;
    /* Smart lincomb recycle Bb8_7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[11], 1);

    /* Computing T2k8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0029363288440026106;
    /* Smart lincomb recycle T2k7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B7 in slot 2 */

    /* Computing B8 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[11], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba8 in slot 3 */
    /* Deallocating Bb8 in slot 12 */

    /* Computing Ba9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0870379046019574;
    /* Smart lincomb recycle Ba9_8 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);

    /* Computing T2k9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.00026219521556040433;
    /* Smart lincomb recycle T2k8 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9425352289352358;
    /* Smart lincomb recycle Bb9_8 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[8], 1);
    /* Deallocating B8 in slot 1 */

    /* Computing B9 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[7], n, memslots[8], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba9 in slot 8 */
    /* Deallocating Bb9 in slot 9 */

    /* Computing T2k11 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0181552414422717;
    /* Smart lincomb recycle T2k9 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B9 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

