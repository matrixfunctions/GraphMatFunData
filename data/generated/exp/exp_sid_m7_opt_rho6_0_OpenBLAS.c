#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_sid_m7_opt_rho6_0(double *A, const size_t n, double *output) {
    size_t max_memalloc = 11;
    size_t max_memslots = 12;

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
    /* Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Bb2 Ba4_2 Bb4_2 Bb8_2 Ba2 B2 Bb3 Bb8_3 Ba4_3 Ba3 B3 Ba7_3 Ba7_4 Ba4 Ba6_3 Bb8_4 T2k3 T2k4 Bb4_3 Bb4 B4 Bb8_5 T2k5 Ba6_4 Ba6_5 Ba7_5 Ba5_2 Ba5_3 Ba5_4 Ba5 Bb5_2 Bb5_3 Bb5_4 Bb5 B5 T2k6 Ba6 Ba7_6 Bb8_6 Bb6_2 Bb6_3 Bb6_4 Bb6_5 Bb6 B6 Ba7 T2k7 Bb8_7 Ba8_2 Ba8_3 Ba8_4 Ba8_5 Ba8_6 Ba8_7 Bb7_2 Bb7_3 Bb7_4 Bb7_5 Bb7_6 Bb7 B7 Ba8 Bb8 B8 T2k8 T2k10 */

    /* Computing Ba3_2 with operation: lincomb */
    coeff1 = -0.05790075880684695;
    coeff2 = 1.050095804216935;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = 2.8512702775691487e-6;
    coeff2 = 8.323858094672198e-6;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 0.9943018399475028;
    coeff2 = 0.9623400714214185;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Ba7_2 with operation: lincomb */
    coeff1 = -0.0006322010382109642;
    coeff2 = 5.966629547283877;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = -0.000951838725229865;
    coeff2 = -0.035990141876148504;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Bb2 with operation: lincomb */
    coeff1 = -0.0674550728089079;
    coeff2 = 1.1639373098746155;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = -0.022273287604033215;
    coeff2 = 0.999419730173801;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = -0.0003564977821414783;
    coeff2 = 0.006790605847297012;
    memcpy(memslots[8], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[8], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[8], n+1);

    /* Computing Bb8_2 with operation: lincomb */
    coeff1 = 0.013739397806079797;
    coeff2 = -0.034167351013704314;
    memcpy(memslots[9], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[9], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[9], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = -0.0674550728089079;
    coeff2 = 1.1639373098746155;
    memcpy(memslots[10], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[10], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[10], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[10], n, memslots[6], n,
            ZERO, memslots[11], n);
    /* Deallocating Ba2 in slot 11 */
    /* Deallocating Bb2 in slot 7 */

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0505773512582828;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.00017439769298246415;
    /* Smart lincomb recycle Bb8_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[9], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.04428302676636615;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.02211891467634549;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[1], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[5], n,
            ZERO, memslots[6], n);
    /* Deallocating Ba3 in slot 2 */
    /* Deallocating Bb3 in slot 6 */

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.3119218262051087;
    /* Smart lincomb recycle Ba7_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01788182758740769;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0013986894432222274;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.3339497680252084e-5;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0006895671006188639;
    /* Smart lincomb recycle Bb8_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.621285209008396;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.22263803654723124;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0863934611339156;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[8], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9986696579293114;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[8], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[7], n, memslots[8], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba4 in slot 8 */
    /* Deallocating Bb4 in slot 9 */

    /* Computing Bb8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.0262187841524036e-5;
    /* Smart lincomb recycle Bb8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.04001562250761491;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0005347112859929923;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.009928519927039425;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0005913130477887621;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = -0.00992506274997938;
    coeff2 = 1.000247094423774;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.010965295621586402;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0007441489467712597;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.3836127966022704e-6;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = -0.00010751494099623147;
    coeff2 = 0.00045481726705312207;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0010096879534240845;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.027269903163941744;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.999827547479686;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[7], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[5], n, memslots[7], n,
            ZERO, memslots[8], n);
    /* Deallocating Ba5 in slot 6 */
    /* Deallocating Bb5 in slot 8 */

    /* Computing T2k6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0017973977825017301;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0000008787568422;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.309674736291994e-5;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.93368482319488e-5;
    /* Smart lincomb recycle Bb8_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -3.9152959471223135e-5;
    coeff2 = 8.334641611255522e-8;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 7.3222974116279275e-9;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.8759928836648602e-10;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.7224589038431254e-12;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.6134727895776298e-13;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[5], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[5], n,
            ZERO, memslots[7], n);
    /* Deallocating Ba6 in slot 3 */
    /* Deallocating Bb6 in slot 6 */

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.1067710016183396;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[4], 1);

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0002723523751212443;
    /* Smart lincomb recycle T2k6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9924325575605532;
    /* Smart lincomb recycle Bb8_6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[9], 1);

    /* Computing Ba8_2 with operation: lincomb */
    coeff1 = 0.20063211995463168;
    coeff2 = 11.100148209740121;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Ba8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.0262688473369668;
    /* Smart lincomb recycle Ba8_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.39898152652406355;
    /* Smart lincomb recycle Ba8_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.01387684466194758;
    /* Smart lincomb recycle Ba8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.002742528074882574;
    /* Smart lincomb recycle Ba8_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 61.75609557189336;
    /* Smart lincomb recycle Ba8_6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb7_2 with operation: lincomb */
    coeff1 = -0.034743137725474355;
    coeff2 = 0.12434316511921449;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.04622622385713797;
    /* Smart lincomb recycle Bb7_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 12 */

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.006600522551149245;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B3 in slot 7 */

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.000256268554613047;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B4 in slot 2 */

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.892795851939188e-5;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B5 in slot 9 */

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.1438201195458555;
    /* Smart lincomb recycle Bb7_6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B6 in slot 8 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[0], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba7 in slot 5 */
    /* Deallocating Bb7 in slot 1 */

    /* Computing Ba8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.1419333121703887;
    /* Smart lincomb recycle Ba8_7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -7.082137587033525e-5;
    /* Smart lincomb recycle Bb8_7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[9], 1);

    /* Computing B8 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[9], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba8 in slot 3 */
    /* Deallocating Bb8 in slot 10 */

    /* Computing T2k8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.012856181363394585;
    /* Smart lincomb recycle T2k7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B7 in slot 2 */

    /* Computing T2k10 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9918580039315373;
    /* Smart lincomb recycle T2k8 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B8 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

