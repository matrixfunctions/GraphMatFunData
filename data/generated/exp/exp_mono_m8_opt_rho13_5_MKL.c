#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_mono_m8_opt_rho13_5(double *A, const size_t n, double *output) {
    size_t max_memalloc = 16;
    size_t max_memslots = 17;

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
    /* Computation order: Ba3_2 Ba6_2 T2k2 Ba7_2 Bb3_2 Ba4_2 Bb9_2 Bb4_2 Ba9_2 Bb8_2 Ba5_2 Bb5_2 Bb6_2 Ba8_2 Bb7_2 B2 Bb3 Bb7_3 Bb6_3 Bb8_3 Ba4_3 Ba3 B3 Bb7_4 Ba5_3 Ba5_4 Bb5_3 Ba7_3 Ba9_3 Ba7_4 Ba8_3 Ba4 Bb6_4 Ba6_3 Ba9_4 T2k3 Bb4_3 Bb9_3 Bb4 B4 T2k4 Bb7_5 Bb9_4 Bb9_5 Ba5 Ba9_5 T2k5 Bb6_5 Ba8_4 Ba8_5 Bb5_4 Bb5 B5 T2k6 Bb9_6 Ba8_6 Ba9_6 Bb7_6 Bb6 Ba6_4 Bb8_4 Bb8_5 Ba6_5 Ba7_5 Ba6 B6 Bb9_7 Ba8_7 Ba9_7 Bb7 T2k7 Ba7_6 Bb8_6 Ba7 Bb8_7 B7 Ba9_8 Bb9_8 Ba8 Bb8 B8 Ba9 T2k9 Bb9 B9 T2k11 */

    /* Computing Ba3_2 with operation: lincomb */
    coeff1 = -0.003584018414055765;
    coeff2 = 0.16376260367565243;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = 4.437908455542293e-5;
    coeff2 = 0.0003258145596094315;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 0.008776543046053789;
    coeff2 = -1.3078408314728238e-5;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Ba7_2 with operation: lincomb */
    coeff1 = -2.8201037195155804e-6;
    coeff2 = 4.8728643050646066e-6;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = 0.30528232877254763;
    coeff2 = 0.23746768414771152;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = 0.0005515398938947843;
    coeff2 = 0.02533967378579724;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Bb9_2 with operation: lincomb */
    coeff1 = -0.06540541282818779;
    coeff2 = 6.846710342127941e-5;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = -0.008872840493334803;
    coeff2 = 0.25571769129339744;
    memcpy(memslots[8], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[8], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[8], n+1);

    /* Computing Ba9_2 with operation: lincomb */
    coeff1 = 0.1346603332119568;
    coeff2 = -5.9198518169123926e-5;
    memcpy(memslots[9], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[9], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[9], n+1);

    /* Computing Bb8_2 with operation: lincomb */
    coeff1 = 1.0002604649982476;
    coeff2 = 0.2422828857323093;
    memcpy(memslots[10], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[10], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[10], n+1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = -0.0012630062267903097;
    coeff2 = -0.0010962696861478022;
    memcpy(memslots[11], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[11], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[11], n+1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 0.8293727566445452;
    coeff2 = 0.1675033881979259;
    memcpy(memslots[12], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[12], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[12], n+1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -0.10275230374603168;
    coeff2 = 0.25648676493721656;
    memcpy(memslots[13], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[13], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[13], n+1);

    /* Computing Ba8_2 with operation: lincomb */
    coeff1 = 1.0002604649982476;
    coeff2 = 0.2422828857323093;
    memcpy(memslots[14], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[14], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[14], n+1);

    /* Computing Bb7_2 with operation: lincomb */
    coeff1 = 0.004650329563627054;
    coeff2 = 0.2598813154789629;
    memcpy(memslots[15], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[15], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[15], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[0], n, memslots[0], n,
            ZERO, memslots[16], n);
    /* Deallocating A in slot 1 */

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0023273769915510135;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.006290014251994485;
    /* Smart lincomb recycle Bb7_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[15], 1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.005628536962635935;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[13], 1);

    /* Computing Bb8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.022001533431922154;
    /* Smart lincomb recycle Bb8_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[10], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.017067262331527198;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[6], 1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.04917059370539749;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[1], 1);

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[5], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba3 in slot 2 */
    /* Deallocating Bb3 in slot 6 */

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.010576806709956765;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[15], 1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.00022839563887310042;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[11], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.013673005502563394;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[11], 1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0030839993650848517;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[12], 1);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.2723646529748017e-6;
    /* Smart lincomb recycle Ba7_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba9_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -8.660244928553117e-6;
    /* Smart lincomb recycle Ba9_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[9], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0061294802470251075;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.022001533431922154;
    /* Smart lincomb recycle Ba8_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[14], 1);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9638185168573307;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[6], 1);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.008329230127484285;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[13], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.00022081927999800417;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba9_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -8.731780915688725e-5;
    /* Smart lincomb recycle Ba9_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -9.364829337290504e-7;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0037256194747560907;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[8], 1);

    /* Computing Bb9_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.738382227960167e-6;
    /* Smart lincomb recycle Bb9_2 */
    cblas_daxpby(n*n, coeff2, memslots[16], 1,
             coeff1, memslots[7], 1);
    /* Deallocating B2 in slot 17 */

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.002052240135641998;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[8], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[6], n, memslots[8], n,
            ZERO, memslots[1], n);
    /* Deallocating Ba4 in slot 7 */
    /* Deallocating Bb4 in slot 9 */

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -4.881187067247219e-6;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0003949349245182783;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[15], 1);

    /* Computing Bb9_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -6.3730674335546956e-6;
    /* Smart lincomb recycle Bb9_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb9_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -5.739005263969626e-6;
    /* Smart lincomb recycle Bb9_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8216219421030205;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[11], 1);

    /* Computing Ba9_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -2.1711606510903095e-5;
    /* Smart lincomb recycle Ba9_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -6.561388347938557e-7;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.00021583757875408033;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[13], 1);

    /* Computing Ba8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.15933029209747282;
    /* Smart lincomb recycle Ba8_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[14], 1);

    /* Computing Ba8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.03184541272812379;
    /* Smart lincomb recycle Ba8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[14], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0031897875289976277;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[12], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.234759791150232e-5;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[12], 1);

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[11], n, memslots[12], n,
            ZERO, memslots[5], n);
    /* Deallocating Ba5 in slot 12 */
    /* Deallocating Bb5 in slot 13 */

    /* Computing T2k6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.0279110807751132e-7;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb9_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -2.157568588204766e-6;
    /* Smart lincomb recycle Bb9_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.007635550247315368;
    /* Smart lincomb recycle Ba8_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[14], 1);

    /* Computing Ba9_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -5.9723490061415135e-6;
    /* Smart lincomb recycle Ba9_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 8.856309878298438e-6;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[15], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.0968604718569746e-7;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[13], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.13039438638797307;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.15933029209747282;
    /* Smart lincomb recycle Bb8_3 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[10], 1);
    /* Deallocating B3 in slot 1 */

    /* Computing Bb8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.03184541272812379;
    /* Smart lincomb recycle Bb8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[10], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8475672106088131;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0688513777099278;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[4], 1);
    /* Deallocating B4 in slot 2 */

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8145974867536132;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[2], 1);

    /* Computing B6 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[13], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba6 in slot 3 */
    /* Deallocating Bb6 in slot 14 */

    /* Computing Bb9_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -3.3636703984819904e-7;
    /* Smart lincomb recycle Bb9_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0008615033049494069;
    /* Smart lincomb recycle Ba8_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[14], 1);

    /* Computing Ba9_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -7.70252692547925e-7;
    /* Smart lincomb recycle Ba9_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 7.553245375703866e-9;
    /* Smart lincomb recycle Bb7_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[15], 1);

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -5.458270135394687e-9;
    /* Smart lincomb recycle T2k6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.12038470670165373;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.007635550247315368;
    /* Smart lincomb recycle Bb8_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[10], 1);
    /* Deallocating B5 in slot 6 */

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9948034699390058;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0008615033049494069;
    /* Smart lincomb recycle Bb8_6 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[10], 1);
    /* Deallocating B6 in slot 1 */

    /* Computing B7 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[15], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba7 in slot 5 */
    /* Deallocating Bb7 in slot 16 */

    /* Computing Ba9_8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -3.1855167508881976e-8;
    /* Smart lincomb recycle Ba9_7 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb9_8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.5592924041261544e-8;
    /* Smart lincomb recycle Bb9_7 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.134369854398072e-5;
    /* Smart lincomb recycle Ba8_7 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[14], 1);

    /* Computing Bb8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.134369854398072e-5;
    /* Smart lincomb recycle Bb8_7 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[10], 1);
    /* Deallocating B7 in slot 1 */

    /* Computing B8 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[14], n, memslots[10], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba8 in slot 15 */
    /* Deallocating Bb8 in slot 11 */

    /* Computing Ba9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0081137407871135;
    /* Smart lincomb recycle Ba9_8 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.06726697857924235;
    /* Smart lincomb recycle T2k7 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb9 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9946441317159572;
    /* Smart lincomb recycle Bb9_8 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);
    /* Deallocating B8 in slot 1 */

    /* Computing B9 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[9], n, memslots[7], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba9 in slot 10 */
    /* Deallocating Bb9 in slot 8 */

    /* Computing T2k11 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9964850858948753;
    /* Smart lincomb recycle T2k9 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B9 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

