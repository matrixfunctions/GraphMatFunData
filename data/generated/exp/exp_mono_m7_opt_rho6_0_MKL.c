#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_mono_m7_opt_rho6_0(double *A, const size_t n, double *output) {
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
    coeff1 = -0.03985250441225717;
    coeff2 = 0.1604494657928954;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = -0.00011974681141880265;
    coeff2 = -9.609608204465964e-5;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = -0.00018687382272586462;
    coeff2 = -0.00010155826657765157;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Ba7_2 with operation: lincomb */
    coeff1 = -4.539313763978297e-6;
    coeff2 = 2.446796588580378e-6;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = 0.3044412670842123;
    coeff2 = 0.4680830251883016;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Bb2 with operation: lincomb */
    coeff1 = 0.19853207576158147;
    coeff2 = 0.4568156227748025;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = -0.01206538482409284;
    coeff2 = -0.007524875465546422;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = -0.011223283513386698;
    coeff2 = 0.5004951105037614;
    memcpy(memslots[8], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[8], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[8], n+1);

    /* Computing Bb8_2 with operation: lincomb */
    coeff1 = 0.9834849370950302;
    coeff2 = 0.40790543324377204;
    memcpy(memslots[9], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[9], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[9], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = 0.19853207576158147;
    coeff2 = 0.4568156227748025;
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
    coeff2 = 0.029166175784732254;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.4202712650569836;
    /* Smart lincomb recycle Bb8_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[9], 1);

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.32310560219905493;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9304976822330905;
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
    coeff2 = 4.215986081751763e-5;
    /* Smart lincomb recycle Ba7_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0061270230013822065;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9639420128969182;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.004160420231960609;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.15990304681757225;
    /* Smart lincomb recycle Bb8_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.0001809250904584009;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.00011452178091543187;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0587301450366299;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[8], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0020530677449521447;
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
    coeff2 = 0.033732264254459236;
    /* Smart lincomb recycle Bb8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[9], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -3.544767352137871e-5;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.13025928619739066;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8471750409761895;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.06884543294512852;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = -0.0014350882470920096;
    coeff2 = -0.0029204003557249366;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.004322247579028231;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.013953282030764238;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8220608344548307;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 0.8268548885673269;
    coeff2 = 0.32692355515774696;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.05385438234818326;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.003348134917199928;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.728916919898747e-5;
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
    coeff2 = -1.1868405154805229e-5;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8150363087806349;
    /* Smart lincomb recycle Ba6_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.1203307382956932;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.008536371457893934;
    /* Smart lincomb recycle Bb8_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -0.10687141220929966;
    coeff2 = 0.49573353466641706;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.09669233528231447;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.008584243736221898;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0002579652548561839;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -6.616194274488254e-8;
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
    coeff2 = 0.9948104309724103;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[4], 1);

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.9506111159992416e-6;
    /* Smart lincomb recycle T2k6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0010322498247622353;
    /* Smart lincomb recycle Bb8_6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[9], 1);

    /* Computing Ba8_2 with operation: lincomb */
    coeff1 = 0.9834849370950302;
    coeff2 = 0.40790543324377204;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Ba8_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.4202712650569836;
    /* Smart lincomb recycle Ba8_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.15990304681757225;
    /* Smart lincomb recycle Ba8_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.033732264254459236;
    /* Smart lincomb recycle Ba8_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.008536371457893934;
    /* Smart lincomb recycle Ba8_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba8_7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0010322498247622353;
    /* Smart lincomb recycle Ba8_6 */
    cblas_daxpby(n*n, coeff2, memslots[7], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb7_2 with operation: lincomb */
    coeff1 = -8.127184169897754e-8;
    coeff2 = 0.5009422778394859;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.10725741908100991;
    /* Smart lincomb recycle Bb7_2 */
    cblas_daxpby(n*n, coeff2, memslots[11], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 12 */

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.011203185515508135;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[6], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B3 in slot 7 */

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0005125194848333935;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B4 in slot 2 */

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.4968901291255363e-5;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B5 in slot 9 */

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.941362084939232e-8;
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
    coeff2 = 4.413245234287325e-5;
    /* Smart lincomb recycle Ba8_7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb8 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.413245234287325e-5;
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
    coeff2 = -1.0864858600469639e-7;
    /* Smart lincomb recycle T2k7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B7 in slot 2 */

    /* Computing T2k10 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.0004223592531323;
    /* Smart lincomb recycle T2k8 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B8 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

