#include<cblas.h>
#include<lapacke.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_ps_opt_rho2_7(double *A, const size_t n, double *output) {
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
    coeff1 = -0.8558107709451741;
    coeff2 = 6.654420725242113;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba6_2 with operation: lincomb */
    coeff1 = 0.0035871003621533335;
    coeff2 = -0.199515898070094;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = -4.433844812389262;
    coeff2 = 0.14174985095473627;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Bb2 with operation: lincomb */
    coeff1 = 2.373110447113857;
    coeff2 = -11.496160069358545;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = -0.7086288697769633;
    coeff2 = 1.5171963077831512;
    memcpy(memslots[5], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[5], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[5], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = -0.7929434423854941;
    coeff2 = 5.6454298259935385;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = 2.373110447113857;
    coeff2 = -11.496160069358545;
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
    coeff2 = -8.805953709057599;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[5], 1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.1464046050697995;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.3717237028855647;
    /* Smart lincomb recycle Ba6_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[2], 1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 14.35374095270583;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 16.186403987022313;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[6], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = 4.521013785174327e-9;
    coeff2 = -0.016780854076753563;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.12296702237640218;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 0.01201228900944982;
    coeff2 = -1.1065314556032247e-12;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -2.0779221004950836e-15;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[7], 1);

    /* Computing Ba7_2 with operation: lincomb */
    coeff1 = -0.09860921274946392;
    coeff2 = 0.2639396617162851;
    memcpy(memslots[9], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[9], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[9], n+1);

    /* Computing Ba7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -4.745046748609644;
    /* Smart lincomb recycle Ba7_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb6_2 with operation: lincomb */
    coeff1 = -0.13628986461873388;
    coeff2 = -0.027077600537154935;
    memcpy(memslots[10], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[10], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[10], n+1);

    /* Computing Bb6_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.004774676777784938;
    /* Smart lincomb recycle Bb6_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[10], 1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = 2.6186718329845604;
    coeff2 = 4.547282965965041;
    memcpy(memslots[11], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[11], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[11], n+1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 5.695287571427549;
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
    coeff2 = 25.592199199005485;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[5], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.6370814945587178;
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
    coeff2 = 0.012069483753455755;
    /* Smart lincomb recycle Bb6_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[10], 1);

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.1554804561626417;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.914091240339738;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.009542640674472407;
    /* Smart lincomb recycle Bb6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[10], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -16.423040024910904;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[3], 1);

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 1.6798059142303008;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -5.249610150102727e-19;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.267589684021346e-27;
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
    coeff2 = 3.7664383775767242e-6;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8690669836434851;
    /* Smart lincomb recycle Bb6_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[10], 1);

    /* Computing Ba6_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -1.0770978480112308;
    /* Smart lincomb recycle Ba6_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 2.1074592605382616;
    /* Smart lincomb recycle Ba6_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[2], 1);

    /* Computing Ba6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.5282957524727887;
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
    coeff2 = 0.04873878019200554;
    /* Smart lincomb recycle T2k6 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[3], 1);

    /* Computing Bb7_2 with operation: lincomb */
    coeff1 = 5.692302453618037;
    coeff2 = -2.4680601957893162;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb7_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.7182801012293478;
    /* Smart lincomb recycle Bb7_2 */
    cblas_daxpby(n*n, coeff2, memslots[8], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 9 */

    /* Computing Bb7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.5030745298888237;
    /* Smart lincomb recycle Bb7_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.145388920252181;
    /* Smart lincomb recycle Ba7_3 */
    cblas_daxpby(n*n, coeff2, memslots[12], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B3 in slot 13 */

    /* Computing Bb7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.2802466342872092;
    /* Smart lincomb recycle Bb7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 3.6815285908979076;
    /* Smart lincomb recycle Ba7_4 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B4 in slot 2 */

    /* Computing Bb7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -7.214590944431316e-5;
    /* Smart lincomb recycle Bb7_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[0], 1);

    /* Computing Ba7_6 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.00016028279238559438;
    /* Smart lincomb recycle Ba7_5 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[9], 1);
    /* Deallocating B5 in slot 6 */

    /* Computing Ba7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -9.362127616405266e-5;
    /* Smart lincomb recycle Ba7_6 */
    cblas_daxpby(n*n, coeff2, memslots[4], 1,
             coeff1, memslots[9], 1);

    /* Computing Bb7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.8753927048889629;
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
    coeff2 = 1.3559720150054728;
    /* Smart lincomb recycle T2k7 */
    cblas_daxpby(n*n, coeff2, memslots[1], 1,
             coeff1, memslots[3], 1);
    /* Deallocating B7 in slot 2 */

    /* Prepare output. */
    memcpy(output, memslots[3], n*n*sizeof(*output));
    free(master_mem);
}

