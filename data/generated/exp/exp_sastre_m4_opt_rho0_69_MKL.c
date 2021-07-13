#include<mkl/mkl.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

/* Code for matrix function evaluation. */
void dexp_sastre_m4_opt_rho0_69(double *A, const size_t n, double *output) {
    size_t max_memalloc = 7;
    size_t max_memslots = 8;

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
    /* Computation order: Bb2 Ba4_2 Bb4_2 Ba2 B2 Ba4_3 Bb4_3 Ba3_2 Ba3 Ba5_2 Ba5_3 Bb5_2 Bb5_3 T2k2 T2k3 Bb3_2 Bb3 B3 Ba4 Bb4 B4 Ba5_4 Ba5 T2k4 Bb5_4 T2k5 Bb5 B5 T2k7 */

    /* Computing Bb2 with operation: lincomb */
    coeff1 = -0.26484760452432327;
    coeff2 = 0.9903047529247974;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba4_2 with operation: lincomb */
    coeff1 = 0.0002727973329217783;
    coeff2 = 0.05658254301026135;
    memcpy(memslots[2], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[2], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[2], n+1);

    /* Computing Bb4_2 with operation: lincomb */
    coeff1 = 0.0045397184361004065;
    coeff2 = 0.0026543192111867;
    memcpy(memslots[3], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[3], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[3], n+1);

    /* Computing Ba2 with operation: lincomb */
    coeff1 = -0.26484760452432327;
    coeff2 = 0.9903047529247974;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing B2 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[1], n,
            ZERO, memslots[5], n);
    /* Deallocating Ba2 in slot 5 */
    /* Deallocating Bb2 in slot 2 */

    /* Computing Ba4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = -0.31616784713314605;
    /* Smart lincomb recycle Ba4_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb4_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.0004235929351940947;
    /* Smart lincomb recycle Bb4_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing Ba3_2 with operation: lincomb */
    coeff1 = -0.3183560628971706;
    coeff2 = 0.9451759356891595;
    memcpy(memslots[1], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[1], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[1], n+1);

    /* Computing Ba3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.08341693362546203;
    /* Smart lincomb recycle Ba3_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[1], 1);

    /* Computing Ba5_2 with operation: lincomb */
    coeff1 = 0.10536695630497932;
    coeff2 = 1.3909408033494417;
    memcpy(memslots[4], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[4], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[4], n+1);

    /* Computing Ba5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.3494801360488722;
    /* Smart lincomb recycle Ba5_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[4], 1);

    /* Computing Bb5_2 with operation: lincomb */
    coeff1 = 0.14340598505310592;
    coeff2 = 0.18419588627720188;
    memcpy(memslots[6], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[6], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[6], n+1);

    /* Computing Bb5_3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.06331721212120305;
    /* Smart lincomb recycle Bb5_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[6], 1);

    /* Computing T2k2 with operation: lincomb */
    coeff1 = 0.9584219043602061;
    coeff2 = 1.0556647430407748;
    memcpy(memslots[7], memslots[0], n*n*sizeof(*master_mem));
    cblas_dscal(n*n, coeff2, memslots[7], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[7], n+1);

    /* Computing T2k3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.4759594481148073;
    /* Smart lincomb recycle T2k2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb3_2 with operation: lincomb */
    coeff1 = 0.041562405017987775;
    coeff2 = -0.30378066944691157;
    /* Smart lincomb recycle A */
    cblas_dscal(n*n, coeff2, memslots[0], 1);
    cblas_daxpby(n, coeff1, &ONE, 0,
             ONE, memslots[0], n+1);

    /* Computing Bb3 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.949941947908575;
    /* Smart lincomb recycle Bb3_2 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[0], 1);
    /* Deallocating B2 in slot 6 */

    /* Computing B3 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[1], n, memslots[0], n,
            ZERO, memslots[5], n);
    /* Deallocating Ba3 in slot 2 */
    /* Deallocating Bb3 in slot 1 */

    /* Computing Ba4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.9736101476608433;
    /* Smart lincomb recycle Ba4_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[2], 1);

    /* Computing Bb4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 4.468866390298052e-5;
    /* Smart lincomb recycle Bb4_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[3], 1);

    /* Computing B4 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[2], n, memslots[3], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba4 in slot 3 */
    /* Deallocating Bb4 in slot 4 */

    /* Computing Ba5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.029870845474797018;
    /* Smart lincomb recycle Ba5_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[4], 1);

    /* Computing Ba5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.4576696798781078;
    /* Smart lincomb recycle Ba5_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[4], 1);

    /* Computing T2k4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.11560289591544066;
    /* Smart lincomb recycle T2k3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5_4 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.016928981566343947;
    /* Smart lincomb recycle Bb5_3 */
    cblas_daxpby(n*n, coeff2, memslots[5], 1,
             coeff1, memslots[6], 1);
    /* Deallocating B3 in slot 6 */

    /* Computing T2k5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 5.152589478327363;
    /* Smart lincomb recycle T2k4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);

    /* Computing Bb5 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.6928011352478407;
    /* Smart lincomb recycle Bb5_4 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[6], 1);
    /* Deallocating B4 in slot 1 */

    /* Computing B5 with operation: mult */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
            ONE, memslots[4], n, memslots[6], n,
            ZERO, memslots[0], n);
    /* Deallocating Ba5 in slot 5 */
    /* Deallocating Bb5 in slot 7 */

    /* Computing T2k7 with operation: lincomb */
    coeff1 = 1.0;
    coeff2 = 0.7070319120970353;
    /* Smart lincomb recycle T2k5 */
    cblas_daxpby(n*n, coeff2, memslots[0], 1,
             coeff1, memslots[7], 1);
    /* Deallocating B5 in slot 1 */

    /* Prepare output. */
    memcpy(output, memslots[7], n*n*sizeof(*output));
    free(master_mem);
}

