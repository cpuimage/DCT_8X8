#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static void DCT(float *data, float *outData) {
    float X07P, X16P, X25P, X34P, X34M, X25M, X16M, X07M;
    float X07P34PP, X16P25PP, X16P25PM, X07P34PM;
    float z1, z2, z3, z4, z5, z11, z13;

    int ctr;

    /* Pass 1: process rows. */

    float *In = data;
    float *Out = outData;
    for (ctr = 7; ctr >= 0; ctr--) {
        X07P = In[0] + In[7];
        X07M = In[0] - In[7];
        X16P = In[1] + In[6];
        X16M = In[1] - In[6];
        X25P = In[2] + In[5];
        X25M = In[2] - In[5];
        X34P = In[3] + In[4];
        X34M = In[3] - In[4];

        /* Even part */

        X07P34PP = X07P + X34P;    /* phase 2 */
        X07P34PM = X07P - X34P;
        X16P25PP = X16P + X25P;
        X16P25PM = X16P - X25P;

        Out[0] = X07P34PP + X16P25PP; /* phase 3 */
        Out[4] = X07P34PP - X16P25PP;

        z1 = (X16P25PM + X07P34PM) * ((float) 0.707106781f); /* c4 */
        Out[2] = X07P34PM + z1;    /* phase 5 */
        Out[6] = X07P34PM - z1;

        /* Odd part */

        X07P34PP = X34M + X25M;    /* phase 2 */
        X16P25PP = X25M + X16M;
        X16P25PM = X16M + X07M;

        /* The rotator is modified from fig 4-8 to avoid extra negations. */
        z5 = (X07P34PP - X16P25PM) * ((float) 0.382683433f); /* c6 */
        z2 = ((float) 0.541196100f) * X07P34PP + z5; /* c2-c6 */
        z4 = ((float) 1.306562965f) * X16P25PM + z5; /* c2+c6 */
        z3 = X16P25PP * ((float) 0.707106781f); /* c4 */

        z11 = X07M + z3;        /* phase 5 */
        z13 = X07M - z3;

        Out[5] = z13 + z2;  /* phase 6 */
        Out[3] = z13 - z2;
        Out[1] = z11 + z4;
        Out[7] = z11 - z4;

        In += 8;     /* advance pointer to next row */
        Out += 8;
    }

    /* Pass 2: process columns. */

    In = outData;

    for (ctr = 8 - 1; ctr >= 0; ctr--) {
        X07P = In[8 * 0] + In[8 * 7];
        X07M = In[8 * 0] - In[8 * 7];
        X16P = In[8 * 1] + In[8 * 6];
        X16M = In[8 * 1] - In[8 * 6];
        X25P = In[8 * 2] + In[8 * 5];
        X25M = In[8 * 2] - In[8 * 5];
        X34P = In[8 * 3] + In[8 * 4];
        X34M = In[8 * 3] - In[8 * 4];

        /* Even part */

        X07P34PP = X07P + X34P;    /* phase 2 */
        X07P34PM = X07P - X34P;
        X16P25PP = X16P + X25P;
        X16P25PM = X16P - X25P;

        In[8 * 0] = X07P34PP + X16P25PP; /* phase 3 */
        In[8 * 4] = X07P34PP - X16P25PP;

        z1 = (X16P25PM + X07P34PM) * ((float) 0.707106781f); /* c4 */
        In[8 * 2] = X07P34PM + z1; /* phase 5 */
        In[8 * 6] = X07P34PM - z1;

        /* Odd part */

        X07P34PP = X34M + X25M;    /* phase 2 */
        X16P25PP = X25M + X16M;
        X16P25PM = X16M + X07M;

        /* The rotator is modified from fig 4-8 to avoid extra negations. */
        z5 = (X07P34PP - X16P25PM) * ((float) 0.382683433f); /* c6 */
        z2 = ((float) 0.541196100f) * X07P34PP + z5; /* c2-c6 */
        z4 = ((float) 1.306562965f) * X16P25PM + z5; /* c2+c6 */
        z3 = X16P25PP * ((float) 0.707106781f); /* c4 */

        z11 = X07M + z3;        /* phase 5 */
        z13 = X07M - z3;

        In[8 * 5] = z13 + z2; /* phase 6 */
        In[8 * 3] = z13 - z2;
        In[8 * 1] = z11 + z4;
        In[8 * 7] = z11 - z4;

        In++;          /* advance pointer to next column */
    }
}


static void IDCT(float *data, float *outData) {
    float Y04P26PP, Y04P26MP, Y04P26MM, Y04P26PM, z4, z3, z2, z1, z0;
    float Y04P, Y04M, Y26M, Y26P;
    float z5, Y53M, Y17P, Y17M, Y53P;

    int ctr;

    /* Pass 1: process rows. */

    float *In = data;

    float *Out = outData;
    for (ctr = 7; ctr >= 0; ctr--) {
        /* Even part */

        Y04P = In[0] + In[4];    /* phase 3 */
        Y04M = In[0] - In[4];

        Y26P = In[2] + In[6];    /* phases 5-3 */
        Y26M = (In[2] - In[6]) * 1.414213562f - Y26P; /* 2*c4 */

        Y04P26PP = Y04P + Y26P;    /* phase 2 */
        Y04P26PM = Y04P - Y26P;
        Y04P26MP = Y04M + Y26M;
        Y04P26MM = Y04M - Y26M;

        /* Odd part */

        Y53P = In[5] + In[3];        /* phase 6 */
        Y53M = In[5] - In[3];
        Y17P = In[1] + In[7];
        Y17M = In[1] - In[7];

        z0 = Y17P + Y53P;        /* phase 5 */
        z4 = (Y17P - Y53P) * 1.414213562f; /* 2*c4 */

        z5 = (Y53M + Y17M) * 1.847759065f; /* 2*c2 */
        Y04P = 1.082392200f * Y17M - z5; /* 2*(c2-c6) */
        Y26M = -2.613125930f * Y53M + z5; /* -2*(c2+c6) */

        z1 = Y26M - z0;    /* phase 2 */
        z2 = z4 - z1;
        z3 = Y04P + z2;

        Out[0] = Y04P26PP + z0;
        Out[7] = Y04P26PP - z0;
        Out[1] = Y04P26MP + z1;
        Out[6] = Y04P26MP - z1;
        Out[2] = Y04P26MM + z2;
        Out[5] = Y04P26MM - z2;
        Out[4] = Y04P26PM + z3;
        Out[3] = Y04P26PM - z3;

        In += 8;     /* advance pointer to next row */
        Out += 8;
    }

    /* Pass 2: process columns. */
    float weight = 1.0f / 64.0f;
    In = outData;
    for (ctr = 8 - 1; ctr >= 0; ctr--) {

        Y04P = In[8 * 0] + In[8 * 4];    /* phase 3 */
        Y04M = In[8 * 0] - In[8 * 4];

        Y26P = In[8 * 2] + In[8 * 6];    /* phases 5-3 */
        Y26M = (In[8 * 2] - In[8 * 6]) * 1.414213562f - Y26P; /* 2*c4 */

        Y04P26PP = Y04P + Y26P;    /* phase 2 */
        Y04P26PM = Y04P - Y26P;
        Y04P26MP = Y04M + Y26M;
        Y04P26MM = Y04M - Y26M;

        /* Odd part */

        Y53P = In[8 * 5] + In[8 * 3];        /* phase 6 */
        Y53M = In[8 * 5] - In[8 * 3];
        Y17P = In[8 * 1] + In[8 * 7];
        Y17M = In[8 * 1] - In[8 * 7];

        z0 = Y17P + Y53P;        /* phase 5 */
        Y04M = (Y17P - Y53P) * 1.414213562f; /* 2*c4 */

        z5 = (Y53M + Y17M) * 1.847759065f; /* 2*c2 */
        Y04P = 1.082392200f * Y17M - z5; /* 2*(c2-c6) */
        Y26M = -2.613125930f * Y53M + z5; /* -2*(c2+c6) */

        z1 = Y26M - z0;    /* phase 2 */
        z2 = Y04M - z1;
        z3 = Y04P + z2;

        In[8 * 0] = weight * (Y04P26PP + z0);
        In[8 * 7] = weight * (Y04P26PP - z0);
        In[8 * 1] = weight * (Y04P26MP + z1);
        In[8 * 6] = weight * (Y04P26MP - z1);
        In[8 * 2] = weight * (Y04P26MM + z2);
        In[8 * 5] = weight * (Y04P26MM - z2);
        In[8 * 4] = weight * (Y04P26PM + z3);
        In[8 * 3] = weight * (Y04P26PM - z3);

        In++;          /* advance pointer to next column */
    }
}

int main(int argc, char *argv[]) {
    printf("dct8x8\n");
    printf("DCT implementation by Thomas G. Lane..\n");
    printf("blog:http://cpuimage.cnblogs.com/\n");
    float fft[64] = {0};
    printf("\n input: \n");
    for (int i = 0; i < 64; i++) {
        fft[i] = i;
        printf("%f \t", fft[i]);
    }
    DCT(fft, fft);
    IDCT(fft, fft);
    printf("\n output: \n");
    for (int i = 0; i < 64; i++) {
        printf("%f \t", fft[i]);
    }
    printf("\n press any key to exit.\n");
    getchar();
    return 0;
}
