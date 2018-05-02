
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "miniz.h"
#include "dct.h"
#include <stdint.h>

int test_miniz(const unsigned char *s_pStr, uLong data_len) {
    int cmp_status;
    uLong src_len = data_len;
    uLong cmp_len = compressBound(src_len);
    uLong uncomp_len = src_len;
    uint8_t *pCmp, *pUncomp;
    printf("miniz.c version: %s\n", MZ_VERSION);
    // Allocate buffers to hold compressed and uncompressed data.
    pCmp = (mz_uint8 *) malloc((size_t) cmp_len);
    pUncomp = (mz_uint8 *) malloc((size_t) src_len);
    if ((!pCmp) || (!pUncomp)) {
        printf("Out of memory!\n");
        return EXIT_FAILURE;
    }
    // Compress the string.
    cmp_status = compress(pCmp, &cmp_len, (const unsigned char *) s_pStr, src_len);
    if (cmp_status != Z_OK) {
        printf("compress() failed!\n");
        free(pCmp);
        free(pUncomp);
        return EXIT_FAILURE;
    }
    printf("Compressed from %u to %u bytes\n", (mz_uint32) src_len, (mz_uint32) cmp_len);
    // Decompress.
    cmp_status = uncompress(pUncomp, &uncomp_len, pCmp, cmp_len);
    if (cmp_status != Z_OK) {
        printf("uncompress failed!\n");
        free(pCmp);
        free(pUncomp);
        return EXIT_FAILURE;
    }
    printf("Decompressed from %u to %u bytes\n", (mz_uint32) cmp_len, (mz_uint32) uncomp_len);
    // Ensure uncompress() returned the expected data.
    if ((uncomp_len != src_len) || (memcmp(pUncomp, s_pStr, (size_t) src_len))) {
        printf("Decompression failed!\n");
        free(pCmp);
        free(pUncomp);
        return EXIT_FAILURE;
    }
    free(pCmp);
    free(pUncomp);
    printf("Success.\n");
    return EXIT_SUCCESS;

}

int test_dct_miniz(float *data, uLong len) {
    uLong nCount = len / 64;
    float *in_data = data;
    for (int i = 0; i < nCount; i++) {
        DCT(in_data, in_data);
        in_data += 64;
    }
    test_miniz((const unsigned char *) data, len * sizeof(float));
    float *out_data = data;
    for (int i = 0; i < nCount; i++) {
        IDCT(out_data, out_data);
        out_data += 64;
    }
}

int main(int argc, char *argv[]) {
    printf("float data loss compression algorithm base DCT 8X8.\n");
    printf("DCT implementation by Thomas G. Lane.\n");
    printf("miniz implementation by Rich Geldreich.\n");
    //http://developer.download.nvidia.com/SDK/9.5/Samples/vidimaging_samples.html#gpgpu_dct
    printf("blog:http://cpuimage.cnblogs.com/\n");
    int is_debug_output = 1;
    const uLong data_len = 8 * 8;// blocksize
    float test_for_miniz[data_len];
    float test_for_dct[data_len];
    for (int i = 0; i < data_len; ++i) {
        test_for_miniz[i] = rand();
    }
    memcpy(test_for_dct, test_for_miniz, data_len * sizeof(float));
    printf("\nonly miniz:\n");
    test_miniz((const unsigned char *) test_for_miniz, data_len * sizeof(float));
    printf("\nwith dct:\n");
    test_dct_miniz(test_for_dct, data_len);

    if (is_debug_output) {
        for (int i = 0; i < data_len; ++i) {
            if (test_for_miniz[i] != test_for_dct[i]) {
                printf("index %d: %f != %f \n", i, test_for_miniz[i], test_for_dct[i]);
            }
        }
    }

    printf("\n press any key to exit.\n");
    return EXIT_SUCCESS;
}

#ifdef __cplusplus
}
#endif