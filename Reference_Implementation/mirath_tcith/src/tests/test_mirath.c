#include <stdio.h>
#include <string.h>

#include "mirath_sign.h"
#include "api.h"
#include "random.h"
#include "parsing.h"

uint32_t N_TESTS = 100;  /* Number of tests. */
#define MSG_LEN 80   /* Message length. */

int test_mirath_keypair(void) {
    printf("\n\ntest_mirath_keypair...\n");

    for (uint32_t i = 0; i < N_TESTS; i++) {
//        printf("\nTest %d of %d...\n", i + 1, N_TESTS);

        uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};
        uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
        uint8_t pk_computed[CRYPTO_PUBLICKEYBYTES] = {0};

        /* Generate public and secret key. */
        mirath_keypair(pk, sk);

        ff_t S[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)];
        ff_t C[mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R)];
        ff_t H[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K)];
        ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];
        ff_t y_computed[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)];

        seed_t seed_pk;

        parse_public_key(seed_pk, y, pk);
        parse_secret_key(S, C, H, sk);

        mirath_tcith_internal_steps_pk(y_computed, S, C, H);

        if (memcmp(y, y_computed, mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)) != 0) {
            printf("y != y_computed\n");
            return -1;
        }

        unparse_public_key(pk_computed, seed_pk, y_computed);

        if (memcmp(pk, pk_computed, CRYPTO_PUBLICKEYBYTES) != 0) {
            printf("pk != pk_computed\n");
            return -1;
        }
    }

    printf("Everything's OK!\n");

    return 0;
}

int test_mirath(void) {
    printf("\n\ntest_mirath...\n");

    for (uint32_t i = 0; i < N_TESTS; i++)
    {
        uint8_t msg[MSG_LEN] = {0};
//        uint8_t msg2[MSG_LEN] = {0};
        uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};
        uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
        uint8_t sig_msg[CRYPTO_BYTES + MSG_LEN];
        // Todo: add message in the signature
//        size_t sig_msg_len = MIRATH_SIGNATURE_BYTES + MSG_LEN;
        size_t sig_msg_len = MIRATH_SIGNATURE_BYTES;
        size_t msg_len = MSG_LEN;
        size_t msg2_len = MSG_LEN;
        size_t pos;
        uint8_t byte;

//        printf("\nTest %d of %d...\n\n", i + 1, N_TESTS);
        printf("Testing [%u]:\t", i);
        printf("%2u%%", 100 * i / N_TESTS);
//        fflush(stdout);
//        printf("\r\x1b[K");

        /* Generate a random message. */
        randombytes(msg, MSG_LEN);

        /* Generate public and secret key. */
        mirath_keypair(pk, sk);

        /* Sign the message.*/
        //mirath_sign(uint8_t *sig_msg, uint8_t *msg, size_t msg_len, uint8_t *sk)
        mirath_sign(sig_msg, msg, MSG_LEN, sk);

        /* Verify the message */
        // Todo: add message in the signature
//        if (mirath_verify(msg2, &msg2_len, sig_msg, pk) != 0)
        if (mirath_verify(msg, &msg2_len, sig_msg, pk) != 0)
        {
            printf("\nError: Verification failed!\n");
            return -1;
        }

        /* Check the message length. */
        if (msg_len != msg2_len)
        {
            printf("\nError: Message lengths don't match!\n");
            return -1;
        }

        // Todo: add message in the signature
        /* Check the message. */
//        if (memcmp(msg, msg2, msg_len) != 0)
//        {
//            printf("Error: Messages don't match!\n");
//            return -1;
//        }

        /* Change one random byte of the signature. */
        randombytes((uint8_t *)&pos, sizeof(pos));
        pos %= sig_msg_len;

        do
        {
            randombytes(&byte, sizeof(byte));
        }
        while (byte == 0);

        sig_msg[pos] ^= byte;
        /* * */

        /* Verify the forged signature. */
        if (mirath_verify(msg, &msg_len, sig_msg, pk) == 0)
        {
            printf("Error: Trivial forgery possible!\n");
            return -1;
        }

//        printf("OK!\n");
        fflush(stdout);
        printf("\r\x1b[K");
    }

    printf("Everything's OK!\n\n");

    return 0;
}

int main(void) {
    int ret;
    printf("\nBytes concerning secret key:\t%d", MIRATH_SECRET_KEY_BYTES);
    printf("\nBytes concerning public key:\t%d", MIRATH_PUBLIC_KEY_BYTES);
    printf("\nBytes concerning signature: \t%d\n\n", MIRATH_SIGNATURE_BYTES);
    ret = test_mirath_keypair();
    ret ^= test_mirath();

    return (ret >> 31);
}