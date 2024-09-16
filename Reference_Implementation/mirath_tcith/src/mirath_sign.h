#ifndef MIRATH_TCITH_MIRATH_SIGN_H
#define MIRATH_TCITH_MIRATH_SIGN_H

#include <stddef.h>
#include <stdint.h>

int mirath_keypair(uint8_t *pk, uint8_t *sk);

int mirath_sign(uint8_t *sig_msg, uint8_t *msg, size_t msg_len, uint8_t *sk);

int mirath_verify(uint8_t *msg, size_t *msg_len, uint8_t *sig_msg, uint8_t *pk);

#endif //MIRATH_TCITH_MIRATH_SIGN_H
