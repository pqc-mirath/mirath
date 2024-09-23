#include <stdint.h>
#include <stdio.h>

#include "random.h"

#ifdef MIRATH_SUPERCOP

/* Nothing to do: 'randombytes' is defined by SUPERCOP. */

#else

#ifdef MIRATH_DETERMINISTIC

#include "prng.h"

static int prng_singleton_up = 0;
static mirath_prng_t prng_singleton;

#else

#ifdef __unix
#include <sys/random.h>
#endif

#ifdef _WIN32
#error "random_bytes on Windows not implemented yet!"
#endif

#endif /* #ifdef MIRATH_DETERMINISTIC */

#endif /* #ifdef MIRATH_SUPERCOP */


#ifndef MIRATH_SUPERCOP

void randombytes(uint8_t *target, size_t length)
{

#ifdef MIRATH_DETERMINISTIC

    if (!prng_singleton_up)
    {
        hash_t salt = {0};

        mirath_prng_init(&prng_singleton, salt, NULL, 0);
        
        prng_singleton_up = 1;
    }
    
    mirath_prng(&prng_singleton, target, length);
    
#else

#if defined (__unix) || defined(__APPLE__)

    FILE* urandom = fopen("/dev/urandom", "r");
    if (urandom == NULL) {
        return;
    }

    if (fread(target, sizeof(uint8_t), length, urandom) != length) {
        return;
    }
    fclose(urandom);

#endif

#if defined(_WIN32) || defined(_WIN64)
#error "randombytes on Windows not implemented yet!"
#endif

#endif /* #ifdef MIRATH_DETERMINISTIC */

}

#endif /* #ifndef MIRATH_SUPERCOP */
