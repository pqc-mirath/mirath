#ifndef MIRATH_TCITH_H
#define MIRATH_TCITH_H

#include <stdint.h>
#include <stdlib.h>
#include "mirath_parameters.h"

typedef uint8_t mirath_tcith_commit_t[2 * MIRATH_SECURITY_BYTES];
typedef mirath_tcith_commit_t mirath_tcith_commit_1_t[MIRATH_PARAM_TAU_1][MIRATH_PARAM_N_1];
typedef mirath_tcith_commit_t mirath_tcith_commit_2_t[MIRATH_PARAM_TAU_2][MIRATH_PARAM_N_2];

void mirath_tcith_commit_set_as_a_grid_list(mirath_tcith_commit_t *seeds[MIRATH_PARAM_TAU],
                                            mirath_tcith_commit_1_t *input_1,
                                            mirath_tcith_commit_2_t *input_2);

void mirath_tcith_commit(mirath_tcith_commit_t commit, const uint8_t *salt, uint8_t e, size_t i, const uint8_t *seed);
size_t mirath_tcith_psi(size_t i, size_t e);

#endif //MIRATH_TCITH_H
