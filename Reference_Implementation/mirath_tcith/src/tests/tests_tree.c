#include <stdio.h>
#include "../mirath_tree.h"
#include "../common/prng.h"

#define NTESTS 1

int main(int argc, char **argv) {
    // TODO
    mirath_tree_t tree;
    mirath_tree_init(&tree);
//    mirath_tree_print_leaves(&tree);

    // Initialize
    uint8_t salt[2 * MIRATH_SECURITY_BYTES] = {0};
    mirath_prng_t prng;
    mirath_prng_init(&prng, NULL, NULL);
    mirath_prng(&prng, salt, 2 * MIRATH_SECURITY_BYTES);

    for (size_t j = 0; j < MIRATH_SECURITY_BYTES; j++) {
        printf("%02X", salt[j]);
    }
    printf("\n");

    printf("Testing mirath_tree_prg(): ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        mirath_prng(&prng, tree.nodes[0], MIRATH_SECURITY_BYTES);
        tree.nonempty[0] = 1;
        mirath_tree_prg(&tree, salt, 0);
        mirath_tree_print_leaves(&tree);
    }
    printf("\n");

    return 0;
}