#include <stdio.h>
#include "mirath_tree.h"
#include "prng.h"

#define NTESTS 1

int main(int argc, char **argv) {
    mirath_tree_t tree;
    mirath_tree_init(&tree);

    // Initialize
    uint8_t salt[2 * MIRATH_SECURITY_BYTES] = {0};
    mirath_prng_t prng;
    mirath_prng_init(&prng, NULL, NULL, 0);
    mirath_prng(&prng, salt, 2 * MIRATH_SECURITY_BYTES);

    printf("Testing mirath_tree: ");
    for (uint32_t i = 0; i < NTESTS; i++) {
        mirath_tree_t partial_tree;
        mirath_tree_init(&partial_tree);

        mirath_prng(&prng, tree.nodes[0], MIRATH_SECURITY_BYTES);
        tree.nonempty[0] = 1;
        mirath_tree_prg(&tree, salt, 0);
////        mirath_tree_print_leaves(&tree);
////        mirath_tree_print(&tree);

        uint8_t revealed[MIRATH_PARAM_TREE_LEAVES] = {0}; // We assume the worst case scenario
        size_t revealed_length;
        const size_t challenge_length = MIRATH_PARAM_TAU;
        uint16_t challenge_list[challenge_length] = {0};
        for (size_t k = 0; k < challenge_length; k++) {
            mirath_prng(&prng, &challenge_list[k], sizeof(uint16_t));
            challenge_list[k] = challenge_list[k] % MIRATH_PARAM_TREE_LEAVES;
        }

        revealed_length = mirath_tree_get_sibling_path(revealed, &tree, challenge_list, challenge_length);

//        printf("revealed length: %ld\n", revealed_length);

//        mirath_tree_print_sibling_path(revealed, revealed_length);

        mirath_tree_get_seeds_from_path(&partial_tree, challenge_list, challenge_length, revealed, revealed_length, salt, 0);
//        mirath_tree_print_leaves(&partial_tree);

        uint32_t first_leaf = MIRATH_PARAM_TREE_NODES - MIRATH_PARAM_TREE_LEAVES;
        const uint8_t zero[MIRATH_SECURITY_BYTES] = {0};
        for (uint32_t j = first_leaf; j < MIRATH_PARAM_TREE_NODES; j++) {
            for (uint32_t k = 0; k < challenge_length; k++) {
                if (j == challenge_list[k]) {
                    if ((memcmp(partial_tree.nodes[j], zero, MIRATH_SECURITY_BYTES) != 0)) {
                        printf("ERROR: partial tree is not zero\n");
                        return -1;
                    }
                }
            }

            if ((memcmp(partial_tree.nodes[j], zero, MIRATH_SECURITY_BYTES) != 0) && memcmp(tree.nodes[j], partial_tree.nodes[j], MIRATH_SECURITY_BYTES) != 0) {
                printf("j: %d\n", j);
                printf("ERROR: trees are not equal\n");
                return -1;
            }
        }

//        for (size_t k = 0; k < challenge_length; k++) {
//            printf("%04u ", challenge_list[k]);
//        }
    }
    printf("OK!\n");

    return 0;
}