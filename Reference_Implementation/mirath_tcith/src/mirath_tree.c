/**
 * @file mirath_tree.c
 * @brief Implementation of tree related functions
 */

#include <stdio.h>
#include <string.h>
#include "mirath_tree.h"
#include "common/hash.h"


/**
 * \fn int mirath_tree_node_contains(size_t* list, size_t len, size_t value)
 * \brief This function determines if value belongs to a list
 *
 * \param[in] list size_t* list of integers
 * \param[in] len size_t length of the input list
 * \param[in] value size_t element to determine if belongs to the given list
 */
int mirath_tree_node_contains(size_t* list, size_t len, size_t value) {
    for (size_t i = 0; i < len; i++) {
        if (list[i] == value) {
            return 1;
        }
    }
    return 0;
}


/**
 * \fn int mirath_tree_node_exists(const mirath_tree_t* tree, size_t i)
 * \brief This function determines if the i-th node exists in the tree
 *
 * \param[in] tree mirath_tree_t* Representation of the tree
 * \param[in] i size_t integer node number
 */
int mirath_tree_node_exists(const mirath_tree_t* tree, size_t i) {
    if (i >= MIRATH_PARAM_TREE_NODES) {
        return 0;
    }
    if (tree->exists[i]) {
        return 1;
    }
    return 0;
}


/**
 * \fn void mirath_tree_init(mirath_tree_t *tree)
 * \brief This function initializes to zero each node of the tree, and sets which nodes exists (resp. does not exist)
 *
 * \param[in,out] tree mirath_tree_t* Representation of the tree
 */
void mirath_tree_init(mirath_tree_t *tree) {
    memset(tree->nodes, 0, MIRATH_PARAM_TREE_NODES * sizeof(mirath_tree_node_t));
    memset(tree->nonempty, 0, MIRATH_PARAM_TREE_NODES);
    memset(tree->exists, 0, MIRATH_PARAM_TREE_NODES);

    memset(&tree->exists[MIRATH_PARAM_TREE_NODES - MIRATH_PARAM_TREE_LEAVES], 1, MIRATH_PARAM_TREE_LEAVES); // Set leaves
    for (int i = MIRATH_PARAM_TREE_NODES - MIRATH_PARAM_TREE_LEAVES; i > 0; i--) {
        if (mirath_tree_node_exists(tree, 2 * i + 1) || mirath_tree_node_exists(tree, 2 * i + 2) ) {
            tree->exists[i] = 1;
        }
    }
    tree->exists[0] = 1;
}


/**
 * \fn void mirath_tree_init(mirath_tree_t *tree)
 * \brief This function makes zero each field of the tree structure
 *
 * \param[in,out] tree mirath_tree_t* Representation of the tree
 */
void mirath_tree_clear(mirath_tree_t *tree) {
    memset(tree->nodes, 0, MIRATH_PARAM_TREE_NODES * sizeof(mirath_tree_node_t));
    memset(tree->nonempty, 0, MIRATH_PARAM_TREE_NODES);
    memset(tree->exists, 0, MIRATH_PARAM_TREE_NODES);
}


/**
 * \fn void mirath_tree_prg(mirath_tree_t *tree, const uint8_t *salt, uint8_t e)
 * \brief This function computes a tree by expanding its master seed.
 *
 * \param[out] tree mirath_tree_t* Representation of the tree with master seed at position tree[0]
 * \param[in] salt uint8_t* Salt used for the signature
 * \param[in] e uint8_t Current iteration
 */
void mirath_tree_prg(mirath_tree_t *tree, const uint8_t *salt, uint8_t e) {
    uint8_t domain_separator = DOMAIN_SEPARATOR_TREE;
    uint8_t digest[2 * MIRATH_SECURITY_BYTES];

    // We only need to calculate children of non-leaf nodes
    size_t non_leaf_nodes = mirath_tree_get_parent(MIRATH_PARAM_TREE_NODES - 1);

    // Next, we proceed by generating the children
    for(size_t i = 0; i <= non_leaf_nodes; i++) {
        if (!tree->nonempty[i]) {
            continue;
        }
        size_t from = i;
        size_t to = i * 2 + 1;

        hash_tree_digest(digest, &domain_separator, salt, &e, tree->nodes[from]);

        if (!tree->nonempty[to]) {
            // left child
            memcpy(tree->nodes[to], digest, MIRATH_SECURITY_BYTES);
            tree->nonempty[to] = 1;
        }

        // The last non-leaf node will only have a left child when there are an odd number of leaves
        if (mirath_tree_node_exists(tree, to + 1) && !tree->nonempty[to + 1]) {
            // right child
            memcpy(tree->nodes[to + 1], &digest[MIRATH_SECURITY_BYTES], MIRATH_SECURITY_BYTES);
            tree->nonempty[to + 1] = 1;
        }
    }
}

/**
 * \fn int mirath_tree_has_right_child(const mirath_tree_t *tree, size_t node)
 * \brief This function determines if the tree has a right child at the input node position.
 *
 * \param[in] tree mirath_tree_t Representation of the tree
 * \param[in] node size_t Node position in the tree
 */
int mirath_tree_has_right_child(const mirath_tree_t *tree, size_t node) {
    return(2 * node + 2 < MIRATH_PARAM_TREE_NODES && mirath_tree_node_exists(tree, node));
}


/**
 * \fn int mirath_tree_is_a_leaf(size_t node)
 * \brief This function determines if the input node position describes a leaves
 *
 * \param[in] node size_t Node position in the tree
 */
int mirath_tree_is_a_leaf(size_t node) {
    return (2 * node + 1 >= MIRATH_PARAM_TREE_NODES);
}

/**
 * \fn size_t mirath_tree_get_parent(size_t node)
 * \brief This function returns the parent position of the input node position
 *
 * \param[in] node size_t Node position in the tree
 */
size_t mirath_tree_get_parent(size_t node) {
    return (node - 2 + (node % 2)) / 2;
}


/**
 * \fn int mirath_tree_has_sibling(const mirath_tree_t* tree, size_t node)
 * \brief This function determines if the tree has a sibling at the input node position.
 *
 * \param[in] tree mirath_tree_t Representation of the tree
 * \param[in] node size_t Node position in the tree
 */
int mirath_tree_has_sibling(const mirath_tree_t* tree, size_t node) {
    if (!mirath_tree_node_exists(tree, node)) {
        return 0;
    }

    if ((node % 2) == 1 && !mirath_tree_node_exists(tree, node + 1)) {
        return 0;
    }

    return 1;
}


/**
 * \fn size_t mirath_tree_get_sibling(size_t node)
 * \brief This function returns the sibling position of the input node position
 *
 * \param[in] node size_t Node position in the tree
 */
size_t mirath_tree_get_sibling(size_t node) {
    if ((node % 2) == 1) {
        if (node + 1 < MIRATH_PARAM_TREE_NODES) {
            return node + 1;
        }
        else {
            // This node does not have siblings");
            return 0;
        }
    }
    else {
        return node - 1;
    }
}


/**
 * \fn void mirath_tree_get_revealed_nodes(size_t *output_revealed, size_t* output_length,
 *                                       const mirath_tree_t *tree,
 *                                       const uint16_t *challenge_list, size_t challenge_length)
 * \brief This function calculates the node positions for a path not revealing the challenges
 *
 * \param[out] output_revealed uint8_t* Node positions of the path that hides the challenges
 * \param[out] output_length size_t* Length of path
 * \param[in] tree mirath_tree_t* Representation of the tree with master seed at position tree[0]
 * \param[in] challenge_list uint16_t* List of challenges (node positions to be hidden)
 * \param[in] challenge_length sizte_t Number of challenges
 */
void mirath_tree_get_revealed_nodes(size_t *output_revealed, size_t* output_length,
                                  const mirath_tree_t *tree,
                                  const uint16_t *challenge_list, const size_t challenge_length) {
    size_t path_length = MIRATH_PARAM_TREE_DEPTH;
    size_t path_sets[path_length][challenge_length];

    // Compute the paths back to the root
    for (size_t i = 0; i < challenge_length; i++) {
        size_t pos = 0;
        size_t node = challenge_list[i] + (MIRATH_PARAM_TREE_NODES - MIRATH_PARAM_TREE_LEAVES);
        path_sets[pos][i] = node;
        pos++;
        while ( (node = mirath_tree_get_parent(node)) != 0 ) {
            path_sets[pos][i] = node;
            pos++;
        }
    }

    /* Determine seeds to reveal */
    size_t position = 0;
    for (size_t d = 0; d < path_length; d++) {
        for (size_t i = 0; i < challenge_length; i++) {
            if (!mirath_tree_has_sibling(tree, path_sets[d][i])) {
                continue;
            }

            size_t sibling = mirath_tree_get_sibling(path_sets[d][i]);
            if (!mirath_tree_node_contains(path_sets[d], challenge_length, sibling )) {
                // Determine the seed to reveal
                while(!mirath_tree_has_right_child(tree, sibling) && !mirath_tree_is_a_leaf(sibling)) {
                    sibling = 2 * sibling + 1;
                }

                // Only reveal if we haven't already
                if (!mirath_tree_node_contains(output_revealed, position, sibling)) {
                    output_revealed[position] = sibling;
                    position++;
                }
            }
        }
    }
    *output_length = position;
}

/**
 * \fn size_t mirath_tree_get_sibling_path(uint8_t* output,
 *                                       const mirath_tree_t* tree,
 *                                       uint16_t *challenge_list, size_t challenge_length)
 * \brief This function calculates the sibling path that hides the challenges, and returns the path length.
 *
 * \param[out] output uint8_t* Sibling path that hides the challenges
 * \param[in] tree mirath_tree_t* Representation of the tree with master seed at position tree[0]
 * \param[in] challenge_list uint16_t* List of challenges (node positions to be hidden)
 * \param[in] challenge_length sizte_t Number of challenges
 */
size_t mirath_tree_get_sibling_path(uint8_t* output,
                                  const mirath_tree_t* tree,
                                  uint16_t *challenge_list, size_t challenge_length) {
    size_t revealed[MIRATH_PARAM_TREE_LEAVES] = {0}; // We assume the worst case scenario
    size_t revealed_length;

    mirath_tree_get_revealed_nodes(revealed, &revealed_length, tree, challenge_list, challenge_length);
    for (size_t i = 0; i < revealed_length; i++) {
        memcpy(&output[i * MIRATH_SECURITY_BYTES], tree->nodes[revealed[i]], MIRATH_SECURITY_BYTES);
    }
    return revealed_length;
}


/**
 * \fn int mirath_tree_get_seeds_from_path(mirath_tree_t* tree,
                                  const uint16_t* challenge_list, size_t challenge_length,
                                  const uint8_t* path, size_t path_length,
                                  const uint8_t *salt, uint8_t e)
 * \brief This function reconstructs the partial tree determined by the sibling path.
 *
 * \param[out] tree mirath_tree_t* Representation of the partial tree determined by the sibling path
 * \param[in] challenge_list uint16_t* List of challenges (node positions to be hidden)
 * \param[in] challenge_length sizte_t Number of challenges
 * \param[in] path uint8_t* Sibling path
 * \param[in] path_length uint8_t Lenght of the sibling path
 * \param[in] salt uint8_t* Salt used for the signature
 * \param[in] e uint8_t Current iteration
 */
int mirath_tree_get_seeds_from_path(mirath_tree_t* tree,
                                  const uint16_t* challenge_list, size_t challenge_length,
                                  const uint8_t* path, size_t path_length,
                                  const uint8_t *salt, uint8_t e) {
    int ret = 0;
    size_t revealed[MIRATH_PARAM_TREE_LEAVES] = {0}; // We assume the worst case scenario
    size_t revealed_length;

    mirath_tree_get_revealed_nodes(revealed, &revealed_length, tree, challenge_list, challenge_length);
    if (revealed_length != path_length) {
        ret = -1;
        goto Exit;
    }
    for (size_t i = 0; i < revealed_length; i++) {
        memcpy(tree->nodes[revealed[i]], &path[i * MIRATH_SECURITY_BYTES], MIRATH_SECURITY_BYTES);
        tree->nonempty[revealed[i]] = 1;
    }

    mirath_tree_prg(tree, salt, e);

Exit:
    memset(revealed, 0, MIRATH_PARAM_TREE_LEAVES * sizeof(size_t));
    return ret;
}


/**
 * \fn void mirath_tree_get_leaves(mirath_tree_leaves_1_t output_1, mirath_tree_leaves_2_t output_2, mirath_tree_t tree)
 * \brief This function returns the leaves as two list of seeds:.
 *
 * \param[out] output mirath_tree_leaves_t Representation of leaves
 * \param[in] tree mirath_tree_t* Representation of the (partial) tree
 */
void mirath_tree_get_leaves(mirath_tree_leaves_t output, mirath_tree_t *tree) {
    size_t first_leaf = MIRATH_PARAM_TREE_NODES - MIRATH_PARAM_TREE_LEAVES;
    for (size_t i = first_leaf; i < MIRATH_PARAM_TREE_NODES; i++) {
        memcpy(&(output[i - first_leaf]), &(tree->nodes[i]), MIRATH_SECURITY_BYTES);
    }
}

/**
 * \fn void mirath_tree_print(mirath_tree_t tree)
 * \brief This function prints the input tree.
 *
 * \param[in] tree mirath_tree_t* Representation of the (partial) tree
 */
void mirath_tree_print(mirath_tree_t *tree) {
    printf("\n");
    for (size_t i = 0; i < MIRATH_PARAM_TREE_NODES; i++) {
        printf("seed %04lu [nonempty: %d; exists: %d]", i, tree->nonempty[i], tree->exists[i]);
        for (size_t j = 0; j < MIRATH_SECURITY_BYTES; j++) {
            printf(" %02X", tree->nodes[i][j]);
        }
        printf("\n");
    }
}


/**
 * \fn void mirath_tree_print_leaves(mirath_tree_t tree)
 * \brief This function prints the leaves of the input tree.
 *
 * \param[in] tree mirath_tree_t* Representation of the (partial) tree
 */
void mirath_tree_print_leaves(mirath_tree_t *tree) {
    printf("\n");
    size_t first_leaf = MIRATH_PARAM_TREE_NODES - MIRATH_PARAM_TREE_LEAVES;
    for (size_t i = first_leaf; i < MIRATH_PARAM_TREE_NODES; i++) {
        printf("               + seed #%04lu: ", i - first_leaf);
        for (size_t j = 0; j < MIRATH_SECURITY_BYTES; j++) {
            printf("%02X", tree->nodes[i][j]);
        }
        printf("\n");
    }
}

/**
 * \fn void mirath_tree_print_sibling_path(const uint8_t *path, size_t length)
 * \brief This function prints the sibling path that hides the challenges.
 *
 * \param[in] path const uint8_t* Sibling path that hides the challenges
 * \param[in] length uint8_t* Sibling path length
 */
void mirath_tree_print_sibling_path(const uint8_t *path, size_t length) {
    printf("\n");
    for (size_t i = 0; i < length; i++) {
        printf("               + seed #%04lu: ", i);
        for (size_t j = 0; j < MIRATH_SECURITY_BYTES; j++) {
            printf("%02X", path[i * MIRATH_SECURITY_BYTES + j]);
        }
        printf("\n");
    }
}
