/**
 * @file mirath_tree.h
 * @brief Header file for mirath_tree.c
 */

#ifndef MIRATH_TREE_H
#define MIRATH_TREE_H

#include <stdint.h>
#include <stdlib.h>
#include "common/fips202.h"
#include "mirath_parameters.h"


typedef uint8_t mirath_tree_node_t[MIRATH_SECURITY_BYTES];

typedef struct mirath_tree_t {
    mirath_tree_node_t nodes[MIRATH_PARAM_TREE_NODES];
    uint8_t nonempty[MIRATH_PARAM_TREE_NODES];
    uint8_t exists[MIRATH_PARAM_TREE_NODES];
} mirath_tree_t;

int mirath_tree_node_contains(size_t* list, size_t len, size_t value);
int mirath_tree_node_exists(const mirath_tree_t *tree, size_t i);
void mirath_tree_prg(mirath_tree_t *tree, const uint8_t *salt, uint8_t e);

int mirath_tree_has_right_child(const mirath_tree_t *tree, size_t node);
int mirath_tree_is_a_leaf(size_t node);
size_t mirath_tree_get_parent(size_t node);
int mirath_tree_has_sibling(const mirath_tree_t* tree, size_t node);
size_t mirath_tree_get_sibling(size_t node);
void mirath_tree_get_revealed_nodes(size_t *output_revealed, size_t* output_length,
                                  const mirath_tree_t *tree,
                                  const uint16_t *challenge_list, size_t challenge_length);
size_t mirath_tree_get_sibling_path(uint8_t* output,
                                  const mirath_tree_t* tree,
                                  uint16_t *challenge_list, size_t challenge_length);

int mirath_tree_get_seeds_from_path(mirath_tree_t* tree,
                                  const uint16_t* challenge_list, size_t challenge_length,
                                  const uint8_t* path, size_t path_length,
                                  const uint8_t *salt, uint8_t e);

typedef mirath_tree_node_t mirath_tree_leaves_t[MIRATH_PARAM_TREE_LEAVES];
void mirath_tree_get_leaves(mirath_tree_leaves_t output, mirath_tree_t *tree);

void mirath_tree_init(mirath_tree_t *tree);
void mirath_tree_clear(mirath_tree_t *tree);
void mirath_tree_print(mirath_tree_t *tree);
void mirath_tree_print_leaves(mirath_tree_t *tree);
void mirath_tree_print_sibling_path(const uint8_t *path, size_t length);

#endif
