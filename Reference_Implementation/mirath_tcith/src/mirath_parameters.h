/**
 * \file parameters.h
 * \brief Parameters of the MIRATH scheme
 */

#ifndef MIRATH_PARAMETER_H
#define MIRATH_PARAMETER_H

#define MIRATH_SECURITY 128			        /**< Expected security level (bits) */
#define MIRATH_SECURITY_BYTES 16		    /**< Expected security level (bytes) */

#define MIRATH_PARAM_N_1 256			    /**< Parameter N_1 of the scheme */
#define MIRATH_PARAM_N_2 128			    /**< Parameter N_2 of the scheme */
#define MIRATH_PARAM_N_1_BITS 8		        /**< Parameter tau_1 of the scheme (bits) */
#define MIRATH_PARAM_N_2_BITS 7		        /**< Parameter tau_2 of the scheme (bits) */
#define MIRATH_PARAM_N_1_BYTES 1		    /**< Parameter tau_1 of the scheme (bytes) */
#define MIRATH_PARAM_N_2_BYTES 1		    /**< Parameter tau_2 of the scheme (bytes) */
#define MIRATH_PARAM_N_1_MASK 0xff		    /**< Parameter tau_1 of the scheme (mask) */
#define MIRATH_PARAM_N_2_MASK 0x7f		    /**< Parameter tau_2 of the scheme (mask) */
#define MIRATH_PARAM_TAU 20			        /**< Parameter tau of the scheme (number of iterations) */
#define MIRATH_PARAM_TAU_1 4			    /**< Parameter tau_1 of the scheme (number of iterations concerning N1) */
#define MIRATH_PARAM_TAU_2 16		        /**< Parameter tau_2 of the scheme (number of iterations concerning N2) */
#define MIRATH_PARAM_RHO 32			        /**< Parameter rho of the scheme (dimension of the extension)*/

// TODO set the correct values
#define MIRATH_PARAM_M 16
#define MIRATH_PARAM_N 16
#define MIRATH_PARAM_K 64
#define MIRATH_PARAM_R 8
#define MIRATH_N_PARTIES 256

#define MIRATH_PARAM_TREE_NODES 7167		/**< Number of nodes in the tree */
#define MIRATH_PARAM_TREE_LEAVES 3072	    /**< Number of leaves in the tree */
#define MIRATH_PARAM_TREE_DEPTH 12		    /**< Depth of the tree */

#define MIRATH_PARAM_CHALLENGE_2_BYTES 18	/**< Number of bytes required to store the second challenge */
#define MIRATH_PARAM_HASH_2_MASK_BYTES 1	/**< Number of last bytes in the second hash to be zero */
#define MIRATH_PARAM_HASH_2_MASK 0xe0	    /**< Mask for the last byte in the second hash to be zero */
#define MIRATH_PARAM_T_OPEN 113		        /**< Maximum sibling path length allowed */

#define DOMAIN_SEPARATOR_MESSAGE 0
#define DOMAIN_SEPARATOR_HASH1 1
#define DOMAIN_SEPARATOR_HASH2 2
#define DOMAIN_SEPARATOR_TREE 3
#define DOMAIN_SEPARATOR_COMMITMENT 4

#endif

