/**
 * \file parameters.h
 * \brief Parameters of the MIRATH scheme
 */

#ifndef MIRATH_PARAMETER_H
#define MIRATH_PARAMETER_H

#if defined(_Ia_)
#define MIRATH_SECURITY 128			        /**< Expected security level (bits) */
#define MIRATH_SECURITY_BYTES 16		    /**< Expected security level (bytes) */

#define MIRATH_PARAM_SALT_BYTES (2 * MIRATH_SECURITY_BYTES)

#define MIRATH_SECRET_KEY_BYTES 32		    /**< Secret key size */
#define MIRATH_PUBLIC_KEY_BYTES 73		    /**< Public key size */
#define MIRATH_SIGNATURE_BYTES 3960		    /**< Signature size */

#define MIRATH_PARAM_Q 16			        /**< Parameter q of the scheme (finite field GF(q^m)) */
#define MIRATH_PARAM_M 16			        /**< Parameter m of the scheme (finite field GF(q^m)) */
#define MIRATH_PARAM_K 143			        /**< Parameter k of the scheme (code dimension) */
#define MIRATH_PARAM_N 16			        /**< Parameter n of the scheme (code length) */
#define MIRATH_PARAM_R 4			        /**< Parameter r of the scheme (rank of vectors) */

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
#define MIRATH_PARAM_RHO 16			        /**< Parameter rho of the scheme (dimension of the extension)*/

#define MIRATH_PARAM_TREE_NODES 7167		/**< Number of nodes in the tree */
#define MIRATH_PARAM_TREE_LEAVES 3072	    /**< Number of leaves in the tree */
#define MIRATH_PARAM_TREE_DEPTH 12		    /**< Depth of the tree */

#define MIRATH_PARAM_CHALLENGE_2_BYTES 18	/**< Number of bytes required to store the second challenge */
#define MIRATH_PARAM_HASH_2_MASK_BYTES 1	/**< Number of last bytes in the second hash to be zero */
#define MIRATH_PARAM_HASH_2_MASK 0xe0	    /**< Mask for the last byte in the second hash to be zero */
#define MIRATH_PARAM_T_OPEN 113		        /**< Maximum sibling path length allowed */

#elif defined(_IIIa_)
#define MIRATH_SECURITY 192			        /**< Expected security level (bits) */
#define MIRATH_SECURITY_BYTES 24		    /**< Expected security level (bytes) */

#define MIRATH_PARAM_SALT_BYTES (2 * MIRATH_SECURITY_BYTES)

#define MIRATH_SECRET_KEY_BYTES 48		    /**< Secret key size */
#define MIRATH_PUBLIC_KEY_BYTES 101		    /**< Public key size */
#define MIRATH_SIGNATURE_BYTES 8816		    /**< Signature size */

#define MIRATH_PARAM_Q 16			        /**< Parameter q of the scheme (finite field GF(q^m)) */
#define MIRATH_PARAM_M 21			        /**< Parameter m of the scheme (finite field GF(q^m)) */
#define MIRATH_PARAM_K 288			        /**< Parameter k of the scheme (code dimension) */
#define MIRATH_PARAM_N 21			        /**< Parameter n of the scheme (code length) */
#define MIRATH_PARAM_R 4			        /**< Parameter r of the scheme (rank of vectors) */

#define MIRATH_PARAM_N_1 256			    /**< Parameter N_1 of the scheme */
#define MIRATH_PARAM_N_2 128			    /**< Parameter N_2 of the scheme */
#define MIRATH_PARAM_N_1_BITS 8		        /**< Parameter tau_1 of the scheme (bits) */
#define MIRATH_PARAM_N_2_BITS 7		        /**< Parameter tau_2 of the scheme (bits) */
#define MIRATH_PARAM_N_1_BYTES 1		    /**< Parameter tau_1 of the scheme (bytes) */
#define MIRATH_PARAM_N_2_BYTES 1		    /**< Parameter tau_2 of the scheme (bytes) */
#define MIRATH_PARAM_N_1_MASK 0xff		    /**< Parameter tau_1 of the scheme (mask) */
#define MIRATH_PARAM_N_2_MASK 0x7f		    /**< Parameter tau_2 of the scheme (mask) */
#define MIRATH_PARAM_TAU 30			        /**< Parameter tau of the scheme (number of iterations) */
#define MIRATH_PARAM_TAU_1 10			    /**< Parameter tau_1 of the scheme (number of iterations concerning N1) */
#define MIRATH_PARAM_TAU_2 20		        /**< Parameter tau_2 of the scheme (number of iterations concerning N2) */
#define MIRATH_PARAM_RHO 24			        /**< Parameter rho of the scheme (dimension of the extension)*/

#define MIRATH_PARAM_TREE_NODES 13311		/**< Number of nodes in the tree */
#define MIRATH_PARAM_TREE_LEAVES 5120	    /**< Number of leaves in the tree */
#define MIRATH_PARAM_TREE_DEPTH 13		    /**< Depth of the tree */

#define MIRATH_PARAM_CHALLENGE_2_BYTES 28	/**< Number of bytes required to store the second challenge */
#define MIRATH_PARAM_HASH_2_MASK_BYTES 1	/**< Number of last bytes in the second hash to be zero */
#define MIRATH_PARAM_HASH_2_MASK 0x80	    /**< Mask for the last byte in the second hash to be zero */
#define MIRATH_PARAM_T_OPEN 178		        /**< Maximum sibling path length allowed */

#elif defined(_Va_)
#define MIRATH_SECURITY 256			        /**< Expected security level (bits) */
#define MIRATH_SECURITY_BYTES 32		    /**< Expected security level (bytes) */

#define MIRATH_PARAM_SALT_BYTES (2 * MIRATH_SECURITY_BYTES)

#define MIRATH_SECRET_KEY_BYTES 64		    /**< Secret key size */
#define MIRATH_PUBLIC_KEY_BYTES 125		    /**< Public key size */
#define MIRATH_SIGNATURE_BYTES 15372		/**< Signature size */

#define MIRATH_PARAM_Q 16			        /**< Parameter q of the scheme (finite field GF(q^m)) */
#define MIRATH_PARAM_M 25			        /**< Parameter m of the scheme (finite field GF(q^m)) */
#define MIRATH_PARAM_K 440			        /**< Parameter k of the scheme (code dimension) */
#define MIRATH_PARAM_N 25			        /**< Parameter n of the scheme (code length) */
#define MIRATH_PARAM_R 4			        /**< Parameter r of the scheme (rank of vectors) */

#define MIRATH_PARAM_N_1 256			    /**< Parameter N_1 of the scheme */
#define MIRATH_PARAM_N_2 128			    /**< Parameter N_2 of the scheme */
#define MIRATH_PARAM_N_1_BITS 8		        /**< Parameter tau_1 of the scheme (bits) */
#define MIRATH_PARAM_N_2_BITS 7		        /**< Parameter tau_2 of the scheme (bits) */
#define MIRATH_PARAM_N_1_BYTES 1		    /**< Parameter tau_1 of the scheme (bytes) */
#define MIRATH_PARAM_N_2_BYTES 1		    /**< Parameter tau_2 of the scheme (bytes) */
#define MIRATH_PARAM_N_1_MASK 0xff		    /**< Parameter tau_1 of the scheme (mask) */
#define MIRATH_PARAM_N_2_MASK 0x7f		    /**< Parameter tau_2 of the scheme (mask) */
#define MIRATH_PARAM_TAU 39			        /**< Parameter tau of the scheme (number of iterations) */
#define MIRATH_PARAM_TAU_1 17			    /**< Parameter tau_1 of the scheme (number of iterations concerning N1) */
#define MIRATH_PARAM_TAU_2 22		        /**< Parameter tau_2 of the scheme (number of iterations concerning N2) */
#define MIRATH_PARAM_RHO 32			        /**< Parameter rho of the scheme (dimension of the extension)*/

#define MIRATH_PARAM_TREE_NODES 15359		/**< Number of nodes in the tree */
#define MIRATH_PARAM_TREE_LEAVES 7168	    /**< Number of leaves in the tree */
#define MIRATH_PARAM_TREE_DEPTH 13		    /**< Depth of the tree */

#define MIRATH_PARAM_CHALLENGE_2_BYTES 37	/**< Number of bytes required to store the second challenge */
#define MIRATH_PARAM_HASH_2_MASK_BYTES 1	/**< Number of last bytes in the second hash to be zero */
#define MIRATH_PARAM_HASH_2_MASK 0xf0	    /**< Mask for the last byte in the second hash to be zero */
#define MIRATH_PARAM_T_OPEN 247		        /**< Maximum sibling path length allowed */
#endif

#define DOMAIN_SEPARATOR_MESSAGE 0
#define DOMAIN_SEPARATOR_HASH1 1
#define DOMAIN_SEPARATOR_HASH2_PARTIAL 2
#define DOMAIN_SEPARATOR_HASH2 3
#define DOMAIN_SEPARATOR_TREE 4
#define DOMAIN_SEPARATOR_COMMITMENT 5

#define MIRATH_VAR_GAMMA (MIRATH_PARAM_RHO * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K))
#define MIRATH_VAR_S (MIRATH_PARAM_M * MIRATH_PARAM_R)
#define MIRATH_VAR_C (MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_BASE_MID (MIRATH_PARAM_M * (MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_E_A (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K)
#define MIRATH_VAR_T (MIRATH_PARAM_M * MIRATH_PARAM_N)

#endif

