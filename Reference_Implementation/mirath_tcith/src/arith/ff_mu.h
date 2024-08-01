#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include "data_type_arith.h"

/// AES modulus
#define MODULUS 0x1B
#define MASK_LSB_PER_BIT ((uint64_t)0x0101010101010101)
#define MASK_MSB_PER_BIT (MASK_LSB_PER_BIT*0x80)
#define MASK_XLSB_PER_BIT (MASK_LSB_PER_BIT*0xFE)

static const uint8_t mirath_ff_mu_mult_base[] __attribute__((aligned(256)))= {
		// row_nr**1, row_nr**2, ..., row_nr**8
        0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
        0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1b,
        0x03,0x06,0x0c,0x18,0x30,0x60,0xc0,0x9b,
        0x04,0x08,0x10,0x20,0x40,0x80,0x1b,0x36,
        0x05,0x0a,0x14,0x28,0x50,0xa0,0x5b,0xb6,
        0x06,0x0c,0x18,0x30,0x60,0xc0,0x9b,0x2d,
        0x07,0x0e,0x1c,0x38,0x70,0xe0,0xdb,0xad,
        0x08,0x10,0x20,0x40,0x80,0x1b,0x36,0x6c,
        0x09,0x12,0x24,0x48,0x90,0x3b,0x76,0xec,
        0x0a,0x14,0x28,0x50,0xa0,0x5b,0xb6,0x77,
        0x0b,0x16,0x2c,0x58,0xb0,0x7b,0xf6,0xf7,
        0x0c,0x18,0x30,0x60,0xc0,0x9b,0x2d,0x5a,
        0x0d,0x1a,0x34,0x68,0xd0,0xbb,0x6d,0xda,
        0x0e,0x1c,0x38,0x70,0xe0,0xdb,0xad,0x41,
        0x0f,0x1e,0x3c,0x78,0xf0,0xfb,0xed,0xc1,
        0x10,0x20,0x40,0x80,0x1b,0x36,0x6c,0xd8,
        0x11,0x22,0x44,0x88,0x0b,0x16,0x2c,0x58,
        0x12,0x24,0x48,0x90,0x3b,0x76,0xec,0xc3,
        0x13,0x26,0x4c,0x98,0x2b,0x56,0xac,0x43,
        0x14,0x28,0x50,0xa0,0x5b,0xb6,0x77,0xee,
        0x15,0x2a,0x54,0xa8,0x4b,0x96,0x37,0x6e,
        0x16,0x2c,0x58,0xb0,0x7b,0xf6,0xf7,0xf5,
        0x17,0x2e,0x5c,0xb8,0x6b,0xd6,0xb7,0x75,
        0x18,0x30,0x60,0xc0,0x9b,0x2d,0x5a,0xb4,
        0x19,0x32,0x64,0xc8,0x8b,0x0d,0x1a,0x34,
        0x1a,0x34,0x68,0xd0,0xbb,0x6d,0xda,0xaf,
        0x1b,0x36,0x6c,0xd8,0xab,0x4d,0x9a,0x2f,
        0x1c,0x38,0x70,0xe0,0xdb,0xad,0x41,0x82,
        0x1d,0x3a,0x74,0xe8,0xcb,0x8d,0x01,0x02,
        0x1e,0x3c,0x78,0xf0,0xfb,0xed,0xc1,0x99,
        0x1f,0x3e,0x7c,0xf8,0xeb,0xcd,0x81,0x19,
        0x20,0x40,0x80,0x1b,0x36,0x6c,0xd8,0xab,
        0x21,0x42,0x84,0x13,0x26,0x4c,0x98,0x2b,
        0x22,0x44,0x88,0x0b,0x16,0x2c,0x58,0xb0,
        0x23,0x46,0x8c,0x03,0x06,0x0c,0x18,0x30,
        0x24,0x48,0x90,0x3b,0x76,0xec,0xc3,0x9d,
        0x25,0x4a,0x94,0x33,0x66,0xcc,0x83,0x1d,
        0x26,0x4c,0x98,0x2b,0x56,0xac,0x43,0x86,
        0x27,0x4e,0x9c,0x23,0x46,0x8c,0x03,0x06,
        0x28,0x50,0xa0,0x5b,0xb6,0x77,0xee,0xc7,
        0x29,0x52,0xa4,0x53,0xa6,0x57,0xae,0x47,
        0x2a,0x54,0xa8,0x4b,0x96,0x37,0x6e,0xdc,
        0x2b,0x56,0xac,0x43,0x86,0x17,0x2e,0x5c,
        0x2c,0x58,0xb0,0x7b,0xf6,0xf7,0xf5,0xf1,
        0x2d,0x5a,0xb4,0x73,0xe6,0xd7,0xb5,0x71,
        0x2e,0x5c,0xb8,0x6b,0xd6,0xb7,0x75,0xea,
        0x2f,0x5e,0xbc,0x63,0xc6,0x97,0x35,0x6a,
        0x30,0x60,0xc0,0x9b,0x2d,0x5a,0xb4,0x73,
        0x31,0x62,0xc4,0x93,0x3d,0x7a,0xf4,0xf3,
        0x32,0x64,0xc8,0x8b,0x0d,0x1a,0x34,0x68,
        0x33,0x66,0xcc,0x83,0x1d,0x3a,0x74,0xe8,
        0x34,0x68,0xd0,0xbb,0x6d,0xda,0xaf,0x45,
        0x35,0x6a,0xd4,0xb3,0x7d,0xfa,0xef,0xc5,
        0x36,0x6c,0xd8,0xab,0x4d,0x9a,0x2f,0x5e,
        0x37,0x6e,0xdc,0xa3,0x5d,0xba,0x6f,0xde,
        0x38,0x70,0xe0,0xdb,0xad,0x41,0x82,0x1f,
        0x39,0x72,0xe4,0xd3,0xbd,0x61,0xc2,0x9f,
        0x3a,0x74,0xe8,0xcb,0x8d,0x01,0x02,0x04,
        0x3b,0x76,0xec,0xc3,0x9d,0x21,0x42,0x84,
        0x3c,0x78,0xf0,0xfb,0xed,0xc1,0x99,0x29,
        0x3d,0x7a,0xf4,0xf3,0xfd,0xe1,0xd9,0xa9,
        0x3e,0x7c,0xf8,0xeb,0xcd,0x81,0x19,0x32,
        0x3f,0x7e,0xfc,0xe3,0xdd,0xa1,0x59,0xb2,
        0x40,0x80,0x1b,0x36,0x6c,0xd8,0xab,0x4d,
        0x41,0x82,0x1f,0x3e,0x7c,0xf8,0xeb,0xcd,
        0x42,0x84,0x13,0x26,0x4c,0x98,0x2b,0x56,
        0x43,0x86,0x17,0x2e,0x5c,0xb8,0x6b,0xd6,
        0x44,0x88,0x0b,0x16,0x2c,0x58,0xb0,0x7b,
        0x45,0x8a,0x0f,0x1e,0x3c,0x78,0xf0,0xfb,
        0x46,0x8c,0x03,0x06,0x0c,0x18,0x30,0x60,
        0x47,0x8e,0x07,0x0e,0x1c,0x38,0x70,0xe0,
        0x48,0x90,0x3b,0x76,0xec,0xc3,0x9d,0x21,
        0x49,0x92,0x3f,0x7e,0xfc,0xe3,0xdd,0xa1,
        0x4a,0x94,0x33,0x66,0xcc,0x83,0x1d,0x3a,
        0x4b,0x96,0x37,0x6e,0xdc,0xa3,0x5d,0xba,
        0x4c,0x98,0x2b,0x56,0xac,0x43,0x86,0x17,
        0x4d,0x9a,0x2f,0x5e,0xbc,0x63,0xc6,0x97,
        0x4e,0x9c,0x23,0x46,0x8c,0x03,0x06,0x0c,
        0x4f,0x9e,0x27,0x4e,0x9c,0x23,0x46,0x8c,
        0x50,0xa0,0x5b,0xb6,0x77,0xee,0xc7,0x95,
        0x51,0xa2,0x5f,0xbe,0x67,0xce,0x87,0x15,
        0x52,0xa4,0x53,0xa6,0x57,0xae,0x47,0x8e,
        0x53,0xa6,0x57,0xae,0x47,0x8e,0x07,0x0e,
        0x54,0xa8,0x4b,0x96,0x37,0x6e,0xdc,0xa3,
        0x55,0xaa,0x4f,0x9e,0x27,0x4e,0x9c,0x23,
        0x56,0xac,0x43,0x86,0x17,0x2e,0x5c,0xb8,
        0x57,0xae,0x47,0x8e,0x07,0x0e,0x1c,0x38,
        0x58,0xb0,0x7b,0xf6,0xf7,0xf5,0xf1,0xf9,
        0x59,0xb2,0x7f,0xfe,0xe7,0xd5,0xb1,0x79,
        0x5a,0xb4,0x73,0xe6,0xd7,0xb5,0x71,0xe2,
        0x5b,0xb6,0x77,0xee,0xc7,0x95,0x31,0x62,
        0x5c,0xb8,0x6b,0xd6,0xb7,0x75,0xea,0xcf,
        0x5d,0xba,0x6f,0xde,0xa7,0x55,0xaa,0x4f,
        0x5e,0xbc,0x63,0xc6,0x97,0x35,0x6a,0xd4,
        0x5f,0xbe,0x67,0xce,0x87,0x15,0x2a,0x54,
        0x60,0xc0,0x9b,0x2d,0x5a,0xb4,0x73,0xe6,
        0x61,0xc2,0x9f,0x25,0x4a,0x94,0x33,0x66,
        0x62,0xc4,0x93,0x3d,0x7a,0xf4,0xf3,0xfd,
        0x63,0xc6,0x97,0x35,0x6a,0xd4,0xb3,0x7d,
        0x64,0xc8,0x8b,0x0d,0x1a,0x34,0x68,0xd0,
        0x65,0xca,0x8f,0x05,0x0a,0x14,0x28,0x50,
        0x66,0xcc,0x83,0x1d,0x3a,0x74,0xe8,0xcb,
        0x67,0xce,0x87,0x15,0x2a,0x54,0xa8,0x4b,
        0x68,0xd0,0xbb,0x6d,0xda,0xaf,0x45,0x8a,
        0x69,0xd2,0xbf,0x65,0xca,0x8f,0x05,0x0a,
        0x6a,0xd4,0xb3,0x7d,0xfa,0xef,0xc5,0x91,
        0x6b,0xd6,0xb7,0x75,0xea,0xcf,0x85,0x11,
        0x6c,0xd8,0xab,0x4d,0x9a,0x2f,0x5e,0xbc,
        0x6d,0xda,0xaf,0x45,0x8a,0x0f,0x1e,0x3c,
        0x6e,0xdc,0xa3,0x5d,0xba,0x6f,0xde,0xa7,
        0x6f,0xde,0xa7,0x55,0xaa,0x4f,0x9e,0x27,
        0x70,0xe0,0xdb,0xad,0x41,0x82,0x1f,0x3e,
        0x71,0xe2,0xdf,0xa5,0x51,0xa2,0x5f,0xbe,
        0x72,0xe4,0xd3,0xbd,0x61,0xc2,0x9f,0x25,
        0x73,0xe6,0xd7,0xb5,0x71,0xe2,0xdf,0xa5,
        0x74,0xe8,0xcb,0x8d,0x01,0x02,0x04,0x08,
        0x75,0xea,0xcf,0x85,0x11,0x22,0x44,0x88,
        0x76,0xec,0xc3,0x9d,0x21,0x42,0x84,0x13,
        0x77,0xee,0xc7,0x95,0x31,0x62,0xc4,0x93,
        0x78,0xf0,0xfb,0xed,0xc1,0x99,0x29,0x52,
        0x79,0xf2,0xff,0xe5,0xd1,0xb9,0x69,0xd2,
        0x7a,0xf4,0xf3,0xfd,0xe1,0xd9,0xa9,0x49,
        0x7b,0xf6,0xf7,0xf5,0xf1,0xf9,0xe9,0xc9,
        0x7c,0xf8,0xeb,0xcd,0x81,0x19,0x32,0x64,
        0x7d,0xfa,0xef,0xc5,0x91,0x39,0x72,0xe4,
        0x7e,0xfc,0xe3,0xdd,0xa1,0x59,0xb2,0x7f,
        0x7f,0xfe,0xe7,0xd5,0xb1,0x79,0xf2,0xff,
        0x80,0x1b,0x36,0x6c,0xd8,0xab,0x4d,0x9a,
        0x81,0x19,0x32,0x64,0xc8,0x8b,0x0d,0x1a,
        0x82,0x1f,0x3e,0x7c,0xf8,0xeb,0xcd,0x81,
        0x83,0x1d,0x3a,0x74,0xe8,0xcb,0x8d,0x01,
        0x84,0x13,0x26,0x4c,0x98,0x2b,0x56,0xac,
        0x85,0x11,0x22,0x44,0x88,0x0b,0x16,0x2c,
        0x86,0x17,0x2e,0x5c,0xb8,0x6b,0xd6,0xb7,
        0x87,0x15,0x2a,0x54,0xa8,0x4b,0x96,0x37,
        0x88,0x0b,0x16,0x2c,0x58,0xb0,0x7b,0xf6,
        0x89,0x09,0x12,0x24,0x48,0x90,0x3b,0x76,
        0x8a,0x0f,0x1e,0x3c,0x78,0xf0,0xfb,0xed,
        0x8b,0x0d,0x1a,0x34,0x68,0xd0,0xbb,0x6d,
        0x8c,0x03,0x06,0x0c,0x18,0x30,0x60,0xc0,
        0x8d,0x01,0x02,0x04,0x08,0x10,0x20,0x40,
        0x8e,0x07,0x0e,0x1c,0x38,0x70,0xe0,0xdb,
        0x8f,0x05,0x0a,0x14,0x28,0x50,0xa0,0x5b,
        0x90,0x3b,0x76,0xec,0xc3,0x9d,0x21,0x42,
        0x91,0x39,0x72,0xe4,0xd3,0xbd,0x61,0xc2,
        0x92,0x3f,0x7e,0xfc,0xe3,0xdd,0xa1,0x59,
        0x93,0x3d,0x7a,0xf4,0xf3,0xfd,0xe1,0xd9,
        0x94,0x33,0x66,0xcc,0x83,0x1d,0x3a,0x74,
        0x95,0x31,0x62,0xc4,0x93,0x3d,0x7a,0xf4,
        0x96,0x37,0x6e,0xdc,0xa3,0x5d,0xba,0x6f,
        0x97,0x35,0x6a,0xd4,0xb3,0x7d,0xfa,0xef,
        0x98,0x2b,0x56,0xac,0x43,0x86,0x17,0x2e,
        0x99,0x29,0x52,0xa4,0x53,0xa6,0x57,0xae,
        0x9a,0x2f,0x5e,0xbc,0x63,0xc6,0x97,0x35,
        0x9b,0x2d,0x5a,0xb4,0x73,0xe6,0xd7,0xb5,
        0x9c,0x23,0x46,0x8c,0x03,0x06,0x0c,0x18,
        0x9d,0x21,0x42,0x84,0x13,0x26,0x4c,0x98,
        0x9e,0x27,0x4e,0x9c,0x23,0x46,0x8c,0x03,
        0x9f,0x25,0x4a,0x94,0x33,0x66,0xcc,0x83,
        0xa0,0x5b,0xb6,0x77,0xee,0xc7,0x95,0x31,
        0xa1,0x59,0xb2,0x7f,0xfe,0xe7,0xd5,0xb1,
        0xa2,0x5f,0xbe,0x67,0xce,0x87,0x15,0x2a,
        0xa3,0x5d,0xba,0x6f,0xde,0xa7,0x55,0xaa,
        0xa4,0x53,0xa6,0x57,0xae,0x47,0x8e,0x07,
        0xa5,0x51,0xa2,0x5f,0xbe,0x67,0xce,0x87,
        0xa6,0x57,0xae,0x47,0x8e,0x07,0x0e,0x1c,
        0xa7,0x55,0xaa,0x4f,0x9e,0x27,0x4e,0x9c,
        0xa8,0x4b,0x96,0x37,0x6e,0xdc,0xa3,0x5d,
        0xa9,0x49,0x92,0x3f,0x7e,0xfc,0xe3,0xdd,
        0xaa,0x4f,0x9e,0x27,0x4e,0x9c,0x23,0x46,
        0xab,0x4d,0x9a,0x2f,0x5e,0xbc,0x63,0xc6,
        0xac,0x43,0x86,0x17,0x2e,0x5c,0xb8,0x6b,
        0xad,0x41,0x82,0x1f,0x3e,0x7c,0xf8,0xeb,
        0xae,0x47,0x8e,0x07,0x0e,0x1c,0x38,0x70,
        0xaf,0x45,0x8a,0x0f,0x1e,0x3c,0x78,0xf0,
        0xb0,0x7b,0xf6,0xf7,0xf5,0xf1,0xf9,0xe9,
        0xb1,0x79,0xf2,0xff,0xe5,0xd1,0xb9,0x69,
        0xb2,0x7f,0xfe,0xe7,0xd5,0xb1,0x79,0xf2,
        0xb3,0x7d,0xfa,0xef,0xc5,0x91,0x39,0x72,
        0xb4,0x73,0xe6,0xd7,0xb5,0x71,0xe2,0xdf,
        0xb5,0x71,0xe2,0xdf,0xa5,0x51,0xa2,0x5f,
        0xb6,0x77,0xee,0xc7,0x95,0x31,0x62,0xc4,
        0xb7,0x75,0xea,0xcf,0x85,0x11,0x22,0x44,
        0xb8,0x6b,0xd6,0xb7,0x75,0xea,0xcf,0x85,
        0xb9,0x69,0xd2,0xbf,0x65,0xca,0x8f,0x05,
        0xba,0x6f,0xde,0xa7,0x55,0xaa,0x4f,0x9e,
        0xbb,0x6d,0xda,0xaf,0x45,0x8a,0x0f,0x1e,
        0xbc,0x63,0xc6,0x97,0x35,0x6a,0xd4,0xb3,
        0xbd,0x61,0xc2,0x9f,0x25,0x4a,0x94,0x33,
        0xbe,0x67,0xce,0x87,0x15,0x2a,0x54,0xa8,
        0xbf,0x65,0xca,0x8f,0x05,0x0a,0x14,0x28,
        0xc0,0x9b,0x2d,0x5a,0xb4,0x73,0xe6,0xd7,
        0xc1,0x99,0x29,0x52,0xa4,0x53,0xa6,0x57,
        0xc2,0x9f,0x25,0x4a,0x94,0x33,0x66,0xcc,
        0xc3,0x9d,0x21,0x42,0x84,0x13,0x26,0x4c,
        0xc4,0x93,0x3d,0x7a,0xf4,0xf3,0xfd,0xe1,
        0xc5,0x91,0x39,0x72,0xe4,0xd3,0xbd,0x61,
        0xc6,0x97,0x35,0x6a,0xd4,0xb3,0x7d,0xfa,
        0xc7,0x95,0x31,0x62,0xc4,0x93,0x3d,0x7a,
        0xc8,0x8b,0x0d,0x1a,0x34,0x68,0xd0,0xbb,
        0xc9,0x89,0x09,0x12,0x24,0x48,0x90,0x3b,
        0xca,0x8f,0x05,0x0a,0x14,0x28,0x50,0xa0,
        0xcb,0x8d,0x01,0x02,0x04,0x08,0x10,0x20,
        0xcc,0x83,0x1d,0x3a,0x74,0xe8,0xcb,0x8d,
        0xcd,0x81,0x19,0x32,0x64,0xc8,0x8b,0x0d,
        0xce,0x87,0x15,0x2a,0x54,0xa8,0x4b,0x96,
        0xcf,0x85,0x11,0x22,0x44,0x88,0x0b,0x16,
        0xd0,0xbb,0x6d,0xda,0xaf,0x45,0x8a,0x0f,
        0xd1,0xb9,0x69,0xd2,0xbf,0x65,0xca,0x8f,
        0xd2,0xbf,0x65,0xca,0x8f,0x05,0x0a,0x14,
        0xd3,0xbd,0x61,0xc2,0x9f,0x25,0x4a,0x94,
        0xd4,0xb3,0x7d,0xfa,0xef,0xc5,0x91,0x39,
        0xd5,0xb1,0x79,0xf2,0xff,0xe5,0xd1,0xb9,
        0xd6,0xb7,0x75,0xea,0xcf,0x85,0x11,0x22,
        0xd7,0xb5,0x71,0xe2,0xdf,0xa5,0x51,0xa2,
        0xd8,0xab,0x4d,0x9a,0x2f,0x5e,0xbc,0x63,
        0xd9,0xa9,0x49,0x92,0x3f,0x7e,0xfc,0xe3,
        0xda,0xaf,0x45,0x8a,0x0f,0x1e,0x3c,0x78,
        0xdb,0xad,0x41,0x82,0x1f,0x3e,0x7c,0xf8,
        0xdc,0xa3,0x5d,0xba,0x6f,0xde,0xa7,0x55,
        0xdd,0xa1,0x59,0xb2,0x7f,0xfe,0xe7,0xd5,
        0xde,0xa7,0x55,0xaa,0x4f,0x9e,0x27,0x4e,
        0xdf,0xa5,0x51,0xa2,0x5f,0xbe,0x67,0xce,
        0xe0,0xdb,0xad,0x41,0x82,0x1f,0x3e,0x7c,
        0xe1,0xd9,0xa9,0x49,0x92,0x3f,0x7e,0xfc,
        0xe2,0xdf,0xa5,0x51,0xa2,0x5f,0xbe,0x67,
        0xe3,0xdd,0xa1,0x59,0xb2,0x7f,0xfe,0xe7,
        0xe4,0xd3,0xbd,0x61,0xc2,0x9f,0x25,0x4a,
        0xe5,0xd1,0xb9,0x69,0xd2,0xbf,0x65,0xca,
        0xe6,0xd7,0xb5,0x71,0xe2,0xdf,0xa5,0x51,
        0xe7,0xd5,0xb1,0x79,0xf2,0xff,0xe5,0xd1,
        0xe8,0xcb,0x8d,0x01,0x02,0x04,0x08,0x10,
        0xe9,0xc9,0x89,0x09,0x12,0x24,0x48,0x90,
        0xea,0xcf,0x85,0x11,0x22,0x44,0x88,0x0b,
        0xeb,0xcd,0x81,0x19,0x32,0x64,0xc8,0x8b,
        0xec,0xc3,0x9d,0x21,0x42,0x84,0x13,0x26,
        0xed,0xc1,0x99,0x29,0x52,0xa4,0x53,0xa6,
        0xee,0xc7,0x95,0x31,0x62,0xc4,0x93,0x3d,
        0xef,0xc5,0x91,0x39,0x72,0xe4,0xd3,0xbd,
        0xf0,0xfb,0xed,0xc1,0x99,0x29,0x52,0xa4,
        0xf1,0xf9,0xe9,0xc9,0x89,0x09,0x12,0x24,
        0xf2,0xff,0xe5,0xd1,0xb9,0x69,0xd2,0xbf,
        0xf3,0xfd,0xe1,0xd9,0xa9,0x49,0x92,0x3f,
        0xf4,0xf3,0xfd,0xe1,0xd9,0xa9,0x49,0x92,
        0xf5,0xf1,0xf9,0xe9,0xc9,0x89,0x09,0x12,
        0xf6,0xf7,0xf5,0xf1,0xf9,0xe9,0xc9,0x89,
        0xf7,0xf5,0xf1,0xf9,0xe9,0xc9,0x89,0x09,
        0xf8,0xeb,0xcd,0x81,0x19,0x32,0x64,0xc8,
        0xf9,0xe9,0xc9,0x89,0x09,0x12,0x24,0x48,
        0xfa,0xef,0xc5,0x91,0x39,0x72,0xe4,0xd3,
        0xfb,0xed,0xc1,0x99,0x29,0x52,0xa4,0x53,
        0xfc,0xe3,0xdd,0xa1,0x59,0xb2,0x7f,0xfe,
        0xfd,0xe1,0xd9,0xa9,0x49,0x92,0x3f,0x7e,
        0xfe,0xe7,0xd5,0xb1,0x79,0xf2,0xff,0xe5,
        0xff,0xe5,0xd1,0xb9,0x69,0xd2,0xbf,0x65,
};

// Warning, getting the inverse of a secret value using this table,
//   would lead to non-constant-time implementation. Only accessing
//   to public positions is allowed.
static const uint8_t mirath_ff_mu_inv_table[256] = {0, 1, 141, 246, 203, 82, 123, 209, 232, 79, 41, 192, 176, 225, 229, 199, 116, 180, 170, 75, 153, 43, 96, 95, 88, 63, 253, 204, 255, 64, 238, 178, 58, 110, 90, 241, 85, 77, 168, 201, 193, 10, 152, 21, 48, 68, 162, 194, 44, 69, 146, 108, 243, 57, 102, 66, 242, 53, 32, 111, 119, 187, 89, 25, 29, 254, 55, 103, 45, 49, 245, 105, 167, 100, 171, 19, 84, 37, 233, 9, 237, 92, 5, 202, 76, 36, 135, 191, 24, 62, 34, 240, 81, 236, 97, 23, 22, 94, 175, 211, 73, 166, 54, 67, 244, 71, 145, 223, 51, 147, 33, 59, 121, 183, 151, 133, 16, 181, 186, 60, 182, 112, 208, 6, 161, 250, 129, 130, 131, 126, 127, 128, 150, 115, 190, 86, 155, 158, 149, 217, 247, 2, 185, 164, 222, 106, 50, 109, 216, 138, 132, 114, 42, 20, 159, 136, 249, 220, 137, 154, 251, 124, 46, 195, 143, 184, 101, 72, 38, 200, 18, 74, 206, 231, 210, 98, 12, 224, 31, 239, 17, 117, 120, 113, 165, 142, 118, 61, 189, 188, 134, 87, 11, 40, 47, 163, 218, 212, 228, 15, 169, 39, 83, 4, 27, 252, 172, 230, 122, 7, 174, 99, 197, 219, 226, 234, 148, 139, 196, 213, 157, 248, 144, 107, 177, 13, 214, 235, 198, 14, 207, 173, 8, 78, 215, 227, 93, 80, 30, 179, 91, 35, 56, 52, 104, 70, 3, 140, 221, 156, 125, 160, 205, 26, 65, 28};

/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a, const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a, const ff_mu_t b) {
    const uint8_t *p = &mirath_ff_mu_mult_base[b * 8];
    uint8_t tmp = 0;
    tmp ^= (a &   1) ? p[0] : 0;
    tmp ^= (a &   2) ? p[1] : 0;
    tmp ^= (a &   4) ? p[2] : 0;
    tmp ^= (a &   8) ? p[3] : 0;
    tmp ^= (a &  16) ? p[4] : 0;
    tmp ^= (a &  32) ? p[5] : 0;
    tmp ^= (a &  64) ? p[6] : 0;
    tmp ^= (a & 128) ? p[7] : 0;
    return tmp;
}

/// \return a^-1
static inline ff_mu_t mirath_ff_mu_inv(const ff_mu_t a) {
    return mirath_ff_mu_inv_table[a];
}

#endif
