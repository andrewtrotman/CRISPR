/*
	ENCODE_KMER_2BIT.CPP
	--------------------
	Copyright (c) 2023 Andrew Trotman
*/
#include "encode_kmer_2bit.h"

/*
	The translation table for converign ascii into bits and bits into ascii
*/
uint64_t encode_kmer_2bit::kmer_encoding_table[256];

