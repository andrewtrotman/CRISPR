/*
	ENCODE_KMER_3BIT.H
	------------------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief Construct and manipulate kmers encoded in 3 bits per base (A=001, C=010, G=100, T=111)
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/
#include <stdint.h>

#include <string>

#include <iostream>

#include "encode_kmer_2bit.h"

/*
	CLASS ENCODE_KMER_3BIT
	----------------------
*/
/*!
	@brief Construct and manipulate kmers encoded in 3 bits per base (A=001, C=010, G=100, T=111)
*/
class encode_kmer_3bit
	{
	private:
		/*
			Fast lookup for encoding encoded kmers
		*/
		static uint64_t kmer_encoding_table[256];

	public:
		/*
			ENCODE_KMER_3BIT::ENCODE_KMER_2BIT()
			------------------------------------
		*/
		/*!
			@brief Constructor
		*/
		encode_kmer_3bit()
			{
			kmer_encoding_table[(size_t)'A'] = 1;
			kmer_encoding_table[(size_t)'C'] = 2;
			kmer_encoding_table[(size_t)'G'] = 4;
			kmer_encoding_table[(size_t)'T'] = 7;

			kmer_encoding_table[(size_t)'a'] = 1;
			kmer_encoding_table[(size_t)'c'] = 2;
			kmer_encoding_table[(size_t)'g'] = 4;
			kmer_encoding_table[(size_t)'t'] = 7;

			kmer_encoding_table[1] = 'A';
			kmer_encoding_table[2] = 'C';
			kmer_encoding_table[4] = 'G';
			kmer_encoding_table[7] = 'T';
			}

		/*
			ENCODE_KMER_3BIT::PACK_20MER()
			------------------------------
		*/
		/*!
			@brief Pack a 20-mer DNA sequence (consisting of just ACGT characters) int a 64-bit integer
			@param sequence [in] The DNA sequence to pack.
			@returns The 20-mer packed into 64-bits
		*/
		static uint64_t pack_20mer(const char *sequence)
			{
			uint64_t packed = 0;

			for (int pos = 0; pos < 20; pos++)
				packed = (packed << 3) | kmer_encoding_table[(size_t)sequence[pos]];

			return packed;
			}

		/*
			ENCODE_KMER_3BIT::UNPACK_20MER()
			--------------------------------
		*/
		/*!
			@brief Turn a packed 20-mer into an ASCII representation (ACGT).
			@param packed_sequence [in] The packed 20-mer.
			@returns The DNA sequence as a 20-character long string.
		*/
		static void unpack_20mer(char *into, uint64_t packed_sequence)
			{
			for (int32_t pos = 19; pos >= 0; pos--)
				*into++ = kmer_encoding_table[(packed_sequence >> (pos * 3)) & 7];

			*into = '\0';
			}

		/*
			ENCODE_KMER_3BIT::UNPACK_20MER()
			--------------------------------
		*/
		/*!
			@brief Turn a packed 20-mer into an ASCII representation (ACGT).
			@param packed_sequence [in] The packed 20-mer.
			@returns The DNA sequence as a 20-character long string.
		*/
		static std::string unpack_20mer(uint64_t packed_sequence)
			{
			std::string sequence;

			for (int32_t pos = 19; pos >= 0; pos--)
				sequence +=kmer_encoding_table[(packed_sequence >> (pos * 3)) & 7];

			return sequence;
			}

		/*
			ENCODE_KMER_3BIT::PACK_20MER_INTO_2BIT()
			----------------------------------------
		*/
		/*!
			@brief Turn a packed 20-mer into a 2-bit packed 20-mer
			@param packed_sequence [in] The packed 20-mer.
			@returns The DNA sequence as a 2-bit encoding.
		*/
		static uint64_t pack_20mer_into_2bit(uint64_t packed_sequence)
			{
			uint64_t re_packed = 0;

			for (size_t pos = 0; pos < 20; pos++)
				re_packed = (re_packed << 2) | encode_kmer_2bit::kmer_encoding_table[kmer_encoding_table[(packed_sequence >> ((19 - pos) * 3)) & 7]];

			return re_packed;
			}
	};
