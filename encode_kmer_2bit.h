/*
	ENCODE_KMER_2BIT.H
	------------------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief Construct and manipulate kmers encoded in 2 bits per base (A=00, C=01, G=10, T=11)
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <stdint.h>

#include <string>

/*
	CLASS ENCODE_KMER_2BIT
	----------------------
*/
/*!
	@brief Construct and manipulate kmers encoded in 2 bits per base (A=00, C=01, G=10, T=11)
*/
class encode_kmer_2bit
	{
	private:
		/*
			Fast lookup for encoding and decoding kmers
		*/
		static uint64_t kmer_encoding_table[256];

	public:
		/*
			ENCODE_KMER_2BIT::ENCODE_KMER_2BIT()
			------------------------------------
		*/
		/*!
			@brief Constructor
		*/
		encode_kmer_2bit()
			{
			kmer_encoding_table[(size_t)'A'] = 0;
			kmer_encoding_table[(size_t)'C'] = 1;
			kmer_encoding_table[(size_t)'G'] = 2;
			kmer_encoding_table[(size_t)'T'] = 3;

			kmer_encoding_table[(size_t)'a'] = 0;
			kmer_encoding_table[(size_t)'c'] = 1;
			kmer_encoding_table[(size_t)'g'] = 2;
			kmer_encoding_table[(size_t)'t'] = 3;

			kmer_encoding_table[0] = 'A';
			kmer_encoding_table[1] = 'C';
			kmer_encoding_table[2] = 'G';
			kmer_encoding_table[3] = 'T';
			}

		/*
			ENCODE_KMER_2BIT::PACK_20MER()
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

			for (size_t pos = 0; pos < 20; pos++)
				packed = (packed << 2) | kmer_encoding_table[(size_t)sequence[pos]];

			return packed;
			}

		/*
			ENCODE_KMER_2BIT::UNPACK_20MER()
			--------------------------------
		*/
		/*!
			@brief Turn a packed 20-mer into an ASCII representation (ACGT).
			@param into [out] The unpaced 20-mer (a sequence of ACTG characters, '\0' terminated).
			@param packed_sequence [in] The packed 20-mer.
			@returns The DNA sequence as a 20-character long string.
		*/
		static void unpack_20mer(char *into, uint64_t packed_sequence)
			{
			for (int32_t pos = 19; pos >= 0; pos--)
				*into++ = kmer_encoding_table[(packed_sequence >> (pos * 2)) & 3];

			*into = '\0';
			}

		/*
			ENCODE_KMER_2BIT::UNPACK_20MER()
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
				sequence += kmer_encoding_table[(packed_sequence >> (pos * 2)) & 3];

			return sequence;
			}

		/*
			ENCODE_KMER_2BIT::REVERSE_COMPLEMENT()
			--------------------------------------
		*/
		/*!
			@brief Compute the reverse complement of a packed DNA sequence.  So "CAT" -> "ATG" (when bases == 3).  Done without unpacking first.
			@param kmer [in] The encoded 20-mer to reverse.
			@param bases [in] The length of the DNA sequence in bases (default=20).
			@returns The reverse complement of the sequence.
		*/
		static uint64_t reverse_complement(uint64_t kmer, uint64_t bases = 20)
			{
			uint64_t result = 0;

			/*
				Compute the complement (A<->T, C<->G), which is a bit-fit because A=00, T=11, C=01, G=10
			*/
			uint64_t complement = ~kmer;

			/*
				Reverse the order of the bases
			*/
			for (uint64_t base = 0; base < bases; base++, complement >>= 2)
				result = (result << 2) | (complement & 0x03);

			return result;
			}
	};
