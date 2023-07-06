/*
	FAST_BINARY_SEARCH_AVX512.H
	---------------------------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief The fast indexed binary search method of finding variants of k-mers in a genome.  But uses AVX512
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <stdlib.h>
#include <stdint.h>
#include <immintrin.h>

#include <vector>
#include <iostream>

#include "fast_binary_search.h"

#define USE_AVX512F 1

/*
	CLASS FAST_BINARY_SEARCH_AVX512
	-------------------------------
*/
/*!
	@brief The fast indexed binary search method of finding variants of k-mers in a genome.  Uses AVX512 instructions
*/
class fast_binary_search_avx512 : public fast_binary_search
	{
	private:
		/*
			AVX_BINARY_SEARCH()
			-------------------
		*/
		/*!
			@brief 8 concurrent binary searches for uint64\_t values using AVX512.
			@param lower [in] pointers to the start of the array (inclusive).
			@param upper [in] pointers to the end of the array (exclusive).
			@param key [in] the 8 values to look for.
			@returns 8 pointers, just like std::lower\_bound().
		*/
		static __m512i avx_binary_search(__m512i lower, __m512i upper, __m512i key);

		/*
			FAST_BINARY_SEARCH::AVX_COMPUTE_INTERSECTION_LIST()
			---------------------------------------------------
		*/
		/*!
			@brief Compute the intersecton between the keys and the genome the has already been indexed.
			@param key [in] the list of 20-mers being seaerched for.
			@param key_length [in] the number of k-mers encoded in key.
			@param positions [out] The vector of match locations (pointers).
		*/
		void avx_compute_intersection_list(const uint64_t *key, size_t key_length, std::vector<sequence_score_pair> &positions) const;

	public:
		/*
			FAST_BINARY_SEARCH::FAST_BINARY_SEARCH()
			----------------------------------------
		*/
		/*!
			@brief Constructor
		*/
		fast_binary_search_avx512() :
			fast_binary_search()
			{
			/*
				Nothing
			*/
			}

		/*
			FAST_BINARY_SEARCH::~FAST_BINARY_SEARCH()
			-----------------------------------------
		*/
		/*!
			@brief Destructor
		*/
		~fast_binary_search_avx512()
			{
			/* Nothing */
			}

		/*
			FAST_BINARY_SEARCH::PROCESS_CHUNK()
			-----------------------------------
		*/
		/*!
			@brief THREAD SAFE.  Search the genome for the given guides in test\_guides
			@param start [in] Where in the test_guides to start searching from.
			@param end [in] Where in the test_guides to start searching to.
			@param test_guides [in] The list of guides to look for - and we only look for those between start and end.
			@param packed_genome_guides [in] The genome to look in.
			@param answer [out] The result set
		*/
		virtual void process_chunk(size_t start, size_t end, std::vector<uint64_t> &test_guides, std::vector<uint64_t> &packed_genome_guides, std::vector<std::vector<sequence_score_pair>> &answer)
			{
			std::vector<uint64_t> variations;

			for (size_t which = start; which < end; which++)
				{
				variations.clear();
				generate_variations(test_guides[which], variations);
				avx_compute_intersection_list(variations.data(), variations.size(), answer[which]);
				}
			}
	};

