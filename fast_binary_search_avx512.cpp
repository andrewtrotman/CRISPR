/*
	FAST_BINARY_SEARCH_AVX512.CPP
	-----------------------------
	Copyright (c) 2023 Andrew Trotman
*/

#include "fast_binary_search_avx512.h"

#if defined(__APPLE__) || defined(__linux__)
	#define __popcnt16 __builtin_popcount
#endif


#ifdef USE_AVX512F
	/*
		FAST_BINARY_SEARCH_AVX512::AVX_BINARY_SEARCH()
		----------------------------------------------
	*/
	inline __m512i fast_binary_search_avx512::avx_binary_search(__m512i lower, __m512i upper, __m512i key)
		{
		__m512i one = _mm512_set1_epi64(1);
		lower = _mm512_srli_epi64(lower, 3);					// Convert from pointers to indexes
		lower = _mm512_sub_epi64(lower, one);					// Necessary because the binary search counts from 1 (not from 0)
		upper = _mm512_srli_epi64(upper, 3);					// Convert from pointers to indexes
		__mmask8 not_finished;

		/*
			Binary search
		*/
		while ((not_finished = _mm512_cmpneq_epi64_mask(_mm512_add_epi64(lower, one), upper)) != 0x00)
			{
			__m512i middle = _mm512_srli_epi64(_mm512_add_epi64(upper, lower), 1);
			__mmask8 results = _mm512_cmplt_epi64_mask(_mm512_i64gather_epi64(middle, nullptr, sizeof(uint64_t)), key);
			lower = _mm512_mask_blend_epi64 (not_finished & results, lower, middle);
			upper = _mm512_mask_blend_epi64 (not_finished & ~results, upper, middle);
			}

		return _mm512_slli_epi64(upper, 3);		// convert back into pointers
		}

	/*
		FAST_BINARY_SEARCH_AVX512::COMPUTE_INTERSECTION_LIST()
		------------------------------------------------------
	*/
	void fast_binary_search_avx512::avx_compute_intersection_list(const uint64_t *key, size_t key_length, std::vector<sequence_score_pair> &positions) const
		{
		positions.push_back(sequence_score_pair(*key, 1));				// push the guide before working on the variants

		const uint64_t *key_end = key + key_length;
		const uint64_t *current_key = key + 1;

		/*
			Do the AVX512 stuff first
		*/
		__m512i one = _mm512_set1_epi64(1);

		while (current_key + 8 < key_end)
			{
			/*
				Get the current key set (that we're looking for)
			*/
			__m512i key_set = _mm512_loadu_epi64(current_key);
			__m512i index_key_set = _mm512_srli_epi64(key_set, (20 - index_width_in_bases) * base_width_in_bits);

			/*
				Compute the start and end of the ranges
			*/
			__m512i start_set = _mm512_i64gather_epi64(index_key_set, &index[0], sizeof(uint64_t *));
			__m512i end_set = _mm512_i64gather_epi64(_mm512_add_epi64(index_key_set, one), &index[0], sizeof(uint64_t *));

			/*
				Do the AVX512-parallel binary search given the pointers
			*/
			__m512i found_pointers = avx_binary_search(start_set, end_set, key_set);

			/*
				At this point found_set is a set of pointers to the data - which may or may not match the key (same as the result of std::lower_bound())
			*/
			__m512i found_values = _mm512_i64gather_epi64(found_pointers, nullptr, 1);		// 1 because found_set is a set of pointers to 64-bit integers
			__mmask8 found_masks = _mm512_cmpeq_epi64_mask(key_set, found_values);

			if (found_masks != 0)
				{
				/*
					I really want to write directly into the positions vector, but it isn't clear how large it needs to be pre-allocated
				*/
				alignas(64) uint64_t *position[8];
				_mm512_mask_compressstoreu_epi64 (&position[0], found_masks, found_pointers);

				switch(__popcnt16(found_masks))
					{
					case 8:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[7], 1));			// FALL THROUGH
					case 7:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[6], 1));			// FALL THROUGH
					case 6:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[5], 1));			// FALL THROUGH
					case 5:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[4], 1));			// FALL THROUGH
					case 4:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[3], 1));			// FALL THROUGH
					case 3:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[2], 1));			// FALL THROUGH
					case 2:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[1], 1));			// FALL THROUGH
					case 1:
						positions.push_back(sequence_score_pair(*(uint64_t *)position[0], 1));			// FALL THROUGH
					}
				}
			current_key += 8;
			}

		/*
			Then fall through to process the last few one at a time (if no AVX512 then do this for all guides)
		*/
		while (current_key < key_end)
			{
			size_t index_key = *current_key >> ((20 - index_width_in_bases) * base_width_in_bits);

			const uint64_t *found = std::lower_bound(index[index_key], index[index_key + 1], *current_key);
			if (*found == *current_key)
				positions.push_back(sequence_score_pair(*found, 1));
			current_key++;
			}
		}
#else

	/*
		FAST_BINARY_SEARCH_AVX512::AVX_BINARY_SEARCH()
		----------------------------------------------
	*/
	__m512i fast_binary_search_avx512::avx_binary_search(__m512i lower, __m512i upper, __m512i key)
		{
		std::cout << "AVX512 not enabled in this build\n";
		exit(0);
		}

	/*
		FAST_BINARY_SEARCH_AVX512::AVX_COMPUTE_INTERSECTION_LIST()
		----------------------------------------------------------
	*/
	void fast_binary_search_avx512::avx_compute_intersection_list(const uint64_t *key, size_t key_length, std::vector<sequence_score_pair> &positions) const
		{
		std::cout << "AVX512 not enabled in this build\n";
		exit(0);
		}

#endif
