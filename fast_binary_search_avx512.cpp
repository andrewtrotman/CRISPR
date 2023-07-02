/*
	FAST_BINARY_SEARCH.CPP
	----------------------
	Copyright (c) 2023 Andrew Trotman
*/



/*
	Define USE_AVX512F to use the AVX512 binary search - which on old AVX512 machines is slightly slower than using scalar operations
*/
#undef USE_AVX512F

namespace fast_binary_search
	{
	/*
		The number of bases to use as the key to the index.
	*/
	constexpr size_t index_width_in_bases = 14;
	constexpr uint64_t in_genome_reverse_bitmap_width_in_bases = 14;
	constexpr uint64_t in_genome_reverse_bitmap_mask = (0x1ULL << (in_genome_reverse_bitmap_width_in_bases * 2)) - 1;
//	constexpr size_t index_width_in_bases = 1;
	constexpr size_t base_width_in_bits = 2;

	/*
		Fast lookup for encoding encoded kmers
	*/
	static std::vector<uint64_t> in_genome_bitmap((1ULL << 32) / 64, 0);  // Bitmap array with 2^32 bits
	static std::vector<uint64_t> in_genome_reverse_bitmap((1ULL << (in_genome_reverse_bitmap_width_in_bases * 2)) / 64, 0);  // Bitmap array

	/*
		The index
	*/
	std::vector<const uint64_t *>index;

	/*
		COMPUTE_INDEX()
		---------------
		data must be sorted
	*/
	void compute_index(const uint64_t *data, size_t length)
		{
		/*
			Allocate space for the index
		*/
		index.resize((size_t)pow(4, index_width_in_bases) + 1);		// +1 because we need a terminator on the end

		/*
			Get the end of each group by passing over the data once.  But first make sure the index starts at the beginning and finishes at the end.
		*/
		index[0] = data;
		index[index.size() - 1] = data + length - 1;
		const uint64_t *end = data + length;
		for (const uint64_t *current = data; current < end - 1; current++)
			index[(*current >> ((20 - index_width_in_bases) * base_width_in_bits)) + 1] = current + 1;

		/*
			Remove nullptrs by setting the start of the current range to the start of the previous range.
		*/
		const uint64_t **index_end = &index[index.size() - 1];
		const uint64_t **previous = &index[0];

		for (const uint64_t **current = &index[0]; current < index_end; current++)
			{
			if (*current == nullptr)
				*current = *previous;
			previous = current;
			}
		}


#ifdef USE_AVX512F
	/*
		AVX_BINARY_SEARCH()
		-------------------
	*/
	__m512i avx_binary_search(__m512i lower, __m512i upper, __m512i key)
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
#endif

	/*
		IN_GENOME_HEAD()
		----------------
	*/
	bool in_genome_head(uint64_t key)
		{
		return in_genome_bitmap[(key >> 8) / 64] & (1ULL << ((key >> 8) % 64));
		}

	/*
		IN_GENOME()
		-----------
	*/
	bool in_genome(uint64_t key)
		{
		return
			(in_genome_bitmap[(key >> 8) / 64] & (1ULL << ((key >> 8) % 64)))
			&&
			(in_genome_reverse_bitmap[(key & in_genome_reverse_bitmap_mask) / 64] & (1ULL << ((key & in_genome_reverse_bitmap_mask) % 64)))
			;
		}

    /*
       BUILD_IN_GENOME_BITMAP()
       ------------------------
     */
	void build_in_genome_bitmap(const std::vector<uint64_t> &integers)
		{
		for (uint64_t num : integers)
			{
			uint64_t index = num >> 8;
			in_genome_bitmap[index / 64] |= (1ULL << (index % 64));  // Set the respective bit to 1

			index = num & in_genome_reverse_bitmap_mask;
			in_genome_reverse_bitmap[index / 64] |= (1ULL << (index % 64));  // Set the respective bit to 1
			}
		}

	/*
		COMPUTE_INTERSECTION_LIST()
		---------------------------
	*/
	void compute_intersection_list(const uint64_t *key, size_t key_length, const uint64_t *data, size_t data_length, std::vector<const uint64_t *> &positions)
		{
		positions.push_back(nullptr);				// push the guide before working on the variants

		const uint64_t *key_end = key + key_length;
		const uint64_t *current_key = key + 1;

		/*
			Do the AVX512 stuff first
		*/
#ifdef USE_AVX512F
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
						positions.push_back(position[7]);			// FALL THROUGH
					case 7:
						positions.push_back(position[6]);			// FALL THROUGH
					case 6:
						positions.push_back(position[5]);			// FALL THROUGH
					case 5:
						positions.push_back(position[4]);			// FALL THROUGH
					case 4:
						positions.push_back(position[3]);			// FALL THROUGH
					case 3:
						positions.push_back(position[2]);			// FALL THROUGH
					case 2:
						positions.push_back(position[1]);			// FALL THROUGH
					case 1:
						positions.push_back(position[0]);			// FALL THROUGH
					}
				}
			current_key += 8;
			}
#endif

		/*
			Then fall through to process the last few one at a time (if no AVX512 then do this for all guides)
		*/
		while (current_key < key_end)
			{
			size_t index_key = *current_key >> ((20 - index_width_in_bases) * base_width_in_bits);
			size_t mers_to_search = index[index_key + 1] - index[index_key];
			if (mers_to_search < 13)				// On my Mac, 13 was the cross-over point between binary and linear search
				{
				/*
					If the length of the list is short then
					Linear search in an ordered list
				*/
				const uint64_t *found;
				for (found = index[index_key]; *found < *current_key; found++)
					{
					/* Nothing */
					}
				if (*found == *current_key)
					positions.push_back(found);
				}
			else
				{
				/*
					If the length of the list is long then
					Binary search in an ordered list
				*/
				const uint64_t *found = std::lower_bound(index[index_key], index[index_key + 1], *current_key);
				if (*found == *current_key)
					positions.push_back(found);
				}
			current_key++;
			}
		}

	/*
		GENERATE_VARIATIONS_BINARY()
		----------------------------
	*/
	void generate_variations_binary(uint64_t sequence, std::vector<uint64_t> &variations, int replacements = 0, int position = 0)
		{
		if (in_genome(sequence))
			variations.push_back(sequence);
		if (replacements == 4)
			return;
		for (uint64_t i = position; i < 20; i++)
			{
			/*
				If the first 16 bases are not in the genome then there is no point in making
				variants in the last 4 positions (because they can't be in the genone either).
			*/
			if (i == 16 && !in_genome_head(sequence))
				return;

			/*
				This optimisation is from Timothy Chappell who points out that by XORing the
				current base with 01, 02, and 03 you automatically get all the variants except the
				original sequence (but not necessarily in a predictable order).  So we don't need to
				mask and replace, we simply XOR with the original sequence and call outselves recursively.
			*/
			uint64_t shifter = (19 - i) * 2;
			generate_variations_binary(sequence ^ (1ULL << shifter), variations, replacements + 1, i + 1);
			generate_variations_binary(sequence ^ (2ULL << shifter), variations, replacements + 1, i + 1);
			generate_variations_binary(sequence ^ (3ULL << shifter), variations, replacements + 1, i + 1);
			}
		}

	/*
		GENERATE_VARIATIONS()
		---------------------
	*/
	void generate_variations(uint64_t sequence, std::vector<uint64_t> &variations, int replacements = 0, int position = 0)
		{
		generate_variations_binary(sequence, variations, replacements, position);
		}
}

