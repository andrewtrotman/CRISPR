/*
	FAST_BINARY_SEARCH.CPP
	----------------------
	Copyright (c) 2023 Andrew Trotman
*/
#include <math.h>

#include "forceinline.h"
#include "fast_binary_search.h"

/*
	FAST_BINARY_SEARCH::MAKE_INDEX()
	--------------------------------
	data must be sorted
*/
void fast_binary_search::make_index(const std::vector<uint64_t> &integers)
	{
	/*
		Allocate space for the index
	*/
	index.resize((size_t)pow(4, index_width_in_bases) + 1);		// +1 because we need a terminator on the end

	/*
		Get the end of each group by passing over the data once.  But first make sure the index starts at the beginning and finishes at the end.
	*/
	index[0] = integers.data();
	index[index.size() - 1] = integers.data() + integers.size() - 1;
	const uint64_t *end = integers.data() + integers.size();
	for (const uint64_t *current = integers.data(); current < end - 1; current++)
		{
		index[(*current >> ((20 - index_width_in_bases) * base_width_in_bits)) + 1] = current + 1;

		uint64_t index = *current >> 8;
		in_genome_bitmap[index / 64] |= (1ULL << (index % 64));  // Set the respective bit to 1

		index = *current & in_genome_reverse_bitmap_mask;
		in_genome_reverse_bitmap[index / 64] |= (1ULL << (index % 64));  // Set the respective bit to 1
		}

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

/*
	FAST_BINARY_SEARCH::COMPUTE_INTERSECTION()
	------------------------------------------
*/
inline void fast_binary_search::compute_intersection(const uint64_t key, std::vector<const uint64_t *> &positions) const
	{
		size_t index_key = key >> ((20 - index_width_in_bases) * base_width_in_bits);
		size_t mers_to_search = index[index_key + 1] - index[index_key];
		if (mers_to_search < 13)				// On my Mac, 13 was the cross-over point between binary and linear search
			{
			/*
				If the length of the list is short then
				Linear search in an ordered list
			*/
			const uint64_t *found;
			for (found = index[index_key]; *found < key; found++)
				{
				/* Nothing */
				}
			if (*found == key)
				positions.push_back(found);
			}
		else
			{
			/*
				If the length of the list is long then
				Binary search in an ordered list
			*/
			const uint64_t *found = std::lower_bound(index[index_key], index[index_key + 1], key);
			if (*found == key)
				positions.push_back(found);
			}
	}

/*
	FAST_BINARY_SEARCH::COMPUTE_INTERSECTION_LIST()
	-----------------------------------------------
*/
void fast_binary_search::compute_intersection_list(const uint64_t *key, size_t key_length, std::vector<const uint64_t *> &positions) const
	{
	positions.push_back(nullptr);				// push the guide before working on the variants

	const uint64_t *key_end = key + key_length;
	const uint64_t *current_key = key + 1;

	while (current_key < key_end)
		{
		compute_intersection(*current_key, positions);
		current_key++;
		}
	}

/*
	FAST_BINARY_SEARCH::GENERATE_VARIATIONS()
	-----------------------------------------
*/
void fast_binary_search::generate_variations(uint64_t sequence, std::vector<uint64_t> &variations, int replacements, int position) const
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
		generate_variations(sequence ^ (1ULL << shifter), variations, replacements + 1, i + 1);
		generate_variations(sequence ^ (2ULL << shifter), variations, replacements + 1, i + 1);
		generate_variations(sequence ^ (3ULL << shifter), variations, replacements + 1, i + 1);
		}
	}
