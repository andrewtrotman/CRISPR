/*
	FAST_BINARY_SEARCH.CPP
	----------------------
	Copyright (c) 2023 Andrew Trotman
*/
#include <math.h>

#include "forceinline.h"
#include "score_mit_local.h"
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
inline double fast_binary_search::compute_intersection(uint64_t guide, uint64_t key) const
	{
	const uint64_t *found;
	size_t index_key = key >> ((20 - index_width_in_bases) * base_width_in_bits);
	size_t mers_to_search = index[index_key + 1] - index[index_key];
	if (mers_to_search < 13)				// On my Mac, 13 was the cross-over point between binary and linear search
		for (found = index[index_key]; *found < key; found++) {/* Nothing */ }// If the length of the list is short then linear search in an ordered list
	else
		found = std::lower_bound(index[index_key], index[index_key + 1], key); // If the length of the list is long then binary search in an ordered list

	if (*found == key)
		{
		double frequency = workload->genome_guide_frequencies[found - &workload->genome_guides[0]];
		return scorer.score(guide, key) * frequency;			// compute the score of the variant.
		}

	return 0.0;		// no match found so there is no effect on the cumulative score.
	}

/*
	FAST_BINARY_SEARCH::COMPUTE_INTERSECTION_LIST()
	-----------------------------------------------
*/
double fast_binary_search::compute_intersection_list(double threshold, const uint64_t *key, size_t key_length) const
	{
	const uint64_t guide	= *key;
	const uint64_t *key_end = key + key_length;
	double score = 0.0;

	for (const uint64_t *current_key = key + 1; current_key < key_end; current_key++)
		{
		score += compute_intersection(guide, *current_key);
		if (score > threshold)
			return 0.0;
		}

	return score;
	}

/*
	FAST_BINARY_SEARCH::GENERATE_VARIATIONS()
	-----------------------------------------
*/
double fast_binary_search::generate_variations(double threshold, uint64_t guide, uint64_t sequence, int replacements, int position) const
	{
	double score = 0;

	if (in_genome(sequence))
		score = compute_intersection(guide, sequence);

	if (replacements == 4)
		return score;

	for (uint64_t i = position; i < 20; i++)
		{
		/*
			If the first 16 bases are not in the genome then there is no point in making
			variants in the last 4 positions (because they can't be in the genone either).
		*/
		if (i == 16 && !in_genome_head(sequence))
			return score;

		/*
			This optimisation is from Timothy Chappell who points out that by XORing the
			current base with 01, 02, and 03 you automatically get all the variants except the
			original sequence (but not necessarily in a predictable order).  So we don't need to
			mask and replace, we simply XOR with the original sequence and call outselves recursively.
		*/
		uint64_t shifter = (19 - i) * 2;
		score += generate_variations(threshold, guide, sequence ^ (1ULL << shifter), replacements + 1, i + 1);
		if (score > threshold)
			throw 1;
		score += generate_variations(threshold, guide, sequence ^ (2ULL << shifter), replacements + 1, i + 1);
		if (score > threshold)
			throw 1;
		score += generate_variations(threshold, guide, sequence ^ (3ULL << shifter), replacements + 1, i + 1);
		if (score > threshold)
			throw 1;
		}

	return score;
	}
