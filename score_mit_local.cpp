/*
	SCORE_MIT_LOCAL.CPP
	-------------------
	Copyright (c) 2023 Andrew Trotman

	The code to compute the MIT score was supplied by Dimitri Perrin.
	For an explination of the matric see: https://dx.doi.org/10.1038/nbt.2647
*/
#include <bitset>
#include <vector>
#include <iostream>
#include <algorithm>

#include "score_mit_local.h"

/*
	MITLOCAL_SCORE
	--------------
	Array of MIT local scores for all possible variants of a k-mer up-to length 20.
*/
double score_mit_local::MITLocal_score[1024 * 1024];

/*
	SCORE_MIT_LOCAL::SCORE_MIT_LOCAL()
	----------------------------------
	Precompute all possible scores and store them in an array indexed by the differance bit-string (which is folded).
*/
score_mit_local::score_mit_local()
	{
	std::vector<uint64_t> positions(20);

	for (uint64_t difflist = 1; difflist < 1024 * 1024; difflist++)
		{
		positions.clear();
		list_differences(positions, difflist);
		MITLocal_score[difflist] = calculate_MITLocal_score(positions.data(), positions.size());
		}
	}

/*
	SCORE_MIT_LOCAL::LIST_DIFFERENCES()
	-----------------------------------
*/
void score_mit_local::list_differences(std::vector<uint64_t> &answer, uint64_t differences)
	{
	uint64_t bit_position[] = {0x000002, 0x000008, 0x000020, 0x000080, 0x000200, 0x000800, 0x002000, 0x008000, 0x020000, 0x080000, 0x000001, 0x000004, 0x000010, 0x000040, 0x000100, 0x000400, 0x001000, 0x004000, 0x010000, 0x040000};

	for (uint64_t bit = 0; bit < 20; bit++)
		if (differences & bit_position[bit])
			answer.push_back(bit);
	}

/*
	SCORE_MIT_LOCAL::UNFOLD()
	-------------------------
	The input bit positions are:
	9 19 8 18 7 17 6 16 5 15 4 14 3 13 2 12 1 11 0 10
*/
uint64_t score_mit_local::unfold(uint64_t folded)
	{
	uint64_t bit_position[] = {0x000002, 0x000008, 0x000020, 0x000080, 0x000200, 0x000800, 0x002000, 0x008000, 0x020000, 0x080000, 0x000001, 0x000004, 0x000010, 0x000040, 0x000100, 0x000400, 0x001000, 0x004000, 0x010000, 0x040000};
	uint64_t answer = 0;

	for (uint64_t bit = 0; bit < 20; bit++)
		if (folded & bit_position[bit])
			answer |= 1 << bit;

	return answer;
	}

/*
	SCORE_MIT_LOCAL::CALCULATE_MITLOCAL_SCORE()
	-------------------------------------------
	Supplied by Dimitri Perrin

	This function calculates the local MIT score based on the positions of mismatches.
	https://dx.doi.org/10.1038/nbt.2647

	@param mismatch_array, An int array that contains the position of mismatches
	@param length, The length of the param `mismatch_array`
	@return The local MIT score
*/
double score_mit_local::calculate_MITLocal_score(uint64_t *mismatch_array, uint64_t length)
	{
	uint64_t i;
	double T1 = 1.0, T2, T3, d = 0.0, score;
	/* Mismatch penalty array */
	double M[] = {0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583};

	/* 1st term */
	for (i = 0; i < length; ++i)
		T1 = T1 * (1.0 - M[mismatch_array[i]]);

	/* 2nd term */
	if (length == 1)
		d = 19.0;
	else
		{
		for (i = 0; i < length - 1; ++i)
			d += mismatch_array[i + 1] - mismatch_array[i];
		d = d / (length - 1);
		}
	T2 = 1.0 / ((19.0 - d) / 19.0 * 4.0 + 1);

	/* 3rd term */
	T3 = 1.0 / (length * length);

	/* Local score */
	score = T1 * T2 * T3 * 100;

	return score;
	}

