/*
	SCORE_MIT_LOCAL.H
	-----------------
	Copyright (c) 2023 Andrew Trotman and Dimitri Perrin.

	The code to compute the MIT score was supplied by Dimitri Perrin.
	For an explination of the matric see: https://dx.doi.org/10.1038/nbt.2647
*/
#pragma once
/*!
	@file
	@brief Compute the MIT score for the difference between k-mers (up to 20-mer) encoded in 2-bits per base.  See https://dx.doi.org/10.1038/nbt.2647
	@author Andrew Trotman,  Dimitri Perrin
	@copyright 2023 Andrew Trotman, Dimitri Perrin
*/
#include <vector>

/*
	CLASS SCORE_MIT_LOCAL
	---------------------
*/
/*!
	@brief Compute MIT local scores
*/
class score_mit_local
	{
	private:
		/*
			MITLOCAL_SCORE
			--------------
			Array of MIT local scores for all possible variants of a k-mer up-to length 20.
		*/
		static double MITLocal_score[1024 * 1024];

		/*
			SCORE_MIT_LOCAL::DIFFERENCES()
			------------------------------
		*/
		/*!
			@brief compute a bitstring of differences between the guide and the variant.
			@param guide [in] the guide.
			@param variant [in] the variant.
			@returns Bitstring representing where differences were found, encoded: 9 19 8 18 7 17 6 16 5 15 4 14 3 13 2 12 1 11 0 10
		*/
		static uint64_t differences(uint64_t guide, uint64_t variant)
			{
			/*
				Compute the "pairs" where there is a difference
			*/
			uint64_t xor_set = guide ^ variant;

			/*
				Take the high bits and the low bits and OR them resulting in only the high bit of the
				pair being set if either bit is set (i.e. there is a change).
			*/
			uint64_t diffs = (((xor_set & 0x5555555555555555) << 1) | (xor_set & 0xAAAAAAAAAAAAAAAA));

			/*
				Take the low 20 bits and fold into it the top 20 bits.  The resultant bitstring will be in the order:
				9 19 8 18 7 17 6 16 5 15 4 14 3 13 2 12 1 11 0 10
			*/
			uint64_t folded = (diffs & 0xFFFFF) | (diffs >> 21);

			return folded;
			}

		/*
			SCORE_MIT_LOCAL::UNFOLD()
			-------------------------
		*/
		/*!
			@brief Turn the folded represenation of the differeces into a bit-encoding of the differences.  That is, put the difference list in order.  So bit 0 represents a difference in the 0th base, etc.
			@param folded [in] The representatio returned by differences().
			@returns Bitstring where bit 0 represents a change in the 0th base, 1 in the first, and so on.
		*/
		static uint64_t unfold(uint64_t folded);

		/*
			SCORE_MIT_LOCAL::LIST_DIFFERENCES()
			-----------------------------------
		*/
		/*!
			@brief Given a bitwise encoded of differences (from differences()), return a vector of integers where each integer is an index of a difference (coutint from 0 of the right hand end).
			@param answer [out] the list of k-mer index positions that are different.
			@param differences [in] the folded bit string of differences.
			@returns A std::vector<uin64_t> if difference locations.
		*/
		static void list_differences(std::vector<uint64_t> &answer, uint64_t differences);

		/*
			SCORE_MIT_LOCAL::CALCULATE_MITLOCAL_SCORE()
			-------------------------------------------
			Supplied by Dimitri Perrin
		*/
		/*!
			@brief Calculate the local MIT score based on the positions of mismatches (see: https://dx.doi.org/10.1038/nbt.2647)
			@param mismatch_array [in]  An int array that contains the position of mismatches
			@param length [in]  The length of the param `mismatch_array`
			@returns The local MIT score
		*/
		static double calculate_MITLocal_score(uint64_t *mismatch_array, uint64_t length);


	public:
		/*
			SCORE_MIT_LOCAL::SCORE_MIT_LOCAL()
			----------------------------------
		*/
		/*!
			@brief Constructor
		*/
		score_mit_local();

		/*
			SCORE_MIT_LOCAL::SCORE()
			------------------------
		*/
		/*!
			@brief Compute the MIT local score for variant with respect to guide.
			@param guide [in] The guide to compare to.
			@param variant [in] The sequence to score with respect to the guide.
			@returns The MIT local score.
		*/
		double score(uint64_t guide, uint64_t variant) const
			{
			return MITLocal_score[differences(guide, variant)];
			}
	};
