/*
	SEQUENCE_SCORE_PAIR.H
	---------------------
	Copyright (c) 2023 Andrew Trotman
*/

#pragma once
/*!
	@file
	@brief Class to save a sequence with a score
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <stdint.h>

/*
	CLASS SEQUENCE_SCORE_PAIR
	-------------------------
	A given 20-mer has a given score
*/
/*!
	@brief Holder for a 20-mer and a score for that 20-mer
*/
class sequence_score_pair
	{
	public:
		uint64_t sequence;
		double score;

	public:
		/*
			SEQUENCE_SCORE_PAIR::SEQUENCE_SCORE_PAIR()
			------------------------------------------
		*/
		/*!
			@brief constructor
		*/
		sequence_score_pair(uint64_t sequence, double score) :
			sequence(sequence),
			score(score)
			{
			/*
				Nothing
			*/
			}
	};


