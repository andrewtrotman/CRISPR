/*
	HAMMING_DISTANCE.H
	------------------
	Copyright (c) 2023 Andrew Trotman
*/
/*!
	@file
	@brief The hamming distance method of finding variants of k-mers in a genome.
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/
#pragma once

#include "score_mit_local.h"

#if defined(__APPLE__) || defined(__linux__)
	#define __popcnt64 __builtin_popcountll
#endif

/*
	CLASS HAMMING_DISTANCE
	----------------------
*/
/*!
	@brief The hamming distance method of finding variants of k-mers in a genome.
*/
class hamming_distance : public finder
	{
	private:
		score_mit_local scorer;

	private:
		/*
			HAMMING_DISTANCE::DISTANCE()
			----------------------------
		*/
		/*!
			@brief Compute the  hamming distance between two integers.
			@param a [in] The first integer.
			@param b [in] The second integer.
			@returns The number of bit positions in which they differ
		*/
		inline uint64_t distance(uint64_t a, uint64_t b)
			{
			return __popcnt64(a ^ b);
			}

		/*
			HAMMING_DISTANCE::COMPUTE_HAMMING_SET()
			---------------------------------------
		*/
		/*!
			@brief Run through the genome looking for matches.
			@param threshold [in] early termination threshold.
			@param twice_the_max_distance [in] Twice the hamming distance looked for (will find all matches where the hamming distance <= twice\_the\_max\_distance / 2).
			@param key [in] The kmer to look for.
			@param workload [in] the queries, data, and frequencies
			@return the local score for the key
		*/
		double compute_hamming_set(double threshold, uint64_t twice_the_max_distance, uint64_t key, job &workload)
			{
			double score = 0;
			const uint64_t *end = workload.genome_guides.data() + workload.genome_guides.size();

			for (const uint64_t *which = workload.genome_guides.data(); which < end; which++)
				if (distance(key, *which) <= twice_the_max_distance)
					{
					double frequency = workload.genome_guide_frequencies[which - &workload.genome_guides[0]];
					score += scorer.score(key, *which) * frequency;
					if (score > threshold)
						return 0.0;
					}

			return score;
			}

	public:
		/*
			HAMMING_DISTANCE::~HAMMING_DISTANCE()
			-------------------------------------
		*/
		virtual ~hamming_distance()
			{
			/* Nothing */
			}

		/*
			HAMMING_DISTANCE::MAKE_INDEX()
			------------------------------
		*/
		/*!
			@brief Do any necessary indexing of the genome before starting the search process (for hamming distance, this does nothing) .
			@param genome [in] The genome being searched.
		*/
		virtual void make_index(const std::vector<uint64_t> &genome)
			{
			/* Nothing */
			}

		/*
			HAMMING_DISTANCE::PROCESS_CHUNK()
			---------------------------------
		*/
		/*!
			@brief THREAD SAFE.  Search the genome for the given guides in test\_guides
			@param workload [in] the queries, data, and frequencies
			@param answer [out] The result set
		*/
		virtual void process_chunk(job &workload)
			{
			uint64_t guide_index;
			uint64_t end = workload.guide.size();
			while ((guide_index = workload.get_next()) < end)
				{
				if (workload.guide_frequencies[guide_index] != 1)
					continue;

				uint64_t guide = workload.guide[guide_index];
				double score = compute_hamming_set(0.75, 8, guide, workload);
				/*
					FIX:
						Write this out to disk using the same code as used in the fast_binary_search version
				*/
				}
			}
	};
