/*
	FAST_BINARY_SEARCH.H
	--------------------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief The fast indexed binary search method of finding variants of k-mers in a genome.
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <stdlib.h>
#include <stdint.h>

#include <vector>
#include <iostream>

#include "job.h"
#include "finder.h"
#include "score_mit_local.h"
#include "encode_kmer_2bit.h"

/*
	CLASS FAST_BINARY_SEARCH
	------------------------
*/
/*!
	@brief The fast indexed binary search method of matching.
*/
class fast_binary_search : public finder
	{
	protected:
		static constexpr size_t index_width_in_bases = 14;			// bases to use as the key in the index
		static constexpr uint64_t in_genome_reverse_bitmap_width_in_bases = 14;
		static constexpr uint64_t in_genome_reverse_bitmap_mask = (0x1ULL << (in_genome_reverse_bitmap_width_in_bases * 2)) - 1;
		static constexpr size_t base_width_in_bits = 2;			// must be 2 (2 bits per base in the uint64_t encoding)

		/*
			Fast lookup for encoding encoded kmers, and the index too.
		*/
		std::vector<uint64_t> in_genome_bitmap; 		 // Bitmap array with 2^32 bits (i.e. 16 bases wide)
		std::vector<uint64_t> in_genome_reverse_bitmap;  // Bitmap array whose size id dependant on in_genome_reverse_bitmap_width_in_bases
		std::vector<const uint64_t *>index;

		/*
			The scoreing mechanism
		*/
		score_mit_local scorer;

	protected:
		/*
			FAST_BINARY_SEARCH::IN_GENOME_HEAD()
			------------------------------------
		*/
		/*!
			@brief Is the prefix of the k-mer anywhere in the genome being searched?
			@param key [in] The k-mer being checked.
			@returns true if the prefix is present, false otherwise.
		*/
		inline bool in_genome_head(uint64_t key) const
			{
			return in_genome_bitmap[(key >> 8) / 64] & (1ULL << ((key >> 8) % 64));
			}

		/*
			FAST_BINARY_SEARCH::IN_GENOME()
			-------------------------------
		*/
		/*!
			@brief Are both  the prefix of the k-mer and the suffix of the k-mer present in the genome?
			@param key [in] The k-mer being checked.
			@returns true if the prefix and suffix are both present, false otherwise.
		*/
		inline bool in_genome(uint64_t key) const
			{
			return
				(in_genome_bitmap[(key >> 8) / 64] & (1ULL << ((key >> 8) % 64)))
				&&
				(in_genome_reverse_bitmap[(key & in_genome_reverse_bitmap_mask) / 64] & (1ULL << ((key & in_genome_reverse_bitmap_mask) % 64)))
				;
			}

		/*
			FAST_BINARY_SEARCH::COMPUTE_INTERSECTION()
			------------------------------------------
		*/
		/*!
			@brief Compute the intersecton between a key and the genome the has already been indexed.
			@param key [in] the 20-mers being searched for.
			@param positions [out] The vector of match locations (pointers).
		*/
		double compute_intersection(uint64_t guide, uint64_t key) const;

		/*
			FAST_BINARY_SEARCH::COMPUTE_INTERSECTION_LIST()
			-----------------------------------------------
		*/
		/*!
			@brief Compute the intersecton between the keys and the genome the has already been indexed.
			@param key [in] the list of 20-mers being searched for.
			@param key_length [in] the number of k-mers encoded in key.
			@param positions [out] The vector of match locations (pointers).
		*/
		double compute_intersection_list(double threshold, const uint64_t *key, size_t key_length) const;

		/*
			FAST_BINARY_SEARCH::GENERATE_VARIATIONS()
			-----------------------------------------
		*/
		/*!
			@brief Generate all 20-mers that are 0 to 4 base mutations away from sequence (and also "likely" to be in the genome begine searched)
			@param threshold [in] early termination threshold.
			@param guide [in] the guide that all the variants are made from.
			@param sequence [in] The sequence being mutated.
			@param replacements [in] The number of mutations that have already been applied to the sequence.
			@param position [in] The position to start mautating from
			@returns The score of the sequence
		*/
		double generate_variations(double threshold, uint64_t guide, uint64_t sequence, int replacements, int position) const;

		/*
			FAST_BINARY_SEARCH::GENERATE_VARIATIONS()
			-----------------------------------------
		*/
		/*!
			@brief Generate all 20-mers that are 0 to 4 base mutations away from sequence (and also "likely" to be in the genome begine searched)
			@param threshold [in] early termination threshold.
			@param sequence [in] The sequence being mutated.
			@param variations [out] The vector of variations.
			@returns The score of the sequence
		*/
		double generate_variations(double threshold, uint64_t sequence) const
			{
			return generate_variations(threshold, sequence, sequence, 0, 0);
			}

	public:
		/*
			FAST_BINARY_SEARCH::FAST_BINARY_SEARCH()
			----------------------------------------
		*/
		/*!
			@brief Constructor
		*/
		fast_binary_search() :
			in_genome_bitmap((1ULL << 32) / 64, 0),
			in_genome_reverse_bitmap((1ULL << (in_genome_reverse_bitmap_width_in_bases * 2)) / 64, 0)
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
		virtual ~fast_binary_search()
			{
			/* Nothing */
			}


		/*
			FAST_BINARY_SEARCH::MAKE_INDEX()
			--------------------------------
			data must be sorted
		*/
		/*!
			@brief Do any necessary indexing of the genome before starting the search process.
			@param genome [in] The genome being searched.
		*/
		virtual void make_index(const std::vector<uint64_t> &genome);

		/*
			FAST_BINARY_SEARCH::PROCESS_CHUNK()
			-----------------------------------
		*/
		/*!
			@brief THREAD SAFE (if answer is not used elsewhere).  Search the genome for the given guides in test\_guides
			@param start [in] Where in the test_guides to start searching from.
			@param end [in] Where in the test_guides to start searching to.
			@param test_guides [in] The list of guides to look for - and we only look for those between start and end.
			@param answer [out] The result set
		*/
		virtual void process_chunk(job &workload)
			{
			constexpr double threshold = 0.75;												// this is the MIT Global score needed to be useful
			constexpr double threshold_sum = (100.0 / threshold) - 100.0;			// the sum of MIT Local scores must be smaller than this to to scores
			char output_buffer[50];

			uint64_t end = workload.guide.size();
			uint64_t guide_index;

			while ((guide_index = workload.get_next()) < end)
				{
				uint64_t guide = workload.guide[guide_index];
				if (workload.genome_guide_frequencies[guide_index] != 1)
					continue;
				try
					{
					double score = generate_variations(threshold_sum, guide);
					if (score != 0.0)
						{
						/*
							Convert the local score to a global score then write to the output file
						*/
						score = 100.0 / (score + 100.0);				// convert to the MIT Global Score
						encode_kmer_2bit::unpack_20mer(output_buffer, guide);
						output_buffer[20] = ' ';
						int bytes = snprintf(&output_buffer[21], &output_buffer[49] - &output_buffer[21], "%2.2f\n", score);
//						auto [string_end, error] = std::to_chars(&output_buffer[21], &output_buffer[49], (int)(score * 100.0));
//						*string_end++ = '\n';
//						*string_end = '\0';

						workload.file_mutex.lock();
//							fwrite(output_buffer, sizeof(char), string_end - output_buffer, workload.output_file);
							fwrite(output_buffer, sizeof(char), bytes + 21, workload.output_file);
						workload.file_mutex.unlock();

						/*
							Update the stats
						*/
						workload.hits++;
						if (score > workload.best_score)
							{
							workload.stats_mutex.lock();
								workload.best_20mer = guide;
								workload.best_score = score;
							workload.stats_mutex.unlock();
							}
						}
					}
				catch (...)
					{
					/*
						This happens if the recursion that generates the variants early terminates
						in which case we do nothing (because there is no work to) other than move
						on to the next guide in the list.
					*/
					}
				}
			}
	};
