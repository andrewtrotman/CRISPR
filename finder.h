
/*
	FINDER.H
	--------
	Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
	@file
	@brief Pure virtual base class for differnt finding techniques
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <vector>

#include "job.h"
#include "encode_kmer_2bit.h"
#include "encode_kmer_3bit.h"
#include "sequence_score_pair.h"

/*
	CLASS FINDER
	------------
*/
/*!
	@brief Pure virtual base class for differnt finding techniques
*/
class finder
	{
	public:
		/*
			FINDER::~FINDER()
			-----------------
		*/
		virtual ~finder()
			{
			/* Nothing */
			}
			
		/*
			FINDER::MAKE_INDEX()
			--------------------
		*/
			/*!
				@brief Do any necessary indexing of the genome before starting the search process (for hamming distance, this does nothing) .
				@param genome [in] The genome being searched.
			*/

		virtual void make_index(const std::vector<uint64_t> &genome) = 0;

		/*
			FINDER::PROCESS_CHUNK()
			-----------------------
		*/
		/*!
			@brief THREAD SAFE.  Search the genome for the given guides in test\_guides
			@param workload [in] the queries, data, and frequencies
			@param answer [out] The result set
		*/
		virtual void process_chunk(job &workload) = 0;



		/*
			FINDER::SAVE_RESULT()
			---------------------
		*/
		/*!
			@brief THREAD SAFE save the result to the disk.
			@param workload [in] The job - including the file handle and the ability to lock it.
			@param kmer [in] The 20-mer to write to the file.
			@param bits_per_guide [in] is this encoding 2 bits per guide or 3 bits per base (fast_binary_search or hamming)
			@param score [in] The score associated with the kmer.
		*/
		void save_result(job &workload, uint64_t guide, size_t bits_per_guide, double score)
			{
			char output_buffer[50];
			/*
				Convert the local score to a global score then write to the output file
			*/
			if (bits_per_guide == 2)
				encode_kmer_2bit::unpack_20mer(output_buffer, guide);
			else
				encode_kmer_3bit::unpack_20mer(output_buffer, guide);

			output_buffer[20] = ' ';
			int bytes = snprintf(&output_buffer[21], &output_buffer[49] - &output_buffer[21], "%2.2f\n", score);
//				auto [string_end, error] = std::to_chars(&output_buffer[21], &output_buffer[49], (int)(score * 100.0));
//				*string_end++ = '\n';
//				*string_end = '\0';

			workload.file_mutex.lock();
//					fwrite(output_buffer, sizeof(char), string_end - output_buffer, workload.output_file);
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
	};
