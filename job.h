/*
	JOB.H
	-----
*/
#pragma once
/*!
	@file
	@brief Class representing the entire workload.
	@author Andrew Trotman
	@copyright 2023 Andrew Trotman
*/

#include <mutex>
#include <atomic>
#include <vector>

/*
	CLASS JOB()
	-----------
*/
/*!
	@brief Class representing the entire workload.
*/
class job
	{
	public:
		std::atomic<uint64_t> current;						// which guide is next to be processed (starts from 0)
		std::vector<uint64_t> guide;							// the vector of guides to perform work on (the things we are looking for)

		std::vector<uint64_t> genome_guides;				// Where we are looking into
		std::vector<uint16_t> genome_guide_frequencies;	// the number of times each guide appears in the genome

		std::mutex file_mutex;									// mutex that controls the file output critical section
		FILE *output_file;										// the FILE handle to the output file (that should be accessed controlled using file_mutex)

		std::mutex stats_mutex;									// mutex that controls the stats (listed here below)
		std::atomic<uint64_t> hits;							// the number of 20-mers that have at least one match
		std::atomic<uint64_t> best_20mer;					// the 20-mer with the highest score so far
		std::atomic<double> best_score;						// the score of the best 20-mer so far

	public:
		/*
			JOB::JOB()
			----------
		*/
		/*!
			@brief Constructor
		*/
		job() :
			current(0),
			output_file(nullptr),
			hits(0),
			best_20mer(0),
			best_score(0)
			{
			/* Nothing */
			}

		/*
			JOB::GET_NEXT()
			---------------
		*/
		/*!
			@brief Get the next value to process
			@returns the next 20-mer to process
		*/
		uint64_t get_next(void)
			{
			uint64_t was;

			do
				was = current;
			while (!std::atomic_compare_exchange_strong(&current, &was, was + 1));

			return was;
			}
	};

